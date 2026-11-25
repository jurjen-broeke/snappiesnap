//! Import every file of the local corpus (`cargo xtask corpus-test`).

#![cfg(feature = "corpus")]
// Test code: the never-crash lints apply to production code only (see clippy.toml).
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use std::path::{Path, PathBuf};

use snappiesnap_doc::{Range, Strand, Topology};
use snappiesnap_io::{Source, import};

fn corpus_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../corpus")
}

fn files(dir: &Path, extension: &str, depth: usize, out: &mut Vec<PathBuf>) {
    assert!(depth < 16);
    for entry in std::fs::read_dir(dir).unwrap() {
        let path = entry.unwrap().path();
        if path
            .file_name()
            .is_some_and(|n| n.to_string_lossy().starts_with('.'))
        {
            continue;
        }
        if path.is_dir() {
            files(&path, extension, depth + 1, out);
        } else if path.extension().is_some_and(|e| e == extension) {
            out.push(path);
        }
    }
}

/// Whether the corpus holds every file a test is pinned to. A corpus is local and partial by
/// nature, so a test whose material is absent says so and skips instead of failing.
fn has(files: &[&str]) -> bool {
    for f in files {
        if !corpus_dir().join(f).exists() {
            println!("no {f} in the corpus: nothing to check");
            return false;
        }
    }
    true
}

fn open(rel: &str) -> snappiesnap_io::Imported {
    let path = corpus_dir().join(rel);
    import(&path.to_string_lossy(), &std::fs::read(&path).unwrap()).unwrap()
}

#[test]
fn every_dna_file_imports_without_losing_features() {
    let mut found = Vec::new();
    files(&corpus_dir().join("dna"), "dna", 0, &mut found);
    if found.is_empty() {
        println!("no .dna files in the corpus: nothing to check");
        return;
    }
    for path in &found {
        let imported = import(&path.to_string_lossy(), &std::fs::read(path).unwrap())
            .unwrap_or_else(|e| panic!("{}: {e}", path.display()));
        let Source::Dna(src) = &imported.source else {
            panic!("{}: not imported as .dna", path.display())
        };
        let in_file = src.file().features().map(|f| f.len()).unwrap_or(0);
        assert_eq!(
            imported.molecule.features().len(),
            in_file,
            "{}: {:?}",
            path.display(),
            imported.warnings
        );
        imported.molecule.validate().unwrap();
        eprintln!(
            "ok  {}: {} bp {:?}, {} features, warnings {:?}",
            path.file_name().unwrap().to_string_lossy(),
            imported.molecule.len(),
            imported.molecule.topology(),
            in_file,
            imported.warnings
        );
    }
}

#[test]
fn coordinates_strands_and_names_are_mapped() {
    if !has(&["dna/edits/mCherry_edit.dna", "dna/roundtrip/pENTR11.dna"]) {
        return;
    }
    let m = open("dna/edits/mCherry_edit.dna").molecule;
    assert_eq!(m.name(), "mCherry_edit");
    assert_eq!(m.topology(), Topology::Linear);
    let start = m
        .features()
        .iter()
        .find(|f| f.name == "start-codon")
        .unwrap();
    assert_eq!(start.segments[0].range, Range::new(0, 3));
    assert_eq!(start.strand, Strand::None);
    assert_eq!(m.residues(Range::new(0, 3)).unwrap(), b"ATG");
    let cds = m.features().iter().find(|f| f.name == "mCherry").unwrap();
    assert_eq!(cds.kind, "CDS");
    assert_eq!(cds.strand, Strand::Forward);
    assert_eq!(cds.segments[0].range, Range::new(0, 711));

    let p = open("dna/roundtrip/pENTR11.dna").molecule;
    assert_eq!(p.topology(), Topology::Circular);
    let lac = p
        .features()
        .iter()
        .find(|f| f.name == "lac UV5 promoter")
        .unwrap();
    assert_eq!(lac.strand, Strand::Reverse);
    // Stored ascending (1926-1932, 1933-1950, 1951-1956); 5'→3' on the reverse strand is
    // descending.
    let starts: Vec<usize> = lac.segments.iter().map(|s| s.range.start).collect();
    assert_eq!(starts, vec![1950, 1932, 1925]);
    assert!(p.notes().fields.contains_key("UUID"));
}

/// Open → save without changes gives back every corpus file byte for byte, through the full
/// import/export path (not just the container).
#[test]
fn unchanged_documents_save_byte_exact() {
    let mut found = Vec::new();
    files(&corpus_dir().join("dna"), "dna", 0, &mut found);
    for path in &found {
        let bytes = std::fs::read(path).unwrap();
        let i = import(&path.to_string_lossy(), &bytes).unwrap();
        let e = snappiesnap_io::export(&i.molecule, Some(&i.source), snappiesnap_io::Format::Dna)
            .unwrap();
        let (out, src) = (e.bytes, e.source);
        assert!(out == bytes, "{}: saved bytes differ", path.display());
        // And again from the source the save returned.
        let again = snappiesnap_io::export(&i.molecule, Some(&src), snappiesnap_io::Format::Dna)
            .unwrap()
            .bytes;
        assert!(again == bytes, "{}: second save differs", path.display());
    }
}

/// Adding the start-codon feature (bp 1-3, misc_feature, no strand) to mCherry.dna produces
/// the same features packet the reference application wrote in mCherry_edit.dna.
#[test]
fn adding_a_feature_matches_the_reference_application() {
    if !has(&["dna/roundtrip/mCherry.dna"]) {
        return;
    }
    use snappiesnap_doc::{NewFeature, Segment};
    use snappiesnap_formats::dna::{DnaFile, tag};

    let mut original = open("dna/roundtrip/mCherry.dna");
    original
        .molecule
        .add_feature(NewFeature {
            color: None,
            name: "start-codon".into(),
            kind: "misc_feature".into(),
            strand: Strand::None,
            segments: vec![Segment {
                range: Range::new(0, 3),
            }],
            qualifiers: vec![],
        })
        .unwrap();
    let ours = snappiesnap_io::export(
        &original.molecule,
        Some(&original.source),
        snappiesnap_io::Format::Dna,
    )
    .unwrap()
    .bytes;
    let reference = std::fs::read(corpus_dir().join("dna/edits/mCherry_edit.dna")).unwrap();
    let ours = DnaFile::parse(&ours).unwrap();
    let reference = DnaFile::parse(&reference).unwrap();
    let ours_text = String::from_utf8_lossy(ours.payload(tag::FEATURES).unwrap()).into_owned();
    let reference_text =
        String::from_utf8_lossy(reference.payload(tag::FEATURES).unwrap()).into_owned();
    // The new feature is written exactly as the reference application wrote it.
    let added = r##"<Feature recentID="1" name="start-codon" type="misc_feature" allowSegmentOverlaps="0" consecutiveTranslationNumbering="1"><Segment range="1-3" color="#a6acb3" type="standard"/></Feature>"##;
    assert!(ours_text.contains(added) && reference_text.contains(added));
    // The only other difference: the newer version that saved the reference file added a
    // derived `translationMW` to the existing CDS (it does that on a plain re-save too). We
    // keep the original feature byte for byte.
    assert_eq!(
        without_attr(&ours_text, "translationMW"),
        without_attr(&reference_text, "translationMW")
    );
    // The sequence packet is untouched by a feature edit, in both files.
    assert_eq!(
        ours.payload(tag::SEQUENCE).unwrap(),
        reference.payload(tag::SEQUENCE).unwrap()
    );
}

/// `text` with every ` name="..."` attribute removed.
fn without_attr(text: &str, name: &str) -> String {
    let needle = format!(" {name}=\"");
    let mut out = String::new();
    let mut rest = text;
    while let Some(i) = rest.find(&needle) {
        out.push_str(&rest[..i]);
        let after = &rest[i + needle.len()..];
        rest = &after[after.find('"').map_or(after.len(), |e| e + 1)..];
    }
    out.push_str(rest);
    out
}

/// Every FASTA file opens and saves back byte for byte; pCas9_GFP's FASTA export has the same
/// sequence as its .dna file.
#[test]
fn fasta_files_round_trip_and_match_their_dna() {
    if !has(&["dna/roundtrip/pCas9_GFP.dna", "fasta/pCas9_GFP_fasta.fa"]) {
        return;
    }
    let mut found = Vec::new();
    files(&corpus_dir().join("fasta"), "fa", 0, &mut found);
    files(&corpus_dir().join("fasta"), "fasta", 0, &mut found);
    assert!(!found.is_empty(), "no FASTA files in corpus/fasta");
    for path in &found {
        let bytes = std::fs::read(path).unwrap();
        let i = import(&path.to_string_lossy(), &bytes).unwrap();
        let out =
            snappiesnap_io::export(&i.molecule, Some(&i.source), snappiesnap_io::Format::Fasta)
                .unwrap();
        assert!(out.bytes == bytes, "{}: saved bytes differ", path.display());
        // The plain parser reproduces SnapGene's export layout too.
        let parsed = snappiesnap_formats::fasta::FastaFile::parse(&bytes).unwrap();
        assert!(
            parsed.to_text().as_bytes() == bytes,
            "{}: layout differs",
            path.display()
        );
    }
    let fasta = open("fasta/pCas9_GFP_fasta.fa").molecule;
    let dna = open("dna/roundtrip/pCas9_GFP.dna").molecule;
    assert_eq!(fasta.seq(), dna.seq());
}

/// (name, type, ranges, strand) per feature, sorted, for comparing files across formats.
fn feature_summary(m: &snappiesnap_doc::Molecule) -> Vec<(String, String, Vec<Range>, Strand)> {
    let mut v: Vec<(String, String, Vec<Range>, Strand)> = m
        .features()
        .iter()
        .map(|f| {
            (
                f.name.clone(),
                f.kind.clone(),
                f.segments.iter().map(|s| s.range).collect(),
                f.strand,
            )
        })
        .collect();
    v.sort_by(|a, b| {
        (&a.0, a.2.first().map(|r| r.start)).cmp(&(&b.0, b.2.first().map(|r| r.start)))
    });
    v
}

/// GenBank exports open byte-exact and match the .dna file they were exported from: same
/// sequence, topology and features. The standard flavour has no direction information for
/// non-coding features, so those may come back without a direction (and only that).
#[test]
fn genbank_files_match_their_dna() {
    if !has(&["dna/roundtrip/pENTR11.dna"]) {
        return;
    }
    let dna = open("dna/roundtrip/pENTR11.dna").molecule;
    for (file, strict) in [
        ("genbank/pENTR11_genbank-snapgene.gbk", true),
        ("genbank/pENTR11_genbank-standard.gb", false),
    ] {
        let bytes = std::fs::read(corpus_dir().join(file)).unwrap();
        let i = import(file, &bytes).unwrap();
        assert!(i.warnings.is_empty(), "{file}: {:?}", i.warnings);
        let out = snappiesnap_io::export(
            &i.molecule,
            Some(&i.source),
            snappiesnap_io::Format::GenBank,
        )
        .unwrap();
        assert!(out.bytes == bytes, "{file}: saved bytes differ");
        let gb = i.molecule;
        assert_eq!(gb.seq(), dna.seq(), "{file}");
        assert_eq!(gb.topology(), dna.topology(), "{file}");
        let (a, b) = (feature_summary(&gb), feature_summary(&dna));
        assert_eq!(a.len(), b.len(), "{file}");
        for (x, y) in a.iter().zip(&b) {
            // The standard flavour can't express direction except for complement(); the
            // reference application imports plain non-coding features without one.
            let same_strand = x.3 == y.3 || (!strict && x.3 == Strand::None);
            assert!(
                x.0 == y.0 && x.1 == y.1 && x.2 == y.2 && same_strand,
                "{file}:\n  genbank {x:?}\n  dna     {y:?}"
            );
        }
    }
}

/// Every .dna file converts to GenBank and back without losing a feature's name, type,
/// ranges or strand.
#[test]
fn dna_converts_to_genbank_and_back() {
    let mut found = Vec::new();
    files(&corpus_dir().join("dna"), "dna", 0, &mut found);
    for path in &found {
        let m = import(&path.to_string_lossy(), &std::fs::read(path).unwrap())
            .unwrap()
            .molecule;
        let gb = snappiesnap_io::export(&m, None, snappiesnap_io::Format::GenBank).unwrap();
        let back = import("x.gb", &gb.bytes).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
        // GenBank has no letter case (some .dna files store lower case), so compare without it.
        assert_eq!(
            back.molecule.seq().as_bytes().to_ascii_uppercase(),
            m.seq().as_bytes().to_ascii_uppercase(),
            "{}",
            path.display()
        );
        assert_eq!(back.molecule.topology(), m.topology());
        assert_eq!(
            feature_summary(&back.molecule),
            feature_summary(&m),
            "{}",
            path.display()
        );
    }
}

/// A GenBank file written by a third-party tool imports the way the reference application
/// imports it: same sequence and features, with the same strands. Two known differences: the
/// reference `.dna` has one feature added after import, and the reference application marks
/// "ori" forward from its own knowledge of common features (the file gives no direction for it;
/// "f1 ori", written the same way, stays without one there too).
#[test]
fn third_party_genbank_imports_like_the_reference_application() {
    if !has(&[
        "dna/edits/third-party-import.dna",
        "genbank/third-party-import.gb",
    ]) {
        return;
    }
    let gb = open("genbank/third-party-import.gb").molecule;
    let dna = open("dna/edits/third-party-import.dna").molecule;
    assert_eq!(
        gb.seq().as_bytes().to_ascii_uppercase(),
        dna.seq().as_bytes().to_ascii_uppercase()
    );
    let ours: Vec<_> = feature_summary(&gb)
        .into_iter()
        .map(|mut f| {
            if f.0 == "ori" {
                f.3 = Strand::Forward;
            }
            f
        })
        .collect();
    let theirs: Vec<_> = feature_summary(&dna)
        .into_iter()
        .filter(|f| f.0 != "This base could be missing")
        .collect();
    assert_eq!(ours, theirs);
}

/// Every CDS with a stored /translation translates to exactly that: segments joined in order,
/// on the feature's strand, from /codon_start, with /transl_table (default 1), the first codon
/// read as a start and the final stop left out.
#[test]
fn cds_translations_match_stored_ones() {
    use snappiesnap_bio::translate::{code, strand_residues, translate};
    let mut files = Vec::new();
    files_with_ext(&corpus_dir().join("genbank"), &mut files);
    let mut checked = 0;
    for path in files {
        let i = import(&path.to_string_lossy(), &std::fs::read(&path).unwrap()).unwrap();
        let m = i.molecule;
        for f in m.features().iter().filter(|f| f.kind == "CDS") {
            let q = |n: &str| {
                f.qualifiers
                    .iter()
                    .find(|q| q.name == n)
                    .and_then(|q| q.values.first())
                    .map(|v| v.replace([' ', '\n'], ""))
            };
            let Some(expected) = q("translation") else {
                continue;
            };
            let start: usize = q("codon_start").and_then(|v| v.parse().ok()).unwrap_or(1);
            let table: u8 = q("transl_table").and_then(|v| v.parse().ok()).unwrap_or(1);
            let mut nt = Vec::new();
            for seg in &f.segments {
                nt.extend(strand_residues(&m, seg.range, f.strand).unwrap());
            }
            // Joins on the reverse strand run 3'->5' along the top strand, segment by segment
            // in feature order, so each piece is already in reading order.
            let nt = nt.get(start.saturating_sub(1)..).unwrap_or_default();
            let protein = translate(nt, code(table).unwrap(), true);
            let protein = String::from_utf8(protein).unwrap();
            let protein = protein.strip_suffix('*').unwrap_or(&protein);
            assert_eq!(protein, expected, "{} {}", path.display(), f.name);
            checked += 1;
        }
    }
    println!("{checked} stored translations checked");
    if checked == 0 {
        println!("no CDS with a stored translation in the corpus: nothing to check");
    }
}

fn files_with_ext(dir: &Path, out: &mut Vec<PathBuf>) {
    for e in std::fs::read_dir(dir).unwrap() {
        let p = e.unwrap().path();
        if p.extension().is_some_and(|x| x == "gb" || x == "gbk") {
            out.push(p);
        }
    }
}

/// The same INSDC record as EMBL (ENA) and GenBank (NCBI) opens as the same document, saves
/// back byte for byte, and an edit rewrites only the edited feature.
#[test]
fn embl_matches_genbank_and_saves_minimal_edits() {
    if !has(&["embl/J01749-pBR322.embl", "genbank-ncbi/J01749-pBR322.gb"]) {
        return;
    }
    use snappiesnap_io::{Format, export};
    let embl = open("embl/J01749-pBR322.embl");
    let gb = open("genbank-ncbi/J01749-pBR322.gb");
    assert_eq!(embl.format, Format::Embl);
    assert_eq!(embl.molecule.seq(), gb.molecule.seq());
    assert_eq!(embl.molecule.topology(), Topology::Circular);
    assert_eq!(
        feature_summary(&embl.molecule),
        feature_summary(&gb.molecule)
    );
    // Unchanged: the original bytes.
    let original = std::fs::read(corpus_dir().join("embl/J01749-pBR322.embl")).unwrap();
    let same = export(&embl.molecule, Some(&embl.source), Format::Embl).unwrap();
    assert!(same.bytes == original);
    // Rename one feature: only its entry changes.
    let mut m = embl.molecule.clone();
    let id = m.features().iter().find(|f| f.kind == "CDS").unwrap().id;
    let mut f = m.feature(id).unwrap().clone();
    f.name = "renamed".into();
    m.update_feature(f).unwrap();
    let out = export(&m, Some(&embl.source), Format::Embl).unwrap();
    let text = String::from_utf8(out.bytes.clone()).unwrap();
    let before: Vec<&str> = std::str::from_utf8(&original).unwrap().lines().collect();
    let after: Vec<&str> = text.lines().collect();
    let changed: Vec<&&str> = after.iter().filter(|l| !before.contains(l)).collect();
    assert_eq!(
        changed,
        vec![&"FT                   /label=renamed"],
        "{changed:?}"
    );
    let again = import("x.embl", &out.bytes).unwrap();
    assert!(
        again
            .molecule
            .features()
            .iter()
            .any(|f| f.name == "renamed")
    );
    // Converted to and from GenBank, features survive.
    let as_gb = export(&embl.molecule, None, Format::GenBank).unwrap();
    let back = import("x.gb", &as_gb.bytes).unwrap();
    let as_embl = export(&back.molecule, None, Format::Embl).unwrap();
    let round = import("x.embl", &as_embl.bytes).unwrap();
    assert_eq!(round.molecule.seq(), embl.molecule.seq());
    assert_eq!(
        feature_summary(&round.molecule).len(),
        feature_summary(&embl.molecule).len()
    );
}

/// NCBI's GFF3 of a record agrees with its GenBank file: same sequence and topology, and every
/// CDS at the same place on the same strand. The annotations-only GFF3 gives the same features
/// as the one with the sequence; unchanged, it saves byte for byte.
#[test]
fn gff3_matches_genbank() {
    if !has(&[
        "genbank-ncbi/J01749-pBR322.gb",
        "gff3/J01749-pBR322-with-fasta.gff3",
    ]) {
        return;
    }
    use snappiesnap_io::{Format, export, gff3_annotations};
    let gff = open("gff3/J01749-pBR322-with-fasta.gff3");
    let gb = open("genbank-ncbi/J01749-pBR322.gb");
    assert_eq!(gff.format, Format::Gff3);
    assert_eq!(gff.molecule.seq(), gb.molecule.seq());
    assert_eq!(gff.molecule.topology(), Topology::Circular);
    let cds = |m: &snappiesnap_doc::Molecule| {
        let mut v: Vec<(Vec<Range>, Strand)> = m
            .features()
            .iter()
            .filter(|f| f.kind == "CDS")
            .map(|f| (f.segments.iter().map(|s| s.range).collect(), f.strand))
            .collect();
        v.sort_by_key(|(r, _)| r.first().map(|r| r.start));
        v
    };
    assert!(!cds(&gb.molecule).is_empty());
    assert_eq!(cds(&gff.molecule), cds(&gb.molecule));
    let bytes = std::fs::read(corpus_dir().join("gff3/J01749-pBR322.gff3")).unwrap();
    let (seqid, features, warnings) =
        gff3_annotations(&bytes, None, gb.molecule.len(), Topology::Circular).unwrap();
    assert_eq!(seqid, "J01749.1");
    assert!(warnings.is_empty(), "{warnings:?}");
    assert_eq!(features.len(), gff.molecule.features().len());
    let original = std::fs::read(corpus_dir().join("gff3/J01749-pBR322-with-fasta.gff3")).unwrap();
    assert!(
        export(&gff.molecule, Some(&gff.source), Format::Gff3)
            .unwrap()
            .bytes
            == original
    );
    eprintln!(
        "gff3: {} features, genbank: {} features, {} CDS",
        gff.molecule.features().len(),
        gb.molecule.features().len(),
        cds(&gb.molecule).len()
    );
}

/// SBOL 2 files from the SBOL test suite: designs with a sequence open, save back byte for
/// byte, and an edit changes only the edited annotation; files without a sequence are a clear
/// error.
#[test]
fn sbol_files_open_and_save_minimal_edits() {
    use snappiesnap_io::{Format, export};
    let mut found = Vec::new();
    files(&corpus_dir().join("sbol"), "xml", 0, &mut found);
    found.sort();
    if found.is_empty() {
        println!("no SBOL files in the corpus: nothing to check");
        return;
    }
    for path in &found {
        let bytes = std::fs::read(path).unwrap();
        let rel = path
            .strip_prefix(corpus_dir())
            .unwrap()
            .display()
            .to_string();
        let i = match import(&path.to_string_lossy(), &bytes) {
            Ok(i) => i,
            Err(e) => {
                // A clear refusal: no sequence to show, or residues a sequence can't hold
                // (gaps); never a silently altered sequence.
                let msg = e.to_string();
                assert!(
                    msg.contains("no sequence") || msg.contains("invalid"),
                    "{rel}: {e}"
                );
                eprintln!("ok  {rel}: refused ({e})");
                continue;
            }
        };
        assert_eq!(i.format, Format::Sbol);
        let same = export(&i.molecule, Some(&i.source), Format::Sbol).unwrap();
        assert!(same.bytes == bytes, "{rel}: unchanged save differs");
        if let Some(f) = i.molecule.features().first() {
            let mut m = i.molecule.clone();
            let mut f = f.clone();
            f.name = "renamed".into();
            m.update_feature(f).unwrap();
            let out = export(&m, Some(&i.source), Format::Sbol).unwrap();
            let text = String::from_utf8(out.bytes.clone()).unwrap();
            let original = String::from_utf8(bytes.clone()).unwrap();
            // Everything before the first changed byte and after the last is untouched.
            let prefix = original
                .bytes()
                .zip(text.bytes())
                .take_while(|(a, b)| a == b)
                .count();
            let suffix = original
                .bytes()
                .rev()
                .zip(text.bytes().rev())
                .take_while(|(a, b)| a == b)
                .count();
            let changed = original.len().saturating_sub(prefix + suffix);
            assert!(
                changed < 1500,
                "{rel}: {changed} bytes changed for one rename"
            );
            let again = import("x.xml", &out.bytes).unwrap();
            assert!(
                again
                    .molecule
                    .features()
                    .iter()
                    .any(|f| f.name == "renamed"),
                "{rel}"
            );
            assert_eq!(again.molecule.seq(), i.molecule.seq(), "{rel}");
            // Saving again from the written source changes nothing.
            assert!(
                export(&m, Some(&out.source), Format::Sbol).unwrap().bytes == out.bytes,
                "{rel}"
            );
        }
        // Written as a new document, it reads back the same.
        let fresh = export(&i.molecule, None, Format::Sbol).unwrap();
        let back = import("x.xml", &fresh.bytes).unwrap();
        assert_eq!(back.molecule.seq(), i.molecule.seq(), "{rel}");
        assert_eq!(
            back.molecule.features().len(),
            i.molecule.features().len(),
            "{rel}"
        );
        eprintln!(
            "ok  {rel}: {} bp {:?}, {} features, warnings {:?}",
            i.molecule.len(),
            i.molecule.topology(),
            i.molecule.features().len(),
            i.warnings
        );
    }
}

/// Matched reads of `.dna` files become molecules: traces with their chromatogram, plain
/// sequences without, all within their stored trimming.
#[test]
fn matched_reads_become_molecules() {
    let mut paths = Vec::new();
    files(&corpus_dir(), "dna", 0, &mut paths);
    let (mut traces, mut plain) = (0, 0);
    for path in &paths {
        let bytes = std::fs::read(path).unwrap();
        let i = import(&path.to_string_lossy(), &bytes).unwrap();
        let Source::Dna(src) = &i.source else {
            continue;
        };
        for r in src.matched_reads() {
            let m = r.molecule.unwrap_or_else(|e| panic!("{}: {e}", r.name));
            match m.trace() {
                Some(t) => {
                    traces += 1;
                    assert_eq!(t.calls().len(), m.len());
                }
                None => plain += 1,
            }
            if let Some((a, b)) = r.trimmed {
                assert!(a < b && b <= m.len());
            }
        }
    }
    eprintln!("matched reads: {traces} traces, {plain} plain sequences");
    if traces == 0 && plain == 0 {
        println!("no .dna file in the corpus carries matched reads: nothing to check");
        return;
    }
    assert!(traces >= 6 && plain >= 1);
}

/// Homology assembly (In-Fusion/HiFi) against a reaction that was actually carried out.
///
/// The corpus holds an entry clone made by In-Fusion from a Gateway entry vector, cut once, and
/// a PCR fragment that its own construction history keeps. Cutting the vector where the product
/// says it was opened and assembling the two must give the product back: same length, two 15 bp
/// junctions, same sequence.
///
/// One base is allowed to differ. An entry vector may carry a single-base variant inside attL1
/// relative to the published site, in which case the product holds a base its own recorded
/// vector does not. Everything else must match exactly.
#[test]
fn infusion_reproduces_an_entry_clone_that_was_built() {
    use snappiesnap_bio::assembly::homology_assembly;
    use snappiesnap_bio::cloning::digest_fragments;
    use snappiesnap_bio::enzymes::builtin;
    use snappiesnap_formats::dna::DnaFile;

    const VECTOR: &str = "genbank/entry-vector.gb";
    const PRODUCT: &str = "dna/entry-clone-infusion.dna";
    if !has(&[VECTOR, PRODUCT]) {
        return;
    }
    let vector = open(VECTOR).molecule;
    let product = open(PRODUCT).molecule;

    // The insert: the PCR fragment the product's own history keeps.
    let bytes = std::fs::read(corpus_dir().join(PRODUCT)).unwrap();
    let file = DnaFile::parse(&bytes).unwrap();
    let history = file
        .history()
        .unwrap()
        .expect("the product records its history");
    let insert = history
        .tree
        .walk()
        .into_iter()
        .filter_map(|n| history.record(n.id)?.residues.clone())
        .find(|r| r.len() == 777)
        .expect("the 777 bp PCR fragment is kept with the history");

    // The vector, opened where the product says it was opened: HindIII cuts it once.
    let hindiii: Vec<_> = builtin()
        .into_iter()
        .filter(|e| e.name == "HindIII")
        .collect();
    assert_eq!(hindiii.len(), 1, "HindIII is in the built-in enzymes");
    let fragments = digest_fragments(&vector, &hindiii);
    assert_eq!(fragments.len(), 1, "HindIII cuts the entry vector once");
    let linear = fragments[0].part.clone();

    let parts = [
        linear,
        snappiesnap_bio::cloning::Part {
            seq: insert,
            features: Vec::new(),
        },
    ];
    let assembly = homology_assembly(&parts, 15, true).expect("the two pieces join");
    assert_eq!(
        assembly.junctions.len(),
        2,
        "a vector and one insert make two junctions"
    );
    for j in &assembly.junctions {
        assert_eq!(j.length, 15, "In-Fusion arms of 15 bp");
    }

    // The same sequence, wherever the product's own origin happens to sit.
    let ours = assembly.part.seq.to_ascii_uppercase();
    let theirs = product.seq().as_bytes().to_ascii_uppercase();
    assert_eq!(ours.len(), theirs.len(), "the product's length");
    let doubled = [ours.as_slice(), ours.as_slice()].concat();
    let at = doubled
        .windows(theirs.len().min(40))
        .position(|w| w == &theirs[..theirs.len().min(40)])
        .expect("the product starts somewhere in ours");
    let rotated: Vec<u8> = doubled[at..at + theirs.len()].to_vec();
    let differences = rotated.iter().zip(&theirs).filter(|(a, b)| a != b).count();
    println!(
        "In-Fusion: {} bp product, {} junction(s), {differences} base(s) differing from the recorded product",
        ours.len(),
        assembly.junctions.len()
    );
    assert!(
        differences <= 1,
        "{differences} bases differ from the recorded product"
    );
}

/// Gateway LR against a reaction that was actually carried out.
///
/// The corpus holds the expression clone of an LR between a destination vector and an entry
/// clone. The destination vector is kept in the product's own construction history, so both
/// inputs are to hand: recombining them must give the product back exactly, with attB1 and
/// attB2 where the att sites were.
#[test]
fn gateway_lr_reproduces_an_expression_clone_that_was_built() {
    use snappiesnap_bio::cloning::{GatewayReaction, gateway};
    use snappiesnap_doc::{Alphabet, Molecule, Seq};

    const PRODUCT: &str = "dna/expression-clone-lr.dna";
    const ENTRY: &str = "dna/entry-clone-tagged.dna";
    if !has(&[PRODUCT, ENTRY]) {
        return;
    }
    let product = open(PRODUCT).molecule;
    let entry = open(ENTRY).molecule;

    // The destination vector: the input of the last step that is not the entry clone.
    let bytes = std::fs::read(corpus_dir().join(PRODUCT)).unwrap();
    let file = snappiesnap_formats::dna::DnaFile::parse(&bytes).unwrap();
    let history = file
        .history()
        .unwrap()
        .expect("the product records its history");
    let destination = history
        .tree
        .walk()
        .into_iter()
        .filter(|n| n.length != entry.len() && n.length > entry.len())
        .filter_map(|n| history.record(n.id)?.residues.clone())
        .map(|r| {
            Molecule::new(
                "destination",
                Seq::new(Alphabet::Dna, r).unwrap(),
                Topology::Circular,
            )
            .unwrap()
        })
        .find(|m| gateway(m, &entry, GatewayReaction::Lr, "p").is_ok())
        .expect("a destination vector in the history recombines with the entry clone");

    let ours = gateway(&destination, &entry, GatewayReaction::Lr, "product").unwrap();
    assert_eq!(ours.len(), product.len(), "the expression clone's length");
    let names: Vec<&str> = ours.features().iter().map(|f| f.name.as_str()).collect();
    assert!(
        names.contains(&"attB1") && names.contains(&"attB2"),
        "the product carries attB1 and attB2: {names:?}"
    );
    // The same sequence, wherever the product's own origin sits.
    let (a, b) = (
        ours.seq().as_bytes().to_ascii_uppercase(),
        product.seq().as_bytes().to_ascii_uppercase(),
    );
    let doubled = [a.as_slice(), a.as_slice()].concat();
    let at = doubled
        .windows(40)
        .position(|w| w == &b[..40])
        .expect("the product starts somewhere in ours");
    assert_eq!(
        &doubled[at..at + b.len()],
        b.as_slice(),
        "the LR product differs from the recorded one"
    );
    println!(
        "Gateway LR: {} bp, identical to the recorded clone",
        ours.len()
    );
}

/// PCR against a reaction that was actually carried out, with tailed oligos.
///
/// `corpus/oligo/pcr-forward.txt` and `pcr-reverse.txt` are the primers of a PCR on one corpus
/// construct whose product went into another; that product is kept in the second construct's
/// construction history. Priming with the two oligo files must give exactly it, tails and all.
#[test]
fn pcr_with_tailed_oligos_matches_its_recorded_product() {
    use snappiesnap_bio::pcr::{DEFAULT_MAX_PRODUCT, pcr};

    const TEMPLATE: &str = "dna/entry-clone-template.dna";
    const CARRIER: &str = "dna/entry-clone-tagged.dna";
    const FORWARD: &str = "oligo/pcr-forward.txt";
    const REVERSE: &str = "oligo/pcr-reverse.txt";
    if !has(&[TEMPLATE, CARRIER, FORWARD, REVERSE]) {
        return;
    }
    let template = open(TEMPLATE).molecule;
    let oligo = |rel: &str| -> Vec<u8> {
        std::fs::read_to_string(corpus_dir().join(rel))
            .unwrap()
            .bytes()
            .filter(u8::is_ascii_alphabetic)
            .map(|b| b.to_ascii_uppercase())
            .collect()
    };
    let (forward, reverse) = (oligo(FORWARD), oligo(REVERSE));

    let products = pcr(&template, [&forward, &reverse], None, DEFAULT_MAX_PRODUCT);
    assert_eq!(products.len(), 1, "one product from these two oligos");
    let ours = products[0].seq.to_ascii_uppercase();

    // The fragment kept with the construct it went into.
    let bytes = std::fs::read(corpus_dir().join(CARRIER)).unwrap();
    let file = snappiesnap_formats::dna::DnaFile::parse(&bytes).unwrap();
    let history = file
        .history()
        .unwrap()
        .expect("the construct records its history");
    let theirs = history
        .tree
        .walk()
        .into_iter()
        .filter_map(|n| history.record(n.id)?.residues.clone())
        .find(|r| r.len() == ours.len())
        .expect("the PCR fragment is kept with the history");
    assert_eq!(
        ours,
        theirs.to_ascii_uppercase(),
        "the PCR product differs from the fragment kept with the construct"
    );
    println!("PCR: {} bp, identical to the recorded fragment", ours.len());
}
