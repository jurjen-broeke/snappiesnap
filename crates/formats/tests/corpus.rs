//! Byte-exact round trips over the local corpus (`cargo xtask corpus-test`).
//!
//! The corpus lives in `corpus/` at the repository root. It is local-only and never committed
//! (AGENTS.md › File fidelity), so these tests only build with the `corpus` feature.

#![cfg(feature = "corpus")]
// Test code: the never-crash lints apply to production code only (see clippy.toml).
#![allow(clippy::unwrap_used, clippy::expect_used, clippy::panic)]

use std::path::{Path, PathBuf};

use snappiesnap_formats::dna::{DnaFile, tag};

fn corpus_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../corpus")
}

/// Every file with `extension` under `dir`, skipping hidden folders (such as SnapGene's
/// `.sglock`).
fn files(dir: &Path, extension: &str, depth: usize, out: &mut Vec<PathBuf>) {
    assert!(depth < 16, "corpus nested too deeply at {}", dir.display());
    for entry in std::fs::read_dir(dir).unwrap() {
        let path = entry.unwrap().path();
        let hidden = path
            .file_name()
            .is_some_and(|n| n.to_string_lossy().starts_with('.'));
        if hidden {
            continue;
        }
        if path.is_dir() {
            files(&path, extension, depth + 1, out);
        } else if path.extension().is_some_and(|e| e == extension) {
            out.push(path);
        }
    }
}

#[test]
fn every_dna_file_round_trips_byte_for_byte() {
    let dir = corpus_dir().join("dna");
    assert!(
        dir.is_dir(),
        "{} is missing: the corpus is local-only, see corpus/README.md",
        dir.display()
    );
    let mut found = Vec::new();
    files(&dir, "dna", 0, &mut found);
    found.sort();
    assert!(!found.is_empty(), "no .dna files under {}", dir.display());
    for path in &found {
        let bytes = std::fs::read(path).unwrap();
        let file = DnaFile::parse(&bytes).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
        let written = file.to_bytes().unwrap();
        assert!(
            written == bytes,
            "{}: written bytes differ from the original",
            path.display()
        );
        eprintln!(
            "ok  {} ({} bytes, {} packets)",
            path.strip_prefix(corpus_dir()).unwrap_or(path).display(),
            bytes.len(),
            file.packets().len()
        );
    }
}

/// Re-serialising every feature from its parsed form (not its raw text) reproduces the
/// reference application's bytes: this pins down attribute order and escaping.
#[test]
fn every_feature_reserialises_exactly() {
    let mut found = Vec::new();
    files(&corpus_dir().join("dna"), "dna", 0, &mut found);
    let mut checked = 0usize;
    for path in &found {
        let file = DnaFile::parse(&std::fs::read(path).unwrap()).unwrap();
        // Some files have no features packet at all (e.g. a FASTA imported and saved).
        let Ok(payload) = file.payload(tag::FEATURES) else {
            continue;
        };
        let payload = String::from_utf8(payload.to_vec()).unwrap();
        let features = file
            .features()
            .unwrap_or_else(|e| panic!("{}: {e}", path.display()));
        assert_eq!(features.to_text(), payload, "{}", path.display());
        for (i, f) in features.features().enumerate() {
            let xml = f.to_xml();
            assert!(
                payload.contains(&xml),
                "{} feature {i} ({}) re-serialises differently:\n{xml}",
                path.display(),
                f.name()
            );
            f.segments()
                .unwrap_or_else(|e| panic!("{} feature {i}: {e}", path.display()));
            checked += 1;
        }
    }
    eprintln!("{checked} features re-serialised exactly");
    assert!(checked > 0);
}

/// Each file's cookie, sequence and notes parse, and match the expectations recorded in
/// corpus/MANIFEST.md (`len=`, `circular`/`linear`, `features=`, `primers=`).
#[test]
fn files_match_their_manifest_expectations() {
    let manifest = std::fs::read_to_string(corpus_dir().join("MANIFEST.md")).unwrap();
    let mut rows = 0usize;
    for line in manifest.lines() {
        let cells: Vec<&str> = line.split('|').map(str::trim).collect();
        let (Some(path), Some(expect)) = (cells.get(1), cells.get(6)) else {
            continue;
        };
        if !path.ends_with(".dna") {
            continue;
        }
        let full = corpus_dir().join(path);
        let file = DnaFile::parse(&std::fs::read(&full).unwrap()).unwrap();
        let cookie = file.cookie().unwrap();
        assert_eq!(cookie.sequence_type, 1, "{path}: not DNA");
        let seq = file.sequence().unwrap();
        assert!(
            seq.residues.iter().all(u8::is_ascii_alphabetic),
            "{path}: non-letter residues"
        );
        let notes = file.notes().unwrap();
        assert!(!notes.is_empty(), "{path}: no notes fields");
        for token in expect.split_whitespace() {
            match token.split_once('=') {
                Some(("len", n)) => {
                    assert_eq!(seq.residues.len().to_string(), n, "{path}: length")
                }
                Some(("features", n)) => assert_eq!(
                    file.features().unwrap().len().to_string(),
                    n,
                    "{path}: features"
                ),
                None if token == "circular" => assert!(seq.is_circular(), "{path}"),
                None if token == "linear" => assert!(!seq.is_circular(), "{path}"),
                _ => {}
            }
        }
        rows += 1;
    }
    eprintln!("{rows} manifest rows checked");
    assert!(rows > 0, "MANIFEST.md has no .dna rows");
}

/// GenBank exports: every file writes back byte for byte, and every feature entry, LOCUS line
/// and ORIGIN block re-serialises exactly as the reference application wrote it.
#[test]
fn genbank_files_reserialise_exactly() {
    use snappiesnap_formats::genbank::{GenBankFile, origin_text};
    let mut found = Vec::new();
    for ext in ["gb", "gbk"] {
        files(&corpus_dir().join("genbank"), ext, 0, &mut found);
    }
    if found.is_empty() {
        println!("no GenBank files in corpus/genbank: nothing to check");
        return;
    }
    let mut checked = 0usize;
    let mut mismatches = Vec::new();
    for path in &found {
        let bytes = std::fs::read(path).unwrap();
        let file = GenBankFile::parse(&bytes).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
        assert!(
            file.to_text().as_bytes() == bytes,
            "{}: text differs",
            path.display()
        );
        // The whole file is compared byte for byte just above. The pieces below are compared
        // with the line endings normalised: a file written with CRLF (older vector records
        // are) keeps them in `raw`, while the writers here join with newlines.
        let text = String::from_utf8(bytes).unwrap().replace("\r\n", "\n");
        for r in &file.records {
            let first = text.lines().next().unwrap_or("");
            // Our layout is checked for standard-width names; longer names are laid out
            // differently by different writers, and edits patch those lines in place.
            if r.locus.name.len() <= 16 && r.locus.to_line() != first {
                mismatches.push(format!(
                    "LOCUS:\n  ours: {}\n  file: {first}",
                    r.locus.to_line()
                ));
            }
            // Some writers pad each ORIGIN line with a trailing space; the file is still
            // written back byte for byte (asserted above), so compare the layout without it.
            let origin: String = text
                .lines()
                .map(|l| format!("{}\n", l.trim_end()))
                .collect();
            if !origin.contains(&origin_text(&r.residues)) {
                mismatches.push(format!("{}: ORIGIN block differs", path.display()));
            }
            for f in &r.features {
                checked += 1;
                let raw = f.raw.as_ref().map(|r| r.replace("\r\n", "\n"));
                if Some(f.entry.to_text()) != raw {
                    mismatches.push(format!(
                        "{} {}:\n--- ours\n{}--- file\n{}",
                        path.display(),
                        f.entry.key,
                        f.entry.to_text(),
                        raw.unwrap_or_default()
                    ));
                }
            }
        }
    }
    eprintln!("{checked} GenBank features checked");
    assert!(mismatches.is_empty(), "{}", mismatches.join("\n"));
}

/// Ids of `node` and its inputs that are hidden or feed only into a hidden step.
fn hidden_ids(
    node: &snappiesnap_formats::dna::history::HistoryNode,
    under: bool,
    out: &mut Vec<u32>,
) {
    let here = under || node.hidden;
    if here {
        out.push(node.id);
    }
    for input in &node.inputs {
        hidden_ids(input, here, out);
    }
}

/// Construction histories: every node's residues are recovered at its declared length (stored,
/// or rebuilt by undoing edits from the construct down), every recorded insertion is found in
/// the rebuilt product at its recorded position, and stored ancestors that are also corpus
/// files match them exactly. No file is named: the corpus is searched.
#[test]
fn histories_rebuild_every_node() {
    use snappiesnap_formats::dna::history::Action;
    let mut paths = Vec::new();
    files(&corpus_dir(), "dna", 0, &mut paths);
    let parsed: Vec<DnaFile> = paths
        .iter()
        .map(|p| DnaFile::parse(&std::fs::read(p).unwrap()).unwrap())
        .collect();
    let corpus_seqs: Vec<Vec<u8>> = parsed
        .iter()
        .filter_map(|f| f.sequence().ok().map(|s| s.residues.to_ascii_uppercase()))
        .collect();
    let (mut files_with, mut nodes, mut stored, mut rebuilt, mut inserts, mut matched) =
        (0, 0, 0, 0, 0, 0);
    let (mut hidden, mut unavailable) = (0, 0);
    for f in &parsed {
        let Some(h) = f.history().unwrap() else {
            continue;
        };
        files_with += 1;
        let root = f.sequence().unwrap().residues;
        // Hidden steps, and the steps that only feed them, may be kept without a record:
        // nothing to rebuild them from.
        let mut under_hidden = Vec::new();
        hidden_ids(&h.tree, false, &mut under_hidden);
        for node in h.tree.walk() {
            if under_hidden.contains(&node.id) && h.record(node.id).is_none() {
                hidden += 1;
                continue;
            }
            // A step the file keeps neither the sequence nor an edit for is one the reference
            // application itself marks unavailable: the chain of its ancestors bottoms out in
            // a sequence nobody stored. Nothing can rebuild it, so it is counted, not failed.
            let Some(r) = h.residues(node.id, &root) else {
                // A record is not always enough to walk backwards: these files record what a
                // step inserted (`redo`) but not what it replaced (`undo` inserts nothing), and
                // the last step's edit may not be recorded at all, so the chain cannot be run
                // back from the final sequence. The reference application marks such ancestors
                // unavailable too. A step whose sequence is stored must always be readable.
                assert!(
                    h.record(node.id).is_none_or(|x| x.residues.is_none()),
                    "node {} ({}) has its sequence stored but can't be read",
                    node.id,
                    node.length
                );
                unavailable += 1;
                continue;
            };
            nodes += 1;
            assert_eq!(r.len(), node.length);
            if h.record(node.id).is_some_and(|x| x.residues.is_some()) {
                stored += 1;
                if corpus_seqs.contains(&r.to_ascii_uppercase()) {
                    matched += 1;
                }
            } else {
                rebuilt += 1;
            }
            // A node made by an insertion: its product holds the inserted residues there.
            if let Some(edit) = h.record(node.id).and_then(|x| x.edit.as_ref())
                && let Some((_, Some(product))) = h.tree.find(node.id)
            {
                let p = h.residues(product.id, &root).unwrap();
                for a in &edit.redo {
                    if let Action::Insert { position, residues } = a {
                        let at = p.get(*position..*position + residues.len()).unwrap();
                        assert!(at.eq_ignore_ascii_case(residues), "insert at {position}");
                        inserts += 1;
                    }
                }
            }
        }
    }
    println!(
        "histories: {files_with} files, {nodes} nodes ({stored} stored, {rebuilt} rebuilt; {hidden} hidden steps without a record, {unavailable} unavailable), {inserts} insertions checked, {matched} stored nodes equal corpus files"
    );
    assert!(files_with > 0 && matched > 0);
}

#[test]
fn every_embl_file_round_trips_byte_for_byte() {
    use snappiesnap_formats::embl::EmblFile;
    let mut found = Vec::new();
    let dir = corpus_dir().join("embl");
    if dir.is_dir() {
        files(&dir, "embl", 0, &mut found);
    }
    found.sort();
    for path in &found {
        let bytes = std::fs::read(path).unwrap();
        let file = EmblFile::parse(&bytes).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
        assert!(
            file.to_text().as_bytes() == bytes,
            "{}: written bytes differ from the original",
            path.display()
        );
        let r = &file.records[0];
        assert_eq!(r.residues.len() as u64, r.id.length, "{}", path.display());
        // Every feature re-serialised from its parsed form reads back the same.
        for f in &r.features {
            let text = snappiesnap_formats::embl::ft_text(&f.entry);
            let again = EmblFile::parse(
                format!("ID   X; SV 1; linear; DNA; STD; SYN; 0 BP.\nFH   Key\n{text}//\n")
                    .as_bytes(),
            )
            .unwrap();
            assert_eq!(
                again.records[0].features[0].entry,
                f.entry,
                "{}",
                path.display()
            );
        }
        eprintln!(
            "ok  {} ({} bytes, {} features)",
            path.strip_prefix(corpus_dir()).unwrap_or(path).display(),
            bytes.len(),
            r.features.len()
        );
    }
}

#[test]
fn every_gff3_file_round_trips_byte_for_byte() {
    use snappiesnap_formats::gff3::Gff3File;
    let mut found = Vec::new();
    let dir = corpus_dir().join("gff3");
    if dir.is_dir() {
        files(&dir, "gff3", 0, &mut found);
    }
    found.sort();
    for path in &found {
        let bytes = std::fs::read(path).unwrap();
        let file = Gff3File::parse(&bytes).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
        assert!(
            file.to_text().as_bytes() == bytes,
            "{}: written bytes differ from the original",
            path.display()
        );
        eprintln!(
            "ok  {} ({} rows, {} sequences)",
            path.strip_prefix(corpus_dir()).unwrap_or(path).display(),
            file.lines.len(),
            file.sequences.len()
        );
    }
}

#[test]
fn sbol_files_parse() {
    use snappiesnap_formats::sbol::SbolDocument;
    let mut found = Vec::new();
    let dir = corpus_dir().join("sbol");
    if dir.is_dir() {
        files(&dir, "xml", 0, &mut found);
    }
    found.sort();
    for path in &found {
        let bytes = std::fs::read(path).unwrap();
        let d = SbolDocument::parse(&bytes).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
        let main = d.main_definition();
        eprintln!(
            "ok  {}: {} definitions, {} sequences; shown: {:?} ({} bp, circular {}, {} annotations, {} locations)",
            path.strip_prefix(corpus_dir()).unwrap_or(path).display(),
            d.component_definitions.len(),
            d.sequences.len(),
            main.and_then(|m| m.display_id.clone()),
            main.and_then(|m| m.sequences.first())
                .and_then(|s| d.sequence(s))
                .map_or(0, |s| s.elements.len()),
            main.is_some_and(|m| m.circular()),
            main.map_or(0, |m| m.annotations.len()),
            main.map_or(0, |m| m
                .annotations
                .iter()
                .map(|a| a.locations.len())
                .sum::<usize>()),
        );
    }
}

/// Traces read the same as Biopython 1.88 reads them (`corpus/abif/biopython-oracle.json`,
/// made by a scratch script with `SeqIO.read(path, "abi")`): bases (Biopython reads the
/// basecaller's calls, we prefer the user-edited ones, so case may differ), qualities, peak
/// positions and the sum of each channel.
#[test]
fn abif_traces_match_biopython() {
    let dir = corpus_dir().join("abif");
    let Ok(bytes) = std::fs::read(dir.join("biopython-oracle.json")) else {
        println!("no corpus/abif/biopython-oracle.json: nothing to check");
        return;
    };
    let Ok(oracle) = serde_json::from_slice::<serde_json::Value>(&bytes) else {
        panic!("corpus/abif/biopython-oracle.json is not JSON");
    };
    let mut checked = 0;
    for (name, want) in oracle.as_object().unwrap() {
        let bytes = std::fs::read(dir.join(name)).unwrap();
        let got = snappiesnap_formats::abif::read_chromatogram(&bytes);
        let seq = want["seq"].as_str().unwrap_or("");
        if seq.is_empty() {
            assert!(got.is_err(), "{name}: no base calls, expected an error");
            continue;
        }
        let got = got.unwrap_or_else(|e| panic!("{name}: {e}"));
        assert!(
            got.bases.eq_ignore_ascii_case(seq.as_bytes()),
            "{name}: bases differ"
        );
        let qual: Vec<u8> = want["qual"]
            .as_array()
            .unwrap()
            .iter()
            .map(|q| q.as_u64().unwrap() as u8)
            .collect();
        assert_eq!(got.quality.as_ref(), Some(&qual), "{name}: qualities");
        let ploc: Vec<u32> = want["ploc2"]
            .as_array()
            .unwrap()
            .iter()
            .map(|q| q.as_u64().unwrap() as u32)
            .collect();
        assert_eq!(got.peaks, ploc, "{name}: peaks");
        assert_eq!(
            got.sample_name.as_deref(),
            want["smpl"].as_str(),
            "{name}: sample name"
        );
        // DATA 9-12 hold G, A, T, C (FWO_ = GATC in every corpus file).
        assert_eq!(want["fwo"], "GATC", "{name}");
        for (tag, slot) in [("data9", 2), ("data10", 0), ("data11", 3), ("data12", 1)] {
            let ch = &got.channels[slot];
            assert_eq!(
                ch.len() as u64,
                want[tag]["len"].as_u64().unwrap(),
                "{name} {tag}"
            );
            let sum: u64 = ch.iter().map(|&v| u64::from(v)).sum();
            assert_eq!(sum, want[tag]["sum"].as_u64().unwrap(), "{name} {tag} sum");
        }
        checked += 1;
    }
    eprintln!("abif: {checked} traces match Biopython");
    assert!(checked >= 5);
}

/// Reads matched to a construct (.dna packets 0x10/0x11, ZTR traces): every listed trace
/// reads, with one peak per base inside its signal, and at every Q30 call the called base's
/// own channel is the highest at its peak (which pins the decoding and the channel order).
#[test]
fn matched_reads_decode() {
    let mut paths = Vec::new();
    files(&corpus_dir(), "dna", 0, &mut paths);
    let (mut traces, mut sequences) = (0, 0);
    let mut q30 = (0usize, 0usize);
    for path in &paths {
        let f = DnaFile::parse(&std::fs::read(path).unwrap()).unwrap();
        for r in f.matched_reads().unwrap() {
            let read = r
                .read
                .unwrap_or_else(|e| panic!("{}: {}: {e}", path.display(), r.name));
            if let Some((a, b)) = r.trimmed {
                assert!(a < b && b <= read.bases.len(), "{} {:?}", r.name, r.trimmed);
            }
            sequences += 1;
            let Some(c) = read.trace else {
                continue;
            };
            assert_eq!(c.peaks.len(), c.bases.len(), "{}", r.name);
            let samples = c.samples();
            assert!(
                c.peaks.iter().all(|&p| (p as usize) < samples),
                "{}",
                r.name
            );
            let quality = c.quality.clone().unwrap();
            for ((b, &p), &q) in c.bases.iter().zip(&c.peaks).zip(&quality) {
                let Some(own) = b"ACGT".iter().position(|x| x == b) else {
                    continue;
                };
                if q < 30 {
                    continue;
                }
                q30.0 += 1;
                let at = |k: usize| c.channels[k][p as usize];
                if (0..4).all(|k| k == own || at(k) <= at(own)) {
                    q30.1 += 1;
                }
            }
            traces += 1;
        }
    }
    eprintln!(
        "matched reads: {sequences} reads, {traces} traces, {}/{} Q30 calls highest in their own channel",
        q30.1, q30.0
    );
    if traces == 0 && sequences == 0 {
        println!("no .dna file in the corpus carries matched reads: nothing to check");
        return;
    }
    assert!(traces >= 6);
    assert_eq!(q30.0, q30.1);
}
