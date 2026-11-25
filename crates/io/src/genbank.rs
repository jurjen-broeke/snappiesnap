//! GenBank ↔ document.
//!
//! The document shows the first record. Feature names come from `/label` (as the reference
//! application writes them), else `/gene`, `/product` or `/locus_tag`. Strands: `complement(...)`
//! is reverse; otherwise an explicit `direction:` note (written by the reference application) is
//! used, else coding and transcript features (CDS, gene, RNAs) are forward and the rest
//! (promoters, origins, ...) have no direction, as the reference application imports them. A segment pair split at the origin of a circular record
//! (`join(3700..3781,1..10)`) becomes one wrapping range, and is written back the same way.

use snappiesnap_doc::{
    Alphabet, DocError, Feature, FeatureId, Molecule, NewFeature, Notes, Qualifier, Range, Segment,
    Seq, Strand, Topology,
};
use snappiesnap_formats::genbank::{
    FeatureEntry, GenBankFile, LocSegment, Location, Locus, Qualifier as GbQualifier, Record,
    format_location,
};

use crate::{ExportOptions, Exported, Format, Imported, IoError, Source};

#[derive(Debug, Clone, PartialEq)]
pub struct GenBankSource {
    original: Vec<u8>,
    file: GenBankFile,
    imported: Molecule,
    /// Document feature id → index in the first record's features.
    feature_map: Vec<(FeatureId, usize)>,
}

pub(crate) fn looks_like_genbank(bytes: &[u8]) -> bool {
    bytes.starts_with(b"LOCUS")
}

const FORWARD_KEYS: &[&str] = &[
    "CDS",
    "gene",
    "mRNA",
    "tRNA",
    "rRNA",
    "ncRNA",
    "misc_RNA",
    "precursor_RNA",
    "prim_transcript",
    "tmRNA",
];

pub(crate) fn import(name: &str, bytes: &[u8]) -> Result<Imported, IoError> {
    let file = GenBankFile::parse(bytes).map_err(|e| IoError::Format(e.to_string()))?;
    let record = file
        .records
        .first()
        .ok_or_else(|| IoError::Format("the file has no records".into()))?;
    let alphabet = if record.locus.protein {
        Alphabet::Protein
    } else if record.locus.molecule.to_ascii_uppercase().contains("RNA") {
        Alphabet::Rna
    } else {
        Alphabet::Dna
    };
    let residues = record.residues.to_ascii_uppercase();
    let seq = Seq::new(alphabet, residues).map_err(DocError::from)?;
    let topology = if record.locus.circular {
        Topology::Circular
    } else {
        Topology::Linear
    };
    let mut m = Molecule::new(name, seq, topology)?;
    let mut warnings = Vec::new();

    let description = definition(&record.header);
    if !description.is_empty() {
        let notes = Notes {
            description,
            ..Notes::default()
        };
        if let Err(e) = m.set_notes(notes) {
            warnings.push(format!("description not shown ({e})"));
        }
    }

    let feature_map = add_features(&record.features, &mut m, topology, &mut warnings);
    if file.records.len() > 1 {
        warnings.push(format!(
            "{} more record(s) in the file are kept but not shown",
            file.records.len() - 1
        ));
    }
    Ok(Imported {
        source: Source::GenBank(Box::new(GenBankSource {
            original: bytes.to_vec(),
            file,
            imported: m.clone(),
            feature_map,
        })),
        molecule: m,
        format: Format::GenBank,
        warnings,
    })
}

/// Show a record's feature table in the document; returns document id → feature index. The
/// record-level source (first, spanning the record) is kept in the file, not shown.
pub(crate) fn add_features(
    features: &[snappiesnap_formats::genbank::Feature],
    m: &mut Molecule,
    topology: Topology,
    warnings: &mut Vec<String>,
) -> Vec<(FeatureId, usize)> {
    let mut feature_map = Vec::new();
    for (index, f) in features.iter().enumerate() {
        // The record-level source (first, spanning the record) describes the whole record; it
        // is kept in the file, not shown. Other `source` features are real annotations.
        if index == 0 && f.entry.key == "source" && spans_record(&f.entry.location, m.len()) {
            continue;
        }
        let added = to_new_feature(&f.entry, m.len(), topology)
            .and_then(|nf| m.add_feature(nf).map_err(|e| e.to_string()));
        match added {
            Ok(id) => feature_map.push((id, index)),
            Err(reason) => warnings.push(format!(
                "feature {} at {} not shown ({reason}); it is kept in the file",
                f.entry.key, f.entry.location
            )),
        }
    }
    feature_map
}

/// Whether a location is exactly `1..len` (the whole record).
fn spans_record(location: &str, len: usize) -> bool {
    Location::parse(location)
        .ok()
        .and_then(|l| l.segments())
        .is_some_and(|s| matches!(s.as_slice(), [one] if one.start == 1 && one.end == len as u64 && !one.reverse))
}

/// The DEFINITION text, continuation lines joined.
fn definition(header: &str) -> String {
    let mut out = String::new();
    let mut inside = false;
    for line in header.lines() {
        if let Some(rest) = line.strip_prefix("DEFINITION") {
            inside = true;
            out.push_str(rest.trim());
        } else if inside && line.starts_with("            ") {
            out.push(' ');
            out.push_str(line.trim());
        } else if inside {
            break;
        }
    }
    let out = out.trim().to_string();
    if out == "." { String::new() } else { out }
}

fn feature_name(e: &FeatureEntry) -> String {
    ["label", "gene", "product", "locus_tag"]
        .iter()
        .find_map(|q| e.qualifier(q))
        .unwrap_or(&e.key)
        .to_string()
}

/// An explicit direction written by the reference application: a `/direction=RIGHT`
/// qualifier, or `direction: RIGHT` in a `/note`.
fn noted_direction(e: &FeatureEntry) -> Option<Strand> {
    let qualifier = e
        .qualifiers
        .iter()
        .filter(|q| q.name == "direction")
        .filter_map(|q| q.value.as_deref())
        .find_map(|v| word_direction(v.trim()));
    qualifier.or_else(|| {
        e.qualifiers
            .iter()
            .filter(|q| q.name == "note")
            .filter_map(|q| q.value.as_deref())
            .find_map(|v| {
                let d = v.split("direction:").nth(1)?.trim_start();
                let word: String = d.chars().take_while(|c| c.is_ascii_alphabetic()).collect();
                word_direction(&word)
            })
    })
}

fn word_direction(word: &str) -> Option<Strand> {
    match word {
        "RIGHT" => Some(Strand::Forward),
        "LEFT" => Some(Strand::Reverse),
        "BOTH" => Some(Strand::Both),
        "NONE" => Some(Strand::None),
        _ => None,
    }
}

fn to_new_feature(e: &FeatureEntry, len: usize, topology: Topology) -> Result<NewFeature, String> {
    let loc = Location::parse(&e.location)?;
    let segs = loc
        .segments()
        .ok_or("location points elsewhere or between bases")?;
    let reverse = segs.iter().all(|s| s.reverse);
    let strand = if reverse && !segs.is_empty() {
        Strand::Reverse
    } else {
        noted_direction(e).unwrap_or_else(|| default_strand(&e.key))
    };
    let mut ranges = Vec::new();
    for s in &segs {
        let (a, b) = (s.start.min(s.end), s.start.max(s.end));
        let r = crate::dna::segment_to_range(a, b, len)?;
        ranges.push(r);
    }
    let ranges = if topology == Topology::Circular {
        merge_origin(ranges, reverse, len)
    } else {
        ranges
    };
    Ok(NewFeature {
        color: None,
        name: feature_name(e),
        kind: e.key.clone(),
        strand,
        segments: ranges.into_iter().map(|range| Segment { range }).collect(),
        qualifiers: e
            .qualifiers
            .iter()
            .filter(|q| q.name != "label" && !is_metadata(q))
            .map(|q| Qualifier {
                name: q.name.clone(),
                values: q.value.iter().cloned().collect(),
            })
            .collect(),
    })
}

/// Display metadata written by the reference application: a `/direction=RIGHT` qualifier or a
/// `/note` holding only `color: #31849b; direction: RIGHT`. Kept in the file, not shown as a
/// qualifier.
fn is_metadata(q: &GbQualifier) -> bool {
    if q.name == "direction" {
        return q
            .value
            .as_deref()
            .is_some_and(|v| word_direction(v.trim()).is_some());
    }
    q.name == "note"
        && q.value.as_deref().is_some_and(|v| {
            !v.is_empty()
                && v.split(';').all(|part| {
                    let part = part.trim();
                    part.strip_prefix("color: #").is_some_and(|hex| {
                        hex.len() == 6 && hex.bytes().all(|b| b.is_ascii_hexdigit())
                    }) || matches!(
                        part,
                        "direction: RIGHT"
                            | "direction: LEFT"
                            | "direction: BOTH"
                            | "direction: NONE"
                    )
                })
        })
}

/// The strand a plain (not complemented) location gets on import without a direction note.
fn default_strand(key: &str) -> Strand {
    if FORWARD_KEYS.contains(&key) {
        Strand::Forward
    } else {
        Strand::None
    }
}

/// A direction note for a new feature whose strand the location and type alone wouldn't give
/// back on import.
fn direction_note(f: &Feature) -> Option<GbQualifier> {
    let word = match f.strand {
        Strand::Reverse => return None,
        s if s == default_strand(&f.kind) => return None,
        Strand::Forward => "RIGHT",
        Strand::Both => "BOTH",
        Strand::None => "NONE",
    };
    Some(GbQualifier {
        name: "note".into(),
        value: Some(format!("direction: {word}")),
        quoted: true,
    })
}

/// Merge consecutive pieces that meet at the origin into one wrapping range.
/// Parts in ascending order on a circular molecule: a part ending at the origin and one
/// starting there are one part across it (formats that can't write a wrapping range split it
/// in two, such as SBOL ranges and some GFF3 writers).
pub(crate) fn join_across_origin(mut ascending: Vec<Range>, len: usize) -> Vec<Range> {
    let joins = ascending.len() > 1
        && ascending
            .first()
            .is_some_and(|f| f.start == 0 && !f.is_wrapping())
        && ascending
            .last()
            .is_some_and(|l| l.end == len && !l.is_wrapping());
    if joins && let (Some(last), Some(first)) = (ascending.pop(), ascending.first_mut()) {
        // The joined part stays first: it starts the feature on the top strand.
        *first = Range::new(last.start, first.end);
    }
    ascending
}

pub(crate) fn merge_origin(ranges: Vec<Range>, reverse: bool, len: usize) -> Vec<Range> {
    let mut out: Vec<Range> = Vec::new();
    for r in ranges {
        if let Some(prev) = out.last_mut() {
            let forward_join = !reverse && prev.end == len && r.start == 0 && !prev.is_wrapping();
            let reverse_join = reverse && prev.start == 0 && r.end == len && !r.is_wrapping();
            if forward_join {
                *prev = Range::new(prev.start, r.end);
                continue;
            }
            if reverse_join {
                *prev = Range::new(r.start, prev.end);
                continue;
            }
        }
        out.push(r);
    }
    out
}

/// A document feature's location as GenBank segments (1-based, biological order).
fn loc_segments(f: &Feature, len: usize) -> Vec<LocSegment> {
    let reverse = f.strand == Strand::Reverse;
    let mut out = Vec::new();
    for s in &f.segments {
        let (first, second) = s.range.segments(len);
        let mut pieces: Vec<std::ops::Range<usize>> =
            std::iter::once(first).chain(second).collect();
        if reverse {
            pieces.reverse();
        }
        for p in pieces {
            out.push(LocSegment {
                start: p.start as u64 + 1,
                end: p.end as u64,
                reverse,
            });
        }
    }
    out
}

fn doc_qualifiers(f: &Feature) -> Vec<GbQualifier> {
    f.qualifiers
        .iter()
        .flat_map(|q| {
            let values: Vec<Option<String>> = if q.values.is_empty() {
                vec![None]
            } else {
                q.values.iter().cloned().map(Some).collect()
            };
            values.into_iter().map(|value| {
                let numeric = value
                    .as_deref()
                    .is_some_and(|v| !v.is_empty() && v.bytes().all(|b| b.is_ascii_digit()));
                GbQualifier {
                    name: q.name.clone(),
                    quoted: value.is_some() && !numeric,
                    value,
                }
            })
        })
        .collect()
}

fn label(name: &str) -> GbQualifier {
    GbQualifier {
        name: "label".into(),
        value: Some(name.to_string()),
        quoted: false,
    }
}

/// A LOCUS name: no spaces, at most 16 characters (the column width).
fn locus_name(name: &str) -> String {
    let n: String = name
        .chars()
        .map(|c| if c.is_whitespace() { '_' } else { c })
        .take(16)
        .collect();
    if n.is_empty() { "Untitled".into() } else { n }
}

pub(crate) fn export(
    m: &Molecule,
    source: Option<&GenBankSource>,
    options: &ExportOptions,
) -> Result<Exported, IoError> {
    let mut warnings = Vec::new();
    if !m.primers().is_empty() {
        warnings.push("primers aren't saved in GenBank files yet".to_string());
    }
    if let Some(src) = source
        && *m == src.imported
    {
        return Ok(Exported {
            bytes: src.original.clone(),
            source: Source::GenBank(Box::new(src.clone())),
            warnings,
        });
    }
    let len = m.len();
    let (mut file, imported, feature_map) = match source {
        Some(src) => (
            src.file.clone(),
            Some(&src.imported),
            src.feature_map.clone(),
        ),
        None => (new_file(m, options), None, Vec::new()),
    };
    let record = file
        .records
        .first_mut()
        .ok_or_else(|| IoError::Format("the file has no records".into()))?;

    // Sequence, length, topology, name.
    record.set_residues(m.seq().as_bytes().to_ascii_lowercase());
    let mut locus = record.locus.clone();
    locus.length = len as u64;
    locus.circular = m.topology() == Topology::Circular;
    if imported.is_some_and(|i| i.name() != m.name()) {
        locus.name = locus_name(m.name());
    }
    record.set_locus(locus);

    let feature_map = sync_features(&mut record.features, m, imported, &feature_map);

    let bytes = file.to_text().into_bytes();
    // The next save starts from what was written.
    let reread = GenBankFile::parse(&bytes).map_err(|e| IoError::Format(e.to_string()))?;
    Ok(Exported {
        source: Source::GenBank(Box::new(GenBankSource {
            original: bytes.clone(),
            file: reread,
            imported: m.clone(),
            feature_map,
        })),
        bytes,
        warnings,
    })
}

/// Bring a feature table up to date with the document: edited features are re-encoded,
/// removed ones dropped and new ones appended; untouched entries keep their raw text. Returns
/// the new document id → feature index map.
pub(crate) fn sync_features(
    features: &mut Vec<snappiesnap_formats::genbank::Feature>,
    m: &Molecule,
    imported: Option<&Molecule>,
    feature_map: &[(FeatureId, usize)],
) -> Vec<(FeatureId, usize)> {
    // Edited and removed features.
    let len = m.len();
    let mut feature_map = feature_map.to_vec();
    let mut keep = vec![true; features.len()];
    for (id, index) in &feature_map {
        let Some(slot) = features.get_mut(*index) else {
            continue;
        };
        let then = imported.and_then(|i| i.feature(*id));
        match (m.feature(*id), then) {
            (None, _) => {
                if let Some(k) = keep.get_mut(*index) {
                    *k = false;
                }
            }
            (Some(now), Some(then)) if now != then => {
                let e = &mut slot.entry;
                e.key = now.kind.clone();
                if now.segments != then.segments || now.strand != then.strand {
                    e.location = format_location(&loc_segments(now, len));
                }
                if now.name != then.name {
                    match e.qualifiers.iter_mut().find(|q| q.name == "label") {
                        Some(q) => *q = label(&now.name),
                        None => e.qualifiers.insert(0, label(&now.name)),
                    }
                }
                if now.qualifiers != then.qualifiers {
                    let label_q = e.qualifiers.iter().find(|q| q.name == "label").cloned();
                    let metadata: Vec<GbQualifier> = e
                        .qualifiers
                        .iter()
                        .filter(|q| is_metadata(q))
                        .cloned()
                        .collect();
                    e.qualifiers = label_q
                        .into_iter()
                        .chain(doc_qualifiers(now))
                        .chain(metadata)
                        .collect();
                }
                if now.strand != then.strand
                    && now.strand != Strand::Reverse
                    && now.strand != default_strand(&now.kind)
                    && noted_direction(e) != Some(now.strand)
                {
                    e.qualifiers.retain(|q| {
                        !(is_metadata(q)
                            && q.value.as_deref().is_some_and(|v| v.contains("direction:")))
                    });
                    e.qualifiers.extend(direction_note(now));
                }
                slot.raw = None;
            }
            _ => {}
        }
    }
    let mut index = 0usize;
    features.retain(|_| {
        let k = keep.get(index).copied().unwrap_or(true);
        index += 1;
        k
    });
    // Re-index the map after removals.
    let mut new_map = Vec::new();
    for (id, old) in &feature_map {
        if keep.get(*old).copied().unwrap_or(false) {
            let shift = keep.iter().take(*old).filter(|k| !**k).count();
            new_map.push((*id, old - shift));
        }
    }
    feature_map = new_map;

    // New features.
    let mapped: Vec<FeatureId> = feature_map.iter().map(|(id, _)| *id).collect();
    for f in m.features().iter().filter(|f| !mapped.contains(&f.id)) {
        let entry = FeatureEntry {
            key: f.kind.clone(),
            location: format_location(&loc_segments(f, len)),
            qualifiers: std::iter::once(label(&f.name))
                .chain(doc_qualifiers(f))
                .chain(direction_note(f))
                .collect(),
        };
        features.push(snappiesnap_formats::genbank::Feature { entry, raw: None });
        feature_map.push((f.id, features.len() - 1));
    }

    feature_map
}

/// A new file for a document that wasn't read from GenBank, laid out like the reference
/// application's export.
fn new_file(m: &Molecule, options: &ExportOptions) -> GenBankFile {
    let description = strip_tags(&m.notes().description);
    let definition = if description.trim().is_empty() {
        ".".to_string()
    } else {
        description
    };
    let header = format!(
        "DEFINITION  {definition}\nACCESSION   .\nVERSION     .\nKEYWORDS    .\nSOURCE      synthetic DNA construct\n  ORGANISM  synthetic DNA construct\n"
    );
    let locus = Locus {
        name: locus_name(m.name()),
        length: m.len() as u64,
        molecule: match m.seq().alphabet() {
            Alphabet::Rna => "RNA".into(),
            Alphabet::Protein => String::new(),
            Alphabet::Dna => "DNA".into(),
        },
        circular: m.topology() == Topology::Circular,
        division: "SYN".into(),
        date: options.date.clone().unwrap_or_default(),
        protein: m.seq().alphabet() == Alphabet::Protein,
    };
    let mut record = Record::new(locus, header);
    record.features.push(snappiesnap_formats::genbank::Feature {
        entry: FeatureEntry {
            key: "source".into(),
            location: format!("1..{}", m.len()),
            qualifiers: vec![
                GbQualifier {
                    name: "mol_type".into(),
                    value: Some("other DNA".into()),
                    quoted: true,
                },
                GbQualifier {
                    name: "organism".into(),
                    value: Some("synthetic DNA construct".into()),
                    quoted: true,
                },
            ],
        },
        raw: None,
    });
    GenBankFile {
        records: vec![record],
        preamble: String::new(),
    }
}

/// Text without `<...>` markup (descriptions may be HTML).
pub(crate) fn strip_tags(s: &str) -> String {
    let mut out = String::new();
    let mut inside = false;
    for c in s.chars() {
        match c {
            '<' => inside = true,
            '>' => inside = false,
            c if !inside => out.push(c),
            _ => {}
        }
    }
    out.replace("&nbsp;", " ").trim().to_string()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn origin_pairs_merge_both_ways() {
        let r = Range::new;
        assert_eq!(
            merge_origin(vec![r(90, 100), r(0, 10)], false, 100),
            vec![r(90, 10)]
        );
        assert_eq!(
            merge_origin(vec![r(0, 10), r(90, 100)], true, 100),
            vec![r(90, 10)]
        );
        assert_eq!(
            merge_origin(vec![r(10, 20), r(30, 40)], false, 100),
            vec![r(10, 20), r(30, 40)]
        );
    }

    #[test]
    fn helpers() {
        assert_eq!(
            strip_tags("<html><body>A &amp; B&nbsp;</body></html>"),
            "A &amp; B"
        );
        assert_eq!(locus_name("my long plasmid name here"), "my_long_plasmid_");
        assert_eq!(
            definition("DEFINITION  One\n            two.\nACCESSION   .\n"),
            "One two."
        );
        assert_eq!(definition("DEFINITION  .\n"), "");
    }

    const PLAIN: &str =
        "LOCUS       pT                       100 bp    DNA     circular SYN 01-JAN-2020
FEATURES             Location/Qualifiers
     source          1..100
                     /mol_type=\"other DNA\"
     promoter        1..10
                     /label=prom
     terminator      11..20
                     /label=term
     CDS             complement(21..40)
                     /label=rev
     protein_bind    41..50
                     /label=site
                     /note=\"color: #31849b; direction: RIGHT\"
     misc_feature    join(91..100,1..5)
                     /label=wrap
     CDS             complement(join(91..100,1..5))
                     /label=wraprev
ORIGIN
        1 acgtacgtac gtacgtacgt acgtacgtac gtacgtacgt acgtacgtac gtacgtacgt
       61 acgtacgtac gtacgtacgt acgtacgtac gtacgtacgt
//
";

    #[test]
    fn import_rules_for_plain_locations() {
        let i = import("pT", PLAIN.as_bytes()).unwrap();
        assert!(i.warnings.is_empty(), "{:?}", i.warnings);
        let f = |name: &str| {
            i.molecule
                .features()
                .iter()
                .find(|f| f.name == name)
                .cloned()
                .unwrap()
        };
        // As the reference application imports plain promoters: no direction.
        assert_eq!(f("prom").strand, Strand::None);
        assert_eq!(f("term").strand, Strand::None);
        assert_eq!(f("rev").strand, Strand::Reverse);
        assert_eq!(f("site").strand, Strand::Forward, "direction note wins");
        let q = PLAIN.replace(
            "                     /label=term",
            "                     /direction=LEFT\n                     /label=term",
        );
        let j = import("pT", q.as_bytes()).unwrap();
        let term = j
            .molecule
            .features()
            .iter()
            .find(|f| f.name == "term")
            .cloned()
            .unwrap();
        assert_eq!(term.strand, Strand::Reverse, "/direction qualifier");
        assert!(term.qualifiers.is_empty());
        assert!(
            f("site").qualifiers.is_empty(),
            "metadata notes aren't shown"
        );
        assert_eq!(
            f("wrap").segments,
            vec![Segment {
                range: Range::new(90, 5)
            }]
        );
        assert_eq!(
            f("wraprev").segments,
            vec![Segment {
                range: Range::new(90, 5)
            }]
        );
        assert_eq!(f("wraprev").strand, Strand::Reverse);
        assert_eq!(i.molecule.features().len(), 6, "source isn't shown");
        // A source feature that isn't the record-level one is a real annotation.
        let extra = PLAIN.replace("     promoter        1..10", "     source          3..9\n                     /label=insert\n     promoter        1..10");
        let j = import("pT", extra.as_bytes()).unwrap();
        assert!(
            j.molecule
                .features()
                .iter()
                .any(|f| f.name == "insert" && f.kind == "source")
        );
    }

    #[test]
    fn new_features_get_direction_notes_only_when_needed() {
        let i = import("pT", PLAIN.as_bytes()).unwrap();
        let out = export(&i.molecule, None, &ExportOptions::default()).unwrap();
        let text = String::from_utf8(out.bytes).unwrap();
        // Forward protein_bind needs a note; the none-strand promoter and terminator don't.
        assert_eq!(text.matches("direction: RIGHT").count(), 1, "{text}");
        assert!(!text.contains("direction: NONE"));
        assert!(
            text.contains("     misc_feature    join(91..100,1..5)\n"),
            "{text}"
        );
        assert!(
            text.contains("     CDS             complement(join(91..100,1..5))\n"),
            "{text}"
        );
        let back = import("x", text.as_bytes()).unwrap().molecule;
        for (a, b) in back.features().iter().zip(i.molecule.features()) {
            assert_eq!(
                (&a.name, &a.segments, a.strand),
                (&b.name, &b.segments, b.strand)
            );
        }
    }

    #[test]
    fn edits_keep_untouched_features_and_metadata() {
        let i = import("pT", PLAIN.as_bytes()).unwrap();
        let Source::GenBank(src) = &i.source else {
            panic!()
        };
        let mut m = i.molecule.clone();
        let mut site = m
            .features()
            .iter()
            .find(|f| f.name == "site")
            .cloned()
            .unwrap();
        site.qualifiers = vec![Qualifier {
            name: "note".into(),
            values: vec!["edited".into()],
        }];
        m.update_feature(site).unwrap();
        let text = String::from_utf8(
            export(&m, Some(src), &ExportOptions::default())
                .unwrap()
                .bytes,
        )
        .unwrap();
        assert!(text.contains("     protein_bind    41..50\n                     /label=site\n                     /note=\"edited\"\n                     /note=\"color: #31849b; direction: RIGHT\"\n"), "{text}");
        assert!(text.contains("     promoter        1..10\n                     /label=prom\n"));
    }
}
