//! GFF3 ↔ document.
//!
//! A GFF3 file opens as a document when its `##FASTA` section holds a sequence for one of its
//! seqids (the first such seqid is shown; other seqids' rows stay in the file). A file with
//! annotations only adds its features to an open document ([`annotations`]).
//!
//! Rows sharing an `ID` are one feature with several parts (a spliced CDS, for example). The
//! name comes from `Name` (unless it merely repeats `protein_id`, as NCBI writes CDS rows),
//! else `gene`, `product`, `locus_tag` or the type. Every other attribute becomes a qualifier,
//! so `ID` and `Parent` survive a save. Types are Sequence Ontology terms; the common ones
//! with an INSDC equivalent are shown with the INSDC key (`origin_of_replication` as
//! `rep_origin`, ...; or NCBI's own `gbkey` attribute when present) and written back as SO
//! terms. Strands: `+` forward, `-` reverse, `.` and `?` none. A region row spanning the whole
//! sequence describes it (`Is_circular=true`) and isn't shown.

use snappiesnap_doc::{
    DocError, Feature, FeatureId, Molecule, NewFeature, Qualifier, Range, Segment, Seq, Strand,
    Topology,
};
use snappiesnap_formats::gff3::{Gff3File, Line, Row};

use crate::genbank::join_across_origin;
use crate::{Exported, Format, Imported, IoError, Source};

#[derive(Debug, Clone, PartialEq)]
pub struct Gff3Source {
    original: Vec<u8>,
    file: Gff3File,
    seqid: String,
    imported: Molecule,
    /// Document feature id → the file lines of its rows.
    feature_map: Vec<(FeatureId, Vec<usize>)>,
}

pub(crate) fn looks_like_gff3(bytes: &[u8]) -> bool {
    bytes.starts_with(b"##gff-version 3")
}

/// Sequence Ontology terms with an INSDC feature key, as listed in the SO's mapping to the
/// INSDC feature table (only where the names differ).
const SO_TO_INSDC: &[(&str, &str)] = &[
    ("origin_of_replication", "rep_origin"),
    ("five_prime_UTR", "5'UTR"),
    ("three_prime_UTR", "3'UTR"),
    ("polyA_signal_sequence", "polyA_signal"),
    ("primer_binding_site", "primer_bind"),
    ("protein_binding_site", "protein_bind"),
    ("binding_site", "misc_binding"),
    ("sequence_feature", "misc_feature"),
    ("mobile_genetic_element", "mobile_element"),
    ("signal_peptide", "sig_peptide"),
    ("mature_protein_region", "mat_peptide"),
    ("recombination_feature", "misc_recomb"),
    ("sequence_difference", "misc_difference"),
    ("sequence_variant", "variation"),
];

fn so_to_insdc(kind: &str) -> String {
    SO_TO_INSDC
        .iter()
        .find(|(so, _)| *so == kind)
        .map_or(kind, |(_, insdc)| insdc)
        .to_string()
}

fn insdc_to_so(kind: &str) -> String {
    SO_TO_INSDC
        .iter()
        .find(|(_, insdc)| *insdc == kind)
        .map_or(kind, |(so, _)| so)
        .to_string()
}

/// Rows that describe the whole sequence rather than annotate it.
fn describes_sequence(row: &Row, len: usize) -> bool {
    matches!(
        row.kind.as_str(),
        "region" | "chromosome" | "contig" | "supercontig" | "databank_entry"
    ) && row.start == 1
        && row.end == len as u64
}

fn feature_name(row: &Row) -> String {
    let name = row
        .attribute("Name")
        .filter(|n| Some(*n) != row.attribute("protein_id"));
    ["label"]
        .iter()
        .filter_map(|t| row.attribute(t))
        .chain(name)
        .chain(
            ["gene", "product", "locus_tag"]
                .iter()
                .filter_map(|t| row.attribute(t)),
        )
        .find(|n| !n.trim().is_empty())
        .map_or_else(|| so_to_insdc(&row.kind), str::to_string)
}

/// Turn a seqid's rows into document features. Returns each feature with the indices of its
/// rows in `lines`, and warnings for rows that can't be shown.
fn rows_to_features(
    lines: &[Line],
    seqid: &str,
    len: usize,
    topology: Topology,
    warnings: &mut Vec<String>,
) -> Vec<(NewFeature, Vec<usize>)> {
    // Group rows by ID (rows without one are features of their own), in file order.
    let mut groups: Vec<(Option<String>, Vec<usize>)> = Vec::new();
    for (i, l) in lines.iter().enumerate() {
        let Line::Feature { row, .. } = l else {
            continue;
        };
        if row.seqid != seqid || describes_sequence(row, len) {
            continue;
        }
        let id = row.attribute("ID").map(str::to_string);
        match id
            .as_ref()
            .and_then(|id| groups.iter_mut().find(|(g, _)| g.as_ref() == Some(id)))
        {
            Some((_, rows)) => rows.push(i),
            None => groups.push((id, vec![i])),
        }
    }
    let row_at = |i: usize| match lines.get(i) {
        Some(Line::Feature { row, .. }) => Some(row),
        _ => None,
    };
    let mut out = Vec::new();
    for (_, indices) in groups {
        let rows: Vec<&Row> = indices.iter().filter_map(|&i| row_at(i)).collect();
        let Some(first) = rows.first() else { continue };
        match to_new_feature(&rows, len, topology) {
            Ok(nf) => out.push((nf, indices)),
            Err(reason) => warnings.push(format!(
                "feature {} at {}..{} not shown ({reason}); it is kept in the file",
                first.kind, first.start, first.end
            )),
        }
    }
    out
}

fn to_new_feature(rows: &[&Row], len: usize, topology: Topology) -> Result<NewFeature, String> {
    let first = rows.first().ok_or("no rows")?;
    let strand = match first.strand {
        '+' => Strand::Forward,
        '-' => Strand::Reverse,
        _ => Strand::None,
    };
    if rows.iter().any(|r| r.strand != first.strand) {
        return Err("its parts are on different strands".into());
    }
    let mut parts: Vec<(u64, u64)> = rows.iter().map(|r| (r.start, r.end)).collect();
    parts.sort_unstable();
    let mut ranges = Vec::new();
    for (start, end) in parts {
        let n = len as u64;
        let range = if end <= n {
            Range::new(start as usize - 1, end as usize)
        } else if topology == Topology::Circular && end - n < start && start <= n {
            // Across the origin: the end is given past the length.
            Range::new(start as usize - 1, (end - n) as usize)
        } else {
            return Err(format!("{start}..{end} is outside the {len} bp sequence"));
        };
        ranges.push(range);
    }
    // Rows are in ascending order here; join two meeting at the origin, then put the parts in
    // biological (5'→3') order: descending on the reverse strand.
    let mut ranges = if topology == Topology::Circular {
        join_across_origin(ranges, len)
    } else {
        ranges
    };
    if strand == Strand::Reverse {
        ranges.reverse();
    }
    let kind = first
        .attribute("gbkey")
        .filter(|k| !matches!(*k, "Src" | "Gene" | "Prot"))
        .map_or_else(|| so_to_insdc(&first.kind), str::to_string);
    Ok(NewFeature {
        color: None,
        name: feature_name(first),
        kind,
        strand,
        segments: ranges.into_iter().map(|range| Segment { range }).collect(),
        qualifiers: first
            .attributes
            .iter()
            .filter(|(t, _)| t != "Name")
            .map(|(t, v)| Qualifier {
                name: t.clone(),
                values: v.clone(),
            })
            .collect(),
    })
}

/// The seqid a document is shown from: the first with a sequence.
fn shown_seqid(file: &Gff3File) -> Option<String> {
    let ids = file.seqids();
    ids.iter()
        .find(|id| file.sequence(id).is_some())
        .cloned()
        .or_else(|| file.sequences.first().map(|s| s.id.clone()))
}

fn is_circular(file: &Gff3File, seqid: &str, len: usize) -> bool {
    file.lines.iter().any(|l| {
        matches!(l, Line::Feature { row, .. }
            if row.seqid == seqid
                && describes_sequence(row, len)
                && row.attribute("Is_circular").is_some_and(|v| v.eq_ignore_ascii_case("true")))
    })
}

pub(crate) fn import(name: &str, bytes: &[u8]) -> Result<Imported, IoError> {
    let file = Gff3File::parse(bytes).map_err(|e| IoError::Format(e.to_string()))?;
    let seqid = shown_seqid(&file).ok_or_else(|| {
        IoError::Format(
            "this GFF3 file has annotations but no sequence (no ##FASTA section): open the sequence, then import the annotations onto it".into(),
        )
    })?;
    let sequence = file
        .sequence(&seqid)
        .ok_or_else(|| IoError::Format(format!("no sequence for {seqid}")))?;
    let residues = sequence.residues.to_ascii_uppercase();
    let alphabet = crate::fasta::guess_alphabet(&residues);
    let seq = Seq::new(alphabet, residues).map_err(DocError::from)?;
    let len = seq.len();
    let topology = if is_circular(&file, &seqid, len) {
        Topology::Circular
    } else {
        Topology::Linear
    };
    let mut m = Molecule::new(name, seq, topology)?;
    let mut warnings = Vec::new();
    let mut feature_map = Vec::new();
    for (nf, rows) in rows_to_features(&file.lines, &seqid, len, topology, &mut warnings) {
        match m.add_feature(nf) {
            Ok(id) => feature_map.push((id, rows)),
            Err(e) => warnings.push(format!(
                "a feature isn't shown ({e}); it is kept in the file"
            )),
        }
    }
    let others = file.seqids().into_iter().filter(|s| *s != seqid).count()
        + file.sequences.iter().filter(|s| s.id != seqid).count();
    if others > 0 {
        warnings.push(format!(
            "the file describes more sequences than {seqid}; they are kept but not shown"
        ));
    }
    Ok(Imported {
        source: Source::Gff3(Box::new(Gff3Source {
            original: bytes.to_vec(),
            file,
            seqid,
            imported: m.clone(),
            feature_map,
        })),
        molecule: m,
        format: Format::Gff3,
        warnings,
    })
}

/// The features of a GFF3 file (annotations only, or with sequences) for a document of `len`
/// residues: those of `seqid`, or of the only seqid, or of the one whose region has that length.
/// Returns the seqid used, the features and warnings.
pub fn annotations(
    bytes: &[u8],
    seqid: Option<&str>,
    len: usize,
    topology: Topology,
) -> Result<(String, Vec<NewFeature>, Vec<String>), IoError> {
    let file = Gff3File::parse(bytes).map_err(|e| IoError::Format(e.to_string()))?;
    let ids = file.seqids();
    let chosen = match seqid {
        Some(s) if ids.iter().any(|i| i == s) => s.to_string(),
        Some(s) => {
            return Err(IoError::Format(format!(
                "the file has no rows for {s:?} (it has: {})",
                ids.join(", ")
            )));
        }
        None => match ids.as_slice() {
            [one] => one.clone(),
            [] => return Err(IoError::Format("the file has no feature rows".into())),
            many => many
                .iter()
                .find(|id| {
                    file.lines.iter().any(|l| {
                        matches!(l, Line::Feature { row, .. }
                            if &row.seqid == *id && describes_sequence(row, len))
                    })
                })
                .cloned()
                .ok_or_else(|| {
                    IoError::Format(format!(
                        "the file annotates several sequences ({}); choose one with `seqid`",
                        many.join(", ")
                    ))
                })?,
        },
    };
    let mut warnings = Vec::new();
    let features = rows_to_features(&file.lines, &chosen, len, topology, &mut warnings)
        .into_iter()
        .map(|(nf, _)| nf)
        .collect();
    Ok((chosen, features, warnings))
}

/// A seqid for a document: its name without whitespace.
fn seqid_for(name: &str) -> String {
    let s: String = name
        .chars()
        .map(|c| {
            if c.is_whitespace() || c.is_control() {
                '_'
            } else {
                c
            }
        })
        .collect();
    if s.is_empty() { "sequence".into() } else { s }
}

/// A document feature as rows: one per part, an end past the length for a part across the
/// origin, ascending order, and phases for CDS parts (in 5'→3' order).
fn feature_rows(
    f: &Feature,
    len: usize,
    seqid: &str,
    source: &str,
    kind: &str,
    id: &str,
) -> Vec<Row> {
    let strand = match f.strand {
        Strand::Forward => '+',
        Strand::Reverse => '-',
        _ => '.',
    };
    let mut attributes = vec![("ID".to_string(), vec![id.to_string()])];
    if !f.name.is_empty() {
        attributes.push(("Name".into(), vec![f.name.clone()]));
    }
    attributes.extend(
        f.qualifiers
            .iter()
            .filter(|q| q.name != "ID" && q.name != "Name" && !q.name.is_empty())
            .map(|q| (q.name.clone(), q.values.clone())),
    );
    // Phases along the coding direction.
    let mut done = 0usize;
    let mut rows: Vec<Row> = f
        .segments
        .iter()
        .map(|s| {
            let r = s.range;
            let (start, end) = if r.is_wrapping() {
                (r.start + 1, len + r.end)
            } else {
                (r.start + 1, r.end)
            };
            let phase = if kind == "CDS" {
                let p = (3 - done % 3) % 3;
                done += r.len_on(len);
                p.to_string()
            } else {
                String::new()
            };
            Row {
                seqid: seqid.to_string(),
                source: source.to_string(),
                kind: kind.to_string(),
                start: start as u64,
                end: end as u64,
                score: String::new(),
                strand,
                phase,
                attributes: attributes.clone(),
            }
        })
        .collect();
    rows.sort_by_key(|r| r.start);
    rows
}

pub(crate) fn export(m: &Molecule, source: Option<&Gff3Source>) -> Result<Exported, IoError> {
    let mut warnings = Vec::new();
    if !m.primers().is_empty() {
        warnings.push("primers aren't saved in GFF3 files".to_string());
    }
    if let Some(src) = source
        && *m == src.imported
    {
        return Ok(Exported {
            bytes: src.original.clone(),
            source: Source::Gff3(Box::new(src.clone())),
            warnings,
        });
    }
    let len = m.len();
    let (mut file, seqid, imported, feature_map) = match source {
        Some(src) => (
            src.file.clone(),
            src.seqid.clone(),
            Some(&src.imported),
            src.feature_map.clone(),
        ),
        None => {
            let seqid = seqid_for(m.name());
            let mut file = Gff3File::new();
            file.lines
                .push(Line::Other(format!("##sequence-region {seqid} 1 {len}\n")));
            file.lines.push(Line::Feature {
                row: region_row(&seqid, len, m.topology()),
                raw: None,
            });
            (file, seqid, None, Vec::new())
        }
    };
    let old_len = imported.map_or(len, Molecule::len);
    file.set_sequence(&seqid, m.seq().as_bytes().to_vec());

    // The sequence's own description: length and topology.
    for l in &mut file.lines {
        match l {
            Line::Other(text) if old_len != len => {
                let words: Vec<&str> = text.split_whitespace().collect();
                if words.first() == Some(&"##sequence-region")
                    && words.get(1) == Some(&seqid.as_str())
                {
                    *text = format!("##sequence-region {seqid} 1 {len}\n");
                }
            }
            Line::Feature { row, raw }
                if row.seqid == seqid && describes_sequence(row, old_len) =>
            {
                let circular = m.topology() == Topology::Circular;
                let was = row
                    .attribute("Is_circular")
                    .is_some_and(|v| v.eq_ignore_ascii_case("true"));
                if row.end != len as u64 || was != circular {
                    row.end = len as u64;
                    row.attributes.retain(|(t, _)| t != "Is_circular");
                    if circular {
                        row.attributes
                            .push(("Is_circular".into(), vec!["true".into()]));
                    }
                    *raw = None;
                }
            }
            _ => {}
        }
    }

    // Features: unchanged rows stay; changed features get new rows in place of their old ones;
    // removed ones go; new ones follow the last row of this sequence.
    let mut replace: Vec<(usize, Vec<Row>)> = Vec::new();
    let mut drop: Vec<usize> = Vec::new();
    for (id, rows) in &feature_map {
        let then = imported.and_then(|i| i.feature(*id));
        match (m.feature(*id), then) {
            (None, _) => drop.extend(rows),
            (Some(now), Some(then)) if now != then => {
                let first = rows.first().and_then(|&i| match file.lines.get(i) {
                    Some(Line::Feature { row, .. }) => Some(row.clone()),
                    _ => None,
                });
                let kind = match &first {
                    Some(r) if now.kind == then.kind => r.kind.clone(),
                    _ => insdc_to_so(&now.kind),
                };
                let source_col = first
                    .as_ref()
                    .map_or_else(String::new, |r| r.source.clone());
                let row_id = now
                    .qualifiers
                    .iter()
                    .find(|q| q.name == "ID")
                    .and_then(|q| q.values.first().cloned())
                    .or_else(|| {
                        first
                            .as_ref()
                            .and_then(|r| r.attribute("ID").map(str::to_string))
                    })
                    .unwrap_or_else(|| format!("feature-{}", id.0));
                let new_rows = feature_rows(now, len, &seqid, &source_col, &kind, &row_id);
                if let Some((&head, rest)) = rows.split_first() {
                    replace.push((head, new_rows));
                    drop.extend(rest);
                }
            }
            _ => {}
        }
    }
    let mapped: Vec<FeatureId> = feature_map.iter().map(|(id, _)| *id).collect();
    let added: Vec<(FeatureId, Vec<Row>)> = m
        .features()
        .iter()
        .filter(|f| !mapped.contains(&f.id))
        .map(|f| {
            let row_id = f
                .qualifiers
                .iter()
                .find(|q| q.name == "ID")
                .and_then(|q| q.values.first().cloned())
                .unwrap_or_else(|| format!("feature-{}", f.id.0));
            (
                f.id,
                feature_rows(
                    f,
                    len,
                    &seqid,
                    "SnappieSnap",
                    &insdc_to_so(&f.kind),
                    &row_id,
                ),
            )
        })
        .collect();
    let insert_after = file
        .lines
        .iter()
        .rposition(|l| matches!(l, Line::Feature { row, .. } if row.seqid == seqid))
        .unwrap_or(file.lines.len().saturating_sub(1));

    // Rebuild the lines, tracking where each kept line goes.
    let mut lines = Vec::new();
    let mut moved: Vec<Option<Vec<usize>>> = vec![None; file.lines.len()];
    let mut new_rows_at: Vec<(FeatureId, Vec<usize>)> = Vec::new();
    for (i, l) in file.lines.iter().enumerate() {
        if let Some((_, rows)) = replace.iter().find(|(head, _)| *head == i) {
            let start = lines.len();
            for r in rows {
                lines.push(Line::Feature {
                    row: r.clone(),
                    raw: None,
                });
            }
            if let Some(slot) = moved.get_mut(i) {
                *slot = Some((start..lines.len()).collect());
            }
        } else if !drop.contains(&i) {
            if let Some(slot) = moved.get_mut(i) {
                *slot = Some(vec![lines.len()]);
            }
            lines.push(l.clone());
        }
        if i == insert_after {
            for (id, rows) in &added {
                let start = lines.len();
                for r in rows {
                    lines.push(Line::Feature {
                        row: r.clone(),
                        raw: None,
                    });
                }
                new_rows_at.push((*id, (start..lines.len()).collect()));
            }
        }
    }
    if file.lines.is_empty() {
        for (id, rows) in &added {
            let start = lines.len();
            lines.extend(rows.iter().map(|r| Line::Feature {
                row: r.clone(),
                raw: None,
            }));
            new_rows_at.push((*id, (start..lines.len()).collect()));
        }
    }
    file.lines = lines;
    let mut new_map: Vec<(FeatureId, Vec<usize>)> = feature_map
        .iter()
        .filter(|(id, _)| m.feature(*id).is_some())
        .filter_map(|(id, rows)| {
            let first = rows.first()?;
            let at = moved.get(*first)?.clone()?;
            // A replaced feature owns all its new rows; an unchanged one keeps its rows.
            let all: Vec<usize> = if replace.iter().any(|(h, _)| h == first) {
                at
            } else {
                rows.iter()
                    .filter_map(|r| moved.get(*r).cloned().flatten())
                    .flatten()
                    .collect()
            };
            Some((*id, all))
        })
        .collect();
    new_map.extend(new_rows_at);

    let bytes = file.to_text().into_bytes();
    let reread = Gff3File::parse(&bytes).map_err(|e| IoError::Format(e.to_string()))?;
    Ok(Exported {
        source: Source::Gff3(Box::new(Gff3Source {
            original: bytes.clone(),
            file: reread,
            seqid,
            imported: m.clone(),
            feature_map: new_map,
        })),
        bytes,
        warnings,
    })
}

fn region_row(seqid: &str, len: usize, topology: Topology) -> Row {
    let mut attributes = vec![("ID".to_string(), vec![seqid.to_string()])];
    if topology == Topology::Circular {
        attributes.push(("Is_circular".into(), vec!["true".into()]));
    }
    Row {
        seqid: seqid.to_string(),
        source: "SnappieSnap".into(),
        kind: "region".into(),
        start: 1,
        end: len as u64,
        score: String::new(),
        strand: '+',
        phase: String::new(),
        attributes,
    }
}

#[cfg(test)]
mod tests {
    #![allow(clippy::unwrap_used, clippy::indexing_slicing)]
    use super::*;

    const FILE: &str = "##gff-version 3
##sequence-region c1 1 40
c1\tlab\tregion\t1\t40\t.\t+\t.\tID=c1;Is_circular=true
c1\tlab\tCDS\t31\t46\t.\t+\t0\tID=wrap;Name=across
c1\tlab\tCDS\t3\t8\t.\t-\t0\tID=split;Name=two parts;gene=g
c1\tlab\tCDS\t12\t20\t.\t-\t0\tID=split;Name=two parts;gene=g
c1\tlab\torigin_of_replication\t21\t28\t.\t.\t.\tID=ori
##FASTA
>c1
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
";

    #[test]
    fn opens_features_across_the_origin_and_in_parts() {
        let i = import("c1", FILE.as_bytes()).unwrap();
        let m = &i.molecule;
        assert_eq!(m.topology(), Topology::Circular);
        assert_eq!(m.features().len(), 3);
        let wrap = m.features().iter().find(|f| f.name == "across").unwrap();
        assert_eq!(wrap.segments[0].range, Range::new(30, 6));
        let split = m.features().iter().find(|f| f.name == "two parts").unwrap();
        assert_eq!(split.strand, Strand::Reverse);
        // 5'→3' on the reverse strand: the later part first.
        assert_eq!(
            split.segments.iter().map(|s| s.range).collect::<Vec<_>>(),
            vec![Range::new(11, 20), Range::new(2, 8)]
        );
        let ori = m
            .features()
            .iter()
            .find(|f| f.name == "rep_origin")
            .unwrap();
        assert_eq!(
            (ori.kind.as_str(), ori.strand),
            ("rep_origin", Strand::None)
        );
    }

    #[test]
    fn saves_unchanged_bytes_and_minimal_edits() {
        let i = import("c1", FILE.as_bytes()).unwrap();
        let src = match &i.source {
            Source::Gff3(s) => s.as_ref().clone(),
            _ => panic!(),
        };
        assert_eq!(
            export(&i.molecule, Some(&src)).unwrap().bytes,
            FILE.as_bytes()
        );
        // Rename the split CDS, remove the origin, add a feature.
        let mut m = i.molecule.clone();
        let mut split = m
            .features()
            .iter()
            .find(|f| f.name == "two parts")
            .unwrap()
            .clone();
        split.name = "renamed".into();
        m.update_feature(split).unwrap();
        let ori = m
            .features()
            .iter()
            .find(|f| f.kind == "rep_origin")
            .unwrap()
            .id;
        m.remove_feature(ori).unwrap();
        m.add_feature(NewFeature {
            color: None,
            name: "new".into(),
            kind: "promoter".into(),
            strand: Strand::Forward,
            segments: vec![Segment {
                range: Range::new(0, 2),
            }],
            qualifiers: vec![],
        })
        .unwrap();
        let out = export(&m, Some(&src)).unwrap();
        let text = String::from_utf8(out.bytes.clone()).unwrap();
        assert!(text.contains("c1\tlab\tCDS\t31\t46\t.\t+\t0\tID=wrap;Name=across\n"));
        assert!(text.contains("c1\tlab\tCDS\t3\t8\t.\t-\t0\tID=split;Name=renamed;gene=g\n"));
        assert!(text.contains("c1\tlab\tCDS\t12\t20\t.\t-\t"));
        assert!(!text.contains("origin_of_replication"));
        assert!(text.contains("c1\tSnappieSnap\tpromoter\t1\t2\t.\t+\t.\tID=feature-"));
        let again = import("c1", text.as_bytes()).unwrap();
        assert_eq!(again.molecule.features().len(), 3);
        // A second save from the written source changes nothing more.
        let src2 = match &out.source {
            Source::Gff3(s) => s.as_ref().clone(),
            _ => panic!(),
        };
        assert_eq!(export(&m, Some(&src2)).unwrap().bytes, out.bytes);
    }

    #[test]
    fn two_rows_meeting_at_the_origin_are_one_part() {
        let file = "##gff-version 3\nc\t.\tregion\t1\t20\t.\t+\t.\tID=c;Is_circular=true\nc\t.\tCDS\t18\t20\t.\t-\t0\tID=x;Name=x\nc\t.\tCDS\t1\t3\t.\t-\t0\tID=x;Name=x\n##FASTA\n>c\nACGTACGTACGTACGTACGT\n";
        let i = import("c", file.as_bytes()).unwrap();
        let f = &i.molecule.features()[0];
        assert_eq!(f.strand, Strand::Reverse);
        assert_eq!(
            f.segments.iter().map(|s| s.range).collect::<Vec<_>>(),
            vec![Range::new(17, 3)]
        );
    }

    #[test]
    fn annotations_only() {
        let rows = "##gff-version 3\nc1\t.\tgene\t2\t5\t.\t+\t.\tID=g1;Name=a\nc2\t.\tgene\t1\t3\t.\t-\t.\tName=b\n";
        assert!(import("x", rows.as_bytes()).is_err());
        assert!(annotations(rows.as_bytes(), None, 40, Topology::Linear).is_err());
        let (id, features, _) =
            annotations(rows.as_bytes(), Some("c2"), 40, Topology::Linear).unwrap();
        assert_eq!((id.as_str(), features.len()), ("c2", 1));
        assert_eq!(features[0].strand, Strand::Reverse);
        // Outside a short sequence: a warning, not a feature.
        let (_, features, warnings) =
            annotations(rows.as_bytes(), Some("c1"), 3, Topology::Linear).unwrap();
        assert!(features.is_empty() && warnings.len() == 1);
    }

    #[test]
    fn new_files_describe_the_sequence() {
        let i = import("c1", FILE.as_bytes()).unwrap();
        let out = export(&i.molecule, None).unwrap();
        let text = String::from_utf8(out.bytes).unwrap();
        assert!(text.starts_with("##gff-version 3\n##sequence-region c1 1 40\nc1\tSnappieSnap\tregion\t1\t40\t.\t+\t.\tID=c1;Is_circular=true\n"));
        let again = import("c1", text.as_bytes()).unwrap();
        assert_eq!(again.molecule.seq(), i.molecule.seq());
        assert_eq!(again.molecule.features().len(), 3);
        let wrap = again
            .molecule
            .features()
            .iter()
            .find(|f| f.name == "across")
            .unwrap();
        assert_eq!(wrap.segments[0].range, Range::new(30, 6));
    }
}
