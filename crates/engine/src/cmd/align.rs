//! Alignment of a sequence (a read, an assembly, another construct) to the open document.

use schemars::JsonSchema;
use serde::Deserialize;
use serde_json::{Value, json};
use snappiesnap_bio::align::{
    Mode, ProteinScoring, Scoring, align_protein_report, align_to_molecule,
};
use snappiesnap_bio::trace::{DEFAULT_CUTOFF, quality_trim};
use snappiesnap_doc::{Alphabet, FeatureId};

use crate::{CmdError, CommandSpec, ErrorCode, Session};

pub(super) fn specs() -> Vec<CommandSpec> {
    vec![
        CommandSpec::query(
            "align_to_reference",
            "Align to This Sequence",
            "Align a sequence (`file`: FASTA, GenBank, .dna, or a sequencing read from .ab1 or .scf; `document`: another document open in a tab, by its id, which is the only way in the web build, where the browser gives the page no file paths; or `residues`) to the open document, on both strands, across the origin of a circular document. `mode`: \"fit\" (default: the whole query somewhere in the document), \"local\" (best matching stretch) or \"global\" (end to end). Reports where it lies, identity, and every substitution, insertion and deletion in the document's coordinates (0-based); `rows` adds the aligned rows. A read from a trace file has its low-quality ends trimmed first (phred's modified Mott trimming, cutoff 0.05; `trim`: false keeps the whole read), reported as `trimmed`. Scores default to NCBI blastn's (2, -3, gaps 5 + 2 per base). Proteins: when the document is a protein, or with `feature` (a CDS of the open DNA, by id: its translation), the query is a protein aligned with BLOSUM62 (gaps 11 + 1 per residue, NCBI BLASTP's defaults); positions are then residues and `positives` counts identical or similar pairs.",
            align_command,
        )
        .menu("Actions/Align to This Sequence…")
        .example(json!({"residues": "ACGTACGTAAAACCCCTT"}))
        .example(json!({"file": "insert.fa", "mode": "local", "rows": true}))
        .example(json!({"file": "read.ab1", "trim": false}))
        .example(json!({"document": 1, "mode": "local"}))
        .example(json!({"feature": 1, "residues": "MKV", "mode": "local"})),
    ]
}

#[derive(Deserialize, JsonSchema, Clone, Copy, Default)]
#[serde(rename_all = "snake_case")]
enum AlignMode {
    #[default]
    Fit,
    Local,
    Global,
}

#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
struct AlignParams {
    /// The sequence to align, from a file (its first record).
    file: Option<String>,
    /// Or the sequence itself.
    residues: Option<String>,
    /// Or a document open in a tab, by its id.
    document: Option<u64>,
    /// Align a protein to the translation of this CDS feature of the open document.
    feature: Option<FeatureId>,
    #[serde(default)]
    mode: AlignMode,
    /// Include the aligned rows.
    #[serde(default)]
    rows: bool,
    /// Trim the low-quality ends of a read from a trace file first (default true).
    trim: Option<bool>,
    /// Score of a match (default 2).
    matched: Option<i32>,
    /// Score of a mismatch (default -3).
    mismatch: Option<i32>,
    /// Cost of opening a gap (default 5).
    gap_open: Option<i32>,
    /// Cost per gap base (default 2).
    gap_extend: Option<i32>,
}

fn align_command(s: &Session, p: AlignParams) -> Result<Value, CmdError> {
    let reference = s.molecule()?;
    let mut trimmed = None;
    let (name, query) = match (&p.file, &p.residues, p.document) {
        (Some(path), None, None) => {
            let bytes = s
                .platform()
                .read(path)
                .map_err(|e| CmdError::new(ErrorCode::Io, format!("cannot read {path}: {e}")))?;
            let m = snappiesnap_io::import(path, &bytes)?.molecule;
            let mut residues = m.seq().as_bytes().to_vec();
            if let Some(t) = m.trace()
                && p.trim.unwrap_or(true)
            {
                let q: Vec<Option<u8>> = t.calls().iter().map(|c| c.quality).collect();
                let keep = quality_trim(&q, DEFAULT_CUTOFF).ok_or_else(|| {
                    CmdError::new(
                        ErrorCode::InvalidInput,
                        format!("no stretch of the read {path} is of good enough quality to align"),
                    )
                    .with_hint("pass \"trim\": false to align the whole read anyway")
                })?;
                residues = residues.get(keep.clone()).unwrap_or(&[]).to_vec();
                trimmed =
                    Some(json!({"start": keep.start, "end": keep.end, "read_length": m.len()}));
            }
            (m.name().to_string(), residues)
        }
        (None, Some(r), None) => {
            let r: Vec<u8> = r.bytes().filter(|b| !b.is_ascii_whitespace()).collect();
            ("sequence".to_string(), r)
        }
        (None, None, Some(id)) => {
            let d = s.document_by_id(id).ok_or_else(|| {
                CmdError::new(ErrorCode::NotFound, format!("document {id} is not open"))
                    .with_hint("the `documents` command lists the open documents and their ids")
            })?;
            let m = &d.molecule;
            let mut residues = m.seq().as_bytes().to_vec();
            if let Some(t) = m.trace()
                && p.trim.unwrap_or(true)
            {
                let q: Vec<Option<u8>> = t.calls().iter().map(|c| c.quality).collect();
                let keep = quality_trim(&q, DEFAULT_CUTOFF).ok_or_else(|| {
                    CmdError::new(
                        ErrorCode::InvalidInput,
                        format!(
                            "no stretch of the read {} is of good enough quality to align",
                            m.name()
                        ),
                    )
                    .with_hint("pass \"trim\": false to align the whole read anyway")
                })?;
                residues = residues.get(keep.clone()).unwrap_or(&[]).to_vec();
                trimmed =
                    Some(json!({"start": keep.start, "end": keep.end, "read_length": m.len()}));
            }
            (m.name().to_string(), residues)
        }
        _ => {
            return Err(CmdError::invalid_params(
                "give the sequence to align as `file`, `document` or `residues` (one of them)",
            ));
        }
    };
    let protein_doc = reference.seq().alphabet() == Alphabet::Protein;
    if protein_doc || p.feature.is_some() {
        return protein(s, &p, &name, &query);
    }
    if let Some(b) = query.iter().find(|b| !b.is_ascii_alphabetic()) {
        return Err(CmdError::invalid_params(format!(
            "the sequence has a character that isn't a base: {:?}",
            char::from(*b)
        )));
    }
    let d = Scoring::default();
    let scoring = Scoring {
        matched: p.matched.unwrap_or(d.matched),
        mismatch: p.mismatch.unwrap_or(d.mismatch),
        gap_open: p.gap_open.unwrap_or(d.gap_open),
        gap_extend: p.gap_extend.unwrap_or(d.gap_extend),
    };
    let limit = 10_000;
    if [
        scoring.matched,
        scoring.mismatch,
        scoring.gap_open,
        scoring.gap_extend,
    ]
    .iter()
    .any(|v| v.abs() > limit)
    {
        return Err(CmdError::invalid_params(format!(
            "scores must be between -{limit} and {limit}"
        )));
    }
    let mode = match p.mode {
        AlignMode::Fit => Mode::Fit,
        AlignMode::Local => Mode::Local,
        AlignMode::Global => Mode::Global,
    };
    let a = align_to_molecule(reference, &query, mode, &scoring).map_err(|e| {
        CmdError::new(ErrorCode::InvalidInput, e)
            .with_hint("align a shorter sequence, or to a shorter document")
    })?;
    let mut v = serde_json::to_value(&a).unwrap_or(Value::Null);
    if let Some(o) = v.as_object_mut() {
        if !p.rows {
            o.remove("reference_row");
            o.remove("query_row");
        }
        o.insert("kind".into(), json!("dna"));
        o.insert("reference_length".into(), json!(reference.len()));
        o.insert("query_name".into(), json!(name));
        o.insert(
            "identity".into(),
            json!((a.identity * 10_000.0).round() / 10_000.0),
        );
        if let Some(t) = trimmed {
            o.insert("trimmed".into(), t);
        }
    }
    Ok(v)
}

/// A protein query against the open protein, or against a CDS translation.
fn protein(s: &Session, p: &AlignParams, name: &str, query: &[u8]) -> Result<Value, CmdError> {
    if p.matched.is_some() || p.mismatch.is_some() {
        return Err(CmdError::invalid_params(
            "proteins are scored with BLOSUM62: only `gap_open` and `gap_extend` can be set",
        ));
    }
    let (reference, reference_name) = match p.feature {
        Some(id) => {
            let m = s.molecule()?;
            let f = m
                .feature(id)
                .ok_or_else(|| CmdError::new(ErrorCode::NotFound, format!("no feature {id}")))?;
            let name = f.name.clone();
            (super::translate::feature_protein(s, id)?, name)
        }
        None => {
            let m = s.molecule()?;
            (m.seq().as_bytes().to_vec(), m.name().to_string())
        }
    };
    let clean = |r: &[u8]| -> Vec<u8> {
        let r: Vec<u8> = r.iter().map(|b| b.to_ascii_uppercase()).collect();
        r.strip_suffix(b"*").map(<[u8]>::to_vec).unwrap_or(r)
    };
    let (reference, query) = (clean(&reference), clean(query));
    if let Some(b) = query
        .iter()
        .find(|b| !(b.is_ascii_alphabetic() || **b == b'*'))
    {
        return Err(CmdError::invalid_params(format!(
            "the protein has a character that isn't an amino acid: {:?}",
            char::from(*b)
        )));
    }
    let d = ProteinScoring::default();
    let scoring = ProteinScoring {
        gap_open: p.gap_open.unwrap_or(d.gap_open),
        gap_extend: p.gap_extend.unwrap_or(d.gap_extend),
    };
    if scoring.gap_open.abs() > 10_000 || scoring.gap_extend.abs() > 10_000 {
        return Err(CmdError::invalid_params(
            "gap costs must be between -10000 and 10000",
        ));
    }
    let mode = match p.mode {
        AlignMode::Fit => Mode::Fit,
        AlignMode::Local => Mode::Local,
        AlignMode::Global => Mode::Global,
    };
    let a = align_protein_report(&reference, &query, mode, &scoring).map_err(|e| {
        CmdError::new(ErrorCode::InvalidInput, e).with_hint("align shorter proteins")
    })?;
    let columns = a.matches + a.mismatches + a.insertions + a.deletions;
    let similarity = a.positives.unwrap_or(0) as f64 / columns.max(1) as f64;
    let mut v = serde_json::to_value(&a).unwrap_or(Value::Null);
    if let Some(o) = v.as_object_mut() {
        o.remove("strand");
        if !p.rows {
            o.remove("reference_row");
            o.remove("query_row");
        }
        o.insert("kind".into(), json!("protein"));
        o.insert("reference_length".into(), json!(reference.len()));
        o.insert("query_name".into(), json!(name));
        o.insert("reference_name".into(), json!(reference_name));
        o.insert(
            "identity".into(),
            json!((a.identity * 10_000.0).round() / 10_000.0),
        );
        o.insert(
            "similarity".into(),
            json!((similarity * 10_000.0).round() / 10_000.0),
        );
    }
    Ok(v)
}
