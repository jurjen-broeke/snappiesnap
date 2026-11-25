//! Checking a construct with sequencing reads (Experiments › Check Reads).

use schemars::JsonSchema;
use serde::Deserialize;
use serde_json::{Value, json};
use snappiesnap_bio::reads::{MAX_READS, Read, check};

use crate::{CmdError, CommandSpec, ErrorCode, Session};

pub(super) fn specs() -> Vec<CommandSpec> {
    vec![
        CommandSpec::query(
            "check_reads",
            "Check Reads",
            "Check a construct with sequencing reads. `document`: the construct (an open document by id; default the active one). Reads: `files` (.ab1, .scf, FASTA, GenBank, .dna; first record), `documents` (open documents by id) and/or, with `matched`, the reads stored with the construct in its .dna file (kept as trimmed there). Each read is trimmed of its low-quality ends when it has qualities (`trim`, default true; phred's modified Mott, cutoff 0.05), aligned to the construct on both strands (the whole read inside it, across the origin of a circle), and kept if at least 80% identical. Returns each read's place, strand and differences with the quality of the read's own call, the stretches no read covers, and every difference pooled over the reads: `confirmed` (all covering reads show it, at least two), `one_read` or `disputed` (other covering reads show the construct's base). Positions are 0-based, half-open; `read_position` is in the untrimmed read, its own orientation. `rows`: give every placed read a `row` as well, its bases as they pair with the construct's top strand, one character per construct position from the read's `start` (`-` for a base the read is missing), with `row_quality` beside it, `row_positions` giving the untrimmed read's base behind each column (its own orientation; null for a gap) and `row_insertions` marking the columns the read has extra bases before; the rows of several reads then line up in one grid.",
            check_reads,
        )
        .example(json!({"files": ["read.ab1"]}))
        .example(json!({"documents": [1], "trim": false}))
        .example(json!({"files": ["read.ab1"], "rows": true})),
        CommandSpec::query(
            "matched_reads",
            "Matched Reads",
            "The sequencing reads matched to a construct and stored in its .dna file (as the reference application keeps them): each read's index, name, length, whether it is a trace, and the stretch kept after trimming (0-based, end exclusive). `document`: an open document's id (default the active one). Check them with check_reads `matched: true`; open one with open_matched_read.",
            list,
        )
        .example(json!({})),
        CommandSpec::session(
            "open_matched_read",
            "Open Matched Read",
            "Open a read matched to a construct (by its `index` in matched_reads) in a tab of its own: a trace opens with its chromatogram. `document`: the construct (default the active one).",
            open_read,
        )
        .needs_document(true)
        .example(json!({"index": 0})),
    ]
}

#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
struct ListParams {
    /// The construct: an open document's id (default the active one).
    document: Option<u64>,
}

fn list(s: &Session, p: ListParams) -> Result<Value, CmdError> {
    let reads = crate::matched::matched_reads(s, p.document)?;
    let out: Vec<Value> = reads
        .iter()
        .enumerate()
        .map(|(i, r)| match &r.molecule {
            Ok(m) => json!({
                "index": i,
                "name": r.name,
                "length": m.len(),
                "trace": m.trace().is_some(),
                "trimmed": r.trimmed,
            }),
            Err(e) => json!({"index": i, "name": r.name, "error": e}),
        })
        .collect();
    Ok(json!({ "reads": out }))
}

#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
struct OpenParams {
    /// The read's index in matched_reads.
    index: usize,
    /// The construct: an open document's id (default the active one).
    document: Option<u64>,
}

fn open_read(s: &mut Session, p: OpenParams) -> Result<Value, CmdError> {
    let reads = crate::matched::matched_reads(s, p.document)?;
    let r = reads.into_iter().nth(p.index).ok_or_else(|| {
        CmdError::new(ErrorCode::NotFound, format!("no matched read {}", p.index))
            .with_hint("matched_reads lists them")
    })?;
    let m = r
        .molecule
        .map_err(|e| CmdError::new(ErrorCode::InvalidInput, format!("{}: {e}", r.name)))?;
    let result = json!({"name": m.name(), "length": m.len(), "trace": m.trace().is_some()});
    let mut doc = crate::Document::new(m);
    doc.dirty = true;
    s.open_document(doc);
    Ok(result)
}

#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
struct Params {
    /// The construct: an open document's id (default the active document).
    document: Option<u64>,
    /// Read files.
    #[serde(default)]
    files: Vec<String>,
    /// Reads that are open documents, by id.
    #[serde(default)]
    documents: Vec<u64>,
    /// Trim low-quality ends of reads that have qualities (default true).
    #[serde(default = "yes")]
    trim: bool,
    /// Also check the reads matched to the construct in its .dna file (kept as trimmed there).
    #[serde(default)]
    matched: bool,
    /// Give each placed read its bases laid out against the construct (`row`).
    #[serde(default)]
    rows: bool,
}

fn yes() -> bool {
    true
}

/// A read from a molecule: its residues, and its qualities if it has a trace.
fn read_of(name: String, m: &snappiesnap_doc::Molecule) -> Read {
    Read {
        name,
        residues: m.seq().as_bytes().to_vec(),
        quality: m
            .trace()
            .map(|t| t.calls().iter().map(|c| c.quality).collect()),
        keep: None,
    }
}

fn check_reads(s: &Session, p: Params) -> Result<Value, CmdError> {
    let construct = match p.document {
        Some(id) => s.document_by_id(id).ok_or_else(|| {
            CmdError::new(ErrorCode::NotFound, format!("document {id} is not open"))
                .with_hint("the `documents` command lists the open documents and their ids")
        })?,
        None => s.document().ok_or_else(CmdError::no_document)?,
    };
    let matched = if p.matched {
        crate::matched::matched_reads(s, p.document)?
    } else {
        Vec::new()
    };
    let count = p
        .files
        .len()
        .saturating_add(p.documents.len())
        .saturating_add(matched.len());
    if count == 0 || count > MAX_READS {
        return Err(CmdError::invalid_params(format!(
            "give 1 to {MAX_READS} reads as `files` or `documents`, or check the construct's matched reads (`matched`)"
        )));
    }
    let mut reads = Vec::with_capacity(count);
    let mut sources = Vec::with_capacity(count);
    for path in &p.files {
        let bytes = s
            .platform()
            .read(path)
            .map_err(|e| CmdError::new(ErrorCode::Io, format!("cannot read {path}: {e}")))?;
        let m = snappiesnap_io::import(path, &bytes)
            .map_err(|e| CmdError::new(ErrorCode::InvalidInput, format!("{path}: {e}")))?
            .molecule;
        reads.push(read_of(m.name().to_string(), &m));
        sources.push(json!({"file": path, "trace": m.trace().is_some()}));
    }
    for id in &p.documents {
        let d = s.document_by_id(*id).ok_or_else(|| {
            CmdError::new(ErrorCode::NotFound, format!("document {id} is not open"))
        })?;
        reads.push(read_of(d.molecule.name().to_string(), &d.molecule));
        sources.push(json!({"document": id, "trace": d.molecule.trace().is_some()}));
    }
    for (i, r) in matched.iter().enumerate() {
        match &r.molecule {
            Ok(m) => {
                let mut read = read_of(r.name.clone(), m);
                read.keep = r.trimmed;
                reads.push(read);
                sources.push(json!({"matched": i, "trace": m.trace().is_some()}));
            }
            Err(e) => {
                return Err(CmdError::new(
                    ErrorCode::InvalidInput,
                    format!("matched read {}: {e}", r.name),
                ));
            }
        }
    }
    let c = check(&construct.molecule, &reads, p.trim, p.rows)
        .map_err(|e| CmdError::new(ErrorCode::InvalidInput, e))?;
    let round = |x: f64| (x * 10_000.0).round() / 10_000.0;
    let reads_json: Vec<Value> = c
        .reads
        .iter()
        .zip(sources)
        .map(|(r, source)| {
            let mut v = serde_json::to_value(r).unwrap_or(Value::Null);
            if let Some(o) = v.as_object_mut() {
                o.insert("identity".into(), json!(round(r.identity)));
                o.insert("source".into(), source);
            }
            v
        })
        .collect();
    let findings: Vec<Value> = c
        .findings
        .iter()
        .map(|f| {
            let mut v = serde_json::to_value(&f.difference).unwrap_or(Value::Null);
            if let Some(o) = v.as_object_mut() {
                o.insert("verdict".into(), json!(f.verdict));
                o.insert(
                    "reads".into(),
                    json!(
                        f.supported_by
                            .iter()
                            .map(
                                |(k, pos, q)| json!({"read": k, "read_position": pos, "quality": q})
                            )
                            .collect::<Vec<_>>()
                    ),
                );
                o.insert(
                    "covered_by".into(),
                    json!(
                        f.covered_by
                            .iter()
                            .map(|(k, pos)| json!({"read": k, "read_position": pos}))
                            .collect::<Vec<_>>()
                    ),
                );
            }
            v
        })
        .collect();
    let len = construct.molecule.len();
    Ok(json!({
        "construct": construct.molecule.name(),
        "length": len,
        "topology": construct.molecule.topology(),
        "reads": reads_json,
        "placed": c.reads.iter().filter(|r| r.error.is_none()).count(),
        "covered": c.covered,
        "covered_twice": c.covered_twice,
        "coverage": round(c.covered as f64 / len.max(1) as f64),
        "uncovered": c.uncovered.iter().map(|(a, b)| json!({"start": a, "end": b})).collect::<Vec<_>>(),
        "findings": findings,
    }))
}
