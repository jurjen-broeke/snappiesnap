//! Multiple sequence alignment of sequences from files, pasted text and the open document.

use schemars::JsonSchema;
use serde::Deserialize;
use serde_json::{Value, json};
use snappiesnap_bio::msa::{Kind, MAX_SEQUENCES, align_many};

use crate::{CmdError, CommandSpec, ErrorCode, Session};

pub(super) fn specs() -> Vec<CommandSpec> {
    vec![CommandSpec::query(
        "align_sequences",
        "Align Multiple Sequences",
        "Multiple sequence alignment (progressive: k-mer distances, UPGMA guide tree, profile alignment with affine gaps; blastn scores for nucleotides, BLOSUM62 for proteins). Sequences come from `files` (every record of a FASTA file; one molecule from other formats), from `documents` (open documents by id, which is the only way in the web build, where the browser gives the page no file paths), from `text` (pasted: FASTA, or plain sequences separated by blank lines), from `sequences` ({name, residues}), and with `include_document` the open document first. `kind` is \"auto\" (default), \"dna\" or \"protein\". Returns the aligned rows (input order, `-` for gaps), a consensus and each row's identity to it.",
        msa_command,
    )
    .needs_document(false)
    .menu("Actions/Align Multiple Sequences…")
    .example(json!({"text": ">a\nACGTTGCAAGGCTTACGATCG\n>b\nACGTTGCAGGCTTACGATCG\n>c\nACGTTGCAAGGCTTTACGATCG"}))
    .example(json!({"include_document": true, "files": ["insert.fa"], "kind": "dna"}))
    .example(json!({"documents": [1, 2], "kind": "dna"}))]
}

#[derive(Deserialize, JsonSchema, Clone, Copy, Default, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
enum KindParam {
    #[default]
    Auto,
    Dna,
    Protein,
}

#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
struct Named {
    #[serde(default)]
    name: String,
    residues: String,
}

#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
struct MsaParams {
    #[serde(default)]
    files: Vec<String>,
    /// Documents open in tabs, by their ids.
    #[serde(default)]
    documents: Vec<u64>,
    /// Pasted sequences: FASTA, or plain sequences separated by blank lines.
    text: Option<String>,
    #[serde(default)]
    sequences: Vec<Named>,
    /// Put the open document first.
    #[serde(default)]
    include_document: bool,
    #[serde(default)]
    kind: KindParam,
}

/// Sequences in pasted text: FASTA records, or blocks of letters separated by blank lines
/// (numbers and spaces, as in copied GenBank or sequence-view text, are dropped).
pub(crate) fn sequences_in_text(text: &str) -> Result<Vec<(String, Vec<u8>)>, CmdError> {
    let text = text.trim_start_matches('\u{feff}').trim();
    if text.starts_with('>') {
        let fasta = snappiesnap_formats::fasta::FastaFile::parse(text.as_bytes())
            .map_err(|e| CmdError::invalid_params(format!("pasted FASTA: {e}")))?;
        return Ok(fasta
            .records
            .into_iter()
            .enumerate()
            .map(|(i, r)| {
                let name = if r.name().is_empty() {
                    format!("sequence {}", i + 1)
                } else {
                    r.name().to_string()
                };
                (name, r.residues)
            })
            .collect());
    }
    let mut out = Vec::new();
    for block in text.split("\n\n").flat_map(|b| b.split("\r\n\r\n")) {
        let residues: Vec<u8> = block
            .bytes()
            .filter(|b| b.is_ascii_alphabetic() || *b == b'*')
            .collect();
        if !residues.is_empty() {
            out.push((format!("sequence {}", out.len() + 1), residues));
        }
    }
    Ok(out)
}

fn msa_command(s: &Session, p: MsaParams) -> Result<Value, CmdError> {
    let mut seqs: Vec<(String, Vec<u8>)> = Vec::new();
    if p.include_document {
        let m = s.molecule()?;
        seqs.push((m.name().to_string(), m.seq().as_bytes().to_vec()));
    }
    for id in &p.documents {
        let d = s.document_by_id(*id).ok_or_else(|| {
            CmdError::new(ErrorCode::NotFound, format!("document {id} is not open"))
                .with_hint("the `documents` command lists the open documents and their ids")
        })?;
        seqs.push((
            d.molecule.name().to_string(),
            d.molecule.seq().as_bytes().to_vec(),
        ));
    }
    for path in &p.files {
        let bytes = s
            .platform()
            .read(path)
            .map_err(|e| CmdError::new(ErrorCode::Io, format!("cannot read {path}: {e}")))?;
        match snappiesnap_formats::fasta::FastaFile::parse(&bytes) {
            Ok(fasta) if fasta.records.len() > 1 => {
                seqs.extend(
                    fasta
                        .records
                        .into_iter()
                        .map(|r| (r.name().to_string(), r.residues)),
                );
            }
            _ => {
                let m = snappiesnap_io::import(path, &bytes)?.molecule;
                seqs.push((m.name().to_string(), m.seq().as_bytes().to_vec()));
            }
        }
    }
    if let Some(text) = &p.text {
        seqs.extend(sequences_in_text(text)?);
    }
    for (i, n) in p.sequences.iter().enumerate() {
        let name = if n.name.trim().is_empty() {
            format!("sequence {}", i + 1)
        } else {
            n.name.clone()
        };
        let residues: Vec<u8> = n
            .residues
            .bytes()
            .filter(|b| !b.is_ascii_whitespace())
            .collect();
        seqs.push((name, residues));
    }
    if seqs.len() < 2 {
        return Err(CmdError::invalid_params(format!(
            "give at least two sequences ({} found)",
            seqs.len()
        )));
    }
    if seqs.len() > MAX_SEQUENCES {
        return Err(CmdError::invalid_params(format!(
            "at most {MAX_SEQUENCES} sequences at once"
        )));
    }
    for (name, r) in &seqs {
        if let Some(b) = r.iter().find(|b| !(b.is_ascii_alphabetic() || **b == b'*')) {
            return Err(CmdError::invalid_params(format!(
                "{name} has a character that isn't a residue: {:?}",
                char::from(*b)
            )));
        }
    }
    let refs: Vec<&[u8]> = seqs.iter().map(|(_, r)| r.as_slice()).collect();
    let kind = match p.kind {
        KindParam::Auto => Kind::guess(&refs),
        KindParam::Dna => Kind::Dna,
        KindParam::Protein => Kind::Protein,
    };
    let msa = align_many(&refs, kind).map_err(|e| {
        CmdError::new(ErrorCode::InvalidInput, e).with_hint("align fewer or shorter sequences")
    })?;
    let rows: Vec<Value> = seqs
        .iter()
        .zip(&msa.rows)
        .zip(&msa.identity)
        .map(|(((name, r), row), id)| {
            json!({
                "name": name,
                "length": r.len(),
                "row": row,
                "identity": (id * 10_000.0).round() / 10_000.0,
            })
        })
        .collect();
    Ok(json!({
        "kind": msa.kind,
        "columns": msa.consensus.len(),
        "identical_columns": msa.identical_columns,
        "consensus": msa.consensus,
        "rows": rows,
    }))
}
