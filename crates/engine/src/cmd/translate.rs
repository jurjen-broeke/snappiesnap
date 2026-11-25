//! Translation, open reading frames and protein properties.

use schemars::JsonSchema;
use serde::Deserialize;
use serde_json::{Value, json};
use snappiesnap_bio::codes::CODES;
use snappiesnap_bio::translate::{
    OrfOptions, code, find_orfs, protein_properties, strand_residues, translate,
};
use snappiesnap_doc::{FeatureId, Range, Strand};

use super::residue_bytes;
use crate::{CmdError, CommandSpec, ErrorCode, Session};

pub(super) fn specs() -> Vec<CommandSpec> {
    vec![
        CommandSpec::query(
            "translate",
            "Translate",
            "Translate a range of the molecule (on the top or bottom strand) or given residues, with an NCBI genetic code (default 1, standard). Ambiguous codons give X unless all their readings agree.",
            translate_command,
        )
        .needs_document(false)
        .example(json!({"range": {"start": 0, "end": 12}}))
        .example(json!({"residues": "ATGGCCTGA", "table": 2})),
        CommandSpec::query(
            "orfs",
            "Find ORFs",
            "Open reading frames (start codon to stop codon) on both strands, longest first; the longest ORF per stop codon. Defaults: at least 75 amino acids, ATG starts, standard code.",
            orfs_command,
        )
        .example(json!({}))
        .example(json!({"min_length": 3, "only_atg": false, "table": 11})),
        CommandSpec::query(
            "protein_properties",
            "Protein Properties",
            "Length, average mass, isoelectric point (Bjellqvist et al. 1993 pKs) and extinction coefficient at 280 nm (Pace et al. 1995) of a protein, given as residues or as a CDS feature of the molecule (translated with its /transl_table and /codon_start).",
            properties_command,
        )
        .needs_document(false)
        .example(json!({"protein": "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEG"}))
        .example(json!({"feature": 1})),
        CommandSpec::query(
            "genetic_codes",
            "Genetic Codes",
            "The NCBI genetic code tables: number and name.",
            codes_command,
        )
        .needs_document(false)
        .example(json!({})),
    ]
}

fn table(id: Option<u8>) -> Result<&'static snappiesnap_bio::codes::GeneticCode, CmdError> {
    let id = id.unwrap_or(1);
    code(id).ok_or_else(|| {
        CmdError::invalid_params(format!("no genetic code {id}"))
            .with_hint("genetic_codes lists them")
    })
}

#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
struct TranslateParams {
    /// A range of the open molecule…
    range: Option<Range>,
    /// …read on this strand ("forward" or "reverse"; default forward)…
    #[serde(default)]
    strand: Option<Strand>,
    /// …or these residues instead.
    residues: Option<String>,
    /// NCBI genetic code (default 1).
    table: Option<u8>,
    /// Read a first codon that can start as Met (default true).
    as_start: Option<bool>,
}

fn translate_command(s: &Session, p: TranslateParams) -> Result<Value, CmdError> {
    let code = table(p.table)?;
    let nt = match (p.residues, p.range) {
        (Some(r), None) => residue_bytes(&r),
        (None, Some(range)) => {
            let m = s.molecule()?;
            let range = range.validate(m.len(), m.topology())?;
            strand_residues(m, range, p.strand.unwrap_or(Strand::Forward))
                .ok_or_else(|| CmdError::invalid_params("range outside the molecule"))?
        }
        _ => {
            return Err(CmdError::invalid_params(
                "give either `range` or `residues`",
            ));
        }
    };
    let protein = translate(&nt, code, p.as_start.unwrap_or(true));
    Ok(json!({
        "protein": String::from_utf8_lossy(&protein),
        "codons": protein.len(),
        "table": code.id,
    }))
}

#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
struct OrfParams {
    /// Shortest ORF in amino acids (default 75).
    min_length: Option<usize>,
    /// Only ATG starts (default true); otherwise every start codon of the code.
    only_atg: Option<bool>,
    /// Search the reverse strand too (default true).
    both_strands: Option<bool>,
    /// NCBI genetic code (default 1).
    table: Option<u8>,
    /// Most ORFs returned (default 100).
    limit: Option<usize>,
}

fn orfs_command(s: &Session, p: OrfParams) -> Result<Value, CmdError> {
    let d = OrfOptions::default();
    let code = table(p.table)?;
    let opts = OrfOptions {
        min_aa: p.min_length.unwrap_or(d.min_aa),
        code: code.id,
        only_atg: p.only_atg.unwrap_or(d.only_atg),
        both_strands: p.both_strands.unwrap_or(d.both_strands),
    };
    let m = s.molecule()?;
    let orfs = find_orfs(m, &opts);
    let total = orfs.len();
    let list: Vec<Value> = orfs
        .iter()
        .take(p.limit.unwrap_or(100).min(10_000))
        .map(
            |o| json!({"range": o.range, "strand": o.strand, "frame": o.frame, "length": o.length}),
        )
        .collect();
    Ok(json!({ "total": total, "orfs": list }))
}

#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
struct PropertiesParams {
    /// One-letter amino acids…
    protein: Option<String>,
    /// …or a CDS feature of the open molecule, by id.
    feature: Option<FeatureId>,
}

/// The protein a CDS feature codes for: segments joined, its strand, /codon_start and
/// /transl_table honoured.
pub(crate) fn feature_protein(s: &Session, id: FeatureId) -> Result<Vec<u8>, CmdError> {
    let m = s.molecule()?;
    let f = m
        .feature(id)
        .ok_or_else(|| CmdError::new(ErrorCode::NotFound, format!("no feature {id}")))?;
    let q = |n: &str| {
        f.qualifiers
            .iter()
            .find(|q| q.name == n)
            .and_then(|q| q.values.first())
            .map(|v| v.trim().to_string())
    };
    let start: usize = q("codon_start").and_then(|v| v.parse().ok()).unwrap_or(1);
    let code = table(q("transl_table").and_then(|v| v.parse().ok()))?;
    let mut nt = Vec::new();
    for seg in &f.segments {
        nt.extend(
            strand_residues(m, seg.range, f.strand)
                .ok_or_else(|| CmdError::invalid_params("feature outside the molecule"))?,
        );
    }
    let nt = nt.get(start.clamp(1, 3) - 1..).unwrap_or_default();
    Ok(translate(nt, code, true))
}

fn properties_command(s: &Session, p: PropertiesParams) -> Result<Value, CmdError> {
    let protein: Vec<u8> = match (p.protein, p.feature) {
        (Some(text), None) => text
            .bytes()
            .filter(|b| !b.is_ascii_whitespace())
            .map(|b| b.to_ascii_uppercase())
            .collect(),
        (None, Some(id)) => feature_protein(s, id)?,
        _ => {
            return Err(CmdError::invalid_params(
                "give either `protein` or `feature`",
            ));
        }
    };
    if protein.len() > 1_000_000 {
        return Err(CmdError::invalid_params("protein too long"));
    }
    let pp = protein_properties(&protein);
    let trimmed = protein.strip_suffix(b"*").unwrap_or(&protein);
    Ok(json!({
        "protein": String::from_utf8_lossy(trimmed),
        "length": pp.length,
        "mass_da": (pp.mass * 100.0).round() / 100.0,
        "unknown_residues": pp.unknown_residues,
        "isoelectric_point": (pp.isoelectric_point * 100.0).round() / 100.0,
        "extinction_280": {"cystines": pp.extinction_oxidised, "reduced": pp.extinction_reduced},
        "internal_stops": trimmed.iter().filter(|&&b| b == b'*').count(),
    }))
}

fn codes_command(_: &Session, _: crate::NoParams) -> Result<Value, CmdError> {
    let list: Vec<Value> = CODES
        .iter()
        .map(|c| json!({"table": c.id, "name": c.name}))
        .collect();
    Ok(json!({ "codes": list }))
}
