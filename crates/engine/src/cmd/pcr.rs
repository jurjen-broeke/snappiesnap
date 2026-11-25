//! PCR simulation and gels.

use schemars::JsonSchema;
use serde::Deserialize;
use serde_json::{Value, json};
use snappiesnap_bio::pcr::{DEFAULT_MAX_PRODUCT, Product, pcr, product_molecule};

use crate::gel::{LaneSpec, bands, primer_residues};
use crate::{CmdError, CommandSpec, Document, ErrorCode, Session};

pub(super) fn specs() -> Vec<CommandSpec> {
    vec![
        CommandSpec::query(
            "pcr",
            "Simulate PCR",
            "The products two primers make from the open molecule (a primer's name, or a sequence). Either primer may prime on either strand; tails and mismatches in the primers become part of the product. Shortest first.",
            pcr_command,
        )
        .example(json!({"forward": "f", "reverse": "r"}))
        .example(json!({"forward": "ACGTACGTAC", "reverse": "GTACGTACGT", "secondary": true, "max_length": 5000})),
        CommandSpec::session(
            "pcr_open",
            "Open PCR Product",
            "Open a PCR product as a new document, with the template's features inside it and the primers as primer_bind features.",
            pcr_open,
        )
        .example(json!({"forward": "f", "reverse": "r", "discard_changes": true})),
        CommandSpec::query(
            "gel",
            "Simulate Gel",
            "Bands of an agarose gel: lanes of a ladder, the uncut molecule, digests, PCR products or given sizes. Migration is log-linear across the separation range of the agarose percentage (Sambrook & Russell 2001); brightness follows DNA mass. `y` runs from 0 (wells) to 1 (front).",
            gel,
        )
        .example(json!({"lanes": [{"kind": "ladder", "name": "1 kb"}, {"kind": "uncut"}]}))
        .example(json!({"agarose": 1.5, "lanes": [{"kind": "digest", "enzymes": ["EcoRI"]}, {"kind": "pcr", "forward": "f", "reverse": "r"}, {"kind": "sizes", "sizes": [500, 1200]}]})),
    ]
}

#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
struct PcrParams {
    /// A primer of this molecule (by name), or a primer sequence 5'->3'.
    forward: String,
    /// The other primer, likewise.
    reverse: String,
    /// Also use weaker (secondary) binding sites, under the session's binding rules
    /// (primer_binding).
    #[serde(default)]
    secondary: bool,
    /// Longest product reported (default 50,000).
    max_length: Option<usize>,
}

fn run_pcr(s: &Session, p: &PcrParams) -> Result<(Vec<Product>, [String; 2]), CmdError> {
    let m = s.molecule()?;
    let bad = |e: String| {
        CmdError::new(ErrorCode::InvalidInput, e)
            .with_hint("name a primer of this molecule or give its sequence (at least 10 bases)")
    };
    let (fname, f) = primer_residues(m, &p.forward).map_err(bad)?;
    let (rname, r) = primer_residues(m, &p.reverse).map_err(bad)?;
    let max = p.max_length.unwrap_or(DEFAULT_MAX_PRODUCT);
    let secondary = p.secondary.then(|| s.primer_binding());
    Ok((pcr(m, [&f, &r], secondary, max), [fname, rname]))
}

fn pcr_command(s: &Session, p: PcrParams) -> Result<Value, CmdError> {
    let (products, names) = run_pcr(s, &p)?;
    let list: Vec<Value> = products
        .iter()
        .map(|x| {
            json!({
                "length": x.length,
                "forward_primer": names.get(x.forward),
                "reverse_primer": names.get(x.reverse),
                "forward_site": x.forward_site,
                "reverse_site": x.reverse_site,
                "template": x.template,
                "residues": String::from_utf8_lossy(&x.seq),
            })
        })
        .collect();
    Ok(json!({ "products": list }))
}

#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
struct OpenParams {
    forward: String,
    reverse: String,
    #[serde(default)]
    secondary: bool,
    /// Which product (0 = shortest; default 0).
    #[serde(default)]
    product: usize,
    /// Accepted for older scripts; the product opens in a new tab, so nothing is replaced.
    #[serde(default)]
    #[allow(dead_code)]
    discard_changes: bool,
}

fn pcr_open(s: &mut Session, p: OpenParams) -> Result<Value, CmdError> {
    let params = PcrParams {
        forward: p.forward,
        reverse: p.reverse,
        secondary: p.secondary,
        max_length: None,
    };
    let (products, names) = run_pcr(s, &params)?;
    let product = products.get(p.product).ok_or_else(|| {
        CmdError::new(
            ErrorCode::NotFound,
            format!("no product {} ({} made)", p.product, products.len()),
        )
    })?;
    let template = s.molecule()?;
    let name = format!("{} PCR {} bp", template.name(), product.length);
    let name: String = name.chars().take(200).collect();
    let mut molecule = product_molecule(template, product, [&names[0], &names[1]], &name)?;
    super::cloning::record(
        &mut molecule,
        "amplifyFragment",
        &format!("primers {} and {}", names[0], names[1]),
        vec![snappiesnap_doc::HistoryNode::input(template)],
        vec![super::cloning::amplified(product, template.len())],
    )?;
    let primers: Vec<(String, Vec<u8>)> = [&params.forward, &params.reverse]
        .iter()
        .filter_map(|t| primer_residues(template, t).ok())
        .collect();
    let primers: Vec<(&str, &[u8])> = primers
        .iter()
        .map(|(n, r)| (n.as_str(), r.as_slice()))
        .collect();
    super::cloning::set_primers(&mut molecule, &primers, true)?;
    let result = json!({ "name": molecule.name(), "length": molecule.len(), "features": molecule.features().len() });
    let mut doc = Document::new(molecule);
    // A new, unsaved result: closing it should ask first.
    doc.dirty = true;
    s.open_document(doc);
    Ok(result)
}

#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
struct GelParams {
    /// Agarose, % (default 1).
    #[serde(default = "one")]
    agarose: f64,
    lanes: Vec<LaneSpec>,
}

pub(super) fn one() -> f64 {
    1.0
}

/// Most lanes on a gel.
pub(super) const MAX_LANES: usize = 50;

fn gel(s: &Session, p: GelParams) -> Result<Value, CmdError> {
    if !p.agarose.is_finite() || !(0.1..=5.0).contains(&p.agarose) {
        return Err(CmdError::invalid_params("agarose must be 0.1-5 %"));
    }
    if p.lanes.is_empty() || p.lanes.len() > MAX_LANES {
        return Err(CmdError::invalid_params("give 1 to 50 lanes"));
    }
    let m = s.molecule()?;
    let mut lanes = Vec::new();
    for spec in &p.lanes {
        let b = bands(m, s.enzymes(), spec, p.agarose)
            .map_err(|e| CmdError::new(ErrorCode::InvalidInput, e))?;
        lanes.push(json!({ "title": spec.title(), "bands": bands_json(&b) }));
    }
    Ok(json!({ "agarose": p.agarose, "lanes": lanes }))
}

/// Bands as reported by the gel commands: positions to 0.001, brightness to 0.01.
pub(super) fn bands_json(bands: &[snappiesnap_bio::gel::Band]) -> Vec<Value> {
    bands
        .iter()
        .map(|band| {
            json!({
                "bp": band.bp,
                "y": (band.y * 1000.0).round() / 1000.0,
                "intensity": (band.intensity * 100.0).round() / 100.0,
                "fragments": band.fragments.iter().map(|f| json!({"bp": f.bp, "label": f.label})).collect::<Vec<_>>(),
            })
        })
        .collect()
}
