//! TA and TOPO cloning: a PCR product into a T-vector or TOPO vector.

use schemars::JsonSchema;
use serde::Deserialize;
use serde_json::{Value, json};
use snappiesnap_bio::topo::{TopoKind, topo_clone};
use snappiesnap_doc::HistoryNode;

use super::cloning::{load_source, open_product, pcr_molecule, record, summary};
use crate::{CmdError, CommandSpec, ErrorCode, Session};

pub(super) fn specs() -> Vec<CommandSpec> {
    vec![
        CommandSpec::session(
            "topo_clone",
            "TA/TOPO Cloning",
            "Put a linear insert (a file, the open document, or the PCR product of two of its primers) into a vector opened at `position` (a circular vector; a linear one is taken as already opened). `kind`: \"ta\" (T-vector or TOPO TA: 3' T overhangs pair with the A that Taq adds, giving T + insert + A, either orientation), \"blunt\" (blunt TOPO, either orientation) or \"directional\" (directional TOPO: only orientations in which the insert starts with CACC). With `preview`, lists the constructs; otherwise opens construct `product` (default 0, forward first).",
            topo_command,
        )
        .menu("Actions/TA-TOPO Cloning…")
        .example(json!({"vector": "dest.gb", "insert": "insert.fa", "position": 10, "preview": true}))
        .example(json!({"vector": "dest.gb", "pcr": {"forward": "f", "reverse": "r"}, "position": 10, "kind": "blunt"})),
    ]
}

/// Which kind of TA/TOPO cloning.
#[derive(Deserialize, JsonSchema, Clone, Copy, Default)]
#[serde(rename_all = "snake_case")]
enum Kind {
    #[default]
    Ta,
    Blunt,
    Directional,
}

#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
struct PcrSpec {
    forward: String,
    reverse: String,
}

#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
struct TopoParams {
    /// Vector file (default: the active document).
    vector: Option<String>,
    /// The vector as an open document's id instead of a file (see `documents`).
    vector_document: Option<u64>,
    /// Insert file (default: the active document).
    insert: Option<String>,
    /// The insert as an open document's id instead of a file.
    insert_document: Option<u64>,
    /// Use the PCR product of these two primers on the insert (a primer's name or a sequence).
    pcr: Option<PcrSpec>,
    /// Where a circular vector is opened: between bases, 0 to its length.
    position: Option<usize>,
    /// "ta" (default), "blunt" or "directional".
    #[serde(default)]
    kind: Kind,
    #[serde(default)]
    preview: bool,
    /// Which construct to open (default 0).
    #[serde(default)]
    product: usize,
}

fn topo_command(s: &mut Session, p: TopoParams) -> Result<Value, CmdError> {
    if p.vector.is_none()
        && p.insert.is_none()
        && p.pcr.is_none()
        && p.vector_document.is_none()
        && p.insert_document.is_none()
    {
        return Err(CmdError::invalid_params(
            "give `vector` or `insert` as a file or an open document's id, or a `pcr` (the \
             other is the active document)",
        ));
    }
    let vector = load_source(s, p.vector.as_deref(), p.vector_document, "the vector")?;
    let source = load_source(s, p.insert.as_deref(), p.insert_document, "the insert")?;
    let insert = match &p.pcr {
        Some(pcr) => pcr_molecule(&source, &pcr.forward, &pcr.reverse)?,
        None => source,
    };
    let kind = match p.kind {
        Kind::Ta => TopoKind::Ta,
        Kind::Blunt => TopoKind::Blunt,
        Kind::Directional => TopoKind::Directional,
    };
    let name: String = format!("{} + {}", vector.name(), insert.name())
        .chars()
        .take(200)
        .collect();
    let constructs = topo_clone(&vector, p.position, &insert, kind, &name).map_err(|e| {
        CmdError::new(ErrorCode::InvalidInput, e).with_hint(
            "the insert must be linear; a circular vector needs `position`; directional inserts start with CACC",
        )
    })?;
    if p.preview {
        let list: Vec<Value> = constructs
            .iter()
            .map(|c| {
                json!({
                    "length": c.molecule.len(),
                    "reversed": c.reversed,
                    "description": c.description,
                    "molecule": c.molecule,
                })
            })
            .collect();
        return Ok(json!({ "constructs": list }));
    }
    let c = constructs.get(p.product).ok_or_else(|| {
        CmdError::new(
            ErrorCode::NotFound,
            format!("no construct {} ({} possible)", p.product, constructs.len()),
        )
    })?;
    let mut m = c.molecule.clone();
    let (operation, insert_summary) = match kind {
        TopoKind::Ta => ("taCloning", "insert"),
        TopoKind::Blunt | TopoKind::Directional => ("topoCloning", "insert"),
    };
    let opened = p.position.unwrap_or(vector.len());
    record(
        &mut m,
        operation,
        &c.description,
        vec![HistoryNode::input(&vector), HistoryNode::input(&insert)],
        vec![
            summary("insertAt", &[("val1", opened.to_string())]),
            summary(
                insert_summary,
                &[
                    ("val1", "0".to_string()),
                    ("val2", insert.len().saturating_sub(1).to_string()),
                ],
            ),
        ],
    )?;
    let mut out = open_product(s, m, true)?;
    if let Some(o) = out.as_object_mut() {
        o.insert("description".into(), json!(c.description));
        o.insert("reversed".into(), json!(c.reversed));
    }
    Ok(out)
}
