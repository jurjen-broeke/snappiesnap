//! Read-only commands.

use schemars::JsonSchema;
use serde::Deserialize;
use serde_json::{Value, json};
use snappiesnap_doc::Range;

use crate::{CmdError, CommandSpec, Session};

pub(super) fn specs() -> Vec<CommandSpec> {
    vec![
        CommandSpec::query(
            "get_residues",
            "Copy Residues",
            "The residues in a range (across the origin if it wraps).",
            get_residues,
        )
        .example(json!({"range": {"start": 10, "end": 2}})),
    ]
}

#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
struct RangeParams {
    range: Range,
}

fn get_residues(s: &Session, p: RangeParams) -> Result<Value, CmdError> {
    let bytes = s.molecule()?.residues(p.range)?;
    Ok(json!({ "residues": String::from_utf8_lossy(&bytes) }))
}
