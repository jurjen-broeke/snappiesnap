//! Never-crash lints (AGENTS.md › Never crash).
//!
//! The workspace manifest must keep the required lint levels, and every workspace crate must
//! opt in with `[lints] workspace = true`.

use std::fs;

use toml::Value;

use crate::util::{self, Result};

const REQUIRED_RUST: &[(&str, &str)] = &[("unsafe_code", "forbid")];

const REQUIRED_CLIPPY: &[&str] = &[
    "unwrap_used",
    "expect_used",
    "panic",
    "unimplemented",
    "todo",
    "unreachable",
    "indexing_slicing",
];

pub fn run() -> Result {
    let mut problems = Vec::new();
    let root = util::root()?;

    let workspace = read_toml(&root.join("Cargo.toml").to_string_lossy())?;
    let lints = workspace.get("workspace").and_then(|w| w.get("lints"));
    for (lint, level) in REQUIRED_RUST {
        if level_of(lints, "rust", lint) != Some(level) {
            problems.push(format!(
                "workspace Cargo.toml must set [workspace.lints.rust] {lint} = \"{level}\""
            ));
        }
    }
    for lint in REQUIRED_CLIPPY {
        if !matches!(level_of(lints, "clippy", lint), Some("deny" | "forbid")) {
            problems.push(format!(
                "workspace Cargo.toml must set [workspace.lints.clippy] {lint} = \"deny\""
            ));
        }
    }

    let meta = util::metadata()?;
    let members = util::workspace_members(&meta);
    for pkg in util::packages(&meta)
        .iter()
        .filter(|p| members.contains(&p.id))
    {
        let manifest = read_toml(pkg.manifest_path)?;
        let opted_in = manifest
            .get("lints")
            .and_then(|l| l.get("workspace"))
            .and_then(Value::as_bool)
            == Some(true);
        if !opted_in {
            problems.push(format!(
                "`{}` must opt in to the workspace lints: add `[lints] workspace = true` to {}",
                pkg.name, pkg.manifest_path
            ));
        }
    }
    util::verdict("lints", problems)
}

fn read_toml(path: &str) -> Result<Value> {
    let text = fs::read_to_string(path).map_err(|e| format!("cannot read {path}: {e}"))?;
    text.parse::<Value>()
        .map_err(|e| format!("cannot parse {path}: {e}"))
}

/// The level of `lint` in a `[workspace.lints.<group>]` table: either `lint = "deny"` or
/// `lint = { level = "deny", ... }`.
fn level_of<'a>(lints: Option<&'a Value>, group: &str, lint: &str) -> Option<&'a str> {
    let entry = lints?.get(group)?.get(lint)?;
    entry
        .as_str()
        .or_else(|| entry.get("level").and_then(Value::as_str))
}
