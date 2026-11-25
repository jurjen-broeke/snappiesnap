//! Shared helpers: workspace root, running cargo, reading `cargo metadata`.

use std::path::{Path, PathBuf};
use std::process::Command;

use serde_json::Value;

pub type Result<T = ()> = std::result::Result<T, String>;

/// The workspace root (the parent of the `xtask` crate).
pub fn root() -> Result<PathBuf> {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .map(Path::to_path_buf)
        .ok_or_else(|| "cannot locate the workspace root".to_string())
}

/// A `cargo` command run from the workspace root.
pub fn cargo() -> Result<Command> {
    let bin = std::env::var("CARGO").unwrap_or_else(|_| "cargo".to_string());
    let mut cmd = Command::new(bin);
    cmd.current_dir(root()?);
    Ok(cmd)
}

/// Run a command with inherited output; `Err` if it can't start or exits non-zero.
pub fn run_status(mut cmd: Command, what: &str) -> Result {
    let status = cmd
        .status()
        .map_err(|e| format!("{what}: could not start: {e}"))?;
    if status.success() {
        Ok(())
    } else {
        Err(format!("{what}: failed ({status})"))
    }
}

/// Parsed `cargo metadata --format-version 1` (including dependencies and the resolve graph).
pub fn metadata() -> Result<Value> {
    let out = cargo()?
        .args(["metadata", "--format-version", "1"])
        .output()
        .map_err(|e| format!("cargo metadata: could not start: {e}"))?;
    if !out.status.success() {
        return Err(format!(
            "cargo metadata failed:\n{}",
            String::from_utf8_lossy(&out.stderr)
        ));
    }
    serde_json::from_slice(&out.stdout).map_err(|e| format!("cannot parse cargo metadata: {e}"))
}

/// A package from the metadata, with only the fields the gates need.
pub struct Package<'a> {
    pub id: &'a str,
    pub name: &'a str,
    pub license: Option<&'a str>,
    pub manifest_path: &'a str,
    pub raw: &'a Value,
}

pub fn packages(meta: &Value) -> Vec<Package<'_>> {
    let Some(list) = meta.get("packages").and_then(Value::as_array) else {
        return Vec::new();
    };
    list.iter()
        .filter_map(|p| {
            Some(Package {
                id: p.get("id")?.as_str()?,
                name: p.get("name")?.as_str()?,
                license: p.get("license").and_then(Value::as_str),
                manifest_path: p.get("manifest_path")?.as_str()?,
                raw: p,
            })
        })
        .collect()
}

/// Ids of the workspace members.
pub fn workspace_members(meta: &Value) -> Vec<&str> {
    meta.get("workspace_members")
        .and_then(Value::as_array)
        .map(|ids| ids.iter().filter_map(Value::as_str).collect())
        .unwrap_or_default()
}

/// Print a gate's verdict; `Err` listing every problem if there are any.
pub fn verdict(gate: &str, problems: Vec<String>) -> Result {
    if problems.is_empty() {
        println!("xtask {gate}: ok");
        Ok(())
    } else {
        Err(format!(
            "xtask {gate}: {} problem(s)\n  - {}",
            problems.len(),
            problems.join("\n  - ")
        ))
    }
}
