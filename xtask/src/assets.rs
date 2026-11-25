//! Asset attribution (AGENTS.md › Assets and brand).
//!
//! Every file under `assets/` must be named in `ATTRIBUTION.md` (path, title, author, source,
//! licence). Licence files that sit next to an asset don't need their own row.

use std::fs;
use std::path::{Path, PathBuf};

use crate::util::{self, Result};

/// Directory depth limit, so a symlink loop can't make the walk run forever.
const MAX_DEPTH: usize = 32;

pub fn run() -> Result {
    let root = util::root()?;
    let attribution_path = root.join("ATTRIBUTION.md");
    let attribution = fs::read_to_string(&attribution_path)
        .map_err(|e| format!("cannot read {}: {e}", attribution_path.display()))?;

    let mut files = Vec::new();
    let assets = root.join("assets");
    if assets.is_dir() {
        walk(&assets, 0, &mut files)?;
    }

    let mut problems = Vec::new();
    for file in files {
        let rel = file
            .strip_prefix(&root)
            .unwrap_or(&file)
            .to_string_lossy()
            .replace('\\', "/");
        if is_exempt(&file) {
            continue;
        }
        if !attribution.contains(&rel) {
            problems.push(format!("`{rel}` has no row in ATTRIBUTION.md"));
        }
    }
    util::verdict("assets", problems)
}

fn walk(dir: &Path, depth: usize, out: &mut Vec<PathBuf>) -> Result {
    if depth > MAX_DEPTH {
        return Err(format!("{} is nested too deeply", dir.display()));
    }
    let entries = fs::read_dir(dir).map_err(|e| format!("cannot list {}: {e}", dir.display()))?;
    for entry in entries {
        let entry = entry.map_err(|e| format!("cannot list {}: {e}", dir.display()))?;
        let path = entry.path();
        if path.is_dir() {
            walk(&path, depth.saturating_add(1), out)?;
        } else {
            out.push(path);
        }
    }
    Ok(())
}

fn is_exempt(path: &Path) -> bool {
    let name = path
        .file_name()
        .map(|n| n.to_string_lossy().to_ascii_uppercase())
        .unwrap_or_default();
    name == ".DS_STORE" || name == ".GITKEEP" || name.starts_with("LICENSE") || name == "README.MD"
}
