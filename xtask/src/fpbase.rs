//! `cargo xtask fpbase`: build the built-in fluorescent-protein set from FPbase.
//!
//! Downloads FPbase's protein list (https://www.fpbase.org/api/proteins/) once with `curl` and a
//! generic User-Agent (cached in `target/fpbase/all.json`). It keeps every protein with a
//! sequence of the 20 standard amino acids, and writes `assets/data/features/fpbase_proteins.json`
//! with the name, the protein sequence, the DOI of the primary reference (the original authors,
//! whom FPbase asks users to credit) and the FPbase page. FPbase data is licensed CC BY-SA 4.0
//! (https://www.fpbase.org/terms/); the file keeps that licence.

use std::fs;

use serde_json::{Value, json};

use crate::util::{self, Result};

const API: &str = "https://www.fpbase.org/api/proteins/?format=json";
const AMINO_ACIDS: &[u8] = b"ACDEFGHIKLMNPQRSTVWY";

pub fn run() -> Result {
    let root = util::root()?;
    let cache = root.join("target/fpbase");
    fs::create_dir_all(&cache).map_err(|e| format!("cannot create {}: {e}", cache.display()))?;
    let path = cache.join("all.json");
    if !path.exists() {
        let status = std::process::Command::new("curl")
            .args(["-sSf", "-A", "SnappieSnap-dev", "-o"])
            .arg(&path)
            .arg(API)
            .status()
            .map_err(|e| format!("cannot run curl: {e}"))?;
        if !status.success() {
            let _ = fs::remove_file(&path);
            return Err(format!("FPbase download failed ({status})"));
        }
    }
    let all: Value = serde_json::from_str(
        &fs::read_to_string(&path).map_err(|e| format!("cannot read {}: {e}", path.display()))?,
    )
    .map_err(|e| format!("FPbase list: {e}"))?;
    let list = all.as_array().ok_or("FPbase list: not an array")?;
    let mut out = Vec::new();
    let mut skipped = 0usize;
    for p in list {
        let text = |k: &str| {
            p.get(k)
                .and_then(Value::as_str)
                .unwrap_or("")
                .trim()
                .to_string()
        };
        let (name, seq) = (text("name"), text("seq").to_ascii_uppercase());
        if name.is_empty()
            || !(50..=2000).contains(&seq.len())
            || !seq.bytes().all(|b| AMINO_ACIDS.contains(&b))
        {
            skipped += 1;
            continue;
        }
        out.push(json!({
            "name": name,
            "protein": seq,
            "doi": text("doi"),
            "genbank": text("genbank"),
            "page": format!("https://www.fpbase.org/protein/{}/", text("slug")),
        }));
    }
    out.sort_by(|a, b| a["name"].as_str().cmp(&b["name"].as_str()));
    let dataset = json!({
        "about": "Fluorescent proteins from FPbase (https://www.fpbase.org), licensed CC BY-SA 4.0 (https://www.fpbase.org/terms/). Credit each protein's original authors: its primary reference is `doi`. Built by `cargo xtask fpbase`; see SOURCES.md.",
        "licence": "CC-BY-SA-4.0",
        "proteins": out,
    });
    let target = root.join("assets/data/features/fpbase_proteins.json");
    let text = serde_json::to_string(&dataset).map_err(|e| e.to_string())?;
    // One protein per line: small diffs when the set is rebuilt.
    let text = text.replace("},{\"name\"", "},\n{\"name\"");
    fs::write(&target, format!("{text}\n"))
        .map_err(|e| format!("cannot write {}: {e}", target.display()))?;
    println!(
        "wrote {} proteins ({skipped} skipped) to {}",
        out.len(),
        target.display()
    );
    util::verdict("fpbase", Vec::new())
}
