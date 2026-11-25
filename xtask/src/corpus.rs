//! `cargo xtask corpus-test`: byte-exact round trips over the local corpus.
//!
//! The corpus (`corpus/`) is local-only and never committed, so this isn't part of `ci`. Run it
//! whenever you touch `formats` or `io` (AGENTS.md › Before you finish a task).

use crate::util::{self, Result};

/// Crates with a `corpus` feature, as `package/feature` pairs.
const CORPUS_FEATURES: &[&str] = &[
    "snappiesnap-formats/corpus",
    "snappiesnap-io/corpus",
    "snappiesnap-bio/corpus",
    "snappiesnap-engine/corpus",
];

pub fn run() -> Result {
    let root = util::root()?;
    if !root.join("corpus").is_dir() {
        return Err(
            "corpus/ is missing; it is local-only (see corpus/README.md in a checkout that has it)"
                .to_string(),
        );
    }
    let mut cmd = util::cargo()?;
    cmd.args(["test", "--features", &CORPUS_FEATURES.join(",")]);
    for feature in CORPUS_FEATURES {
        if let Some((package, _)) = feature.split_once('/') {
            cmd.args(["-p", package]);
        }
    }
    cmd.args(["--test", "corpus", "--", "--nocapture"]);
    util::run_status(cmd, "corpus tests")?;
    println!("xtask corpus-test: ok");
    Ok(())
}
