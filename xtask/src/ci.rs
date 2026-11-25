//! `cargo xtask ci`: the full quality gate (AGENTS.md › Before you finish a task).
//!
//! Every step runs even if an earlier one fails, so one run reports every problem.

use crate::util::{self, Result};
use crate::{assets, hygiene, layers, licenses, lints, wasm};

type Step = (&'static str, fn() -> Result);

const STEPS: &[Step] = &[
    ("fmt", fmt),
    ("clippy", clippy),
    ("test", test),
    ("layers", layers::run),
    ("lints", lints::run),
    ("licenses", licenses::run),
    ("assets", assets::run),
    ("hygiene", hygiene::run),
    ("wasm", wasm::run),
];

pub fn run() -> Result {
    let mut failed = Vec::new();
    for (name, step) in STEPS {
        println!("\n=== xtask ci: {name} ===");
        if let Err(msg) = step() {
            eprintln!("{msg}");
            failed.push(*name);
        }
    }
    println!();
    if failed.is_empty() {
        println!("xtask ci: all {} steps passed", STEPS.len());
        Ok(())
    } else {
        Err(format!("xtask ci: failed steps: {}", failed.join(", ")))
    }
}

fn fmt() -> Result {
    let mut cmd = util::cargo()?;
    cmd.args(["fmt", "--all", "--check"]);
    util::run_status(cmd, "cargo fmt --check")
}

fn clippy() -> Result {
    let mut cmd = util::cargo()?;
    cmd.args([
        "clippy",
        "--workspace",
        "--all-targets",
        // Also lints code behind features, such as the local-only corpus tests.
        "--all-features",
        "--",
        "-D",
        "warnings",
    ]);
    util::run_status(cmd, "cargo clippy")
}

fn test() -> Result {
    let mut cmd = util::cargo()?;
    cmd.args(["test", "--workspace"]);
    util::run_status(cmd, "cargo test")
}
