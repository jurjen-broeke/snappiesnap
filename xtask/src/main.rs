//! `cargo xtask <task>`: SnappieSnap developer tasks and quality gates.
//!
//! Every gate reports problems as text and a non-zero exit code; none of them panic.

mod assets;
mod ci;
mod corpus;
mod enzymes;
mod features;
mod fpbase;
mod hygiene;
mod icons;
mod layers;
mod licenses;
mod lints;
mod util;
mod wasm;
mod web;

use std::process::ExitCode;

const USAGE: &str = "usage: cargo xtask <task>

tasks:
  ci        run every gate below plus fmt, clippy and tests
  layers    crate layering (and no UI crates below L5)
  lints     every crate opts in to the workspace never-crash lints
  licenses  every dependency has a permissively licensed dependency graph
  assets    every file under assets/ has a row in ATTRIBUTION.md
  hygiene   local-only folders (plan/, log/, corpus/) stay out of git
  wasm      L0-L5 crates build for wasm32-unknown-unknown

local only (not part of ci):
  web          build the browser version (trunk) into dist/web and pack it for a server
  corpus-test  byte-exact round trips over the local corpus/ (never committed)
  enzymes      rebuild the built-in enzyme table from a REBASE withrefm file
  features     rebuild assets/data/features/insdc_features.json from public NCBI records
  fpbase       rebuild assets/data/features/fpbase_proteins.json from FPbase (CC BY-SA 4.0)
  icons        render the icon pack in assets/brand/icons/ from assets/brand/snappiesnap.svg";

fn main() -> ExitCode {
    let task = std::env::args().nth(1);
    let result = match task.as_deref() {
        Some("ci") => ci::run(),
        Some("layers") => layers::run(),
        Some("lints") => lints::run(),
        Some("licenses") => licenses::run(),
        Some("assets") => assets::run(),
        Some("hygiene") => hygiene::run(),
        Some("wasm") => wasm::run(),
        Some("web") => web::run(),
        Some("corpus-test") => corpus::run(),
        Some("enzymes") => enzymes::run(std::env::args().nth(2).as_deref()),
        Some("features") => features::run(),
        Some("fpbase") => fpbase::run(),
        Some("icons") => icons::run(),
        None | Some("help" | "--help" | "-h") => {
            println!("{USAGE}");
            return ExitCode::SUCCESS;
        }
        Some(other) => Err(format!("unknown task `{other}`\n\n{USAGE}")),
    };
    match result {
        Ok(()) => ExitCode::SUCCESS,
        Err(msg) => {
            eprintln!("error: {msg}");
            ExitCode::FAILURE
        }
    }
}
