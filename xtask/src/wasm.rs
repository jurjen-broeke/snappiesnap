//! Never break WASM (AGENTS.md › Architecture rules): every L0–L5 crate, and the web app, must
//! build for `wasm32-unknown-unknown`. Clippy runs on that build too, so code that only exists
//! on the web is held to the same never-panic lints as the rest.

use crate::layers;
use crate::util::{self, Result};

const TARGET: &str = "wasm32-unknown-unknown";

/// The browser app (`apps/snappiesnap-web`, built with trunk; see docs/web.md).
const WEB_APP: &str = "snappiesnap-web";

pub fn run() -> Result {
    let mut cmd = util::cargo()?;
    cmd.args(["clippy", "--target", TARGET]);
    for name in layers::core_crates().into_iter().chain([WEB_APP]) {
        cmd.args(["-p", name]);
    }
    cmd.args(["--", "-D", "warnings"]);
    // Separate target dir so wasm checks don't invalidate the native build cache.
    if let Ok(root) = util::root() {
        cmd.env("CARGO_TARGET_DIR", root.join("target").join("wasm-check"));
    }
    util::run_status(cmd, "wasm check")?;
    println!("xtask wasm: ok");
    Ok(())
}
