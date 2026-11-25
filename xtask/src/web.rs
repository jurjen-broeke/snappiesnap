//! Build the browser version as static files, and pack them for a web server.
//!
//! This is the same build the release workflow's `web` job runs, done locally so a feature can
//! be tried in the browser as soon as it lands: `trunk build
//! --release` writes `dist/web`, and the folder is packed as
//! `dist/SnappieSnap-<version>-web.tar.gz` to copy to the server (see docs/web.md).
//!
//! It needs `trunk` (`cargo install trunk`) and the `wasm32-unknown-unknown` target; both are
//! reported by name when they are missing, because neither can be installed from here.

use std::path::Path;
use std::process::Command;

use crate::util::{self, Result};

/// Where the web app's `Trunk.toml` lives, relative to the workspace root.
const WEB_APP_DIR: &str = "apps/snappiesnap-web";
/// Where trunk writes the built files (`dist` in that `Trunk.toml`).
const DIST: &str = "dist/web";

pub fn run() -> Result {
    let root = util::root()?;
    if which("trunk").is_none() {
        return Err(
            "trunk is not installed: `cargo install trunk` (it builds the browser version; \
             see docs/web.md)"
                .to_string(),
        );
    }
    let version = version()?;

    let mut build = Command::new("trunk");
    build.current_dir(root.join(WEB_APP_DIR));
    build.args(["build", "--release", "--locked"]);
    // An absolute target dir of its own, inside the ignored `target/`. Trunk runs cargo from
    // the app's folder, so a relative CARGO_TARGET_DIR inherited from the caller (an agent's
    // `target/agent-<name>`) would otherwise leave a build tree inside `apps/`.
    build.env("CARGO_TARGET_DIR", root.join("target").join("web-build"));
    util::run_status(build, "trunk build --release")?;

    let dist = root.join(DIST);
    if !dist.join("index.html").is_file() {
        return Err(format!(
            "trunk reported success but {DIST}/index.html is missing"
        ));
    }

    // One archive of the folder, as the CI job makes it: unpack it anywhere a web server
    // serves static files.
    let archive = format!("SnappieSnap-{version}-web.tar.gz");
    let mut tar = Command::new("tar");
    tar.current_dir(&root);
    tar.args(["-czf", &format!("dist/{archive}"), "-C", "dist", "web"]);
    util::run_status(tar, "tar")?;

    let size = std::fs::metadata(root.join("dist").join(&archive))
        .map(|m| m.len())
        .unwrap_or(0);
    println!(
        "xtask web: dist/{archive} ({:.1} MB), unpacked in {DIST}/",
        size as f64 / 1_000_000.0
    );
    Ok(())
}

/// The workspace version, for the archive's name.
fn version() -> Result<String> {
    let metadata = util::metadata()?;
    metadata
        .get("packages")
        .and_then(|p| p.as_array())
        .and_then(|packages| {
            packages
                .iter()
                .find(|p| p.get("name").and_then(|n| n.as_str()) == Some("snappiesnap"))
        })
        .and_then(|p| p.get("version"))
        .and_then(|v| v.as_str())
        .map(str::to_string)
        .ok_or_else(|| "cannot read the workspace version from cargo metadata".to_string())
}

/// Where a program is on PATH, if it is there at all.
fn which(program: &str) -> Option<std::path::PathBuf> {
    let path = std::env::var_os("PATH")?;
    std::env::split_paths(&path)
        .map(|dir| dir.join(program))
        .find(|p| is_program(p))
}

fn is_program(path: &Path) -> bool {
    path.is_file()
}
