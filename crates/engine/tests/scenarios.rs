//! Replays every scenario in `tests/scenarios/*.jsonl` (repository root) in a fresh session.
//! Every step must succeed, including its `expect` / `expect_error` checks.

// Test code: the never-crash lints apply to production code only (see clippy.toml).
#![allow(clippy::unwrap_used, clippy::expect_used, clippy::panic)]

use std::io;
use std::path::PathBuf;
use std::sync::Arc;

use snappiesnap_engine::{DirEntry, NativeFs, Platform, Session, script};

fn scenario_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/scenarios")
}

/// The real file system (scenarios read files under `tests/scenarios/data`), but with a data
/// folder of the scenario's own, so saved enzyme sets, libraries and caches never reach the
/// user's real data folder or leak from one scenario into the next.
struct ScenarioFs {
    data: PathBuf,
}

impl Platform for ScenarioFs {
    fn read(&self, path: &str) -> io::Result<Vec<u8>> {
        NativeFs.read(path)
    }
    fn write(&self, path: &str, bytes: &[u8]) -> io::Result<()> {
        NativeFs.write(path, bytes)
    }
    fn today(&self) -> Option<(i32, u32, u32)> {
        NativeFs.today()
    }
    fn exists(&self, path: &str) -> bool {
        NativeFs.exists(path)
    }
    fn rename(&self, from: &str, to: &str) -> io::Result<()> {
        NativeFs.rename(from, to)
    }
    fn create_dir(&self, path: &str) -> io::Result<()> {
        NativeFs.create_dir(path)
    }
    fn list_dir(&self, path: &str) -> io::Result<Vec<DirEntry>> {
        NativeFs.list_dir(path)
    }
    fn data_dir(&self) -> Option<String> {
        Some(self.data.to_string_lossy().into_owned())
    }
}

#[test]
fn every_scenario_passes() {
    let mut files: Vec<PathBuf> = std::fs::read_dir(scenario_dir())
        .expect("tests/scenarios exists")
        .map(|e| e.unwrap().path())
        .filter(|p| p.extension().is_some_and(|x| x == "jsonl"))
        .collect();
    files.sort();
    assert!(!files.is_empty(), "no scenarios found");

    let data_root =
        std::env::temp_dir().join(format!("snappiesnap-scenarios-{}", std::process::id()));
    let mut failures = Vec::new();
    for file in &files {
        let text = std::fs::read_to_string(file).unwrap();
        let name = file.file_name().unwrap().to_string_lossy().into_owned();
        let data = data_root.join(&name);
        let _ = std::fs::remove_dir_all(&data);
        std::fs::create_dir_all(&data).unwrap();
        let mut session = Session::new().with_platform(Arc::new(ScenarioFs { data }));
        let summary = script::run(&mut session, &text, |step| {
            if let Err(e) = &step.result {
                failures.push(format!("{name}:{}: {e}\n    {}", step.line, step.source));
            }
            true
        });
        assert!(summary.steps > 0, "{name} has no steps");
    }
    let _ = std::fs::remove_dir_all(&data_root);
    assert!(
        failures.is_empty(),
        "scenario failures:\n{}",
        failures.join("\n")
    );
}
