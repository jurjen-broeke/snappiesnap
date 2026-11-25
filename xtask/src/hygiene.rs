//! Local-only folders never get committed (AGENTS.md › Committing).
//!
//! `.gitignore` must exclude `plan/`, `log/` and `corpus/`, and, in a git checkout, git must not
//! be tracking anything inside them.

use std::fs;
use std::process::Command;

use crate::util::{self, Result};

const LOCAL_ONLY: &[&str] = &["plan", "log", "corpus"];

pub fn run() -> Result {
    let root = util::root()?;
    let mut problems = Vec::new();

    let gitignore = fs::read_to_string(root.join(".gitignore")).unwrap_or_default();
    let ignored: Vec<&str> = gitignore
        .lines()
        .map(|l| l.trim().trim_start_matches('/').trim_end_matches('/'))
        .collect();
    for dir in LOCAL_ONLY {
        if !ignored.contains(dir) {
            problems.push(format!("`.gitignore` must contain `/{dir}/`"));
        }
    }

    // Outside a git checkout (or without git) there is nothing tracked to check.
    let tracked = Command::new("git")
        .current_dir(&root)
        .args(["ls-files", "--"])
        .args(LOCAL_ONLY)
        .output();
    if let Ok(out) = tracked
        && out.status.success()
    {
        for path in String::from_utf8_lossy(&out.stdout).lines() {
            problems.push(format!(
                "`{path}` is tracked by git but is local-only; run `git rm --cached` on it"
            ));
        }
    }
    util::verdict("hygiene", problems)
}
