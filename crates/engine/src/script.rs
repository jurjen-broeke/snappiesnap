//! Command scripts: JSON lines replayed against a session (`snappiesnap-cli run --script`,
//! scenario tests, and later `--record` sessions).
//!
//! Each non-blank line that isn't a comment (`#` or `//`) is one step:
//!
//! - `{"command": "insert", "params": {...}}` runs a command (same as `engine.execute`);
//! - `{"method": "doc.get", "params": {...}}` calls a protocol method.
//!
//! A step may also check its outcome, which turns a script into a test:
//!
//! - `"expect": {"path": "/molecule/length", "eq": 10}` (or a list of such checks) inspects the
//!   result by JSON pointer with `eq`, `len` or `contains`;
//! - `"expect_error": "invalid_input"` requires the step to fail with that error code.
//!
//! A failed check is an `assertion_failed` error. `"note"` holds free text.

use serde_json::Value;

use crate::{CmdError, ErrorCode, Session, api};

/// Longest script we accept, in bytes.
pub const MAX_SCRIPT: usize = 64 * 1024 * 1024;

/// What one step did.
#[derive(Debug)]
pub struct Step {
    /// 1-based line number in the script.
    pub line: usize,
    /// The step as written.
    pub source: String,
    pub result: Result<Value, CmdError>,
}

/// Totals for a script run.
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub struct Summary {
    pub steps: usize,
    pub failed: usize,
    /// True if the run stopped early because `on_step` returned false.
    pub stopped: bool,
}

/// Run `text` step by step. `on_step` sees every step and returns whether to continue.
pub fn run(session: &mut Session, text: &str, mut on_step: impl FnMut(&Step) -> bool) -> Summary {
    let mut summary = Summary::default();
    if text.len() > MAX_SCRIPT {
        let step = Step {
            line: 0,
            source: String::new(),
            result: Err(CmdError::invalid_params(format!(
                "script is larger than {MAX_SCRIPT} bytes"
            ))),
        };
        summary.failed = 1;
        summary.stopped = !on_step(&step);
        return summary;
    }
    for (index, raw) in text.lines().enumerate() {
        let source = raw.trim();
        if source.is_empty() || source.starts_with('#') || source.starts_with("//") {
            continue;
        }
        let step = Step {
            line: index.saturating_add(1),
            source: source.to_string(),
            result: run_step(session, source),
        };
        summary.steps = summary.steps.saturating_add(1);
        if step.result.is_err() {
            summary.failed = summary.failed.saturating_add(1);
        }
        if !on_step(&step) {
            summary.stopped = true;
            break;
        }
    }
    summary
}

/// Keys a step may have. Anything else is a typo and an error.
const STEP_KEYS: &[&str] = &[
    "command",
    "method",
    "params",
    "expect",
    "expect_error",
    "note",
];

fn run_step(session: &mut Session, source: &str) -> Result<Value, CmdError> {
    let step: Value = serde_json::from_str(source)
        .map_err(|e| CmdError::invalid_params(format!("not a JSON object: {e}")))?;
    let Some(fields) = step.as_object() else {
        return Err(CmdError::invalid_params("a step must be a JSON object"));
    };
    if let Some(key) = fields.keys().find(|k| !STEP_KEYS.contains(&k.as_str())) {
        return Err(CmdError::invalid_params(format!(
            "unknown step key {key:?}; allowed: {}",
            STEP_KEYS.join(", ")
        )));
    }
    let params = step.get("params").cloned().unwrap_or(Value::Null);
    let result = match (
        step.get("command").and_then(Value::as_str),
        step.get("method").and_then(Value::as_str),
    ) {
        (Some(command), None) => session.execute(command, params),
        (None, Some(method)) => api::call(session, method, params),
        _ => {
            return Err(CmdError::new(
                ErrorCode::InvalidParams,
                "each step needs exactly one of \"command\" or \"method\"",
            )
            .with_hint(
                r#"for example {"command": "insert", "params": {"at": 0, "residues": "ACGT"}}"#,
            ));
        }
    };
    match step.get("expect_error") {
        Some(expected) => check_error(result, expected),
        None => {
            let value = result?;
            if let Some(expect) = step.get("expect") {
                check_expectations(&value, expect)?;
            }
            Ok(value)
        }
    }
}

/// `"expect_error": "invalid_input"`: the step must fail with that error code.
fn check_error(result: Result<Value, CmdError>, expected: &Value) -> Result<Value, CmdError> {
    let Some(expected) = expected.as_str() else {
        return Err(CmdError::invalid_params(
            "\"expect_error\" must be an error code such as \"invalid_input\"",
        ));
    };
    match result {
        Ok(value) => Err(assertion(format!(
            "expected error {expected:?}, but the step succeeded with {value}"
        ))),
        Err(e) => {
            let code = serde_json::to_value(e.code).unwrap_or(Value::Null);
            if code == expected {
                Ok(serde_json::json!({ "expected_error": e }))
            } else {
                Err(assertion(format!(
                    "expected error {expected:?}, got {code}: {e}"
                )))
            }
        }
    }
}

/// `"expect"`: one check or a list of checks, each `{"path": <JSON pointer>, <op>: <value>}`
/// with `eq` (equal), `len` (length of an array, string or object) or `contains` (substring or
/// array element).
fn check_expectations(value: &Value, expect: &Value) -> Result<(), CmdError> {
    let checks = match expect {
        Value::Array(list) => list.iter().collect(),
        other => vec![other],
    };
    for check in checks {
        let path = check.get("path").and_then(Value::as_str).unwrap_or("");
        let actual = value
            .pointer(path)
            .ok_or_else(|| assertion(format!("nothing at {path:?} in {value}")))?;
        if let Some(want) = check.get("eq")
            && actual != want
        {
            return Err(assertion(format!("{path:?} is {actual}, expected {want}")));
        }
        if let Some(want) = check.get("len") {
            let len = match actual {
                Value::Array(a) => Some(a.len()),
                Value::String(s) => Some(s.chars().count()),
                Value::Object(o) => Some(o.len()),
                _ => None,
            };
            if len.map(Value::from) != Some(want.clone()) {
                return Err(assertion(format!(
                    "{path:?} has length {len:?}, expected {want}"
                )));
            }
        }
        if let Some(want) = check.get("contains") {
            let ok = match (actual, want) {
                (Value::String(s), Value::String(w)) => s.contains(w.as_str()),
                (Value::Array(a), w) => a.contains(w),
                _ => false,
            };
            if !ok {
                return Err(assertion(format!(
                    "{path:?} is {actual}, which doesn't contain {want}"
                )));
            }
        }
        if ["eq", "len", "contains"]
            .iter()
            .all(|op| check.get(op).is_none())
        {
            return Err(CmdError::invalid_params(
                "an expectation needs one of \"eq\", \"len\" or \"contains\"",
            ));
        }
    }
    Ok(())
}

fn assertion(message: String) -> CmdError {
    CmdError::new(ErrorCode::AssertionFailed, message)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn runs_commands_and_methods_and_reports_failures() {
        let mut s = Session::new();
        let script = r#"
            # build a tiny plasmid
            {"command": "new", "params": {"residues": "ACGT", "topology": "circular"}}
            // comments of both kinds are skipped
            {"command": "insert", "params": {"at": 4, "residues": "GG"}}
            {"method": "doc.get", "params": {"residues": false}}
            {"command": "insert", "params": {"at": 99, "residues": "A"}}
            not json
            {"command": "x", "method": "y"}
        "#;
        let mut lines = Vec::new();
        let summary = run(&mut s, script, |step| {
            lines.push((step.line, step.result.is_ok()));
            true
        });
        assert_eq!(
            summary,
            Summary {
                steps: 6,
                failed: 3,
                stopped: false
            }
        );
        assert_eq!(lines[0], (3, true));
        assert_eq!(s.molecule().map(|m| m.len()), Ok(6));
    }

    #[test]
    fn stops_when_asked() {
        let mut s = Session::new();
        let summary = run(
            &mut s,
            "{\"command\": \"undo\"}\n{\"command\": \"undo\"}",
            |_| false,
        );
        assert_eq!(summary.steps, 1);
        assert!(summary.stopped);
    }

    fn outcome(script: &str) -> Vec<Result<Value, CmdError>> {
        let mut s = Session::new();
        let mut out = Vec::new();
        run(&mut s, script, |step| {
            out.push(step.result.clone());
            true
        });
        out
    }

    #[test]
    fn expectations_pass_and_fail() {
        let r = outcome(
            r#"
            {"command": "new", "params": {"residues": "ACGT"}, "expect": {"path": "/length", "eq": 4}}
            {"method": "doc.get", "expect": [{"path": "/molecule/seq/residues", "contains": "CG"}, {"path": "/molecule/features", "len": 0}]}
            {"command": "insert", "params": {"at": 9, "residues": "A"}, "expect_error": "invalid_input"}
            {"method": "doc.get", "expect": {"path": "/molecule/length", "eq": 5}}
            {"command": "insert", "params": {"at": 0, "residues": "A"}, "expect_error": "invalid_input"}
            {"command": "insert", "params": {"at": 99, "residues": "A"}, "expect_error": "not_found"}
            {"method": "doc.get", "expect": {"path": "/nope", "eq": 1}}
            {"method": "doc.get", "expect": {"path": "/open"}}
            {"method": "doc.get", "expected": {"path": "/open", "eq": true}}
        "#,
        );
        let codes: Vec<Option<ErrorCode>> =
            r.iter().map(|r| r.as_ref().err().map(|e| e.code)).collect();
        use ErrorCode::*;
        assert_eq!(
            codes,
            vec![
                None,
                None,
                None,
                Some(AssertionFailed),
                Some(AssertionFailed),
                Some(AssertionFailed),
                Some(AssertionFailed),
                Some(InvalidParams),
                Some(InvalidParams),
            ]
        );
    }
}
