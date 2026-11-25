//! The CLI binary, run as a separate process.

// Test code: the never-crash lints apply to production code only (see clippy.toml).
#![allow(clippy::unwrap_used, clippy::expect_used, clippy::indexing_slicing)]

use std::io::Write;
use std::process::{Command, Output, Stdio};

use serde_json::Value;

fn cli(args: &[&str], stdin: &str) -> Output {
    let mut child = Command::new(env!("CARGO_BIN_EXE_snappiesnap-cli"))
        .args(args)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("start snappiesnap-cli");
    child
        .stdin
        .take()
        .unwrap()
        .write_all(stdin.as_bytes())
        .unwrap();
    child.wait_with_output().unwrap()
}

fn json_lines(out: &Output) -> Vec<Value> {
    String::from_utf8_lossy(&out.stdout)
        .lines()
        .map(|l| serde_json::from_str(l).unwrap())
        .collect()
}

#[test]
fn version_help_and_usage_errors() {
    let out = cli(&["--version"], "");
    assert!(out.status.success());
    assert!(String::from_utf8_lossy(&out.stdout).starts_with("snappiesnap-cli "));
    assert!(cli(&["help"], "").status.success());
    assert_eq!(cli(&["bogus"], "").status.code(), Some(2));
    assert_eq!(cli(&["info"], "").status.code(), Some(2));
    assert_eq!(
        cli(&["info", "/nonexistent/x.dna"], "").status.code(),
        Some(1)
    );
    assert_eq!(cli(&["convert", "a"], "").status.code(), Some(2));
    assert_eq!(
        cli(&["convert", "/nonexistent.dna", "x.gb"], "")
            .status
            .code(),
        Some(1)
    );
    assert_eq!(cli(&["run"], "").status.code(), Some(2));
    assert_eq!(cli(&["describe"], "").status.code(), Some(2));
}

#[test]
fn run_a_script_from_stdin() {
    let script = r#"
        {"command": "new", "params": {"name": "pCli", "residues": "ACGTACGT", "topology": "circular"}}
        {"command": "insert", "params": {"at": 0, "residues": "TT"}}
        {"method": "doc.get", "params": {"residues": false}}
    "#;
    let out = cli(&["run", "--script", "-"], script);
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    let lines = json_lines(&out);
    assert_eq!(lines.len(), 3);
    assert_eq!(lines[2]["result"]["molecule"]["length"], 10);
}

#[test]
fn a_failing_step_stops_the_script_unless_keep_going() {
    let script = "{\"command\": \"undo\"}\n{\"command\": \"new\"}\n";
    let out = cli(&["run", "--script", "-"], script);
    assert_eq!(out.status.code(), Some(1));
    assert_eq!(json_lines(&out).len(), 1);
    assert!(String::from_utf8_lossy(&out.stderr).contains("--keep-going"));

    let out = cli(&["run", "--script", "-", "--keep-going", "--quiet"], script);
    assert_eq!(out.status.code(), Some(1));
    let lines = json_lines(&out);
    assert_eq!(lines.len(), 1, "quiet prints only the failure");
    assert_eq!(lines[0]["error"]["code"], "no_document");
}

#[test]
fn missing_script_file_is_an_error() {
    let out = cli(&["run", "--script", "/nonexistent/x.jsonl"], "");
    assert_eq!(out.status.code(), Some(1));
    assert!(String::from_utf8_lossy(&out.stderr).contains("cannot read script"));
}

#[test]
fn mcp_over_pipes() {
    let input = concat!(
        r#"{"jsonrpc":"2.0","id":0,"method":"initialize","params":{"protocolVersion":"2025-06-18"}}"#,
        "\n",
        r#"{"jsonrpc":"2.0","method":"notifications/initialized"}"#,
        "\n",
        r#"{"jsonrpc":"2.0","id":1,"method":"tools/call","params":{"name":"new","arguments":{"residues":"ACGT"}}}"#,
        "\n",
    );
    let out = cli(&["mcp"], input);
    assert!(out.status.success());
    let lines = json_lines(&out);
    assert_eq!(
        lines.len(),
        2,
        "one response per request, none for the notification"
    );
    assert_eq!(lines[0]["result"]["serverInfo"]["name"], "snappiesnap");
    assert_eq!(lines[1]["result"]["structuredContent"]["length"], 4);
}

#[test]
fn commands_and_describe() {
    let out = cli(&["commands"], "");
    let v: Value = serde_json::from_slice(&out.stdout).unwrap();
    assert!(v["commands"].as_array().unwrap().len() >= 15);
    let out = cli(&["describe", "insert"], "");
    let v: Value = serde_json::from_slice(&out.stdout).unwrap();
    assert_eq!(v["command"]["id"], "insert");
    assert_eq!(cli(&["describe", "nope"], "").status.code(), Some(1));
}

#[test]
fn run_and_export_an_image() {
    let out = std::env::temp_dir().join(format!("snappiesnap-cli-{}.png", std::process::id()));
    let out_s = out.to_string_lossy().into_owned();
    let script = "{\"command\": \"new\", \"params\": {\"residues\": \"ACGTACGTACGTACGTACGT\", \"topology\": \"circular\"}}\n";
    let r = cli(
        &[
            "run", "--script", "-", "--quiet", "--export", &out_s, "--view", "sequence", "--size",
            "320x240",
        ],
        script,
    );
    assert!(r.status.success(), "{}", String::from_utf8_lossy(&r.stderr));
    let bytes = std::fs::read(&out).unwrap();
    let _ = std::fs::remove_file(&out);
    assert_eq!(&bytes[1..4], b"PNG");
    // Width and height in the IHDR chunk.
    assert_eq!(u32::from_be_bytes(bytes[16..20].try_into().unwrap()), 320);
    assert_eq!(u32::from_be_bytes(bytes[20..24].try_into().unwrap()), 240);
    for bad in [&["--size", "huge"][..], &["--view", "nope"], &["--theme"]] {
        let mut args = vec!["run", "--script", "-", "--export", &out_s];
        args.extend_from_slice(bad);
        assert_eq!(cli(&args, script).status.code(), Some(2), "{bad:?}");
    }
    assert_eq!(
        cli(&["render", "/nonexistent.dna", &out_s], "")
            .status
            .code(),
        Some(1)
    );
    assert_eq!(cli(&["render", "only-one-arg"], "").status.code(), Some(2));
}

#[test]
fn convert_between_formats() {
    let dir = std::env::temp_dir().join(format!("snappiesnap-convert-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let p = |name: &str| dir.join(name).to_string_lossy().into_owned();
    // Built with serde_json, never by formatting a path into JSON text: a Windows path is full
    // of backslashes, and `C:\Users\...` written straight into a JSON string is invalid JSON
    // (`\U` is not an escape), so the script would not parse there.
    let script = [
        serde_json::json!({
            "command": "new",
            "params": {
                "name": "pConv",
                "residues": "ACGTACGTACGTACGTACGT",
                "topology": "circular",
            }
        }),
        serde_json::json!({
            "command": "add_feature",
            "params": {
                "feature": {
                    "name": "f",
                    "kind": "CDS",
                    "strand": "forward",
                    "segments": [{ "range": { "start": 18, "end": 3 } }],
                }
            }
        }),
        serde_json::json!({ "command": "save", "params": { "path": p("a.dna") } }),
    ]
    .iter()
    .map(|step| serde_json::to_string(step).unwrap() + "\n")
    .collect::<String>();
    let r = cli(&["run", "--script", "-", "--quiet"], &script);
    assert!(r.status.success(), "{}", String::from_utf8_lossy(&r.stderr));
    let r = cli(&["convert", &p("a.dna"), &p("b.gb")], "");
    assert!(r.status.success(), "{}", String::from_utf8_lossy(&r.stderr));
    let gb = std::fs::read_to_string(p("b.gb")).unwrap();
    // Documents opened from a file are named after it.
    assert!(gb.starts_with("LOCUS       a "), "{gb}");
    assert!(
        gb.contains("     CDS             join(19..20,1..3)\n"),
        "{gb}"
    );
    let r = cli(&["convert", &p("b.gb"), &p("c.fa")], "");
    assert!(String::from_utf8_lossy(&r.stderr).contains("warning: FASTA stores only the sequence"));
    assert_eq!(
        std::fs::read_to_string(p("c.fa")).unwrap(),
        ">b  (20 bp)\nACGTACGTACGTACGTACGT\n"
    );
    let r = cli(&["convert", &p("b.gb"), &p("d.dna")], "");
    assert!(r.status.success());
    let info: Value = serde_json::from_slice(&cli(&["info", &p("d.dna")], "").stdout).unwrap();
    assert_eq!(info["features"][0]["segments"][0], "19..3");
    let _ = std::fs::remove_dir_all(&dir);
}
