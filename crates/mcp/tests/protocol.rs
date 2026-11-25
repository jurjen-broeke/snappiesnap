//! The MCP server: handshake, tool coverage, tool calls, and hostile input.

// Test code: the never-crash lints apply to production code only (see clippy.toml).
#![allow(clippy::unwrap_used, clippy::indexing_slicing)]

use std::collections::BTreeSet;

use proptest::prelude::*;
use serde_json::{Value, json};
use snappiesnap_engine::{MemoryFs, Session};
use snappiesnap_mcp::{Server, serve};

/// In-memory files, so fuzzed `open`/`save` tool calls never touch the real file system.
fn server() -> Server {
    Server::new(Session::new().with_platform(std::sync::Arc::new(MemoryFs::new())))
}

fn request(s: &mut Server, id: u64, method: &str, params: Value) -> Value {
    let line = json!({"jsonrpc": "2.0", "id": id, "method": method, "params": params});
    let out = s.handle_line(&line.to_string()).unwrap();
    serde_json::from_str(&out).unwrap()
}

fn call(s: &mut Server, name: &str, arguments: Value) -> Value {
    request(
        s,
        1,
        "tools/call",
        json!({"name": name, "arguments": arguments}),
    )["result"]
        .clone()
}

#[test]
fn handshake() {
    let mut s = server();
    let r = request(
        &mut s,
        0,
        "initialize",
        json!({"protocolVersion": "2025-03-26", "capabilities": {},
               "clientInfo": {"name": "test", "version": "0"}}),
    );
    assert_eq!(r["result"]["protocolVersion"], "2025-03-26");
    assert_eq!(r["result"]["serverInfo"]["name"], "snappiesnap");
    assert!(r["result"]["capabilities"]["tools"].is_object());
    // Unknown versions get our newest.
    let r = request(
        &mut s,
        1,
        "initialize",
        json!({"protocolVersion": "1999-01-01"}),
    );
    assert_eq!(r["result"]["protocolVersion"], "2025-06-18");
    // Notifications get no response.
    let note = json!({"jsonrpc": "2.0", "method": "notifications/initialized"});
    assert!(s.handle_line(&note.to_string()).is_none());
    assert_eq!(request(&mut s, 2, "ping", json!({}))["result"], json!({}));
}

/// Every registered command is a tool with an object schema, and tool names are unique.
#[test]
fn every_command_is_a_tool() {
    let mut s = server();
    let tools = request(&mut s, 1, "tools/list", json!({}))["result"]["tools"].clone();
    let tools = tools.as_array().unwrap();
    let names: BTreeSet<&str> = tools.iter().map(|t| t["name"].as_str().unwrap()).collect();
    assert_eq!(names.len(), tools.len(), "duplicate tool names");
    for c in s.session().registry().iter() {
        let tool = tools.iter().find(|t| t["name"] == c.id);
        let tool = tool.unwrap_or_else(|| panic!("command {} has no tool", c.id));
        assert_eq!(tool["inputSchema"]["type"], "object", "{}", c.id);
        assert!(!tool["description"].as_str().unwrap().is_empty());
    }
    for name in [
        "doc_get",
        "commands_list",
        "commands_describe",
        "log_tail",
        "app_version",
    ] {
        assert!(names.contains(name), "{name}");
    }
    for name in &names {
        assert!(
            name.bytes().all(|b| b.is_ascii_alphanumeric() || b == b'_'),
            "{name}"
        );
    }
}

#[test]
fn agents_can_build_and_inspect_a_plasmid() {
    let mut s = server();
    let r = call(
        &mut s,
        "new",
        json!({"name": "pAgent", "residues": "ACGTACGTAC", "topology": "circular"}),
    );
    assert_eq!(r["isError"], false);
    assert_eq!(r["structuredContent"]["length"], 10);

    call(&mut s, "insert", json!({"at": 10, "residues": "GG"}));
    call(
        &mut s,
        "add_feature",
        json!({"feature": {"name": "site", "segments": [{"range": {"start": 11, "end": 2}}]}}),
    );
    let doc = call(&mut s, "doc_get", json!({}))["structuredContent"].clone();
    assert_eq!(doc["molecule"]["length"], 12);
    assert_eq!(doc["molecule"]["features"][0]["name"], "site");
    assert_eq!(doc["can_undo"], true);

    // A failed edit is a tool error with a structured payload, not a protocol error.
    let r = call(&mut s, "insert", json!({"at": 99, "residues": "A"}));
    assert_eq!(r["isError"], true);
    assert_eq!(r["structuredContent"]["error"]["code"], "invalid_input");
    let log = call(&mut s, "log_tail", json!({"n": 5}))["structuredContent"].clone();
    assert_eq!(log["entries"][0]["source"], "insert");
}

#[test]
fn protocol_errors() {
    let mut s = server();
    let err = |out: Option<String>| -> i64 {
        let v: Value = serde_json::from_str(&out.unwrap()).unwrap();
        v["error"]["code"].as_i64().unwrap()
    };
    assert_eq!(err(s.handle_line("{not json")), -32700);
    assert_eq!(err(s.handle_line("[1,2]")), -32600);
    assert_eq!(err(s.handle_line(r#"{"jsonrpc":"2.0","id":1}"#)), -32600);
    assert_eq!(
        err(s.handle_line(r#"{"jsonrpc":"2.0","id":1,"method":"nope"}"#)),
        -32601
    );
    assert_eq!(
        err(s.handle_line(
            r#"{"jsonrpc":"2.0","id":1,"method":"tools/call","params":{"name":"nope"}}"#
        )),
        -32602
    );
    assert!(s.handle_line("   ").is_none());
}

#[test]
fn serve_handles_a_stream_with_garbage_and_overlong_lines() {
    let mut s = server();
    let long = "x".repeat(snappiesnap_mcp::MAX_LINE + 10);
    let input = format!(
        "{}\n\u{0}garbage\n{long}\n{}\n",
        json!({"jsonrpc":"2.0","id":1,"method":"ping"}),
        json!({"jsonrpc":"2.0","id":2,"method":"tools/call","params":{"name":"app_version"}}),
    );
    let mut bytes = input.into_bytes();
    bytes.extend_from_slice(b"\xff\xfe\n");
    let mut out = Vec::new();
    serve(&mut s, bytes.as_slice(), &mut out).unwrap();
    let lines: Vec<Value> = String::from_utf8(out)
        .unwrap()
        .lines()
        .map(|l| serde_json::from_str(l).unwrap())
        .collect();
    assert_eq!(lines.len(), 5);
    assert_eq!(lines[0]["id"], 1);
    assert_eq!(lines[1]["error"]["code"], -32700);
    assert!(
        lines[2]["error"]["message"]
            .as_str()
            .unwrap()
            .contains("longer than")
    );
    assert_eq!(
        lines[3]["result"]["structuredContent"]["name"],
        "SnappieSnap"
    );
    assert_eq!(lines[4]["error"]["code"], -32700);
}

proptest! {
    /// Arbitrary input lines never panic, and every response is valid JSON-RPC.
    #[test]
    fn arbitrary_lines_never_panic(line in ".{0,200}") {
        let mut s = server();
        if let Some(out) = s.handle_line(&line) {
            let v: Value = serde_json::from_str(&out).unwrap();
            prop_assert_eq!(&v["jsonrpc"], "2.0");
        }
    }

    /// Tool calls with arbitrary JSON arguments never panic.
    #[test]
    fn arbitrary_tool_arguments_never_panic(args in arb_json(), tool in 0usize..40) {
        let mut s = server();
        call(&mut s, "new", json!({"residues": "ACGTACGT", "topology": "circular"}));
        let tools = s.tools();
        let name = tools[tool % tools.len()]["name"].as_str().unwrap().to_string();
        let r = call(&mut s, &name, args);
        prop_assert!(r["isError"].is_boolean());
        if r["isError"] == true {
            prop_assert_ne!(&r["structuredContent"]["error"]["code"], "internal");
        }
    }
}

fn arb_json() -> impl Strategy<Value = Value> {
    let leaf = prop_oneof![
        Just(Value::Null),
        any::<bool>().prop_map(Value::from),
        any::<i64>().prop_map(Value::from),
        any::<f64>()
            .prop_map(|f| serde_json::Number::from_f64(f).map_or(Value::Null, Value::Number)),
        ".{0,12}".prop_map(Value::from),
        Just(json!({"start": 6, "end": 2})),
    ];
    leaf.prop_recursive(3, 24, 6, |inner| {
        prop_oneof![
            proptest::collection::vec(inner.clone(), 0..4).prop_map(Value::Array),
            proptest::collection::btree_map(
                prop_oneof![
                    Just("at".to_string()),
                    Just("range".to_string()),
                    Just("residues".to_string()),
                    Just("name".to_string()),
                    Just("feature".to_string()),
                    Just("id".to_string()),
                    Just("origin".to_string()),
                    Just("topology".to_string()),
                    ".{0,6}".prop_map(|s| s),
                ],
                inner,
                0..4
            )
            .prop_map(|m| Value::Object(m.into_iter().collect())),
        ]
    })
}
