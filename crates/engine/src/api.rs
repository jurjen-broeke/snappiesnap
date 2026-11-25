//! The control protocol: methods shared by MCP (`snappiesnap-cli mcp`) and the GUI control
//! channel (`--control PORT`). See `docs/control-protocol.md`.
//!
//! [`call`] never panics: `engine.execute` goes through the guarded [`Session::execute`], and
//! every other method runs under its own `catch_unwind`.

use std::panic::{AssertUnwindSafe, catch_unwind};

use schemars::JsonSchema;
use serde::Deserialize;
use serde::de::DeserializeOwned;
use serde_json::{Value, json};

use crate::{CmdError, ErrorCode, NoParams, Session};

/// Version of the method set below; bumped on breaking changes.
pub const PROTOCOL_VERSION: u32 = 1;

/// A protocol method: name, description and params schema.
pub struct Method {
    pub name: &'static str,
    pub description: &'static str,
    pub params_schema: fn() -> Value,
}

/// Every protocol method except the per-command tools.
pub const METHODS: &[Method] = &[
    Method {
        name: "engine.execute",
        description: "Run a command by id with JSON params (see commands.list).",
        params_schema: schema::<ExecuteParams>,
    },
    Method {
        name: "commands.list",
        description: "List every command: id, label, menu path, shortcut, effect.",
        params_schema: schema::<NoParams>,
    },
    Method {
        name: "commands.describe",
        description: "Describe one command, including the JSON Schema of its params.",
        params_schema: schema::<DescribeParams>,
    },
    Method {
        name: "doc.get",
        description: "The open document as JSON: molecule (sequence, features, primers, notes), path, unsaved changes, selection, undo/redo availability.",
        params_schema: schema::<DocGetParams>,
    },
    Method {
        name: "log.tail",
        description: "The most recent warnings and errors, oldest first.",
        params_schema: schema::<LogTailParams>,
    },
    Method {
        name: "app.version",
        description: "Application name, version and protocol version.",
        params_schema: schema::<NoParams>,
    },
];

fn schema<P: JsonSchema>() -> Value {
    serde_json::to_value(schemars::schema_for!(P)).unwrap_or(Value::Null)
}

#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
struct ExecuteParams {
    /// Command id, for example `insert`.
    command: String,
    /// The command's params object.
    #[serde(default)]
    params: Value,
}

#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
struct DescribeParams {
    /// Command id.
    id: String,
}

#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
struct DocGetParams {
    /// Include the residues (set false for large molecules when only annotations matter).
    #[serde(default = "yes")]
    residues: bool,
}

fn yes() -> bool {
    true
}

#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
struct LogTailParams {
    /// How many entries (at most 1000).
    #[serde(default = "default_tail")]
    n: usize,
}

fn default_tail() -> usize {
    50
}

/// Handle one protocol call.
pub fn call(session: &mut Session, method: &str, params: Value) -> Result<Value, CmdError> {
    if method == "engine.execute" {
        let p: ExecuteParams = parse(method, params)?;
        return session.execute(&p.command, p.params);
    }
    let result = catch_unwind(AssertUnwindSafe(|| dispatch(session, method, params)));
    result.unwrap_or_else(|_| {
        let err = CmdError::new(
            ErrorCode::Internal,
            format!("internal error in protocol method {method:?}"),
        )
        .with_hint("this is a bug in SnappieSnap; please report it");
        session.log(crate::Level::Error, method, err.message.clone());
        Err(err)
    })
}

fn dispatch(session: &mut Session, method: &str, params: Value) -> Result<Value, CmdError> {
    match method {
        "commands.list" => {
            let _: NoParams = parse(method, params)?;
            let list: Vec<_> = session.registry().iter().map(|c| c.info(false)).collect();
            Ok(json!({ "commands": list }))
        }
        "commands.describe" => {
            let p: DescribeParams = parse(method, params)?;
            let spec = session.registry().get(&p.id).ok_or_else(|| {
                CmdError::new(
                    ErrorCode::UnknownCommand,
                    format!("unknown command {:?}", p.id),
                )
                .with_hint("`commands.list` lists every command")
            })?;
            Ok(json!({ "command": spec.info(true) }))
        }
        "doc.get" => {
            let p: DocGetParams = parse(method, params)?;
            Ok(doc_get(session, p.residues))
        }
        "log.tail" => {
            let p: LogTailParams = parse(method, params)?;
            Ok(json!({ "entries": session.log_tail(p.n.min(1000)) }))
        }
        "app.version" => {
            let _: NoParams = parse(method, params)?;
            Ok(json!({
                "name": "SnappieSnap",
                "version": env!("CARGO_PKG_VERSION"),
                "protocol": PROTOCOL_VERSION,
            }))
        }
        _ => Err(CmdError::new(
            ErrorCode::UnknownMethod,
            format!("unknown method {method:?}"),
        )
        .with_hint(format!(
            "methods: {}",
            METHODS
                .iter()
                .map(|m| m.name)
                .collect::<Vec<_>>()
                .join(", ")
        ))),
    }
}

fn doc_get(session: &Session, residues: bool) -> Value {
    let Some(doc) = session.document() else {
        return json!({ "open": false });
    };
    let mut molecule = serde_json::to_value(&doc.molecule).unwrap_or(Value::Null);
    if !residues && let Some(seq) = molecule.get_mut("seq").and_then(Value::as_object_mut) {
        seq.remove("residues");
    }
    // The construction history holds copies of every ancestor: summarised here, opened with
    // `history_open`.
    // The trace's signal is tens of thousands of numbers: summarised here, read with
    // `trace_info`.
    if let Some(m) = molecule.as_object_mut() {
        m.remove("history");
        m.remove("trace");
    }
    let trace = doc.molecule.trace().map(|t| {
        json!({
            "sample_name": t.sample_name(),
            "samples": t.samples(),
            "reversed": t.is_reversed(),
        })
    });
    let construction = doc.molecule.history().map(|h| h.steps());
    json!({
        "open": true,
        "path": doc.path,
        "dirty": doc.dirty,
        "selection": doc.selection,
        "can_undo": session.can_undo(),
        "can_redo": session.can_redo(),
        "history": session.history(),
        "construction": construction,
        "trace": trace,
        "molecule": molecule,
    })
}

fn parse<P: DeserializeOwned>(method: &str, params: Value) -> Result<P, CmdError> {
    let params = if params.is_null() {
        Value::Object(serde_json::Map::new())
    } else {
        params
    };
    serde_json::from_value(params).map_err(|e| CmdError::invalid_params(format!("{method}: {e}")))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn methods_work_and_fail_cleanly() {
        let mut s = Session::new();
        assert_eq!(call(&mut s, "doc.get", Value::Null).unwrap()["open"], false);
        call(
            &mut s,
            "engine.execute",
            json!({"command": "new", "params": {"residues": "ACGT"}}),
        )
        .unwrap();
        let doc = call(&mut s, "doc.get", json!({"residues": false})).unwrap();
        assert_eq!(doc["molecule"]["length"], 4);
        assert!(doc["molecule"]["seq"].get("residues").is_none());

        let list = call(&mut s, "commands.list", json!({})).unwrap();
        assert!(list["commands"].as_array().unwrap().len() >= 15);
        let d = call(&mut s, "commands.describe", json!({"id": "insert"})).unwrap();
        assert_eq!(d["command"]["params"]["type"], "object");

        for (m, p) in [
            ("nope", json!({})),
            ("commands.describe", json!({"id": "nope"})),
            ("commands.list", json!({"x": 1})),
            ("engine.execute", json!({})),
            ("log.tail", json!({"n": -1})),
        ] {
            assert!(call(&mut s, m, p).is_err(), "{m}");
        }
        assert_eq!(
            call(&mut s, "app.version", json!({})).unwrap()["protocol"],
            1
        );
    }

    #[test]
    fn every_method_has_an_object_schema() {
        for m in METHODS {
            assert_eq!((m.params_schema)()["type"], "object", "{}", m.name);
        }
    }
}
