//! SnappieSnap `mcp` (layer L5).
//!
//! A Model Context Protocol server over newline-delimited JSON-RPC 2.0 (the MCP stdio
//! transport). Tools are generated from the engine's command registry, so every command is
//! available to agents without extra wiring; the protocol methods (`doc.get`, `log.tail`, ...)
//! are tools too, with `.` replaced by `_`.
//!
//! The server never panics and never writes anything but protocol messages to its output:
//! malformed input gets a JSON-RPC error, and tool failures are returned as tool results with
//! `isError: true` and a structured error.
//!
//! Layering is enforced by `cargo xtask layers`; see `AGENTS.md`.

use std::io::{self, BufRead, Write};

#[cfg(not(target_arch = "wasm32"))]
mod attach;
#[cfg(target_arch = "wasm32")]
mod attach {
    pub fn ui_tools() -> Vec<serde_json::Value> {
        Vec::new()
    }
}

use serde_json::{Map, Value, json};
use snappiesnap_engine::{CmdError, Session, api};

/// MCP protocol revisions we speak, newest first.
pub const PROTOCOL_VERSIONS: &[&str] = &["2025-06-18", "2025-03-26", "2024-11-05"];

/// Longest request line we accept; longer lines are discarded with an error.
pub const MAX_LINE: usize = 16 * 1024 * 1024;

// JSON-RPC error codes.
const PARSE_ERROR: i64 = -32700;
const INVALID_REQUEST: i64 = -32600;
const METHOD_NOT_FOUND: i64 = -32601;
const INVALID_PARAMS: i64 = -32602;

const INSTRUCTIONS: &str = "SnappieSnap is a molecular biology sequence editor. Every action is a \
tool. Start with `new` (or check `doc_get`), then edit with tools such as `insert`, `delete`, \
`add_feature`. Positions are 0-based; ranges are half-open {start, end} and wrap past the origin \
of circular molecules when start > end. Failed tools change nothing and return a structured \
error with a hint. `doc_get` shows the whole document; `log_tail` shows recent problems.";

const ATTACHED: &str = " This server is attached to the running SnappieSnap app: tools act on \
its open document, and `ui_state`, `ui_set`, `ui_render` (an offscreen image of the window) and \
`ui_screenshot` show and drive the window.";

/// An MCP server around one [`Session`], or attached to the running app's control channel.
pub struct Server {
    session: Session,
    #[cfg(not(target_arch = "wasm32"))]
    remote: Option<attach::Remote>,
}

impl Server {
    pub fn new(session: Session) -> Self {
        Server {
            session,
            #[cfg(not(target_arch = "wasm32"))]
            remote: None,
        }
    }

    /// A server whose tools act on the running app with `--control PORT` (on 127.0.0.1).
    /// `session` provides the tool list (the command registry); the app does the work.
    #[cfg(not(target_arch = "wasm32"))]
    pub fn attach(session: Session, port: u16) -> io::Result<Self> {
        Ok(Server {
            session,
            remote: Some(attach::Remote::connect(port)?),
        })
    }

    fn attached(&self) -> bool {
        #[cfg(not(target_arch = "wasm32"))]
        {
            self.remote.is_some()
        }
        #[cfg(target_arch = "wasm32")]
        {
            false
        }
    }

    pub fn session(&self) -> &Session {
        &self.session
    }

    /// Every tool: one per command, then one per protocol method (except `engine.execute`,
    /// which the per-command tools replace).
    pub fn tools(&self) -> Vec<Value> {
        let commands = self.session.registry().iter().map(|c| {
            let mut description = c.description.to_string();
            if let Some(menu) = c.menu {
                description.push_str(&format!(" (menu: {menu})"));
            }
            json!({
                "name": c.id,
                "title": c.label,
                "description": description,
                "inputSchema": c.params_schema(),
            })
        });
        let methods = api::METHODS
            .iter()
            .filter(|m| m.name != "engine.execute")
            .map(|m| {
                json!({
                    "name": tool_name(m.name),
                    "title": m.name,
                    "description": m.description,
                    "inputSchema": (m.params_schema)(),
                })
            });
        let mut tools: Vec<Value> = commands.chain(methods).collect();
        if self.attached() {
            tools.extend(attach::ui_tools());
        }
        tools
    }

    /// Handle one request line. Returns the response line, or `None` for notifications.
    pub fn handle_line(&mut self, line: &str) -> Option<String> {
        let line = line.trim();
        if line.is_empty() {
            return None;
        }
        let msg: Value = match serde_json::from_str(line) {
            Ok(v) => v,
            Err(e) => {
                return Some(error_response(
                    Value::Null,
                    PARSE_ERROR,
                    &format!("parse error: {e}"),
                ));
            }
        };
        let Some(obj) = msg.as_object() else {
            return Some(error_response(
                Value::Null,
                INVALID_REQUEST,
                "expected a JSON-RPC request object (batches are not supported)",
            ));
        };
        let id = obj.get("id").cloned();
        let Some(method) = obj.get("method").and_then(Value::as_str) else {
            // A response from the client (we send no requests) or a malformed message.
            return id.map(|id| error_response(id, INVALID_REQUEST, "missing \"method\""));
        };
        let params = obj.get("params").cloned().unwrap_or(Value::Null);
        let result = self.dispatch(method, params);
        // Notifications (no id) get no response, whatever happened.
        let id = id?;
        Some(match result {
            Ok(value) => json!({ "jsonrpc": "2.0", "id": id, "result": value }).to_string(),
            Err((code, message)) => error_response(id, code, &message),
        })
    }

    fn dispatch(&mut self, method: &str, params: Value) -> Result<Value, (i64, String)> {
        match method {
            "initialize" => Ok(initialize(&params, self.attached())),
            "ping" => Ok(json!({})),
            "tools/list" => Ok(json!({ "tools": self.tools() })),
            "tools/call" => self.call_tool(&params),
            m if m.starts_with("notifications/") => Ok(Value::Null),
            other => Err((METHOD_NOT_FOUND, format!("method not found: {other}"))),
        }
    }

    fn call_tool(&mut self, params: &Value) -> Result<Value, (i64, String)> {
        let name = params
            .get("name")
            .and_then(Value::as_str)
            .ok_or((INVALID_PARAMS, "tools/call needs a \"name\"".to_string()))?;
        let arguments = params.get("arguments").cloned().unwrap_or(Value::Null);
        #[cfg(not(target_arch = "wasm32"))]
        if let Some(remote) = &mut self.remote {
            let (method, params) = if self.session.registry().get(name).is_some() {
                (
                    "engine.execute".to_string(),
                    json!({ "command": name, "params": arguments }),
                )
            } else if let Some(m) = api::METHODS
                .iter()
                .find(|m| m.name != "engine.execute" && tool_name(m.name) == name)
            {
                (m.name.to_string(), arguments)
            } else if let Some(m) = attach::UI_METHODS.iter().find(|m| tool_name(m) == name) {
                ((*m).to_string(), arguments)
            } else {
                return Err((INVALID_PARAMS, format!("unknown tool: {name}")));
            };
            return Ok(match remote.call(&method, params) {
                Ok(value) => tool_result(Ok(value)),
                Err(error) => error_result(error),
            });
        }
        let outcome = if self.session.registry().get(name).is_some() {
            self.session.execute(name, arguments)
        } else if let Some(m) = api::METHODS
            .iter()
            .find(|m| m.name != "engine.execute" && tool_name(m.name) == name)
        {
            api::call(&mut self.session, m.name, arguments)
        } else {
            return Err((INVALID_PARAMS, format!("unknown tool: {name}")));
        };
        Ok(tool_result(outcome))
    }
}

/// `doc.get` → `doc_get`: MCP tool names are safest without dots.
fn tool_name(method: &str) -> String {
    method.replace('.', "_")
}

fn initialize(params: &Value, attached: bool) -> Value {
    let requested = params.get("protocolVersion").and_then(Value::as_str);
    let version = requested
        .filter(|v| PROTOCOL_VERSIONS.contains(v))
        .or(PROTOCOL_VERSIONS.first().copied())
        .unwrap_or_default();
    json!({
        "protocolVersion": version,
        "capabilities": { "tools": { "listChanged": false } },
        "serverInfo": {
            "name": "snappiesnap",
            "title": "SnappieSnap",
            "version": env!("CARGO_PKG_VERSION"),
        },
        "instructions": if attached {
            format!("{INSTRUCTIONS}{ATTACHED}")
        } else {
            INSTRUCTIONS.to_string()
        },
    })
}

fn tool_result(outcome: Result<Value, CmdError>) -> Value {
    match outcome {
        Ok(value) => {
            let structured = if value.is_object() {
                value
            } else {
                let mut m = Map::new();
                m.insert("value".into(), value);
                Value::Object(m)
            };
            json!({
                "content": [{ "type": "text", "text": structured.to_string() }],
                "structuredContent": structured,
                "isError": false,
            })
        }
        Err(e) => {
            let error = json!({ "error": e });
            json!({
                "content": [{ "type": "text", "text": error.to_string() }],
                "structuredContent": error,
                "isError": true,
            })
        }
    }
}

/// A failed tool call with a structured error (`{code, message, hint?}`). Used by `--attach`,
/// which the web build doesn't have.
#[cfg_attr(target_arch = "wasm32", allow(dead_code))]
fn error_result(error: Value) -> Value {
    let error = json!({ "error": error });
    json!({
        "content": [{ "type": "text", "text": error.to_string() }],
        "structuredContent": error,
        "isError": true,
    })
}

fn error_response(id: Value, code: i64, message: &str) -> String {
    json!({
        "jsonrpc": "2.0",
        "id": id,
        "error": { "code": code, "message": message },
    })
    .to_string()
}

/// Serve requests from `input` until it closes, writing one response line per request.
pub fn serve(
    server: &mut Server,
    mut input: impl BufRead,
    mut output: impl Write,
) -> io::Result<()> {
    let mut buf = Vec::new();
    loop {
        buf.clear();
        let response = match read_line_bounded(&mut input, &mut buf)? {
            Line::Eof => return Ok(()),
            Line::TooLong => Some(error_response(
                Value::Null,
                PARSE_ERROR,
                &format!("request longer than {MAX_LINE} bytes"),
            )),
            Line::Complete => match std::str::from_utf8(&buf) {
                Ok(text) => server.handle_line(text),
                Err(_) => Some(error_response(
                    Value::Null,
                    PARSE_ERROR,
                    "request is not valid UTF-8",
                )),
            },
        };
        if let Some(line) = response {
            output.write_all(line.as_bytes())?;
            output.write_all(b"\n")?;
            output.flush()?;
        }
    }
}

enum Line {
    Complete,
    TooLong,
    Eof,
}

/// Read up to and including `\n`, keeping at most [`MAX_LINE`] bytes; the rest of an
/// over-long line is skipped.
fn read_line_bounded(input: &mut impl BufRead, buf: &mut Vec<u8>) -> io::Result<Line> {
    let mut too_long = false;
    let mut read_any = false;
    loop {
        let available = match input.fill_buf() {
            Ok(b) => b,
            Err(e) if e.kind() == io::ErrorKind::Interrupted => continue,
            Err(e) => return Err(e),
        };
        if available.is_empty() {
            return Ok(match (read_any, too_long) {
                (false, _) => Line::Eof,
                (true, true) => Line::TooLong,
                (true, false) => Line::Complete,
            });
        }
        read_any = true;
        let (chunk, found_newline) = match available.iter().position(|&b| b == b'\n') {
            Some(i) => (available.get(..i).unwrap_or_default(), true),
            None => (available, false),
        };
        if !too_long {
            if buf.len().saturating_add(chunk.len()) > MAX_LINE {
                too_long = true;
                buf.clear();
            } else {
                buf.extend_from_slice(chunk);
            }
        }
        let used = chunk.len().saturating_add(usize::from(found_newline));
        input.consume(used);
        if found_newline {
            return Ok(if too_long {
                Line::TooLong
            } else {
                Line::Complete
            });
        }
    }
}
