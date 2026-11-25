//! Attaching to the running app: tool calls become requests on its JSON-lines control channel
//! (`snappiesnap --control PORT`, 127.0.0.1 only; see `docs/control-protocol.md`).

use std::io::{self, BufRead, BufReader, Read, Write};
use std::net::{Ipv4Addr, SocketAddr, TcpStream};
use std::time::Duration;

use serde_json::{Value, json};

/// Longest a request may take (a screenshot waits for the next frame).
const TIMEOUT: Duration = Duration::from_secs(60);
/// Longest response line accepted.
const MAX_RESPONSE: usize = 256 * 1024 * 1024;

/// The window's methods, available only when attached.
pub const UI_METHODS: [&str; 4] = ["ui.state", "ui.set", "ui.render", "ui.screenshot"];

/// Tool definitions for [`UI_METHODS`].
pub fn ui_tools() -> Vec<Value> {
    vec![
        json!({
            "name": "ui_state",
            "title": "ui.state",
            "description": "The window's state: the view, theme, map style and layers, the gel, the enzyme set, map zoom, and what is open.",
            "inputSchema": {"type": "object", "properties": {}, "additionalProperties": false},
        }),
        json!({
            "name": "ui_set",
            "title": "ui.set",
            "description": "Change the window: any of view, theme, map_style, folder, sidebar, the layers (features, enzymes, primers, positions, ruler, orfs), translation, genetic_code, enzyme_set, map_zoom, compare, gel, history. Returns the new state (see docs/control-protocol.md).",
            "inputSchema": {"type": "object", "additionalProperties": true},
        }),
        json!({
            "name": "ui_render",
            "title": "ui.render",
            "description": "Draw the window offscreen to a PNG (works with the screen locked or the window hidden).",
            "inputSchema": {
                "type": "object",
                "properties": {
                    "path": {"type": "string"},
                    "width": {"type": "integer", "minimum": 1},
                    "height": {"type": "integer", "minimum": 1},
                    "scale": {"type": "number"},
                },
                "required": ["path"],
            },
        }),
        json!({
            "name": "ui_screenshot",
            "title": "ui.screenshot",
            "description": "A PNG of the window as it is on screen (needs the window visible; ui_render doesn't).",
            "inputSchema": {
                "type": "object",
                "properties": {"path": {"type": "string"}},
                "required": ["path"],
            },
        }),
    ]
}

/// A connection to the app's control channel.
pub struct Remote {
    reader: BufReader<TcpStream>,
    writer: TcpStream,
    next_id: u64,
}

impl Remote {
    pub fn connect(port: u16) -> io::Result<Remote> {
        let address = SocketAddr::from((Ipv4Addr::LOCALHOST, port));
        let stream = TcpStream::connect_timeout(&address, Duration::from_secs(5)).map_err(|e| {
            io::Error::new(
                e.kind(),
                format!(
                    "cannot reach SnappieSnap on 127.0.0.1:{port} ({e}); start it with `snappiesnap --control {port}`"
                ),
            )
        })?;
        stream.set_read_timeout(Some(TIMEOUT))?;
        stream.set_write_timeout(Some(TIMEOUT))?;
        Ok(Remote {
            reader: BufReader::new(stream.try_clone()?),
            writer: stream,
            next_id: 1,
        })
    }

    /// Send one request and wait for its answer: the result, or a structured error
    /// (`{code, message, hint?}`).
    pub fn call(&mut self, method: &str, params: Value) -> Result<Value, Value> {
        let id = self.next_id;
        self.next_id = self.next_id.wrapping_add(1);
        let request = json!({ "id": id, "method": method, "params": params });
        let lost = |e: &dyn std::fmt::Display| {
            json!({
                "code": "io",
                "message": format!("lost the connection to SnappieSnap: {e}"),
                "hint": "is the app still running with --control?",
            })
        };
        self.writer
            .write_all(format!("{request}\n").as_bytes())
            .and_then(|()| self.writer.flush())
            .map_err(|e| lost(&e))?;
        // Answers come in order; skip anything that isn't ours.
        loop {
            let mut line = Vec::new();
            let n = Read::by_ref(&mut self.reader)
                .take(MAX_RESPONSE as u64)
                .read_until(b'\n', &mut line)
                .map_err(|e| lost(&e))?;
            if n == 0 {
                return Err(lost(&"the app closed the connection"));
            }
            let Ok(answer) = serde_json::from_slice::<Value>(&line) else {
                return Err(lost(&"an answer that isn't JSON"));
            };
            if answer.get("id").and_then(Value::as_u64) != Some(id) {
                continue;
            }
            if let Some(error) = answer.get("error") {
                return Err(error.clone());
            }
            return Ok(answer.get("result").cloned().unwrap_or(Value::Null));
        }
    }
}

#[cfg(test)]
mod tests {
    use std::io::{BufRead, BufReader, Write};
    use std::net::TcpListener;

    use serde_json::{Value, json};
    use snappiesnap_engine::Session;

    use crate::Server;

    /// A stand-in for the app's control channel: echoes engine.execute, refuses ui.state, then
    /// closes after `answers` requests.
    fn fake_app(answers: usize) -> (u16, std::thread::JoinHandle<Vec<Value>>) {
        let listener = TcpListener::bind("127.0.0.1:0").unwrap();
        let port = listener.local_addr().unwrap().port();
        let handle = std::thread::spawn(move || {
            let (stream, _) = listener.accept().unwrap();
            let mut reader = BufReader::new(stream.try_clone().unwrap());
            let mut writer = stream;
            let mut seen = Vec::new();
            for _ in 0..answers {
                let mut line = String::new();
                if reader.read_line(&mut line).unwrap() == 0 {
                    break;
                }
                let req: Value = serde_json::from_str(&line).unwrap();
                let id = req["id"].clone();
                let answer = match req["method"].as_str() {
                    Some("ui.state") => {
                        json!({"id": id, "error": {"code": "not_enabled", "message": "no window"}})
                    }
                    _ => json!({"id": id, "result": {"echo": req.clone()}}),
                };
                // Something unrelated first: the client must skip it.
                writeln!(writer, "{}", json!({"id": 999_999, "result": {}})).unwrap();
                writeln!(writer, "{answer}").unwrap();
                seen.push(req);
            }
            seen
        });
        (port, handle)
    }

    fn call(server: &mut Server, name: &str, args: Value) -> Value {
        let line = json!({"jsonrpc": "2.0", "id": 1, "method": "tools/call", "params": {"name": name, "arguments": args}}).to_string();
        let out: Value = serde_json::from_str(&server.handle_line(&line).unwrap()).unwrap();
        out["result"].clone()
    }

    #[test]
    fn tool_calls_go_to_the_running_app() {
        let (port, app) = fake_app(3);
        let mut server = Server::attach(Session::new(), port).unwrap();
        // The window's tools are listed only when attached.
        let names: Vec<String> = server
            .tools()
            .iter()
            .map(|t| t["name"].as_str().unwrap().to_string())
            .collect();
        assert!(names.contains(&"ui_render".to_string()));
        assert!(names.contains(&"insert".to_string()));
        assert!(
            !Server::new(Session::new())
                .tools()
                .iter()
                .any(|t| t["name"] == "ui_state")
        );

        // A command becomes engine.execute; a protocol method keeps its name.
        let r = call(&mut server, "insert", json!({"at": 0, "residues": "A"}));
        assert_eq!(r["isError"], false);
        assert_eq!(r["structuredContent"]["echo"]["method"], "engine.execute");
        assert_eq!(
            r["structuredContent"]["echo"]["params"]["command"],
            "insert"
        );
        let r = call(&mut server, "doc_get", json!({}));
        assert_eq!(r["structuredContent"]["echo"]["method"], "doc.get");
        // The app's errors come back structured.
        let r = call(&mut server, "ui_state", json!({}));
        assert_eq!(r["isError"], true);
        assert_eq!(r["structuredContent"]["error"]["code"], "not_enabled");
        let seen = app.join().unwrap();
        assert_eq!(seen.len(), 3);
        // The app went away: an error, not a hang or a crash.
        let r = call(&mut server, "doc_get", json!({}));
        assert_eq!(r["isError"], true);
        assert_eq!(r["structuredContent"]["error"]["code"], "io");
        // Unknown tools are still refused.
        let line =
            json!({"jsonrpc": "2.0", "id": 2, "method": "tools/call", "params": {"name": "nope"}})
                .to_string();
        assert!(server.handle_line(&line).unwrap().contains("unknown tool"));
    }

    #[test]
    fn nothing_listening_is_a_clear_error() {
        let port = {
            let l = TcpListener::bind("127.0.0.1:0").unwrap();
            l.local_addr().unwrap().port()
        };
        let e = Server::attach(Session::new(), port).err().unwrap();
        assert!(e.to_string().contains("--control"), "{e}");
    }
}
