//! The JSON control channel (`--control PORT`), native only.
//!
//! JSON lines on `127.0.0.1:PORT`: each request `{"id": .., "method": "..", "params": {..}}`
//! gets one response `{"id": .., "result": ..}` or `{"id": .., "error": {code, message, hint}}`.
//! Connections are served on background threads; every request is handed to the UI thread, which
//! owns the session, so the app behaves exactly as if a user had clicked. Methods are listed in
//! `docs/control-protocol.md`.

use std::io::{self, BufRead, BufReader, Write};
use std::net::{TcpListener, TcpStream};
use std::sync::mpsc::{self, Receiver, Sender};
use std::time::Duration;

use eframe::egui;
use serde_json::{Value, json};
use snappiesnap_engine::{CmdError, ErrorCode};

/// Longest request line accepted.
const MAX_LINE: usize = 16 * 1024 * 1024;
/// How long a connection waits for the UI to answer.
const REPLY_TIMEOUT: Duration = Duration::from_secs(60);

/// A request waiting for the UI thread.
pub struct Request {
    pub method: String,
    pub params: Value,
    pub reply: Sender<Result<Value, CmdError>>,
}

pub struct ControlServer {
    pub port: u16,
    requests: Receiver<Request>,
}

impl ControlServer {
    /// Listen on `127.0.0.1:port` (0 picks a free port). Only the local machine can connect.
    pub fn start(port: u16, ctx: egui::Context) -> io::Result<ControlServer> {
        let listener = TcpListener::bind(("127.0.0.1", port))?;
        let port = listener.local_addr()?.port();
        let (tx, rx) = mpsc::channel();
        std::thread::Builder::new()
            .name("control-accept".into())
            .spawn(move || {
                for stream in listener.incoming().flatten() {
                    let tx = tx.clone();
                    let ctx = ctx.clone();
                    let _ = std::thread::Builder::new()
                        .name("control-conn".into())
                        .spawn(move || serve(stream, &tx, &ctx));
                }
            })?;
        Ok(ControlServer { port, requests: rx })
    }

    /// The next waiting request, if any.
    pub fn try_recv(&self) -> Option<Request> {
        self.requests.try_recv().ok()
    }
}

fn serve(stream: TcpStream, tx: &Sender<Request>, ctx: &egui::Context) {
    let Ok(write_half) = stream.try_clone() else {
        return;
    };
    let mut out = io::BufWriter::new(write_half);
    let mut input = BufReader::new(stream);
    let mut buf = Vec::new();
    loop {
        buf.clear();
        let response = match read_line(&mut input, &mut buf) {
            Ok(Line::Eof) | Err(_) => return,
            Ok(Line::TooLong) => error_line(
                &Value::Null,
                &CmdError::invalid_params(format!("request longer than {MAX_LINE} bytes")),
            ),
            Ok(Line::Complete) => match handle(&buf, tx, ctx) {
                Some(line) => line,
                None => continue,
            },
        };
        if out
            .write_all(response.as_bytes())
            .and_then(|()| out.write_all(b"\n"))
            .and_then(|()| out.flush())
            .is_err()
        {
            return;
        }
    }
}

/// One request line → one response line (`None` for blank lines).
fn handle(line: &[u8], tx: &Sender<Request>, ctx: &egui::Context) -> Option<String> {
    let text = String::from_utf8_lossy(line);
    let text = text.trim();
    if text.is_empty() {
        return None;
    }
    let request: Value = match serde_json::from_str(text) {
        Ok(v) => v,
        Err(e) => {
            return Some(error_line(
                &Value::Null,
                &CmdError::invalid_params(format!("not JSON: {e}")),
            ));
        }
    };
    let id = request.get("id").cloned().unwrap_or(Value::Null);
    let Some(method) = request.get("method").and_then(Value::as_str) else {
        return Some(error_line(
            &id,
            &CmdError::invalid_params("a request needs a \"method\""),
        ));
    };
    let (reply_tx, reply_rx) = mpsc::channel();
    let sent = tx.send(Request {
        method: method.to_string(),
        params: request.get("params").cloned().unwrap_or(Value::Null),
        reply: reply_tx,
    });
    if sent.is_err() {
        return Some(error_line(
            &id,
            &CmdError::new(ErrorCode::Internal, "the app is shutting down"),
        ));
    }
    ctx.request_repaint();
    Some(match reply_rx.recv_timeout(REPLY_TIMEOUT) {
        Ok(Ok(result)) => json!({ "id": id, "result": result }).to_string(),
        Ok(Err(e)) => error_line(&id, &e),
        Err(_) => error_line(
            &id,
            &CmdError::new(ErrorCode::Internal, "the app didn't answer in time")
                .with_hint("is the window minimised or the screen locked? try ui.render"),
        ),
    })
}

fn error_line(id: &Value, e: &CmdError) -> String {
    json!({ "id": id, "error": e }).to_string()
}

enum Line {
    Complete,
    TooLong,
    Eof,
}

fn read_line(input: &mut impl BufRead, buf: &mut Vec<u8>) -> io::Result<Line> {
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
        let (chunk, newline) = match available.iter().position(|&b| b == b'\n') {
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
        let used = chunk.len().saturating_add(usize::from(newline));
        input.consume(used);
        if newline {
            return Ok(if too_long {
                Line::TooLong
            } else {
                Line::Complete
            });
        }
    }
}

/// Write an egui screenshot or render as a PNG file.
pub fn save_png(path: &str, image: &egui::ColorImage) -> Result<Value, CmdError> {
    crate::offscreen::write_png(path, image)
        .map_err(|e| CmdError::new(ErrorCode::Io, format!("cannot write {path}: {e}")))?;
    let [width, height] = image.size;
    Ok(json!({ "path": path, "width": width, "height": height }))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn requests_reach_the_ui_and_answers_come_back() {
        let ctx = egui::Context::default();
        let server = ControlServer::start(0, ctx).unwrap();
        let port = server.port;
        let client = std::thread::spawn(move || {
            let mut stream = TcpStream::connect(("127.0.0.1", port)).unwrap();
            stream
                .write_all(
                    b"{\"id\":1,\"method\":\"app.version\"}\nnot json\n{\"id\":2}\n\n{\"id\":3,\"method\":\"x\"}\n",
                )
                .unwrap();
            let mut lines = BufReader::new(stream).lines();
            (0..4)
                .map(|_| serde_json::from_str::<Value>(&lines.next().unwrap().unwrap()).unwrap())
                .collect::<Vec<_>>()
        });
        // Play the UI thread: answer the two requests that have a method.
        let mut answered = 0;
        while answered < 2 {
            if let Some(req) = server.try_recv() {
                let result = if req.method == "app.version" {
                    Ok(json!({"name": "test"}))
                } else {
                    Err(CmdError::new(ErrorCode::UnknownMethod, "nope"))
                };
                req.reply.send(result).unwrap();
                answered += 1;
            }
            std::thread::sleep(Duration::from_millis(5));
        }
        let lines = client.join().unwrap();
        assert_eq!(lines[0], json!({"id": 1, "result": {"name": "test"}}));
        assert_eq!(lines[1]["error"]["code"], "invalid_params");
        assert_eq!(lines[2]["id"], 2);
        assert_eq!(lines[3]["error"]["code"], "unknown_method");
    }
}
