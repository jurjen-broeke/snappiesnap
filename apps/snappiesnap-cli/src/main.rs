//! `snappiesnap-cli`: headless SnappieSnap.
//!
//! Every failure is a message on stderr and a non-zero exit code; nothing panics. As a last
//! resort the whole run is wrapped in `catch_unwind` so even a bug ends in an error message.

use std::io::{self, BufReader, Read, Write};
use std::panic::{AssertUnwindSafe, catch_unwind};
use std::process::ExitCode;

use serde_json::{Value, json};
use snappiesnap_engine::io::{self as sio, Source};
use snappiesnap_engine::{NativeFs, Platform, Session, api, script};
use snappiesnap_mcp::Server;
use snappiesnap_ui_egui::state::UiPatch;
use snappiesnap_ui_egui::{App, offscreen};

const USAGE: &str = "usage: snappiesnap-cli <command> [args]

commands:
  mcp                          serve the MCP protocol on stdin/stdout (for agents)
      [--attach PORT]          act on the running app started with `snappiesnap --control PORT`
  run --script <file|->        run a JSON-lines command script ('-' reads stdin)
      [--keep-going]           continue after a failing step
      [--quiet]                print only failures and the summary
      [--export out.png]       afterwards, render the document as the app shows it
  render <file> <out.png>      render a file as the app shows it, without a window
                               (out.svg: the map as vector graphics, text kept as text)
  image options (render, run --export):
      [--view map|sequence|features|primers|enzymes|gel|history] [--map auto|circular|linear]
      [--theme light|dark] [--size 1280x820] [--scale 1.0]
  convert <in> <out>           convert between .dna, GenBank (.gb/.gbk), EMBL (.embl), GFF3 (.gff3), SBOL 2 (.sbol/.xml) and FASTA (.fa)
  info <file>                  summarise a file (.dna, GenBank, EMBL, GFF3, SBOL, FASTA) as JSON
  commands                     list every command
  describe <command-id>        show one command with its parameter schema
  --version                    print the version
  help                         print this help


See docs/control-protocol.md for the protocol and script format.";

/// Exit code for usage errors and unavailable features.
const USAGE_ERROR: u8 = 2;

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().skip(1).collect();
    match catch_unwind(AssertUnwindSafe(|| run(&args))) {
        Ok(code) => code,
        Err(_) => {
            eprintln!("error: internal error (this is a bug in SnappieSnap; please report it)");
            ExitCode::FAILURE
        }
    }
}

fn run(args: &[String]) -> ExitCode {
    let mut args = args.iter().map(String::as_str);
    match args.next() {
        Some("--version" | "-V") => {
            println!("snappiesnap-cli {}", env!("CARGO_PKG_VERSION"));
            ExitCode::SUCCESS
        }
        None | Some("help" | "--help" | "-h") => {
            println!("{USAGE}");
            ExitCode::SUCCESS
        }
        Some("mcp") => match (args.next(), args.next(), args.next()) {
            (None, _, _) => mcp(None),
            (Some("--attach"), Some(port), None) => match port.parse::<u16>() {
                Ok(p) if p > 0 => mcp(Some(p)),
                _ => usage_error("--attach needs a port number (1-65535)"),
            },
            _ => usage_error("`mcp` takes only --attach PORT"),
        },
        Some("run") => run_script(args.collect()),
        Some("commands") => protocol("commands.list", json!({})),
        Some("describe") => match args.next() {
            Some(id) => protocol("commands.describe", json!({ "id": id })),
            None => usage_error("`describe` needs a command id"),
        },
        Some("render") => {
            let rest: Vec<&str> = args.collect();
            match rest.as_slice() {
                [file, out, flags @ ..] => render_file(file, out, flags),
                _ => usage_error("`render` needs <file> <out.png>"),
            }
        }
        Some("info") => match args.next() {
            Some(path) => info(path),
            None => usage_error("`info` needs a file"),
        },
        Some("convert") => {
            let rest: Vec<&str> = args.collect();
            match rest.as_slice() {
                [input, output] => convert(input, output),
                _ => usage_error("`convert` needs <in> <out>"),
            }
        }
        Some(other) => usage_error(&format!("unknown command `{other}`")),
    }
}

fn usage_error(message: &str) -> ExitCode {
    eprintln!("error: {message}\n\n{USAGE}");
    ExitCode::from(USAGE_ERROR)
}

/// A session with the standard commands and the cached REBASE enzymes, if downloaded before.
fn session() -> Session {
    let mut s = Session::new();
    s.load_enzyme_cache();
    s.load_enzyme_sets();
    s
}

fn mcp(attach: Option<u16>) -> ExitCode {
    let mut server = match attach {
        None => Server::new(session()),
        Some(port) => match Server::attach(session(), port) {
            Ok(s) => s,
            Err(e) => {
                eprintln!("error: {e}");
                return ExitCode::FAILURE;
            }
        },
    };
    let stdin = io::stdin();
    let stdout = io::stdout();
    match snappiesnap_mcp::serve(&mut server, stdin.lock(), stdout.lock()) {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("error: MCP transport failed: {e}");
            ExitCode::FAILURE
        }
    }
}

/// Run one protocol method and print its result as pretty JSON.
fn protocol(method: &str, params: Value) -> ExitCode {
    let mut session = session();
    match api::call(&mut session, method, params) {
        Ok(v) => {
            println!(
                "{}",
                serde_json::to_string_pretty(&v).unwrap_or_else(|_| v.to_string())
            );
            ExitCode::SUCCESS
        }
        Err(e) => {
            eprintln!("error: {e}");
            ExitCode::FAILURE
        }
    }
}

fn run_script(args: Vec<&str>) -> ExitCode {
    let mut path = None;
    let mut keep_going = false;
    let mut quiet = false;
    let mut export_to = None;
    let mut image_flags = Vec::new();
    let mut it = args.into_iter();
    while let Some(arg) = it.next() {
        match arg {
            "--script" => path = it.next(),
            "--keep-going" => keep_going = true,
            "--quiet" => quiet = true,
            "--export" => export_to = it.next(),
            flag @ ("--view" | "--map" | "--theme" | "--size" | "--scale") => {
                image_flags.push(flag);
                if let Some(v) = it.next() {
                    image_flags.push(v);
                }
            }
            other => return usage_error(&format!("unknown option `{other}` for `run`")),
        }
    }
    let image = match ImageOptions::parse(&image_flags) {
        Ok(o) => o,
        Err(e) => return usage_error(&e),
    };
    let Some(path) = path else {
        return usage_error("`run` needs --script <file|->");
    };
    let text = match read_script(path) {
        Ok(t) => t,
        Err(e) => {
            eprintln!("error: cannot read script {path}: {e}");
            return ExitCode::FAILURE;
        }
    };

    let mut session = session();
    let stdout = io::stdout();
    let mut out = stdout.lock();
    let summary = script::run(&mut session, &text, |step| {
        let ok = step.result.is_ok();
        if !quiet || !ok {
            let line = match &step.result {
                Ok(result) => json!({ "line": step.line, "ok": true, "result": result }),
                Err(error) => json!({ "line": step.line, "ok": false, "error": error }),
            };
            // A closed stdout must not stop us from reporting the summary on stderr.
            let _ = writeln!(out, "{line}");
        }
        ok || keep_going
    });
    eprintln!(
        "{} step(s), {} failed{}",
        summary.steps,
        summary.failed,
        if summary.stopped {
            " (stopped at the first failure; use --keep-going to continue)"
        } else {
            ""
        }
    );
    if let Some(out) = export_to {
        if let Err(e) = export(session, out, &image) {
            eprintln!("error: {e}");
            return ExitCode::FAILURE;
        }
        eprintln!("rendered {out}");
    }
    if summary.failed == 0 {
        ExitCode::SUCCESS
    } else {
        ExitCode::FAILURE
    }
}

/// How to render an image (`render`, `run --export`).
struct ImageOptions {
    patch: serde_json::Map<String, Value>,
    width: u32,
    height: u32,
    scale: f32,
}

impl ImageOptions {
    fn parse(flags: &[&str]) -> Result<Self, String> {
        let mut o = ImageOptions {
            patch: serde_json::Map::from_iter([("theme".to_string(), json!("light"))]),
            width: 1280,
            height: 820,
            scale: 1.0,
        };
        let mut it = flags.iter();
        while let Some(flag) = it.next() {
            let value = it.next().ok_or(format!("{flag} needs a value"))?;
            match *flag {
                "--view" => {
                    o.patch.insert("view".into(), json!(value));
                }
                "--map" => {
                    o.patch.insert("map_style".into(), json!(value));
                }
                "--theme" => {
                    o.patch.insert("theme".into(), json!(value));
                }
                "--size" => {
                    let (w, h) = value
                        .split_once('x')
                        .ok_or("--size takes WIDTHxHEIGHT, e.g. 1280x820")?;
                    o.width = w.parse().map_err(|_| format!("--size: bad width {w:?}"))?;
                    o.height = h.parse().map_err(|_| format!("--size: bad height {h:?}"))?;
                }
                "--scale" => {
                    o.scale = value
                        .parse()
                        .map_err(|_| format!("--scale: bad number {value:?}"))?;
                }
                other => return Err(format!("unknown option `{other}`")),
            }
        }
        // Check the view/map/theme names now, for a clear error.
        serde_json::from_value::<UiPatch>(Value::Object(o.patch.clone()))
            .map_err(|e| e.to_string())?;
        Ok(o)
    }
}

/// Render the session's document as the app would show it.
fn export(session: Session, out: &str, o: &ImageOptions) -> Result<(), String> {
    let mut app = App::headless(session);
    let patch: UiPatch =
        serde_json::from_value(Value::Object(o.patch.clone())).map_err(|e| e.to_string())?;
    app.set_ui(patch);
    // A vector map: the map alone, as SVG.
    if out.to_ascii_lowercase().ends_with(".svg") {
        if o.patch
            .get("view")
            .and_then(Value::as_str)
            .is_some_and(|v| v != "map")
        {
            return Err("SVG output is for the map view (leave out --view, or use .png)".into());
        }
        app.save_map(out, o.width, o.height)?;
        return Ok(());
    }
    let image = app.render(o.width, o.height, o.scale)?;
    offscreen::write_png(out, &image).map_err(|e| format!("cannot write {out}: {e}"))
}

/// Open `input` and save it as `output` (format from the extension), reporting warnings.
fn convert(input: &str, output: &str) -> ExitCode {
    let mut session = Session::new();
    let opened = match session.execute("open", json!({ "path": input })) {
        Ok(v) => v,
        Err(e) => {
            eprintln!("error: {e}");
            return ExitCode::FAILURE;
        }
    };
    match session.execute("save", json!({ "path": output })) {
        Ok(saved) => {
            let warnings = |v: &Value| -> Vec<String> {
                v.get("warnings")
                    .and_then(Value::as_array)
                    .map(|a| {
                        a.iter()
                            .filter_map(Value::as_str)
                            .map(str::to_string)
                            .collect()
                    })
                    .unwrap_or_default()
            };
            for w in warnings(&opened).into_iter().chain(warnings(&saved)) {
                eprintln!("warning: {w}");
            }
            let format = |v: &Value| {
                v.get("format")
                    .and_then(Value::as_str)
                    .unwrap_or("?")
                    .to_string()
            };
            eprintln!(
                "converted {input} ({}) to {output} ({})",
                format(&opened),
                format(&saved)
            );
            ExitCode::SUCCESS
        }
        Err(e) => {
            eprintln!("error: {e}");
            ExitCode::FAILURE
        }
    }
}

fn render_file(file: &str, out: &str, flags: &[&str]) -> ExitCode {
    let image = match ImageOptions::parse(flags) {
        Ok(o) => o,
        Err(e) => return usage_error(&e),
    };
    let mut session = session();
    if let Err(e) = session.execute("open", json!({ "path": file })) {
        eprintln!("error: {e}");
        return ExitCode::FAILURE;
    }
    match export(session, out, &image) {
        Ok(()) => {
            eprintln!("rendered {out}");
            ExitCode::SUCCESS
        }
        Err(e) => {
            eprintln!("error: {e}");
            ExitCode::FAILURE
        }
    }
}

fn read_script(path: &str) -> io::Result<String> {
    let limit = u64::try_from(script::MAX_SCRIPT)
        .unwrap_or(u64::MAX)
        .saturating_add(1);
    let mut text = String::new();
    if path == "-" {
        io::stdin().lock().take(limit).read_to_string(&mut text)?;
    } else {
        let file = std::fs::File::open(path)?;
        BufReader::new(file).take(limit).read_to_string(&mut text)?;
    }
    Ok(text)
}

/// Print a JSON summary of a file: what the document shows, plus the raw packet layout.
fn info(path: &str) -> ExitCode {
    let bytes = match NativeFs.read(path) {
        Ok(b) => b,
        Err(e) => {
            eprintln!("error: cannot read {path}: {e}");
            return ExitCode::FAILURE;
        }
    };
    let imported = match sio::import(path, &bytes) {
        Ok(i) => i,
        Err(e) => {
            eprintln!("error: {path}: {e}");
            return ExitCode::FAILURE;
        }
    };
    let m = &imported.molecule;
    let features: Vec<Value> = m
        .features()
        .iter()
        .map(|f| {
            let segments: Vec<String> = f
                .segments
                .iter()
                .map(|s| format!("{}..{}", s.range.start.saturating_add(1), s.range.end))
                .collect();
            json!({ "name": f.name, "type": f.kind, "strand": f.strand, "segments": segments })
        })
        .collect();
    let dna = match &imported.source {
        Source::Dna(src) => {
            let file = src.file();
            let packets: Vec<Value> = file
                .packets()
                .iter()
                .map(|p| json!({ "tag": format!("0x{:02X}", p.tag), "bytes": p.payload.len() }))
                .collect();
            let cookie = file.cookie().ok().map(|c| {
                json!({
                    "sequence_type": c.sequence_type,
                    "export_version": c.export_version,
                    "import_version": c.import_version,
                })
            });
            Some(json!({ "cookie": cookie, "packets": packets }))
        }
        _ => None,
    };
    let out = json!({
        "file": path,
        "format": imported.format.name(),
        "bytes": bytes.len(),
        "name": m.name(),
        "length": m.len(),
        "topology": m.topology(),
        "features": features,
        "description": m.notes().description,
        "warnings": imported.warnings,
        "dna": dna,
    });
    println!(
        "{}",
        serde_json::to_string_pretty(&out).unwrap_or_else(|_| out.to_string())
    );
    ExitCode::SUCCESS
}
