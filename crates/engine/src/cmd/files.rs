//! File management for the folder panel and agents: list, rename/move, copy, new folder.
//! Nothing is ever overwritten, and nothing is deleted.

use schemars::JsonSchema;
use serde::Deserialize;
use serde_json::{Value, json};

use crate::{CmdError, CommandSpec, ErrorCode, Session};

pub(super) fn specs() -> Vec<CommandSpec> {
    vec![
        CommandSpec::query(
            "list_folder",
            "List Folder",
            "The entries of a folder (folders first; hidden entries left out).",
            list_folder,
        )
        .needs_document(false)
        .example(json!({"path": "."})),
        CommandSpec::query(
            "file_info",
            "File Info",
            "What a sequence file holds, without opening it as a document: `name`, `length`, `topology`, `alphabet`, and how many `features` and `primers`. With `molecule: true` the whole molecule comes too, so a map can be drawn of a construct before choosing it (large for a chromosome-sized file; `doc.get` reads an open document in ranges instead).",
            file_info,
        )
        .needs_document(false)
        .example(json!({"path": "example.dna"}))
        .example(json!({"path": "example.dna", "molecule": true})),
        CommandSpec::session(
            "file_rename",
            "Rename or Move File",
            "Rename or move a file or folder. Never replaces an existing file. If the open document is renamed, it follows.",
            file_rename,
        )
        .example(json!({"from": "example.dna", "to": "renamed.dna"})),
        CommandSpec::session(
            "file_copy",
            "Duplicate File",
            "Copy a file. Never replaces an existing file.",
            file_copy,
        )
        .example(json!({"from": "example.dna", "to": "example copy.dna"})),
        CommandSpec::session(
            "folder_create",
            "New Folder",
            "Create a folder.",
            folder_create,
        )
        .example(json!({"path": "new folder"})),
    ]
}

fn io_error(what: &str, e: &std::io::Error) -> CmdError {
    if e.kind() == std::io::ErrorKind::AlreadyExists {
        CmdError::new(ErrorCode::InvalidInput, format!("{what}: {e}"))
            .with_hint("choose a name that isn't taken")
    } else {
        CmdError::new(ErrorCode::Io, format!("{what}: {e}"))
    }
}

fn non_empty(path: &str, name: &str) -> Result<(), CmdError> {
    if path.trim().is_empty() {
        return Err(CmdError::invalid_params(format!("{name} can't be empty")));
    }
    Ok(())
}

#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
struct PathParams {
    path: String,
}

#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
struct FileInfoParams {
    /// The sequence file to look into.
    path: String,
    /// Also return the whole molecule, for a map of the construct (default false).
    #[serde(default)]
    molecule: bool,
}

/// Read a sequence file and describe it, without opening it as a document. The file is only
/// read: nothing is kept, and the open document is untouched.
fn file_info(s: &Session, p: FileInfoParams) -> Result<Value, CmdError> {
    if p.path.trim().is_empty() {
        return Err(CmdError::invalid_params("`path` can't be empty"));
    }
    let bytes = s
        .platform()
        .read(&p.path)
        .map_err(|e| io_error(&format!("cannot read {}", p.path), &e))?;
    let m = snappiesnap_io::import(&p.path, &bytes)?.molecule;
    let mut out = json!({
        "path": p.path,
        "name": m.name(),
        "length": m.len(),
        "topology": match m.topology() {
            snappiesnap_doc::Topology::Circular => "circular",
            snappiesnap_doc::Topology::Linear => "linear",
        },
        "alphabet": m.seq().alphabet(),
        "features": m.features().len(),
        "primers": m.primers().len(),
    });
    if p.molecule
        && let Some(o) = out.as_object_mut()
    {
        o.insert(
            "molecule".into(),
            serde_json::to_value(&m).unwrap_or(Value::Null),
        );
    }
    Ok(out)
}

fn list_folder(s: &Session, p: PathParams) -> Result<Value, CmdError> {
    let entries = s
        .platform()
        .list_dir(&p.path)
        .map_err(|e| io_error(&format!("cannot list {}", p.path), &e))?;
    Ok(json!({ "entries": entries }))
}

#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
struct FromTo {
    from: String,
    to: String,
}

fn file_rename(s: &mut Session, p: FromTo) -> Result<Value, CmdError> {
    non_empty(&p.from, "from")?;
    non_empty(&p.to, "to")?;
    if p.from == p.to {
        return Ok(json!({ "path": p.to }));
    }
    s.platform()
        .rename(&p.from, &p.to)
        .map_err(|e| io_error(&format!("cannot rename {} to {}", p.from, p.to), &e))?;
    // Keep the open document pointing at its file (also when a parent folder moved).
    if let Some(doc) = s.document_mut()
        && let Some(path) = doc.path.as_mut()
    {
        if *path == p.from {
            *path = p.to.clone();
        } else if let Some(rest) = path.strip_prefix(&format!("{}/", p.from)) {
            *path = format!("{}/{rest}", p.to);
        }
    }
    Ok(json!({ "path": p.to }))
}

fn file_copy(s: &mut Session, p: FromTo) -> Result<Value, CmdError> {
    non_empty(&p.to, "to")?;
    if s.platform().exists(&p.to) {
        return Err(
            CmdError::new(ErrorCode::InvalidInput, format!("{} already exists", p.to))
                .with_hint("choose a name that isn't taken"),
        );
    }
    let bytes = s
        .platform()
        .read(&p.from)
        .map_err(|e| io_error(&format!("cannot read {}", p.from), &e))?;
    s.platform()
        .write(&p.to, &bytes)
        .map_err(|e| io_error(&format!("cannot write {}", p.to), &e))?;
    Ok(json!({ "path": p.to, "bytes": bytes.len() }))
}

fn folder_create(s: &mut Session, p: PathParams) -> Result<Value, CmdError> {
    non_empty(&p.path, "path")?;
    s.platform()
        .create_dir(&p.path)
        .map_err(|e| io_error(&format!("cannot create {}", p.path), &e))?;
    Ok(json!({ "path": p.path }))
}
