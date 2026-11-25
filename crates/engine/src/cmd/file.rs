//! Creating, opening, switching between and closing documents. Several documents can be open
//! (one per tab); commands work on the active one.

use schemars::JsonSchema;
use serde::Deserialize;
use serde_json::{Value, json};
use snappiesnap_doc::{Alphabet, Molecule, Seq, Topology};

use super::{check_discard, residue_bytes};
use crate::NoParams;
use crate::{CmdError, CommandSpec, Document, ErrorCode, Level, Session};

pub(super) fn specs() -> Vec<CommandSpec> {
    vec![
        CommandSpec::session(
            "new",
            "New Sequence",
            "Create a new document in a new tab and make it active.",
            new,
        )
        .menu("File/New")
        .shortcut("Cmd+N")
        .example(json!({"name": "pDemo", "residues": "ACGTACGT", "topology": "circular"})),
        CommandSpec::session(
            "open",
            "Open",
            "Open a file (SnapGene .dna, GenBank, EMBL, GFF3 with a ##FASTA section, SBOL 2, FASTA) in a new tab, or switch to its tab if it is already open. Anything the document can't show is kept in the file for saving and listed in `warnings`.",
            open,
        )
        .menu("File/Open…")
        .shortcut("Cmd+O")
        .example(json!({"path": "example.dna"})),
        CommandSpec::session(
            "save",
            "Save",
            "Save the document (to its own path, or to `path`; the format follows the extension: .dna, .gb/.gbk, .fa/.fasta). Content read from a file in the same format is written back byte for byte; only edits are re-encoded. What the format can't hold is listed in `warnings`. With `without_history`, writes a copy to `path` without the construction history (the open document keeps its history and its own file).",
            save,
        )
        .menu("File/Save")
        .shortcut("Cmd+S")
        .needs_document(true)
        .example(json!({"path": "copy.dna"}))
        .example(json!({"path": "lean.dna", "without_history": true})),
        CommandSpec::session(
            "close",
            "Close",
            "Close the active document (refused if it has unsaved changes, unless `discard_changes`); the tab next to it becomes active.",
            close,
        )
        .menu("File/Close")
        .shortcut("Cmd+W")
        .needs_document(true)
        .example(json!({"discard_changes": true})),
        CommandSpec::query(
            "documents",
            "Open Documents",
            "The open documents in tab order: index, name, path, unsaved changes, and which is active.",
            documents,
        )
        .needs_document(false)
        .example(json!({})),
        CommandSpec::session(
            "switch_document",
            "Switch Document",
            "Make the open document at `index` (from `documents`) the active one.",
            switch_document,
        )
        .example(json!({"index": 0})),
    ]
}

#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
struct NewParams {
    /// Name of the molecule.
    #[serde(default = "default_name")]
    name: String,
    /// Initial residues (IUPAC letters; whitespace is ignored).
    #[serde(default)]
    residues: String,
    #[serde(default = "default_alphabet")]
    alphabet: Alphabet,
    #[serde(default)]
    topology: Topology,
    /// Accepted for older scripts; documents open in new tabs, so nothing is replaced.
    #[serde(default)]
    #[allow(dead_code)]
    discard_changes: bool,
}

fn default_name() -> String {
    "Untitled".to_string()
}

fn default_alphabet() -> Alphabet {
    Alphabet::Dna
}

fn new(s: &mut Session, p: NewParams) -> Result<Value, CmdError> {
    let seq = Seq::new(p.alphabet, residue_bytes(&p.residues))?;
    let molecule = Molecule::new(&p.name, seq, p.topology)?;
    let result = json!({ "name": molecule.name(), "length": molecule.len() });
    s.open_document(Document::new(molecule));
    Ok(result)
}

#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
struct CloseParams {
    /// Close even if the document has unsaved changes.
    #[serde(default)]
    discard_changes: bool,
}

fn close(s: &mut Session, p: CloseParams) -> Result<Value, CmdError> {
    check_discard(s, p.discard_changes, "closing")?;
    s.close_document();
    Ok(json!({}))
}

#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
struct OpenParams {
    /// Path of the file to open.
    path: String,
    /// Accepted for older scripts; documents open in new tabs, so nothing is replaced.
    #[serde(default)]
    #[allow(dead_code)]
    discard_changes: bool,
}

fn open(s: &mut Session, p: OpenParams) -> Result<Value, CmdError> {
    if let Some(index) = s.document_with_path(&p.path) {
        s.switch_document(index)?;
        let m = s.molecule()?;
        return Ok(json!({
            "name": m.name(),
            "length": m.len(),
            "topology": m.topology(),
            "features": m.features().len(),
            "already_open": true,
            "warnings": [],
        }));
    }
    let bytes = s
        .platform()
        .read(&p.path)
        .map_err(|e| CmdError::new(ErrorCode::Io, format!("cannot read {}: {e}", p.path)))?;
    let imported = snappiesnap_io::import(&p.path, &bytes)?;
    for w in &imported.warnings {
        s.log(Level::Warn, "open", format!("{}: {w}", p.path));
    }
    let m = &imported.molecule;
    let matched = match &imported.source {
        snappiesnap_io::Source::Dna(src) => src.matched_reads().len(),
        _ => 0,
    };
    let result = json!({
        "name": m.name(),
        "format": imported.format.name(),
        "length": m.len(),
        "topology": m.topology(),
        "features": m.features().len(),
        "matched_reads": matched,
        "warnings": imported.warnings,
    });
    s.open_document(Document {
        molecule: imported.molecule,
        path: Some(p.path),
        dirty: false,
        selection: None,
        source: Some(imported.source),
    });
    Ok(result)
}

#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
struct SaveParams {
    /// Where to save; defaults to the path the document was opened from or last saved to.
    #[serde(default)]
    path: Option<String>,
    /// Write a copy to `path` without the construction history (which holds every ancestor's
    /// sequence); the open document is left as it is.
    #[serde(default)]
    without_history: bool,
}

fn save(s: &mut Session, p: SaveParams) -> Result<Value, CmdError> {
    let doc = s.document().ok_or_else(CmdError::no_document)?;
    if p.without_history {
        return save_copy_without_history(s, p.path);
    }
    let path = p.path.or_else(|| doc.path.clone()).ok_or_else(|| {
        CmdError::invalid_params("this document has never been saved; pass a path")
            .with_hint(r#"for example {"path": "/path/to/plasmid.dna"}"#)
    })?;
    let format = snappiesnap_io::format_for_path(&path).ok_or_else(|| {
        CmdError::unsupported(format!(
            "saving {path}: only .dna files can be written so far"
        ))
    })?;
    let options = snappiesnap_io::ExportOptions {
        date: s
            .platform()
            .today()
            .and_then(|(y, m, d)| snappiesnap_io::genbank_date(y, m, d)),
    };
    let exported =
        snappiesnap_io::export_with(&doc.molecule, doc.source.as_ref(), format, &options)?;
    s.platform()
        .write(&path, &exported.bytes)
        .map_err(|e| CmdError::new(ErrorCode::Io, format!("cannot write {path}: {e}")))?;
    for w in &exported.warnings {
        s.log(Level::Warn, "save", format!("{path}: {w}"));
    }
    if let Some(doc) = s.document_mut() {
        doc.path = Some(path.clone());
        doc.dirty = false;
        doc.source = Some(exported.source);
    }
    Ok(json!({
        "path": path,
        "format": format.name(),
        "bytes": exported.bytes.len(),
        "warnings": exported.warnings,
    }))
}

/// A copy of the document without its construction history, at `path` (not the document's own
/// file, which keeps the history).
fn save_copy_without_history(s: &mut Session, path: Option<String>) -> Result<Value, CmdError> {
    let doc = s.document().ok_or_else(CmdError::no_document)?;
    let path = path.ok_or_else(|| {
        CmdError::invalid_params("give the `path` of the copy")
            .with_hint("the open document keeps its history and its own file")
    })?;
    if doc.path.as_deref() == Some(path.as_str()) {
        return Err(CmdError::invalid_params(
            "that is the document's own file: save the copy elsewhere, or use remove_history and save",
        ));
    }
    let format = snappiesnap_io::format_for_path(&path).ok_or_else(|| {
        CmdError::unsupported(format!(
            "saving {path}: only .dna, GenBank, EMBL, GFF3, SBOL (.sbol, .xml) and FASTA files can be written"
        ))
    })?;
    let mut copy = doc.molecule.clone();
    copy.take_history();
    let options = snappiesnap_io::ExportOptions {
        date: s
            .platform()
            .today()
            .and_then(|(y, m, d)| snappiesnap_io::genbank_date(y, m, d)),
    };
    let exported = snappiesnap_io::export_with(&copy, doc.source.as_ref(), format, &options)?;
    s.platform()
        .write(&path, &exported.bytes)
        .map_err(|e| CmdError::new(ErrorCode::Io, format!("cannot write {path}: {e}")))?;
    Ok(json!({
        "path": path,
        "format": format.name(),
        "bytes": exported.bytes.len(),
        "warnings": exported.warnings,
        "copy": true,
    }))
}

fn documents(s: &Session, _: NoParams) -> Result<Value, CmdError> {
    Ok(json!({ "documents": s.documents() }))
}

#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
struct SwitchParams {
    /// Tab index, from `documents`.
    index: usize,
}

fn switch_document(s: &mut Session, p: SwitchParams) -> Result<Value, CmdError> {
    s.switch_document(p.index)?;
    let m = s.molecule()?;
    Ok(json!({ "index": p.index, "name": m.name(), "length": m.len() }))
}
