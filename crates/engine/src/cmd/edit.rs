//! Sequence edits, topology, selection and undo/redo.

use schemars::JsonSchema;
use serde::Deserialize;
use serde_json::{Value, json};
use snappiesnap_doc::{EditReport, Range, Topology};

use super::residue_bytes;
use crate::{CmdError, CommandSpec, NoParams, Session};

pub(super) fn specs() -> Vec<CommandSpec> {
    vec![
            CommandSpec::session("undo", "Undo", "Undo the last edit.", undo)
                .menu("Edit/Undo")
                .shortcut("Cmd+Z")
                .needs_document(true)
    .example(json!({})),
            CommandSpec::session("redo", "Redo", "Redo the last undone edit.", redo)
                .menu("Edit/Redo")
                .shortcut("Cmd+Shift+Z")
                .needs_document(true)
    .example(json!({})),
            CommandSpec::edit(
                "insert",
                "Insert Residues",
                "Insert residues before position `at` (0-based; `at` = length appends). Features containing the insertion point grow.",
                insert,
            )
            .menu("Edit/Insert Residues…")
    .example(json!({"at": 2, "residues": "GG"})),
            CommandSpec::edit(
                "delete",
                "Delete Residues",
                "Delete a range. Features shrink; features left empty are removed and listed.",
                delete,
            )
            .menu("Edit/Delete Residues…")
    .example(json!({"range": {"start": 10, "end": 2}})),
            CommandSpec::edit(
                "replace",
                "Replace Residues",
                "Replace a range with new residues.",
                replace,
            )
            .menu("Edit/Replace Residues…")
    .example(json!({"range": {"start": 1, "end": 3}, "residues": "TT"})),
            CommandSpec::edit(
                "reverse_complement",
                "Reverse Complement",
                "Reverse-complement the whole molecule; features move with their residues.",
                reverse_complement,
            )
            .menu("Actions/Reverse Complement")
    .example(json!({})),
            CommandSpec::edit(
                "rotate_origin",
                "Set Origin",
                "Make position `origin` the new position 0 of a circular molecule.",
                rotate_origin,
            )
            .menu("Actions/Set Origin…")
    .example(json!({"origin": 3})),
            CommandSpec::edit(
                "set_topology",
                "Set Topology",
                "Make the molecule linear or circular.",
                set_topology,
            )
            .menu("Actions/Topology")
    .example(json!({"topology": "circular"})),
            CommandSpec::edit("rename", "Rename", "Rename the molecule.", rename)
                .menu("File/Rename…")
    .example(json!({"name": "pRenamed"})),
            CommandSpec::session(
                "select",
                "Select",
                "Set or clear the selection (view state, not undoable).",
                select,
            )
            .menu("Edit/Select Range…")
            .needs_document(true)
    .example(json!({"range": {"start": 10, "end": 2}})),
        ]
}

fn report(s: &Session, r: EditReport) -> Result<Value, CmdError> {
    Ok(json!({
        "length": s.molecule()?.len(),
        "removed_features": r.removed_features,
    }))
}

fn undo(s: &mut Session, _: NoParams) -> Result<Value, CmdError> {
    s.undo()?;
    Ok(json!({ "length": s.molecule()?.len() }))
}

fn redo(s: &mut Session, _: NoParams) -> Result<Value, CmdError> {
    s.redo()?;
    Ok(json!({ "length": s.molecule()?.len() }))
}

#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
struct InsertParams {
    /// 0-based position to insert before.
    at: usize,
    /// IUPAC residues to insert (whitespace is ignored).
    residues: String,
}

fn insert(s: &mut Session, p: InsertParams) -> Result<Value, CmdError> {
    let r = s
        .molecule_mut()?
        .insert(p.at, &residue_bytes(&p.residues))?;
    report(s, r)
}

#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
struct RangeParams {
    /// Half-open range; `start > end` wraps past the origin of a circular molecule.
    range: Range,
}

fn delete(s: &mut Session, p: RangeParams) -> Result<Value, CmdError> {
    let r = s.molecule_mut()?.delete(p.range)?;
    report(s, r)
}

#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
struct ReplaceParams {
    range: Range,
    /// Replacement residues (whitespace is ignored).
    residues: String,
}

fn replace(s: &mut Session, p: ReplaceParams) -> Result<Value, CmdError> {
    let r = s
        .molecule_mut()?
        .replace(p.range, &residue_bytes(&p.residues))?;
    report(s, r)
}

fn reverse_complement(s: &mut Session, _: NoParams) -> Result<Value, CmdError> {
    s.molecule_mut()?.reverse_complement()?;
    Ok(json!({}))
}

#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
struct OriginParams {
    /// The position that becomes the new origin.
    origin: usize,
}

fn rotate_origin(s: &mut Session, p: OriginParams) -> Result<Value, CmdError> {
    s.molecule_mut()?.rotate_origin(p.origin)?;
    Ok(json!({}))
}

#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
struct TopologyParams {
    topology: Topology,
}

fn set_topology(s: &mut Session, p: TopologyParams) -> Result<Value, CmdError> {
    s.molecule_mut()?.set_topology(p.topology)?;
    Ok(json!({}))
}

#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
struct RenameParams {
    name: String,
}

fn rename(s: &mut Session, p: RenameParams) -> Result<Value, CmdError> {
    s.molecule_mut()?.rename(&p.name)?;
    Ok(json!({}))
}

#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
struct SelectParams {
    /// The range to select, or null to clear the selection.
    range: Option<Range>,
}

fn select(s: &mut Session, p: SelectParams) -> Result<Value, CmdError> {
    let m = s.molecule()?;
    let range = p
        .range
        .map(|r| r.validate(m.len(), m.topology()))
        .transpose()?;
    if let Some(doc) = s.document_mut() {
        doc.selection = range;
    }
    Ok(json!({ "selection": range }))
}
