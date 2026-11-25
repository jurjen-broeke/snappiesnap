//! Construction history: how the open construct was made, and its earlier steps as molecules.

use schemars::JsonSchema;
use serde::Deserialize;
use serde_json::{Value, json};

use crate::{CmdError, CommandSpec, Document, ErrorCode, NoParams, Session};

pub(super) fn specs() -> Vec<CommandSpec> {
    vec![
        CommandSpec::query(
            "history_tree",
            "Construction History",
            "How the open construct was made (read from a .dna file, or recorded by cloning, assembly, mutagenesis and PCR commands): each step's operation, details, name, length and the steps it was made from. `available` says whether a copy of the step's molecule is kept, so it can be opened. `steps` is null if there is no history.",
            history_tree,
        )
        .example(json!({})),
        CommandSpec::session(
            "history_open",
            "Open History Step",
            "Open a step of the construction history (by `id` from history_tree) in a new tab, with the features, primers and notes stored with it.",
            history_open,
        )
        .example(json!({"id": 1})),
        CommandSpec::edit(
            "remove_history",
            "Remove Construction History",
            "Remove the construction history from the open document (undoable); the next save writes the file without it, which can make .dna files much smaller.",
            remove_history,
        )
        .menu("Edit/Remove Construction History")
        .example(json!({})),
    ]
}

fn remove_history(s: &mut Session, _: NoParams) -> Result<Value, CmdError> {
    let m = s.molecule_mut()?;
    let removed = m.take_history().map_or(0, |h| {
        let mut n = 0usize;
        h.walk(&mut |_| n += 1);
        n
    });
    if removed == 0 {
        return Err(CmdError::new(
            ErrorCode::NotFound,
            "the document has no construction history",
        ));
    }
    Ok(json!({ "removed_steps": removed }))
}

fn history_tree(s: &Session, _: NoParams) -> Result<Value, CmdError> {
    match s.molecule()?.history() {
        None => Ok(json!({ "steps": null, "reason": "the document has no construction history" })),
        Some(h) => Ok(json!({ "steps": h.steps() })),
    }
}

#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
struct OpenParams {
    /// Step id, from history_tree.
    id: u32,
}

fn history_open(s: &mut Session, p: OpenParams) -> Result<Value, CmdError> {
    let history = s.molecule()?.history().ok_or_else(|| {
        CmdError::new(
            ErrorCode::NotFound,
            "the document has no construction history",
        )
    })?;
    let mut node = history.find(p.id).cloned().ok_or_else(|| {
        CmdError::new(ErrorCode::NotFound, format!("no history step {}", p.id))
            .with_hint("history_tree lists the steps and their ids")
    })?;
    if p.id == history.id {
        return Err(CmdError::invalid_params(
            "that step is the open document itself",
        ));
    }
    let mut m = node.molecule.take().map(|m| *m).ok_or_else(|| {
        CmdError::new(
            ErrorCode::NotFound,
            format!("the sequence of step {} isn't kept", p.id),
        )
    })?;
    // The step keeps the history that led to it.
    if !node.operation.is_empty() || !node.inputs.is_empty() {
        m.set_history(Some(node))?;
    }
    let result = json!({
        "name": m.name(),
        "length": m.len(),
        "topology": m.topology(),
        "features": m.features().len(),
    });
    let mut doc = Document::new(m);
    // A molecule from the history, not saved as a file of its own yet.
    doc.dirty = true;
    s.open_document(doc);
    Ok(result)
}
