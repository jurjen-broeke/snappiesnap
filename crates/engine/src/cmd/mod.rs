//! The standard commands, one module per area (AGENTS.md › Parallel agents).

mod align;
mod annotations;
mod binding;
mod blast;
mod clipboard;

/// The text a copy puts on the clipboard: the window uses it for constructs that are not
/// the open document, so a copy is the same text wherever it is made.
pub use clipboard::copy_text;
mod cloning;
mod edit;
mod enzyme_sets;
mod enzymes;
mod experiment;
mod feature;
mod file;
mod files;
mod find;
mod history;
mod library;
mod msa;
mod navigate;
mod notes;
mod oligos;
mod pcr;
mod primer;
mod query;
mod reads;
mod sample;
mod topo;
mod trace;
mod translate;

use crate::{CmdError, CommandSpec};

/// Every standard command, in menu order.
pub(crate) fn all() -> Vec<CommandSpec> {
    let mut v = Vec::new();
    v.extend(file::specs());
    v.extend(files::specs());
    v.extend(edit::specs());
    v.extend(navigate::specs());
    v.extend(feature::specs());
    v.extend(primer::specs());
    v.extend(binding::specs());
    v.extend(oligos::specs());
    v.extend(query::specs());
    v.extend(enzymes::specs());
    v.extend(enzyme_sets::specs());
    v.extend(pcr::specs());
    v.extend(experiment::specs());
    v.extend(translate::specs());
    v.extend(cloning::specs());
    v.extend(history::specs());
    v.extend(topo::specs());
    v.extend(align::specs());
    v.extend(msa::specs());
    v.extend(reads::specs());
    v.extend(sample::specs());
    v.extend(trace::specs());
    v.extend(clipboard::specs());
    v.extend(find::specs());
    v.extend(notes::specs());
    v.extend(library::specs());
    v.extend(annotations::specs());
    v.extend(blast::specs());
    v
}

/// Residues from a parameter: ASCII whitespace (as in pasted sequences) is ignored.
pub(crate) fn residue_bytes(text: &str) -> Vec<u8> {
    text.bytes().filter(|b| !b.is_ascii_whitespace()).collect()
}

/// Refuse to drop unsaved changes unless the caller says so.
pub(crate) fn check_discard(
    session: &crate::Session,
    discard_changes: bool,
    action: &str,
) -> Result<(), CmdError> {
    if session.document().is_some_and(|d| d.dirty) && !discard_changes {
        return Err(CmdError::new(
            crate::ErrorCode::NotEnabled,
            format!("the open document has unsaved changes; {action} would lose them"),
        )
        .with_hint("save first, or pass \"discard_changes\": true"));
    }
    Ok(())
}
