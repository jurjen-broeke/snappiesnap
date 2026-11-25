//! The molecule's notes: its description and named fields (author, organism, ...).

use std::collections::BTreeMap;

use schemars::JsonSchema;
use serde::Deserialize;
use serde_json::{Value, json};

use crate::{CmdError, CommandSpec, Session};

pub(super) fn specs() -> Vec<CommandSpec> {
    vec![
        CommandSpec::edit(
            "set_notes",
            "Edit Notes",
            "Change the molecule's notes: `description` replaces the description; `fields` sets named fields (a null value removes one). What isn't given stays. Field names are kept as given (as in .dna files: CreatedBy, Organism, ...).",
            set_notes,
        )
        .menu("Edit/Notes…")
        .example(json!({"description": "a test construct"}))
        .example(json!({"fields": {"Organism": "E. coli", "Old": null}})),
    ]
}

#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
struct NotesParams {
    description: Option<String>,
    /// Fields to set; a null value removes the field.
    #[serde(default)]
    fields: BTreeMap<String, Option<String>>,
}

fn set_notes(s: &mut Session, p: NotesParams) -> Result<Value, CmdError> {
    let m = s.molecule_mut()?;
    let mut notes = m.notes().clone();
    if let Some(d) = p.description {
        notes.description = d;
    }
    for (name, value) in p.fields {
        let name = name.trim().to_string();
        if name.is_empty() {
            return Err(CmdError::invalid_params("a field name can't be empty"));
        }
        if name == "Description" {
            return Err(CmdError::invalid_params(
                "the description is set with `description`",
            ));
        }
        if !name
            .chars()
            .all(|c| c.is_ascii_alphanumeric() || c == '_' || c == '-')
        {
            return Err(CmdError::invalid_params(format!(
                "field name {name:?}: use letters, digits, - and _"
            )));
        }
        match value {
            Some(v) => {
                notes.fields.insert(name, v);
            }
            None => {
                notes.fields.remove(&name);
            }
        }
    }
    m.set_notes(notes)?;
    let n = m.notes();
    Ok(json!({ "description": n.description, "fields": n.fields }))
}
