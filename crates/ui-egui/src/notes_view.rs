//! The Notes tab: the molecule's description and named fields, edited and saved through the
//! `set_notes` command (so the edit is undoable and lands in the file on save).

use std::collections::BTreeMap;

use eframe::egui::{self, RichText, Ui};
use serde_json::{Value, json};
use snappiesnap_engine::Session;

use crate::theme::Tokens;

#[derive(Debug, Default)]
pub struct NotesEdit {
    /// The document and revision the form was filled from.
    key: Option<(Option<u64>, u64)>,
    /// The description as shown (HTML tags taken out).
    description: String,
    original_description: String,
    fields: Vec<(String, String)>,
    original: BTreeMap<String, String>,
    new_name: String,
    new_value: String,
}

/// Descriptions in `.dna` files are HTML: the text without tags and with the common entities
/// decoded.
pub fn plain(html: &str) -> String {
    let mut out = String::new();
    let mut in_tag = false;
    for c in html.chars() {
        match c {
            '<' => in_tag = true,
            '>' if in_tag => in_tag = false,
            _ if !in_tag => out.push(c),
            _ => {}
        }
    }
    out.replace("&nbsp;", " ")
        .replace("&lt;", "<")
        .replace("&gt;", ">")
        .replace("&quot;", "\"")
        .replace("&amp;", "&")
        .trim()
        .to_string()
}

impl NotesEdit {
    fn load(&mut self, session: &Session) {
        let key = (session.active_document_id(), session.revision());
        if self.key == Some(key) {
            return;
        }
        self.key = Some(key);
        let Some(doc) = session.document() else {
            return;
        };
        let n = doc.molecule.notes();
        self.description = plain(&n.description);
        self.original_description = self.description.clone();
        self.fields = n
            .fields
            .iter()
            .map(|(k, v)| (k.clone(), v.clone()))
            .collect();
        self.original = n.fields.clone();
    }

    /// The `set_notes` params for what changed (None if nothing did).
    fn changes(&self) -> Option<Value> {
        let mut p = serde_json::Map::new();
        if self.description != self.original_description {
            p.insert("description".into(), json!(self.description));
        }
        let mut fields = serde_json::Map::new();
        let now: BTreeMap<&str, &str> = self
            .fields
            .iter()
            .map(|(k, v)| (k.as_str(), v.as_str()))
            .collect();
        for (k, v) in &self.original {
            match now.get(k.as_str()) {
                None => {
                    fields.insert(k.clone(), Value::Null);
                }
                Some(n) if n != v => {
                    fields.insert(k.clone(), json!(n));
                }
                _ => {}
            }
        }
        for (k, v) in &now {
            if !self.original.contains_key(*k) {
                fields.insert((*k).to_string(), json!(v));
            }
        }
        if !fields.is_empty() {
            p.insert("fields".into(), Value::Object(fields));
        }
        (!p.is_empty()).then_some(Value::Object(p))
    }
}

/// Draw the tab; returns `set_notes` params when the user saves.
pub fn show(ui: &mut Ui, t: &Tokens, session: &Session, state: &mut NotesEdit) -> Option<Value> {
    state.load(session);
    let mut save = None;
    ui.label(RichText::new("Description").strong());
    ui.add(
        egui::TextEdit::multiline(&mut state.description)
            .desired_rows(4)
            .desired_width(f32::INFINITY),
    );
    ui.add_space(t.spacing * 2.0);
    ui.label(RichText::new("Fields").strong());
    let mut remove = None;
    egui::Grid::new("notes-fields")
        .num_columns(3)
        .spacing(egui::vec2(t.spacing * 2.0, t.spacing))
        .show(ui, |ui| {
            for (i, (name, value)) in state.fields.iter_mut().enumerate() {
                ui.label(name.as_str());
                ui.add(egui::TextEdit::singleline(value).desired_width(360.0));
                if ui.small_button("Remove").clicked() {
                    remove = Some(i);
                }
                ui.end_row();
            }
            ui.add(
                egui::TextEdit::singleline(&mut state.new_name)
                    .hint_text("new field, e.g. Organism")
                    .desired_width(160.0),
            );
            ui.add(egui::TextEdit::singleline(&mut state.new_value).desired_width(360.0));
            let name = state.new_name.trim().to_string();
            let ok = !name.is_empty()
                && name
                    .chars()
                    .all(|c| c.is_ascii_alphanumeric() || c == '_' || c == '-')
                && !state.fields.iter().any(|(k, _)| *k == name);
            if ui
                .add_enabled(ok, egui::Button::new("Add"))
                .on_disabled_hover_text("A new name of letters, digits, - and _")
                .clicked()
            {
                state
                    .fields
                    .push((name, std::mem::take(&mut state.new_value)));
                state.new_name.clear();
            }
            ui.end_row();
        });
    if let Some(i) = remove
        && i < state.fields.len()
    {
        state.fields.remove(i);
    }
    ui.add_space(t.spacing * 2.0);
    let changes = state.changes();
    ui.horizontal(|ui| {
        if ui
            .add_enabled(changes.is_some(), egui::Button::new("Save notes"))
            .on_hover_text("Apply the changes to the document (undoable; saved with the file)")
            .clicked()
        {
            save = changes.clone();
        }
        if ui
            .add_enabled(changes.is_some(), egui::Button::new("Revert"))
            .clicked()
        {
            state.key = None;
        }
    });
    ui.label(
        RichText::new("Fields such as Created and UUID come from the file; changes to them are saved as typed.")
            .small()
            .color(t.text_muted),
    );
    save
}

#[cfg(test)]
mod tests {
    use snappiesnap_engine::Document;
    use snappiesnap_engine::doc::{Alphabet, Molecule, Notes, Seq, Topology};

    use super::*;

    #[test]
    fn edits_become_set_notes_params() {
        assert_eq!(
            plain("<html><body>a &amp; b&nbsp;c</body></html>"),
            "a & b c"
        );
        let mut m = Molecule::new(
            "m",
            Seq::new(Alphabet::Dna, "ACGT").unwrap(),
            Topology::Linear,
        )
        .unwrap();
        let mut notes = Notes {
            description: "<html><body>old</body></html>".into(),
            ..Notes::default()
        };
        notes.fields.insert("Organism".into(), "E. coli".into());
        notes.fields.insert("UUID".into(), "abc".into());
        m.set_notes(notes).unwrap();
        let mut session = Session::new();
        session.open_document(Document::new(m));
        let mut e = NotesEdit::default();
        e.load(&session);
        assert_eq!(e.description, "old");
        assert_eq!(e.changes(), None, "loading changes nothing");
        e.description = "new".into();
        e.fields.retain(|(k, _)| k != "UUID");
        e.fields.push(("CreatedBy".into(), "Lab".into()));
        let p = e.changes().unwrap();
        assert_eq!(
            p,
            json!({"description": "new", "fields": {"UUID": null, "CreatedBy": "Lab"}})
        );
        session.execute("set_notes", p).unwrap();
        let n = session.document().unwrap().molecule.notes().clone();
        assert_eq!(n.description, "new");
        assert_eq!(n.fields.get("CreatedBy").map(String::as_str), Some("Lab"));
        assert!(!n.fields.contains_key("UUID"));
    }
}
