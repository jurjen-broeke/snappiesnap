//! The Primers tab's "Lab oligos" section: the lab's oligo library (imported once from a list,
//! kept between sessions) searched for the oligos that bind the open construct, so the ones
//! ticked can be added as its primers. Uses `oligo_library_find`; changes go back to the shell
//! as [`OligoAction`]s.

use std::collections::BTreeSet;

use eframe::egui::{self, RichText, Sense, Ui};
use serde_json::{Value, json};
use snappiesnap_engine::Session;
use snappiesnap_engine::doc::Range;

use crate::theme::Tokens;

/// What the user asked for in the section.
#[derive(Debug, Clone, PartialEq)]
pub enum OligoAction {
    /// Pick a list file and import it into the library.
    Import,
    /// Pick a folder of oligo files (one primer per file, named after the file) and import
    /// all of them.
    ImportFolder,
    /// Add these library oligos (by name) to the construct as primers.
    Add(Vec<String>),
    /// Select a binding site.
    Select(Range),
    /// Empty the library.
    RemoveAll,
}

/// What a search depends on: session revision, weaker sites, library generation.
type SearchKey = (u64, bool, u64);

#[derive(Default)]
pub struct OligoViewState {
    /// Also list weaker sites (a 10-base 3' match with one mismatch).
    pub weaker: bool,
    /// Oligos ticked to be added, by name.
    pub ticked: BTreeSet<String>,
    /// The search, and what it was for.
    found: Option<(SearchKey, Result<Value, String>)>,
    /// Bumped whenever the library changes, so the search runs again.
    pub generation: u64,
    /// "Remove all" was clicked once and waits for confirmation.
    confirm_remove: bool,
}

impl OligoViewState {
    fn found(&mut self, session: &mut Session) -> Result<Value, String> {
        let key = (session.revision(), self.weaker, self.generation);
        if self.found.as_ref().is_none_or(|(k, _)| *k != key) {
            let result = session
                .execute(
                    "oligo_library_find",
                    json!({ "primary_only": !self.weaker }),
                )
                .map_err(|e| e.message);
            self.found = Some((key, result));
        }
        self.found
            .as_ref()
            .map_or(Ok(Value::Null), |(_, r)| r.clone())
    }
}

/// The section's contents (the caller puts it in a collapsible header).
pub fn show(
    ui: &mut Ui,
    t: &Tokens,
    session: &mut Session,
    state: &mut OligoViewState,
) -> Option<OligoAction> {
    let mut action = None;
    let len = session.document().map_or(0, |d| d.molecule.len());
    let found = match state.found(session) {
        Ok(v) => v,
        Err(e) => {
            ui.label(RichText::new(e).color(t.error));
            if ui.button("Import oligo list\u{2026}").clicked() {
                action = Some(OligoAction::Import);
            }
            return action;
        }
    };
    let total = found.get("library").and_then(Value::as_u64).unwrap_or(0);
    if total == 0 {
        ui.label(
            RichText::new(
                "Import your lab's oligo list once (name and sequence separated by a tab, one \
                 per line, as SnapGene or a spreadsheet exports it; or FASTA). A folder works \
                 too, where each file holds one primer's sequence and is named after it. The \
                 library is kept between sessions, and each construct then shows which of the \
                 oligos bind it.",
            )
            .color(t.text_muted),
        );
        ui.horizontal(|ui| {
            if ui.button("Import oligo list\u{2026}").clicked() {
                action = Some(OligoAction::Import);
            }
            if ui
                .button("Import folder\u{2026}")
                .on_hover_text(
                    "Every oligo list, FASTA and single-primer file in the folder. A file \
                     holding nothing but a sequence is named after the file.",
                )
                .clicked()
            {
                action = Some(OligoAction::ImportFolder);
            }
        });
        return action;
    }
    let oligos = found
        .get("oligos")
        .and_then(Value::as_array)
        .cloned()
        .unwrap_or_default();
    ui.horizontal(|ui| {
        ui.label(
            RichText::new(format!(
                "{} of the {total} lab oligos bind this construct",
                oligos.len()
            ))
            .color(t.text_muted),
        );
        ui.checkbox(&mut state.weaker, "Include weaker sites")
            .on_hover_text("Also sites where only the last 10 bases pair, with up to one mismatch. By default an oligo must pair exactly over at least 15 bases (or its whole length).");
    });
    ui.horizontal(|ui| {
        if ui.button("Import more\u{2026}").clicked() {
            action = Some(OligoAction::Import);
        }
        if ui
            .button("Import folder\u{2026}")
            .on_hover_text(
                "Every oligo list, FASTA and single-primer file in the folder. A file holding \
                 nothing but a sequence is named after the file.",
            )
            .clicked()
        {
            action = Some(OligoAction::ImportFolder);
        }
        if state.confirm_remove {
            ui.label(
                RichText::new(format!("Remove all {total} oligos from the library?"))
                    .color(t.error),
            );
            if ui.button("Remove").clicked() {
                action = Some(OligoAction::RemoveAll);
                state.confirm_remove = false;
            }
            if ui.button("Keep").clicked() {
                state.confirm_remove = false;
            }
        } else if ui
            .button("Empty library\u{2026}")
            .on_hover_text("Remove every oligo from the lab's library (asks first)")
            .clicked()
        {
            state.confirm_remove = true;
        }
    });
    if oligos.is_empty() {
        return action;
    }
    ui.add_space(t.spacing);
    // Names that no longer bind (after an edit) aren't added by mistake.
    state.ticked.retain(|n| {
        oligos
            .iter()
            .any(|o| o.get("name").and_then(Value::as_str) == Some(n))
    });
    egui::ScrollArea::vertical()
        .id_salt("lab-oligos")
        .max_height(260.0)
        .show(ui, |ui| {
            egui::Grid::new("lab-oligos-grid")
                .striped(true)
                .spacing(egui::vec2(t.spacing * 2.0, t.spacing * 0.5))
                .show(ui, |ui| {
                    for h in ["", "Name", "Length", "Tm", "Binds at", ""] {
                        ui.label(RichText::new(h).strong());
                    }
                    ui.end_row();
                    for o in &oligos {
                        let text = |k: &str| o.get(k).and_then(Value::as_str).unwrap_or_default();
                        let name = text("name").to_string();
                        let in_molecule =
                            o.get("in_molecule").and_then(Value::as_bool).unwrap_or(false);
                        let mut on = state.ticked.contains(&name);
                        if ui
                            .add_enabled(!in_molecule, egui::Checkbox::without_text(&mut on))
                            .changed()
                        {
                            if on {
                                state.ticked.insert(name.clone());
                            } else {
                                state.ticked.remove(&name);
                            }
                        }
                        let description = text("description");
                        let r = ui.label(&name);
                        let hover = if description.is_empty() {
                            text("residues").to_string()
                        } else {
                            format!("{}\n{description}", text("residues"))
                        };
                        r.on_hover_text(hover);
                        ui.label(format!(
                            "{} nt",
                            o.get("length").and_then(Value::as_u64).unwrap_or(0)
                        ));
                        ui.label(
                            o.get("tm")
                                .and_then(Value::as_f64)
                                .map_or("–".to_string(), |v| format!("{v:.1} °C")),
                        );
                        ui.horizontal(|ui| {
                            let sites = o.get("sites").and_then(Value::as_array);
                            for site in sites.into_iter().flatten().take(4) {
                                let Some(range) = site
                                    .get("range")
                                    .and_then(|r| serde_json::from_value::<Range>(r.clone()).ok())
                                else {
                                    continue;
                                };
                                let reverse = site.get("strand").and_then(Value::as_str)
                                    == Some("reverse");
                                let mismatches =
                                    site.get("mismatches").and_then(Value::as_u64).unwrap_or(0);
                                let end = if range.end == 0 { len } else { range.end };
                                let mut label = format!(
                                    "{} {}..{}",
                                    if reverse { "rev" } else { "fwd" },
                                    range.start.saturating_add(1),
                                    end
                                );
                                if mismatches > 0 {
                                    label.push_str(&format!(" ({mismatches} mm)"));
                                }
                                let r = ui
                                    .add(egui::Label::new(label).sense(Sense::click()))
                                    .on_hover_text(if reverse {
                                        "Reverse: binds the top strand, reads leftwards. Click to select."
                                    } else {
                                        "Forward: reads rightwards. Click to select."
                                    });
                                if r.clicked() {
                                    action = Some(OligoAction::Select(range));
                                }
                            }
                            if sites.is_some_and(|s| s.len() > 4) {
                                ui.label(RichText::new("\u{2026}").color(t.text_muted));
                            }
                        });
                        if in_molecule {
                            ui.label(RichText::new("a primer already").small().color(t.text_muted));
                        } else {
                            ui.label("");
                        }
                        ui.end_row();
                    }
                });
        });
    ui.add_space(t.spacing);
    let n = state.ticked.len();
    let label = match n {
        0 => "Add as primers".to_string(),
        1 => "Add 1 oligo as a primer".to_string(),
        _ => format!("Add {n} oligos as primers"),
    };
    if ui
        .add_enabled(n > 0, egui::Button::new(label))
        .on_disabled_hover_text("Tick the oligos to add")
        .clicked()
    {
        action = Some(OligoAction::Add(state.ticked.iter().cloned().collect()));
        state.ticked.clear();
    }
    action
}

#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use snappiesnap_engine::doc::{Alphabet, Molecule, Seq, Topology};
    use snappiesnap_engine::{Document, MemoryFs};

    use super::*;

    /// Empty library, then a library with binding oligos; set `SNAPPIESNAP_OLIGOS_PNG=/path.png`
    /// to look at the second.
    #[test]
    fn draws_and_searches() {
        let fs = Arc::new(MemoryFs::new().with_file(
            "list.txt",
            b"M13F\tGTAAAACGACGGCCAGT\nT7\tTAATACGACTCACTATAGGG\tpromoter\n".to_vec(),
        ));
        let mut session = Session::new().with_platform(fs);
        let residues = format!("{0}TAATACGACTCACTATAGGG{0}", "C".repeat(40));
        let m = Molecule::new(
            "pT",
            Seq::new(Alphabet::Dna, residues).unwrap(),
            Topology::Circular,
        )
        .unwrap();
        session.open_document(Document::new(m));
        let mut state = OligoViewState::default();
        let t = Tokens::dark();
        let render = |session: &mut Session, state: &mut OligoViewState| {
            let mut action = None;
            let image = crate::offscreen::render(
                800,
                300,
                1.0,
                |ctx| t.apply(ctx),
                |ui| {
                    action = show(ui, &t, session, state);
                },
            )
            .unwrap();
            (image, action)
        };
        render(&mut session, &mut state);
        assert_eq!(
            state.found.as_ref().unwrap().1.as_ref().unwrap()["library"],
            0
        );
        session
            .execute("oligo_library_import", json!({"path": "list.txt"}))
            .unwrap();
        state.generation += 1;
        state.ticked.insert("T7".into());
        let (image, _) = render(&mut session, &mut state);
        if let Ok(path) = std::env::var("SNAPPIESNAP_OLIGOS_PNG") {
            crate::offscreen::write_png(&path, &image).unwrap();
        }
        let found = state.found.as_ref().unwrap().1.as_ref().unwrap();
        assert_eq!(found["oligos"].as_array().unwrap().len(), 1);
        assert_eq!(found["oligos"][0]["name"], "T7");
        assert!(state.ticked.contains("T7"));
    }
}
