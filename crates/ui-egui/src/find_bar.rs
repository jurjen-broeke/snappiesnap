//! The find bar of the sequence view: a query, the number of matches, and stepping through them.
//! Matches come from the `find` command, so agents find exactly what the user finds.

use eframe::egui::{self, Key, RichText, Ui};
use serde_json::{Value, json};
use snappiesnap_engine::Session;
use snappiesnap_engine::doc::Range;

use crate::theme::Tokens;

#[derive(Debug, Default, Clone)]
pub struct FindState {
    pub open: bool,
    pub query: String,
    pub mismatches: usize,
    /// Matches: range and what it is ("sequence", "reverse", feature or primer name).
    pub matches: Vec<(Range, String)>,
    pub truncated: bool,
    pub error: Option<String>,
    /// The match shown (selected).
    pub current: Option<usize>,
    /// What the matches were computed for: document, revision, query, mismatches.
    key: Option<(Option<u64>, u64, String, usize)>,
    /// Put the cursor in the text field on the next frame.
    pub focus: bool,
}

/// What the bar asks the app to do.
pub enum FindAction {
    /// Select this match and scroll to it.
    Show(Range),
}

impl FindState {
    /// Open the bar (or focus it if it's open).
    pub fn open(&mut self) {
        self.open = true;
        self.focus = true;
    }

    /// Recompute the matches if the document or the query changed.
    fn refresh(&mut self, session: &mut Session) {
        let key = (
            session.active_document_id(),
            session.revision(),
            self.query.clone(),
            self.mismatches,
        );
        if self.key.as_ref() == Some(&key) {
            return;
        }
        self.key = Some(key);
        self.matches.clear();
        self.error = None;
        self.truncated = false;
        self.current = None;
        if self.query.trim().is_empty() {
            return;
        }
        match session.execute(
            "find",
            json!({ "query": self.query, "mismatches": self.mismatches }),
        ) {
            Ok(v) => {
                self.truncated = v.get("truncated").and_then(Value::as_bool) == Some(true);
                for m in v
                    .get("matches")
                    .and_then(Value::as_array)
                    .into_iter()
                    .flatten()
                {
                    let Some(range) = m
                        .get("range")
                        .and_then(|r| serde_json::from_value::<Range>(r.clone()).ok())
                    else {
                        continue;
                    };
                    let label = match m.get("kind").and_then(Value::as_str) {
                        Some("sequence") => {
                            let strand = m.get("strand").and_then(Value::as_str).unwrap_or("");
                            let mm = m.get("mismatches").and_then(Value::as_u64).unwrap_or(0);
                            let mut l = if strand == "reverse" {
                                "sequence, bottom strand".to_string()
                            } else {
                                "sequence".to_string()
                            };
                            if mm > 0 {
                                l.push_str(&format!(", {mm} mismatch(es)"));
                            }
                            l
                        }
                        Some(kind) => format!(
                            "{kind} {}",
                            m.get("name").and_then(Value::as_str).unwrap_or("")
                        ),
                        None => String::new(),
                    };
                    self.matches.push((range, label));
                }
            }
            Err(e) => self.error = Some(e.message),
        }
    }

    /// Step to the next (or previous) match.
    fn step(&mut self, forward: bool) -> Option<FindAction> {
        let n = self.matches.len();
        if n == 0 {
            return None;
        }
        let next = match (self.current, forward) {
            (None, true) => 0,
            (None, false) => n - 1,
            (Some(c), true) => (c + 1) % n,
            (Some(c), false) => (c + n - 1) % n,
        };
        self.current = Some(next);
        self.matches.get(next).map(|(r, _)| FindAction::Show(*r))
    }
}

/// Draw the bar (when open); returns the match to show.
pub fn show(
    ui: &mut Ui,
    t: &Tokens,
    session: &mut Session,
    state: &mut FindState,
) -> Option<FindAction> {
    if !state.open {
        return None;
    }
    let mut action = None;
    let before = (state.query.clone(), state.mismatches);
    ui.horizontal(|ui| {
        ui.label(RichText::new("Find").color(t.text_muted));
        let field = ui.add(
            egui::TextEdit::singleline(&mut state.query)
                .hint_text("sequence (IUPAC), or a feature or primer name")
                .desired_width(320.0),
        );
        if state.focus {
            field.request_focus();
            state.focus = false;
        }
        let enter = field.lost_focus() && ui.input(|i| i.key_pressed(Key::Enter));
        let back = ui.input(|i| i.modifiers.shift);
        egui::ComboBox::from_id_salt("find-mismatches")
            .width(110.0)
            .selected_text(match state.mismatches {
                0 => "exact".to_string(),
                n => format!("{n} mismatch(es)"),
            })
            .show_ui(ui, |ui| {
                ui.selectable_value(&mut state.mismatches, 0, "exact");
                for n in 1..=3 {
                    ui.selectable_value(&mut state.mismatches, n, format!("{n} mismatch(es)"));
                }
            });
        if (state.query.clone(), state.mismatches) != before {
            state.key = None;
        }
        state.refresh(session);
        let n = state.matches.len();
        let count = match (state.current, n) {
            (_, 0) if state.query.trim().is_empty() => String::new(),
            (_, 0) => "no matches".to_string(),
            (Some(c), _) => format!("{} of {n}{}", c + 1, if state.truncated { "+" } else { "" }),
            (None, _) => format!("{n}{} matches", if state.truncated { "+" } else { "" }),
        };
        ui.label(RichText::new(count).color(t.text_muted));
        if ui
            .add_enabled(n > 0, egui::Button::new("<"))
            .on_hover_text("Previous match (Shift+Enter)")
            .clicked()
        {
            action = state.step(false);
        }
        if ui
            .add_enabled(n > 0, egui::Button::new(">"))
            .on_hover_text("Next match (Enter)")
            .clicked()
        {
            action = state.step(true);
        }
        if enter {
            action = state.step(!back);
            state.focus = true;
        }
        if let Some(label) = state
            .current
            .and_then(|c| state.matches.get(c))
            .map(|(_, l)| l.clone())
        {
            ui.label(RichText::new(label).color(t.text_muted));
        }
        if ui
            .small_button("\u{d7}")
            .on_hover_text("Close (Esc)")
            .clicked()
            || ui.input(|i| i.key_pressed(Key::Escape))
        {
            state.open = false;
        }
    });
    if let Some(e) = &state.error {
        ui.label(RichText::new(e).color(t.error));
    }
    action
}

#[cfg(test)]
mod tests {
    use snappiesnap_engine::Document;
    use snappiesnap_engine::doc::{Alphabet, Molecule, Seq, Topology};

    use super::*;

    #[test]
    fn matches_step_round() {
        let mut session = Session::new();
        session.open_document(Document::new(
            Molecule::new(
                "p",
                Seq::new(Alphabet::Dna, "GAATTCAAAGAATTCAAA").unwrap(),
                Topology::Circular,
            )
            .unwrap(),
        ));
        let mut f = FindState {
            query: "GAATTC".into(),
            ..Default::default()
        };
        f.refresh(&mut session);
        assert_eq!(f.matches.len(), 2);
        let starts: Vec<usize> = (0..3)
            .filter_map(|_| f.step(true).map(|FindAction::Show(r)| r.start))
            .collect();
        assert_eq!(starts, vec![0, 9, 0]);
        assert!(matches!(f.step(false), Some(FindAction::Show(r)) if r.start == 9));
        // A new query starts over; nothing found means nothing to step to.
        f.query = "CCCCCCCC".into();
        f.refresh(&mut session);
        assert!(f.matches.is_empty() && f.step(true).is_none());
    }
}
