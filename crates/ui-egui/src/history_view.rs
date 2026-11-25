//! The construction history of the open document ("how this construct was made"), as a tree
//! of steps; any step whose sequence the file holds can be opened in a tab of its own.

use eframe::egui::{self, RichText, Ui};
use snappiesnap_engine::Session;
use snappiesnap_engine::doc::HistoryStep as Step;

use crate::theme::Tokens;

/// The steps of the active document, recomputed only when the document changes.
#[derive(Default)]
pub struct HistoryCache {
    key: Option<(Option<u64>, u64)>,
    steps: Option<Step>,
}

impl HistoryCache {
    fn get(&mut self, session: &Session) -> Option<&Step> {
        let key = (session.active_document_id(), session.revision());
        if self.key != Some(key) {
            self.key = Some(key);
            self.steps = session
                .document()
                .and_then(|d| d.molecule.history())
                .map(|h| h.steps());
        }
        self.steps.as_ref()
    }
}

/// Draw the tree; returns the step the user asked to open.
pub fn show(ui: &mut Ui, t: &Tokens, session: &Session, cache: &mut HistoryCache) -> Option<u32> {
    let mut open = None;
    ui.label(RichText::new("How this construct was made").strong());
    match cache.get(session) {
        None => {
            ui.label(
                RichText::new("No construction history. Cloning, assembly, mutagenesis and PCR products record theirs, and .dna files can carry one.")
                    .color(t.text_muted),
            );
        }
        Some(root) => {
            ui.label(
                RichText::new("Each step lists the molecules it was made from. Open shows a step as a molecule of its own.")
                    .color(t.text_muted),
            );
            ui.add_space(t.spacing);
            step(ui, t, root, 0, true, &mut open);
        }
    }
    open
}

fn describe(s: &Step) -> String {
    let op = if s.operation.is_empty() {
        "Starting molecule".to_string()
    } else {
        s.operation.clone()
    };
    let name = s.name.strip_suffix(".dna").unwrap_or(&s.name);
    format!(
        "{op}: {name} ({} bp{})",
        s.length,
        if s.circular { ", circular" } else { "" }
    )
}

fn step(ui: &mut Ui, t: &Tokens, s: &Step, depth: usize, root: bool, open: &mut Option<u32>) {
    let row = |ui: &mut Ui, open: &mut Option<u32>| {
        if root {
            ui.label(RichText::new("this document").color(t.text_muted));
        } else if s.available {
            if ui
                .small_button("Open")
                .on_hover_text("Open this step in a new tab")
                .clicked()
            {
                *open = Some(s.id);
            }
        } else {
            ui.label(RichText::new("sequence not stored").color(t.text_muted));
        }
    };
    let details = |ui: &mut Ui| {
        if !s.details.is_empty() {
            ui.label(RichText::new(&s.details).color(t.text_muted));
        }
    };
    if s.inputs.is_empty() {
        ui.horizontal(|ui| {
            ui.add_space(18.0);
            ui.label(RichText::new(describe(s)).color(t.text));
            row(ui, open);
        });
        return;
    }
    egui::CollapsingHeader::new(RichText::new(describe(s)).color(t.text))
        .id_salt(("history-step", s.id))
        .default_open(depth < 2)
        .show(ui, |ui| {
            ui.horizontal(|ui| row(ui, open));
            details(ui);
            for input in &s.inputs {
                step(ui, t, input, depth + 1, false, open);
            }
        });
}
