//! Multiple sequence alignment in the Align tab: sequences from this document, the other open
//! tabs, files and pasted text, aligned by the `align_sequences` command and shown in blocks
//! with a consensus.

use eframe::egui::{self, RichText, Ui, text::LayoutJob};
use serde_json::{Value, json};

use crate::theme::Tokens;

/// What the view asks the app to do.
pub enum MsaAction {
    /// Run `align_sequences` with these params.
    Align(Value),
    /// Let the user choose sequence files to add, wherever the build gets files from.
    AddFiles,
}

#[derive(Debug)]
pub struct MsaState {
    /// Include the open document.
    pub document: bool,
    /// Other documents open in tabs, by id.
    pub documents: Vec<u64>,
    /// Files chosen from the open folder, or loaded with the file picker.
    pub files: Vec<String>,
    /// Pasted sequences: FASTA, or plain sequences separated by blank lines.
    pub text: String,
    /// The last answer of `align_sequences`, or its error.
    pub result: Option<Result<Value, String>>,
}

impl Default for MsaState {
    fn default() -> Self {
        MsaState {
            document: true,
            documents: Vec::new(),
            files: Vec::new(),
            text: String::new(),
            result: None,
        }
    }
}

impl MsaState {
    pub fn params(&self) -> Value {
        let mut p = json!({
            "include_document": self.document,
            "documents": self.documents,
            "files": self.files,
        });
        if let Some(o) = p.as_object_mut()
            && !self.text.trim().is_empty()
        {
            o.insert("text".into(), json!(self.text));
        }
        p
    }

    /// How many sequences the inputs hold, roughly (for enabling the button).
    fn count(&self) -> usize {
        let pasted = if self.text.trim_start().starts_with('>') {
            self.text.matches('>').count()
        } else {
            self.text
                .split("\n\n")
                .filter(|b| b.chars().any(|c| c.is_ascii_alphabetic()))
                .count()
        };
        usize::from(self.document) + self.documents.len() + self.files.len() + pasted
    }
}

/// The inputs and the result; returns what the user asked for. Side by side on a wide page:
/// the alignment itself is wide, and under the inputs it left most of the page empty.
pub fn show(
    ui: &mut Ui,
    t: &Tokens,
    state: &mut MsaState,
    document: &str,
    files: &[String],
    tabs: &[(u64, String)],
) -> Option<MsaAction> {
    let full = ui.available_width();
    let height = ui.available_height();
    let Some(left) = crate::two_columns(full) else {
        let mut action = inputs(ui, t, state, document, files, tabs);
        match &state.result {
            None => {}
            Some(Err(e)) => {
                ui.label(RichText::new(e).color(t.error));
            }
            Some(Ok(r)) => result(ui, t, r),
        }
        if action.is_none() {
            action = None;
        }
        return action;
    };
    let mut action = None;
    ui.horizontal_top(|ui| {
        ui.allocate_ui_with_layout(
            egui::vec2(left, height),
            egui::Layout::top_down(egui::Align::Min),
            |ui| {
                ui.set_max_width(left);
                egui::ScrollArea::vertical()
                    .id_salt("msa-inputs")
                    .max_width(left)
                    .show(ui, |ui| {
                        ui.set_max_width(left);
                        action = inputs(ui, t, state, document, files, tabs);
                    });
            },
        );
        ui.separator();
        egui::ScrollArea::vertical()
            .id_salt("msa-result")
            .show(ui, |ui| {
                ui.vertical(|ui| match &state.result {
                    None => {
                        ui.label(RichText::new("The alignment appears here.").color(t.text_muted));
                    }
                    Some(Err(e)) => {
                        ui.label(RichText::new(e).color(t.error));
                    }
                    Some(Ok(r)) => result(ui, t, r),
                });
            });
    });
    action
}

/// What to align: this document, the other tabs, files and pasted text.
fn inputs(
    ui: &mut Ui,
    t: &Tokens,
    state: &mut MsaState,
    document: &str,
    files: &[String],
    tabs: &[(u64, String)],
) -> Option<MsaAction> {
    let mut action = None;
    ui.label(
        RichText::new(
            "Align several sequences to each other: this document, the other open tabs, files \
             and pasted sequences (FASTA, or plain sequences separated by blank lines).",
        )
        .color(t.text_muted),
    );
    ui.add_space(t.spacing);
    ui.checkbox(&mut state.document, format!("This document ({document})"));
    if !tabs.is_empty() {
        ui.label(RichText::new("Other open tabs").strong());
        egui::ScrollArea::vertical()
            .id_salt("msa-tabs")
            .max_height(120.0)
            .show(ui, |ui| {
                for (id, name) in tabs {
                    let mut on = state.documents.contains(id);
                    if ui.checkbox(&mut on, name).changed() {
                        if on {
                            state.documents.push(*id);
                        } else {
                            state.documents.retain(|x| x != id);
                        }
                    }
                }
            });
    }
    // Tabs that were closed since they were ticked.
    state
        .documents
        .retain(|id| tabs.iter().any(|(i, _)| i == id));
    // Files from the open folder, and any loaded with the picker (the only way in the browser,
    // which gives the page no folder to list).
    let loaded: Vec<String> = state
        .files
        .iter()
        .filter(|f| !files.contains(f))
        .cloned()
        .collect();
    ui.horizontal(|ui| {
        ui.label(RichText::new("Files").strong());
        if ui
            .button("Add files\u{2026}")
            .on_hover_text("Choose one or more sequence files to align")
            .clicked()
        {
            action = Some(MsaAction::AddFiles);
        }
    });
    if files.is_empty() && loaded.is_empty() {
        let hint = if tabs.is_empty() {
            "None chosen yet: use Add files\u{2026}, or paste sequences below."
        } else {
            "None chosen yet: use Add files\u{2026}, or tick the tabs above."
        };
        ui.label(RichText::new(hint).color(t.text_muted));
    } else {
        egui::ScrollArea::vertical()
            .id_salt("msa-files")
            .max_height(140.0)
            .show(ui, |ui| {
                for f in files.iter().chain(loaded.iter()) {
                    let mut on = state.files.contains(f);
                    if ui.checkbox(&mut on, crate::sidebar::base_name(f)).changed() {
                        if on {
                            state.files.push(f.clone());
                        } else {
                            state.files.retain(|x| x != f);
                        }
                    }
                }
            });
    }
    ui.label(RichText::new("Pasted sequences").strong());
    ui.add(
        egui::TextEdit::multiline(&mut state.text)
            .hint_text(">name\nACGT\u{2026}\n>other\nACGT\u{2026}")
            .desired_rows(4)
            .desired_width(560.0)
            .font(egui::TextStyle::Monospace),
    );
    ui.add_space(t.spacing);
    let n = state.count();
    ui.horizontal(|ui| {
        if ui
            .add_enabled(n >= 2, egui::Button::new("Align"))
            .on_disabled_hover_text("Choose at least two sequences")
            .clicked()
        {
            action = Some(MsaAction::Align(state.params()));
        }
        ui.label(RichText::new(format!("{n} sequence(s)")).color(t.text_muted));
    });
    ui.add_space(t.spacing);
    action
}

fn result(ui: &mut Ui, t: &Tokens, r: &Value) {
    let empty = Vec::new();
    let rows = r.get("rows").and_then(Value::as_array).unwrap_or(&empty);
    let columns = r.get("columns").and_then(Value::as_u64).unwrap_or(0);
    let identical = r
        .get("identical_columns")
        .and_then(Value::as_u64)
        .unwrap_or(0);
    let kind = if r.get("kind").and_then(Value::as_str) == Some("protein") {
        "proteins"
    } else {
        "nucleotides"
    };
    ui.label(
        RichText::new(format!(
            "{} sequences ({kind}), {columns} columns; {identical} identical in all",
            rows.len()
        ))
        .strong(),
    );
    egui::Grid::new("msa-identity")
        .num_columns(3)
        .spacing(egui::vec2(t.spacing * 3.0, 2.0))
        .show(ui, |ui| {
            for row in rows {
                let name = row.get("name").and_then(Value::as_str).unwrap_or("");
                let len = row.get("length").and_then(Value::as_u64).unwrap_or(0);
                let id = row.get("identity").and_then(Value::as_f64).unwrap_or(0.0) * 100.0;
                ui.label(name);
                ui.label(RichText::new(format!("{len}")).color(t.text_muted));
                ui.label(RichText::new(format!("{id:.1}% to consensus")).color(t.text_muted));
                ui.end_row();
            }
        });
    ui.add_space(t.spacing);
    let names: Vec<&str> = rows
        .iter()
        .map(|row| row.get("name").and_then(Value::as_str).unwrap_or(""))
        .collect();
    let lines: Vec<Vec<u8>> = rows
        .iter()
        .map(|row| {
            row.get("row")
                .and_then(Value::as_str)
                .unwrap_or("")
                .as_bytes()
                .to_vec()
        })
        .collect();
    let consensus = r
        .get("consensus")
        .and_then(Value::as_str)
        .unwrap_or("")
        .as_bytes();
    blocks(ui, t, &names, &lines, consensus);
}

/// The alignment in blocks of 60 columns: one line per sequence (name, residues, position of its
/// last residue), then the consensus. Residues that differ from the consensus stand out.
fn blocks(ui: &mut Ui, t: &Tokens, names: &[&str], rows: &[Vec<u8>], consensus: &[u8]) {
    const NAME: usize = 18;
    let width = consensus.len();
    let font = egui::FontId::monospace(12.0);
    let (line, cw) = ui.fonts_mut(|f| (f.row_height(&font), f.glyph_width(&font, 'A').max(1.0)));
    // As many residues to a block as the page is wide enough for.
    let per_block = crate::alignment_columns(ui.available_width() - cw * (NAME + 12) as f32, cw);
    let blocks = width.div_ceil(per_block);
    let fmt = |color| egui::TextFormat {
        font_id: font.clone(),
        color,
        ..Default::default()
    };
    let label = |n: &str| {
        let short: String = n.chars().take(NAME).collect();
        format!("{short:<w$}  ", w = NAME)
    };
    egui::ScrollArea::vertical().id_salt("msa-rows").show_rows(
        ui,
        line * (rows.len() + 2) as f32,
        blocks,
        |ui, range| {
            for b in range {
                let from = b * per_block;
                let to = (from + per_block).min(width);
                let cons = consensus.get(from..to).unwrap_or_default();
                let mut job = LayoutJob::default();
                for (name, row) in names.iter().zip(rows) {
                    job.append(&label(name), 0.0, fmt(t.text_muted));
                    let part = row.get(from..to).unwrap_or_default();
                    // Runs of the same colour.
                    let mut k = 0;
                    while k < part.len() {
                        let class = |i: usize| {
                            let (b, c) = (part.get(i).copied(), cons.get(i).copied());
                            match (b, c) {
                                (Some(b'-'), _) => 0,
                                (Some(b), Some(c)) if b.eq_ignore_ascii_case(&c) => 1,
                                _ => 2,
                            }
                        };
                        let cl = class(k);
                        let mut e = k + 1;
                        while e < part.len() && class(e) == cl {
                            e += 1;
                        }
                        let text = String::from_utf8_lossy(part.get(k..e).unwrap_or_default());
                        let color = match cl {
                            0 => t.text_muted,
                            1 => t.text,
                            _ => t.error,
                        };
                        job.append(&text, 0.0, fmt(color));
                        k = e;
                    }
                    let done = row
                        .get(..to)
                        .map_or(0, |s| s.iter().filter(|c| **c != b'-').count());
                    job.append(&format!("  {done}\n"), 0.0, fmt(t.text_muted));
                }
                job.append(&label("consensus"), 0.0, fmt(t.accent));
                job.append(&String::from_utf8_lossy(cons), 0.0, fmt(t.accent));
                ui.label(job);
                ui.add_space(line);
            }
        },
    );
}

#[cfg(test)]
mod tests {
    use snappiesnap_engine::doc::{Alphabet, Molecule, Seq, Topology};
    use snappiesnap_engine::{Document, Session};

    use super::*;

    /// Three sequences pasted with the document: aligned and drawn. Set
    /// `SNAPPIESNAP_MSA_PNG=/path.png` to look at it.
    #[test]
    fn multiple_alignment_draws() {
        let mut session = Session::new();
        let base = "ATGGCTAGCAAAGGAGAAGAACTTTTCACTGGAGTTGTCCCAATTCTTGTTGAATTAGATGGTGATGTTAATGG";
        session.open_document(Document::new(
            Molecule::new(
                "gfp start",
                Seq::new(Alphabet::Dna, base).unwrap(),
                Topology::Linear,
            )
            .unwrap(),
        ));
        let mut state = MsaState {
            text: format!(
                ">del\n{}\n>ins\n{}\n>sub\n{}\n",
                base.replacen("CTTTTCACT", "CTTTCACT", 1),
                base.replacen("GGAGTTGTC", "GGAGTTTGTC", 1),
                base.replacen("AATTCTTG", "AATACTTG", 1)
            ),
            ..Default::default()
        };
        assert_eq!(state.count(), 4);
        let r = session
            .execute("align_sequences", state.params())
            .map_err(|e| e.message);
        let rows = r.as_ref().unwrap()["rows"].as_array().unwrap().clone();
        assert_eq!(rows.len(), 4);
        assert_eq!(rows[0]["name"], "gfp start");
        state.result = Some(r);
        let t = Tokens::light();
        let image = crate::offscreen::render(
            1500,
            800,
            1.0,
            |ctx| t.apply(ctx),
            |ui| {
                egui::Frame::new().fill(t.panel).show(ui, |ui| {
                    let _ = show(ui, &t, &mut state, "gfp start", &[], &[]);
                });
            },
        )
        .unwrap();
        if let Ok(path) = std::env::var("SNAPPIESNAP_MSA_PNG") {
            crate::offscreen::write_png(&path, &image).unwrap();
        }
    }

    /// Sequences come from the other open tabs as well, which is the only way in the browser:
    /// ticking a tab counts it, puts its id in the params, and aligning it works.
    #[test]
    fn open_tabs_are_sequences_too() {
        let mut session = Session::new();
        let base = "ATGGCTAGCAAAGGAGAAGAACTTTTCACTGGAGTTGTCCCAATTCTTGTTGAATTAGATGGTGATGTTAATGG";
        let doc = |name: &str, residues: &str| {
            Document::new(
                Molecule::new(
                    name,
                    Seq::new(Alphabet::Dna, residues).expect("test sequence"),
                    Topology::Linear,
                )
                .expect("test molecule"),
            )
        };
        session.open_document(doc("gfp start", base));
        session.open_document(doc("gfp sub", &base.replacen("AATTCTTG", "AATACTTG", 1)));
        let tabs = [(2u64, "gfp sub".to_string())];
        let mut state = MsaState::default();
        assert_eq!(state.count(), 1, "the document alone is not an alignment");
        state.documents.push(2);
        assert_eq!(state.count(), 2);
        assert_eq!(state.params()["documents"], json!([2]));
        let r = session
            .execute("align_sequences", state.params())
            .map_err(|e| e.message);
        let r = r.expect("two open documents align");
        let rows = r["rows"].as_array().expect("aligned rows");
        assert_eq!(rows.len(), 2);
        assert_eq!(rows[1]["name"], "gfp sub");
        // A tab that was closed since it was ticked is dropped, not sent as a missing id.
        state.result = Some(Ok(r));
        let t = Tokens::dark();
        crate::offscreen::render(
            900,
            600,
            1.0,
            |ctx| t.apply(ctx),
            |ui| {
                egui::Frame::new().fill(t.panel).show(ui, |ui| {
                    let _ = show(ui, &t, &mut state, "gfp start", &[], &[]);
                });
            },
        )
        .expect("the panel draws");
        assert!(state.documents.is_empty());
        assert_eq!(tabs.len(), 1);
    }
}
