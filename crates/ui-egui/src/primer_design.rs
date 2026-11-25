//! The primer design panel in the Primers tab: design a pair for the selected region (through
//! the `design_primers` command), review the candidates and keep one.

use eframe::egui::{self, RichText, Sense, Ui, vec2};
use serde_json::{Value, json};
use snappiesnap_engine::Session;
use snappiesnap_engine::doc::Range;

use crate::theme::Tokens;

/// What the user asked for in the panel.
#[derive(Debug, Clone, PartialEq)]
pub enum DesignAction {
    /// Run `design_primers` with these params.
    Design(Value),
    /// Keep a pair: (name, residues) for the forward and reverse primer.
    Add([(String, String); 2]),
    /// Select a range of the molecule.
    Select(Range),
    /// Show the PCR with these primers (names) on the gel.
    PcrGel(String, String),
    /// Open the PCR product as a new document.
    PcrOpen(String, String),
}

/// A designed primer as shown.
#[derive(Debug, Clone, PartialEq)]
pub struct ShownPrimer {
    pub residues: String,
    pub range: Range,
    pub tm: f64,
    pub gc: f64,
    pub warnings: Vec<String>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct ShownPair {
    pub forward: ShownPrimer,
    pub reverse: ShownPrimer,
    pub product: u64,
    pub warnings: Vec<String>,
}

pub struct DesignState {
    pub target_tm: f64,
    pub min_len: usize,
    pub max_len: usize,
    pub search: usize,
    /// The region the results are for, and the results (or the error).
    pub results: Option<(Range, Result<Vec<ShownPair>, String>)>,
    /// The primer pair chosen for a PCR (names).
    pub pcr: (String, String),
}

impl Default for DesignState {
    fn default() -> Self {
        DesignState {
            target_tm: 60.0,
            min_len: 18,
            max_len: 30,
            search: 0,
            results: None,
            pcr: (String::new(), String::new()),
        }
    }
}

impl DesignState {
    /// The `design_primers` params for `region`.
    pub fn params(&self, region: Range) -> Value {
        json!({
            "range": region,
            "target_tm": self.target_tm,
            "min_len": self.min_len,
            "max_len": self.max_len.max(self.min_len),
            "search": self.search,
            "count": 5,
        })
    }

    /// Keep the command's answer for `region`.
    pub fn set_results(&mut self, region: Range, answer: Result<Value, String>) {
        let parsed = answer.map(|v| {
            v.get("pairs")
                .and_then(Value::as_array)
                .map(|pairs| pairs.iter().filter_map(parse_pair).collect())
                .unwrap_or_default()
        });
        self.results = Some((region, parsed));
    }
}

fn parse_primer(v: &Value) -> Option<ShownPrimer> {
    let range: Range = serde_json::from_value(v.get("range")?.clone()).ok()?;
    Some(ShownPrimer {
        residues: v.get("residues")?.as_str()?.to_string(),
        range,
        tm: v.get("tm")?.as_f64()?,
        gc: v.get("gc_percent")?.as_f64()?,
        warnings: strings(v.get("warnings")),
    })
}

fn parse_pair(v: &Value) -> Option<ShownPair> {
    Some(ShownPair {
        forward: parse_primer(v.get("forward")?)?,
        reverse: parse_primer(v.get("reverse")?)?,
        product: v.get("product")?.as_u64()?,
        warnings: strings(v.get("warnings")),
    })
}

fn strings(v: Option<&Value>) -> Vec<String> {
    v.and_then(Value::as_array)
        .map(|a| {
            a.iter()
                .filter_map(|s| s.as_str().map(str::to_string))
                .collect()
        })
        .unwrap_or_default()
}

/// 1-based inclusive end of a range.
fn end1(r: Range, len: usize) -> usize {
    if r.end == 0 { len } else { r.end }
}

/// The "Simulate PCR" row: pick two primers of the molecule, then show the product on the gel or
/// open it.
pub fn pcr_row(
    ui: &mut Ui,
    t: &Tokens,
    session: &Session,
    state: &mut DesignState,
) -> Option<DesignAction> {
    let doc = session.document()?;
    let primers: Vec<String> = doc
        .molecule
        .primers()
        .iter()
        .map(|p| p.name.clone())
        .collect();
    if primers.is_empty() {
        return None;
    }
    let mut action = None;
    ui.horizontal(|ui| {
        ui.label(RichText::new("Simulate PCR").strong());
        for (k, (id, value)) in [("pcr-f", &mut state.pcr.0), ("pcr-r", &mut state.pcr.1)]
            .into_iter()
            .enumerate()
        {
            if !primers.contains(value) {
                // Forward defaults to the first primer, reverse to the second.
                *value = primers
                    .get(k)
                    .or(primers.first())
                    .cloned()
                    .unwrap_or_default();
            }
            egui::ComboBox::from_id_salt(id)
                .selected_text(value.as_str())
                .show_ui(ui, |ui| {
                    for p in &primers {
                        ui.selectable_value(value, p.clone(), p);
                    }
                });
        }
        let (f, r) = state.pcr.clone();
        if ui.button("Show on gel").clicked() {
            action = Some(DesignAction::PcrGel(f.clone(), r.clone()));
        }
        if ui.button("Open product").clicked() {
            action = Some(DesignAction::PcrOpen(f, r));
        }
        ui.label(
            RichText::new("Either primer may prime on either strand.")
                .small()
                .color(t.text_muted),
        );
    });
    action
}

pub fn show(
    ui: &mut Ui,
    t: &Tokens,
    session: &Session,
    state: &mut DesignState,
) -> Option<DesignAction> {
    let mut action = None;
    let doc = session.document()?;
    let len = doc.molecule.len();
    let selection = doc.selection.filter(|s| !s.is_empty());
    ui.label(RichText::new("Design primers").strong());
    ui.horizontal_wrapped(|ui| {
        match selection {
            Some(r) => {
                ui.label(format!(
                    "Amplify {}..{} ({} bp)",
                    r.start + 1,
                    end1(r, len),
                    r.len_on(len)
                ));
            }
            None => {
                ui.label(
                    RichText::new("Select the region to amplify on the map or in the sequence.")
                        .color(t.text_muted),
                );
            }
        }
        ui.separator();
        ui.label("Tm");
        ui.add(
            egui::DragValue::new(&mut state.target_tm)
                .range(40.0..=80.0)
                .speed(0.5)
                .suffix(" \u{b0}C"),
        );
        ui.label("length");
        ui.add(egui::DragValue::new(&mut state.min_len).range(10..=60));
        ui.label("to");
        ui.add(egui::DragValue::new(&mut state.max_len).range(10..=60));
        ui.label("outside by up to");
        ui.add(
            egui::DragValue::new(&mut state.search)
                .range(0..=500)
                .suffix(" bp"),
        );
        let design = ui.add_enabled(selection.is_some(), egui::Button::new("Design"));
        if design.clicked()
            && let Some(r) = selection
        {
            action = Some(DesignAction::Design(state.params(r)));
        }
    });
    let Some((region, results)) = &state.results else {
        return action;
    };
    ui.add_space(t.spacing);
    let pairs = match results {
        Ok(p) if !p.is_empty() => p,
        Ok(_) => {
            ui.label(RichText::new("No primer pairs found.").color(t.text_muted));
            return action;
        }
        Err(e) => {
            ui.label(RichText::new(e).color(t.error));
            return action;
        }
    };
    ui.label(
        RichText::new(format!(
            "Best pairs for {}..{}, best first. Click a primer to select where it binds.",
            region.start + 1,
            end1(*region, len)
        ))
        .small()
        .color(t.text_muted),
    );
    egui::Grid::new("designed")
        .striped(true)
        .spacing(vec2(t.spacing * 2.0, t.spacing * 0.6))
        .show(ui, |ui| {
            for h in [
                "",
                "Forward (5' to 3')",
                "Tm",
                "GC",
                "Reverse (5' to 3')",
                "Tm",
                "GC",
                "Product",
                "",
            ] {
                ui.label(RichText::new(h).strong());
            }
            ui.end_row();
            for (i, p) in pairs.iter().enumerate() {
                let mut notes: Vec<String> = Vec::new();
                notes.extend(p.forward.warnings.iter().map(|w| format!("forward: {w}")));
                notes.extend(p.reverse.warnings.iter().map(|w| format!("reverse: {w}")));
                notes.extend(p.warnings.iter().cloned());
                let marker = if notes.is_empty() {
                    RichText::new(format!("{}", i + 1))
                } else {
                    RichText::new(format!("{} !", i + 1)).color(t.error)
                };
                let r = ui.label(marker);
                if !notes.is_empty() {
                    r.on_hover_text(notes.join("\n"));
                }
                for primer in [&p.forward, &p.reverse] {
                    let seq = ui.add(
                        egui::Label::new(RichText::new(&primer.residues).monospace())
                            .sense(Sense::click()),
                    );
                    let hover = format!(
                        "{} nt, binds {}..{}{}",
                        primer.residues.len(),
                        primer.range.start + 1,
                        end1(primer.range, len),
                        if primer.warnings.is_empty() {
                            String::new()
                        } else {
                            format!("\n{}", primer.warnings.join("\n"))
                        }
                    );
                    if seq.on_hover_text(hover).clicked() {
                        action = Some(DesignAction::Select(primer.range));
                    }
                    ui.label(format!("{:.1} \u{b0}C", primer.tm));
                    ui.label(format!("{:.0} %", primer.gc));
                }
                ui.label(format!("{} bp", p.product));
                if ui.button("Add pair").clicked() {
                    action = Some(DesignAction::Add([
                        (
                            format!("F{}", p.forward.range.start + 1),
                            p.forward.residues.clone(),
                        ),
                        (
                            format!("R{}", end1(p.reverse.range, len)),
                            p.reverse.residues.clone(),
                        ),
                    ]));
                }
                ui.end_row();
            }
        });
    action
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_command_answers() {
        let mut s = DesignState::default();
        let answer = json!({"pairs": [{
            "forward": {"residues": "ACGT", "range": {"start": 1, "end": 5}, "tm": 60.1, "gc_percent": 50.0, "warnings": []},
            "reverse": {"residues": "TTGG", "range": {"start": 9, "end": 13}, "tm": 59.0, "gc_percent": 50.0, "warnings": ["x"]},
            "product": 12, "penalty": 1.0, "warnings": []
        }, {"broken": true}]});
        s.set_results(Range::new(1, 13), Ok(answer));
        let (_, r) = s.results.as_ref().unwrap();
        let pairs = r.as_ref().unwrap();
        assert_eq!(pairs.len(), 1);
        assert_eq!(pairs[0].reverse.warnings, vec!["x".to_string()]);
        assert_eq!(s.params(Range::new(1, 13))["max_len"], 30);
    }

    /// The panel with real results draws. Set `SNAPPIESNAP_DESIGN_PNG=/path.png` to look.
    #[test]
    fn panel_with_results_draws() {
        use snappiesnap_engine::Document;
        use snappiesnap_engine::doc::{Alphabet, Molecule, Seq, Topology};
        let mut x: u32 = 99;
        let seq: String = (0..600)
            .map(|_| {
                x = x.wrapping_mul(1_103_515_245).wrapping_add(12345);
                b"ACGT"[((x >> 16) % 4) as usize] as char
            })
            .collect();
        let m = Molecule::new(
            "m",
            Seq::new(Alphabet::Dna, seq).unwrap(),
            Topology::Circular,
        )
        .unwrap();
        let mut session = Session::new();
        session.open_document(Document::new(m));
        let region = Range::new(100, 500);
        session
            .execute("select", json!({ "range": region }))
            .unwrap();
        let mut state = DesignState {
            search: 30,
            ..DesignState::default()
        };
        let answer = session
            .execute("design_primers", state.params(region))
            .map_err(|e| e.message);
        state.set_results(region, answer);
        assert!(matches!(&state.results, Some((_, Ok(p))) if !p.is_empty()));
        let t = Tokens::light();
        let image = crate::offscreen::render(
            1100,
            300,
            1.0,
            |ctx| t.apply(ctx),
            |ui| {
                egui::Frame::new().fill(t.panel).show(ui, |ui| {
                    show(ui, &t, &session, &mut state);
                });
            },
        )
        .unwrap();
        if let Ok(path) = std::env::var("SNAPPIESNAP_DESIGN_PNG") {
            crate::offscreen::write_png(&path, &image).unwrap();
        }
    }
}
