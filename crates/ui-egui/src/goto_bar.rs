//! The Go to bar of the map and sequence views (Cmd/Ctrl+L): type a position or range as
//! people write them ("1200", "100-250", "4000-200" across the origin of a circle) and it is
//! selected. Parsing and selecting are the `go_to` command's, so agents get the same result.

use eframe::egui::{self, Key, RichText, Ui};
use serde_json::{Value, json};
use snappiesnap_engine::Session;
use snappiesnap_engine::doc::Range;

use crate::theme::Tokens;

#[derive(Debug, Default, Clone)]
pub struct GoToState {
    pub open: bool,
    pub text: String,
    pub error: Option<String>,
    /// Put the cursor in the text field on the next frame.
    focus: bool,
}

impl GoToState {
    /// Open the bar (or focus it if it's open), with its text selected for typing over.
    pub fn open(&mut self) {
        self.open = true;
        self.focus = true;
        self.error = None;
    }
}

/// The bar, if open. Returns the range selected when the user pressed Enter or Go.
pub fn show(ui: &mut Ui, t: &Tokens, session: &mut Session, st: &mut GoToState) -> Option<Range> {
    if !st.open {
        return None;
    }
    let mut selected = None;
    ui.horizontal(|ui| {
        ui.label(RichText::new("Go to").color(t.text_muted));
        let r = ui.add(
            egui::TextEdit::singleline(&mut st.text)
                .hint_text("1200, or 100-250")
                .desired_width(180.0),
        );
        if std::mem::take(&mut st.focus) {
            r.request_focus();
        }
        let enter = r.lost_focus() && ui.input(|i| i.key_pressed(Key::Enter));
        let go = ui.button("Go").clicked();
        if (enter || go) && !st.text.trim().is_empty() {
            match session.execute("go_to", json!({ "location": st.text })) {
                Ok(v) => {
                    selected = v
                        .get("selection")
                        .and_then(|r| serde_json::from_value::<Range>(r.clone()).ok());
                    st.error = None;
                    st.open = false;
                }
                Err(e) => {
                    st.error = Some(e.message);
                    st.focus = true;
                }
            }
        }
        let len = session.document().map_or(0, |d| d.molecule.len());
        match &st.error {
            Some(e) => ui.label(RichText::new(e).color(t.error)),
            None => ui.label(
                RichText::new(format!("1 to {len}; a range selects it")).color(t.text_muted),
            ),
        };
        if ui
            .small_button("\u{d7}")
            .on_hover_text("Close (Esc)")
            .clicked()
            || (r.has_focus() && ui.input(|i| i.key_pressed(Key::Escape)))
        {
            st.open = false;
        }
    });
    selected
}

/// The status bar's summary of a `sequence_info` result: GC %, Tm and weight.
pub fn summary(info: &Value) -> Vec<String> {
    let mut parts = Vec::new();
    if let Some(gc) = info.get("gc_percent").and_then(Value::as_f64) {
        parts.push(format!("GC {gc:.1}%"));
    }
    if let Some(tm) = info.get("tm").and_then(Value::as_f64) {
        parts.push(format!("Tm {tm:.1} \u{b0}C"));
    }
    if let Some(w) = info.get("weight") {
        let single = info.get("whole").and_then(Value::as_bool) != Some(true)
            && info
                .get("length")
                .and_then(Value::as_u64)
                .is_some_and(|n| n <= 200);
        let (key, label) = if single {
            ("single_strand", "ss")
        } else {
            ("double_strand", "ds")
        };
        if let Some(da) = w.get(key).and_then(Value::as_f64) {
            let approx = if w.get("approximate").and_then(Value::as_bool) == Some(true) {
                "\u{2248}"
            } else {
                ""
            };
            parts.push(format!("{approx}{} ({label})", mass(da)));
        }
    }
    if let Some(da) = info.get("mass").and_then(Value::as_f64) {
        parts.push(mass(da));
    }
    parts
}

/// A mass in Da, kDa or MDa.
fn mass(da: f64) -> String {
    if da >= 1e6 {
        format!("{:.2} MDa", da / 1e6)
    } else if da >= 1e3 {
        format!("{:.1} kDa", da / 1e3)
    } else {
        format!("{da:.0} Da")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn summaries() {
        let oligo = json!({"length": 20, "whole": false, "gc_percent": 50.0, "tm": 61.24,
                           "weight": {"single_strand": 6113.0, "double_strand": 12300.0, "approximate": false}});
        assert_eq!(
            summary(&oligo),
            ["GC 50.0%", "Tm 61.2 \u{b0}C", "6.1 kDa (ss)"]
        );
        let plasmid = json!({"length": 4361, "whole": true, "gc_percent": 53.75,
                             "weight": {"single_strand": 1.3e6, "double_strand": 2_694_000.0, "approximate": true}});
        assert_eq!(summary(&plasmid), ["GC 53.8%", "\u{2248}2.69 MDa (ds)"]);
        assert_eq!(summary(&json!({"mass": 512.6})), ["513 Da"]);
        assert!(summary(&json!({})).is_empty());
    }
}
