//! The Align view: align a sequence (a sequencing read, an assembly, another construct) to the
//! open document through the `align_to_reference` command, list the differences, show them in
//! the sequence, annotate them, and read the aligned rows.

use eframe::egui::{self, RichText, Ui, text::LayoutJob};
use serde_json::{Value, json};
use snappiesnap_engine::Session;

use crate::theme::Tokens;

/// Where the sequence to align comes from.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum Source {
    #[default]
    File,
    /// Another document open in a tab: the only way in the web build, where the browser gives
    /// the page no file paths.
    Tab,
    Paste,
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum Mode {
    #[default]
    Fit,
    Local,
    Global,
}

impl Mode {
    fn id(self) -> &'static str {
        match self {
            Mode::Fit => "fit",
            Mode::Local => "local",
            Mode::Global => "global",
        }
    }
}

#[derive(Debug, Default)]
pub struct AlignState {
    /// A file to align, when `source` is [`Source::File`].
    pub file: String,
    pub pasted: String,
    /// A document open in a tab, when `source` is [`Source::Tab`].
    pub document: Option<u64>,
    pub source: Source,
    pub mode: Mode,
    /// Align a protein to the translation of this CDS feature (DNA documents).
    pub cds: Option<u64>,
    /// The last answer of `align_to_reference`, or its error.
    pub result: Option<Result<Value, String>>,
    /// Multiple alignment instead of aligning one sequence to this document.
    pub multiple: bool,
    pub msa: crate::msa_view::MsaState,
}

/// What the user asked for.
pub enum AlignAction {
    /// Run `align_to_reference` with these params.
    Align(Value),
    /// Run `align_sequences` with these params.
    AlignMany(Value),
    /// Let the user choose sequence files, wherever the build gets files from. `multiple` for
    /// the multiple alignment's list; otherwise one file for the pairwise alignment.
    AddFiles { multiple: bool },
    /// Select this range of the document and show it in the sequence view.
    Show { start: usize, end: usize },
    /// Add this feature (an `add_feature` params object).
    Annotate(Value),
    /// Open the sequence view's "Align with" panel (reads under the sequence, with traces).
    AlignWith,
}

impl AlignState {
    fn params(&self) -> Value {
        let mut p = json!({ "mode": self.mode.id(), "rows": true });
        if let Some(o) = p.as_object_mut() {
            match self.source {
                Source::Paste => o.insert("residues".into(), json!(self.pasted)),
                Source::Tab => o.insert("document".into(), json!(self.document)),
                Source::File => o.insert("file".into(), json!(self.file)),
            };
            if let Some(id) = self.cds {
                o.insert("feature".into(), json!(id));
            }
        }
        p
    }
}

fn num(v: &Value, k: &str) -> usize {
    v.get(k)
        .and_then(Value::as_u64)
        .and_then(|n| usize::try_from(n).ok())
        .unwrap_or(0)
}

fn text<'a>(v: &'a Value, k: &str) -> &'a str {
    v.get(k).and_then(Value::as_str).unwrap_or("")
}

/// A difference in words, its range in the document, and a feature name for it.
fn describe(d: &Value, len: usize) -> (String, String, usize, usize, String) {
    let wrap = |p: usize| if len == 0 { 0 } else { p % len };
    match text(d, "kind") {
        "substitution" => {
            let p = num(d, "position");
            let (r, q) = (
                text(d, "reference").to_uppercase(),
                text(d, "query").to_uppercase(),
            );
            (
                "Substitution".into(),
                format!("{r} \u{2192} {q}"),
                p,
                wrap(p + 1),
                format!("{r}{}{q} (alignment)", p + 1),
            )
        }
        "insertion" => {
            let p = num(d, "position");
            let bases = text(d, "bases").to_uppercase();
            // Show the bases either side of the insertion point.
            (
                "Insertion".into(),
                format!("+{bases}"),
                wrap(p + len.saturating_sub(1)),
                wrap(p + 1),
                format!(
                    "ins {bases} after {} (alignment)",
                    if p == 0 { len } else { p }
                ),
            )
        }
        _ => {
            let (s, n) = (num(d, "start"), num(d, "length"));
            let bases = text(d, "bases").to_uppercase();
            (
                "Deletion".into(),
                format!("\u{2212}{bases}"),
                s,
                wrap(s + n),
                format!("del {bases} at {} (alignment)", s + 1),
            )
        }
    }
}

pub fn show(
    ui: &mut Ui,
    t: &Tokens,
    session: &Session,
    state: &mut AlignState,
    files: &[String],
) -> Option<AlignAction> {
    let doc = session.document()?;
    // Every other open document, as a sequence to align.
    let active = session.active_document_id();
    let tabs: Vec<(u64, String)> = session
        .documents()
        .into_iter()
        .filter(|d| Some(d.id) != active)
        .map(|d| (d.id, d.name))
        .collect();
    ui.horizontal(|ui| {
        ui.selectable_value(&mut state.multiple, false, "To this sequence")
            .on_hover_text("Align one sequence (a read, an assembly, a protein) to this document");
        ui.selectable_value(&mut state.multiple, true, "Multiple")
            .on_hover_text("Align several sequences to each other");
    });
    if ui
        .button("Align reads under the sequence\u{2026}")
        .on_hover_text("Open the sequence view with the Align with panel: reads and sequences aligned under this one, with their traces")
        .clicked()
    {
        return Some(AlignAction::AlignWith);
    }
    ui.add_space(t.spacing);
    if state.multiple {
        return crate::msa_view::show(ui, t, &mut state.msa, doc.molecule.name(), files, &tabs)
            .map(|a| match a {
                crate::msa_view::MsaAction::Align(p) => AlignAction::AlignMany(p),
                crate::msa_view::MsaAction::AddFiles => AlignAction::AddFiles { multiple: true },
            });
    }
    let height = ui.available_height();
    let Some(left) = crate::two_columns(ui.available_width()) else {
        let action = pairwise(ui, t, session, state, files, &tabs);
        return action.or_else(|| pairwise_result(ui, t, session, state));
    };
    let mut action = None;
    ui.horizontal_top(|ui| {
        ui.allocate_ui_with_layout(
            egui::vec2(left, height),
            egui::Layout::top_down(egui::Align::Min),
            |ui| {
                ui.set_max_width(left);
                egui::ScrollArea::vertical()
                    .id_salt("align-inputs-column")
                    .max_width(left)
                    .show(ui, |ui| {
                        ui.set_max_width(left);
                        action = pairwise(ui, t, session, state, files, &tabs);
                    });
            },
        );
        ui.separator();
        egui::ScrollArea::vertical()
            .id_salt("align-result-column")
            .show(ui, |ui| {
                ui.vertical(|ui| {
                    if state.result.is_none() {
                        ui.label(RichText::new("The alignment appears here.").color(t.text_muted));
                    }
                    if let Some(a) = pairwise_result(ui, t, session, state) {
                        action = Some(a);
                    }
                });
            });
    });
    action
}

fn pairwise(
    ui: &mut Ui,
    t: &Tokens,
    session: &Session,
    state: &mut AlignState,
    files: &[String],
    tabs: &[(u64, String)],
) -> Option<AlignAction> {
    let doc = session.document()?;
    let m = &doc.molecule;
    let len = m.len();
    let protein_doc = m.seq().alphabet() == snappiesnap_engine::doc::Alphabet::Protein;
    let cds: Vec<(u64, String)> = if protein_doc {
        Vec::new()
    } else {
        m.features()
            .iter()
            .filter(|f| f.kind == "CDS")
            .map(|f| (f.id.0, f.name.clone()))
            .collect()
    };
    if state
        .cds
        .is_some_and(|id| !cds.iter().any(|(c, _)| *c == id))
    {
        state.cds = None;
    }
    let mut action = None;
    let intro = if protein_doc {
        format!("Align a protein to {} ({len} aa), with BLOSUM62.", m.name())
    } else {
        format!(
            "Align a sequence to {} ({len} bp): a sequencing read, an assembly or another construct (both strands are tried), or a protein to the translation of one of its CDS features.",
            m.name()
        )
    };
    ui.label(RichText::new(intro).color(t.text_muted));
    ui.add_space(t.spacing);
    egui::Grid::new("align-inputs")
        .num_columns(2)
        .spacing(egui::vec2(t.spacing * 2.0, t.spacing))
        .show(ui, |ui| {
            ui.label("Sequence");
            ui.horizontal(|ui| {
                ui.selectable_value(&mut state.source, Source::File, "From a file");
                ui.selectable_value(&mut state.source, Source::Tab, "From a tab")
                    .on_hover_text("Another construct or read already open");
                ui.selectable_value(&mut state.source, Source::Paste, "Paste");
            });
            ui.end_row();
            ui.label("");
            match state.source {
                Source::Paste => {
                    ui.add(
                        egui::TextEdit::multiline(&mut state.pasted)
                            .hint_text("ACGT\u{2026}")
                            .desired_rows(3)
                            .desired_width(520.0)
                            .font(egui::TextStyle::Monospace),
                    );
                }
                Source::Tab => {
                    if tabs.is_empty() {
                        ui.label(
                            RichText::new("No other document is open.").color(t.text_muted),
                        );
                    } else {
                        let shown = state
                            .document
                            .and_then(|id| tabs.iter().find(|(i, _)| *i == id))
                            .map_or("Choose\u{2026}", |(_, n)| n.as_str())
                            .to_string();
                        egui::ComboBox::from_id_salt("align-tab")
                            .width(320.0)
                            .selected_text(shown)
                            .show_ui(ui, |ui| {
                                for (id, name) in tabs {
                                    ui.selectable_value(&mut state.document, Some(*id), name);
                                }
                            });
                    }
                }
                Source::File => {
                    ui.horizontal(|ui| {
                        crate::sidebar::file_picker_pub(
                            ui,
                            t,
                            "align-file",
                            &mut state.file,
                            files,
                        );
                        if ui
                            .button("Choose\u{2026}")
                            .on_hover_text("Pick a sequence file to align")
                            .clicked()
                        {
                            action = Some(AlignAction::AddFiles { multiple: false });
                        }
                    });
                }
            }
            ui.end_row();
            if !protein_doc {
                ui.label("Against");
                let shown = match state.cds.and_then(|id| cds.iter().find(|(c, _)| *c == id)) {
                    Some((_, n)) => format!("Translation of CDS {n}"),
                    None => "This DNA".to_string(),
                };
                egui::ComboBox::from_id_salt("align-against")
                    .width(280.0)
                    .selected_text(shown)
                    .show_ui(ui, |ui| {
                        ui.selectable_value(&mut state.cds, None, "This DNA");
                        for (id, n) in &cds {
                            ui.selectable_value(
                                &mut state.cds,
                                Some(*id),
                                format!("Translation of CDS {n} (align a protein)"),
                            );
                        }
                    });
                ui.end_row();
            }
            ui.label("Mode");
            ui.horizontal(|ui| {
                for (mode, label, tip) in [
                    (Mode::Fit, "Fit", "the whole sequence, somewhere in this document (reads, inserts, assemblies)"),
                    (Mode::Local, "Local", "only the best matching stretch of each"),
                    (Mode::Global, "Global", "both end to end"),
                ] {
                    ui.selectable_value(&mut state.mode, mode, label)
                        .on_hover_text(tip);
                }
            });
            ui.end_row();
        });
    if state
        .document
        .is_some_and(|id| !tabs.iter().any(|(i, _)| *i == id))
    {
        state.document = None;
    }
    let ready = match state.source {
        Source::Paste => state.pasted.chars().any(|c| c.is_ascii_alphabetic()),
        Source::Tab => state.document.is_some(),
        Source::File => !state.file.is_empty(),
    };
    ui.add_space(t.spacing);
    if ui.add_enabled(ready, egui::Button::new("Align")).clicked() {
        action = Some(AlignAction::Align(state.params()));
    }
    ui.add_space(t.spacing);
    action
}

/// The alignment itself, after `pairwise` has run it.
fn pairwise_result(
    ui: &mut Ui,
    t: &Tokens,
    session: &Session,
    state: &AlignState,
) -> Option<AlignAction> {
    let doc = session.document()?;
    let m = &doc.molecule;
    let len = m.len();
    match &state.result {
        None => None,
        Some(Err(e)) => {
            ui.label(RichText::new(e).color(t.error));
            None
        }
        Some(Ok(r)) => {
            // Positions of a protein against a CDS map to that CDS's codons.
            let codons = r
                .get("kind")
                .and_then(Value::as_str)
                .filter(|k| *k == "protein")
                .and(state.cds)
                .and_then(|id| m.features().iter().find(|f| f.id.0 == id))
                .map(|f| cds_positions(m, f));
            let to_doc = |a: usize, b: usize| -> Option<(usize, usize)> {
                match &codons {
                    None => Some((a, b)),
                    Some(pos) => codon_range(pos, a, b, len),
                }
            };
            result(ui, t, r, &to_doc)
        }
    }
}

fn result(
    ui: &mut Ui,
    t: &Tokens,
    r: &Value,
    to_doc: &dyn Fn(usize, usize) -> Option<(usize, usize)>,
) -> Option<AlignAction> {
    let mut action = None;
    let protein = text(r, "kind") == "protein";
    let len = num(r, "reference_length");
    let units = if protein { "residue(s)" } else { "base(s)" };
    let identity = r.get("identity").and_then(Value::as_f64).unwrap_or(0.0) * 100.0;
    let (start, end, covered) = (num(r, "start"), num(r, "end"), num(r, "length"));
    let shown_end = if end == 0 { len } else { end };
    let summary = if protein {
        let similarity = r.get("similarity").and_then(Value::as_f64).unwrap_or(0.0) * 100.0;
        format!(
            "{} ({} aa) aligns to residues {}..{shown_end} of {} ({len} aa): {identity:.2}% identical, {similarity:.2}% similar",
            text(r, "query_name"),
            num(r, "query_length"),
            start + 1,
            text(r, "reference_name"),
        )
    } else {
        let strand = if text(r, "strand") == "reverse" {
            "reverse"
        } else {
            "forward"
        };
        format!(
            "{} ({} bp) aligns to {}..{shown_end} on the {strand} strand ({covered} bp): {identity:.2}% identical",
            text(r, "query_name"),
            num(r, "query_length"),
            start + 1,
        )
    };
    ui.label(RichText::new(summary).strong());
    let (qs, qe, ql) = (
        num(r, "query_start"),
        num(r, "query_end"),
        num(r, "query_length"),
    );
    if qs > 0 || qe < ql {
        ui.label(
            RichText::new(format!(
                "Only {} {}..{qe} of the sequence align.",
                if protein { "residues" } else { "bases" },
                qs + 1
            ))
            .color(t.text_muted),
        );
    }
    ui.label(format!(
        "{} substitution(s), {} inserted {units}, {} deleted {units}",
        num(r, "mismatches"),
        num(r, "insertions"),
        num(r, "deletions")
    ));
    ui.add_space(t.spacing);
    let empty = Vec::new();
    let diffs = r
        .get("differences")
        .and_then(Value::as_array)
        .unwrap_or(&empty);
    if r.get("differences_truncated").and_then(Value::as_bool) == Some(true) {
        ui.label(RichText::new("Only the first differences are listed.").color(t.text_muted));
    }
    if !diffs.is_empty() {
        egui::ScrollArea::vertical()
            .id_salt("align-diffs")
            .max_height(220.0)
            .show(ui, |ui| {
                egui::Grid::new("align-diff-table")
                    .num_columns(4)
                    .striped(true)
                    .spacing(egui::vec2(t.spacing * 3.0, t.spacing * 0.5))
                    .show(ui, |ui| {
                        for h in ["Position", "Difference", "Change", ""] {
                            ui.label(RichText::new(h).strong());
                        }
                        ui.end_row();
                        for d in diffs.iter().take(2000) {
                            let (kind, change, s, e, name) = describe(d, len);
                            let place = to_doc(s, e);
                            let at = match text(d, "kind") {
                                "deletion" => num(d, "start") + 1,
                                _ => num(d, "position") + 1,
                            };
                            ui.label(at.to_string());
                            ui.label(kind);
                            ui.label(RichText::new(change).monospace());
                            ui.horizontal(|ui| {
                                let Some((s, e)) = place else {
                                    return;
                                };
                                if ui
                                    .small_button("Show")
                                    .on_hover_text("Select it in the sequence")
                                    .clicked()
                                {
                                    action = Some(AlignAction::Show { start: s, end: e });
                                }
                                if ui
                                    .small_button("Annotate")
                                    .on_hover_text("Add a feature here noting the difference")
                                    .clicked()
                                {
                                    action = Some(AlignAction::Annotate(json!({
                                        "feature": {
                                            "name": name,
                                            "kind": "misc_difference",
                                            "segments": [{"range": {"start": s, "end": e}}],
                                        }
                                    })));
                                }
                            });
                            ui.end_row();
                        }
                    });
            });
        ui.add_space(t.spacing);
    }
    let (rr, qr) = (text(r, "reference_row"), text(r, "query_row"));
    if !rr.is_empty() {
        ui.label(RichText::new("Alignment").strong());
        rows(ui, t, rr, qr, start, len);
    }
    action
}

/// The document positions of a CDS's bases, in reading order (its strand), after
/// `/codon_start`.
fn cds_positions(
    m: &snappiesnap_engine::doc::Molecule,
    f: &snappiesnap_engine::doc::Feature,
) -> Vec<usize> {
    use snappiesnap_engine::doc::Strand;
    let len = m.len();
    let mut out = Vec::new();
    if len == 0 {
        return out;
    }
    for seg in &f.segments {
        let n = seg.range.len_on(len).min(len);
        let mut stretch: Vec<usize> = (0..n).map(|k| (seg.range.start + k) % len).collect();
        if f.strand == Strand::Reverse {
            stretch.reverse();
        }
        out.extend(stretch);
    }
    let skip = f
        .qualifiers
        .iter()
        .find(|q| q.name == "codon_start")
        .and_then(|q| q.values.first())
        .and_then(|v| v.trim().parse::<usize>().ok())
        .unwrap_or(1)
        .clamp(1, 3)
        - 1;
    out.drain(..skip.min(out.len()));
    out
}

/// The document range of residues `a..b` of a CDS (`positions` from [`cds_positions`]); an
/// insertion (`a == b`) shows the codons either side.
fn codon_range(positions: &[usize], a: usize, b: usize, len: usize) -> Option<(usize, usize)> {
    let (a, b) = if a == b {
        (a.saturating_sub(1), a + 1)
    } else {
        (a, b)
    };
    let first = *positions.get(a * 3)?;
    let last = *positions.get((b * 3).min(positions.len()).checked_sub(1)?)?;
    // Along the top strand: forward CDSs read up, reverse ones down.
    let (lo, hi) = if positions.len() > 1
        && positions
            .first()
            .zip(positions.get(1))
            .is_some_and(|(x, y)| (*y + len - *x) % len.max(1) == 1)
    {
        (first, last)
    } else {
        (last, first)
    };
    Some((lo, (hi + 1) % len.max(1)))
}

/// The aligned rows in blocks of 60 columns: document (with positions), marks, sequence.
fn rows(ui: &mut Ui, t: &Tokens, reference: &str, query: &str, start: usize, len: usize) {
    let font = egui::FontId::monospace(12.0);
    let (line, cw) = ui.fonts_mut(|f| (f.row_height(&font), f.glyph_width(&font, 'A').max(1.0)));
    // The position gutter takes about eight characters.
    let width = crate::alignment_columns(ui.available_width() - cw * 10.0, cw);
    let r: Vec<char> = reference.chars().collect();
    let q: Vec<char> = query.chars().collect();
    let blocks = r.len().div_ceil(width);
    // Document position at the start of each block.
    let mut starts = Vec::with_capacity(blocks);
    let mut pos = start;
    for b in 0..blocks {
        starts.push(pos);
        let n = r
            .get(b * width..((b + 1) * width).min(r.len()))
            .map_or(0, |s| s.iter().filter(|c| **c != '-').count());
        pos += n;
    }
    egui::ScrollArea::vertical()
        .id_salt("align-rows")
        .show_rows(ui, line * 4.0, blocks, |ui, range| {
            for b in range {
                let from = b * width;
                let to = ((b + 1) * width).min(r.len());
                let (Some(rs), Some(qs)) = (r.get(from..to), q.get(from..to.min(q.len()))) else {
                    continue;
                };
                let p = starts.get(b).copied().unwrap_or(0);
                let shown = if len == 0 { 0 } else { p % len } + 1;
                let mut job = LayoutJob::default();
                let fmt = |color| egui::TextFormat {
                    font_id: font.clone(),
                    color,
                    ..Default::default()
                };
                job.append(&format!("{shown:>8}  "), 0.0, fmt(t.text_muted));
                job.append(&rs.iter().collect::<String>(), 0.0, fmt(t.text));
                job.append("\n          ", 0.0, fmt(t.text));
                for (a, c) in rs.iter().zip(qs) {
                    let same = a.eq_ignore_ascii_case(c);
                    job.append(if same { "|" } else { " " }, 0.0, fmt(t.text_muted));
                }
                job.append("\n          ", 0.0, fmt(t.text));
                for (a, c) in rs.iter().zip(qs) {
                    let same = a.eq_ignore_ascii_case(c);
                    job.append(
                        &c.to_string(),
                        0.0,
                        fmt(if same { t.text } else { t.error }),
                    );
                }
                ui.label(job);
                ui.add_space(line * 0.6);
            }
        });
}

#[cfg(test)]
mod tests {
    use snappiesnap_engine::Document;
    use snappiesnap_engine::doc::{Alphabet, Molecule, Seq, Topology};

    use super::*;

    /// Each source makes its own params, and a tab aligns like a file.
    #[test]
    fn a_sequence_to_align_comes_from_a_file_a_tab_or_pasted_text() {
        let mut state = AlignState {
            file: "read.ab1".into(),
            pasted: "ACGT".into(),
            document: Some(3),
            ..Default::default()
        };
        assert_eq!(state.params()["file"], "read.ab1");
        state.source = Source::Tab;
        let p = state.params();
        assert_eq!(p["document"], 3);
        assert!(p.get("file").is_none() && p.get("residues").is_none());
        state.source = Source::Paste;
        assert_eq!(state.params()["residues"], "ACGT");

        // A tab aligns to the active document the same way a file does.
        let plasmid = "ATGCGTACGTTAGCCGATCGGATCCTTAGGCATGCAAGCTTGACCTGAAGTCCGATACGGTACCGAGCTCGAATTCACTGGCCGTCGTTTTACAACGTCGTGACTGGGAAAACCCTGGCG";
        let doc = |name: &str, residues: &str, topology| {
            Document::new(
                Molecule::new(
                    name,
                    Seq::new(Alphabet::Dna, residues).expect("test sequence"),
                    topology,
                )
                .expect("test molecule"),
            )
        };
        let mut session = Session::new();
        session.open_document(doc("pTest", plasmid, Topology::Circular));
        let part = plasmid.get(10..90).expect("part of the plasmid");
        session.open_document(doc("insert", part, Topology::Linear));
        // The second tab is active now, so align back to the plasmid from its own tab.
        let answer = session
            .execute(
                "align_to_reference",
                json!({"document": 1, "mode": "local"}),
            )
            .expect("an open tab aligns");
        assert_eq!(answer["identity"], 1.0);
        // An id that is not open is refused, not ignored.
        let e = session
            .execute("align_to_reference", json!({"document": 99}))
            .expect_err("an id that is not open");
        assert!(e.message.contains("99"), "{}", e.message);

        // The panel offering the open tabs, with no folder to list. Set
        // `SNAPPIESNAP_ALIGN_TAB_PNG=/path.png` to look at it.
        let mut state = AlignState {
            source: Source::Tab,
            document: Some(1),
            ..Default::default()
        };
        let t = Tokens::dark();
        let image = crate::offscreen::render(
            900,
            420,
            1.0,
            |ctx| t.apply(ctx),
            |ui| {
                egui::Frame::new().fill(t.panel).show(ui, |ui| {
                    let _ = show(ui, &t, &session, &mut state, &[]);
                });
            },
        )
        .expect("the panel draws");
        if let Ok(path) = std::env::var("SNAPPIESNAP_ALIGN_TAB_PNG") {
            crate::offscreen::write_png(&path, &image).expect("writing the png");
        }
        // The chosen tab is still open, so it stays chosen.
        assert_eq!(state.document, Some(1));
    }

    /// A read with one substitution, one insertion and one deletion, aligned and drawn. Set
    /// `SNAPPIESNAP_ALIGN_PNG=/path.png` to look at it.
    #[test]
    fn alignment_draws_with_its_differences() {
        let plasmid = "ATGCGTACGTTAGCCGATCGGATCCTTAGGCATGCAAGCTTGACCTGAAGTCCGATACGGTACCGAGCTCGAATTCACTGGCCGTCGTTTTACAACGTCGTGACTGGGAAAACCCTGGCG";
        let mut session = Session::new();
        session.open_document(Document::new(
            Molecule::new(
                "pTest",
                Seq::new(Alphabet::Dna, plasmid).unwrap(),
                Topology::Circular,
            )
            .unwrap(),
        ));
        let mut read = plasmid.as_bytes()[10..110].to_vec();
        read[20] = b'A'; // plasmid 30 is C
        read.insert(50, b'T');
        read.remove(80);
        let mut state = AlignState {
            source: Source::Paste,
            pasted: String::from_utf8(read).unwrap(),
            ..Default::default()
        };
        let answer = session
            .execute("align_to_reference", state.params())
            .map_err(|e| e.message);
        let r = answer.clone().unwrap();
        assert_eq!(
            (
                num(&r, "mismatches"),
                num(&r, "insertions"),
                num(&r, "deletions")
            ),
            (1, 1, 1)
        );
        let diffs = r["differences"].as_array().unwrap();
        let (kind, change, s, e, name) = describe(&diffs[0], plasmid.len());
        assert_eq!(
            (kind.as_str(), change.as_str(), s, e),
            ("Substitution", "C \u{2192} A", 30, 31)
        );
        assert_eq!(name, "C31A (alignment)");
        state.result = Some(answer);
        let t = Tokens::light();
        let mut asked = None;
        let image = crate::offscreen::render(
            1400,
            700,
            1.0,
            |ctx| t.apply(ctx),
            |ui| {
                egui::Frame::new().fill(t.panel).show(ui, |ui| {
                    asked = Some(show(ui, &t, &session, &mut state, &[]).is_some());
                });
            },
        )
        .unwrap();
        assert_eq!(asked, Some(false));
        if let Ok(path) = std::env::var("SNAPPIESNAP_ALIGN_PNG") {
            crate::offscreen::write_png(&path, &image).unwrap();
        }
    }

    /// A protein against the translation of a CDS: its differences map to codons, on either
    /// strand.
    #[test]
    fn protein_differences_map_to_codons() {
        use snappiesnap_engine::doc::{NewFeature, Range, Segment, Strand};
        // 6 bases, then ATG GCC AAA TGG TAA (M A K W *), then 6 bases.
        let dna = "GGGGGGATGGCCAAATGGTAAGGGGGG";
        for strand in [Strand::Forward, Strand::Reverse] {
            let seq = match strand {
                Strand::Forward => dna.to_string(),
                _ => dna
                    .bytes()
                    .rev()
                    .map(|b| match b {
                        b'A' => 'T',
                        b'T' => 'A',
                        b'C' => 'G',
                        _ => 'C',
                    })
                    .collect(),
            };
            let mut m = Molecule::new(
                "p",
                Seq::new(Alphabet::Dna, seq.as_str()).unwrap(),
                Topology::Linear,
            )
            .unwrap();
            let range = match strand {
                Strand::Forward => Range::new(6, 21),
                _ => Range::new(6, 21),
            };
            let id = m
                .add_feature(NewFeature {
                    color: None,
                    name: "gene".into(),
                    kind: "CDS".into(),
                    strand,
                    segments: vec![Segment { range }],
                    qualifiers: Vec::new(),
                })
                .unwrap();
            let mut session = Session::new();
            session.open_document(Document::new(m));
            // K (residue 3) changed to R.
            let state = AlignState {
                source: Source::Paste,
                pasted: "MARW".into(),
                cds: Some(id.0),
                mode: Mode::Global,
                ..Default::default()
            };
            let r = session
                .execute("align_to_reference", state.params())
                .unwrap();
            assert_eq!(r["kind"], "protein");
            assert_eq!(r["reference_length"], 4);
            let d = &r["differences"][0];
            assert_eq!(
                (d["position"].as_u64(), d["reference"].as_str()),
                (Some(2), Some("K"))
            );
            let (_, _, s, e, name) = describe(d, 4);
            assert_eq!(name, "K3R (alignment)");
            let doc = session.document().unwrap();
            let f = doc.molecule.features().first().unwrap();
            let pos = cds_positions(&doc.molecule, f);
            // The AAA codon: bases 12..15 forward; on the reverse strand its image, 12..15 too
            // (the CDS sits in the middle of a 27 bp molecule).
            assert_eq!(
                codon_range(&pos, s, e, doc.molecule.len()),
                Some((12, 15)),
                "{strand:?}"
            );
        }
    }
}
