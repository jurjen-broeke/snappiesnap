//! The Features tab: a sortable, searchable table of the molecule's features, and an editor for
//! a feature's name, type, direction, parts and qualifiers. Changes go through the
//! `add_feature`, `update_feature` and `remove_feature` commands (undoable, saved with the file).

use eframe::egui::{self, CornerRadius, RichText, Sense, Ui, vec2};
use serde_json::{Value, json};
use snappiesnap_engine::Session;
use snappiesnap_engine::doc::{Feature, FeatureId, Range, Strand};

use crate::theme::Tokens;

/// Common INSDC feature types, offered in the editor (any other type can be typed).
pub const KINDS: &[&str] = &[
    "CDS",
    "gene",
    "promoter",
    "terminator",
    "rep_origin",
    "enhancer",
    "RBS",
    "polyA_signal",
    "LTR",
    "intron",
    "exon",
    "5'UTR",
    "3'UTR",
    "mRNA",
    "ncRNA",
    "misc_RNA",
    "sig_peptide",
    "mat_peptide",
    "protein_bind",
    "primer_bind",
    "misc_binding",
    "misc_recomb",
    "mobile_element",
    "oriT",
    "repeat_region",
    "regulatory",
    "stem_loop",
    "variation",
    "misc_difference",
    "misc_feature",
];

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum Column {
    #[default]
    Location,
    Name,
    Size,
    Direction,
    Kind,
}

/// The feature being edited (a new one when `id` is `None`). Positions are 1-based, as shown.
#[derive(Debug, Clone, Default)]
pub struct Editor {
    pub id: Option<FeatureId>,
    pub name: String,
    pub kind: String,
    pub strand: Strand,
    /// Parts as (start, end), 1-based inclusive; end < start goes across the origin.
    pub parts: Vec<(String, String)>,
    pub qualifiers: Vec<(String, String)>,
    /// The feature's own colour as `#rrggbb` (`None`: automatic, by name and type).
    pub color: Option<String>,
    pub error: Option<String>,
}

#[derive(Debug, Default)]
pub struct FeaturesState {
    pub search: String,
    pub sort: Column,
    pub descending: bool,
    pub editor: Option<Editor>,
    /// What the feature library holds (or why there is none), once asked.
    pub library: Option<Result<Value, String>>,
    /// Features found by detection, each with whether it is ticked to be added.
    pub detected: Option<Result<Vec<(Value, bool)>, String>>,
    /// The folder last used to build the library.
    pub library_folder: String,
    /// A feature of the user's own design, being typed: name, type, sequence.
    pub design: (String, String, String),
    /// The terms of the vector-library download, while the user is asked to accept them.
    pub download_terms: Option<String>,
}

/// What the tab asks the app to run.
pub enum FeaturesAction {
    Run(&'static str, Value),
    /// Choose a folder (native dialog) and build the library from it.
    PickLibraryFolder,
    /// Choose a GFF3 file (native dialog) and add its features to the document.
    PickAnnotations,
}

/// The 1-based (start, end) of a 0-based half-open range on a molecule of `len`.
fn shown(r: Range, len: usize) -> (usize, usize) {
    (r.start + 1, if r.end == 0 { len } else { r.end })
}

/// The stretch a feature covers along the top strand, from its first part to its last.
fn span(f: &Feature) -> Option<Range> {
    let (first, last) = (f.segments.first()?, f.segments.last()?);
    Some(if f.strand == Strand::Reverse {
        Range::new(last.range.start, first.range.end)
    } else {
        Range::new(first.range.start, last.range.end)
    })
}

fn direction(s: Strand) -> &'static str {
    match s {
        Strand::Forward => "forward",
        Strand::Reverse => "reverse",
        Strand::Both => "both",
        Strand::None => "none",
    }
}

impl Editor {
    /// An editor filled from a feature.
    pub fn of(f: &Feature, len: usize) -> Editor {
        Editor {
            id: Some(f.id),
            name: f.name.clone(),
            kind: f.kind.clone(),
            strand: f.strand,
            parts: f
                .segments
                .iter()
                .map(|s| {
                    let (a, b) = shown(s.range, len);
                    (a.to_string(), b.to_string())
                })
                .collect(),
            qualifiers: f
                .qualifiers
                .iter()
                .flat_map(|q| q.values.iter().map(move |v| (q.name.clone(), v.clone())))
                .collect(),
            color: f.color.clone(),
            error: None,
        }
    }

    /// A new feature over `range`.
    pub fn new_over(range: Range, len: usize) -> Editor {
        let (a, b) = shown(range, len);
        Editor {
            id: None,
            name: "New feature".into(),
            kind: "misc_feature".into(),
            strand: Strand::Forward,
            parts: vec![(a.to_string(), b.to_string())],
            qualifiers: Vec::new(),
            color: None,
            error: None,
        }
    }

    /// The `add_feature` / `update_feature` params, starting from the feature as it is (so
    /// what the editor doesn't show is kept). Errors are in words.
    fn params(
        &self,
        existing: Option<&Feature>,
        len: usize,
    ) -> Result<(&'static str, Value), String> {
        let mut segments = Vec::new();
        for (i, (a, b)) in self.parts.iter().enumerate() {
            let parse = |t: &str| {
                t.trim()
                    .parse::<usize>()
                    .ok()
                    .filter(|&n| n >= 1 && n <= len)
            };
            let (Some(a), Some(b)) = (parse(a), parse(b)) else {
                return Err(format!(
                    "part {}: positions must be whole numbers from 1 to {len}",
                    i + 1
                ));
            };
            // 1-based inclusive to 0-based half-open; b < a goes across the origin.
            let end = if b == len { 0 } else { b };
            let start = a - 1;
            let end = if end == 0 && start == 0 { len } else { end };
            segments.push(json!({"range": {"start": start, "end": end}}));
        }
        if segments.is_empty() {
            return Err("a feature needs at least one part".into());
        }
        let mut qualifiers: Vec<(String, Vec<String>)> = Vec::new();
        for (name, value) in &self.qualifiers {
            let name = name.trim();
            if name.is_empty() {
                continue;
            }
            match qualifiers.iter_mut().find(|(n, _)| n == name) {
                Some((_, values)) => values.push(value.clone()),
                None => qualifiers.push((name.to_string(), vec![value.clone()])),
            }
        }
        let strand = match self.strand {
            Strand::Forward => "forward",
            Strand::Reverse => "reverse",
            Strand::Both => "both",
            Strand::None => "none",
        };
        let mut feature = existing
            .and_then(|f| serde_json::to_value(f).ok())
            .unwrap_or_else(|| json!({}));
        if let Some(o) = feature.as_object_mut() {
            o.insert("name".into(), json!(self.name));
            o.insert("kind".into(), json!(self.kind.trim()));
            o.insert("strand".into(), json!(strand));
            match &self.color {
                Some(c) => o.insert("color".into(), json!(c)),
                None => o.remove("color"),
            };
            o.insert("segments".into(), Value::Array(segments));
            o.insert(
                "qualifiers".into(),
                Value::Array(
                    qualifiers
                        .into_iter()
                        .map(|(n, v)| json!({"name": n, "values": v}))
                        .collect(),
                ),
            );
        }
        Ok(match self.id {
            Some(_) => ("update_feature", json!({ "feature": feature })),
            None => ("add_feature", json!({ "feature": feature })),
        })
    }
}

fn editor_ui(ui: &mut Ui, t: &Tokens, e: &mut Editor) -> (bool, bool) {
    let (mut save, mut cancel) = (false, false);
    egui::Frame::new()
        .stroke(egui::Stroke::new(1.0, t.border))
        .corner_radius(t.radius)
        .inner_margin(8)
        .show(ui, |ui| {
            ui.label(
                RichText::new(if e.id.is_some() {
                    "Edit feature"
                } else {
                    "New feature"
                })
                .strong(),
            );
            egui::Grid::new("feature-editor")
                .num_columns(2)
                .spacing(vec2(t.spacing * 2.0, t.spacing))
                .show(ui, |ui| {
                    ui.label("Name");
                    ui.add(egui::TextEdit::singleline(&mut e.name).desired_width(260.0));
                    ui.end_row();
                    ui.label("Type");
                    ui.horizontal(|ui| {
                        ui.add(egui::TextEdit::singleline(&mut e.kind).desired_width(160.0));
                        egui::ComboBox::from_id_salt("feature-kind")
                            .width(140.0)
                            .selected_text("common types")
                            .show_ui(ui, |ui| {
                                for k in KINDS {
                                    if ui.selectable_label(e.kind == *k, *k).clicked() {
                                        e.kind = (*k).to_string();
                                    }
                                }
                            });
                    });
                    ui.end_row();
                    ui.label("Colour");
                    ui.horizontal(|ui| {
                        let automatic = t.automatic_color(e.name.trim(), e.kind.trim());
                        let current = e.color.as_deref().and_then(crate::theme::parse_hex);
                        match crate::color_picker::color_picker(
                            ui,
                            t,
                            "feature-colour",
                            current,
                            automatic,
                        ) {
                            Some(crate::color_picker::ColorChoice::Set(c)) => {
                                e.color = Some(crate::theme::hex(c));
                            }
                            Some(crate::color_picker::ColorChoice::Default) => e.color = None,
                            None => {}
                        }
                        let note = if e.color.is_some() {
                            "its own colour, kept wherever the feature goes"
                        } else {
                            "automatic, from its name and type"
                        };
                        ui.label(RichText::new(note).small().color(t.text_muted));
                    });
                    ui.end_row();
                    ui.label("Direction");
                    ui.horizontal(|ui| {
                        for (s, label) in [
                            (Strand::Forward, "forward"),
                            (Strand::Reverse, "reverse"),
                            (Strand::Both, "both"),
                            (Strand::None, "none"),
                        ] {
                            ui.selectable_value(&mut e.strand, s, label);
                        }
                    });
                    ui.end_row();
                    ui.label("Parts");
                    ui.vertical(|ui| {
                        let mut remove = None;
                        let n = e.parts.len();
                        for (i, (a, b)) in e.parts.iter_mut().enumerate() {
                            ui.horizontal(|ui| {
                                ui.add(egui::TextEdit::singleline(a).desired_width(80.0));
                                ui.label("..");
                                ui.add(egui::TextEdit::singleline(b).desired_width(80.0));
                                if n > 1 && ui.small_button("Remove").clicked() {
                                    remove = Some(i);
                                }
                            });
                        }
                        if let Some(i) = remove
                            && i < e.parts.len()
                        {
                            e.parts.remove(i);
                        }
                        if ui
                            .small_button("Add a part")
                            .on_hover_text("For a feature in pieces (a join), listed 5' to 3' along it")
                            .clicked()
                        {
                            e.parts.push((String::new(), String::new()));
                        }
                    });
                    ui.end_row();
                    ui.label("Qualifiers");
                    ui.vertical(|ui| {
                        let mut remove = None;
                        for (i, (name, value)) in e.qualifiers.iter_mut().enumerate() {
                            ui.horizontal(|ui| {
                                ui.add(
                                    egui::TextEdit::singleline(name)
                                        .hint_text("note")
                                        .desired_width(110.0),
                                );
                                ui.add(egui::TextEdit::singleline(value).desired_width(320.0));
                                if ui.small_button("Remove").clicked() {
                                    remove = Some(i);
                                }
                            });
                        }
                        if let Some(i) = remove
                            && i < e.qualifiers.len()
                        {
                            e.qualifiers.remove(i);
                        }
                        if ui.small_button("Add a qualifier").clicked() {
                            e.qualifiers.push(("note".into(), String::new()));
                        }
                    });
                    ui.end_row();
                });
            if let Some(err) = &e.error {
                ui.label(RichText::new(err).color(t.error));
            }
            ui.horizontal(|ui| {
                save = ui.button("Save").clicked();
                cancel = ui.button("Cancel").clicked();
            });
            ui.label(
                RichText::new("Positions are 1-based and inclusive; an end before its start goes across the origin.")
                    .small()
                    .color(t.text_muted),
            );
        });
    (save, cancel)
}

/// What a CDS codes for: length, mass, pI, extinction coefficient and the sequence.
fn protein_card(ui: &mut Ui, t: &Tokens, session: &mut Session, id: FeatureId) {
    let Ok(p) = session.execute("protein_properties", json!({ "feature": id })) else {
        ui.label(RichText::new("The protein couldn't be worked out").color(t.error));
        return;
    };
    let num = |k: &str| p.get(k).and_then(Value::as_f64).unwrap_or(0.0);
    let seq = p.get("protein").and_then(Value::as_str).unwrap_or("");
    let shown: String = seq.chars().take(80).collect();
    let more = if seq.chars().count() > 80 {
        "\u{2026}"
    } else {
        ""
    };
    ui.label(
        RichText::new(format!(
            "{} aa \u{b7} {:.1} kDa \u{b7} pI {:.2} \u{b7} \u{3b5}280 {}",
            p.get("length").and_then(Value::as_u64).unwrap_or(0),
            num("mass_da") / 1000.0,
            num("isoelectric_point"),
            p.get("extinction_280")
                .and_then(|e| e.get("cystines"))
                .and_then(Value::as_u64)
                .unwrap_or(0),
        ))
        .strong(),
    );
    let stops = p.get("internal_stops").and_then(Value::as_u64).unwrap_or(0);
    if stops > 0 {
        ui.label(RichText::new(format!("! {stops} internal stop codon(s)")).color(t.error));
    }
    ui.label(RichText::new(format!("{shown}{more}")).monospace());
}

/// "Detect features": the library learned from the user's files, and what it finds here.
fn detection(
    ui: &mut Ui,
    t: &Tokens,
    session: &mut Session,
    state: &mut FeaturesState,
) -> Option<FeaturesAction> {
    let mut action = None;
    let mut clear = false;
    egui::CollapsingHeader::new(RichText::new("Detect known features").strong())
        .id_salt("detect-features")
        .default_open(false)
        .show(ui, |ui| {
            ui.label(
                RichText::new("Finds copies (either strand, up to 2% of bases differing) of known elements in this sequence: the built-in features (elements from public GenBank/RefSeq records, and fluorescent proteins from FPbase, CC BY-SA 4.0, found by translation), and a library SnappieSnap learns from a folder of your own .dna and GenBank files. Your library stays on this computer.")
                    .color(t.text_muted),
            );
            if state.library.is_none() {
                state.library = Some(
                    session
                        .execute("feature_library", json!({}))
                        .map_err(|e| e.message),
                );
            }
            ui.horizontal(|ui| {
                match &state.library {
                    Some(Ok(v)) => {
                        let n = v.get("count").and_then(Value::as_u64).unwrap_or(0);
                        let files = v.get("files").and_then(Value::as_u64).unwrap_or(0);
                        let folder = v.get("folder").and_then(Value::as_str).unwrap_or("");
                        let builtin = v.get("builtin").and_then(Value::as_u64).unwrap_or(0);
                        let mine = v.get("mine").and_then(Value::as_u64).unwrap_or(0);
                        let downloaded = v.get("downloaded").and_then(Value::as_u64).unwrap_or(0);
                        let proteins = v.get("proteins").and_then(Value::as_u64).unwrap_or(0);
                        let mut parts = vec![format!("{builtin} built-in elements"), format!("{proteins} FPbase proteins")];
                        if n > 0 {
                            parts.push(format!("{n} from {files} of your files"));
                        }
                        if mine > 0 {
                            parts.push(format!("{mine} added by you"));
                        }
                        if downloaded > 0 {
                            parts.push(format!("{downloaded} downloaded"));
                        }
                        let data = v.get("data_folder").and_then(Value::as_str).unwrap_or("");
                        let mut where_ = format!("Kept in {data}");
                        if !folder.is_empty() {
                            where_.push_str(&format!("\nLearned from {folder}"));
                        }
                        ui.label(parts.join(" + ")).on_hover_text(where_);
                    }
                    _ => {
                        ui.label(RichText::new("Built-in features only (no library of yours yet)").color(t.text_muted));
                    }
                }
                if ui
                    .button("Build from a folder\u{2026}")
                    .on_hover_text("Read every annotated sequence file in a folder and its subfolders")
                    .clicked()
                {
                    action = Some(FeaturesAction::PickLibraryFolder);
                }
                let downloaded = matches!(&state.library, Some(Ok(v)) if v.get("downloaded").and_then(Value::as_u64).unwrap_or(0) > 0);
                if downloaded {
                    if ui
                        .button("Delete downloaded library")
                        .on_hover_text("Remove the vector library downloaded from GenBank from this computer")
                        .clicked()
                    {
                        action = Some(FeaturesAction::Run("feature_library_delete_download", json!({})));
                    }
                } else if ui
                    .button("Download vector library\u{2026}")
                    .on_hover_text("Build a library from common vector records in GenBank (asks you to accept its terms first)")
                    .clicked()
                {
                    // The command states its own terms when asked without consent.
                    state.download_terms = Some(
                        session
                            .execute("feature_library_download", json!({}))
                            .err()
                            .map(|e| e.message)
                            .unwrap_or_default(),
                    );
                }
                if ui.button("Detect").clicked() {
                    state.detected = Some(
                        session
                            .execute("detect_features", json!({}))
                            .map(|v| {
                                v.get("features")
                                    .and_then(Value::as_array)
                                    .cloned()
                                    .unwrap_or_default()
                                    .into_iter()
                                    .map(|f| (f, true))
                                    .collect()
                            })
                            .map_err(|e| e.message),
                    );
                }
            });
            // The user's own entries, each removable.
            if let Some(Ok(v)) = &state.library {
                let mine: Vec<(String, String, u64, String)> = v
                    .get("features")
                    .and_then(Value::as_array)
                    .map(|list| {
                        list.iter()
                            .filter(|f| f.get("mine").and_then(Value::as_bool) == Some(true))
                            .map(|f| {
                                let text = |k: &str| f.get(k).and_then(Value::as_str).unwrap_or("").to_string();
                                (text("name"), text("kind"), f.get("length").and_then(Value::as_u64).unwrap_or(0), text("source"))
                            })
                            .collect()
                    })
                    .unwrap_or_default();
                if !mine.is_empty() {
                    egui::CollapsingHeader::new(format!("My features ({})", mine.len()))
                        .id_salt("my-features")
                        .show(ui, |ui| {
                            egui::Grid::new("my-features-grid")
                                .num_columns(5)
                                .spacing(vec2(t.spacing * 2.0, 2.0))
                                .show(ui, |ui| {
                                    for (name, kind, length, source) in &mine {
                                        ui.label(name);
                                        ui.label(kind);
                                        ui.label(format!("{length} bp"));
                                        ui.label(RichText::new(source).color(t.text_muted));
                                        if ui.small_button("Remove").clicked() {
                                            action = Some(FeaturesAction::Run(
                                                "feature_library_remove",
                                                json!({ "name": name }),
                                            ));
                                        }
                                        ui.end_row();
                                    }
                                });
                        });
                }
            }
            egui::CollapsingHeader::new("Design a feature")
                .id_salt("design-feature")
                .show(ui, |ui| {
                    let (name, kind, residues) = &mut state.design;
                    if kind.is_empty() {
                        *kind = "misc_feature".into();
                    }
                    ui.horizontal(|ui| {
                        ui.add(egui::TextEdit::singleline(name).hint_text("Name").desired_width(160.0));
                        egui::ComboBox::from_id_salt("design-kind")
                            .selected_text(kind.as_str())
                            .show_ui(ui, |ui| {
                                for k in KINDS {
                                    ui.selectable_value(kind, (*k).to_string(), *k);
                                }
                            });
                    });
                    ui.add(
                        egui::TextEdit::multiline(residues)
                            .hint_text("Sequence, 5' to 3' (A, C, G, T and IUPAC codes such as N)")
                            .desired_rows(2)
                            .desired_width(f32::INFINITY),
                    );
                    let ready = !name.trim().is_empty() && !residues.trim().is_empty();
                    if ui
                        .add_enabled(ready, egui::Button::new("Add to my library"))
                        .on_hover_text("It will be found in other sequences by Detect")
                        .clicked()
                    {
                        action = Some(FeaturesAction::Run(
                            "feature_library_add",
                            json!({ "name": name.trim(), "kind": kind, "residues": residues }),
                        ));
                        name.clear();
                        residues.clear();
                    }
                });
            match &mut state.detected {
                None => {}
                Some(Err(e)) => {
                    ui.label(RichText::new(e.as_str()).color(t.error));
                }
                Some(Ok(found)) if found.is_empty() => {
                    ui.label(
                        RichText::new("Nothing new found: every known feature here is already annotated (or none occurs).")
                            .color(t.text_muted),
                    );
                }
                Some(Ok(found)) => {
                    egui::ScrollArea::vertical()
                        .id_salt("detected")
                        .max_height(200.0)
                        .show(ui, |ui| {
                            egui::Grid::new("detected-grid")
                                .num_columns(6)
                                .spacing(vec2(t.spacing * 2.0, 2.0))
                                .show(ui, |ui| {
                                    for (f, on) in found.iter_mut() {
                                        let text = |k: &str| {
                                            f.get(k).and_then(Value::as_str).unwrap_or("").to_string()
                                        };
                                        ui.checkbox(on, text("name"));
                                        let range = f
                                            .get("range")
                                            .and_then(|r| serde_json::from_value::<Range>(r.clone()).ok());
                                        let len = session.document().map_or(0, |d| d.molecule.len());
                                        ui.label(range.map_or(String::new(), |r| {
                                            let (a, b) = shown(r, len);
                                            format!("{a}..{b}")
                                        }));
                                        ui.label(text("strand"));
                                        ui.label(text("kind"));
                                        let differ = f.get("mismatches").and_then(Value::as_u64).unwrap_or(0);
                                        if differ == 0 {
                                            ui.label("exact");
                                        } else {
                                            let identity = f.get("identity").and_then(Value::as_f64).unwrap_or(0.0);
                                            let n = |k: &str| f.get(k).and_then(Value::as_u64).unwrap_or(0);
                                            let mut kinds = Vec::new();
                                            for (k, label) in [("substitutions", "substituted"), ("insertions", "inserted"), ("deletions", "missing")] {
                                                if n(k) > 0 {
                                                    kinds.push(format!("{} {label}", n(k)));
                                                }
                                            }
                                            let unit = if text("compared") == "protein" { "amino acid(s)" } else { "base(s)" };
                                            ui.label(RichText::new(format!("{identity}%")).color(t.accent))
                                                .on_hover_text(format!(
                                                    "{differ} {unit} differ from the known sequence: {}",
                                                    kinds.join(", ")
                                                ));
                                        }
                                        ui.label(RichText::new(format!("from {}", text("source"))).color(t.text_muted));
                                        ui.end_row();
                                    }
                                });
                        });
                    ui.horizontal(|ui| {
                        let names: Vec<String> = found
                            .iter()
                            .filter(|(_, on)| *on)
                            .filter_map(|(f, _)| f.get("name").and_then(Value::as_str).map(str::to_string))
                            .collect();
                        if ui
                            .add_enabled(!names.is_empty(), egui::Button::new(format!("Add {} feature(s)", names.len())))
                            .clicked()
                        {
                            action = Some(FeaturesAction::Run(
                                "add_detected_features",
                                json!({ "names": names }),
                            ));
                            clear = true;
                        }
                    });
                }
            }
        });
    if clear {
        state.detected = None;
    }
    // The download's terms, to accept or decline.
    let mut close = false;
    if let Some(terms) = &state.download_terms {
        egui::Window::new("Download vector library")
            .collapsible(false)
            .resizable(false)
            .anchor(egui::Align2::CENTER_CENTER, vec2(0.0, 0.0))
            .show(ui.ctx(), |ui| {
                ui.set_max_width(460.0);
                ui.label(terms.as_str());
                ui.add_space(t.spacing);
                ui.horizontal(|ui| {
                    if ui.button("Accept and download").clicked() {
                        action = Some(FeaturesAction::Run(
                            "feature_library_download",
                            json!({ "accept": true }),
                        ));
                        close = true;
                    }
                    if ui.button("Cancel").clicked() {
                        close = true;
                    }
                });
            });
    }
    if close {
        state.download_terms = None;
    }
    action
}

/// Draw the tab; returns a command to run.
pub fn show(
    ui: &mut Ui,
    t: &Tokens,
    session: &mut Session,
    state: &mut FeaturesState,
) -> Option<FeaturesAction> {
    let doc = session.document()?;
    let features = doc.molecule.features().to_vec();
    let len = doc.molecule.len();
    let selection = doc.selection.filter(|s| !s.is_empty());
    let mut action = None;

    // The editor, when open.
    if let Some(e) = &mut state.editor {
        let (save, cancel) = editor_ui(ui, t, e);
        if cancel {
            state.editor = None;
        } else if save {
            let existing = e.id.and_then(|id| features.iter().find(|f| f.id == id));
            match e.params(existing, len) {
                Ok((id, params)) => match session.execute(id, params) {
                    Ok(_) => state.editor = None,
                    Err(err) => e.error = Some(err.message),
                },
                Err(err) => e.error = Some(err),
            }
        }
        ui.add_space(t.spacing);
    }

    ui.horizontal(|ui| {
        if ui
            .add_enabled(
                selection.is_some(),
                egui::Button::new("New feature from selection"),
            )
            .on_disabled_hover_text("Select a stretch of the sequence first")
            .clicked()
            && let Some(sel) = selection
        {
            state.editor = Some(Editor::new_over(sel, len));
        }
        if ui
            .button("Import annotations\u{2026}")
            .on_hover_text("Add the features of a GFF3 file (for example from NCBI or a genome browser) to this sequence")
            .clicked()
        {
            action = Some(FeaturesAction::PickAnnotations);
        }
        ui.add(
            egui::TextEdit::singleline(&mut state.search)
                .hint_text("Find features (name or type)")
                .desired_width(220.0),
        );
    });
    ui.add_space(t.spacing);
    if let Some(a) = detection(ui, t, session, state) {
        action = Some(a);
    }
    ui.add_space(t.spacing);

    let needle = state.search.trim().to_lowercase();
    let mut rows: Vec<&Feature> = features
        .iter()
        .filter(|f| {
            needle.is_empty()
                || f.name.to_lowercase().contains(&needle)
                || f.kind.to_lowercase().contains(&needle)
        })
        .collect();
    let size = |f: &Feature| -> usize { f.segments.iter().map(|s| s.range.len_on(len)).sum() };
    rows.sort_by(|a, b| {
        let o = match state.sort {
            Column::Location => span(a).map(|r| r.start).cmp(&span(b).map(|r| r.start)),
            Column::Name => a.name.to_lowercase().cmp(&b.name.to_lowercase()),
            Column::Size => size(a).cmp(&size(b)),
            Column::Direction => direction(a.strand).cmp(direction(b.strand)),
            Column::Kind => a.kind.cmp(&b.kind),
        }
        .then_with(|| a.id.0.cmp(&b.id.0));
        if state.descending { o.reverse() } else { o }
    });
    ui.label(
        RichText::new(if rows.len() == features.len() {
            format!("{} features", features.len())
        } else {
            format!("{} of {} features", rows.len(), features.len())
        })
        .color(t.text_muted),
    );

    egui::ScrollArea::vertical()
        .id_salt("features-table")
        .show(ui, |ui| {
            egui::Grid::new("features")
                .striped(true)
                .num_columns(6)
                .spacing(vec2(t.spacing * 3.0, t.spacing))
                .show(ui, |ui| {
                    for (column, label) in [
                        (Column::Name, "Feature"),
                        (Column::Location, "Location"),
                        (Column::Size, "Size"),
                        (Column::Direction, "Direction"),
                        (Column::Kind, "Type"),
                    ] {
                        let on = state.sort == column;
                        let clicked = ui
                            .horizontal(|ui| {
                                let r = ui
                                    .add(
                                        egui::Label::new(RichText::new(label).strong())
                                            .sense(Sense::click()),
                                    )
                                    .on_hover_text("Sort by this column (again: the other way)");
                                if on {
                                    // A small triangle: up when ascending, down when descending.
                                    let (rect, _) =
                                        ui.allocate_exact_size(vec2(8.0, 8.0), Sense::hover());
                                    let c = rect.center();
                                    let d = if state.descending { 1.0 } else { -1.0 };
                                    ui.painter().add(egui::Shape::convex_polygon(
                                        vec![
                                            egui::pos2(c.x - 3.5, c.y - 2.0 * d),
                                            egui::pos2(c.x + 3.5, c.y - 2.0 * d),
                                            egui::pos2(c.x, c.y + 2.5 * d),
                                        ],
                                        t.text_muted,
                                        egui::Stroke::NONE,
                                    ));
                                }
                                r.clicked()
                            })
                            .inner;
                        if clicked {
                            if on {
                                state.descending = !state.descending;
                            } else {
                                state.sort = column;
                                state.descending = false;
                            }
                        }
                    }
                    ui.label("");
                    ui.end_row();
                    for f in &rows {
                        ui.horizontal(|ui| {
                            let (r, _) = ui.allocate_exact_size(vec2(10.0, 10.0), Sense::hover());
                            ui.painter()
                                .rect_filled(r, CornerRadius::same(2), t.color_of(f));
                            let name = ui
                                .add(egui::Label::new(&f.name).sense(Sense::click()))
                                .on_hover_text("Click to select it; double-click to edit");
                            if name.double_clicked() {
                                state.editor = Some(Editor::of(f, len));
                            } else if name.clicked()
                                && let Some(s) = span(f)
                            {
                                action = Some(FeaturesAction::Run("select", json!({ "range": s })));
                            }
                        });
                        let location = span(f).map_or(String::new(), |s| {
                            let (a, b) = shown(s, len);
                            if f.segments.len() > 1 {
                                format!("{a}..{b} ({} parts)", f.segments.len())
                            } else {
                                format!("{a}..{b}")
                            }
                        });
                        let parts: Vec<String> = f
                            .segments
                            .iter()
                            .map(|s| {
                                let (a, b) = shown(s.range, len);
                                format!("{a}..{b}")
                            })
                            .collect();
                        ui.label(location).on_hover_text(parts.join(", "));
                        ui.label(format!("{} bp", size(f)));
                        ui.label(direction(f.strand));
                        let kind = ui.label(&f.kind);
                        if f.kind == "CDS" {
                            // The protein it codes for, worked out only while hovered.
                            kind.on_hover_ui(|ui| protein_card(ui, t, session, f.id));
                        }
                        ui.horizontal(|ui| {
                            if ui.small_button("Edit").clicked() {
                                state.editor = Some(Editor::of(f, len));
                            }
                            if ui
                                .add_enabled(f.segments.len() == 1, egui::Button::new("Keep").small())
                                .on_hover_text("Add to your feature library, so Detect finds it in other sequences")
                                .clicked()
                            {
                                action = Some(FeaturesAction::Run(
                                    "feature_library_add",
                                    json!({ "feature": f.id }),
                                ));
                            }
                            if ui.small_button("Delete").clicked() {
                                action = Some(FeaturesAction::Run(
                                    "remove_feature",
                                    json!({ "id": f.id }),
                                ));
                            }
                        });
                        ui.end_row();
                    }
                });
            if features.is_empty() {
                ui.label(RichText::new("No features").color(t.text_muted));
            }
        });
    action
}

#[cfg(test)]
mod tests {
    use snappiesnap_engine::Document;
    use snappiesnap_engine::doc::{
        Alphabet, Molecule, NewFeature, Qualifier, Segment, Seq, Topology,
    };

    use super::*;

    fn session() -> Session {
        let mut m = Molecule::new(
            "p",
            Seq::new(Alphabet::Dna, "ACGT".repeat(25)).unwrap(),
            Topology::Circular,
        )
        .unwrap();
        m.add_feature(NewFeature {
            color: Some("#ff0000".into()),
            name: "gene".into(),
            kind: "CDS".into(),
            strand: Strand::Reverse,
            segments: vec![
                Segment {
                    range: Range::new(60, 70),
                },
                Segment {
                    range: Range::new(10, 20),
                },
            ],
            qualifiers: vec![Qualifier {
                name: "note".into(),
                values: vec!["a".into(), "b".into()],
            }],
        })
        .unwrap();
        let mut s = Session::new();
        s.open_document(Document::new(m));
        s
    }

    #[test]
    fn the_editor_round_trips_a_feature() {
        let mut s = session();
        let f = s.document().unwrap().molecule.features()[0].clone();
        let mut e = Editor::of(&f, 100);
        // Parts 1-based, in the feature's order; both note values.
        assert_eq!(
            e.parts,
            vec![("61".into(), "70".into()), ("11".into(), "20".into())]
        );
        assert_eq!(e.qualifiers.len(), 2);
        // Unchanged: saving changes nothing (the colour is kept).
        let (id, p) = e.params(Some(&f), 100).unwrap();
        assert_eq!(id, "update_feature");
        s.execute(id, p).unwrap();
        assert_eq!(s.document().unwrap().molecule.features()[0], f);
        // Edited: name, type, direction, one part across the origin, a new qualifier.
        e.name = "renamed".into();
        e.kind = "gene".into();
        e.strand = Strand::Forward;
        e.parts = vec![("95".into(), "5".into())];
        e.qualifiers.push(("gene".into(), "abc".into()));
        let (id, p) = e.params(Some(&f), 100).unwrap();
        s.execute(id, p).unwrap();
        let g = s.document().unwrap().molecule.features()[0].clone();
        assert_eq!(
            (g.name.as_str(), g.kind.as_str(), g.strand),
            ("renamed", "gene", Strand::Forward)
        );
        assert_eq!(g.segments[0].range, Range::new(94, 5));
        assert_eq!(g.qualifiers.len(), 2);
        assert_eq!(g.color.as_deref(), Some("#ff0000"));
        // A colour picked in the editor is saved; "Default" removes it (undoably).
        e.parts = vec![("95".into(), "5".into())];
        e.color = Some("#00ff00".into());
        let (id, p) = e.params(Some(&g), 100).unwrap();
        s.execute(id, p).unwrap();
        let green = s.document().unwrap().molecule.features()[0].clone();
        assert_eq!(green.color.as_deref(), Some("#00ff00"));
        e.color = None;
        let (id, p) = e.params(Some(&green), 100).unwrap();
        s.execute(id, p).unwrap();
        assert_eq!(s.document().unwrap().molecule.features()[0].color, None);
        s.execute("undo", serde_json::json!({})).unwrap();
        assert_eq!(
            s.document().unwrap().molecule.features()[0]
                .color
                .as_deref(),
            Some("#00ff00")
        );
        // Bad positions are explained, not sent.
        e.parts = vec![("0".into(), "5".into())];
        assert!(e.params(Some(&f), 100).unwrap_err().contains("1 to 100"));
        e.parts.clear();
        assert!(e.params(Some(&f), 100).is_err());
        // A new feature over a selection, including the whole molecule.
        let n = Editor::new_over(Range::new(0, 0), 100);
        assert_eq!(n.parts, vec![("1".into(), "100".into())]);
        let (id, p) = n.params(None, 100).unwrap();
        assert_eq!(id, "add_feature");
        s.execute(id, p).unwrap();
        assert_eq!(s.document().unwrap().molecule.features().len(), 2);
    }

    #[test]
    fn the_tab_draws() {
        let mut s = session();
        let mut state = FeaturesState {
            editor: Some(Editor::new_over(Range::new(3, 9), 100)),
            ..Default::default()
        };
        let t = Tokens::light();
        let image = crate::offscreen::render(
            900,
            600,
            1.0,
            |ctx| t.apply(ctx),
            |ui| {
                let _ = show(ui, &t, &mut s, &mut state);
            },
        )
        .unwrap();
        if let Ok(path) = std::env::var("SNAPPIESNAP_FEATURES_PNG") {
            crate::offscreen::write_png(&path, &image).unwrap();
        }
    }
}
