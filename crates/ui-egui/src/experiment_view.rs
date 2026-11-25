//! Experiments, opened from the toolbar's Experiments menu in tabs of their own (beside the
//! document tabs), so they are bound to no construct.
//! For now one kind, Simulate Gel: one gel across constructs. Step 1 chooses and orders the
//! constructs (the open documents to start with; one lane each, lane 1 being the ladder).
//! Step 2 sets each lane's method (uncut, digest or PCR), lane by lane or for all lanes at
//! once. Lanes run with `snappiesnap_engine::gel::experiment`, the same code as the
//! `gel_experiment` command and as each construct's own gel. The experiment is UI state
//! ([`Experiment`]), so agents can read and set it with `ui.state` / `ui.set`.

use eframe::egui::{self, RichText, Ui};
use snappiesnap_engine::Session;
use snappiesnap_engine::bio::gel::LADDERS;
use snappiesnap_engine::gel::{LaneSpec, RunLane, experiment};

use crate::gel_view::draw_gel;
use crate::state::{Experiment, ExperimentKind, ExperimentStep};
use crate::theme::Tokens;

/// The method being set up for the selected lane, before it is applied.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
enum Kind {
    #[default]
    Uncut,
    Digest,
    Pcr,
}

/// Inputs of the method editor, and the gel's cache.
#[derive(Default)]
pub struct ExperimentInputs {
    kind: Kind,
    enzymes: String,
    forward: String,
    reverse: String,
    /// The lane the editor was last filled from.
    loaded: Option<usize>,
    cache: Option<(String, Vec<RunLane>)>,
    /// Diagnostic Digest: single enzymes only.
    singles_only: bool,
    /// The enzyme set the app is showing (what the lab holds), set each frame by the app.
    /// When there is one, digests are suggested from it, so every suggestion can be run.
    pub enzyme_pool: Option<String>,
    /// Diagnostic Digest: also search enzymes the set hasn't got (marked, and listed after
    /// the ones it has).
    include_out_of_stock: bool,
    /// What the shown suggestions were asked for, and the `diagnostic_digest` answer (None
    /// while it is being computed).
    suggest_key: Option<String>,
    suggestions: Option<Result<serde_json::Value, String>>,
    /// Check Reads.
    pub reads: crate::reads_view::ReadsInputs,
    /// Cloning.
    pub cloning: crate::cloning_experiment::CloningInputs,
}

impl ExperimentInputs {
    /// The answer of `diagnostic_digest` to the last [`ExperimentAction::Suggest`].
    pub fn set_suggestions(&mut self, answer: Result<serde_json::Value, String>) {
        self.suggestions = Some(answer);
    }

    /// Fill the editor from a lane's method (a lane without one keeps the last method, so it
    /// can be applied lane after lane).
    fn load(&mut self, method: Option<&LaneSpec>) {
        match method {
            Some(LaneSpec::Uncut) => self.kind = Kind::Uncut,
            Some(LaneSpec::Digest { enzymes }) => {
                self.kind = Kind::Digest;
                self.enzymes = enzymes.join(", ");
            }
            Some(LaneSpec::Pcr { forward, reverse }) => {
                self.kind = Kind::Pcr;
                self.forward.clone_from(forward);
                self.reverse.clone_from(reverse);
            }
            _ => {}
        }
    }

    /// The method in the editor, if it is complete.
    fn method(&self) -> Option<LaneSpec> {
        match self.kind {
            Kind::Uncut => Some(LaneSpec::Uncut),
            Kind::Digest => {
                let enzymes: Vec<String> = self
                    .enzymes
                    .split([',', ' ', '+'])
                    .map(str::trim)
                    .filter(|s| !s.is_empty())
                    .map(str::to_string)
                    .collect();
                (!enzymes.is_empty()).then_some(LaneSpec::Digest { enzymes })
            }
            Kind::Pcr => {
                (!self.forward.trim().is_empty() && !self.reverse.trim().is_empty()).then(|| {
                    LaneSpec::Pcr {
                        forward: self.forward.trim().to_string(),
                        reverse: self.reverse.trim().to_string(),
                    }
                })
            }
        }
    }

    /// The gel's lanes, recomputed only when a document, the enzymes or the setup changed.
    fn lanes(&mut self, session: &Session, exp: &Experiment) -> &[RunLane] {
        let open: Vec<u64> = session.documents().iter().map(|d| d.id).collect();
        let key = format!(
            "{}|{open:?}|{}",
            session.revision(),
            serde_json::to_string(&(&exp.ladder, exp.agarose, &exp.lanes)).unwrap_or_default()
        );
        if self.cache.as_ref().is_none_or(|(k, _)| *k != key) {
            let lanes = experiment(session, &exp.ladder, exp.agarose, &exp.lanes);
            self.cache = Some((key, lanes));
        }
        self.cache
            .as_ref()
            .map_or(&[], |(_, lanes)| lanes.as_slice())
    }
}

/// What the user asked for in an experiment tab.
#[derive(Debug, Clone, PartialEq)]
pub enum ExperimentAction {
    /// Save the gel as an image.
    SaveImage,
    /// Save the experiment to its file (Save As when it has none).
    Save,
    /// Save the experiment to a new file.
    SaveAs,
    /// Run `diagnostic_digest` with these params; give the answer to
    /// [`ExperimentInputs::set_suggestions`].
    Suggest(serde_json::Value),
    /// Check Reads: pick read files to add.
    AddReads,
    /// Check Reads: run `check_reads` with these params; give the answer to
    /// [`crate::reads_view::ReadsInputs::set_result`].
    CheckReads(serde_json::Value),
    /// Run a command (opening a read, or opening a cloning product).
    Run(&'static str, serde_json::Value),
    /// Cloning: read a file with `file_info`; the answer goes to `CloningInputs::set_info`.
    CloningInfo(String),
    /// Cloning: run the reaction as a preview; the answer goes to `CloningInputs::set_preview`.
    CloningPreview(&'static str, serde_json::Value),
    /// Cloning: open the sample constructs and set the experiment up to clone them.
    CloningSamples,
}

/// An experiment tab: its name (editable) with Save buttons, then the current step.
pub fn show(
    ui: &mut Ui,
    t: &Tokens,
    session: &Session,
    exp: &mut Experiment,
    inputs: &mut ExperimentInputs,
    folder: Option<&str>,
    pane: &mut crate::cloning_experiment::Pane,
) -> Option<ExperimentAction> {
    let mut action = None;
    ui.horizontal(|ui| {
        ui.label(RichText::new(exp.kind.label()).color(t.text_muted));
        let r = ui.add(
            egui::TextEdit::singleline(&mut exp.name)
                .font(egui::TextStyle::Heading)
                .char_limit(Experiment::MAX_NAME)
                .desired_width(320.0),
        );
        if r.lost_focus() && exp.name.trim().is_empty() {
            exp.name = exp.kind.label().to_string();
        }
        let has_content = !exp.lanes.is_empty() || exp.reference.is_some() || !exp.reads.is_empty();
        if !matches!(
            exp.step,
            ExperimentStep::NotStarted | ExperimentStep::Simulation
        ) && has_content
        {
            ui.separator();
            let tip = match &exp.path {
                Some(p) => format!("Save to {p}"),
                None => "Save the experiment and its constructs in a .sslab file".to_string(),
            };
            if ui
                .add_enabled(
                    exp.unsaved() || exp.path.is_none(),
                    egui::Button::new("Save"),
                )
                .on_hover_text(tip)
                .clicked()
            {
                action = Some(ExperimentAction::Save);
            }
            if ui
                .button("Save As\u{2026}")
                .on_hover_text("Save the experiment in a new .sslab file")
                .clicked()
            {
                action = Some(ExperimentAction::SaveAs);
            }
        }
    });
    ui.add_space(t.spacing);
    if exp.kind == ExperimentKind::Cloning {
        use crate::cloning_experiment::CloningAction;
        // The enzymes the lab holds reach the cloning fields too, so a reaction designed here
        // is one that can be run.
        inputs.cloning.pool = inputs.enzyme_pool.clone();
        return match crate::cloning_experiment::show(
            ui,
            t,
            session,
            exp,
            &mut inputs.cloning,
            folder,
            pane,
        ) {
            Some(CloningAction::Info(path)) => Some(ExperimentAction::CloningInfo(path)),
            Some(CloningAction::Preview(c, p)) => Some(ExperimentAction::CloningPreview(c, p)),
            Some(CloningAction::Open(c, p)) => Some(ExperimentAction::Run(c, p)),
            Some(CloningAction::Samples) => Some(ExperimentAction::CloningSamples),
            None => action,
        };
    }
    if session.documents().is_empty() {
        ui.label(
            RichText::new("Open the constructs for this experiment first (Open in the toolbar).")
                .color(t.text_muted),
        );
        return action;
    }
    if exp.kind == ExperimentKind::Reads && exp.step != ExperimentStep::NotStarted {
        use crate::reads_view::{ReadsAction, results, setup};
        let a = if exp.step == ExperimentStep::Results {
            results(ui, t, session, exp, &mut inputs.reads)
        } else {
            setup(ui, t, session, exp, &mut inputs.reads, folder)
        };
        return match a {
            Some(ReadsAction::AddFiles) => Some(ExperimentAction::AddReads),
            Some(ReadsAction::Check(p)) => Some(ExperimentAction::CheckReads(p)),
            Some(ReadsAction::Run(command, p)) => Some(ExperimentAction::Run(command, p)),
            None => action,
        };
    }
    match exp.step {
        ExperimentStep::NotStarted => not_started(ui, t, session, exp),
        ExperimentStep::Simulation => simulation(ui, t, session, exp),
        ExperimentStep::Constructs => constructs(ui, t, session, exp),
        ExperimentStep::Results => {}
        ExperimentStep::Enzymes => {
            if let Some(params) = enzymes_step(ui, t, session, exp, inputs) {
                action = Some(ExperimentAction::Suggest(params));
            }
        }
        ExperimentStep::Lanes => {
            if lanes(ui, t, session, exp, inputs) {
                action = Some(ExperimentAction::SaveImage);
            }
        }
        // The cloning steps are drawn above, before this match.
        ExperimentStep::Technique
        | ExperimentStep::Vector
        | ExperimentStep::Insert
        | ExperimentStep::Product => {}
    }
    action
}

/// The open documents, as (id, name).
fn open_documents(session: &Session) -> Vec<(u64, String)> {
    session
        .documents()
        .into_iter()
        .map(|d| (d.id, d.name))
        .collect()
}

fn name_of(session: &Session, id: u64) -> Option<String> {
    session
        .document_by_id(id)
        .map(|d| d.molecule.name().to_string())
}

/// Step 1 of Simulate Gel: which gel to simulate. Every later step follows the choice, so the
/// kinds that end in a gel share one menu item.
fn simulation(ui: &mut Ui, t: &Tokens, session: &Session, exp: &mut Experiment) {
    steps(ui, t, exp.kind, 1);
    ui.label(RichText::new("What would you like to simulate?").strong());
    ui.add_space(t.spacing);
    for kind in ExperimentKind::SIMULATIONS {
        let chosen = exp.kind == kind;
        let r = ui.selectable_label(chosen, RichText::new(kind.simulation_label()).strong());
        if r.clicked() {
            exp.set_simulation(kind);
        }
        ui.label(RichText::new(description_of(kind)).color(t.text_muted));
        ui.add_space(t.spacing);
    }
    // Coming back to this step keeps the lanes already chosen; starting fresh takes the
    // documents open now.
    if exp.lanes.is_empty() {
        if ui.button("Choose constructs").clicked() {
            exp.start(open_documents(session).into_iter().map(|(id, _)| id));
        }
    } else if ui.button("Continue: choose constructs").clicked() {
        exp.step = crate::state::ExperimentStep::Constructs;
    }
}

/// What each simulation does, in the first step and on the intro screen.
fn description_of(kind: ExperimentKind) -> &'static str {
    match kind {
        ExperimentKind::Gel => {
            "Run several constructs side by side on one gel, each lane with its own method: \
             uncut, a restriction digest or a PCR."
        }
        ExperimentKind::Diagnostic => {
            "Find restriction enzymes, alone or in pairs, whose digests tell your constructs \
             apart on a gel: a clone from its empty vector, one insert orientation from the \
             other. With one construct, find a digest whose pattern confirms it."
        }
        ExperimentKind::Reads => {
            "Check a construct with sequencing reads: which stretches the reads cover, how \
             each read lies, and every difference from the construct."
        }
        ExperimentKind::Cloning => {
            "Clone two constructs together, step by step, and open the product as a construct \
             of its own."
        }
    }
}

fn not_started(ui: &mut Ui, t: &Tokens, session: &Session, exp: &mut Experiment) {
    let (title, text) = match exp.kind {
        ExperimentKind::Gel => (
            "Gel experiment",
            "Run several constructs side by side on one gel, each lane with its own method: \
             uncut, a restriction digest or a PCR. Start with the constructs you have open; \
             you can add, remove and reorder them next.",
        ),
        ExperimentKind::Reads => (
            "Check reads",
            "Check a construct with sequencing reads: which stretches the reads cover, how \
             each read lies, and every difference from the construct, pooled over the reads \
             (confirmed by several, disputed, or seen in one), with the traces to judge it. \
             Start with the construct you have open; choose the read files next.",
        ),
        ExperimentKind::Cloning => ("Cloning", description_of(ExperimentKind::Cloning)),
        ExperimentKind::Diagnostic => (
            "Diagnostic digest",
            "Find restriction enzymes, alone or in pairs, whose digests tell your constructs \
             apart on a gel: a clone from its empty vector, one insert orientation from the \
             other. With one construct, find a digest whose pattern confirms it. Start with \
             the constructs you have open; you can add, remove and reorder them next.",
        ),
    };
    ui.label(RichText::new(title).strong());
    ui.label(RichText::new(text).color(t.text_muted));
    ui.add_space(t.spacing);
    if ui.button("Choose constructs").clicked() {
        exp.start(open_documents(session).into_iter().map(|(id, _)| id));
    }
}

/// Step 2: the constructs, one lane each, in lane order.
fn constructs(ui: &mut Ui, t: &Tokens, session: &Session, exp: &mut Experiment) {
    steps(ui, t, exp.kind, 2);
    ui.label(
        RichText::new(format!(
            "Each construct runs in its own lane, top of the list first; lane 1 is the ladder. \
             {} lanes in all.",
            exp.lanes.len() + 1
        ))
        .color(t.text_muted),
    );
    ui.add_space(t.spacing);
    let open = open_documents(session);
    ui.vertical(|ui| {
        ui.set_max_width(520.0);
        ui.horizontal(|ui| {
            ui.label(RichText::new("Lane 1").color(t.text_muted));
            ui.label(format!("Ladder ({})", exp.ladder));
        });
        let mut action: Option<(usize, i8)> = None;
        let count = exp.lanes.len();
        for (i, lane) in exp.lanes.iter().enumerate() {
            ui.horizontal(|ui| {
                ui.label(RichText::new(format!("Lane {}", i + 2)).color(t.text_muted));
                match name_of(session, lane.document) {
                    Some(name) => ui.add(egui::Label::new(name).truncate()),
                    None => ui.label(RichText::new("closed construct").color(t.error)),
                };
                ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                    if ui.small_button("Remove").clicked() {
                        action = Some((i, 0));
                    }
                    if ui
                        .add_enabled(i + 1 < count, egui::Button::new("Down").small())
                        .on_hover_text("Move one lane later")
                        .clicked()
                    {
                        action = Some((i, 1));
                    }
                    if ui
                        .add_enabled(i > 0, egui::Button::new("Up").small())
                        .on_hover_text("Move one lane earlier")
                        .clicked()
                    {
                        action = Some((i, -1));
                    }
                });
            });
        }
        match action {
            Some((i, 0)) => exp.remove(i),
            Some((i, d)) => exp.move_lane(i, d < 0),
            None => {}
        }
        if exp.lanes.is_empty() {
            ui.label(RichText::new("No constructs yet: add one below.").color(t.text_muted));
        }
        ui.add_space(t.spacing);
        ui.horizontal(|ui| {
            let full = exp.lanes.len() >= Experiment::MAX_LANES;
            ui.add_enabled_ui(!full, |ui| {
                egui::ComboBox::from_id_salt("experiment-add")
                    .selected_text("Add a construct\u{2026}")
                    .width(200.0)
                    .show_ui(ui, |ui| {
                        for (id, name) in &open {
                            if ui.selectable_label(false, name).clicked() {
                                exp.add(*id);
                            }
                        }
                    });
            });
            if ui
                .button("Use the open documents")
                .on_hover_text("Start the list again from the documents open now")
                .clicked()
            {
                let keep = (exp.agarose, exp.ladder.clone());
                exp.start(open.iter().map(|(id, _)| *id));
                (exp.agarose, exp.ladder) = keep;
            }
        });
        ui.add_space(t.spacing * 2.0);
        let diagnostic = exp.kind == ExperimentKind::Diagnostic;
        let next = if diagnostic {
            "Continue: find enzymes"
        } else {
            "Continue: set up lanes"
        };
        ui.horizontal(|ui| {
            if exp.kind.chooses_simulation()
                && ui
                    .button("Back")
                    .on_hover_text("Back to the simulation; the constructs are kept")
                    .clicked()
            {
                exp.step = ExperimentStep::Simulation;
            }
            if ui
                .add_enabled(!exp.lanes.is_empty(), egui::Button::new(next))
                .clicked()
            {
                exp.step = if diagnostic {
                    ExperimentStep::Enzymes
                } else {
                    ExperimentStep::Lanes
                };
                exp.selected = exp.selected.or(Some(0));
            }
        });
    });
}

/// Step 3 (4 for a diagnostic digest): a method per lane, and the gel.
fn lanes(
    ui: &mut Ui,
    t: &Tokens,
    session: &Session,
    exp: &mut Experiment,
    inputs: &mut ExperimentInputs,
) -> bool {
    let mut save = false;
    let diagnostic = exp.kind == ExperimentKind::Diagnostic;
    steps(ui, t, exp.kind, if diagnostic { 4 } else { 3 });
    if exp.selected.is_some_and(|s| s >= exp.lanes.len()) {
        exp.selected = None;
    }
    if inputs.loaded != exp.selected {
        let method = exp
            .selected
            .and_then(|s| exp.lanes.get(s))
            .and_then(|l| l.method.as_ref());
        inputs.load(method);
        inputs.loaded = exp.selected;
    }
    let run = inputs.lanes(session, exp).to_vec();
    ui.horizontal_top(|ui| {
        ui.vertical(|ui| {
            ui.set_width(300.0);
            ui.horizontal(|ui| {
                if diagnostic {
                    if ui.button("Back to enzymes").clicked() {
                        exp.step = ExperimentStep::Enzymes;
                    }
                } else if ui.button("Back to constructs").clicked() {
                    exp.step = ExperimentStep::Constructs;
                }
                save = ui
                    .button("Save image\u{2026}")
                    .on_hover_text("Save this gel as a PNG or JPEG image, as shown or inverted")
                    .clicked();
                if ui
                    .button("Start over")
                    .on_hover_text("Start this experiment again from the documents open now (lane methods are cleared)")
                    .clicked()
                {
                    exp.start(open_documents(session).into_iter().map(|(id, _)| id));
                }
            });
            ui.add_space(t.spacing);
            ui.horizontal(|ui| {
                ui.label("Agarose");
                ui.add(
                    egui::DragValue::new(&mut exp.agarose)
                        .range(0.3..=2.0)
                        .speed(0.05)
                        .fixed_decimals(1)
                        .suffix(" %"),
                );
                ui.label("Ladder");
                egui::ComboBox::from_id_salt("experiment-ladder")
                    .selected_text(&exp.ladder)
                    .show_ui(ui, |ui| {
                        for (name, _) in LADDERS {
                            ui.selectable_value(&mut exp.ladder, name.to_string(), *name);
                        }
                    });
            });
            ui.add_space(t.spacing);
            ui.label(RichText::new("Lanes").strong());
            ui.label(RichText::new(format!("1  Ladder ({})", exp.ladder)).color(t.text_muted));
            for (i, lane) in exp.lanes.iter().enumerate() {
                let name = name_of(session, lane.document)
                    .unwrap_or_else(|| "closed construct".to_string());
                let method = lane
                    .method
                    .as_ref()
                    .map_or("no method yet".to_string(), LaneSpec::title);
                let error = run.get(i + 1).and_then(|(_, r)| r.as_ref().err());
                let mut text = RichText::new(format!("{}  {name}: {method}", i + 2));
                if error.is_some() {
                    text = text.color(t.error);
                } else if lane.method.is_none() {
                    text = text.color(t.text_muted);
                }
                let mut r = ui.selectable_label(exp.selected == Some(i), text);
                if let Some(e) = error {
                    r = r.on_hover_text(e.as_str());
                }
                if r.clicked() {
                    exp.selected = Some(i);
                }
            }
            ui.add_space(t.spacing);
            ui.separator();
            method_editor(ui, t, session, exp, inputs);
            ui.add_space(t.spacing * 2.0);
            ui.label(
                RichText::new(
                    "Click a lane in the list or on the gel to set it up. Each lane runs as on \
                     its construct's own gel: migration log-linear over the agarose's \
                     separation range (Sambrook & Russell 2001), brightness by DNA mass.",
                )
                .small()
                .color(t.text_muted),
            );
        });
        ui.separator();
        let selected = exp.selected.map_or(Some(usize::MAX), |s| Some(s + 1));
        if let Some(lane) = draw_gel(ui, t, &run, selected)
            && let Some(i) = lane.checked_sub(1)
        {
            exp.selected = Some(i);
        }
    });
    save
}

/// The method of the selected lane: uncut, a digest or a PCR, applied to that lane or to all.
fn method_editor(
    ui: &mut Ui,
    t: &Tokens,
    session: &Session,
    exp: &mut Experiment,
    inputs: &mut ExperimentInputs,
) {
    let Some(i) = exp.selected else {
        ui.label(RichText::new("Select a lane to set its method.").color(t.text_muted));
        return;
    };
    let Some(document) = exp.lanes.get(i).map(|l| l.document) else {
        return;
    };
    let Some(doc) = session.document_by_id(document) else {
        ui.label(
            RichText::new("This lane's construct is no longer open. Remove it in step 1.")
                .color(t.error),
        );
        return;
    };
    let m = &doc.molecule;
    ui.label(RichText::new(format!("Lane {}: {}", i + 2, m.name())).strong());
    ui.horizontal(|ui| {
        for (kind, label) in [
            (Kind::Uncut, "Uncut"),
            (Kind::Digest, "Digest"),
            (Kind::Pcr, "PCR"),
        ] {
            ui.selectable_value(&mut inputs.kind, kind, label);
        }
    });
    match inputs.kind {
        Kind::Uncut => {
            ui.label(
                RichText::new("The construct itself (supercoiled if circular).")
                    .small()
                    .color(t.text_muted),
            );
        }
        Kind::Digest => {
            ui.add(
                egui::TextEdit::singleline(&mut inputs.enzymes)
                    .hint_text("EcoRI, BamHI")
                    .desired_width(260.0),
            );
        }
        Kind::Pcr => {
            let primers: Vec<String> = m.primers().iter().map(|p| p.name.clone()).collect();
            if primers.is_empty() {
                ui.label(
                    RichText::new("This construct has no primers. Type primer names or sequences:")
                        .small()
                        .color(t.text_muted),
                );
            }
            for (label, id, value) in [
                ("Forward", "exp-fwd", &mut inputs.forward),
                ("Reverse", "exp-rev", &mut inputs.reverse),
            ] {
                ui.horizontal(|ui| {
                    ui.add_sized(
                        [56.0, 18.0],
                        egui::Label::new(RichText::new(label).color(t.text_muted)),
                    );
                    ui.add(egui::TextEdit::singleline(value).desired_width(140.0));
                    if !primers.is_empty() {
                        egui::ComboBox::from_id_salt(id)
                            .width(60.0)
                            .selected_text("")
                            .show_ui(ui, |ui| {
                                for p in &primers {
                                    ui.selectable_value(value, p.clone(), p);
                                }
                            });
                    }
                });
            }
        }
    }
    let method = inputs.method();
    ui.horizontal(|ui| {
        if ui
            .add_enabled(method.is_some(), egui::Button::new(format!("Apply to lane {}", i + 2)))
            .clicked()
            && let (Some(method), Some(lane)) = (&method, exp.lanes.get_mut(i))
        {
            lane.method = Some(method.clone());
        }
        if ui
            .add_enabled(method.is_some(), egui::Button::new("Apply to all lanes"))
            .on_hover_text("Give every construct this method (primers are looked up by name in each construct)")
            .clicked()
            && let Some(method) = &method
        {
            exp.apply_to_all(method);
        }
    });
}

/// "1 Constructs  ›  2 Lanes", with the current step highlighted.
/// Diagnostic Digest, step 2: enzymes whose digests tell the constructs apart, best first.
/// Returns params for `diagnostic_digest` when the suggestions must be (re)computed.
fn enzymes_step(
    ui: &mut Ui,
    t: &Tokens,
    session: &Session,
    exp: &mut Experiment,
    inputs: &mut ExperimentInputs,
) -> Option<serde_json::Value> {
    steps(ui, t, exp.kind, 3);
    let mut documents: Vec<u64> = Vec::new();
    for l in &exp.lanes {
        if !documents.contains(&l.document) && session.document_by_id(l.document).is_some() {
            documents.push(l.document);
        }
    }
    let single = documents.len() == 1;
    ui.label(
        RichText::new(if single {
            "Digests that give this construct a clear pattern of bands, to confirm it on a gel. \
             Best first."
        } else {
            "Digests that tell these constructs apart on a gel: every two differ by a clearly \
             visible band. Best first."
        })
        .color(t.text_muted),
    );
    ui.add_space(t.spacing);
    let mut back = false;
    ui.horizontal(|ui| {
        back = ui.button("Back to constructs").clicked();
        ui.separator();
        ui.label("Agarose");
        ui.add(
            egui::DragValue::new(&mut exp.agarose)
                .range(0.3..=3.0)
                .speed(0.05)
                .suffix(" %"),
        );
        let mut doubles = !inputs.singles_only;
        if ui
            .checkbox(&mut doubles, "Also double digests")
            .on_hover_text("Pairs of enzymes used together; they need a buffer both work in")
            .changed()
        {
            inputs.singles_only = !doubles;
        }
    });
    if let Some(pool) = inputs.enzyme_pool.clone() {
        ui.horizontal(|ui| {
            ui.checkbox(
                &mut inputs.include_out_of_stock,
                "Include enzymes I haven't got",
            )
            .on_hover_text(format!(
                "Off: only the enzymes of \"{pool}\", so every digest here is one you can \
                     run. On: every enzyme is searched, and a digest needing one that set \
                     hasn't got is marked and listed after the rest."
            ));
        });
    }
    if back {
        exp.step = ExperimentStep::Constructs;
        return None;
    }
    ui.add_space(t.spacing);
    let key = format!(
        "{}|{documents:?}|{}|{}|{:?}|{}",
        session.revision(),
        exp.agarose,
        inputs.singles_only,
        inputs.enzyme_pool,
        inputs.include_out_of_stock
    );
    let mut request = None;
    if inputs.suggest_key.as_deref() != Some(key.as_str()) {
        inputs.suggest_key = Some(key);
        inputs.suggestions = None;
        let mut params = serde_json::json!({
            "documents": documents,
            "agarose": exp.agarose,
            "max_enzymes": if inputs.singles_only { 1 } else { 2 },
            "limit": 20,
        });
        if let (Some(o), Some(pool)) = (params.as_object_mut(), &inputs.enzyme_pool) {
            o.insert("set".into(), serde_json::json!(pool));
            o.insert(
                "include_out_of_stock".into(),
                serde_json::json!(inputs.include_out_of_stock),
            );
        }
        request = Some(params);
    }
    let names: Vec<String> = documents
        .iter()
        .map(|d| name_of(session, *d).unwrap_or_default())
        .collect();
    let mut chosen: Option<Vec<String>> = None;
    match &inputs.suggestions {
        None => {
            ui.label(RichText::new("Searching\u{2026}").color(t.text_muted));
        }
        Some(Err(e)) => {
            ui.label(RichText::new(e).color(t.error));
        }
        Some(Ok(v)) => {
            if let Some(name) = v.get("set").and_then(serde_json::Value::as_str) {
                let held = v
                    .get("in_stock_searched")
                    .and_then(serde_json::Value::as_u64)
                    .unwrap_or(0);
                ui.label(
                    RichText::new(if inputs.include_out_of_stock {
                        format!("Every enzyme searched; \"{name}\" holds {held} of them.")
                    } else {
                        format!("Searching the {held} enzymes of \"{name}\".")
                    })
                    .small()
                    .color(t.text_muted),
                );
                let absent: Vec<&str> = v
                    .get("missing_from_database")
                    .and_then(serde_json::Value::as_array)
                    .into_iter()
                    .flatten()
                    .filter_map(serde_json::Value::as_str)
                    .collect();
                if !absent.is_empty() {
                    ui.label(
                        RichText::new(format!(
                            "Not in the enzyme database, so not searched: {}. Download REBASE \
                             (Preferences > Enzymes) to use them.",
                            absent.join(", ")
                        ))
                        .small()
                        .color(t.warning),
                    );
                }
            } else {
                ui.label(
                    RichText::new(
                        "Every commercial enzyme is searched. Choose an enzyme set of your own \
                         (Preferences > Enzymes) to be offered only digests you can run.",
                    )
                    .small()
                    .color(t.text_muted),
                );
            }
            ui.add_space(t.spacing);
            let list = v
                .get("suggestions")
                .and_then(serde_json::Value::as_array)
                .cloned()
                .unwrap_or_default();
            if list.is_empty() {
                ui.label(
                    RichText::new(
                        "No digest tells these constructs apart on a gel. If they differ by only \
                         a few bases, sequence them instead (Align tab).",
                    )
                    .color(t.text_muted),
                );
            }
            egui::ScrollArea::vertical()
                .id_salt("diagnostic-suggestions")
                .show(ui, |ui| {
                    for (i, sug) in list.iter().enumerate() {
                        let enzymes: Vec<String> = sug
                            .get("enzymes")
                            .and_then(|e| serde_json::from_value(e.clone()).ok())
                            .unwrap_or_default();
                        let alternatives: Vec<String> = sug
                            .get("alternatives")
                            .and_then(serde_json::Value::as_array)
                            .into_iter()
                            .flatten()
                            .filter_map(serde_json::Value::as_array)
                            .flatten()
                            .filter_map(serde_json::Value::as_str)
                            .map(str::to_string)
                            .collect();
                        ui.horizontal(|ui| {
                            if ui
                                .button("Use")
                                .on_hover_text(
                                    "Digest every lane with these enzymes and show the gel",
                                )
                                .clicked()
                            {
                                chosen = Some(enzymes.clone());
                            }
                            let r = ui.label(RichText::new(enzymes.join(" + ")).strong());
                            if !alternatives.is_empty() {
                                r.on_hover_text(format!(
                                    "Cut the same way: {}",
                                    alternatives.join(", ")
                                ));
                            }
                            if i == 0 {
                                ui.label(RichText::new("best").color(t.accent_text).small());
                            }
                            if sug.get("in_stock") == Some(&serde_json::json!(false)) {
                                ui.label(RichText::new("not in your set").color(t.warning).small())
                                    .on_hover_text(
                                        "You'd have to buy or borrow an enzyme for this digest.",
                                    );
                            }
                        });
                        let lanes = sug
                            .get("lanes")
                            .and_then(serde_json::Value::as_array)
                            .cloned()
                            .unwrap_or_default();
                        for (lane, name) in lanes.iter().zip(&names) {
                            let sizes: Vec<String> = lane
                                .get("fragments")
                                .and_then(serde_json::Value::as_array)
                                .into_iter()
                                .flatten()
                                .filter_map(serde_json::Value::as_u64)
                                .map(group_digits)
                                .collect();
                            ui.label(
                                RichText::new(format!("    {name}: {} bp", sizes.join(" + ")))
                                    .color(t.text_muted),
                            );
                        }
                        for note in sug
                            .get("notes")
                            .and_then(serde_json::Value::as_array)
                            .into_iter()
                            .flatten()
                            .filter_map(serde_json::Value::as_str)
                        {
                            ui.label(
                                RichText::new(format!("    {note}"))
                                    .small()
                                    .color(t.text_muted),
                            );
                        }
                        ui.add_space(t.spacing);
                    }
                });
        }
    }
    if let Some(enzymes) = chosen {
        for lane in &mut exp.lanes {
            lane.method = Some(LaneSpec::Digest {
                enzymes: enzymes.clone(),
            });
        }
        exp.step = ExperimentStep::Lanes;
        exp.selected = Some(0);
    }
    request
}

/// 4124 as "4,124".
fn group_digits(n: u64) -> String {
    let s = n.to_string();
    let mut out = String::new();
    for (i, c) in s.chars().enumerate() {
        if i > 0 && (s.len() - i).is_multiple_of(3) {
            out.push(',');
        }
        out.push(c);
    }
    out
}

fn steps(ui: &mut Ui, t: &Tokens, kind: ExperimentKind, current: u8) {
    let names: &[&str] = match kind {
        ExperimentKind::Gel => &["Simulation", "Constructs", "Lanes"],
        ExperimentKind::Diagnostic => &["Simulation", "Constructs", "Enzymes", "Gel"],
        ExperimentKind::Reads => &["Construct and reads", "Results"],
        // Cloning draws its own breadcrumb, which follows the technique.
        ExperimentKind::Cloning => &[],
    };
    ui.horizontal(|ui| {
        for (i, label) in names.iter().enumerate() {
            let n = u8::try_from(i + 1).unwrap_or(u8::MAX);
            if i > 0 {
                ui.label(RichText::new("\u{203a}").color(t.text_muted));
            }
            let text = RichText::new(format!("{n}  {label}"));
            ui.label(if n == current {
                text.strong().color(t.text)
            } else {
                text.color(t.text_muted)
            });
        }
    });
    ui.add_space(t.spacing);
}

#[cfg(test)]
mod tests {
    use snappiesnap_engine::Document;
    use snappiesnap_engine::doc::{Alphabet, Molecule, Seq, Topology};

    use super::*;

    #[test]
    fn editor_methods() {
        let mut inputs = ExperimentInputs::default();
        assert_eq!(inputs.method(), Some(LaneSpec::Uncut));
        inputs.kind = Kind::Digest;
        assert_eq!(inputs.method(), None);
        inputs.enzymes = "EcoRI + BamHI,".into();
        let digest = inputs.method().unwrap();
        assert_eq!(digest.title(), "EcoRI + BamHI");
        inputs.kind = Kind::Pcr;
        inputs.forward = " f ".into();
        assert_eq!(inputs.method(), None);
        inputs.reverse = "r".into();
        assert_eq!(
            inputs.method(),
            Some(LaneSpec::Pcr {
                forward: "f".into(),
                reverse: "r".into()
            })
        );
        // Loading a lane fills the editor; a lane without a method keeps it.
        inputs.load(Some(&digest));
        assert_eq!(
            (inputs.kind, inputs.enzymes.as_str()),
            (Kind::Digest, "EcoRI, BamHI")
        );
        inputs.load(None);
        assert_eq!(inputs.kind, Kind::Digest);
    }

    /// Both steps draw with two constructs; set `SNAPPIESNAP_EXPERIMENT_PNG=/path.png` to look
    /// at them.
    #[test]
    fn draws_both_steps() {
        let mut session =
            Session::new().with_platform(std::sync::Arc::new(snappiesnap_engine::MemoryFs::new()));
        let alpha = format!("{0}GAATTC{1}GGATCC{0}", "A".repeat(800), "C".repeat(1500));
        let beta = "ACGTTGCA".repeat(400);
        for (name, residues) in [("pAlpha", alpha), ("pBeta", beta)] {
            let m = Molecule::new(
                name,
                Seq::new(Alphabet::Dna, residues).unwrap(),
                Topology::Circular,
            )
            .unwrap();
            session.open_document(Document::new(m));
        }
        let mut exp = Experiment::default();
        let mut inputs = ExperimentInputs::default();
        let t = Tokens::dark();
        let render = |exp: &mut Experiment, inputs: &mut ExperimentInputs, step: &str| {
            let image = crate::offscreen::render(
                1000,
                560,
                1.0,
                |ctx| t.apply(ctx),
                |ui| {
                    egui::Frame::new().fill(t.panel).show(ui, |ui| {
                        let mut pane = crate::cloning_experiment::Pane {
                            sequence: false,
                            layers: crate::map::Layers::default(),
                            style: crate::state::MapStyle::Auto,
                            translation: crate::state::Translation::Orfs,
                            genetic_code: 1,
                            tails: crate::state::TailStyle::default(),
                            single_line: false,
                        };
                        show(ui, &t, &session, exp, inputs, None, &mut pane);
                    });
                },
            )
            .unwrap();
            if let Ok(path) = std::env::var("SNAPPIESNAP_EXPERIMENT_PNG") {
                let path = path.replace(".png", &format!("_{step}.png"));
                crate::offscreen::write_png(&path, &image).unwrap();
            }
        };
        render(&mut exp, &mut inputs, "start");
        exp.start(session.documents().iter().map(|d| d.id));
        exp.add(exp.lanes[0].document);
        render(&mut exp, &mut inputs, "constructs");
        exp.step = ExperimentStep::Lanes;
        exp.selected = Some(0);
        exp.apply_to_all(&LaneSpec::Digest {
            enzymes: vec!["EcoRI".into(), "BamHI".into()],
        });
        exp.lanes[2].method = Some(LaneSpec::Uncut);
        render(&mut exp, &mut inputs, "lanes");
        let run = inputs.lanes(&session, &exp).to_vec();
        assert_eq!(run.len(), 4);
        assert_eq!(run[1].0, "pAlpha: EcoRI + BamHI");
        assert_eq!(run[2].1.as_ref().unwrap().len(), 1);
        assert!(run[3].0.ends_with("Uncut"));
        assert_eq!(inputs.loaded, Some(0));
        assert_eq!(inputs.kind, Kind::Digest);
    }
    /// The Diagnostic Digest step asks for the enzymes the lab holds, and only widens the
    /// search when told to. Set `SNAPPIESNAP_DIGEST_PNG=/path.png` to look at it.
    #[test]
    fn diagnostic_digest_searches_the_lab_set() {
        let mut session =
            Session::new().with_platform(std::sync::Arc::new(snappiesnap_engine::MemoryFs::new()));
        let residues = format!("{0}GAATTC{1}GGATCC{0}", "A".repeat(800), "C".repeat(1500));
        let m = Molecule::new(
            "pAlpha",
            Seq::new(Alphabet::Dna, residues).unwrap(),
            Topology::Circular,
        )
        .unwrap();
        session.open_document(Document::new(m));
        let mut exp = Experiment {
            kind: ExperimentKind::Diagnostic,
            ..Default::default()
        };
        exp.start(session.documents().iter().map(|d| d.id));
        exp.step = ExperimentStep::Enzymes;
        let t = Tokens::dark();
        let mut inputs = ExperimentInputs::default();
        let ask = |exp: &mut Experiment, inputs: &mut ExperimentInputs| {
            let mut params = None;
            let image = crate::offscreen::render(
                1000,
                560,
                1.0,
                |ctx| t.apply(ctx),
                |ui| {
                    egui::Frame::new().fill(t.panel).show(ui, |ui| {
                        // Three passes settle the layout; the question is asked on the first.
                        if let Some(p) = enzymes_step(ui, &t, &session, exp, inputs) {
                            params = Some(p);
                        }
                    });
                },
            )
            .unwrap();
            (params, image)
        };

        // With no set chosen, nothing is restricted: the params say nothing about a set.
        let (params, _) = ask(&mut exp, &mut inputs);
        let params = params.expect("the step asks for suggestions");
        assert!(params.get("set").is_none(), "{params}");
        assert!(params.get("include_out_of_stock").is_none(), "{params}");

        // With one, the search is the set's enzymes, and every suggestion is runnable.
        inputs.enzyme_pool = Some("Lab freezer".into());
        inputs.suggest_key = None;
        let (params, image) = ask(&mut exp, &mut inputs);
        let params = params.expect("a new pool is a new question");
        assert_eq!(params["set"], "Lab freezer");
        assert_eq!(params["include_out_of_stock"], serde_json::json!(false));
        if let Ok(path) = std::env::var("SNAPPIESNAP_DIGEST_PNG") {
            crate::offscreen::write_png(&path, &image).unwrap();
        }

        // Asking again with the same question doesn't re-run the search.
        let (again, _) = ask(&mut exp, &mut inputs);
        assert!(again.is_none(), "{again:?}");

        // Widening is a different question, and says so.
        inputs.include_out_of_stock = true;
        let (params, _) = ask(&mut exp, &mut inputs);
        let params = params.expect("widening re-runs the search");
        assert_eq!(params["include_out_of_stock"], serde_json::json!(true));
        assert_eq!(params["set"], "Lab freezer");
    }
}
