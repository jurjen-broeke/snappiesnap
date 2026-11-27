//! UI state that isn't document state. Serialisable so the control channel can read and set it
//! (`ui.state`, `ui.set`).

use eframe::egui;
use serde::{Deserialize, Serialize};

use crate::theme::ThemeChoice;

/// The main views, as tabs above the canvas.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum View {
    #[default]
    Map,
    Sequence,
    Features,
    Primers,
    Enzymes,
    Gel,
    Align,
    History,
    Notes,
    /// The chromatogram of a sequencing read (shown only for reads).
    Trace,
}

impl View {
    pub const ALL: [View; 10] = [
        View::Map,
        View::Sequence,
        View::Trace,
        View::Features,
        View::Primers,
        View::Enzymes,
        View::Gel,
        View::Align,
        View::History,
        View::Notes,
    ];

    pub fn label(self) -> &'static str {
        match self {
            View::Map => "Map",
            View::Sequence => "Sequence",
            View::Features => "Features",
            View::Primers => "Primers",
            View::Enzymes => "Enzymes",
            View::Gel => "Gel",
            View::Align => "Align",
            View::History => "History",
            View::Notes => "Notes",
            View::Trace => "Trace",
        }
    }
}

/// The simulated gels: one for restriction digests, one for PCR, each with its own agarose,
/// ladder and lanes; `mode` says which one the Gel tab shows.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct GelSetup {
    pub mode: GelMode,
    pub digest: GelLanes,
    pub pcr: GelLanes,
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum GelMode {
    #[default]
    Digest,
    Pcr,
}

/// One gel: agarose %, the ladder in the first lane, and the other lanes.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct GelLanes {
    pub agarose: f64,
    pub ladder: String,
    pub lanes: Vec<snappiesnap_engine::gel::LaneSpec>,
}

impl Default for GelLanes {
    fn default() -> Self {
        GelLanes {
            agarose: 1.0,
            ladder: "1 kb".to_string(),
            lanes: Vec::new(),
        }
    }
}

impl Default for GelSetup {
    fn default() -> Self {
        GelSetup {
            mode: GelMode::Digest,
            digest: GelLanes {
                lanes: vec![snappiesnap_engine::gel::LaneSpec::Uncut],
                ..GelLanes::default()
            },
            // PCR products are usually small: a denser gel and a finer ladder.
            pcr: GelLanes {
                agarose: 1.5,
                ladder: "100 bp".to_string(),
                lanes: Vec::new(),
            },
        }
    }
}

impl GelSetup {
    /// The gel the Gel tab shows.
    pub fn current(&mut self) -> &mut GelLanes {
        match self.mode {
            GelMode::Digest => &mut self.digest,
            GelMode::Pcr => &mut self.pcr,
        }
    }
}

/// An experiment, opened from the toolbar's Experiments menu in a tab of its own beside the
/// document tabs: it spans documents, so it belongs to none of them. For now a gel experiment:
/// one gel across constructs. The user first chooses and orders the constructs (one lane each;
/// lane 1 is the ladder), then sets a method per lane.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct Experiment {
    /// The tab's name.
    pub name: String,
    pub kind: ExperimentKind,
    pub step: ExperimentStep,
    pub agarose: f64,
    pub ladder: String,
    /// One lane per construct, in lane order (lanes 2, 3, ... of the gel).
    pub lanes: Vec<snappiesnap_engine::gel::ExperimentLane>,
    /// The lane being set up (an index into `lanes`).
    pub selected: Option<usize>,
    /// Check Reads: the construct (an open document by tab id) and the read files.
    pub reference: Option<u64>,
    pub reads: Vec<String>,
    /// Check Reads: reads that are open in tabs, by id. Not saved with the experiment: ids
    /// mean nothing in another session, and a wrong one would check the wrong read.
    #[serde(skip)]
    pub read_documents: Vec<u64>,
    /// Check Reads: also check the reads stored with the construct in its .dna file.
    pub matched: bool,
    /// Cloning: the technique, and the two sides of the reaction.
    #[serde(default)]
    pub technique: CloningTechnique,
    #[serde(default)]
    pub vector: CloningSide,
    #[serde(default)]
    pub insert: CloningSide,
    /// Assembly: parts after the second, in order.
    #[serde(default)]
    pub extra: Vec<CloningSide>,
    /// The file the experiment was opened from or last saved to.
    pub path: Option<String>,
    /// [`Experiment::content_key`] when it was last opened or saved, to tell unsaved changes.
    #[serde(skip)]
    pub saved: Option<String>,
}

/// The kinds of experiment the toolbar's Experiments menu offers (more designs will follow).
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum ExperimentKind {
    /// One gel across constructs.
    #[default]
    Gel,
    /// Enzymes whose digests tell constructs apart, then their gel.
    Diagnostic,
    /// A construct checked with sequencing reads.
    Reads,
    /// Two constructs cloned together, step by step.
    Cloning,
}

impl ExperimentKind {
    /// Every kind, in menu order.
    pub const ALL: [ExperimentKind; 4] = [
        ExperimentKind::Gel,
        ExperimentKind::Diagnostic,
        ExperimentKind::Reads,
        ExperimentKind::Cloning,
    ];

    /// The kinds the Experiments menu lists. Both gel kinds live behind Simulate Gel, which
    /// asks in its first step which gel to simulate, so [`ExperimentKind::Diagnostic`] has no
    /// menu item of its own.
    pub const MENU: [ExperimentKind; 3] = [
        ExperimentKind::Gel,
        ExperimentKind::Cloning,
        ExperimentKind::Reads,
    ];

    /// The kinds Simulate Gel offers in its first step, in the order shown.
    pub const SIMULATIONS: [ExperimentKind; 2] = [ExperimentKind::Gel, ExperimentKind::Diagnostic];

    /// Whether this kind is one of [`ExperimentKind::SIMULATIONS`], so it begins by asking
    /// which gel to simulate.
    pub fn chooses_simulation(self) -> bool {
        Self::SIMULATIONS.contains(&self)
    }

    pub fn label(self) -> &'static str {
        match self {
            ExperimentKind::Gel => "Simulate Gel",
            ExperimentKind::Diagnostic => "Diagnostic Digest",
            ExperimentKind::Reads => "Check Reads",
            ExperimentKind::Cloning => "Cloning",
        }
    }

    /// One line for the menu.
    pub fn description(self) -> &'static str {
        match self {
            ExperimentKind::Gel => {
                "Run constructs side by side on one gel, or find the digest that tells them apart"
            }
            ExperimentKind::Diagnostic => {
                "Find enzymes whose digests tell your constructs apart on a gel, such as a clone from its empty vector"
            }
            ExperimentKind::Reads => {
                "Check a construct with sequencing reads: what they cover and every difference, with the traces"
            }
            ExperimentKind::Cloning => {
                "Clone two constructs together step by step: the technique, the vector, the insert, the product"
            }
        }
    }

    /// What a new experiment of this kind is called in its tab (numbered when taken).
    pub fn tab_base(self) -> &'static str {
        match self {
            ExperimentKind::Gel => "Gel experiment",
            ExperimentKind::Diagnostic => "Diagnostic digest",
            ExperimentKind::Reads => "Read check",
            ExperimentKind::Cloning => "Cloning",
        }
    }

    /// The name of this simulation in the first step, where both gel kinds are offered.
    pub fn simulation_label(self) -> &'static str {
        match self {
            ExperimentKind::Gel => "Gel of several constructs",
            ExperimentKind::Diagnostic => "Diagnostic digest",
            ExperimentKind::Reads => "Check reads",
            ExperimentKind::Cloning => "Cloning",
        }
    }
}

/// A cloning technique the engine simulates, offered in step 1 of a cloning experiment.
/// Each knows what its two sides are called and which command runs it, so the later steps
/// follow from the choice.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum CloningTechnique {
    #[default]
    Restriction,
    GatewayLr,
    GatewayBp,
    Gibson,
    Infusion,
    GoldenGate,
    TopoTa,
    TopoBlunt,
    TopoDirectional,
    Mutagenesis,
}

impl CloningTechnique {
    /// Every technique, in the order step 1 lists them.
    pub const ALL: [CloningTechnique; 10] = [
        CloningTechnique::Restriction,
        CloningTechnique::GatewayLr,
        CloningTechnique::GatewayBp,
        CloningTechnique::Gibson,
        CloningTechnique::Infusion,
        CloningTechnique::GoldenGate,
        CloningTechnique::TopoTa,
        CloningTechnique::TopoBlunt,
        CloningTechnique::TopoDirectional,
        CloningTechnique::Mutagenesis,
    ];

    pub fn label(self) -> &'static str {
        match self {
            CloningTechnique::Restriction => "Restriction cloning",
            CloningTechnique::GatewayLr => "Gateway LR",
            CloningTechnique::GatewayBp => "Gateway BP",
            CloningTechnique::Gibson => "Gibson / HiFi",
            CloningTechnique::Infusion => "In-Fusion",
            CloningTechnique::GoldenGate => "Golden Gate",
            CloningTechnique::TopoTa => "TA / TOPO TA",
            CloningTechnique::TopoBlunt => "Blunt TOPO",
            CloningTechnique::TopoDirectional => "Directional TOPO",
            CloningTechnique::Mutagenesis => "Site-directed mutagenesis",
        }
    }

    pub fn description(self) -> &'static str {
        match self {
            CloningTechnique::Restriction => {
                "Cut both constructs with restriction enzymes and ligate the fragments whose ends fit."
            }
            CloningTechnique::GatewayLr => {
                "An entry clone (attL) into a destination vector (attR), giving an expression clone with attB sites."
            }
            CloningTechnique::GatewayBp => {
                "An attB insert into a donor vector (attP), giving an entry clone with attL sites."
            }
            CloningTechnique::Gibson => {
                "Join pieces whose ends match over at least 20 bp, in one isothermal reaction."
            }
            CloningTechnique::Infusion => "Join pieces whose ends match over at least 15 bp.",
            CloningTechnique::GoldenGate => {
                "Cut both parts with one Type IIS enzyme (BsaI by default) and chain the fragments by their overhangs."
            }
            CloningTechnique::TopoTa => "A PCR product with 3' A overhangs into a T-vector.",
            CloningTechnique::TopoBlunt => "A blunt-ended insert into a blunt TOPO vector.",
            CloningTechnique::TopoDirectional => {
                "A CACC-tailed insert into a directional TOPO vector, which goes in one way round."
            }
            CloningTechnique::Mutagenesis => {
                "Copy one plasmid with back-to-back primers to change, insert or delete bases. It has no second construct: step 3 takes the primers."
            }
        }
    }

    /// What step 2 is called for this technique.
    pub fn vector_label(self) -> &'static str {
        match self {
            CloningTechnique::GatewayLr => "Destination vector",
            CloningTechnique::GatewayBp => "Donor vector",
            CloningTechnique::Gibson
            | CloningTechnique::Infusion
            | CloningTechnique::GoldenGate => "First part (the vector)",
            CloningTechnique::Mutagenesis => "Plasmid",
            _ => "Vector",
        }
    }

    /// What step 3 is called for this technique.
    pub fn insert_label(self) -> &'static str {
        match self {
            CloningTechnique::GatewayLr => "Entry clone",
            CloningTechnique::GatewayBp => "attB insert",
            CloningTechnique::Gibson
            | CloningTechnique::Infusion
            | CloningTechnique::GoldenGate => "Second part",
            CloningTechnique::Mutagenesis => "Primers",
            _ => "Insert",
        }
    }

    /// Whether step 3 takes primers instead of a second construct.
    pub fn insert_is_primers(self) -> bool {
        self == CloningTechnique::Mutagenesis
    }

    /// Whether this technique is one of the assemblies, whose parts may be more than two.
    pub fn is_assembly(self) -> bool {
        matches!(
            self,
            CloningTechnique::Gibson | CloningTechnique::Infusion | CloningTechnique::GoldenGate
        )
    }

    /// Whether a side may be cut or PCR-amplified before the reaction.
    pub fn prepares_sides(self) -> bool {
        self.is_assembly()
    }

    /// The command that runs it.
    pub fn command(self) -> &'static str {
        match self {
            CloningTechnique::Restriction => "restriction_clone",
            CloningTechnique::GatewayLr | CloningTechnique::GatewayBp => "gateway_clone",
            CloningTechnique::Gibson
            | CloningTechnique::Infusion
            | CloningTechnique::GoldenGate => "assemble",
            CloningTechnique::TopoTa
            | CloningTechnique::TopoBlunt
            | CloningTechnique::TopoDirectional => "topo_clone",
            CloningTechnique::Mutagenesis => "mutagenesis",
        }
    }
}

/// How a construct of a cloning experiment is used.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum SideUse {
    /// As it is in the file.
    #[default]
    Whole,
    /// Cut with enzymes first, keeping the largest fragment.
    Cut,
    /// The PCR product of two primers.
    Pcr,
}

/// One side of a cloning experiment: which construct, and how it is prepared.
#[derive(Debug, Clone, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct CloningSide {
    /// The file it is read from; empty until one is chosen.
    #[serde(default)]
    pub file: String,
    /// An open document's tab id, when the construct is one of the open tabs rather than a
    /// file: this is how a construct is chosen in the web build, which has no file paths of
    /// its own, and how an edited or unsaved tab is cloned from. Never saved in a `.sslab`
    /// file, where a tab id would mean nothing: `file` is kept alongside it when the document
    /// has one, so a reopened experiment falls back to the file.
    #[serde(default, skip)]
    pub document: Option<u64>,
    #[serde(default)]
    pub using: SideUse,
    /// Enzymes, as typed (comma, space or + separated).
    #[serde(default)]
    pub enzymes: String,
    /// Primers: a name of the construct's own, or a sequence. For mutagenesis these are the
    /// two back-to-back primers.
    #[serde(default)]
    pub forward: String,
    #[serde(default)]
    pub reverse: String,
    /// TOPO: where a circular vector is opened (1-based, as typed; empty: already linear).
    #[serde(default)]
    pub position: String,
}

impl CloningSide {
    /// Whether a construct has been chosen.
    pub fn is_set(&self) -> bool {
        self.document.is_some() || !self.file.trim().is_empty()
    }

    /// What identifies this construct while the app runs: the open document, or the file.
    /// An open document wins, because that is what the user picked and it holds their edits.
    pub fn source_key(&self) -> String {
        match self.document {
            Some(id) => format!("doc:{id}"),
            None => self.file.trim().to_string(),
        }
    }

    /// The file's last path element, for showing.
    pub fn file_name(&self) -> &str {
        self.file.rsplit(['/', '\\']).next().unwrap_or(&self.file)
    }
}

/// Where the experiment is.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum ExperimentStep {
    /// Nothing set up yet.
    #[default]
    NotStarted,
    /// Choosing which simulation to run (the first step of Simulate Gel).
    Simulation,
    /// Choosing and ordering the constructs.
    Constructs,
    /// Choosing enzymes from suggestions (Diagnostic Digest).
    Enzymes,
    /// The reads against the construct (Check Reads).
    Results,
    /// Setting the method of each lane.
    Lanes,
    /// Cloning: choosing the technique (step 1).
    Technique,
    /// Cloning: the vector and how it is prepared (step 2).
    Vector,
    /// Cloning: the insert, or the primers, and how it is prepared (step 3).
    Insert,
    /// Cloning: the product, and opening it as a construct (step 4).
    Product,
}

impl Default for Experiment {
    fn default() -> Self {
        Experiment {
            name: "Gel experiment".to_string(),
            kind: ExperimentKind::Gel,
            step: ExperimentStep::NotStarted,
            agarose: 1.0,
            ladder: "1 kb".to_string(),
            lanes: Vec::new(),
            selected: None,
            reference: None,
            reads: Vec::new(),
            read_documents: Vec::new(),
            matched: false,
            technique: CloningTechnique::default(),
            vector: CloningSide::default(),
            insert: CloningSide::default(),
            extra: Vec::new(),
            path: None,
            saved: None,
        }
    }
}

impl Experiment {
    /// Most constructs (the ladder takes one more lane of the 50).
    pub const MAX_LANES: usize = 49;

    /// Longest experiment name, in characters.
    pub const MAX_NAME: usize = 80;

    /// A new experiment of `kind` called `name`, with these documents (tab ids) as the
    /// constructs, in this order.
    pub fn new(kind: ExperimentKind, name: &str, documents: impl IntoIterator<Item = u64>) -> Self {
        let mut e = Experiment {
            name: name.to_string(),
            kind,
            ..Experiment::default()
        };
        if kind == ExperimentKind::Cloning {
            // Cloning is bound to no construct: it starts by choosing the technique, and its
            // two constructs are chosen from files in the steps after it.
            e.step = ExperimentStep::Technique;
            return e;
        }
        e.start(documents);
        // Simulate Gel asks first which gel to simulate; the constructs it seeded are kept
        // and shown in the next step.
        if kind.chooses_simulation() {
            e.step = ExperimentStep::Simulation;
        }
        e
    }

    /// What a saved file holds: name, gel and lanes.
    pub fn content_key(&self) -> String {
        serde_json::to_string(&(
            &self.name,
            self.agarose,
            &self.ladder,
            &self.lanes,
            self.reference,
            &self.reads,
            self.matched,
        ))
        .unwrap_or_default()
    }

    /// Whether what is shown differs from the file (or was never saved).
    pub fn unsaved(&self) -> bool {
        self.saved.as_deref() != Some(self.content_key().as_str())
    }

    /// The params of `experiment_save` for this experiment at `path`.
    pub fn save_params(&self, path: &str) -> serde_json::Value {
        if self.kind == ExperimentKind::Reads {
            return serde_json::json!({
                "path": path,
                "name": self.name,
                "kind": "reads",
                "reference": self.reference,
                "reads": self.reads,
                "matched": self.matched,
            });
        }
        let mut v = serde_json::json!({
            "path": path,
            "name": self.name,
            "agarose": self.agarose,
            "ladder": self.ladder,
            "lanes": self.lanes,
        });
        if self.kind == ExperimentKind::Diagnostic
            && let Some(o) = v.as_object_mut()
        {
            o.insert("design".into(), "diagnostic_digest".into());
        }
        v
    }

    /// An experiment from the result of `experiment_open`, at its step of setting methods.
    pub fn from_opened(path: &str, opened: &serde_json::Value) -> Option<Experiment> {
        let lanes: Vec<snappiesnap_engine::gel::ExperimentLane> =
            serde_json::from_value(opened.get("lanes")?.clone()).ok()?;
        let mut e = Experiment {
            name: opened
                .get("name")?
                .as_str()?
                .chars()
                .take(Self::MAX_NAME)
                .collect(),
            kind: match (
                opened.get("kind").and_then(serde_json::Value::as_str),
                opened.get("design").and_then(serde_json::Value::as_str),
            ) {
                (Some("reads"), _) => ExperimentKind::Reads,
                (_, Some("diagnostic_digest")) => ExperimentKind::Diagnostic,
                _ => ExperimentKind::Gel,
            },
            step: ExperimentStep::Lanes,
            agarose: opened.get("agarose")?.as_f64()?,
            ladder: opened.get("ladder")?.as_str()?.to_string(),
            lanes: lanes.into_iter().take(Self::MAX_LANES).collect(),
            selected: None,
            reference: opened.get("reference").and_then(serde_json::Value::as_u64),
            reads: opened
                .get("reads")
                .and_then(|r| serde_json::from_value(r.clone()).ok())
                .unwrap_or_default(),
            read_documents: Vec::new(),
            matched: opened.get("matched").and_then(serde_json::Value::as_bool) == Some(true),
            technique: CloningTechnique::default(),
            vector: CloningSide::default(),
            insert: CloningSide::default(),
            extra: Vec::new(),
            path: Some(path.to_string()),
            saved: None,
        };
        if e.kind == ExperimentKind::Reads {
            e.step = ExperimentStep::Results;
        }
        e.saved = Some(e.content_key());
        Some(e)
    }

    /// Choose which gel to simulate, in the first step. The tab keeps a name of the user's
    /// own; one still left as the default follows the choice.
    pub fn set_simulation(&mut self, kind: ExperimentKind) {
        if self.kind == kind {
            return;
        }
        let was_default = self
            .name
            .trim()
            .strip_prefix(self.kind.tab_base())
            .is_some_and(|rest| {
                rest.trim().is_empty() || rest.trim().chars().all(|c| c.is_ascii_digit())
            });
        if was_default {
            let number = self
                .name
                .trim()
                .strip_prefix(self.kind.tab_base())
                .map(str::trim)
                .filter(|rest| !rest.is_empty());
            self.name = match number {
                Some(n) => format!("{} {n}", kind.tab_base()),
                None => kind.tab_base().to_string(),
            };
        }
        self.kind = kind;
    }

    /// Start over with these documents (tab ids) as the constructs, in this order.
    /// Check Reads keeps its read files and takes the first document as its construct (unless
    /// it has one) instead of making lanes.
    pub fn start(&mut self, documents: impl IntoIterator<Item = u64>) {
        let mut documents = documents.into_iter();
        let reads = self.kind == ExperimentKind::Reads;
        let reference = if reads {
            self.reference.or_else(|| documents.next())
        } else {
            None
        };
        let lanes = if reads {
            Vec::new()
        } else {
            documents
                .take(Self::MAX_LANES)
                .map(|document| snappiesnap_engine::gel::ExperimentLane {
                    document,
                    method: None,
                })
                .collect()
        };
        *self = Experiment {
            name: std::mem::take(&mut self.name),
            kind: self.kind,
            step: ExperimentStep::Constructs,
            agarose: self.agarose,
            ladder: self.ladder.clone(),
            lanes,
            selected: None,
            reference,
            reads: std::mem::take(&mut self.reads),
            read_documents: std::mem::take(&mut self.read_documents),
            // A new read check includes the construct's matched reads, if it has any.
            matched: self.matched || reads,
            technique: self.technique,
            vector: std::mem::take(&mut self.vector),
            insert: std::mem::take(&mut self.insert),
            extra: std::mem::take(&mut self.extra),
            path: self.path.take(),
            saved: self.saved.take(),
        };
    }

    /// Add a construct as the last lane (the same construct may run in several lanes).
    pub fn add(&mut self, document: u64) {
        if self.lanes.len() < Self::MAX_LANES {
            self.lanes.push(snappiesnap_engine::gel::ExperimentLane {
                document,
                method: None,
            });
        }
    }

    /// Move lane `i` one place up (towards lane 1) or down; its method moves with it.
    pub fn move_lane(&mut self, i: usize, up: bool) {
        let j = if up {
            i.checked_sub(1)
        } else {
            i.checked_add(1)
        };
        if let Some(j) = j
            && j < self.lanes.len()
            && i < self.lanes.len()
        {
            self.lanes.swap(i, j);
            if self.selected == Some(i) {
                self.selected = Some(j);
            } else if self.selected == Some(j) {
                self.selected = Some(i);
            }
        }
    }

    pub fn remove(&mut self, i: usize) {
        if i < self.lanes.len() {
            self.lanes.remove(i);
        }
        self.selected = match self.selected {
            Some(s) if s > i => Some(s - 1),
            Some(s) if s == i || s >= self.lanes.len() => None,
            other => other,
        };
    }

    /// Give every lane this method.
    pub fn apply_to_all(&mut self, method: &snappiesnap_engine::gel::LaneSpec) {
        for lane in &mut self.lanes {
            lane.method = Some(method.clone());
        }
    }

    /// Keep the values in range (for `ui.set`).
    pub fn clamped(mut self) -> Self {
        self.lanes.truncate(Self::MAX_LANES);
        if !self.agarose.is_finite() {
            self.agarose = 1.0;
        }
        self.agarose = self.agarose.clamp(0.3, 2.0);
        self.selected = self.selected.filter(|&s| s < self.lanes.len());
        self.name = self.name.trim().chars().take(Self::MAX_NAME).collect();
        if self.name.is_empty() {
            self.name = self.kind.label().to_string();
        }
        self
    }
}

/// The background of the map and sequence views, per theme, as `#rrggbb` (`None`: the
/// theme's own). Set in Preferences.
#[derive(Debug, Clone, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct CanvasColors {
    pub light: Option<String>,
    pub dark: Option<String>,
}

impl CanvasColors {
    /// The colour for the light or dark theme, if one is set.
    pub fn for_theme(&self, dark: bool) -> Option<egui::Color32> {
        let c = if dark { &self.dark } else { &self.light };
        c.as_deref().and_then(crate::theme::parse_hex)
    }

    /// Set (or clear, with `None`) the colour for one theme.
    pub fn set(&mut self, dark: bool, color: Option<egui::Color32>) {
        let value = color.map(crate::theme::hex);
        if dark {
            self.dark = value;
        } else {
            self.light = value;
        }
    }

    /// Drop values that aren't colours (from `ui.set` or a hand-edited settings file).
    pub fn clamped(self) -> Self {
        let keep =
            |c: Option<String>| c.and_then(|c| crate::theme::parse_hex(&c).map(crate::theme::hex));
        CanvasColors {
            light: keep(self.light),
            dark: keep(self.dark),
        }
    }
}

/// Most experiment tabs open at once.
pub const MAX_EXPERIMENTS: usize = 20;

impl UiState {
    /// The experiment tab being shown, if any.
    pub fn experiment(&mut self) -> Option<&mut Experiment> {
        self.active_experiment
            .and_then(|i| self.experiments.get_mut(i))
    }

    /// Open a new experiment of `kind` in a tab of its own, with these documents (tab ids) as
    /// the constructs, and show it. Named "Gel experiment", "Gel experiment 2", ...
    pub fn new_experiment(
        &mut self,
        kind: ExperimentKind,
        documents: impl IntoIterator<Item = u64>,
    ) {
        if self.experiments.len() >= MAX_EXPERIMENTS {
            self.status = Some(format!(
                "At most {MAX_EXPERIMENTS} experiments can be open; close one first"
            ));
            return;
        }
        let base = kind.tab_base();
        let name = (1..)
            .map(|n| {
                if n == 1 {
                    base.to_string()
                } else {
                    format!("{base} {n}")
                }
            })
            .find(|name| !self.experiments.iter().any(|e| e.name == *name))
            .unwrap_or_else(|| base.to_string());
        self.experiments
            .push(Experiment::new(kind, &name, documents));
        self.active_experiment = Some(self.experiments.len() - 1);
    }

    /// Close an experiment tab. The next tab (or the previous) is shown if this one was.
    pub fn close_experiment(&mut self, i: usize) {
        if i >= self.experiments.len() {
            return;
        }
        self.experiments.remove(i);
        self.active_experiment = match self.active_experiment {
            Some(a) if a == i => {
                (!self.experiments.is_empty()).then(|| i.min(self.experiments.len() - 1))
            }
            Some(a) if a > i => Some(a - 1),
            other => other,
        };
    }
}

/// How the sequence view shows a primer's 5' tail (the part that doesn't pair).
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum TailStyle {
    /// A small loop marked with the tail's length.
    #[default]
    Simplified,
    /// The tail's bases written out over the sequence.
    Full,
}

/// How the map view draws the molecule.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum MapStyle {
    /// Circular for circular molecules, linear for linear ones.
    #[default]
    Auto,
    Circular,
    Linear,
}

/// What the ORF layer shows in the sequence view: ORFs, or the translation of one frame or a
/// set of frames. (The map always shows ORFs.)
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Translation {
    #[default]
    #[serde(alias = "none")]
    Orfs,
    #[serde(alias = "first")]
    Plus1,
    Plus2,
    Plus3,
    Minus1,
    Minus2,
    Minus3,
    /// Frames +1, +2, +3.
    #[serde(alias = "three")]
    Top,
    /// Frames -1, -2, -3.
    Bottom,
    /// All six frames.
    Six,
}

impl Translation {
    /// Every choice, in menu order, with its label; groups are separated by `None`.
    pub const MENU: [Option<(Translation, &'static str)>; 14] = [
        Some((Translation::Orfs, "ORFs only")),
        None,
        Some((Translation::Plus1, "Frame +1")),
        Some((Translation::Plus2, "Frame +2")),
        Some((Translation::Plus3, "Frame +3")),
        None,
        Some((Translation::Minus1, "Frame \u{2212}1")),
        Some((Translation::Minus2, "Frame \u{2212}2")),
        Some((Translation::Minus3, "Frame \u{2212}3")),
        None,
        Some((Translation::Top, "Top strand: 3 frames")),
        Some((Translation::Bottom, "Bottom strand: 3 frames")),
        Some((Translation::Six, "All 6 frames")),
        None,
    ];

    /// The frames translated: 1..=3 forward, -1..=-3 reverse (none for ORFs only).
    pub fn frames(self) -> &'static [i8] {
        match self {
            Translation::Orfs => &[],
            Translation::Plus1 => &[1],
            Translation::Plus2 => &[2],
            Translation::Plus3 => &[3],
            Translation::Minus1 => &[-1],
            Translation::Minus2 => &[-2],
            Translation::Minus3 => &[-3],
            Translation::Top => &[1, 2, 3],
            Translation::Bottom => &[-1, -2, -3],
            Translation::Six => &[1, 2, 3, -1, -2, -3],
        }
    }

    pub fn label(self) -> &'static str {
        Translation::MENU
            .iter()
            .flatten()
            .find(|(t, _)| *t == self)
            .map_or("", |(_, l)| l)
    }

    /// Whether the chosen frames are on the bottom strand.
    pub fn reverse(self) -> bool {
        matches!(
            self,
            Translation::Minus1 | Translation::Minus2 | Translation::Minus3 | Translation::Bottom
        )
    }
}

/// Largest map zoom.
pub const MAX_ZOOM: f32 = 50.0;

/// Which restriction enzymes the map and the sequence view show.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct EnzymeSet {
    /// Most sites an enzyme may have (1 = unique cutters; at most [`EnzymeSet::MAX_CUTS`]).
    pub max_cuts: usize,
    /// Fewest specific bases in the site (0 = any).
    pub min_length: usize,
    pub commercial_only: bool,
    /// Draw the enzymes from this enzyme set of the user's own (such as the lab's stock)
    /// instead of the whole database; the filters above still apply, except `commercial_only`.
    pub pool: Option<String>,
    /// Only these enzymes (when not empty, the filters above and `pool` don't apply).
    pub only: Vec<String>,
}

impl Default for EnzymeSet {
    /// Unique 6+ cutters that are sold.
    fn default() -> Self {
        EnzymeSet {
            max_cuts: 1,
            min_length: 6,
            commercial_only: true,
            pool: None,
            only: Vec::new(),
        }
    }
}

impl EnzymeSet {
    pub const MAX_CUTS: usize = 3;

    /// Keep the values in range.
    pub fn clamped(mut self) -> Self {
        self.max_cuts = self.max_cuts.clamp(1, Self::MAX_CUTS);
        self.min_length = self.min_length.min(12);
        self.only.truncate(500);
        for n in &mut self.only {
            *n = n.trim().chars().take(40).collect();
        }
        self.only.retain(|n| !n.is_empty());
        self.pool = self
            .pool
            .map(|p| {
                p.trim()
                    .chars()
                    .take(snappiesnap_engine::enzymes::MAX_SET_NAME)
                    .collect()
            })
            .filter(|p: &String| !p.is_empty());
        self
    }

    /// The engine query for this set.
    pub fn query(&self) -> snappiesnap_engine::enzymes::SiteQuery {
        snappiesnap_engine::enzymes::SiteQuery {
            enzymes: (!self.only.is_empty()).then(|| self.only.clone()),
            set: self.pool.clone(),
            min_length: self.min_length,
            max_cuts: self.max_cuts,
            commercial_only: self.commercial_only,
        }
    }

    /// In words: "unique 6+ cutters", "chosen enzymes (EcoRI, BamHI)", ...
    pub fn describe(&self) -> String {
        if !self.only.is_empty() {
            let mut names = self
                .only
                .iter()
                .take(4)
                .cloned()
                .collect::<Vec<_>>()
                .join(", ");
            if self.only.len() > 4 {
                names.push_str(", \u{2026}");
            }
            return format!("chosen enzymes ({names})");
        }
        let cuts = match self.max_cuts {
            1 => "unique",
            2 => "1-2 site",
            _ => "1-3 site",
        };
        let length = if self.min_length == 0 {
            String::new()
        } else {
            format!(" {}+", self.min_length)
        };
        match &self.pool {
            Some(set) => format!("{cuts}{length} cutters from {set}"),
            None => {
                let sold = if self.commercial_only { ", sold" } else { "" };
                format!("{cuts}{length} cutters{sold}")
            }
        }
    }
}

/// How the History tab draws the construction history.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum HistoryLayout {
    /// Indented lines.
    #[default]
    List,
    /// A tree of cards with maps.
    Diagram,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct UiState {
    pub view: View,
    /// Whether the Cloning experiment shows its construct as a sequence rather than a map.
    /// The layers (features, enzymes, primers, ORFs, ruler) are the same ones the document
    /// views use, so a toggle means the same thing wherever it is pressed.
    pub cloning_sequence: bool,
    /// Whether the sequence view shows the whole molecule on one line, scrolled sideways,
    /// instead of wrapping it into rows. Some work (reading along an alignment, following a
    /// primer past a mismatch) is easier to follow when nothing wraps.
    pub sequence_single_line: bool,
    /// Whether the sequence view shows its "Align with" panel: reads and other sequences
    /// aligned under the document, in a pane below the sequence.
    pub align_panel: bool,
    /// What each document is aligned with (the "Align with" lists).
    pub align_with: Vec<AlignWith>,
    pub theme: ThemeChoice,
    pub map_style: MapStyle,
    /// The folder shown in the folder panel.
    pub folder: Option<String>,
    /// Whether the folder panel is shown.
    pub sidebar: bool,
    /// Map layers: features, restriction sites, primer binding sites, positions in site
    /// labels, and the ruler (ticks).
    pub features: bool,
    pub enzymes: bool,
    pub primers: bool,
    pub positions: bool,
    pub ruler: bool,
    /// Open reading frames on the map.
    pub orfs: bool,
    /// The simulated gel: agarose, ladder and lanes.
    pub gel: GelSetup,
    /// Open experiments, each in a tab after the document tabs.
    pub experiments: Vec<Experiment>,
    /// The experiment tab shown instead of the active document, if any.
    pub active_experiment: Option<usize>,
    /// The History tab as a list or a diagram.
    pub history: HistoryLayout,
    /// What the ORF layer shows in the sequence view.
    pub translation: Translation,
    /// Which enzymes the map and the sequence view show.
    pub enzyme_set: EnzymeSet,
    /// The background of the map and sequence views (Preferences).
    pub map_background: CanvasColors,
    /// The Preferences window is open, and the category it shows.
    pub preferences: bool,
    pub preferences_page: crate::preferences::PrefsPage,
    /// Interface size (1 = normal), as set in Preferences or with Ctrl/Cmd + and -.
    pub ui_scale: f32,
    /// Size of the residues in the sequence view and the Align with pane, in points.
    #[serde(default = "sequence_font_default")]
    pub sequence_font: f32,
    /// Reopen the window at the size it had when the app last ran.
    pub restore_window: bool,
    /// Map zoom (1 = the whole molecule fits).
    pub map_zoom: f32,
    /// NCBI genetic code of the frame translations.
    pub genetic_code: u8,
    /// How primer tails show in the sequence view (Preferences > Primers).
    pub primer_tails: TailStyle,
    /// Chromatogram width: points per sample point of the signal.
    pub trace_zoom: f32,
    /// Chromatogram height: 1 = the highest peak of the read fills the view.
    pub trace_gain: f32,
    /// The last status message shown in the status bar.
    #[serde(skip_deserializing)]
    pub status: Option<String>,
}

impl Default for UiState {
    fn default() -> Self {
        UiState {
            view: View::default(),
            cloning_sequence: false,
            sequence_single_line: false,
            align_panel: false,
            align_with: Vec::new(),
            theme: ThemeChoice::default(),
            map_style: MapStyle::default(),
            folder: None,
            sidebar: true,
            features: true,
            enzymes: true,
            primers: true,
            positions: true,
            ruler: true,
            orfs: false,
            gel: GelSetup::default(),
            experiments: Vec::new(),
            active_experiment: None,
            history: HistoryLayout::default(),
            translation: Translation::default(),
            enzyme_set: EnzymeSet::default(),
            map_background: CanvasColors::default(),
            preferences: false,
            preferences_page: crate::preferences::PrefsPage::default(),
            ui_scale: 1.0,
            sequence_font: SEQUENCE_FONT_DEFAULT,
            restore_window: true,
            map_zoom: 1.0,
            genetic_code: 1,
            primer_tails: TailStyle::default(),
            trace_zoom: TRACE_ZOOM_DEFAULT,
            trace_gain: 1.0,
            status: None,
        }
    }
}

/// Most sequences one document can be aligned with: the most `check_reads` takes.
pub const MAX_ALIGN_WITH: usize = snappiesnap_engine::bio::reads::MAX_READS;

/// Where a sequence aligned under a document comes from.
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum AlignSource {
    /// A file (a trace, FASTA, GenBank, .dna).
    File(String),
    /// A document open in a tab, by id.
    Document(u64),
    /// A read stored with the document in its .dna file, by its index there.
    Matched(usize),
}

/// One sequence of an "Align with" list.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct AlignEntry {
    pub source: AlignSource,
    /// Ticked: aligned and shown in the pane.
    #[serde(default = "yes")]
    pub on: bool,
    /// Its trace is folded away, leaving one line of bases.
    #[serde(default)]
    pub collapsed: bool,
}

fn yes() -> bool {
    true
}

/// The sequences one document is aligned with, in the order they are shown.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct AlignWith {
    /// The document's id.
    pub document: u64,
    #[serde(default)]
    pub entries: Vec<AlignEntry>,
}

impl UiState {
    /// The "Align with" list of document `id`, made if it has none.
    pub fn align_list(&mut self, id: u64) -> Option<&mut Vec<AlignEntry>> {
        if !self.align_with.iter().any(|l| l.document == id) {
            self.align_with.push(AlignWith {
                document: id,
                entries: Vec::new(),
            });
        }
        self.align_with
            .iter_mut()
            .find(|l| l.document == id)
            .map(|l| &mut l.entries)
    }

    /// The "Align with" list of document `id`, if it has one.
    pub fn align_entries(&self, id: u64) -> &[AlignEntry] {
        self.align_with
            .iter()
            .find(|l| l.document == id)
            .map_or(&[], |l| l.entries.as_slice())
    }

    /// Add a sequence to document `id`'s list once, ticked; false when the list is full.
    pub fn align_add(&mut self, id: u64, source: AlignSource) -> bool {
        let Some(list) = self.align_list(id) else {
            return false;
        };
        if list.iter().any(|e| e.source == source) {
            return true;
        }
        if list.len() >= MAX_ALIGN_WITH {
            return false;
        }
        list.push(AlignEntry {
            source,
            on: true,
            collapsed: false,
        });
        true
    }
}

/// Chromatogram zoom: points per sample point (traces have about 12 samples per base).
pub const TRACE_ZOOM_DEFAULT: f32 = 1.5;
pub const TRACE_ZOOM_RANGE: std::ops::RangeInclusive<f32> = 0.1..=12.0;
pub const TRACE_GAIN_RANGE: std::ops::RangeInclusive<f32> = 0.25..=16.0;

/// A chromatogram zoom within range (the default for a value that isn't a number).
pub fn clamp_trace_zoom(z: f32) -> f32 {
    if z.is_finite() {
        z.clamp(*TRACE_ZOOM_RANGE.start(), *TRACE_ZOOM_RANGE.end())
    } else {
        TRACE_ZOOM_DEFAULT
    }
}

/// A chromatogram height within range (1 for a value that isn't a number).
pub fn clamp_trace_gain(g: f32) -> f32 {
    if g.is_finite() {
        g.clamp(*TRACE_GAIN_RANGE.start(), *TRACE_GAIN_RANGE.end())
    } else {
        1.0
    }
}

/// The sequence font size until the user picks another, in points.
pub const SEQUENCE_FONT_DEFAULT: f32 = 13.0;

fn sequence_font_default() -> f32 {
    SEQUENCE_FONT_DEFAULT
}

/// A sequence font size kept within Preferences' range, in whole points (the default for a
/// value that isn't a number).
pub fn clamp_sequence_font(z: f32) -> f32 {
    let range = crate::preferences::SEQUENCE_FONT_RANGE;
    if z.is_finite() {
        z.round().max(*range.start()).min(*range.end())
    } else {
        SEQUENCE_FONT_DEFAULT
    }
}

/// An interface size kept within Preferences' range (1 for a value that isn't a number).
pub fn clamp_ui_scale(z: f32) -> f32 {
    let range = crate::preferences::UI_SCALE_RANGE;
    if z.is_finite() {
        z.clamp(*range.start(), *range.end())
    } else {
        1.0
    }
}

/// A partial update for `ui.set`: every field is optional.
#[derive(Debug, Clone, Default, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct UiPatch {
    pub view: Option<View>,
    pub theme: Option<ThemeChoice>,
    pub map_style: Option<MapStyle>,
    pub folder: Option<String>,
    pub sidebar: Option<bool>,
    pub features: Option<bool>,
    pub enzymes: Option<bool>,
    pub primers: Option<bool>,
    pub positions: Option<bool>,
    pub ruler: Option<bool>,
    pub orfs: Option<bool>,
    pub gel: Option<GelSetup>,
    pub experiments: Option<Vec<Experiment>>,
    /// Show this experiment tab (index into `experiments`; negative: the active document).
    pub active_experiment: Option<i64>,
    pub history: Option<HistoryLayout>,
    pub translation: Option<Translation>,
    pub genetic_code: Option<u8>,
    pub enzyme_set: Option<EnzymeSet>,
    pub map_background: Option<CanvasColors>,
    pub preferences: Option<bool>,
    pub preferences_page: Option<crate::preferences::PrefsPage>,
    pub ui_scale: Option<f32>,
    pub sequence_font: Option<f32>,
    pub restore_window: Option<bool>,
    pub map_zoom: Option<f32>,
    pub trace_zoom: Option<f32>,
    pub trace_gain: Option<f32>,
    pub primer_tails: Option<TailStyle>,
    /// The sequence view on one line, scrolled sideways, instead of wrapped into rows.
    pub sequence_single_line: Option<bool>,
    /// Show or hide the sequence view's "Align with" panel.
    pub align_panel: Option<bool>,
    /// Replace the "Align with" lists (at most [`MAX_ALIGN_WITH`] sequences each).
    pub align_with: Option<Vec<AlignWith>>,
    /// Show this tab's document beside the active one on the map (tab index; negative: off).
    pub compare: Option<i64>,
}

impl UiState {
    pub fn apply(&mut self, patch: UiPatch) {
        if let Some(v) = patch.view {
            self.view = v;
        }
        if let Some(t) = patch.theme {
            self.theme = t;
        }
        if let Some(m) = patch.map_style {
            self.map_style = m;
        }
        if let Some(f) = patch.folder {
            self.folder = Some(f);
        }
        if let Some(b) = patch.sidebar {
            self.sidebar = b;
        }
        if let Some(b) = patch.sequence_single_line {
            self.sequence_single_line = b;
        }
        if let Some(b) = patch.align_panel {
            self.align_panel = b;
        }
        if let Some(lists) = patch.align_with {
            self.align_with = lists
                .into_iter()
                .take(MAX_EXPERIMENTS * 5)
                .map(|mut l| {
                    l.entries.truncate(MAX_ALIGN_WITH);
                    l
                })
                .collect();
        }
        if let Some(g) = patch.gel {
            self.gel = g;
        }
        if let Some(list) = patch.experiments {
            self.experiments = list
                .into_iter()
                .take(MAX_EXPERIMENTS)
                .map(Experiment::clamped)
                .collect();
            self.active_experiment = self
                .active_experiment
                .filter(|&i| i < self.experiments.len());
        }
        if let Some(i) = patch.active_experiment {
            self.active_experiment = usize::try_from(i)
                .ok()
                .filter(|&i| i < self.experiments.len());
        }
        if let Some(h) = patch.history {
            self.history = h;
        }
        if let Some(z) = patch.map_zoom
            && z.is_finite()
        {
            self.map_zoom = z.clamp(1.0, MAX_ZOOM);
        }
        if let Some(t) = patch.primer_tails {
            self.primer_tails = t;
        }
        if let Some(z) = patch.trace_zoom {
            self.trace_zoom = clamp_trace_zoom(z);
        }
        if let Some(g) = patch.trace_gain {
            self.trace_gain = clamp_trace_gain(g);
        }
        if let Some(e) = patch.enzyme_set {
            self.enzyme_set = e.clamped();
        }
        if let Some(c) = patch.map_background {
            self.map_background = c.clamped();
        }
        if let Some(p) = patch.preferences_page {
            self.preferences_page = p;
        }
        if let Some(z) = patch.ui_scale {
            self.ui_scale = clamp_ui_scale(z);
        }
        if let Some(z) = patch.sequence_font {
            self.sequence_font = clamp_sequence_font(z);
        }
        if let Some(tr) = patch.translation {
            self.translation = tr;
        }
        if let Some(c) = patch.genetic_code
            && snappiesnap_engine::bio::translate::code(c).is_some()
        {
            self.genetic_code = c;
        }
        for (value, field) in [
            (patch.features, &mut self.features),
            (patch.enzymes, &mut self.enzymes),
            (patch.primers, &mut self.primers),
            (patch.positions, &mut self.positions),
            (patch.ruler, &mut self.ruler),
            (patch.orfs, &mut self.orfs),
            (patch.preferences, &mut self.preferences),
            (patch.restore_window, &mut self.restore_window),
        ] {
            if let Some(b) = value {
                *field = b;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn experiment_lanes_keep_their_order_and_methods() {
        use snappiesnap_engine::gel::LaneSpec;
        let mut e = Experiment::new(ExperimentKind::Gel, "Run 1", [7, 8, 9]);
        // Simulate Gel opens on its first step, which asks what to simulate; the constructs
        // are already seeded behind it.
        assert_eq!(
            (e.step, e.name.as_str()),
            (ExperimentStep::Simulation, "Run 1")
        );
        e.step = ExperimentStep::Constructs;
        let docs = |e: &Experiment| e.lanes.iter().map(|l| l.document).collect::<Vec<_>>();
        e.lanes[2].method = Some(LaneSpec::Uncut);
        e.selected = Some(2);
        e.move_lane(2, true);
        assert_eq!(docs(&e), vec![7, 9, 8]);
        assert_eq!(
            (e.selected, &e.lanes[1].method),
            (Some(1), &Some(LaneSpec::Uncut))
        );
        // Moves past either end do nothing.
        e.move_lane(0, true);
        e.move_lane(2, false);
        e.move_lane(99, false);
        assert_eq!(docs(&e), vec![7, 9, 8]);
        e.add(7);
        assert_eq!(docs(&e), vec![7, 9, 8, 7]);
        e.remove(0);
        assert_eq!((docs(&e), e.selected), (vec![9, 8, 7], Some(0)));
        e.remove(0);
        assert_eq!(e.selected, None);
        e.remove(99);
        e.apply_to_all(&LaneSpec::Digest {
            enzymes: vec!["EcoRI".into()],
        });
        assert!(e.lanes.iter().all(|l| l.method.is_some()));
        for _ in 0..100 {
            e.add(1);
        }
        assert_eq!(e.lanes.len(), Experiment::MAX_LANES);
        let wild = Experiment {
            agarose: f64::NAN,
            selected: Some(500),
            ..e.clone()
        }
        .clamped();
        assert_eq!((wild.agarose, wild.selected), (1.0, None));
        // ui.set takes the experiments as JSON.
        let patch: UiPatch = serde_json::from_str(
            r#"{"experiments": [{"name": "  ", "step": "lanes", "lanes": [{"document": 3, "method": {"kind": "uncut"}}], "selected": 0}], "active_experiment": 0}"#,
        )
        .unwrap();
        let mut s = UiState::default();
        s.apply(patch);
        let e = s.experiment().unwrap();
        assert_eq!((e.step, e.lanes.len()), (ExperimentStep::Lanes, 1));
        assert_eq!(e.name, "Simulate Gel");
        s.apply(serde_json::from_str(r#"{"active_experiment": -1}"#).unwrap());
        assert_eq!(s.active_experiment, None);
        s.apply(serde_json::from_str(r#"{"active_experiment": 9}"#).unwrap());
        assert_eq!(s.active_experiment, None);
    }

    #[test]
    fn map_background_per_theme() {
        let mut c = CanvasColors::default();
        assert_eq!(c.for_theme(true), None);
        c.set(true, Some(egui::Color32::from_rgb(1, 2, 3)));
        assert_eq!(c.dark.as_deref(), Some("#010203"));
        assert_eq!(c.for_theme(true), Some(egui::Color32::from_rgb(1, 2, 3)));
        assert_eq!(c.for_theme(false), None);
        c.set(true, None);
        assert_eq!(c, CanvasColors::default());
        let mut s = UiState::default();
        s.apply(
            serde_json::from_str(
                r##"{"map_background": {"light": "#FFFFEE", "dark": "not a colour"}}"##,
            )
            .unwrap(),
        );
        assert_eq!(s.map_background.light.as_deref(), Some("#ffffee"));
        assert_eq!(s.map_background.dark, None);
    }

    #[test]
    fn experiment_tabs_open_and_close() {
        let mut s = UiState::default();
        s.new_experiment(ExperimentKind::Gel, [1, 2]);
        s.new_experiment(ExperimentKind::Gel, [2]);
        s.new_experiment(ExperimentKind::Gel, []);
        let names = |s: &UiState| {
            s.experiments
                .iter()
                .map(|e| e.name.clone())
                .collect::<Vec<_>>()
        };
        assert_eq!(
            names(&s),
            vec!["Gel experiment", "Gel experiment 2", "Gel experiment 3"]
        );
        assert_eq!(s.active_experiment, Some(2));
        // Closing the shown tab shows its neighbour; closing an earlier one keeps the shown one.
        s.close_experiment(2);
        assert_eq!(s.active_experiment, Some(1));
        s.close_experiment(0);
        assert_eq!(
            (s.active_experiment, names(&s)),
            (Some(0), vec!["Gel experiment 2".to_string()])
        );
        // A freed name is used again.
        s.new_experiment(ExperimentKind::Gel, [1]);
        assert_eq!(s.experiments[1].name, "Gel experiment");
        s.close_experiment(99);
        s.close_experiment(0);
        s.close_experiment(0);
        assert_eq!((s.experiments.len(), s.active_experiment), (0, None));
        for _ in 0..MAX_EXPERIMENTS + 3 {
            s.new_experiment(ExperimentKind::Gel, [1]);
        }
        assert_eq!(s.experiments.len(), MAX_EXPERIMENTS);
        assert!(s.status.is_some());
    }

    #[test]
    fn patches_and_serde() {
        let mut s = UiState::default();
        let patch: UiPatch =
            serde_json::from_str(r#"{"view": "sequence", "theme": "dark"}"#).unwrap();
        s.apply(patch);
        assert_eq!(s.view, View::Sequence);
        assert_eq!(s.theme, ThemeChoice::Dark);
        assert!(serde_json::from_str::<UiPatch>(r#"{"view": "nope"}"#).is_err());
        assert!(serde_json::from_str::<UiPatch>(r#"{"zoom": 2}"#).is_err());
        s.apply(serde_json::from_str(r#"{"primers": false, "ruler": false}"#).unwrap());
        assert!(!s.primers && !s.ruler && s.features && s.enzymes && s.positions);
        let json = serde_json::to_value(&s).unwrap();
        assert_eq!(json["view"], "sequence");
    }

    #[test]
    fn opened_experiments_track_unsaved_changes() {
        let opened = serde_json::json!({
            "name": "Digest check", "agarose": 2.0, "ladder": "100 bp",
            "lanes": [{"document": 3, "method": {"kind": "uncut"}}, {"document": 4, "method": null}],
        });
        let mut e = Experiment::from_opened("/x/run.sslab", &opened).unwrap();
        assert_eq!(e.step, ExperimentStep::Lanes);
        assert_eq!(e.lanes.len(), 2);
        assert!(!e.unsaved());
        assert_eq!(e.save_params("/y.sslab")["lanes"][0]["document"], 3);
        e.agarose = 1.0;
        assert!(e.unsaved());
        // Starting over keeps the file it belongs to.
        e.start([7]);
        assert_eq!(e.path.as_deref(), Some("/x/run.sslab"));
        assert!(Experiment::from_opened("p", &serde_json::json!({"name": 1})).is_none());
        // A diagnostic digest saves its design and reopens as one.
        let mut d = Experiment::new(ExperimentKind::Diagnostic, "Check", [3]);
        d.lanes[0].method = Some(snappiesnap_engine::gel::LaneSpec::Uncut);
        let params = d.save_params("/c.sslab");
        assert_eq!(params["design"], "diagnostic_digest");
        let mut reopened = params.clone();
        reopened["lanes"] = serde_json::json!([{"document": 3, "method": {"kind": "uncut"}}]);
        let back = Experiment::from_opened("/c.sslab", &reopened).unwrap();
        assert_eq!(back.kind, ExperimentKind::Diagnostic);
        assert!(e.save_params("/y.sslab").get("design").is_none());
        assert!(Experiment::default().unsaved());
    }
    /// Simulate Gel is one menu item covering both gel kinds: its first step chooses which,
    /// and a tab name still left as the default follows the choice.
    #[test]
    fn simulate_gel_covers_both_gel_kinds() {
        assert_eq!(
            ExperimentKind::MENU,
            [
                ExperimentKind::Gel,
                ExperimentKind::Cloning,
                ExperimentKind::Reads
            ]
        );
        assert!(ExperimentKind::Gel.chooses_simulation());
        assert!(ExperimentKind::Diagnostic.chooses_simulation());
        assert!(!ExperimentKind::Reads.chooses_simulation());
        // Check Reads keeps its own item and its own first step.
        let reads = Experiment::new(ExperimentKind::Reads, "Read check", [1]);
        assert_eq!(reads.step, ExperimentStep::Constructs);

        let mut e = Experiment::new(ExperimentKind::Gel, "Gel experiment", [1, 2]);
        assert_eq!(e.step, ExperimentStep::Simulation);
        e.set_simulation(ExperimentKind::Diagnostic);
        assert_eq!(e.kind, ExperimentKind::Diagnostic);
        // The default name followed the choice, and the constructs were kept.
        assert_eq!(e.name, "Diagnostic digest");
        assert_eq!(e.lanes.len(), 2);
        // A numbered default follows it too.
        let mut numbered = Experiment::new(ExperimentKind::Gel, "Gel experiment 3", [1]);
        numbered.set_simulation(ExperimentKind::Diagnostic);
        assert_eq!(numbered.name, "Diagnostic digest 3");
        // A name of the user's own is left alone.
        let mut mine = Experiment::new(ExperimentKind::Gel, "Monday clones", [1]);
        mine.set_simulation(ExperimentKind::Diagnostic);
        assert_eq!(mine.name, "Monday clones");
    }
}
