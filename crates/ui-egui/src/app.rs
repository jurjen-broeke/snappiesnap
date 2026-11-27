//! The application shell: toolbar, view tabs, status bar and the central view. Thin by design:
//! every action is an engine command run through [`App::run`], exactly like the control channel
//! and MCP.

use eframe::egui::{
    self, Align, CornerRadius, Key, KeyboardShortcut, Layout, Modifiers, Rect, RichText, Sense,
    Stroke, Theme, Vec2, vec2,
};
use serde_json::{Value, json};
#[cfg(not(target_arch = "wasm32"))]
use snappiesnap_engine::{CmdError, api};
use snappiesnap_engine::{ErrorCode, Level, Session};

use crate::icons::{self, Icon};
use crate::map::{self, MapAction};
use crate::seqview::{self, SeqAction, SeqViewState};
use crate::sidebar::{self, Sidebar, SidebarAction};
#[cfg(not(target_arch = "wasm32"))]
use crate::state::UiPatch;
use crate::state::{MapStyle, Translation, UiState, View};
use crate::theme::{ThemeChoice, Tokens};

/// What the app should do at startup.
#[derive(Debug, Clone, Default)]
pub struct LaunchOptions {
    /// Start the JSON control channel on this port (0 = any free port).
    pub control_port: Option<u16>,
    /// Record every successful edit to this file as a replayable command script.
    pub record: Option<String>,
    /// Load and save settings (theme, view, folder panel) between sessions.
    pub remember_settings: bool,
}

/// A command that needs confirmation because it would drop unsaved changes.
struct Confirm {
    id: String,
    params: Value,
    question: String,
}

pub struct App {
    session: Session,
    state: UiState,
    tokens: Tokens,
    applied: Option<Tokens>,
    #[cfg(not(target_arch = "wasm32"))]
    control: Option<crate::control::ControlServer>,
    confirm: Option<Confirm>,
    title: String,
    seq_view: SeqViewState,
    trace_view: crate::trace_view::TraceViewState,
    sidebar: Sidebar,
    /// Where settings are remembered (desktop app only).
    #[cfg(not(target_arch = "wasm32"))]
    settings_path: Option<std::path::PathBuf>,
    #[cfg(not(target_arch = "wasm32"))]
    saved_settings: Option<crate::settings::Settings>,
    /// Where `--record` writes the session's commands.
    #[cfg(not(target_arch = "wasm32"))]
    recorder: Option<std::io::BufWriter<std::fs::File>>,
    /// The feature being renamed in the features table, and the text so far.
    /// Restriction sites for the map, per document revision.
    site_cache: crate::enzymes::SiteCache,
    enzyme_view: crate::enzyme_view::EnzymeViewState,
    primer_design: crate::primer_design::DesignState,
    gel_inputs: crate::gel_view::GelInputs,
    align: crate::align_view::AlignState,
    /// The sequence view's "Align with" pane: the last check and what drawing it needs.
    align_with: crate::align_with::AlignWithView,
    find: crate::find_bar::FindState,
    /// The Go to bar of the map and sequence views.
    goto: crate::goto_bar::GoToState,
    /// The status bar's composition summary, and what it was computed for (document,
    /// revision, selection).
    info: Option<(InfoKey, Vec<String>)>,
    features: crate::features_view::FeaturesState,
    /// Map pan (with `state.map_zoom`).
    map_pan: Vec2,
    /// The document shown beside the active one on the map (tab id), and its marks.
    compare: Option<u64>,
    compare_cache: crate::enzymes::SiteCache,
    notes_edit: crate::notes_view::NotesEdit,
    /// The primer being edited in the Primers tab: id, name, sequence, description.
    primer_edit: Option<(u64, String, String, String)>,
    /// How each open document was being viewed, by tab id; the active one's lives in the
    /// fields above.
    tab_states: std::collections::HashMap<u64, TabUi>,
    current_tab: Option<u64>,
    gel_cache: crate::gel_view::GelCache,
    /// The shown experiment tab's method editor and cache, and which tab it was filled for
    /// (experiments span documents, so this
    /// is not kept per tab).
    experiment_inputs: (Option<usize>, crate::experiment_view::ExperimentInputs),
    /// The interface size last given to egui, to tell a keyboard zoom from a Preferences change.
    applied_scale: Option<f32>,
    /// The window's size to remember, and the last size seen with when it was first seen (it
    /// is remembered once it has stayed the same for a moment, not on every resize step).
    #[cfg(not(target_arch = "wasm32"))]
    window: Option<crate::settings::WindowSize>,
    #[cfg(not(target_arch = "wasm32"))]
    window_seen: Option<(crate::settings::WindowSize, f64)>,
    /// A gel the user asked to save as an image, while asking whether to invert it.
    gel_export: Option<GelSource>,
    /// The Primers tab's lab-oligo section (the library is the lab's, not one document's).
    oligo_view: crate::oligo_view::OligoViewState,
    history_cache: crate::history_view::HistoryCache,
    /// An update check the user started, and the last one's answer.
    update_check: Option<std::sync::mpsc::Receiver<crate::updates::CheckResult>>,
    update_answer: Option<String>,
    /// A REBASE download running in the background.
    #[cfg(not(target_arch = "wasm32"))]
    download: Option<std::sync::mpsc::Receiver<Result<String, String>>>,
    /// Web build: files the browser read for us, waiting to be opened.
    #[cfg_attr(not(target_arch = "wasm32"), allow(dead_code))]
    web_inbox: crate::web_files::Inbox,
    /// Enzyme set files the browser has read, waiting to be imported.
    web_sets: crate::web_files::Inbox,
    /// Sequence files the browser read for the Align view or the "Align with" list, and which
    /// of them they are for.
    web_align: crate::web_files::Inbox,
    web_align_for: AlignFor,
    /// Read files the browser read for Check Reads, and the experiment they belong to.
    web_reads: crate::web_files::Inbox,
    web_reads_experiment: usize,
    /// Web build: the Download prompt that stands in for a Save dialog.
    #[cfg(target_arch = "wasm32")]
    web_save: Option<crate::web_files::SavePrompt>,
    /// Web build: the context, for the file picker to wake the app when files arrive.
    #[cfg(target_arch = "wasm32")]
    ctx: egui::Context,
    /// Web build: whether anything is unsaved, for the browser to ask before the page closes.
    #[cfg(target_arch = "wasm32")]
    web_unsaved: std::sync::Arc<std::sync::atomic::AtomicBool>,
}

impl App {
    pub fn new(ctx: &egui::Context, session: Session, options: &LaunchOptions) -> Self {
        #[cfg_attr(target_arch = "wasm32", allow(unused_mut))]
        let mut app = App {
            session,
            state: UiState::default(),
            tokens: Tokens::dark(),
            applied: None,
            #[cfg(not(target_arch = "wasm32"))]
            control: None,
            confirm: None,
            title: String::new(),
            seq_view: SeqViewState::default(),
            trace_view: Default::default(),
            sidebar: Sidebar::default(),
            #[cfg(not(target_arch = "wasm32"))]
            settings_path: None,
            #[cfg(not(target_arch = "wasm32"))]
            saved_settings: None,
            #[cfg(not(target_arch = "wasm32"))]
            recorder: None,
            site_cache: Default::default(),
            enzyme_view: Default::default(),
            primer_design: Default::default(),
            gel_inputs: Default::default(),
            align: Default::default(),
            align_with: Default::default(),
            find: Default::default(),
            goto: Default::default(),
            info: None,
            features: Default::default(),
            map_pan: Vec2::ZERO,
            compare: None,
            compare_cache: Default::default(),
            notes_edit: Default::default(),
            primer_edit: None,
            tab_states: Default::default(),
            current_tab: None,
            gel_cache: Default::default(),
            experiment_inputs: Default::default(),
            applied_scale: None,
            #[cfg(not(target_arch = "wasm32"))]
            window: None,
            #[cfg(not(target_arch = "wasm32"))]
            window_seen: None,
            gel_export: None,
            oligo_view: Default::default(),
            history_cache: Default::default(),
            update_check: None,
            update_answer: None,
            #[cfg(not(target_arch = "wasm32"))]
            download: None,
            web_inbox: Default::default(),
            web_sets: Default::default(),
            web_align: Default::default(),
            web_align_for: AlignFor::Pairwise,
            web_reads: Default::default(),
            web_reads_experiment: 0,
            #[cfg(target_arch = "wasm32")]
            web_save: None,
            #[cfg(target_arch = "wasm32")]
            ctx: ctx.clone(),
            #[cfg(target_arch = "wasm32")]
            web_unsaved: Default::default(),
        };
        #[cfg(not(target_arch = "wasm32"))]
        if options.remember_settings
            && let Some(path) = crate::settings::default_path()
        {
            let (s, problem) = crate::settings::load(&path);
            app.state.theme = s.theme;
            app.state.view = s.view;
            app.state.map_style = s.map_style;
            app.state.folder = s.folder.clone();
            app.state.sidebar = s.sidebar;
            app.state.features = s.features;
            app.state.enzymes = s.enzymes;
            app.state.primers = s.primers;
            app.state.positions = s.positions;
            app.state.ruler = s.ruler;
            app.state.orfs = s.orfs;
            app.state.enzyme_set = s.enzyme_set.clone().clamped();
            app.enzyme_view.from = s.enzyme_list.clone();
            app.state.map_background = s.map_background.clone();
            app.state.preferences_page = s.preferences_page;
            app.state.ui_scale = crate::state::clamp_ui_scale(s.ui_scale);
            app.state.sequence_font = crate::state::clamp_sequence_font(s.sequence_font);
            app.state.restore_window = s.restore_window;
            app.session.set_primer_binding(s.primer_binding);
            app.state.translation = s.translation;
            app.state.primer_tails = s.primer_tails;
            if snappiesnap_engine::bio::translate::code(s.genetic_code).is_some() {
                app.state.genetic_code = s.genetic_code;
            }
            app.window = s.window;
            app.sidebar.expanded = s.expanded.iter().cloned().collect();
            if let Some(p) = problem {
                app.status_error(p);
            }
            app.saved_settings = Some(s);
            app.settings_path = Some(path);
        }
        app.follow_document_folder();
        #[cfg(not(target_arch = "wasm32"))]
        if let Some(path) = &options.record {
            match std::fs::File::create(path) {
                Ok(file) => {
                    app.recorder = Some(std::io::BufWriter::new(file));
                    app.session.set_journaling(true);
                }
                Err(e) => app.status_error(format!("cannot record to {path}: {e}")),
            }
        }
        #[cfg(not(target_arch = "wasm32"))]
        if let Some(port) = options.control_port {
            match crate::control::ControlServer::start(port, ctx.clone()) {
                Ok(server) => {
                    // Agents read the port from stderr when they asked for port 0.
                    eprintln!("control channel listening on 127.0.0.1:{}", server.port);
                    app.control = Some(server);
                }
                Err(e) => app.status_error(format!("control channel not started: {e}")),
            }
        }
        #[cfg(target_arch = "wasm32")]
        {
            let _ = options;
            app.state.sidebar = false;
            if let Err(e) = crate::web::warn_before_leaving(app.web_unsaved.clone()) {
                app.status_error(format!(
                    "closing this page won't ask about unsaved work: {e}"
                ));
            }
        }
        app
    }

    /// An app without a window, for offscreen rendering (CLI export, tests).
    pub fn headless(session: Session) -> Self {
        App::new(
            &egui::Context::default(),
            session,
            &LaunchOptions::default(),
        )
    }

    pub fn session(&self) -> &Session {
        &self.session
    }

    /// The current view, theme and map style.
    pub fn ui_settings(&self) -> &UiState {
        &self.state
    }

    /// Change the view, theme or map style (as `ui.set` does).
    pub fn set_ui(&mut self, patch: crate::state::UiPatch) {
        self.state.apply(patch);
    }

    /// Show the open document's folder if no folder is chosen yet.
    /// Show the open document's folder when it isn't inside the panel's folder.
    fn follow_document_folder(&mut self) {
        let Some(path) = self.session.document().and_then(|d| d.path.clone()) else {
            return;
        };
        let inside = self
            .state
            .folder
            .as_deref()
            .is_some_and(|f| sidebar::is_within(&path, f));
        if !inside && let Some(folder) = sidebar::parent_folder(&path) {
            self.state.folder = Some(folder);
        }
    }

    fn sidebar_panel(&mut self, ui: &mut egui::Ui) {
        let t = self.tokens;
        let platform = std::sync::Arc::clone(&self.platform_handle());
        let current = self.session.document().and_then(|d| d.path.clone());
        let action = sidebar::show(
            ui,
            &t,
            platform.as_ref(),
            self.state.folder.as_deref(),
            current.as_deref(),
            &mut self.sidebar,
        );
        match action {
            Some(SidebarAction::Open(path)) => {
                if sidebar::is_experiment_file(&path) {
                    self.open_experiment(&path);
                } else if current.as_deref() != Some(path.as_str()) {
                    self.run("open", json!({ "path": path }));
                }
            }
            Some(SidebarAction::ChooseFolder) => {
                #[cfg(not(target_arch = "wasm32"))]
                if let Some(dir) = rfd::FileDialog::new().pick_folder() {
                    self.state.folder = Some(dir.to_string_lossy().into_owned());
                    self.sidebar.refresh();
                }
            }
            Some(SidebarAction::Rename { from, to }) => {
                self.file_command("file_rename", json!({ "from": from, "to": to }));
            }
            Some(SidebarAction::Copy { from, to }) => {
                self.file_command("file_copy", json!({ "from": from, "to": to }));
            }
            Some(SidebarAction::Move { from, into }) => {
                let to = sidebar::join(&into, sidebar::base_name(&from));
                if self
                    .file_command("file_rename", json!({ "from": from, "to": to }))
                    .is_some()
                {
                    self.sidebar.expanded.insert(into);
                }
            }
            Some(SidebarAction::NewFolder(path)) => {
                self.file_command("folder_create", json!({ "path": path }));
            }
            Some(SidebarAction::Reveal(path)) => self.reveal_path(&path),
            Some(SidebarAction::CopyPath(path)) => ui.ctx().copy_text(path),
            Some(SidebarAction::Hide) => self.state.sidebar = false,
            None => {}
        }
    }

    /// Show a path in the platform's file manager (desktop only).
    fn reveal_path(&mut self, path: &str) {
        #[cfg(not(target_arch = "wasm32"))]
        if let Err(e) = reveal(path) {
            self.status_error(format!("cannot show {path}: {e}"));
        }
        #[cfg(target_arch = "wasm32")]
        let _ = path;
    }

    /// Run a file-management command and re-read the folder listings.
    fn file_command(&mut self, id: &str, params: Value) -> Option<Value> {
        let result = self.run(id, params);
        self.sidebar.refresh();
        result
    }

    /// Open files dropped onto the window from the OS (a dropped folder is shown in the panel).
    fn handle_dropped_files(&mut self, ctx: &egui::Context) {
        // In the browser a dropped file is read in the background; it is opened from the
        // inbox once it arrives (`take_web_files`).
        #[cfg(target_arch = "wasm32")]
        {
            let files = ctx.input(|i| i.raw.dropped_files.clone());
            let names: Vec<String> = files
                .iter()
                .map(|f| f.path().to_string_lossy().into_owned())
                .collect();
            if !files.is_empty() && self.drops_align_with(ctx, &names) {
                self.web_align_for = AlignFor::AlignWith;
                crate::web::read_dropped(&self.web_align, files, ctx);
            } else if !files.is_empty() {
                crate::web::read_dropped(&self.web_inbox, files, ctx);
            }
        }
        #[cfg(not(target_arch = "wasm32"))]
        self.handle_dropped_paths(ctx);
    }

    /// Open files dropped from the OS, by path (desktop).
    #[cfg(not(target_arch = "wasm32"))]
    fn handle_dropped_paths(&mut self, ctx: &egui::Context) {
        let dropped: Vec<String> = ctx.input(|i| {
            i.raw
                .dropped_files
                .iter()
                .map(|f| f.path().to_string_lossy().into_owned())
                .collect()
        });
        if dropped.is_empty() {
            return;
        }
        if self.drops_align_with(ctx, &dropped) {
            self.align_with_files(dropped);
            return;
        }
        let platform = self.session.platform_arc();
        if let Some(file) = dropped.iter().find(|p| sidebar::is_experiment_file(p)) {
            let file = file.clone();
            self.open_experiment(&file);
        } else if let Some(file) = dropped.iter().find(|p| sidebar::is_sequence_file(p)) {
            self.run("open", json!({ "path": file }));
        } else if let Some(dir) = dropped.iter().find(|p| platform.list_dir(p).is_ok()) {
            self.state.folder = Some(dir.clone());
            self.state.sidebar = true;
            self.sidebar.refresh();
        } else {
            self.status_error(
                "drop .dna, GenBank, EMBL, GFF3, SBOL or FASTA files, or a folder".to_string(),
            );
        }
    }

    /// Settings to remember, from the current state.
    #[cfg_attr(target_arch = "wasm32", allow(dead_code))]
    fn settings(&self) -> crate::settings::Settings {
        crate::settings::Settings {
            theme: self.state.theme,
            view: self.state.view,
            map_style: self.state.map_style,
            folder: self.state.folder.clone(),
            sidebar: self.state.sidebar,
            expanded: self.sidebar.expanded.iter().cloned().collect(),
            features: self.state.features,
            enzymes: self.state.enzymes,
            primers: self.state.primers,
            positions: self.state.positions,
            ruler: self.state.ruler,
            orfs: self.state.orfs,
            enzyme_set: self.state.enzyme_set.clone(),
            enzyme_list: self.enzyme_view.from.clone(),
            map_background: self.state.map_background.clone(),
            preferences_page: self.state.preferences_page,
            ui_scale: self.state.ui_scale,
            sequence_font: self.state.sequence_font,
            restore_window: self.state.restore_window,
            translation: self.state.translation,
            genetic_code: self.state.genetic_code,
            primer_tails: self.state.primer_tails,
            primer_binding: *self.session.primer_binding(),
            #[cfg(not(target_arch = "wasm32"))]
            window: self.window,
            #[cfg(target_arch = "wasm32")]
            window: None,
        }
    }

    /// Save settings when they changed (native only, when enabled).
    #[cfg(not(target_arch = "wasm32"))]
    fn save_settings(&mut self) {
        let Some(path) = self.settings_path.clone() else {
            return;
        };
        let current = self.settings();
        if self.saved_settings.as_ref() == Some(&current) {
            return;
        }
        if let Err(e) = crate::settings::save(&path, &current) {
            self.settings_path = None;
            self.status_error(format!("settings not saved: {e}"));
        }
        self.saved_settings = Some(current);
    }

    fn platform_handle(&self) -> std::sync::Arc<dyn snappiesnap_engine::Platform> {
        self.session.platform_arc()
    }

    /// Append journaled commands to the `--record` file.
    #[cfg(not(target_arch = "wasm32"))]
    fn write_recording(&mut self) {
        use std::io::Write;
        let entries = self.session.take_journal();
        let Some(out) = self.recorder.as_mut() else {
            return;
        };
        if entries.is_empty() {
            return;
        }
        let mut result = Ok(());
        for e in entries {
            result = result.and_then(|()| {
                writeln!(
                    out,
                    "{}",
                    json!({ "command": e.command, "params": e.params })
                )
            });
        }
        if let Err(e) = result.and_then(|()| out.flush()) {
            self.recorder = None;
            self.session.set_journaling(false);
            self.status_error(format!("recording stopped: {e}"));
        }
    }

    /// Run a command; failures go to the status bar and the log, never a panic.
    /// The tokens of the map and sequence views: the theme's, with the background chosen in
    /// Preferences for this theme.
    fn canvas_tokens(&self) -> Tokens {
        let mut t = self.tokens;
        if let Some(c) = self.state.map_background.for_theme(t.dark) {
            t.canvas = c;
        }
        t
    }

    /// The Preferences window, and what it asks the app to do.
    fn preferences(&mut self, ctx: &egui::Context) {
        use crate::preferences::PrefsAction;
        #[cfg(not(target_arch = "wasm32"))]
        let downloading = self.download.is_some();
        #[cfg(target_arch = "wasm32")]
        let downloading = false;
        let mut binding = *self.session.primer_binding();
        let extra = crate::preferences::Extra {
            primer_binding: &mut binding,
            enzymes: self.session.enzymes(),
            enzyme_list: &mut self.enzyme_view.from,
            downloading,
            checking: self.update_check.is_some(),
            update: self.update_answer.as_deref(),
        };
        let action = crate::preferences::window(ctx, &self.tokens, &mut self.state, extra);
        if binding != *self.session.primer_binding() {
            self.run(
                "primer_binding",
                serde_json::to_value(binding).unwrap_or_default(),
            );
        }
        match action {
            None => {}
            Some(PrefsAction::CheckForUpdates) => {
                self.update_answer = None;
                self.update_check = Some(crate::updates::start(self.platform_handle()));
            }
            Some(PrefsAction::UseBuiltinEnzymes) => {
                self.run("enzymes_reset", json!({}));
            }
            Some(PrefsAction::RemoveEnzymeSet(name)) => {
                if self
                    .run("enzyme_set_delete", json!({ "name": name }))
                    .is_some()
                {
                    let same = |n: &str| n.eq_ignore_ascii_case(&name);
                    if self.state.enzyme_set.pool.as_deref().is_some_and(same) {
                        self.state.enzyme_set.pool = None;
                    }
                    if matches!(&self.enzyme_view.from, Some(crate::enzyme_view::ListFrom::Set(n)) if same(n))
                    {
                        self.enzyme_view.from = None;
                    }
                    self.state.status = Some(format!("Enzyme set \u{201c}{name}\u{201d} removed"));
                }
            }
            #[cfg(not(target_arch = "wasm32"))]
            Some(PrefsAction::ChooseFolder) => {
                if let Some(dir) = rfd::FileDialog::new().pick_folder() {
                    self.state.folder = Some(dir.to_string_lossy().into_owned());
                    self.sidebar.refresh();
                }
            }
            #[cfg(not(target_arch = "wasm32"))]
            Some(PrefsAction::ImportEnzymeSet) => self.import_enzyme_set(),
            #[cfg(not(target_arch = "wasm32"))]
            Some(PrefsAction::DownloadRebase) if self.download.is_none() => {
                self.download = Some(crate::enzyme_view::start_download(&self.session));
            }
            #[cfg(not(target_arch = "wasm32"))]
            Some(PrefsAction::DownloadRebase) => {}
            #[cfg(target_arch = "wasm32")]
            Some(
                PrefsAction::ChooseFolder
                | PrefsAction::ImportEnzymeSet
                | PrefsAction::DownloadRebase,
            ) => {}
        }
    }

    /// Take the answer of a running update check, once it is there.
    fn poll_update_check(&mut self, ctx: &egui::Context) {
        let Some(rx) = &self.update_check else {
            return;
        };
        match rx.try_recv() {
            Ok(result) => {
                self.update_answer = Some(result.message());
                self.update_check = None;
            }
            Err(std::sync::mpsc::TryRecvError::Empty) => {
                ctx.request_repaint_after(std::time::Duration::from_millis(200));
            }
            Err(std::sync::mpsc::TryRecvError::Disconnected) => {
                self.update_answer = Some("Update check failed.".into());
                self.update_check = None;
            }
        }
    }

    /// Follow the interface size: a change in Preferences (or `ui.set`) goes to egui, and a
    /// zoom with the keyboard (Ctrl/Cmd + and -) comes back into Preferences.
    fn sync_scale(&mut self, ctx: &egui::Context) {
        let zoom = ctx.zoom_factor();
        if self.applied_scale == Some(self.state.ui_scale) && zoom != self.state.ui_scale {
            self.state.ui_scale = crate::state::clamp_ui_scale(zoom);
        }
        if zoom != self.state.ui_scale {
            ctx.set_zoom_factor(self.state.ui_scale);
        }
        self.applied_scale = Some(self.state.ui_scale);
    }

    /// Note the window's size once it has settled, to reopen at it next time. Sizes are in
    /// points before the interface size; while maximized the size from before is kept.
    #[cfg(not(target_arch = "wasm32"))]
    fn track_window(&mut self, ctx: &egui::Context) {
        let (rect, maximized, hidden, now) = ctx.input(|i| {
            let v = i.viewport();
            (
                v.inner_rect,
                v.maximized.unwrap_or(false),
                v.minimized == Some(true) || v.fullscreen == Some(true),
                i.time,
            )
        });
        let Some(rect) = rect.filter(|_| !hidden) else {
            return;
        };
        let zoom = ctx.zoom_factor();
        let seen = match (maximized, self.window) {
            (true, Some(w)) => crate::settings::WindowSize {
                maximized: true,
                ..w
            },
            _ => crate::settings::WindowSize {
                width: (rect.width() * zoom).round(),
                height: (rect.height() * zoom).round(),
                maximized,
            },
        };
        match self.window_seen {
            Some((last, since)) if last == seen => {
                if now - since >= 1.0 && self.window != Some(seen) && seen.usable().is_some() {
                    self.window = Some(seen);
                }
            }
            _ => self.window_seen = Some((seen, now)),
        }
    }

    /// A gel as an image: the open document's digest or PCR gel, or an experiment tab's gel.
    /// Returns the image and a file name to suggest.
    #[cfg_attr(target_arch = "wasm32", allow(dead_code))]
    fn gel_image(
        &mut self,
        source: GelSource,
        invert: bool,
    ) -> Result<(egui::ColorImage, String), String> {
        let (title, name, lanes) = match source {
            GelSource::Gel => {
                let doc = self
                    .session
                    .document()
                    .ok_or("open a document to save its gel")?;
                let name = doc.molecule.name().to_string();
                let kind = match self.state.gel.mode {
                    crate::state::GelMode::Digest => "restriction digest",
                    crate::state::GelMode::Pcr => "PCR",
                };
                let setup = self.state.gel.current().clone();
                let lanes = self.gel_cache.lanes(&self.session, &setup).clone();
                (
                    format!("{name}: {kind} gel, {:.1} % agarose", setup.agarose),
                    format!("{name} {kind} gel"),
                    lanes,
                )
            }
            GelSource::Experiment(i) => {
                let exp = self
                    .state
                    .experiments
                    .get(i)
                    .ok_or("that experiment is no longer open")?;
                let lanes = snappiesnap_engine::gel::experiment(
                    &self.session,
                    &exp.ladder,
                    exp.agarose,
                    &exp.lanes,
                );
                (
                    format!("{}, {:.1} % agarose", exp.name, exp.agarose),
                    exp.name.clone(),
                    lanes,
                )
            }
        };
        let image = crate::gel_export::gel_image(&title, &lanes, invert)?;
        // A file name without folder separators or other characters file systems refuse.
        let file: String = name
            .chars()
            .map(|c| if "/\\:*?\"<>|".contains(c) { '-' } else { c })
            .collect();
        Ok((image, file))
    }

    /// Ask whether to invert the gel, then where to save it.
    #[cfg(not(target_arch = "wasm32"))]
    fn gel_export_dialog(&mut self, ctx: &egui::Context) {
        let Some(source) = self.gel_export else {
            return;
        };
        let t = self.tokens;
        let mut answer: Option<Option<bool>> = None;
        egui::Window::new("Save gel image")
            .collapsible(false)
            .resizable(false)
            .anchor(egui::Align2::CENTER_CENTER, vec2(0.0, 0.0))
            .show(ctx, |ui| {
                ui.label("Invert the image before saving?");
                ui.label(
                    RichText::new(
                        "As shown: light bands on a dark gel, like a photo on a transilluminator. \
                         Inverted: dark bands on white, for printing and lab notebooks.",
                    )
                    .small()
                    .color(t.text_muted),
                );
                ui.add_space(t.spacing);
                ui.horizontal(|ui| {
                    if ui.button("As shown").clicked() {
                        answer = Some(Some(false));
                    }
                    if ui.button("Inverted").clicked() {
                        answer = Some(Some(true));
                    }
                    if ui.button("Cancel").clicked() {
                        answer = Some(None);
                    }
                });
            });
        let Some(answer) = answer else {
            return;
        };
        self.gel_export = None;
        let Some(invert) = answer else {
            return;
        };
        let (image, name) = match self.gel_image(source, invert) {
            Ok(x) => x,
            Err(e) => {
                self.status_error(format!("the gel image couldn't be made: {e}"));
                return;
            }
        };
        let suffix = if invert { " (inverted)" } else { "" };
        let mut dialog = rfd::FileDialog::new()
            .set_title("Save Gel Image")
            .add_filter("PNG image", &["png"])
            .add_filter("JPEG image", &["jpg", "jpeg"])
            .set_file_name(format!("{name}{suffix}.png"));
        if let Some(folder) = &self.state.folder {
            dialog = dialog.set_directory(folder);
        }
        let Some(file) = dialog.save_file() else {
            return;
        };
        let mut path = file.to_string_lossy().into_owned();
        // No extension given: PNG.
        if file.extension().is_none() {
            path.push_str(".png");
        }
        match crate::gel_export::save(&path, &image) {
            Ok(()) => self.state.status = Some(format!("Gel saved to {path}")),
            Err(e) => self.status_error(e),
        }
    }

    /// The whole map of the active document, as on screen (layers, style, enzymes and
    /// background) but unzoomed and without a selection, written to `path`: SVG (vector, text
    /// kept as text) or PNG (twice `width`×`height` pixels, for print), by its extension.
    #[cfg(not(target_arch = "wasm32"))]
    pub fn save_map(&mut self, path: &str, width: u32, height: u32) -> Result<Value, String> {
        let lower = path.to_ascii_lowercase();
        let svg = lower.ends_with(".svg");
        if !svg && !lower.ends_with(".png") {
            return Err("save the map as .svg or .png".into());
        }
        let bytes = self.map_file(svg, width, height)?;
        std::fs::write(path, bytes).map_err(|e| format!("cannot write {path}: {e}"))?;
        Ok(json!({
            "path": path,
            "format": if svg { "svg" } else { "png" },
            "width": width,
            "height": height,
        }))
    }

    /// The map file [`App::save_map`] writes: SVG text or PNG bytes.
    pub fn map_file(&mut self, svg: bool, width: u32, height: u32) -> Result<Vec<u8>, String> {
        if !(200..=4000).contains(&width) || !(200..=4000).contains(&height) {
            return Err("width and height must be 200-4000 points".into());
        }
        // Without a window there is no system theme: "system" keeps the current darkness.
        self.tokens = Tokens::for_choice(self.state.theme, self.tokens.dark);
        let t = self.canvas_tokens();
        let layers = self.layers();
        let style = self.state.map_style;
        let query = self.enzyme_query();
        let marks = self.site_cache.get(&self.session, &query);
        let doc = self.session.document().ok_or("no document is open")?;
        let m = doc.molecule.clone();
        let draw = |ui: &mut egui::Ui| {
            let rect = ui.max_rect();
            ui.painter().rect_filled(rect, 0.0, t.canvas);
            let _ = map::show(ui, &t, &m, None, style, &marks, layers);
        };
        if svg {
            let shapes = crate::offscreen::shapes(width, height, |ctx| t.apply(ctx), draw)?;
            let text = crate::svg::to_svg(
                &shapes,
                vec2(width as f32, height as f32),
                t.canvas,
                &format!("{} map", m.name()),
            );
            Ok(text.into_bytes())
        } else {
            let image =
                crate::offscreen::render(width * 2, height * 2, 2.0, |ctx| t.apply(ctx), draw)?;
            crate::offscreen::png_bytes(&image)
        }
    }

    /// The default size of an exported map, in points: square for circular maps, wide for
    /// linear ones.
    fn map_export_size(&self) -> (u32, u32) {
        let circular = match self.state.map_style {
            MapStyle::Auto => self.session.document().is_some_and(|d| {
                d.molecule.topology() == snappiesnap_engine::doc::Topology::Circular
            }),
            MapStyle::Circular => true,
            MapStyle::Linear => false,
        };
        if circular { (1000, 1000) } else { (1400, 500) }
    }

    /// Ask where to save the map, then save it.
    fn save_map_dialog(&mut self) {
        #[cfg(not(target_arch = "wasm32"))]
        {
            let name = self
                .session
                .document()
                .map(|d| format!("{} map.svg", d.molecule.name()))
                .unwrap_or_else(|| "map.svg".into());
            if let Some(path) = rfd::FileDialog::new()
                .add_filter("SVG (vector, for papers and posters)", &["svg"])
                .add_filter("PNG image", &["png"])
                .set_file_name(name)
                .save_file()
            {
                let mut path = path.to_string_lossy().into_owned();
                let lower = path.to_ascii_lowercase();
                if !lower.ends_with(".svg") && !lower.ends_with(".png") {
                    path.push_str(".svg");
                }
                let (w, h) = self.map_export_size();
                match self.save_map(&path, w, h) {
                    Ok(_) => self.state.status = Some(format!("Saved the map to {path}")),
                    Err(e) => self.status_error(e),
                }
            }
        }
        #[cfg(target_arch = "wasm32")]
        self.web_save_prompt(crate::web_files::SaveTarget::Map);
    }

    /// `ui.export_gel`: save the open document's gel, or an experiment tab's, as an image.
    #[cfg(not(target_arch = "wasm32"))]
    fn export_gel_request(&mut self, params: &Value) -> Result<Value, CmdError> {
        let path = params.get("path").and_then(Value::as_str).ok_or_else(|| {
            CmdError::invalid_params(
                "ui.export_gel needs {\"path\": \"/tmp/gel.png\"} (.png or .jpg)",
            )
        })?;
        let invert = params
            .get("invert")
            .and_then(Value::as_bool)
            .unwrap_or(false);
        let source = match params.get("experiment") {
            None | Some(Value::Null) => GelSource::Gel,
            Some(v) => GelSource::Experiment(
                v.as_u64()
                    .and_then(|i| usize::try_from(i).ok())
                    .ok_or_else(|| CmdError::invalid_params("experiment is a tab index"))?,
            ),
        };
        let (image, _) = self
            .gel_image(source, invert)
            .map_err(|e| CmdError::new(ErrorCode::NotFound, e))?;
        crate::gel_export::save(path, &image).map_err(|e| CmdError::new(ErrorCode::Io, e))?;
        let [w, h] = image.size;
        Ok(json!({ "path": path, "width": w, "height": h, "inverted": invert }))
    }

    /// Carry out what the Primers tab's lab-oligo section asked for.
    fn oligo_action(&mut self, action: crate::oligo_view::OligoAction) {
        use crate::oligo_view::OligoAction;
        match action {
            #[cfg(not(target_arch = "wasm32"))]
            OligoAction::Import | OligoAction::ImportFolder => {
                let folder = action == OligoAction::ImportFolder;
                let mut dialog = rfd::FileDialog::new().set_title(if folder {
                    "Import Oligo Folder"
                } else {
                    "Import Oligo List"
                });
                if !folder {
                    dialog = dialog
                        .add_filter("Oligo list (text, FASTA)", &["txt", "tsv", "fasta", "fa"])
                        .add_filter("All files", &["*"]);
                }
                if let Some(open) = &self.state.folder {
                    dialog = dialog.set_directory(open);
                }
                let picked = if folder {
                    dialog.pick_folder()
                } else {
                    dialog.pick_file()
                };
                let Some(picked) = picked else {
                    return;
                };
                let path = picked.to_string_lossy().into_owned();
                if let Some(v) = self.run("oligo_library_import", json!({ "path": path })) {
                    self.state.status = Some(Self::oligo_import_status(&v));
                }
                self.oligo_view.generation = self.oligo_view.generation.wrapping_add(1);
            }
            OligoAction::Add(names) => {
                if let Some(v) = self.run("oligo_library_add", json!({ "names": names })) {
                    let n = v.get("added").and_then(Value::as_array).map_or(0, Vec::len);
                    self.state.status = Some(format!("{n} lab oligos added as primers"));
                }
            }
            OligoAction::Select(range) => {
                self.run("select", json!({ "range": range }));
            }
            OligoAction::RemoveAll => {
                self.run("oligo_library_remove", json!({ "all": true }));
                self.oligo_view.generation = self.oligo_view.generation.wrapping_add(1);
            }
            #[allow(unreachable_patterns)]
            _ => {}
        }
    }

    /// The engine query for the enzymes the map shows. An enzyme set that no longer exists
    /// (deleted, or its file lost) falls back to the whole database.
    fn enzyme_query(&self) -> snappiesnap_engine::enzymes::SiteQuery {
        let mut q = self.state.enzyme_set.query();
        if q.set
            .as_ref()
            .is_some_and(|n| self.session.enzymes().set(n).is_none())
        {
            q.set = None;
        }
        q
    }

    /// What an oligo import did, for the status bar.
    #[cfg(not(target_arch = "wasm32"))]
    fn oligo_import_status(v: &Value) -> String {
        let count = |k: &str| v.get(k).and_then(Value::as_u64).unwrap_or(0);
        let list = |k: &str| v.get(k).and_then(Value::as_array).map_or(0, Vec::len);
        let mut status = String::from("Lab oligos: ");
        if let Some(files) = v.get("files").and_then(Value::as_u64) {
            status.push_str(&format!("{files} files read, "));
        }
        status.push_str(&format!(
            "{} added, {} already in the library ({} in all)",
            count("added"),
            count("duplicates"),
            count("total")
        ));
        let conflicts = list("conflicts");
        if conflicts > 0 {
            status.push_str(&format!(
                "; {conflicts} names already used for another sequence were kept as they were"
            ));
        }
        let skipped = count("skipped_count");
        if skipped > 0 {
            status.push_str(&format!("; {skipped} lines weren't oligos"));
        }
        let unreadable = list("unreadable");
        if unreadable > 0 {
            status.push_str(&format!("; {unreadable} files held no oligos"));
        }
        status
    }

    /// Open the sample vector and insert and set experiment `i` up to clone them: a first
    /// cloning experiment for someone with nothing of their own open, and the only way to try
    /// one in the browser, where a page has no files to pick from. Unsaved work is never
    /// replaced: `sample_open` refuses, and the error reaches the status bar like any other.
    fn open_sample_cloning(&mut self, i: usize) {
        let Some(v) = self.run("sample_open", json!({ "which": "cloning" })) else {
            return;
        };
        let ids: Vec<u64> = v
            .get("documents")
            .and_then(Value::as_array)
            .map(|d| {
                d.iter()
                    .filter_map(|x| x.get("id").and_then(Value::as_u64))
                    .collect()
            })
            .unwrap_or_default();
        // The two the command just opened are the last two tabs, vector first.
        let (Some(&insert), Some(&vector)) = (ids.last(), ids.get(ids.len().wrapping_sub(2)))
        else {
            self.status_error("the sample constructs did not open".to_string());
            return;
        };
        let Some(exp) = self.state.experiments.get_mut(i) else {
            return;
        };
        exp.technique = crate::state::CloningTechnique::Restriction;
        exp.vector = crate::state::CloningSide {
            document: Some(vector),
            enzymes: "EcoRI, BamHI".into(),
            ..Default::default()
        };
        exp.insert = crate::state::CloningSide {
            document: Some(insert),
            ..Default::default()
        };
        exp.extra.clear();
        exp.step = crate::state::ExperimentStep::Vector;
        self.experiment_inputs.1.cloning = Default::default();
        self.state.status = Some(
            "Opened pSample and demoInsert: step through the vector and the insert to the \
             product."
                .to_string(),
        );
    }

    /// Ask for an enzyme set file, import it and show its enzymes on the map.
    #[cfg(not(target_arch = "wasm32"))]
    fn import_enzyme_set(&mut self) {
        let mut dialog = rfd::FileDialog::new()
            .set_title("Import Enzyme Set")
            .add_filter("Enzyme set (text)", &["txt"])
            .add_filter("All files", &["*"]);
        if let Some(folder) = &self.state.folder {
            dialog = dialog.set_directory(folder);
        }
        let Some(file) = dialog.pick_file() else {
            return;
        };
        self.use_enzyme_set_file(&file.to_string_lossy());
    }

    /// Import an enzyme set from a file already in reach (on disk, or stored in the session's
    /// file system from the browser's picker) and show its enzymes on the map.
    fn use_enzyme_set_file(&mut self, path: &str) {
        let Some(v) = self.run("enzyme_set_import", json!({ "path": path })) else {
            return;
        };
        let name = v.get("name").and_then(Value::as_str).unwrap_or_default();
        let count = v.get("count").and_then(Value::as_u64).unwrap_or(0);
        let missing = v
            .get("missing")
            .and_then(Value::as_array)
            .map_or(0, Vec::len);
        self.state.enzyme_set.pool = Some(name.to_string());
        self.state.enzyme_set.only.clear();
        let mut status =
            format!("Enzyme set \u{201c}{name}\u{201d}: {count} enzymes, now shown on the map");
        if missing > 0 {
            let hint = if self.session.enzymes().is_builtin() {
                "; download REBASE in the Enzymes tab to include them"
            } else {
                ""
            };
            status.push_str(&format!(
                ". {missing} of them aren't in the enzyme database{hint}"
            ));
        }
        self.state.status = Some(status);
    }

    /// Keep the enzymes ticked in the Enzymes view as a set of the user's own, and show it
    /// on the map. Saving under a name that exists replaces that set, which is how a set is
    /// edited. On the web the sets live in memory, so they last the visit.
    fn save_enzyme_set(&mut self, name: String, enzymes: Vec<String>) {
        if enzymes.is_empty() {
            self.status_error("Tick the enzymes the set should hold before saving it.".to_string());
            return;
        }
        let count = enzymes.len();
        let Some(v) = self.run(
            "enzyme_set_create",
            json!({ "name": name, "enzymes": enzymes }),
        ) else {
            return;
        };
        let replaced = v.get("replaced").and_then(Value::as_bool).unwrap_or(false);
        self.state.enzyme_set.pool = Some(name.clone());
        self.state.enzyme_set.only.clear();
        let kept = if cfg!(target_arch = "wasm32") {
            ", kept for this visit"
        } else {
            ""
        };
        self.state.status = Some(format!(
            "Enzyme set \u{201c}{name}\u{201d} {}: {count} enzymes, now shown on the map{kept}",
            if replaced { "replaced" } else { "saved" }
        ));
    }

    /// Ask the browser for an enzyme set file. The bytes arrive in [`Self::web_sets`] and are
    /// imported on the next frame, like any other file the page is given.
    fn pick_enzyme_set_file(&mut self) {
        #[cfg(target_arch = "wasm32")]
        if let Err(e) = crate::web::pick_files(&self.web_sets, ".txt,.tsv,.csv", false, &self.ctx) {
            self.status_error(format!("the file picker didn't open: {e}"));
        }
    }

    /// Choose sequence files for the Align view: one for the pairwise alignment, or any number
    /// for the multiple alignment. The desktop opens its own picker; the browser hands the
    /// bytes to [`Self::web_align`], which the next frame copies into the session's file system,
    /// so both builds end up aligning a path.
    fn pick_align_files(&mut self, target: AlignFor) {
        let multiple = target != AlignFor::Pairwise;
        #[cfg(not(target_arch = "wasm32"))]
        {
            let dialog = rfd::FileDialog::new().add_filter(
                "Sequences",
                &[
                    "dna", "gb", "gbk", "genbank", "embl", "emb", "fa", "fasta", "fna", "fas",
                    "ab1", "abi", "scf",
                ],
            );
            let chosen: Vec<String> = if multiple {
                dialog.pick_files().unwrap_or_default()
            } else {
                dialog.pick_file().into_iter().collect()
            }
            .into_iter()
            .map(|p| p.to_string_lossy().into_owned())
            .collect();
            self.add_align_files(target, chosen);
        }
        #[cfg(target_arch = "wasm32")]
        {
            self.web_align_for = target;
            if let Err(e) = crate::web::pick_files(
                &self.web_align,
                crate::web_files::OPEN_ACCEPT,
                multiple,
                &self.ctx,
            ) {
                self.status_error(format!("the file picker didn't open: {e}"));
            }
        }
    }

    /// Put chosen files where the Align view will use them: ticked in the multiple alignment's
    /// list, or as the one file of the pairwise alignment.
    fn add_align_files(&mut self, target: AlignFor, paths: Vec<String>) {
        if target == AlignFor::AlignWith {
            self.align_with_files(paths);
        } else if target == AlignFor::Multiple {
            for path in paths {
                if !self.align.msa.files.contains(&path) {
                    self.align.msa.files.push(path);
                }
            }
        } else if let Some(path) = paths.into_iter().next() {
            self.align.file = path;
        }
    }

    /// Web build: take the sequence files the browser read for the Align view, copy each into
    /// the session's file system and use it, exactly as a file on disk is used.
    #[cfg_attr(not(target_arch = "wasm32"), allow(dead_code))]
    fn take_web_align_files(&mut self) {
        let target = self.web_align_for;
        let mut paths = Vec::new();
        for file in self.web_align.take() {
            if !sidebar::is_sequence_file(&file.name) {
                self.status_error(format!(
                    "{}: align .dna, GenBank, EMBL, FASTA or trace files",
                    file.name
                ));
                continue;
            }
            match file.bytes.and_then(|bytes| {
                crate::web_files::store(self.session.platform(), &file.name, &bytes)
            }) {
                Ok(path) => paths.push(path),
                Err(e) => self.status_error(e),
            }
        }
        self.add_align_files(target, paths);
    }

    /// Web build: take the read files the browser read for a Check Reads experiment, copy each
    /// into the session's file system and add it to that experiment, as a path on disk is.
    #[cfg_attr(not(target_arch = "wasm32"), allow(dead_code))]
    fn take_web_reads(&mut self) {
        let which = self.web_reads_experiment;
        for file in self.web_reads.take() {
            if !sidebar::is_sequence_file(&file.name) {
                self.status_error(format!(
                    "{}: reads are .ab1, .scf, FASTA or GenBank files",
                    file.name
                ));
                continue;
            }
            let stored = file.bytes.and_then(|bytes| {
                crate::web_files::store(self.session.platform(), &file.name, &bytes)
            });
            match stored {
                Ok(path) => {
                    if let Some(exp) = self.state.experiments.get_mut(which) {
                        crate::reads_view::add_read(exp, path);
                    }
                }
                Err(e) => self.status_error(e),
            }
        }
    }

    /// Web build: import the enzyme set files the browser has read since the last frame. The
    /// file is copied into the session's file system first, so the ordinary import command
    /// reads it by path, exactly as on the desktop.
    #[cfg_attr(not(target_arch = "wasm32"), allow(dead_code))]
    fn take_web_enzyme_sets(&mut self) {
        for file in self.web_sets.take() {
            let stored = file.bytes.and_then(|bytes| {
                crate::web_files::store(self.session.platform(), &file.name, &bytes)
            });
            match stored {
                Ok(path) => self.use_enzyme_set_file(&path),
                Err(e) => self.status_error(e),
            }
        }
    }

    fn run(&mut self, id: &str, params: Value) -> Option<Value> {
        match self.session.execute(id, params.clone()) {
            Ok(v) => {
                self.state.status = None;
                if matches!(id, "open" | "save") {
                    self.sidebar.refresh();
                    self.follow_document_folder();
                }
                if id == "open"
                    && let Some(n) = v
                        .get("matched_reads")
                        .and_then(Value::as_u64)
                        .filter(|n| *n > 0)
                {
                    let name = v.get("name").and_then(Value::as_str).unwrap_or("This file");
                    self.state.status = Some(format!(
                        "{name} has {n} matched sequencing reads: check them with Experiments \u{203a} Check Reads"
                    ));
                }
                Some(v)
            }
            Err(e) if e.code == ErrorCode::NotEnabled && e.message.contains("unsaved changes") => {
                let mut params = params;
                if let Some(obj) = params.as_object_mut() {
                    obj.insert("discard_changes".into(), Value::Bool(true));
                }
                self.confirm = Some(Confirm {
                    id: id.to_string(),
                    params,
                    question: format!("{}\nDiscard the unsaved changes?", e.message),
                });
                None
            }
            Err(e) => {
                self.status_error(e.to_string());
                None
            }
        }
    }

    fn status_error(&mut self, message: String) {
        self.session.log(Level::Warn, "ui", message.clone());
        self.state.status = Some(message);
    }

    #[cfg(not(target_arch = "wasm32"))]
    fn ui_state(&self) -> Value {
        json!({
            "state": self.state,
            "dark": self.tokens.dark,
            "title": self.title,
            "document_open": self.session.document().is_some(),
            "compare": self.compare.and_then(|id| self.session.documents().into_iter().find(|d| d.id == id).map(|d| d.index)),
        })
    }

    /// Answer waiting control-channel requests.
    #[cfg(not(target_arch = "wasm32"))]
    fn handle_control(&mut self, ctx: &egui::Context) {
        while let Some(req) = self.control.as_ref().and_then(|c| c.try_recv()) {
            let result = match req.method.as_str() {
                "ui.state" => Ok(self.ui_state()),
                "ui.set" => match serde_json::from_value::<UiPatch>(normalise(req.params)) {
                    Ok(patch) => {
                        if let Some(c) = patch.compare {
                            self.compare = usize::try_from(c)
                                .ok()
                                .and_then(|i| self.session.documents().into_iter().nth(i))
                                .map(|d| d.id);
                        }
                        self.state.apply(patch);
                        ctx.request_repaint();
                        Ok(self.ui_state())
                    }
                    Err(e) => Err(CmdError::invalid_params(format!("ui.set: {e}"))),
                },
                "ui.render" => self.render_request(&req.params),
                "ui.export_gel" => self.export_gel_request(&req.params),
                "ui.export_map" => {
                    let path = req.params.get("path").and_then(Value::as_str);
                    let (dw, dh) = self.map_export_size();
                    let size = |k: &str, d: u32| {
                        req.params
                            .get(k)
                            .and_then(Value::as_u64)
                            .and_then(|v| u32::try_from(v).ok())
                            .unwrap_or(d)
                    };
                    match path {
                        Some(p) => {
                            let p = p.to_string();
                            self.save_map(&p, size("width", dw), size("height", dh))
                                .map_err(|e| CmdError::new(ErrorCode::InvalidInput, e))
                        }
                        None => Err(CmdError::invalid_params(
                            "ui.export_map needs {\"path\": \"/tmp/map.svg\"} (.svg or .png)",
                        )),
                    }
                }
                m @ ("ui.open_experiment" | "ui.save_experiment") => {
                    self.experiment_request(m, &req.params)
                }
                "ui.screenshot" => {
                    let path = req
                        .params
                        .get("path")
                        .and_then(Value::as_str)
                        .map(str::to_string);
                    match path {
                        Some(path) => {
                            let reply = req.reply;
                            ctx.request_screenshot(move |image| {
                                let _ = reply.send(crate::control::save_png(&path, &image));
                            });
                            ctx.request_repaint();
                            continue;
                        }
                        None => Err(CmdError::invalid_params(
                            "ui.screenshot needs {\"path\": \"/tmp/shot.png\"}",
                        )),
                    }
                }
                method => api::call(&mut self.session, method, req.params),
            };
            let _ = req.reply.send(result);
            ctx.request_repaint();
        }
    }

    /// `ui.render`: draw the window offscreen into a PNG (works with a locked screen).
    #[cfg(not(target_arch = "wasm32"))]
    fn render_request(&mut self, params: &Value) -> Result<Value, CmdError> {
        let path = params.get("path").and_then(Value::as_str).ok_or_else(|| {
            CmdError::invalid_params("ui.render needs {\"path\": \"/tmp/render.png\"}")
        })?;
        let dim =
            |key: &str, default: u64| params.get(key).and_then(Value::as_u64).unwrap_or(default);
        let width = u32::try_from(dim("width", 1280)).unwrap_or(u32::MAX);
        let height = u32::try_from(dim("height", 820)).unwrap_or(u32::MAX);
        let scale = params.get("scale").and_then(Value::as_f64).unwrap_or(1.0) as f32;
        let image = self
            .render(width, height, scale)
            .map_err(|e| CmdError::invalid_params(format!("ui.render: {e}")))?;
        crate::control::save_png(path, &image)
    }

    /// Render the window as it would look, without a display.
    pub fn render(
        &mut self,
        width: u32,
        height: u32,
        scale: f32,
    ) -> Result<egui::ColorImage, String> {
        // Without a window there is no system theme: "system" keeps the current darkness.
        let tokens = Tokens::for_choice(self.state.theme, self.tokens.dark);
        self.tokens = tokens;
        crate::offscreen::render(
            width,
            height,
            scale,
            |ctx| tokens.apply(ctx),
            |ui| self.draw(ui),
        )
    }

    fn shortcuts(&mut self, ctx: &egui::Context) {
        let pressed = |m: Modifiers, k: Key| {
            ctx.input_mut(|i| i.consume_shortcut(&KeyboardShortcut::new(m, k)))
        };
        if pressed(Modifiers::COMMAND | Modifiers::SHIFT, Key::Z) {
            self.run("redo", json!({}));
        } else if pressed(Modifiers::COMMAND, Key::Z) {
            self.run("undo", json!({}));
        }
        if pressed(Modifiers::COMMAND, Key::N) {
            self.run("new", json!({}));
        }
        if pressed(Modifiers::COMMAND, Key::O) {
            self.open_dialog();
        }
        if pressed(Modifiers::COMMAND | Modifiers::SHIFT, Key::S) {
            self.save_dialog();
        } else if pressed(Modifiers::COMMAND, Key::S) {
            self.save();
        }
    }

    fn open_dialog(&mut self) {
        #[cfg(not(target_arch = "wasm32"))]
        if let Some(path) = rfd::FileDialog::new()
            .add_filter(
                "Sequence files",
                &[
                    "dna", "gb", "gbk", "genbank", "embl", "emb", "gff3", "gff", "sbol", "xml",
                    "fa", "fasta", "fna", "fas", "ab1", "abi", "scf", "sslab", "ssexp",
                ],
            )
            .add_filter("SnapGene", &["dna"])
            .add_filter("GenBank", &["gb", "gbk", "genbank"])
            .add_filter("EMBL", &["embl", "emb"])
            .add_filter("GFF3 with sequence", &["gff3", "gff"])
            .add_filter("SBOL 2", &["sbol", "xml", "rdf"])
            .add_filter("FASTA", &["fa", "fasta", "fna", "fas"])
            .add_filter("Sequencing traces", &["ab1", "abi", "scf"])
            .add_filter("SnappieSnap experiments", &["sslab", "ssexp"])
            .pick_file()
        {
            let path = path.to_string_lossy().into_owned();
            if sidebar::is_experiment_file(&path) {
                self.open_experiment(&path);
            } else {
                self.run("open", json!({ "path": path }));
            }
        }
        #[cfg(target_arch = "wasm32")]
        if let Err(e) = crate::web::pick_files(
            &self.web_inbox,
            crate::web_files::OPEN_ACCEPT,
            true,
            &self.ctx,
        ) {
            self.status_error(format!("the file picker didn't open: {e}"));
        }
    }

    /// Web build: open the files the browser has read since the last frame. Each is copied
    /// into the in-memory file system and opened from there, like a file on disk.
    #[cfg_attr(not(target_arch = "wasm32"), allow(dead_code))]
    fn take_web_files(&mut self) {
        for file in self.web_inbox.take() {
            let experiment = sidebar::is_experiment_file(&file.name);
            if !experiment && !sidebar::is_sequence_file(&file.name) {
                self.status_error(format!(
                    "{}: open .dna, GenBank, EMBL, GFF3, SBOL, FASTA, trace or .sslab files",
                    file.name
                ));
                continue;
            }
            let stored = file.bytes.and_then(|bytes| {
                crate::web_files::store(self.session.platform(), &file.name, &bytes)
            });
            match stored {
                Ok(path) if experiment => self.open_experiment(&path),
                Ok(path) => {
                    self.run("open", json!({ "path": path }));
                }
                Err(e) => self.status_error(e),
            }
        }
    }

    /// Web build: hand a file the app saved in the in-memory file system to the browser as a
    /// download.
    #[cfg(target_arch = "wasm32")]
    fn download_saved(&mut self, path: &str) {
        let name = crate::web_files::safe_name(path);
        let result = self
            .session
            .platform()
            .read(path)
            .map_err(|e| e.to_string())
            .and_then(|bytes| crate::web::download(&name, &bytes));
        match result {
            Ok(()) => self.state.status = Some(format!("Downloaded {name}")),
            Err(e) => self.status_error(format!("{name} couldn't be downloaded: {e}")),
        }
    }

    /// Open an experiment file in a tab of its own (or show it, if it is open already). Errors
    /// and warnings go to the status bar.
    fn open_experiment(&mut self, path: &str) {
        match self.try_open_experiment(path) {
            Ok((_, warnings)) if !warnings.is_empty() => self.status_error(warnings.join("; ")),
            Ok(_) => {}
            Err(e) => self.status_error(e),
        }
    }

    /// The experiment tab's index, and the warnings of opening its constructs.
    fn try_open_experiment(&mut self, path: &str) -> Result<(usize, Vec<String>), String> {
        if let Some(i) = self
            .state
            .experiments
            .iter()
            .position(|e| e.path.as_deref() == Some(path))
        {
            self.state.active_experiment = Some(i);
            return Ok((i, Vec::new()));
        }
        if self.state.experiments.len() >= crate::state::MAX_EXPERIMENTS {
            return Err("Too many experiment tabs: close one first".into());
        }
        let opened = self
            .session
            .execute("experiment_open", json!({ "path": path }))
            .map_err(|e| e.to_string())?;
        let exp = crate::state::Experiment::from_opened(path, &opened)
            .ok_or_else(|| format!("{path}: the experiment couldn't be shown"))?;
        self.state.experiments.push(exp);
        let i = self.state.experiments.len() - 1;
        self.state.active_experiment = Some(i);
        self.state.status = None;
        self.sidebar.refresh();
        Ok((i, strings(&opened, "warnings")))
    }

    fn save_experiment_dialog(&mut self, i: usize) {
        #[cfg(not(target_arch = "wasm32"))]
        {
            let Some(exp) = self.state.experiments.get(i) else {
                return;
            };
            let name = format!("{}.sslab", exp.name.trim());
            if let Some(path) = rfd::FileDialog::new()
                .add_filter("SnappieSnap experiment", &["sslab"])
                .set_file_name(name)
                .save_file()
            {
                let mut path = path.to_string_lossy().into_owned();
                if !sidebar::is_experiment_file(&path) {
                    path.push_str(".sslab");
                }
                self.save_experiment(i, &path);
            }
        }
        #[cfg(target_arch = "wasm32")]
        self.web_save_prompt(crate::web_files::SaveTarget::Experiment(i));
    }

    fn save_experiment(&mut self, i: usize, path: &str) {
        match self.try_save_experiment(i, path) {
            Ok(warnings) if warnings.is_empty() => {
                self.state.status = Some(format!("Saved the experiment to {path}"));
                #[cfg(target_arch = "wasm32")]
                self.download_saved(path);
            }
            Ok(warnings) => {
                #[cfg(target_arch = "wasm32")]
                self.download_saved(path);
                self.state.status = Some(format!("Saved to {path}. {}", warnings.join("; ")));
            }
            Err(e) => self.status_error(e),
        }
    }

    /// Save experiment `i` to `path`; the warnings of saving its constructs.
    fn try_save_experiment(&mut self, i: usize, path: &str) -> Result<Vec<String>, String> {
        let exp = self
            .state
            .experiments
            .get(i)
            .ok_or_else(|| format!("no experiment tab {i}"))?;
        let (params, key) = (exp.save_params(path), exp.content_key());
        let r = self
            .session
            .execute("experiment_save", params)
            .map_err(|e| e.to_string())?;
        if let Some(exp) = self.state.experiments.get_mut(i) {
            exp.path = Some(path.to_string());
            exp.saved = Some(key);
        }
        self.sidebar.refresh();
        Ok(strings(&r, "warnings"))
    }

    /// `ui.open_experiment {path}` and `ui.save_experiment {experiment?, path?}`: the same as the
    /// Experiments menu and the experiment tab's Save button.
    #[cfg(not(target_arch = "wasm32"))]
    fn experiment_request(&mut self, method: &str, params: &Value) -> Result<Value, CmdError> {
        let path = params
            .get("path")
            .and_then(Value::as_str)
            .map(str::to_string);
        if method == "ui.open_experiment" {
            let path = path.ok_or_else(|| {
                CmdError::invalid_params("ui.open_experiment needs {\"path\": \"run.sslab\"}")
            })?;
            let (i, warnings) = self
                .try_open_experiment(&path)
                .map_err(|e| CmdError::new(ErrorCode::InvalidInput, e))?;
            return Ok(json!({ "experiment": i, "warnings": warnings }));
        }
        let i = match params.get("experiment") {
            Some(v) => v
                .as_u64()
                .and_then(|i| usize::try_from(i).ok())
                .ok_or_else(|| CmdError::invalid_params("`experiment` is a tab index"))?,
            None => self.state.active_experiment.ok_or_else(|| {
                CmdError::invalid_params("no experiment tab is shown: give `experiment`")
            })?,
        };
        let path = path
            .or_else(|| self.state.experiments.get(i).and_then(|e| e.path.clone()))
            .ok_or_else(|| {
                CmdError::invalid_params("this experiment has never been saved: give `path`")
            })?;
        let warnings = self
            .try_save_experiment(i, &path)
            .map_err(|e| CmdError::new(ErrorCode::InvalidInput, e))?;
        Ok(json!({ "experiment": i, "path": path, "warnings": warnings }))
    }

    fn save(&mut self) {
        let has_path = self.session.document().is_some_and(|d| d.path.is_some());
        // Trace files are read only: an edited read is saved as a new file.
        let edited_trace = self.session.document().is_some_and(|d| {
            d.dirty
                && d.path.as_deref().is_some_and(|p| {
                    matches!(
                        snappiesnap_engine::io::format_for_path(p),
                        Some(
                            snappiesnap_engine::io::Format::Abif
                                | snappiesnap_engine::io::Format::Scf
                        )
                    )
                })
        });
        if has_path && !edited_trace {
            #[cfg_attr(not(target_arch = "wasm32"), allow(unused_variables))]
            if let Some(saved) = self.run("save", json!({})) {
                #[cfg(target_arch = "wasm32")]
                if let Some(path) = saved.get("path").and_then(Value::as_str) {
                    let path = path.to_string();
                    self.download_saved(&path);
                }
            }
        } else {
            self.save_dialog();
        }
    }

    fn save_dialog(&mut self) {
        #[cfg(not(target_arch = "wasm32"))]
        {
            let name = self
                .session
                .document()
                .map(|d| format!("{}.dna", d.molecule.name()))
                .unwrap_or_else(|| "Untitled.dna".into());
            if let Some(path) = rfd::FileDialog::new()
                .add_filter("SnapGene", &["dna"])
                .add_filter("GenBank", &["gb", "gbk"])
                .add_filter("EMBL", &["embl"])
                .add_filter("GFF3", &["gff3"])
                .add_filter("SBOL 2", &["sbol", "xml"])
                .add_filter("FASTA", &["fa", "fasta"])
                .set_file_name(name)
                .save_file()
            {
                self.run("save", json!({ "path": path.to_string_lossy() }));
            }
        }
        #[cfg(target_arch = "wasm32")]
        self.web_save_prompt(crate::web_files::SaveTarget::Document);
    }

    /// The window title: the open document (• when it has unsaved changes), then SnappieSnap
    /// and its version, so it is always clear which build is running.
    fn window_title(&self) -> String {
        let app = format!("SnappieSnap {}", env!("CARGO_PKG_VERSION"));
        match self.session.document() {
            Some(d) => {
                let name = d
                    .path
                    .as_deref()
                    .map(|p| p.rsplit(['/', '\\']).next().unwrap_or(p).to_string())
                    .unwrap_or_else(|| d.molecule.name().to_string());
                format!("{}{} — {app}", name, if d.dirty { " •" } else { "" })
            }
            None => app,
        }
    }

    fn toolbar(&mut self, ui: &mut egui::Ui) {
        let t = self.tokens;
        // While an experiment tab is shown, the buttons that act on a construct are off.
        let in_experiment = self.state.active_experiment.is_some();
        let has_doc = self.session.document().is_some() && !in_experiment;
        ui.horizontal(|ui| {
            ui.spacing_mut().item_spacing.x = 2.0;
            // The folder panel lists folders on disk, which a web page can't see.
            if cfg!(not(target_arch = "wasm32")) {
                if tool_button(ui, &t, Icon::Sidebar, "Files", true) {
                    self.state.sidebar = !self.state.sidebar;
                }
                separator(ui, &t);
            }
            if tool_button(ui, &t, Icon::Open, "Open", true) {
                self.open_dialog();
            }
            if tool_button(ui, &t, Icon::Save, "Save", has_doc) {
                self.save();
            }
            separator(ui, &t);
            if tool_button(
                ui,
                &t,
                Icon::Undo,
                "Undo",
                has_doc && self.session.can_undo(),
            ) {
                self.run("undo", json!({}));
            }
            if tool_button(
                ui,
                &t,
                Icon::Redo,
                "Redo",
                has_doc && self.session.can_redo(),
            ) {
                self.run("redo", json!({}));
            }
            separator(ui, &t);
            let nucleotide = has_doc
                && self
                    .session
                    .document()
                    .is_some_and(|d| d.molecule.seq().alphabet().is_nucleotide());
            if tool_button(ui, &t, Icon::ReverseComplement, "Rev. comp.", nucleotide) {
                self.run("reverse_complement", json!({}));
            }
            separator(ui, &t);
            // Experiments work across the open constructs, not on the active one.
            let experiments = tool_button_response(ui, &t, Icon::Experiments, "Experiments", true)
                .on_hover_text("Experiments across your open constructs");
            let mut chosen = None;
            let mut open_file = false;
            egui::Popup::menu(&experiments).show(|ui| {
                ui.set_min_width(280.0);
                for kind in crate::state::ExperimentKind::MENU {
                    let r = ui.add(egui::Button::new(kind.label()).min_size(vec2(270.0, 0.0)));
                    ui.label(
                        RichText::new(kind.description())
                            .small()
                            .color(t.text_muted),
                    );
                    if r.clicked() {
                        chosen = Some(kind);
                        ui.close();
                    }
                }
                ui.separator();
                if ui
                    .add(egui::Button::new("Open Experiment\u{2026}").min_size(vec2(270.0, 0.0)))
                    .clicked()
                {
                    open_file = true;
                    ui.close();
                }
                ui.label(
                    RichText::new("Reopen an experiment saved in a .sslab file")
                        .small()
                        .color(t.text_muted),
                );
            });
            #[cfg(target_arch = "wasm32")]
            if open_file
                && let Err(e) =
                    crate::web::pick_files(&self.web_inbox, ".sslab,.ssexp", false, &self.ctx)
            {
                self.status_error(format!("the file picker didn't open: {e}"));
            }
            #[cfg(not(target_arch = "wasm32"))]
            if open_file
                && let Some(path) = rfd::FileDialog::new()
                    .add_filter("SnappieSnap experiment", &["sslab", "ssexp"])
                    .pick_file()
            {
                self.open_experiment(&path.to_string_lossy());
            }
            if let Some(kind) = chosen {
                let mut open: Vec<u64> = self.session.documents().iter().map(|d| d.id).collect();
                // Check Reads takes the first document as its construct: the active one.
                if kind == crate::state::ExperimentKind::Reads
                    && let Some(active) = self.session.active_document_id()
                {
                    open.retain(|id| *id != active);
                    open.insert(0, active);
                }
                self.state.new_experiment(kind, open);
            }
            ui.with_layout(Layout::right_to_left(Align::Center), |ui| {
                let label = match self.state.theme {
                    ThemeChoice::System => "System",
                    ThemeChoice::Light => "Light",
                    ThemeChoice::Dark => "Dark",
                };
                if tool_button(ui, &t, Icon::Theme, label, true) {
                    self.state.theme = match self.state.theme {
                        ThemeChoice::System => ThemeChoice::Light,
                        ThemeChoice::Light => ThemeChoice::Dark,
                        ThemeChoice::Dark => ThemeChoice::System,
                    };
                }
                if tool_button(ui, &t, Icon::Preferences, "Preferences", true) {
                    self.state.preferences = !self.state.preferences;
                }
                if let Some(d) = self.session.document().filter(|_| !in_experiment) {
                    ui.add_space(t.spacing * 2.0);
                    ui.label(RichText::new(format!("{} bp", d.molecule.len())).color(t.text_muted));
                    // Long names are shortened rather than drawn over the toolbar buttons.
                    ui.add(egui::Label::new(RichText::new(d.molecule.name()).strong()).truncate())
                        .on_hover_text(d.molecule.name());
                }
            });
        });
    }

    /// When the active document changes, keep how the old one was viewed and restore how the
    /// new one was (a document seen for the first time starts on its map).
    fn sync_tab(&mut self) {
        let active = self.session.active_document_id();
        if active == self.current_tab {
            return;
        }
        // Which enzyme list the Enzymes view shows is a setting, not a per-document view.
        let enzyme_list = self.enzyme_view.from.clone();
        // Going to another document (switching, opening, creating) shows that document, not
        // the experiment tab.
        if self.current_tab.is_some() {
            self.state.active_experiment = None;
        }
        if let Some(old) = self.current_tab {
            let state = TabUi {
                view: self.state.view,
                map_style: self.state.map_style,
                gel: std::mem::take(&mut self.state.gel),
                seq_view: std::mem::take(&mut self.seq_view),
                enzyme_view: std::mem::take(&mut self.enzyme_view),
                primer_design: std::mem::take(&mut self.primer_design),
                gel_inputs: std::mem::take(&mut self.gel_inputs),
                align: std::mem::take(&mut self.align),
                find: std::mem::take(&mut self.find),
            };
            self.tab_states.insert(old, state);
        }
        let open: std::collections::HashSet<u64> =
            self.session.documents().iter().map(|d| d.id).collect();
        self.tab_states.retain(|id, _| open.contains(id));
        let first = self.current_tab.is_none();
        self.current_tab = active;
        let is_read = self
            .session
            .document()
            .is_some_and(|d| d.molecule.trace().is_some());
        let next = match active.and_then(|id| self.tab_states.remove(&id)) {
            Some(state) => state,
            // A sequencing read opens on its chromatogram.
            None if is_read => TabUi {
                view: View::Trace,
                ..TabUi::fresh()
            },
            // The first document keeps the view remembered from last time.
            None if first => return,
            None => TabUi::fresh(),
        };
        self.state.view = next.view;
        self.state.map_style = next.map_style;
        self.state.gel = next.gel;
        self.seq_view = next.seq_view;
        self.enzyme_view = next.enzyme_view;
        self.enzyme_view.from = enzyme_list;
        self.primer_design = next.primer_design;
        self.gel_inputs = next.gel_inputs;
        self.align = next.align;
        self.find = next.find;
        self.primer_edit = None;
        self.features.editor = None;
    }

    /// One tab per open document: its name, a dot when it has unsaved changes (or was never
    /// saved), and a close button. Clicking a tab makes it the active document.
    fn document_tabs(&mut self, ui: &mut egui::Ui) {
        let t = self.tokens;
        let docs = self.session.documents();
        let shown_experiment = self.state.active_experiment;
        let experiments: Vec<String> = self
            .state
            .experiments
            .iter()
            .map(|e| e.name.clone())
            .collect();
        let mut switch = None;
        let mut close = None;
        let mut new = false;
        let mut show_experiment = None;
        let mut close_experiment = None;
        egui::ScrollArea::horizontal()
            .id_salt("document-tabs")
            .show(ui, |ui| {
                ui.horizontal(|ui| {
                    ui.spacing_mut().item_spacing.x = 2.0;
                    for d in &docs {
                        let unsaved = d.dirty || d.path.is_none();
                        // The shown tab joins the toolbar below it (an experiment tab may be
                        // shown instead of the active document).
                        let active = d.active && shown_experiment.is_none();
                        let fill = if active { t.chrome } else { t.panel };
                        let stroke = if active {
                            Stroke::new(1.0, t.border)
                        } else {
                            Stroke::NONE
                        };
                        egui::Frame::new()
                            .fill(fill)
                            .stroke(stroke)
                            .corner_radius(CornerRadius {
                                nw: 6,
                                ne: 6,
                                sw: 0,
                                se: 0,
                            })
                            .inner_margin(egui::Margin::symmetric(10, 5))
                            .show(ui, |ui| {
                                ui.horizontal(|ui| {
                                    ui.spacing_mut().item_spacing.x = 6.0;
                                    let name: String = if d.name.chars().count() > 32 {
                                        d.name.chars().take(31).collect::<String>() + "\u{2026}"
                                    } else {
                                        d.name.clone()
                                    };
                                    let color = if active { t.text } else { t.text_muted };
                                    let label = ui.add(
                                        egui::Label::new(RichText::new(name).color(color))
                                            .sense(Sense::click()),
                                    );
                                    let hover = match (&d.path, d.dirty) {
                                        (None, _) => format!("{}\nNot saved yet", d.name),
                                        (Some(p), true) => format!("{p}\nUnsaved changes"),
                                        (Some(p), false) => p.clone(),
                                    };
                                    if label.on_hover_text(hover).clicked() {
                                        switch = Some(d.index);
                                    }
                                    if unsaved {
                                        let (r, _) =
                                            ui.allocate_exact_size(vec2(8.0, 8.0), Sense::hover());
                                        ui.painter().circle_filled(r.center(), 3.5, t.accent);
                                    }
                                    let x = ui
                                        .add(
                                            egui::Button::new(
                                                RichText::new("\u{d7}").color(t.text_muted),
                                            )
                                            .frame(false)
                                            .small(),
                                        )
                                        .on_hover_text("Close");
                                    if x.clicked() {
                                        close = Some(d.index);
                                    }
                                });
                            });
                    }
                    // Experiment tabs, after the documents: each spans the open constructs.
                    for (i, name) in experiments.iter().enumerate() {
                        let active = shown_experiment == Some(i);
                        egui::Frame::new()
                            .fill(if active { t.chrome } else { t.panel })
                            .stroke(if active {
                                Stroke::new(1.0, t.border)
                            } else {
                                Stroke::NONE
                            })
                            .corner_radius(CornerRadius {
                                nw: 6,
                                ne: 6,
                                sw: 0,
                                se: 0,
                            })
                            .inner_margin(egui::Margin::symmetric(10, 5))
                            .show(ui, |ui| {
                                ui.horizontal(|ui| {
                                    ui.spacing_mut().item_spacing.x = 6.0;
                                    let color = if active { t.text } else { t.text_muted };
                                    let (r, _) =
                                        ui.allocate_exact_size(vec2(14.0, 14.0), Sense::hover());
                                    icons::paint(ui.painter(), r, Icon::Experiments, color);
                                    let short: String = if name.chars().count() > 32 {
                                        name.chars().take(31).collect::<String>() + "\u{2026}"
                                    } else {
                                        name.clone()
                                    };
                                    let label = ui.add(
                                        egui::Label::new(RichText::new(short).color(color))
                                            .sense(Sense::click()),
                                    );
                                    if label
                                        .on_hover_text(format!(
                                            "{name}\nExperiment across the open constructs"
                                        ))
                                        .clicked()
                                    {
                                        show_experiment = Some(i);
                                    }
                                    let x = ui
                                        .add(
                                            egui::Button::new(
                                                RichText::new("\u{d7}").color(t.text_muted),
                                            )
                                            .frame(false)
                                            .small(),
                                        )
                                        .on_hover_text("Close the experiment (it isn't saved)");
                                    if x.clicked() {
                                        close_experiment = Some(i);
                                    }
                                });
                            });
                    }
                    let plus = ui
                        .add(
                            egui::Button::new(RichText::new("+").color(t.text_muted))
                                .frame(false)
                                .min_size(vec2(26.0, 24.0)),
                        )
                        .on_hover_text("New sequence (Cmd+N)");
                    if plus.clicked() {
                        new = true;
                    }
                });
            });
        if new {
            self.run("new", json!({}));
        }
        if let Some(i) = switch {
            self.state.active_experiment = None;
            self.run("switch_document", json!({ "index": i }));
        }
        if let Some(i) = show_experiment {
            self.state.active_experiment = Some(i);
        }
        if let Some(i) = close_experiment {
            self.state.close_experiment(i);
        }
        if let Some(i) = close {
            // Close acts on the active document: switch first. Unsaved work asks first.
            if self.run("switch_document", json!({ "index": i })).is_some() {
                self.run("close", json!({}));
            }
        }
    }

    fn tabs(&mut self, ui: &mut egui::Ui) {
        let t = self.tokens;
        let is_read = self
            .session
            .document()
            .is_some_and(|d| d.molecule.trace().is_some());
        ui.horizontal(|ui| {
            for view in View::ALL {
                if view == View::Trace && !is_read {
                    continue;
                }
                let selected = self.state.view == view;
                let text = RichText::new(view.label()).color(if selected {
                    t.accent_text
                } else {
                    t.text_muted
                });
                let r = ui.add(egui::Button::selectable(selected, text));
                if r.clicked() {
                    self.state.view = view;
                }
            }
        });
    }

    /// The composition of the selection (or the whole molecule) for the status bar, recomputed
    /// only when the document, its revision or the selection changed.
    fn refresh_info(&mut self) {
        let Some(doc) = self.session.document() else {
            self.info = None;
            return;
        };
        let key = (
            self.session.active_document_id(),
            self.session.revision(),
            doc.selection,
        );
        if self.info.as_ref().is_some_and(|(k, _)| *k == key) {
            return;
        }
        let parts = self
            .session
            .execute("sequence_info", json!({}))
            .map(|v| crate::goto_bar::summary(&v))
            .unwrap_or_default();
        self.info = Some((key, parts));
    }

    fn status_bar(&mut self, ui: &mut egui::Ui) {
        let t = self.tokens;
        let mut open_goto = false;
        ui.horizontal(|ui| {
            if let Some(msg) = &self.state.status {
                ui.label(RichText::new(msg).color(t.error));
            }
            ui.with_layout(Layout::right_to_left(Align::Center), |ui| {
                if let Some(e) = self
                    .state
                    .active_experiment
                    .and_then(|i| self.state.experiments.get(i))
                {
                    let set = e.lanes.iter().filter(|l| l.method.is_some()).count();
                    let saved = match (&e.path, e.unsaved()) {
                        (None, _) => "not saved".to_string(),
                        (Some(_), true) => "unsaved changes".to_string(),
                        (Some(p), false) => {
                            format!("saved as {}", crate::sidebar::base_name(p))
                        }
                    };
                    if e.kind == crate::state::ExperimentKind::Cloning {
                        let chosen = match (e.vector.is_set(), e.insert.is_set()) {
                            (false, _) => "no construct chosen yet".to_string(),
                            (true, false) => {
                                crate::cloning_experiment::side_label(&self.session, &e.vector)
                            }
                            (true, true) => format!(
                                "{} and {}",
                                crate::cloning_experiment::side_label(&self.session, &e.vector),
                                crate::cloning_experiment::side_label(&self.session, &e.insert)
                            ),
                        };
                        let text = format!(
                            "{}  \u{b7}  {}  \u{b7}  {chosen}",
                            e.kind.label(),
                            e.technique.label()
                        );
                        ui.label(RichText::new(text).color(t.text_muted));
                        return;
                    }
                    if e.kind == crate::state::ExperimentKind::Reads {
                        let matched = if e.matched {
                            " and the construct's matched reads"
                        } else {
                            ""
                        };
                        let text = format!(
                            "{}  \u{b7}  {} read files{matched}  \u{b7}  {saved}",
                            e.kind.label(),
                            e.reads.len()
                        );
                        ui.label(RichText::new(text).color(t.text_muted));
                        return;
                    }
                    let constructs = e
                        .lanes
                        .iter()
                        .map(|l| l.document)
                        .collect::<std::collections::BTreeSet<_>>()
                        .len();
                    let text = format!(
                        "{}  \u{b7}  {constructs} constructs  \u{b7}  {set} of {} lanes set up  \u{b7}  {saved}",
                        e.kind.label(),
                        e.lanes.len()
                    );
                    ui.label(RichText::new(text).color(t.text_muted));
                    return;
                }
                let Some(d) = self.session.document() else {
                    ui.label(RichText::new("No document").color(t.text_muted));
                    return;
                };
                let m = &d.molecule;
                let unit = if m.seq().alphabet().is_nucleotide() {
                    "bp"
                } else {
                    "aa"
                };
                let mut parts = vec![
                    format!("{} {unit}", m.len()),
                    format!("{:?}", m.topology()).to_lowercase(),
                    format!("{} features", m.features().len()),
                ];
                if let Some(sel) = d.selection.filter(|r| r.len_on(m.len()) > 0) {
                    parts.push(format!(
                        "selection {}..{} ({} {unit})",
                        sel.start.saturating_add(1),
                        sel.end,
                        sel.len_on(m.len())
                    ));
                }
                if let Some((_, info)) = &self.info {
                    parts.extend(info.iter().cloned());
                }
                if d.dirty {
                    parts.push("unsaved changes".into());
                }
                ui.label(RichText::new(parts.join("  ·  ")).color(t.text_muted));
                open_goto = ui
                    .small_button("Go to\u{2026}")
                    .on_hover_text("Select a position or range by its numbers (Cmd/Ctrl+L)")
                    .clicked();
            });
        });
        if open_goto {
            if !matches!(self.state.view, View::Map | View::Sequence) {
                self.state.view = View::Sequence;
            }
            self.goto.open();
        }
    }

    fn central(&mut self, ui: &mut egui::Ui) {
        // An experiment tab spans documents: it shows instead of the document views.
        if let Some(i) = self.state.active_experiment {
            let t = self.tokens;
            if self.experiment_inputs.0 != Some(i) {
                self.experiment_inputs = (Some(i), Default::default());
            }
            self.experiment_inputs.1.cloning.seq.font_size = self.state.sequence_font;
            use crate::experiment_view::ExperimentAction;
            let folder = self.state.folder.clone();
            // What the lab holds: the enzyme set the map is drawn from, so a diagnostic digest
            // is suggested from the enzymes at hand. A set that no longer exists is no pool.
            self.experiment_inputs.1.enzyme_pool = self
                .state
                .enzyme_set
                .pool
                .clone()
                .filter(|n| self.session.enzymes().set(n).is_some());
            // The construct pane of the Cloning experiment shows what a document tab shows,
            // with the same layers: taken by value here and written back after, because the
            // experiments themselves are borrowed from the same state.
            let mut pane = crate::cloning_experiment::Pane {
                sequence: self.state.cloning_sequence,
                single_line: self.state.sequence_single_line,
                layers: self.layers(),
                style: self.state.map_style,
                translation: self.state.translation,
                genetic_code: self.state.genetic_code,
                tails: self.state.primer_tails,
            };
            let action = self.state.experiments.get_mut(i).and_then(|exp| {
                crate::experiment_view::show(
                    ui,
                    &t,
                    &self.session,
                    exp,
                    &mut self.experiment_inputs.1,
                    folder.as_deref(),
                    &mut pane,
                )
            });
            self.state.cloning_sequence = pane.sequence;
            self.state.sequence_single_line = pane.single_line;
            self.state.features = pane.layers.features;
            self.state.enzymes = pane.layers.enzymes;
            self.state.primers = pane.layers.primers;
            self.state.orfs = pane.layers.orfs;
            self.state.ruler = pane.layers.ruler;
            self.state.positions = pane.layers.positions;
            self.state.map_style = pane.style;
            self.state.translation = pane.translation;
            match action {
                Some(ExperimentAction::SaveImage) => {
                    self.gel_export = Some(GelSource::Experiment(i));
                }
                Some(ExperimentAction::Save) => {
                    match self.state.experiments.get(i).and_then(|e| e.path.clone()) {
                        Some(path) => self.save_experiment(i, &path),
                        None => self.save_experiment_dialog(i),
                    }
                }
                Some(ExperimentAction::SaveAs) => self.save_experiment_dialog(i),
                Some(ExperimentAction::AddReads) => {
                    #[cfg(not(target_arch = "wasm32"))]
                    if let Some(files) = rfd::FileDialog::new()
                        .add_filter("Reads", crate::reads_view::READ_EXTENSIONS)
                        .pick_files()
                        && let Some(exp) = self.state.experiments.get_mut(i)
                    {
                        for f in files {
                            crate::reads_view::add_read(exp, f.to_string_lossy().into_owned());
                        }
                    }
                    // The browser reads the files itself; they arrive on a later frame.
                    #[cfg(target_arch = "wasm32")]
                    {
                        self.web_reads_experiment = i;
                        if let Err(e) = crate::web::pick_files(
                            &self.web_reads,
                            ".ab1,.abi,.scf,.fa,.fasta,.seq,.gb,.gbk",
                            true,
                            &self.ctx,
                        ) {
                            self.status_error(format!("the file picker didn't open: {e}"));
                        }
                    }
                }
                Some(ExperimentAction::Run(command, params)) => {
                    if self.run(command, params).is_some() {
                        // Show the read just opened.
                        self.state.active_experiment = None;
                    }
                }
                Some(ExperimentAction::CheckReads(params)) => {
                    let answer = self
                        .session
                        .execute("check_reads", params)
                        .map_err(|e| e.message);
                    self.experiment_inputs.1.reads.set_result(answer);
                }
                Some(ExperimentAction::CloningInfo(path)) => {
                    let answer = self
                        .session
                        .execute("file_info", json!({"path": path, "molecule": true}))
                        .map_err(|e| e.message);
                    self.experiment_inputs.1.cloning.set_info(path, answer);
                }
                Some(ExperimentAction::CloningSamples) => self.open_sample_cloning(i),
                Some(ExperimentAction::CloningPreview(command, params)) => {
                    let answer = self.session.execute(command, params).map_err(|e| e.message);
                    self.experiment_inputs.1.cloning.set_preview(answer);
                }
                Some(ExperimentAction::Suggest(params)) => {
                    let answer = self
                        .session
                        .execute("diagnostic_digest", params)
                        .map_err(|e| e.message);
                    self.experiment_inputs.1.set_suggestions(answer);
                }
                None => {}
            }
            return;
        }
        if self.session.document().is_none() {
            self.empty_state(ui);
            return;
        }
        match self.state.view {
            View::Map => self.map_view(ui),
            View::Sequence => self.sequence_view(ui),
            View::Features => self.features_table(ui),
            View::Primers => self.primers_table(ui),
            View::Enzymes => self.enzymes_view(ui),
            View::Align => self.align_view(ui),
            View::Gel => {
                let t = self.tokens;
                if crate::gel_view::show(
                    ui,
                    &t,
                    &self.session,
                    &mut self.state.gel,
                    &mut self.gel_inputs,
                    &mut self.gel_cache,
                ) {
                    self.gel_export = Some(GelSource::Gel);
                }
            }
            View::History => self.history(ui),
            View::Notes => self.notes_tab(ui),
            View::Trace => self.trace_tab(ui),
        }
    }

    fn trace_tab(&mut self, ui: &mut egui::Ui) {
        let t = self.canvas_tokens();
        let Some(doc) = self.session.document() else {
            return;
        };
        if doc.molecule.trace().is_none() {
            // Another document became active while the Trace tab was shown.
            self.sequence_view(ui);
            return;
        }
        let actions = crate::trace_view::show(
            ui,
            &t,
            &doc.molecule,
            doc.selection,
            &mut self.state.trace_zoom,
            &mut self.state.trace_gain,
            &mut self.trace_view,
        );
        for action in actions {
            match action {
                crate::trace_view::TraceAction::Select(range) => {
                    self.run("select", json!({ "range": range }));
                }
                crate::trace_view::TraceAction::Trim => {
                    if let Some(r) = self.run("trace_trim", json!({})) {
                        let a = r.get("removed_start").and_then(Value::as_u64).unwrap_or(0);
                        let b = r.get("removed_end").and_then(Value::as_u64).unwrap_or(0);
                        self.state.status = Some(format!(
                            "Trimmed {a} bases from the start and {b} from the end (Undo restores them)"
                        ));
                    }
                }
            }
        }
    }

    /// The map layers that are switched on (the sequence view shows the same ones).
    fn layers(&self) -> map::Layers {
        map::Layers {
            features: self.state.features,
            enzymes: self.state.enzymes,
            primers: self.state.primers,
            positions: self.state.positions,
            ruler: self.state.ruler,
            orfs: self.state.orfs,
            labels: true,
        }
    }

    /// Layer toggles down the left edge of the map and the sequence view. Right-clicking the
    /// ORF button picks the reading-frame translations of the sequence view.
    fn layer_column(&mut self, ui: &mut egui::Ui) {
        let t = self.tokens;
        ui.vertical(|ui| {
            ui.spacing_mut().item_spacing.y = 4.0;
            if layer_toggle(ui, &t, Icon::Features, "Features", self.state.features).clicked() {
                self.state.features = !self.state.features;
            }
            // Enzymes, with a menu picking which ones.
            let mut set_action = None;
            if let Some(pool) = &self.state.enzyme_set.pool
                && self.session.enzymes().set(pool).is_none()
            {
                self.state.enzyme_set.pool = None;
            }
            let db = self.session.enzymes();
            let enzymes_tip = format!(
                "Restriction sites: {} ({}, {} enzymes). The arrow or a right-click picks which enzymes",
                self.state.enzyme_set.describe(),
                db.version
                    .as_ref()
                    .map_or(db.source.clone(), |v| format!("{} {v}", db.source)),
                db.enzymes.len()
            );
            let unknown: Vec<String> = self
                .state
                .enzyme_set
                .only
                .iter()
                .filter(|n| db.get(n).is_none())
                .cloned()
                .collect();
            ui.horizontal(|ui| {
                ui.spacing_mut().item_spacing.x = 2.0;
                let enzymes =
                    layer_toggle(ui, &t, Icon::Enzymes, &enzymes_tip, self.state.enzymes);
                if enzymes.clicked() {
                    self.state.enzymes = !self.state.enzymes;
                }
                let chevron = chevron(ui, &t).on_hover_text("Which enzymes to show");
                let state = &mut self.state;
                egui::Popup::menu(&chevron)
                    .close_behavior(egui::PopupCloseBehavior::CloseOnClickOutside)
                    .show(|ui| {
                        set_action = set_action.take().or(enzyme_menu(ui, &t, state, db, &unknown));
                    });
                enzymes.context_menu(|ui| {
                    set_action = set_action.take().or(enzyme_menu(ui, &t, state, db, &unknown));
                });
            });
            #[cfg(not(target_arch = "wasm32"))]
            if matches!(set_action, Some(EnzymeMenuAction::Import)) {
                self.import_enzyme_set();
            }
            if let Some(EnzymeMenuAction::Remove(name)) = set_action
                && self
                    .run("enzyme_set_delete", json!({ "name": name }))
                    .is_some()
            {
                self.state.enzyme_set.pool = None;
                self.state.status = Some(format!("Enzyme set \u{201c}{name}\u{201d} removed"));
            }
            for (icon, tip, on) in [
                (
                    Icon::Primers,
                    "Primer binding sites".to_string(),
                    &mut self.state.primers,
                ),
                (
                    Icon::Positions,
                    "Cut positions in enzyme labels".to_string(),
                    &mut self.state.positions,
                ),
                (
                    Icon::Ruler,
                    "Ruler (coordinates)".to_string(),
                    &mut self.state.ruler,
                ),
            ] {
                if layer_toggle(ui, &t, icon, &tip, *on).clicked() {
                    *on = !*on;
                }
            }
            // Only the sequence view wraps, so only it can stop wrapping.
            if self.state.view == View::Sequence
                && layer_toggle(
                    ui,
                    &t,
                    Icon::SingleLine,
                    "The whole sequence on one line, scrolled sideways",
                    self.state.sequence_single_line,
                )
                .clicked()
            {
                self.state.sequence_single_line = !self.state.sequence_single_line;
                // The reads line up under the one line: wrapping the sequence closes them.
                if !self.state.sequence_single_line {
                    self.state.align_panel = false;
                }
            }
            if self.state.view == View::Sequence
                && layer_toggle(
                    ui,
                    &t,
                    Icon::AlignWith,
                    "Align with: reads and other sequences aligned under this one, with their traces",
                    self.state.align_panel,
                )
                .clicked()
            {
                self.state.align_panel = !self.state.align_panel;
            }
            let shows = if self.state.translation == Translation::Orfs {
                "ORFs".to_string()
            } else {
                format!("translation, {}", self.state.translation.label())
            };
            let has_selection = self
                .session
                .document()
                .and_then(|d| d.selection)
                .is_some_and(|r| !r.is_empty());
            let mut menu_action = None;
            ui.horizontal(|ui| {
                ui.spacing_mut().item_spacing.x = 2.0;
                let orfs = layer_toggle(
                    ui,
                    &t,
                    Icon::Orfs,
                    &format!(
                        "ORFs and translations (the sequence view shows {shows}; the map, ORFs of 75+ codons). The arrow or a right-click picks what it shows",
                    ),
                    self.state.orfs,
                );
                if orfs.clicked() {
                    self.state.orfs = !self.state.orfs;
                }
                let chevron = chevron(ui, &t);
                let chevron = chevron.on_hover_text("What the ORF button shows");
                let state = &mut self.state;
                egui::Popup::menu(&chevron).show(|ui| {
                    menu_action = translation_menu(ui, &t, state, has_selection);
                });
                orfs.context_menu(|ui| {
                    menu_action = translation_menu(ui, &t, state, has_selection);
                });
            });
            if let Some(TranslationAction::AddFeature) = menu_action {
                self.add_translated_feature();
            }
        });
    }

    /// A CDS feature over the selection, on the strand of the chosen frames, translated with the
    /// chosen genetic code.
    fn add_translated_feature(&mut self) {
        let Some(range) = self.session.document().and_then(|d| d.selection) else {
            return;
        };
        let mut qualifiers = Vec::new();
        if self.state.genetic_code != 1 {
            qualifiers.push(
                json!({"name": "transl_table", "values": [self.state.genetic_code.to_string()]}),
            );
        }
        let strand = if self.state.translation.reverse() {
            "reverse"
        } else {
            "forward"
        };
        self.run(
            "add_feature",
            json!({"feature": {
                "name": "translated region",
                "kind": "CDS",
                "strand": strand,
                "segments": [{"range": range}],
                "qualifiers": qualifiers,
            }}),
        );
    }

    fn map_view(&mut self, ui: &mut egui::Ui) {
        let t = self.canvas_tokens();
        let mut export_map = false;
        // Find works on the sequence: Cmd+F switches to it.
        if ui.input(|i| i.modifiers.command && i.key_pressed(egui::Key::F)) {
            self.state.view = View::Sequence;
            self.find.open();
        }
        if ui.input(|i| i.modifiers.command && i.key_pressed(egui::Key::L)) {
            self.goto.open();
        }
        // The selection shows on the map; the sequence view scrolls to it when shown.
        if let Some(range) = crate::goto_bar::show(ui, &t, &mut self.session, &mut self.goto) {
            self.seq_view.scroll_to = Some(range.start);
        }
        // A compared document that was closed, or became the active one, is dropped.
        let active = self.session.active_document_id();
        if self
            .compare
            .is_some_and(|id| Some(id) == active || self.session.document_by_id(id).is_none())
        {
            self.compare = None;
        }
        let others: Vec<(u64, String)> = self
            .session
            .documents()
            .into_iter()
            .filter(|d| !d.active)
            .map(|d| (d.id, d.name))
            .collect();
        ui.horizontal(|ui| {
            // Zoom.
            if ui
                .small_button("\u{2212}")
                .on_hover_text("Zoom out")
                .clicked()
            {
                self.state.map_zoom = (self.state.map_zoom / 1.5).max(1.0);
            }
            ui.label(
                RichText::new(format!("{:.0}%", self.state.map_zoom * 100.0)).color(t.text_muted),
            )
            .on_hover_text("Zoom: pinch or Cmd/Ctrl+scroll over the map; drag to move it");
            if ui.small_button("+").on_hover_text("Zoom in").clicked() {
                self.state.map_zoom = (self.state.map_zoom * 1.5).min(crate::state::MAX_ZOOM);
            }
            if ui
                .add_enabled(self.state.map_zoom > 1.0, egui::Button::new("Fit").small())
                .clicked()
            {
                self.state.map_zoom = 1.0;
                self.map_pan = Vec2::ZERO;
            }
            if ui
                .small_button("Export\u{2026}")
                .on_hover_text("Save the whole map as SVG (vector, editable text) or PNG")
                .clicked()
            {
                export_map = true;
            }
            if !others.is_empty() {
                ui.separator();
                ui.label(RichText::new("Compare with").color(t.text_muted));
                let shown = self
                    .compare
                    .and_then(|id| others.iter().find(|(o, _)| *o == id))
                    .map_or("nothing".to_string(), |(_, n)| n.clone());
                egui::ComboBox::from_id_salt("map-compare")
                    .width(220.0)
                    .selected_text(shown)
                    .show_ui(ui, |ui| {
                        ui.selectable_value(&mut self.compare, None, "nothing");
                        for (id, name) in &others {
                            ui.selectable_value(&mut self.compare, Some(*id), name);
                        }
                    });
            }
            ui.with_layout(Layout::right_to_left(Align::Center), |ui| {
                for (style, label) in [
                    (MapStyle::Linear, "Linear"),
                    (MapStyle::Circular, "Circular"),
                    (MapStyle::Auto, "Auto"),
                ] {
                    let selected = self.state.map_style == style;
                    let text = RichText::new(label).color(if selected {
                        t.accent_text
                    } else {
                        t.text_muted
                    });
                    if ui.add(egui::Button::selectable(selected, text)).clicked() {
                        self.state.map_style = style;
                    }
                }
            });
        });
        if export_map {
            self.save_map_dialog();
        }
        let layers = self.layers();
        let style = self.state.map_style;
        let query = self.enzyme_query();
        let mut action = None;
        ui.horizontal_top(|ui| {
            self.layer_column(ui);
            let area = ui.available_rect_before_wrap();
            ui.allocate_rect(area, Sense::hover());
            let Some(doc) = self.session.document() else {
                return;
            };
            let circular = match style {
                MapStyle::Auto => {
                    doc.molecule.topology() == snappiesnap_engine::doc::Topology::Circular
                }
                MapStyle::Circular => true,
                MapStyle::Linear => false,
            };
            let compared = self
                .compare
                .and_then(|id| self.session.document_by_id(id).map(|d| (id, d)));
            // Side by side: the active document left, the compared one right (same zoom).
            let (left, right) = match compared {
                Some(_) => {
                    let w = (area.width() - 12.0) / 2.0;
                    (
                        Rect::from_min_size(area.min, vec2(w, area.height())),
                        Some(Rect::from_min_size(
                            area.min + vec2(w + 12.0, 0.0),
                            vec2(w, area.height()),
                        )),
                    )
                }
                None => (area, None),
            };
            let marks = self.site_cache.get(&self.session, &query);
            let zoom = &mut self.state.map_zoom;
            let pan = &mut self.map_pan;
            action = zoomed(ui, "map-main", left, zoom, pan, circular, |ui| {
                map::show(ui, &t, &doc.molecule, doc.selection, style, &marks, layers)
            });
            if let (Some((id, other)), Some(right)) = (compared, right) {
                ui.painter().line_segment(
                    [
                        egui::pos2(right.left() - 6.0, right.top()),
                        egui::pos2(right.left() - 6.0, right.bottom()),
                    ],
                    Stroke::new(1.0, t.border),
                );
                let marks =
                    self.compare_cache
                        .get_for(&self.session, Some(id), Some(other), &query);
                let other_circular = match style {
                    MapStyle::Auto => {
                        other.molecule.topology() == snappiesnap_engine::doc::Topology::Circular
                    }
                    MapStyle::Circular => true,
                    MapStyle::Linear => false,
                };
                // The compared map follows the main one's zoom; it is for looking, not editing.
                let _ = zoomed(ui, "map-compare", right, zoom, pan, other_circular, |ui| {
                    map::show(
                        ui,
                        &t,
                        &other.molecule,
                        other.selection,
                        style,
                        &marks,
                        layers,
                    )
                });
            }
        });
        // Choosing an enzyme or a primer is for the cloning experiment; here a map click
        // selects.
        if let Some(MapAction::Select(range)) = action {
            self.run("select", json!({ "range": range }));
        }
    }

    fn sequence_view(&mut self, ui: &mut egui::Ui) {
        let t = self.canvas_tokens();
        if self.session.document().is_none() {
            return;
        }
        let revision = self.session.revision();
        // The ORF layer shows ORFs or translated frames, as picked in its menu.
        let mut layers = self.layers();
        layers.orfs = self.state.orfs && self.state.translation == Translation::Orfs;
        let frames: &[i8] = if self.state.orfs {
            self.state.translation.frames()
        } else {
            &[]
        };
        let table = self.state.genetic_code;
        if ui.input(|i| i.modifiers.command && i.key_pressed(egui::Key::F)) {
            self.find.open();
        }
        if ui.input(|i| i.modifiers.command && i.key_pressed(egui::Key::L)) {
            self.goto.open();
        }
        if let Some(range) = crate::goto_bar::show(ui, &t, &mut self.session, &mut self.goto) {
            self.seq_view.scroll_to = Some(range.start);
        }
        if let Some(crate::find_bar::FindAction::Show(range)) =
            crate::find_bar::show(ui, &t, &mut self.session, &mut self.find)
            && self.run("select", json!({ "range": range })).is_some()
        {
            self.seq_view.scroll_to = Some(range.start);
        }
        let found: Vec<snappiesnap_engine::doc::Range> = if self.find.open {
            self.find.matches.iter().map(|(r, _)| *r).collect()
        } else {
            Vec::new()
        };
        let mut actions = Vec::new();
        // With Align with open, the view is split in fixed parts: the sequence on one line in
        // the top half, the reads on the same columns in the bottom half, and one navigator
        // under both. The parts are laid out here, not sized by what they contain, so nothing
        // can creep from frame to frame.
        // However Align with was opened (toggle, menu, a dropped file, the control channel),
        // the sequence goes on one line for it.
        if self.state.align_panel {
            self.state.sequence_single_line = true;
        }
        let aligning = self.state.align_panel;
        self.seq_view.navigator_outside = aligning;
        self.seq_view.font_size = self.state.sequence_font;
        if !aligning {
            self.seq_view.set_reads(Vec::new());
        }
        let full = ui.available_rect_before_wrap();
        let mut seq_rect = full;
        let mut lower: Option<(Rect, Rect)> = None;
        if aligning
            && let Some(id) = self.session.active_document_id()
            && let Some(doc) = self.session.document()
        {
            let entries = self.state.align_list(id).cloned().unwrap_or_default();
            self.seq_view.set_reads(crate::align_with::read_spans(
                &self.align_with,
                &entries,
                doc.molecule.len(),
            ));
            let nav_h = self
                .seq_view
                .navigator_height(&doc.molecule, revision, layers.features)
                .min(full.height());
            let (top, bottom) = split_for_alignment(full, nav_h);
            seq_rect = top;
            lower = Some(bottom);
        }
        ui.scope_builder(
            egui::UiBuilder::new()
                .max_rect(seq_rect)
                .layout(Layout::left_to_right(Align::Min)),
            |ui| {
                ui.shrink_clip_rect(seq_rect);
                if self.state.align_panel {
                    self.align_with_panel(ui);
                }
                self.layer_column(ui);
                let query = self.enzyme_query();
                let marks = self.site_cache.get(&self.session, &query);
                let Some(doc) = self.session.document() else {
                    return;
                };
                // The sequence takes all the room right of the toggles.
                let size = ui.available_size();
                ui.allocate_ui_with_layout(size, Layout::top_down(Align::Min), |ui| {
                    ui.set_min_size(size);
                    actions = seqview::show(
                        ui,
                        &t,
                        &doc.molecule,
                        doc.selection,
                        &mut self.seq_view,
                        revision,
                        &marks,
                        layers,
                        frames,
                        table,
                        &found,
                        self.state.primer_tails,
                        self.state.sequence_single_line,
                    );
                });
            },
        );
        if let Some((pane_rect, nav_rect)) = lower {
            ui.scope_builder(egui::UiBuilder::new().max_rect(pane_rect), |ui| {
                ui.shrink_clip_rect(pane_rect);
                self.align_with_pane(ui);
            });
            ui.scope_builder(egui::UiBuilder::new().max_rect(nav_rect), |ui| {
                ui.shrink_clip_rect(nav_rect);
                if let Some(doc) = self.session.document() {
                    self.seq_view.show_navigator(
                        ui,
                        &t,
                        &doc.molecule,
                        layers.features,
                        nav_rect.height(),
                        doc.selection,
                    );
                }
            });
            // A swipe over the reads moves the sequence, and the reads with it.
            self.seq_view.scroll_by += std::mem::take(&mut self.align_with.scroll_by);
        }
        ui.advance_cursor_after_rect(full);
        for action in actions {
            match action {
                SeqAction::Select(range) => {
                    self.run("select", json!({ "range": range }));
                }
                // Choosing an enzyme is for the cloning experiment; a document tab ignores it.
                SeqAction::SiteChosen(_) => {}
                SeqAction::Insert { at, residues } => {
                    let n = residues.len();
                    if self
                        .run("insert", json!({ "at": at, "residues": residues }))
                        .is_some()
                    {
                        self.run(
                            "select",
                            json!({ "range": { "start": at + n, "end": at + n } }),
                        );
                    }
                }
                SeqAction::Replace { range, residues } => {
                    let n = residues.len();
                    if self
                        .run("replace", json!({ "range": range, "residues": residues }))
                        .is_some()
                    {
                        let at = if range.is_wrapping() {
                            0
                        } else {
                            range.start + n
                        };
                        self.run("select", json!({ "range": { "start": at, "end": at } }));
                    }
                }
                SeqAction::Delete(range) => {
                    if self.run("delete", json!({ "range": range })).is_some() {
                        let at = if range.is_wrapping() { 0 } else { range.start };
                        self.run("select", json!({ "range": { "start": at, "end": at } }));
                    }
                }
                SeqAction::Copy(format) => {
                    self.copy_selection(ui.ctx(), format);
                }
                SeqAction::Blast(against) => {
                    if let Some(r) = self.run("blast_search", json!({ "against": against })) {
                        let url = r.get("url").and_then(Value::as_str).unwrap_or("");
                        let program = r.get("program").and_then(Value::as_str).unwrap_or("BLAST");
                        let n = r.get("query_length").and_then(Value::as_u64).unwrap_or(0);
                        match crate::browser::open(ui.ctx(), url) {
                            Ok(()) => {
                                self.state.status = Some(format!(
                                    "Opened NCBI {program} in your browser with {n} residues; press BLAST there to search"
                                ));
                            }
                            Err(e) => self.status_error(e),
                        }
                    }
                }
                SeqAction::Cut => {
                    let range = self.session.document().and_then(|d| d.selection);
                    if self.copy_selection(ui.ctx(), "plain")
                        && let Some(range) = range
                        && self.run("delete", json!({ "range": range })).is_some()
                    {
                        let at = if range.is_wrapping() { 0 } else { range.start };
                        self.run("select", json!({ "range": { "start": at, "end": at } }));
                    }
                }
                SeqAction::Paste(text) => {
                    if let Some(r) = self.run("paste", json!({ "text": text })) {
                        let n = r.get("pasted").and_then(Value::as_u64).unwrap_or(0);
                        let notes: Vec<&str> = r
                            .get("notes")
                            .and_then(Value::as_array)
                            .map(|a| a.iter().filter_map(Value::as_str).collect())
                            .unwrap_or_default();
                        let mut msg = format!("Pasted {n} residues");
                        if !notes.is_empty() {
                            msg.push_str(&format!(" ({})", notes.join("; ")));
                        }
                        self.state.status = Some(msg);
                    }
                }
            }
        }
    }

    /// The "Align with" pane under the sequence, in the room it is given.
    fn align_with_pane(&mut self, ui: &mut egui::Ui) {
        let Some(id) = self.session.active_document_id() else {
            return;
        };
        let t = self.tokens;
        let mut entries = self
            .state
            .align_list(id)
            .map(std::mem::take)
            .unwrap_or_default();
        let action = egui::Frame::new()
            .fill(t.panel)
            .inner_margin(egui::Margin::symmetric(0, 4))
            .show(ui, |ui| {
                ui.set_min_size(ui.available_size());
                crate::align_with::pane(
                    ui,
                    &t,
                    &self.session,
                    id,
                    &mut entries,
                    &mut self.align_with,
                    self.seq_view.line,
                )
            })
            .inner;
        if let Some(list) = self.state.align_list(id) {
            *list = entries;
        }
        self.align_with_action(action);
    }

    /// The "Align with" list, left of the sequence.
    fn align_with_panel(&mut self, ui: &mut egui::Ui) {
        let Some(id) = self.session.active_document_id() else {
            return;
        };
        let t = self.tokens;
        let dragging = ui.ctx().input(|i| !i.raw.hovered_files.is_empty());
        let mut entries = self
            .state
            .align_list(id)
            .map(std::mem::take)
            .unwrap_or_default();
        let action = crate::align_with::panel(
            ui,
            &t,
            &self.session,
            id,
            &mut entries,
            &mut self.align_with,
            dragging,
        );
        if let Some(list) = self.state.align_list(id) {
            *list = entries;
        }
        self.align_with_action(action);
    }

    fn align_with_action(&mut self, action: Option<crate::align_with::AlignWithAction>) {
        use crate::align_with::AlignWithAction;
        match action {
            Some(AlignWithAction::AddFiles) => self.pick_align_files(AlignFor::AlignWith),
            Some(AlignWithAction::Check(params)) => {
                let answer = self
                    .session
                    .execute("check_reads", params)
                    .map_err(|e| e.message);
                self.align_with.set_result(answer);
            }
            Some(AlignWithAction::Select(p)) => {
                if self
                    .run("select", json!({ "range": { "start": p, "end": p + 1 } }))
                    .is_some()
                {
                    self.seq_view.scroll_to = Some(p);
                }
            }
            Some(AlignWithAction::Close) => self.state.align_panel = false,
            None => {}
        }
    }

    /// Add files to the active document's "Align with" list.
    fn align_with_files(&mut self, paths: Vec<String>) {
        let Some(id) = self.session.active_document_id() else {
            return;
        };
        self.state.align_panel = true;
        for path in paths {
            if !self
                .state
                .align_add(id, crate::state::AlignSource::File(path))
            {
                self.status_error(format!(
                    "at most {} sequences can be aligned with one document",
                    crate::state::MAX_ALIGN_WITH
                ));
                break;
            }
        }
    }

    /// Whether files dropped now go to the "Align with" list rather than opening in tabs: the
    /// panel is shown and they were dropped on it, or they are all traces (which are aligned
    /// far more often than opened).
    fn drops_align_with(&self, ctx: &egui::Context, names: &[String]) -> bool {
        let on_panel = ctx
            .input(|i| i.pointer.latest_pos())
            .zip(self.align_with.panel_rect)
            .is_some_and(|(p, r)| r.contains(p));
        let trace = |n: &String| {
            n.rsplit_once('.').is_some_and(|(_, x)| {
                matches!(x.to_ascii_lowercase().as_str(), "ab1" | "abi" | "scf")
            })
        };
        self.state.align_panel
            && self.state.view == View::Sequence
            && self.state.active_experiment.is_none()
            && self.session.document().is_some()
            && !names.is_empty()
            && names.iter().all(|n| sidebar::is_sequence_file(n))
            && (on_panel || names.iter().all(trace))
    }

    /// Put the selection on the clipboard (`format` as `copy_selection` takes it).
    fn copy_selection(&mut self, ctx: &egui::Context, format: &str) -> bool {
        match self
            .session
            .execute("copy_selection", json!({ "format": format }))
        {
            Ok(v) => {
                let text = v.get("text").and_then(Value::as_str).unwrap_or("");
                let n = v.get("length").and_then(Value::as_u64).unwrap_or(0);
                ctx.copy_text(text.to_string());
                self.state.status = Some(format!("Copied {n} residues"));
                true
            }
            Err(e) => {
                self.status_error(e.message);
                false
            }
        }
    }

    fn empty_state(&mut self, ui: &mut egui::Ui) {
        let t = self.tokens;
        ui.vertical_centered(|ui| {
            ui.add_space(ui.available_height() * 0.3);
            ui.label(RichText::new("No document open").size(20.0));
            ui.add_space(t.spacing);
            ui.horizontal(|ui| {
                let w = 220.0;
                ui.add_space((ui.available_width() - w).max(0.0) / 2.0);
                if ui.button("New sequence").clicked() {
                    self.run("new", json!({}));
                }
                if ui.button("Open…").clicked() {
                    self.open_dialog();
                }
            });
        });
    }

    fn features_table(&mut self, ui: &mut egui::Ui) {
        let t = self.tokens;
        let action = crate::features_view::show(ui, &t, &mut self.session, &mut self.features);
        match action {
            Some(crate::features_view::FeaturesAction::Run(id, params)) => {
                let result = self.run(id, params);
                if let Some(r) = &result
                    && id == "feature_library_add"
                {
                    let name = r.get("name").and_then(Value::as_str).unwrap_or("");
                    self.state.status = Some(format!("Added {name} to your feature library"));
                } else if let Some(r) = &result
                    && id == "feature_library_download"
                {
                    let n = r.get("features").and_then(Value::as_u64).unwrap_or(0);
                    let records = r.get("records").and_then(Value::as_u64).unwrap_or(0);
                    self.state.status = Some(format!(
                        "Vector library: {n} features from {records} GenBank records, for use on this computer only"
                    ));
                }
                if id.starts_with("feature_library") {
                    self.features.library = None;
                    self.features.detected = None;
                }
            }
            Some(crate::features_view::FeaturesAction::PickAnnotations) => {
                #[cfg(not(target_arch = "wasm32"))]
                if let Some(path) = rfd::FileDialog::new()
                    .add_filter("GFF3", &["gff3", "gff"])
                    .pick_file()
                    && let Some(r) = self.run(
                        "import_annotations",
                        json!({ "path": path.to_string_lossy() }),
                    )
                {
                    let added = r.get("added").and_then(Value::as_u64).unwrap_or(0);
                    let skipped = r.get("skipped").and_then(Value::as_u64).unwrap_or(0);
                    let warned = r
                        .get("warnings")
                        .and_then(Value::as_array)
                        .map_or(0, Vec::len);
                    let mut msg = format!("Imported {added} feature(s)");
                    if skipped > 0 {
                        msg.push_str(&format!(", {skipped} already present"));
                    }
                    if warned > 0 {
                        msg.push_str(&format!(", {warned} outside the sequence (see the log)"));
                    }
                    self.state.status = Some(msg);
                }
            }
            Some(crate::features_view::FeaturesAction::PickLibraryFolder) => {
                #[cfg(not(target_arch = "wasm32"))]
                {
                    let mut dialog = rfd::FileDialog::new();
                    if let Some(folder) = &self.state.folder {
                        dialog = dialog.set_directory(folder);
                    }
                    if let Some(folder) = dialog.pick_folder() {
                        let folder = folder.to_string_lossy().into_owned();
                        if let Some(r) =
                            self.run("feature_library_build", json!({ "folder": folder }))
                        {
                            let n = r.get("features").and_then(Value::as_u64).unwrap_or(0);
                            let files = r.get("files").and_then(Value::as_u64).unwrap_or(0);
                            self.state.status =
                                Some(format!("Feature library: {n} features from {files} files"));
                        }
                        self.features.library = None;
                        self.features.detected = None;
                    }
                }
            }
            None => {}
        }
    }

    fn enzymes_view(&mut self, ui: &mut egui::Ui) {
        let t = self.tokens;
        #[cfg(not(target_arch = "wasm32"))]
        {
            self.enzyme_view.downloading = self.download.is_some();
        }
        let action = crate::enzyme_view::show(
            ui,
            &t,
            &self.session,
            &mut self.enzyme_view,
            self.state.enzyme_set.pool.as_deref(),
        );
        match action {
            Some(crate::enzyme_view::EnzymeAction::Select(range)) => {
                self.run("select", json!({ "range": range }));
            }
            Some(crate::enzyme_view::EnzymeAction::UseBuiltin) => {
                self.run("enzymes_reset", json!({}));
            }
            Some(crate::enzyme_view::EnzymeAction::SaveSet(name, enzymes)) => {
                self.save_enzyme_set(name, enzymes);
            }
            Some(crate::enzyme_view::EnzymeAction::ImportSet) => self.pick_enzyme_set_file(),
            Some(crate::enzyme_view::EnzymeAction::ShowOnGel(enzymes)) => {
                let gel = &mut self.state.gel;
                gel.mode = crate::state::GelMode::Digest;
                gel.digest
                    .lanes
                    .push(snappiesnap_engine::gel::LaneSpec::Digest { enzymes });
                self.state.view = View::Gel;
            }
            #[cfg(not(target_arch = "wasm32"))]
            Some(crate::enzyme_view::EnzymeAction::Download) if self.download.is_none() => {
                self.download = Some(crate::enzyme_view::start_download(&self.session));
            }
            _ => {}
        }
    }

    /// Use REBASE once a background download has finished.
    #[cfg(not(target_arch = "wasm32"))]
    fn poll_download(&mut self, ctx: &egui::Context) {
        let Some(rx) = &self.download else {
            return;
        };
        match rx.try_recv() {
            Ok(Ok(path)) => {
                self.download = None;
                if let Some(v) = self.run("enzymes_load", json!({ "path": path })) {
                    let version = v.get("version").and_then(Value::as_str).unwrap_or("?");
                    self.state.status = Some(format!(
                        "REBASE {version} loaded: {} enzymes",
                        v.get("count").and_then(Value::as_u64).unwrap_or(0)
                    ));
                }
            }
            Ok(Err(e)) => {
                self.download = None;
                self.status_error(e);
            }
            Err(std::sync::mpsc::TryRecvError::Empty) => {
                ctx.request_repaint_after(std::time::Duration::from_millis(200));
            }
            Err(std::sync::mpsc::TryRecvError::Disconnected) => {
                self.download = None;
                self.status_error("the download stopped unexpectedly".to_string());
            }
        }
    }

    fn align_view(&mut self, ui: &mut egui::Ui) {
        use crate::align_view::AlignAction;
        let t = self.tokens;
        let files: Vec<String> = self
            .state
            .folder
            .as_deref()
            .and_then(|f| self.session.platform().list_dir(f).ok())
            .unwrap_or_default()
            .into_iter()
            .filter(|e| !e.is_dir && crate::sidebar::is_sequence_file(&e.name))
            .map(|e| e.path)
            .collect();
        let action = egui::ScrollArea::vertical()
            .id_salt("align-view")
            .show(ui, |ui| {
                crate::align_view::show(ui, &t, &self.session, &mut self.align, &files)
            })
            .inner;
        match action {
            Some(AlignAction::Align(params)) => {
                let answer = self
                    .session
                    .execute("align_to_reference", params)
                    .map_err(|e| e.message);
                self.align.result = Some(answer);
            }
            Some(AlignAction::AlignMany(params)) => {
                let answer = self
                    .session
                    .execute("align_sequences", params)
                    .map_err(|e| e.message);
                self.align.msa.result = Some(answer);
            }
            Some(AlignAction::Show { start, end }) => {
                if self
                    .run("select", json!({ "range": { "start": start, "end": end } }))
                    .is_some()
                {
                    self.state.view = View::Sequence;
                }
            }
            Some(AlignAction::AddFiles { multiple }) => self.pick_align_files(if multiple {
                AlignFor::Multiple
            } else {
                AlignFor::Pairwise
            }),
            Some(AlignAction::AlignWith) => {
                self.state.align_panel = true;
                self.state.view = View::Sequence;
            }
            Some(AlignAction::Annotate(params)) => {
                self.run("add_feature", params);
            }
            None => {}
        }
    }

    fn primers_table(&mut self, ui: &mut egui::Ui) {
        let t = self.tokens;
        use crate::primer_design::DesignAction;
        let pcr = crate::primer_design::pcr_row(ui, &t, &self.session, &mut self.primer_design);
        if pcr.is_some() {
            ui.separator();
        }
        let design = crate::primer_design::show(ui, &t, &self.session, &mut self.primer_design);
        match pcr.or(design) {
            Some(DesignAction::Design(params)) => {
                let region = params
                    .get("range")
                    .and_then(|r| serde_json::from_value(r.clone()).ok())
                    .unwrap_or_default();
                let answer = self
                    .session
                    .execute("design_primers", params)
                    .map_err(|e| e.message);
                self.primer_design.set_results(region, answer);
            }
            Some(DesignAction::Add(pair)) => {
                for (name, residues) in pair {
                    self.run("add_primer", json!({ "name": name, "residues": residues }));
                }
            }
            Some(DesignAction::Select(range)) => {
                self.run("select", json!({ "range": range }));
            }
            Some(DesignAction::PcrGel(forward, reverse)) => {
                let gel = &mut self.state.gel;
                gel.mode = crate::state::GelMode::Pcr;
                gel.pcr
                    .lanes
                    .push(snappiesnap_engine::gel::LaneSpec::Pcr { forward, reverse });
                self.state.view = View::Gel;
            }
            Some(DesignAction::PcrOpen(forward, reverse)) => {
                // The product opens in a new tab, which starts on its map.
                self.run(
                    "pcr_open",
                    json!({ "forward": forward, "reverse": reverse }),
                );
            }
            None => {}
        }
        ui.separator();
        let mut oligo_action = None;
        egui::CollapsingHeader::new(RichText::new("Lab oligos").strong())
            .id_salt("lab-oligos")
            .default_open(true)
            .show(ui, |ui| {
                oligo_action =
                    crate::oligo_view::show(ui, &t, &mut self.session, &mut self.oligo_view);
            });
        if let Some(a) = oligo_action {
            self.oligo_action(a);
        }
        ui.separator();
        let query = self.enzyme_query();
        let marks = self.site_cache.get(&self.session, &query).primers;
        let Some(d) = self.session.document() else {
            return;
        };
        if d.molecule.primers().is_empty() {
            ui.label(RichText::new("No primers yet").color(t.text_muted));
            return;
        }
        let len = d.molecule.len();
        let mut select = None;
        let mut edit = None;
        let mut delete = None;
        let mut apply = None;
        if let Some((id, name, seq, desc)) = &mut self.primer_edit {
            egui::Frame::new()
                .stroke(Stroke::new(1.0, t.border))
                .corner_radius(t.radius)
                .inner_margin(8)
                .show(ui, |ui| {
                    egui::Grid::new("primer-edit")
                        .num_columns(2)
                        .spacing(vec2(t.spacing * 2.0, t.spacing))
                        .show(ui, |ui| {
                            ui.label("Name");
                            ui.add(egui::TextEdit::singleline(name).desired_width(240.0));
                            ui.end_row();
                            ui.label("Sequence (5' to 3')");
                            ui.add(
                                egui::TextEdit::singleline(seq)
                                    .font(egui::TextStyle::Monospace)
                                    .desired_width(480.0),
                            );
                            ui.end_row();
                            ui.label("Description");
                            ui.add(egui::TextEdit::singleline(desc).desired_width(480.0));
                            ui.end_row();
                        });
                    ui.horizontal(|ui| {
                        if ui.button("Save").clicked() {
                            apply = Some(json!({
                                "id": *id, "name": name.clone(), "residues": seq.clone(), "description": desc.clone(),
                            }));
                        }
                        if ui.button("Cancel").clicked() {
                            edit = Some(None);
                        }
                    });
                });
            ui.add_space(t.spacing);
        }
        egui::ScrollArea::both().show(ui, |ui| {
            egui::Grid::new("primers")
                .striped(true)
                .spacing(vec2(t.spacing * 3.0, t.spacing))
                .show(ui, |ui| {
                    for h in [
                        "Name",
                        "Sequence (5' to 3')",
                        "Length",
                        "Tm",
                        "Structure",
                        "Binds at",
                        "",
                    ] {
                        ui.label(RichText::new(h).strong());
                    }
                    ui.end_row();
                    for (i, p) in d.molecule.primers().iter().enumerate() {
                        ui.label(&p.name);
                        ui.monospace(p.seq.as_str());
                        ui.label(format!("{} nt", p.seq.len()));
                        let full = snappiesnap_engine::bio::tm::tm(
                            p.seq.as_bytes(),
                            &snappiesnap_engine::bio::tm::TmConditions::default(),
                        );
                        ui.label(full.map_or("-".to_string(), |t| format!("{t:.1} \u{b0}C")))
                            .on_hover_text(
                                "Whole primer with its exact complement (100 nM, 100 mM Na+)",
                            );
                        let (hp, dimer, checks) = crate::enzymes::primer_structure(
                            p.seq.as_bytes(),
                            &snappiesnap_engine::bio::tm::TmConditions::default(),
                        );
                        let fmt = crate::enzymes::format_dg;
                        let detail = format!(
                            "Hairpin \u{394}G {}, self-dimer \u{394}G {} kcal/mol (37 \u{b0}C){}",
                            fmt(hp),
                            fmt(dimer),
                            checks
                                .iter()
                                .map(|w| format!("\n! {w}"))
                                .collect::<String>()
                        );
                        let summary = if checks.is_empty() {
                            RichText::new("OK").color(t.text_muted)
                        } else {
                            RichText::new(format!("! {}", checks.join("; "))).color(t.error)
                        };
                        ui.add(egui::Label::new(summary).truncate())
                            .on_hover_text(detail);
                        ui.horizontal_wrapped(|ui| {
                            ui.set_min_width(260.0);
                            let mut any = false;
                            for m in marks.iter().filter(|m| m.primer == i) {
                                any = true;
                                let end = if m.range.end == 0 { len } else { m.range.end };
                                let arrow = if m.forward { "fwd" } else { "rev" };
                                let text = format!("{arrow} {}..{}", m.range.start + 1, end);
                                let rich = if m.primary {
                                    RichText::new(text).color(t.text)
                                } else {
                                    RichText::new(format!("{text}*")).color(t.text_muted)
                                };
                                let r = ui.add(egui::Label::new(rich).sense(Sense::click()));
                                if r.on_hover_text(m.detail(len)).clicked() {
                                    select = Some(m.range);
                                }
                            }
                            if !any {
                                ui.label(RichText::new("doesn't bind").color(t.text_muted));
                            }
                        });
                        ui.horizontal(|ui| {
                            if ui.small_button("Edit").clicked() {
                                edit = Some(Some((
                                    p.id.0,
                                    p.name.clone(),
                                    p.seq.as_str().to_string(),
                                    crate::notes_view::plain(&p.description),
                                )));
                            }
                            if ui.small_button("Delete").clicked() {
                                delete = Some(p.id.0);
                            }
                        });
                        ui.end_row();
                    }
                });
            ui.label(
                RichText::new(
                    "* secondary site (shorter or with a mismatch). Click a site to select it.",
                )
                .small()
                .color(t.text_muted),
            );
        });
        if let Some(range) = select {
            self.run("select", json!({ "range": range }));
        }
        if let Some(e) = edit {
            self.primer_edit = e;
        }
        if let Some(id) = delete {
            self.run("remove_primer", json!({ "id": id }));
            self.primer_edit = None;
        }
        if let Some(mut params) = apply {
            // An untouched description stays as the file has it (HTML).
            if let (Some(o), Some(d)) = (params.as_object_mut(), self.session.document()) {
                let id = o.get("id").and_then(Value::as_u64);
                let old = d
                    .molecule
                    .primers()
                    .iter()
                    .find(|p| Some(p.id.0) == id)
                    .map(|p| crate::notes_view::plain(&p.description));
                if old.as_deref() == o.get("description").and_then(Value::as_str) {
                    o.remove("description");
                }
            }
            if self.run("edit_primer", params).is_some() {
                self.primer_edit = None;
            }
        }
    }

    /// Ask where to save a copy without the construction history, and save it.
    fn save_copy_without_history(&mut self) {
        #[cfg(not(target_arch = "wasm32"))]
        {
            let name = self
                .session
                .document()
                .map(|d| format!("{} (no history).dna", d.molecule.name()))
                .unwrap_or_else(|| "Untitled.dna".into());
            if let Some(path) = rfd::FileDialog::new()
                .add_filter("SnapGene", &["dna"])
                .set_file_name(name)
                .save_file()
                && let Some(r) = self.run(
                    "save",
                    json!({ "path": path.to_string_lossy(), "without_history": true }),
                )
            {
                let bytes = r.get("bytes").and_then(Value::as_u64).unwrap_or(0);
                self.state.status = Some(format!(
                    "Saved a copy without history ({} kB)",
                    bytes.div_ceil(1000)
                ));
            }
        }
        #[cfg(target_arch = "wasm32")]
        self.web_save_prompt(crate::web_files::SaveTarget::DocumentWithoutHistory);
    }

    fn notes_tab(&mut self, ui: &mut egui::Ui) {
        let t = self.tokens;
        let save = egui::ScrollArea::vertical()
            .id_salt("notes")
            .show(ui, |ui| {
                crate::notes_view::show(ui, &t, &self.session, &mut self.notes_edit)
            })
            .inner;
        if let Some(params) = save {
            self.run("set_notes", params);
        }
    }

    fn history(&mut self, ui: &mut egui::Ui) {
        use crate::state::HistoryLayout;
        let t = self.tokens;
        let mut open = None;
        let has_history = self
            .session
            .document()
            .is_some_and(|d| d.molecule.history().is_some());
        let mut lean = false;
        let mut remove = false;
        if has_history {
            ui.horizontal(|ui| {
                ui.label(RichText::new("Show as").color(t.text_muted));
                ui.selectable_value(&mut self.state.history, HistoryLayout::List, "List");
                ui.selectable_value(&mut self.state.history, HistoryLayout::Diagram, "Diagram");
                ui.separator();
                lean = ui
                    .button("Save a copy without history\u{2026}")
                    .on_hover_text("A .dna file without the history, which holds every ancestor's sequence. This document keeps its history and its own file.")
                    .clicked();
                remove = ui
                    .button("Remove history")
                    .on_hover_text("Take the history out of this document (undoable); the next save writes the file without it")
                    .clicked();
            });
            ui.add_space(t.spacing);
        }
        if remove {
            self.run("remove_history", json!({}));
        }
        if lean {
            self.save_copy_without_history();
        }
        if has_history && self.state.history == HistoryLayout::Diagram {
            egui::ScrollArea::both()
                .id_salt("history-diagram")
                .show(ui, |ui| {
                    ui.label(RichText::new("How this construct was made").strong());
                    ui.label(
                        RichText::new("The construct is at the top, with the molecules it was made from below it. Open shows a step as a molecule of its own.")
                            .color(t.text_muted),
                    );
                    ui.add_space(t.spacing);
                    let doc_id = self.session.active_document_id();
                    if let Some(d) = self.session.document()
                        && let Some(h) = d.molecule.history()
                        && let Some(crate::history_diagram::Action::Open(id)) =
                            crate::history_diagram::show(ui, &t, doc_id, &d.molecule, h)
                    {
                        open = Some(id);
                    }
                    ui.add_space(t.spacing * 3.0);
                    self.session_edits(ui);
                });
        } else {
            egui::ScrollArea::vertical()
                .id_salt("history")
                .show(ui, |ui| {
                    open =
                        crate::history_view::show(ui, &t, &self.session, &mut self.history_cache);
                    ui.add_space(t.spacing * 3.0);
                    self.session_edits(ui);
                });
        }
        if let Some(id) = open {
            self.run("history_open", json!({ "id": id }));
        }
    }

    /// Edits made in this session (undo history) and recent problems.
    fn session_edits(&mut self, ui: &mut egui::Ui) {
        let t = self.tokens;
        let h = self.session.history();
        let mut steps: Option<(&str, usize)> = None;
        ui.label(RichText::new("Edits in this session").strong());
        ui.label(
            RichText::new("Click an entry to go back (or forward) to the state after it.")
                .color(t.text_muted),
        );
        ui.add_space(t.spacing);
        let row = |ui: &mut egui::Ui, text: String, current: bool, undone: bool| {
            let mut rt = RichText::new(text);
            if current {
                rt = rt.strong().color(t.accent);
            } else if undone {
                rt = rt.color(t.text_muted).italics();
            }
            ui.add(egui::Label::new(rt).sense(Sense::click())).clicked()
        };
        if row(ui, "Opened / created".into(), h.undo.is_empty(), false) {
            steps = Some(("undo", h.undo.len()));
        }
        for (i, label) in h.undo.iter().enumerate() {
            let current = i + 1 == h.undo.len();
            if row(ui, format!("{}. {label}", i + 1), current, false) {
                steps = Some(("undo", h.undo.len() - (i + 1)));
            }
        }
        for (j, label) in h.redo.iter().enumerate() {
            if row(
                ui,
                format!("{}. {label} (undone)", h.undo.len() + j + 1),
                false,
                true,
            ) {
                steps = Some(("redo", j + 1));
            }
        }
        ui.add_space(t.spacing * 2.0);
        ui.label(RichText::new("Recent problems").strong());
        for e in self.session.log_tail(20).iter().rev() {
            ui.label(RichText::new(format!("{}: {}", e.source, e.message)).color(t.text_muted));
        }
        if let Some((id, n)) = steps {
            for _ in 0..n {
                if self.run(id, json!({})).is_none() {
                    break;
                }
            }
        }
    }

    /// Web build: ask for a file name to download `target` under, named after what is saved.
    #[cfg(target_arch = "wasm32")]
    fn web_save_prompt(&mut self, target: crate::web_files::SaveTarget) {
        use crate::web_files::{SavePrompt, SaveTarget};
        let doc = self
            .session
            .document()
            .map(|d| d.molecule.name().to_string());
        let name = match target {
            SaveTarget::Document => doc.unwrap_or_else(|| "Untitled".into()),
            SaveTarget::DocumentWithoutHistory => {
                format!("{} (no history)", doc.unwrap_or_else(|| "Untitled".into()))
            }
            SaveTarget::Map => format!("{} map", doc.unwrap_or_else(|| "Untitled".into())),
            SaveTarget::Experiment(i) => self
                .state
                .experiments
                .get(i)
                .map(|e| e.name.trim().to_string())
                .unwrap_or_else(|| "Experiment".into()),
        };
        self.web_save = Some(SavePrompt::new(target, &name));
    }

    /// Web build: the Download prompt, and the download once a name is chosen.
    #[cfg(target_arch = "wasm32")]
    fn web_save_dialog(&mut self, ctx: &egui::Context) {
        use crate::web_files::SaveTarget;
        let Some(prompt) = self.web_save.as_mut() else {
            return;
        };
        match crate::web_files::prompt(ctx, &self.tokens, prompt) {
            None => {}
            Some(false) => self.web_save = None,
            Some(true) => {
                let target = prompt.target;
                let name = target.file_name(&prompt.name);
                self.web_save = None;
                match target {
                    SaveTarget::Document => {
                        if self.run("save", json!({ "path": name })).is_some() {
                            self.download_saved(&name);
                        }
                    }
                    SaveTarget::DocumentWithoutHistory => {
                        if self
                            .run("save", json!({ "path": name, "without_history": true }))
                            .is_some()
                        {
                            self.download_saved(&name);
                        }
                    }
                    SaveTarget::Experiment(i) => self.save_experiment(i, &name),
                    SaveTarget::Map => {
                        let (w, h) = self.map_export_size();
                        let svg = name.to_ascii_lowercase().ends_with(".svg");
                        let result = self
                            .map_file(svg, w, h)
                            .and_then(|bytes| crate::web::download(&name, &bytes));
                        match result {
                            Ok(()) => self.state.status = Some(format!("Downloaded {name}")),
                            Err(e) => self.status_error(format!("the map couldn't be saved: {e}")),
                        }
                    }
                }
            }
        }
    }

    fn confirm_dialog(&mut self, ctx: &egui::Context) {
        let Some(c) = &self.confirm else {
            return;
        };
        let question = c.question.clone();
        let mut answer = None;
        egui::Window::new("Unsaved changes")
            .collapsible(false)
            .resizable(false)
            .anchor(egui::Align2::CENTER_CENTER, vec2(0.0, 0.0))
            .show(ctx, |ui| {
                ui.label(question);
                ui.horizontal(|ui| {
                    if ui.button("Discard").clicked() {
                        answer = Some(true);
                    }
                    if ui.button("Cancel").clicked() {
                        answer = Some(false);
                    }
                });
            });
        match answer {
            Some(true) => {
                if let Some(c) = self.confirm.take() {
                    self.run(&c.id, c.params);
                }
            }
            Some(false) => self.confirm = None,
            None => {}
        }
    }
}

/// Which gel to save as an image.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum GelSource {
    /// The open document's gel in the Gel tab.
    Gel,
    /// An experiment tab's gel.
    Experiment(usize),
}

/// Accept `null` for "no changes".
#[cfg(not(target_arch = "wasm32"))]
fn normalise(params: Value) -> Value {
    if params.is_null() { json!({}) } else { params }
}

/// A toolbar button: icon above a small label. Returns true when clicked.
fn tool_button(ui: &mut egui::Ui, t: &Tokens, icon: Icon, label: &str, enabled: bool) -> bool {
    tool_button_response(ui, t, icon, label, enabled).clicked()
}

/// [`tool_button`], returning its response (for a button that opens a menu).
fn tool_button_response(
    ui: &mut egui::Ui,
    t: &Tokens,
    icon: Icon,
    label: &str,
    enabled: bool,
) -> egui::Response {
    // Wide enough for the label (about 6 px per character at 11 px).
    let size = vec2(
        (label.chars().count() as f32 * 6.0 + 8.0).max(58.0),
        t.toolbar_icon + 22.0,
    );
    let sense = if enabled {
        Sense::click()
    } else {
        Sense::hover()
    };
    let (rect, response) = ui.allocate_exact_size(size, sense);
    let color = if enabled {
        t.text
    } else {
        t.text_muted.gamma_multiply(0.5)
    };
    if enabled && response.hovered() {
        ui.painter()
            .rect_filled(rect, CornerRadius::same(t.radius), t.border);
    }
    let icon_rect = egui::Rect::from_center_size(
        rect.center_top() + vec2(0.0, 4.0 + t.toolbar_icon / 2.0),
        vec2(t.toolbar_icon, t.toolbar_icon),
    );
    icons::paint(ui.painter(), icon_rect, icon, color);
    ui.painter().text(
        rect.center_bottom() - vec2(0.0, 4.0),
        egui::Align2::CENTER_BOTTOM,
        label,
        egui::FontId::proportional(11.0),
        color,
    );
    response
}

/// A square on/off button with an icon (map layers); returns whether it was clicked.
/// Draw `draw` into `area` at `zoom` and `pan`: the map gets a larger rectangle (in both
/// directions when circular, only along the sequence when linear) clipped to `area`. Pinch or
/// Cmd/Ctrl+scroll over the area zooms at the pointer; dragging moves the map.
pub(crate) fn zoomed<R>(
    ui: &mut egui::Ui,
    id: &str,
    area: Rect,
    zoom: &mut f32,
    pan: &mut Vec2,
    circular: bool,
    draw: impl FnOnce(&mut egui::Ui) -> R,
) -> R {
    let response = ui.interact(area, ui.make_persistent_id(id), Sense::drag());
    let old = (*zoom).clamp(1.0, crate::state::MAX_ZOOM);
    if response.hovered() {
        let factor = ui.input(|i| i.zoom_delta());
        if (factor - 1.0).abs() > f32::EPSILON
            && let Some(p) = ui.input(|i| i.pointer.hover_pos())
        {
            let new = (old * factor).clamp(1.0, crate::state::MAX_ZOOM);
            // Keep the point under the pointer where it is.
            let center = area.center() + *pan;
            let rel = (p - center) / old;
            let moved = p - rel * new - area.center();
            *pan = if circular { moved } else { vec2(moved.x, 0.0) };
            *zoom = new;
        }
    }
    if response.dragged() && *zoom > 1.0 {
        *pan += response.drag_delta();
    }
    let z = (*zoom).clamp(1.0, crate::state::MAX_ZOOM);
    let size = if circular {
        area.size() * z
    } else {
        vec2(area.width() * z, area.height())
    };
    // The map always covers the area.
    let room = (size - area.size()) / 2.0;
    pan.x = pan.x.clamp(-room.x, room.x);
    pan.y = if circular {
        pan.y.clamp(-room.y, room.y)
    } else {
        0.0
    };
    let rect = Rect::from_center_size(area.center() + *pan, size);
    let mut child = ui.new_child(egui::UiBuilder::new().max_rect(rect).id_salt(id));
    child.set_clip_rect(area.intersect(ui.clip_rect()));
    draw(&mut child)
}

/// The small arrow next to a layer button that opens its menu.
fn chevron(ui: &mut egui::Ui, t: &Tokens) -> egui::Response {
    let (rect, response) = ui.allocate_exact_size(vec2(14.0, 30.0), Sense::click());
    if response.hovered() {
        ui.painter()
            .rect_filled(rect, CornerRadius::same(t.radius), t.border);
    }
    let c = rect.center();
    ui.painter().add(egui::Shape::line(
        vec![
            egui::pos2(c.x - 3.5, c.y - 1.5),
            egui::pos2(c.x, c.y + 2.0),
            egui::pos2(c.x + 3.5, c.y - 1.5),
        ],
        Stroke::new(1.4, t.text_muted),
    ));
    response
}

/// What the enzyme menu asked for, besides changing the settings.
enum EnzymeMenuAction {
    /// Import an enzyme set from a file.
    #[cfg_attr(target_arch = "wasm32", allow(dead_code))]
    Import,
    /// Delete the user's enzyme set of this name.
    Remove(String),
}

/// The enzyme menu: which enzymes to draw from (the whole database or one of the user's enzyme
/// sets), how often they may cut, their shortest site, whether they must be sold, or a list of
/// chosen enzymes (which overrides the rest). `unknown` lists chosen names the database
/// doesn't have.
fn enzyme_menu(
    ui: &mut egui::Ui,
    t: &Tokens,
    state: &mut UiState,
    db: &snappiesnap_engine::EnzymeDb,
    unknown: &[String],
) -> Option<EnzymeMenuAction> {
    use crate::state::EnzymeSet;
    ui.set_min_width(260.0);
    ui.spacing_mut().item_spacing.y = 0.0;
    let chosen = !state.enzyme_set.only.is_empty();
    let set = &mut state.enzyme_set;
    let mut action = None;
    ui.label(RichText::new("Enzymes from").strong());
    ui.add_space(3.0);
    ui.add_enabled_ui(!chosen, |ui| {
        let database = match &db.version {
            Some(v) => format!("{} {v}", db.source),
            None => db.source.clone(),
        };
        let all = format!("All enzymes ({database}, {})", db.enzymes.len());
        if check_row(ui, t, &all, !chosen && set.pool.is_none()) {
            set.pool = None;
        }
        for user in &db.sets {
            let missing = user.missing(db);
            let label = if missing.is_empty() {
                format!("{} ({})", user.name, user.enzymes.len())
            } else {
                format!(
                    "{} ({} of {} known)",
                    user.name,
                    user.enzymes.len().saturating_sub(missing.len()),
                    user.enzymes.len()
                )
            };
            let on = !chosen
                && set
                    .pool
                    .as_ref()
                    .is_some_and(|p| p.eq_ignore_ascii_case(&user.name));
            if check_row(ui, t, &label, on) {
                set.pool = Some(user.name.clone());
            }
        }
    });
    // The chosen set's enzymes that the database lacks, and how to get them.
    if !chosen && let Some(user) = set.pool.as_ref().and_then(|p| db.set(p)) {
        let missing = user.missing(db);
        if !missing.is_empty() {
            let mut names = missing
                .iter()
                .take(6)
                .copied()
                .collect::<Vec<_>>()
                .join(", ");
            if missing.len() > 6 {
                names.push_str(", \u{2026}");
            }
            let fix = if db.is_builtin() {
                "Download REBASE in the Enzymes tab to include them."
            } else {
                "Check their spelling in the set file."
            };
            ui.add_space(3.0);
            ui.add(
                egui::Label::new(
                    RichText::new(format!("Not in the enzyme database: {names}. {fix}"))
                        .color(t.error),
                )
                .wrap(),
            );
        }
    }
    ui.add_space(4.0);
    ui.horizontal(|ui| {
        #[cfg(not(target_arch = "wasm32"))]
        if ui
            .button("Import enzyme set\u{2026}")
            .on_hover_text(
                "A text file of enzyme names, such as a set exported from another program",
            )
            .clicked()
        {
            action = Some(EnzymeMenuAction::Import);
            ui.close();
        }
        if let Some(name) = set.pool.clone()
            && db.set(&name).is_some()
            && ui
                .button("Remove set")
                .on_hover_text(format!("Delete the enzyme set \u{201c}{name}\u{201d}"))
                .clicked()
        {
            action = Some(EnzymeMenuAction::Remove(name));
            ui.close();
        }
    });
    ui.add_space(3.0);
    ui.separator();
    ui.add_space(3.0);
    ui.add_enabled_ui(!chosen, |ui| {
        for (cuts, label) in [
            (1, "Unique cutters"),
            (2, "Cutting once or twice"),
            (3, "Cutting up to 3 times"),
        ] {
            if check_row(ui, t, label, !chosen && set.max_cuts == cuts) {
                set.max_cuts = cuts;
            }
        }
        ui.add_space(3.0);
        ui.separator();
        ui.add_space(3.0);
        for (length, label) in [
            (0, "Any site length"),
            (4, "Sites of 4+ bases"),
            (5, "Sites of 5+ bases"),
            (6, "Sites of 6+ bases"),
            (7, "Sites of 7+ bases"),
            (8, "Sites of 8+ bases"),
        ] {
            if check_row(ui, t, label, !chosen && set.min_length == length) {
                set.min_length = length;
            }
        }
        ui.add_space(3.0);
        ui.separator();
        ui.add_space(3.0);
        // A set says what is at hand, so whether an enzyme is sold doesn't matter then.
        ui.add_enabled_ui(set.pool.is_none(), |ui| {
            if check_row(
                ui,
                t,
                "Commercially available only",
                !chosen && set.pool.is_none() && set.commercial_only,
            ) {
                set.commercial_only = !set.commercial_only;
            }
        });
    });
    ui.add_space(3.0);
    ui.separator();
    ui.spacing_mut().item_spacing.y = 4.0;
    ui.label(RichText::new("Only these enzymes").strong());
    // The text being typed is kept apart from the list, so commas and spaces survive.
    let id = ui.make_persistent_id("enzyme-only-text");
    let mut text: String = ui
        .data_mut(|d| d.get_temp::<String>(id))
        .unwrap_or_else(|| set.only.join(", "));
    let edit = ui.add(
        egui::TextEdit::singleline(&mut text)
            .hint_text("EcoRI, BamHI, \u{2026}")
            .desired_width(240.0),
    );
    if edit.changed() {
        let only: Vec<String> = text
            .split([',', ' ', ';'])
            .map(str::trim)
            .filter(|n| !n.is_empty())
            .map(str::to_string)
            .collect();
        *set = EnzymeSet {
            only,
            ..set.clone()
        }
        .clamped();
    }
    ui.data_mut(|d| d.insert_temp(id, text.clone()));
    if !unknown.is_empty() {
        ui.label(
            RichText::new(format!("Not in the enzyme list: {}", unknown.join(", "))).color(t.error),
        );
    }
    ui.horizontal(|ui| {
        if chosen && ui.button("Show the filters again").clicked() {
            set.only.clear();
            ui.data_mut(|d| d.insert_temp(id, String::new()));
        }
        if ui
            .button("Default")
            .on_hover_text("Unique 6+ cutters (sold ones, unless an enzyme set is chosen)")
            .clicked()
        {
            *set = EnzymeSet {
                pool: set.pool.clone(),
                ..EnzymeSet::default()
            };
            ui.data_mut(|d| d.insert_temp(id, String::new()));
        }
    });
    action
}

/// What the translation menu asked for, besides changing the settings.
enum TranslationAction {
    AddFeature,
}

/// One row of a pick-one menu: a check mark when chosen.
fn check_row(ui: &mut egui::Ui, t: &Tokens, label: &str, on: bool) -> bool {
    let h = 24.0;
    let w = ui.available_width().max(220.0);
    let (rect, response) = ui.allocate_exact_size(vec2(w, h), Sense::click());
    if response.hovered() {
        ui.painter()
            .rect_filled(rect, CornerRadius::same(t.radius), t.selection);
    }
    if on {
        let c = egui::pos2(rect.left() + 12.0, rect.center().y);
        ui.painter().add(egui::Shape::line(
            vec![
                egui::pos2(c.x - 4.0, c.y),
                egui::pos2(c.x - 1.0, c.y + 3.5),
                egui::pos2(c.x + 5.0, c.y - 4.0),
            ],
            Stroke::new(1.8, t.accent),
        ));
    }
    ui.painter().text(
        egui::pos2(rect.left() + 26.0, rect.center().y),
        egui::Align2::LEFT_CENTER,
        label,
        egui::FontId::proportional(13.0),
        t.text,
    );
    response.clicked()
}

/// The ORF/translation menu: one choice from the list (ORFs only, one frame, three frames, all
/// six), then adding a translated feature and the genetic code.
fn translation_menu(
    ui: &mut egui::Ui,
    t: &Tokens,
    state: &mut UiState,
    has_selection: bool,
) -> Option<TranslationAction> {
    let mut action = None;
    ui.set_min_width(240.0);
    ui.spacing_mut().item_spacing.y = 0.0;
    for item in Translation::MENU {
        match item {
            Some((value, label)) => {
                if check_row(ui, t, label, state.translation == value) {
                    state.translation = value;
                    state.orfs = true;
                    ui.close();
                }
            }
            None => {
                ui.add_space(3.0);
                ui.separator();
                ui.add_space(3.0);
            }
        }
    }
    ui.spacing_mut().item_spacing.y = 4.0;
    if ui
        .add_enabled(
            has_selection,
            egui::Button::new("Add translated feature").frame(false),
        )
        .on_hover_text("A CDS over the selection, on the strand of the chosen frames")
        .on_disabled_hover_text("Select a stretch of the sequence first")
        .clicked()
    {
        action = Some(TranslationAction::AddFeature);
        ui.close();
    }
    let current = snappiesnap_engine::bio::translate::code(state.genetic_code)
        .map_or("Standard".to_string(), |c| format!("{} {}", c.id, c.name));
    ui.menu_button(format!("Genetic code: {current}"), |ui| {
        egui::ScrollArea::vertical()
            .max_height(360.0)
            .show(ui, |ui| {
                for c in snappiesnap_engine::bio::codes::CODES {
                    if check_row(
                        ui,
                        t,
                        &format!("{} {}", c.id, c.name),
                        state.genetic_code == c.id,
                    ) {
                        state.genetic_code = c.id;
                        ui.close();
                    }
                }
            });
    });
    action
}

pub(crate) fn layer_toggle(
    ui: &mut egui::Ui,
    t: &Tokens,
    icon: Icon,
    tip: &str,
    on: bool,
) -> egui::Response {
    let (rect, response) = ui.allocate_exact_size(vec2(30.0, 30.0), Sense::click());
    let fill = if on {
        Some(t.accent)
    } else if response.hovered() {
        Some(t.border)
    } else {
        None
    };
    if let Some(fill) = fill {
        ui.painter()
            .rect_filled(rect, CornerRadius::same(t.radius), fill);
    }
    let color = if on { t.accent_text } else { t.text_muted };
    icons::paint(ui.painter(), rect.shrink(7.0), icon, color);
    let state = if on { "shown" } else { "hidden" };
    response.on_hover_text(format!("{tip} ({state})"))
}

fn separator(ui: &mut egui::Ui, t: &Tokens) {
    let (rect, _) = ui.allocate_exact_size(vec2(9.0, t.toolbar_icon + 18.0), Sense::hover());
    ui.painter().line_segment(
        [rect.center_top(), rect.center_bottom()],
        Stroke::new(1.0, t.border),
    );
}

impl eframe::App for App {
    fn logic(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        #[cfg(not(target_arch = "wasm32"))]
        self.handle_control(ctx);
        #[cfg(not(target_arch = "wasm32"))]
        self.write_recording();
        #[cfg(not(target_arch = "wasm32"))]
        self.save_settings();
        #[cfg(not(target_arch = "wasm32"))]
        self.poll_download(ctx);
        self.handle_dropped_files(ctx);
        #[cfg(target_arch = "wasm32")]
        {
            self.take_web_files();
            self.take_web_enzyme_sets();
            self.take_web_align_files();
            self.take_web_reads();
            let unsaved = self.session.documents().iter().any(|d| d.dirty)
                || self.state.experiments.iter().any(|e| e.unsaved());
            self.web_unsaved
                .store(unsaved, std::sync::atomic::Ordering::Relaxed);
        }
        self.poll_update_check(ctx);
        self.sync_scale(ctx);
        #[cfg(not(target_arch = "wasm32"))]
        self.track_window(ctx);
        let system_dark = ctx.system_theme().is_none_or(|t| t == Theme::Dark);
        self.tokens = Tokens::for_choice(self.state.theme, system_dark);
        if self.applied != Some(self.tokens) {
            self.tokens.apply(ctx);
            self.applied = Some(self.tokens);
        }
        let title = self.window_title();
        if title != self.title {
            ctx.send_viewport_cmd(egui::ViewportCommand::Title(title.clone()));
            #[cfg(target_arch = "wasm32")]
            crate::web::set_title(&title);
            self.title = title;
        }
    }

    fn ui(&mut self, ui: &mut egui::Ui, _frame: &mut eframe::Frame) {
        self.draw(ui);
    }
}

impl App {
    /// Draw the whole window. Used by the native window and by offscreen rendering, so both
    /// show exactly the same thing.
    pub fn draw(&mut self, ui: &mut egui::Ui) {
        self.sync_tab();
        let ctx = ui.ctx().clone();
        self.shortcuts(&ctx);
        let t = self.tokens;
        let chrome = egui::Frame::new()
            .fill(t.chrome)
            .inner_margin(egui::Margin::symmetric(8, 4));
        // Document tabs on top: the toolbar and views below act on the active one.
        egui::Panel::top("documents")
            .frame(egui::Frame::new().fill(t.panel).inner_margin(egui::Margin {
                left: 8,
                right: 8,
                top: 6,
                bottom: 0,
            }))
            .show(ui, |ui| self.document_tabs(ui));
        egui::Panel::top("toolbar")
            .frame(chrome)
            .show(ui, |ui| self.toolbar(ui));
        self.refresh_info();
        egui::Panel::bottom("status")
            .frame(chrome)
            .show(ui, |ui| self.status_bar(ui));
        if self.state.sidebar && cfg!(not(target_arch = "wasm32")) {
            egui::Panel::left("files")
                .resizable(true)
                .default_size(250.0)
                .min_size(160.0)
                .frame(
                    egui::Frame::new()
                        .fill(t.chrome)
                        .inner_margin(egui::Margin::same(10)),
                )
                .show(ui, |ui| self.sidebar_panel(ui));
        }
        // The views (Map, Sequence, ...) belong to a construct, not to an experiment.
        if self.state.active_experiment.is_none() {
            egui::Panel::top("tabs")
                .frame(
                    egui::Frame::new()
                        .fill(t.panel)
                        .inner_margin(egui::Margin::symmetric(8, 6)),
                )
                .show(ui, |ui| self.tabs(ui));
        }
        // The Map and Sequence views take the background chosen in Preferences.
        let canvas = if self.state.active_experiment.is_none()
            && matches!(self.state.view, View::Map | View::Sequence)
        {
            self.canvas_tokens().canvas
        } else {
            t.canvas
        };
        egui::CentralPanel::default()
            .frame(
                egui::Frame::new()
                    .fill(canvas)
                    .inner_margin(egui::Margin::same(12)),
            )
            .show(ui, |ui| {
                // Widget state (scroll positions and the like) is kept per tab.
                let tab = self.current_tab.unwrap_or(0);
                let experiment = self.state.active_experiment;
                ui.push_id(("tab", tab, experiment), |ui| self.central(ui));
            });
        self.confirm_dialog(&ctx);
        #[cfg(target_arch = "wasm32")]
        self.web_save_dialog(&ctx);
        self.preferences(&ctx);
        #[cfg(not(target_arch = "wasm32"))]
        self.gel_export_dialog(&ctx);
        if ctx.input(|i| !i.raw.hovered_files.is_empty()) {
            let screen = ctx.content_rect();
            let painter = ctx.layer_painter(egui::LayerId::new(
                egui::Order::Foreground,
                egui::Id::unique("drop-overlay"),
            ));
            painter.rect_filled(screen, 0.0, t.selection);
            painter.text(
                screen.center(),
                egui::Align2::CENTER_CENTER,
                "Drop to open",
                egui::FontId::proportional(22.0),
                t.text,
            );
        }
    }
}

/// What a status-bar composition summary was computed for.
type InfoKey = (Option<u64>, u64, Option<snappiesnap_engine::doc::Range>);

/// The strings of array `key` in a command result.
fn strings(v: &Value, key: &str) -> Vec<String> {
    v.get(key)
        .and_then(Value::as_array)
        .map(|a| {
            a.iter()
                .filter_map(Value::as_str)
                .map(str::to_string)
                .collect()
        })
        .unwrap_or_default()
}

/// What chosen sequence files are for.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum AlignFor {
    /// The one file of the Align view's pairwise alignment.
    Pairwise,
    /// The Align view's multiple alignment list.
    Multiple,
    /// The sequence view's "Align with" list of the active document.
    AlignWith,
}

/// How one open document is being viewed (kept per tab while it is in the background).
struct TabUi {
    view: View,
    map_style: MapStyle,
    gel: crate::state::GelSetup,
    seq_view: SeqViewState,
    enzyme_view: crate::enzyme_view::EnzymeViewState,
    primer_design: crate::primer_design::DesignState,
    gel_inputs: crate::gel_view::GelInputs,
    align: crate::align_view::AlignState,
    find: crate::find_bar::FindState,
}

impl TabUi {
    /// A document seen for the first time: its map, and empty panels.
    fn fresh() -> Self {
        TabUi {
            view: View::Map,
            map_style: MapStyle::Auto,
            gel: Default::default(),
            seq_view: Default::default(),
            enzyme_view: Default::default(),
            primer_design: Default::default(),
            gel_inputs: Default::default(),
            align: Default::default(),
            find: Default::default(),
        }
    }
}

/// Split the sequence view's room for Align with: the sequence in the top half, the reads in
/// the bottom half, and a strip of `nav_h` under both for the navigator. Whole points, so the
/// parts don't shift with rounding; returns (sequence, (reads, navigator)).
fn split_for_alignment(full: Rect, nav_h: f32) -> (Rect, (Rect, Rect)) {
    let nav_h = if nav_h.is_finite() {
        nav_h.max(0.0).min(full.height().max(0.0))
    } else {
        0.0
    };
    let nav_top = (full.bottom() - nav_h).max(full.top());
    let split = (full.top() + ((nav_top - full.top()) / 2.0).floor())
        .max(full.top())
        .min(nav_top);
    let top = Rect::from_min_max(full.min, egui::pos2(full.right(), split));
    let reads = Rect::from_min_max(
        egui::pos2(full.left(), split),
        egui::pos2(full.right(), nav_top),
    );
    let nav = Rect::from_min_max(egui::pos2(full.left(), nav_top), full.max);
    (top, (reads, nav))
}

#[cfg(all(test, not(target_arch = "wasm32")))]
mod tests {
    use snappiesnap_engine::script;

    use super::*;

    /// Edits made in the app with `--record` replay to the same document.
    #[test]
    fn recorded_sessions_replay_exactly() {
        let path =
            std::env::temp_dir().join(format!("snappiesnap-record-{}.jsonl", std::process::id()));
        let options = LaunchOptions {
            control_port: None,
            record: Some(path.to_string_lossy().into_owned()),
            remember_settings: false,
        };
        let mut app = App::new(&egui::Context::default(), Session::new(), &options);
        app.run(
            "new",
            json!({"residues": "ACGTACGTAC", "topology": "circular"}),
        );
        app.run("insert", json!({"at": 4, "residues": "GG"}));
        app.run(
            "add_feature",
            json!({"feature": {"name": "f", "segments": [{"range": {"start": 9, "end": 2}}]}}),
        );
        app.run("insert", json!({"at": 99, "residues": "A"})); // fails: not recorded
        app.run("undo", json!({}));
        app.run("redo", json!({}));
        app.write_recording();
        drop(app.recorder.take());

        let text = std::fs::read_to_string(&path).unwrap();
        let _ = std::fs::remove_file(&path);
        assert_eq!(text.lines().count(), 5, "{text}");
        let mut replay = Session::new();
        let summary = script::run(&mut replay, &text, |_| true);
        assert_eq!(summary.failed, 0);
        assert_eq!(
            replay.document().map(|d| &d.molecule),
            app.session().document().map(|d| &d.molecule)
        );
    }

    /// Web build: files the browser hands over are opened from memory, the same name twice
    /// opens a second copy, and unreadable or unknown files are reported, not opened.
    #[test]
    fn browser_files_open_from_memory() {
        use crate::web_files::Incoming;
        let platform = std::sync::Arc::new(snappiesnap_engine::MemoryFs::new());
        let mut app = App::headless(Session::new().with_platform(platform.clone()));
        let fasta = b">pWeb\nACGTACGTACGTAAAA\n".to_vec();
        for (name, bytes) in [
            ("C:\\lab\\pWeb.fa", Ok(fasta.clone())),
            ("pWeb.fa", Ok(fasta.clone())),
            ("notes.docx", Ok(b"PK".to_vec())),
            (
                "broken.gb",
                Err("the browser couldn't read broken.gb".to_string()),
            ),
        ] {
            app.web_inbox.push(Incoming {
                name: name.into(),
                bytes,
            });
        }
        app.take_web_files();
        let paths: Vec<_> = app
            .session()
            .documents()
            .iter()
            .filter_map(|d| d.path.clone())
            .collect();
        assert_eq!(paths, ["pWeb.fa", "pWeb (2).fa"]);
        assert_eq!(app.session().document().map(|d| d.molecule.len()), Some(16));
        assert!(platform.get("notes.docx").is_none());
        let log = app.session().log_tail(10);
        assert!(
            log.iter().any(|e| e.message.contains("notes.docx")),
            "{log:?}"
        );
        assert!(
            log.iter().any(|e| e.message.contains("broken.gb")),
            "{log:?}"
        );
        assert!(app.web_inbox.take().is_empty());
    }

    #[test]
    fn headless_render_draws_every_view() {
        let mut session = Session::new();
        session
            .execute(
                "new",
                json!({"residues": "ACGT".repeat(50), "topology": "circular"}),
            )
            .unwrap();
        let mut app = App::headless(session);
        for view in View::ALL {
            app.set_ui(crate::state::UiPatch {
                view: Some(view),
                ..Default::default()
            });
            let image = app.render(400, 300, 1.0).unwrap();
            assert_eq!(image.size, [400, 300], "{view:?}");
        }
    }

    /// The map is saved as SVG (vector, labels as text) and as PNG, at the default size for
    /// its shape.
    #[test]
    fn map_export() {
        let mut session = Session::new();
        session
            .execute(
                "new",
                json!({"name": "pExport", "residues": format!("{}GAATTC{}", "A".repeat(500), "T".repeat(500)), "topology": "circular"}),
            )
            .unwrap();
        session
            .execute(
                "add_feature",
                json!({"feature": {"name": "my gene", "kind": "CDS", "strand": "forward", "segments": [{"range": {"start": 100, "end": 400}}], "qualifiers": []}}),
            )
            .unwrap();
        let mut app = App::headless(session);
        assert_eq!(app.map_export_size(), (1000, 1000));
        let dir = std::env::temp_dir().join(format!("snappiesnap-map-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let svg = dir.join("m.svg").to_string_lossy().into_owned();
        let r = app.save_map(&svg, 1000, 1000).unwrap();
        assert_eq!(r["format"], "svg");
        let text = std::fs::read_to_string(&svg).unwrap();
        assert!(text.contains("<title>pExport map</title>"));
        assert!(text.contains(">pExport</text>"), "the name in the middle");
        assert!(text.contains(">my gene</text>"), "the feature label");
        assert!(text.contains("EcoRI"), "the enzyme label");
        let png = dir.join("m.png").to_string_lossy().into_owned();
        app.save_map(&png, 400, 300).unwrap();
        assert!(std::fs::read(&png).unwrap().starts_with(b"\x89PNG"));
        assert!(
            app.save_map(&dir.join("m.gif").to_string_lossy(), 400, 300)
                .is_err()
        );
        assert!(app.save_map(&svg, 10, 10).is_err());
        let _ = std::fs::remove_dir_all(&dir);
    }

    /// The Go to bar draws in both views, and its error stays in the bar. Set
    /// `SNAPPIESNAP_GOTO_PNG=/path.png` to look at it.
    #[test]
    fn go_to_bar() {
        let mut session = Session::new();
        session
            .execute("new", json!({"name": "a", "residues": "ACGT".repeat(50)}))
            .unwrap();
        let mut app = App::headless(session);
        for view in [View::Map, View::Sequence] {
            app.state.view = view;
            app.goto.open();
            app.goto.text = "999".into();
            app.goto.error = Some("positions run from 1 to 200, not 999".into());
            let image = app.render(900, 400, 1.0).unwrap();
            if let Ok(path) = std::env::var("SNAPPIESNAP_GOTO_PNG") {
                crate::offscreen::write_png(&format!("{path}.{view:?}.png"), &image).unwrap();
            }
            assert!(app.goto.open);
        }
        // The summary follows the selection.
        app.session
            .execute("go_to", json!({"location": "1-20"}))
            .unwrap();
        app.render(400, 300, 1.0).unwrap();
        let (_, parts) = app.info.clone().unwrap();
        assert_eq!(parts.first().map(String::as_str), Some("GC 50.0%"));
        assert!(parts.iter().any(|p| p.starts_with("Tm ")), "{parts:?}");
    }

    /// The Preferences window draws, and the background it sets applies to the Map and
    /// Sequence views only. Set `SNAPPIESNAP_PREFERENCES_PNG=/path.png` to look at it.
    #[test]
    fn map_background_from_preferences() {
        let mut session = Session::new();
        session
            .execute("new", json!({"name": "a", "residues": "ACGT".repeat(50)}))
            .unwrap();
        let mut app = App::headless(session);
        app.state.theme = ThemeChoice::Dark;
        assert_eq!(app.canvas_tokens().canvas, Tokens::dark().canvas);
        let bg = egui::Color32::from_rgb(0x20, 0x30, 0x40);
        app.state.map_background.set(true, Some(bg));
        app.state.preferences = true;
        let image = app.render(900, 600, 1.0).unwrap();
        if let Ok(path) = std::env::var("SNAPPIESNAP_PREFERENCES_PNG") {
            crate::offscreen::write_png(&path, &image).unwrap();
        }
        assert_eq!(app.canvas_tokens().canvas, bg);
        // The chosen dark colour leaves the light theme alone, and is remembered.
        app.state.theme = ThemeChoice::Light;
        app.render(200, 150, 1.0).unwrap();
        assert_eq!(app.canvas_tokens().canvas, Tokens::light().canvas);
        assert_eq!(
            app.settings().map_background.dark.as_deref(),
            Some("#203040")
        );
    }

    /// The title names the document and the version that is running.
    #[test]
    fn title_shows_the_version() {
        let mut app = App::headless(Session::new());
        assert_eq!(
            app.window_title(),
            format!("SnappieSnap {}", env!("CARGO_PKG_VERSION"))
        );
        app.run("new", json!({"name": "pT", "residues": "ACGT"}));
        let title = app.window_title();
        assert!(title.starts_with("pT"), "{title}");
        assert!(title.ends_with(&format!("— SnappieSnap {}", env!("CARGO_PKG_VERSION"))));
    }

    /// The Gel tab's gel and an experiment's gel both become images, as shown and inverted.
    #[test]
    fn gels_become_images() {
        let mut session = Session::new();
        session
            .execute(
                "new",
                json!({"name": "pG", "residues": format!("{0}GAATTC{0}", "A".repeat(600)), "topology": "circular"}),
            )
            .unwrap();
        let mut app = App::headless(session);
        app.render(400, 300, 1.0).unwrap();
        let (shown, name) = app.gel_image(GelSource::Gel, false).unwrap();
        assert_eq!(name, "pG restriction digest gel");
        let (inverted, _) = app.gel_image(GelSource::Gel, true).unwrap();
        assert_eq!(shown.size, inverted.size);
        assert_ne!(shown.pixels[0], inverted.pixels[0]);
        let id = app.session.documents()[0].id;
        app.state
            .new_experiment(crate::state::ExperimentKind::Gel, [id, id]);
        app.state.experiments[0].name = "Run 1/2".into();
        let (_, name) = app.gel_image(GelSource::Experiment(0), false).unwrap();
        assert_eq!(name, "Run 1-2");
        assert!(app.gel_image(GelSource::Experiment(5), false).is_err());
        // The question before saving draws; set `SNAPPIESNAP_GEL_DIALOG_PNG=/path.png` to look.
        app.gel_export = Some(GelSource::Gel);
        app.state.active_experiment = None;
        app.state.view = View::Gel;
        let image = app.render(900, 600, 1.0).unwrap();
        if let Ok(path) = std::env::var("SNAPPIESNAP_GEL_DIALOG_PNG") {
            crate::offscreen::write_png(&path, &image).unwrap();
        }
        assert_eq!(app.gel_export, Some(GelSource::Gel), "waits for an answer");
    }

    #[test]
    fn each_tab_remembers_how_it_was_viewed() {
        let mut session = Session::new();
        session
            .execute("new", json!({"name": "a", "residues": "ACGT"}))
            .unwrap();
        let mut app = App::headless(session);
        let set = |app: &mut App, view| {
            app.set_ui(crate::state::UiPatch {
                view: Some(view),
                ..Default::default()
            });
            app.render(200, 150, 1.0).unwrap();
        };
        set(&mut app, View::Sequence);
        app.enzyme_view.filter = "EcoRI".into();
        app.enzyme_view.from = Some(crate::enzyme_view::ListFrom::All);
        app.run("new", json!({"name": "b", "residues": "GGGG"}));
        app.render(200, 150, 1.0).unwrap();
        // A new tab starts on its map, with fresh panels, but the same enzyme list.
        assert_eq!(app.state.view, View::Map);
        assert!(app.enzyme_view.filter.is_empty());
        assert_eq!(
            app.enzyme_view.from,
            Some(crate::enzyme_view::ListFrom::All)
        );
        assert_eq!(
            app.settings().enzyme_list,
            Some(crate::enzyme_view::ListFrom::All)
        );
        set(&mut app, View::Gel);
        // An experiment tab spans documents: going to another document shows that document,
        // and the experiment tab stays open.
        app.state
            .new_experiment(crate::state::ExperimentKind::Gel, [1, 2]);
        app.render(200, 150, 1.0).unwrap();
        assert_eq!(app.state.active_experiment, Some(0));
        app.run("switch_document", json!({"index": 0}));
        app.render(200, 150, 1.0).unwrap();
        assert_eq!(app.state.view, View::Sequence);
        assert_eq!(app.enzyme_view.filter, "EcoRI");
        assert_eq!(app.state.active_experiment, None);
        assert_eq!(app.state.experiments[0].lanes.len(), 2);
        app.run("switch_document", json!({"index": 1}));
        app.render(200, 150, 1.0).unwrap();
        assert_eq!(app.state.view, View::Gel);
    }
}

/// Show a file or folder in the platform's file manager.
#[cfg(not(target_arch = "wasm32"))]
fn reveal(path: &str) -> std::io::Result<()> {
    use std::process::Command;
    if cfg!(target_os = "macos") {
        Command::new("open").arg("-R").arg(path).spawn()?;
    } else if cfg!(target_os = "windows") {
        Command::new("explorer")
            .arg(format!("/select,{path}"))
            .spawn()?;
    } else {
        let dir = sidebar::parent_folder(path).unwrap_or_else(|| path.to_string());
        Command::new("xdg-open").arg(dir).spawn()?;
    }
    Ok(())
}

#[cfg(test)]
mod translation_menu_tests {
    use super::*;

    /// The menu draws with every choice; set `SNAPPIESNAP_MENU_PNG=/path.png` to look at it.
    #[test]
    fn translation_menu_draws() {
        let t = Tokens::dark();
        let mut state = UiState {
            translation: Translation::Plus2,
            ..UiState::default()
        };
        let image = crate::offscreen::render(
            320,
            520,
            1.0,
            |ctx| t.apply(ctx),
            |ui| {
                egui::Frame::new()
                    .fill(t.panel)
                    .inner_margin(8)
                    .show(ui, |ui| {
                        let _ = translation_menu(ui, &t, &mut state, false);
                    });
            },
        )
        .unwrap();
        assert_eq!(state.translation, Translation::Plus2);
        if let Ok(path) = std::env::var("SNAPPIESNAP_MENU_PNG") {
            crate::offscreen::write_png(&path, &image).unwrap();
        }
        // Old settings still load.
        let old: UiState = serde_json::from_str(r#"{"translation": "three"}"#).unwrap();
        assert_eq!(old.translation, Translation::Top);
        assert_eq!(Translation::Bottom.frames(), &[-1, -2, -3]);
        assert!(Translation::Minus2.reverse() && !Translation::Six.reverse());
    }

    /// The enzyme menu draws; the set describes itself and turns into the engine query. Set
    /// `SNAPPIESNAP_ENZYME_MENU_PNG=/path.png` to look at it.
    #[test]
    fn enzyme_menu_draws() {
        use crate::state::EnzymeSet;
        let t = Tokens::dark();
        let mut state = UiState::default();
        let mut db = snappiesnap_engine::EnzymeDb::default();
        db.sets.push(snappiesnap_engine::enzymes::UserSet {
            name: "Lab stock".into(),
            enzymes: vec!["EcoRI".into(), "XcmI".into()],
        });
        state.enzyme_set.pool = Some("Lab stock".into());
        let image = crate::offscreen::render(
            320,
            720,
            1.0,
            |ctx| t.apply(ctx),
            |ui| {
                egui::Frame::new()
                    .fill(t.panel)
                    .inner_margin(8)
                    .show(ui, |ui| {
                        enzyme_menu(ui, &t, &mut state, &db, &["Foo1".to_string()]);
                    });
            },
        )
        .unwrap();
        if let Ok(path) = std::env::var("SNAPPIESNAP_ENZYME_MENU_PNG") {
            crate::offscreen::write_png(&path, &image).unwrap();
        }
        assert_eq!(
            state.enzyme_set.describe(),
            "unique 6+ cutters from Lab stock"
        );
        assert_eq!(state.enzyme_set.query().set.as_deref(), Some("Lab stock"));
        let d = EnzymeSet::default();
        assert_eq!(d.describe(), "unique 6+ cutters, sold");
        let q = d.query();
        assert_eq!((q.max_cuts, q.min_length, q.enzymes), (1, 6, None));
        let chosen = EnzymeSet {
            only: vec!["EcoRI".into(), " BamHI ".into(), "".into()],
            max_cuts: 99,
            ..EnzymeSet::default()
        }
        .clamped();
        assert_eq!(chosen.only, vec!["EcoRI", "BamHI"]);
        assert_eq!(chosen.max_cuts, EnzymeSet::MAX_CUTS);
        assert_eq!(chosen.describe(), "chosen enzymes (EcoRI, BamHI)");
        assert_eq!(
            chosen.query().enzymes,
            Some(vec!["EcoRI".to_string(), "BamHI".to_string()])
        );
    }

    /// Zoom stays within its limits, and the pan keeps the map covering the area.
    #[test]
    fn zoom_and_pan_are_kept_in_range() {
        let t = Tokens::light();
        let mut zoom = 500.0f32;
        let mut pan = vec2(1.0e6, -1.0e6);
        let mut drawn = None;
        let _ = crate::offscreen::render(
            400,
            300,
            1.0,
            |ctx| t.apply(ctx),
            |ui| {
                let area = Rect::from_min_size(egui::pos2(0.0, 0.0), vec2(400.0, 300.0));
                drawn = Some(zoomed(ui, "z", area, &mut zoom, &mut pan, true, |ui| {
                    ui.max_rect()
                }));
            },
        )
        .unwrap();
        let rect = drawn.unwrap();
        // 50x at most: the map is 50 times the area, and still covers it.
        assert!((rect.width() - 400.0 * crate::state::MAX_ZOOM).abs() < 1.0);
        assert!(rect.left() <= 0.0 && rect.right() >= 400.0);
        assert!(rect.top() <= 0.0 && rect.bottom() >= 300.0);
        // Linear maps zoom along the sequence only.
        let mut zoom = 3.0f32;
        let mut pan = vec2(0.0, 80.0);
        let mut drawn = None;
        let _ = crate::offscreen::render(
            400,
            300,
            1.0,
            |ctx| t.apply(ctx),
            |ui| {
                let area = Rect::from_min_size(egui::pos2(0.0, 0.0), vec2(400.0, 300.0));
                drawn = Some(zoomed(ui, "z", area, &mut zoom, &mut pan, false, |ui| {
                    ui.max_rect()
                }));
            },
        )
        .unwrap();
        let rect = drawn.unwrap();
        assert!((rect.width() - 1200.0).abs() < 1.0 && (rect.height() - 300.0).abs() < 1.0);
        assert_eq!(pan.y, 0.0);
    }
    /// Simulate Gel covers both gel kinds behind one menu item: every step of each draws,
    /// starting with the step that asks what to simulate.
    #[test]
    fn simulate_gel_steps_draw() {
        use crate::state::{ExperimentKind, ExperimentStep};
        let mut session = Session::new();
        session
            .execute(
                "new",
                json!({"name": "pClone", "residues": "ACGT".repeat(500), "topology": "circular"}),
            )
            .unwrap();
        let mut app = App::headless(session);
        let ids: Vec<u64> = app.session.documents().iter().map(|d| d.id).collect();
        app.state.new_experiment(ExperimentKind::Gel, ids);
        for kind in ExperimentKind::SIMULATIONS {
            let steps = [
                ExperimentStep::Simulation,
                ExperimentStep::Constructs,
                if kind == ExperimentKind::Diagnostic {
                    ExperimentStep::Enzymes
                } else {
                    ExperimentStep::Lanes
                },
            ];
            for step in steps {
                let Some(e) = app.state.experiments.first_mut() else {
                    panic!("the experiment is open");
                };
                e.kind = kind;
                e.step = step;
                e.selected = Some(0);
                let image = app.render(900, 620, 1.0).unwrap();
                assert_eq!(image.size, [900, 620], "{kind:?} {step:?}");
            }
        }
    }
    /// Every step of a cloning experiment draws, and the two constructs are read from the
    /// open folder without any document being open.
    #[test]
    fn cloning_steps_draw() {
        use crate::state::{ExperimentKind, ExperimentStep};
        let fs = snappiesnap_engine::MemoryFs::new()
            .with_file("lab/vector.gb", VECTOR_GB.as_bytes().to_vec())
            .with_file(
                "lab/insert.fa",
                b">insert\nAAGCTTCCGGAATTCGGCCAAGCTT\n".to_vec(),
            );
        let mut session = Session::new().with_platform(std::sync::Arc::new(fs));
        // The lab's own enzyme set: the enzyme fields offer what it holds (HindIII cuts the
        // vector; SpeI doesn't cut it, EcoRI isn't in the set).
        session
            .set_enzyme_sets(vec![snappiesnap_engine::enzymes::UserSet {
                name: "Fridge".into(),
                enzymes: vec!["HindIII".into(), "SpeI".into()],
            }])
            .expect("sets kept");
        let mut app = App::headless(session);
        app.state.folder = Some("lab".to_string());
        app.state.enzyme_set.pool = Some("Fridge".to_string());
        app.state.new_experiment(ExperimentKind::Cloning, []);
        assert_eq!(
            app.state.experiments.first().map(|e| e.step),
            Some(ExperimentStep::Technique)
        );
        for step in [
            ExperimentStep::Technique,
            ExperimentStep::Vector,
            ExperimentStep::Insert,
            ExperimentStep::Product,
        ] {
            let Some(e) = app.state.experiments.first_mut() else {
                panic!("the experiment is open");
            };
            e.step = step;
            e.vector.file = "lab/vector.gb".into();
            e.vector.enzymes = "HindIII".into();
            e.insert.file = "lab/insert.fa".into();
            // Drawing asks for the file and the preview; the next frames show the answers.
            for _ in 0..3 {
                let image = app.render(900, 700, 1.0).unwrap();
                assert_eq!(image.size, [900, 700], "{step:?}");
                if let Ok(path) = std::env::var("SNAPPIESNAP_CLONING_PNG") {
                    let path = path.replace(".png", &format!("_{step:?}.png"));
                    crate::offscreen::write_png(&path, &image).unwrap();
                }
            }
        }
        // Enzymes the lab's set hasn't got, and a name no database knows, are said so under
        // the field (drawn, so the notes can't quietly clip).
        if let Some(e) = app.state.experiments.first_mut() {
            e.step = ExperimentStep::Vector;
            e.vector.enzymes = "HindIII, EcoRI, Nonesuch".into();
        }
        for _ in 0..3 {
            let image = app.render(900, 700, 1.0).unwrap();
            if let Ok(path) = std::env::var("SNAPPIESNAP_CLONING_PNG") {
                crate::offscreen::write_png(&path.replace(".png", "_notes.png"), &image).unwrap();
            }
        }
        if let Some(e) = app.state.experiments.first_mut() {
            e.step = ExperimentStep::Product;
            e.vector.enzymes = "HindIII".into();
        }
        // The product of the reaction is opened as a construct of its own.
        assert!(app.session.documents().is_empty());
        let Some((command, params)) = crate::cloning_experiment::call(
            app.state.experiments.first().expect("the experiment"),
            false,
        ) else {
            panic!("the reaction is ready to run");
        };
        app.run(command, params);
        assert_eq!(app.session.documents().len(), 1);
    }

    /// The web build's case: no folder at all, the constructs chosen from the open tabs.
    /// Set `SNAPPIESNAP_CLONING_TABS_PNG=/path.png` to look at it.
    #[test]
    fn cloning_uses_the_open_tabs_when_there_is_no_folder() {
        use crate::state::{ExperimentKind, ExperimentStep};
        let session =
            Session::new().with_platform(std::sync::Arc::new(snappiesnap_engine::MemoryFs::new()));
        let mut app = App::headless(session);
        // As the browser does it: files are opened into tabs, and there is no folder.
        app.run(
            "new",
            json!({"name": "pVec", "residues": "GGGGGAAGCTTGGGGGGGGGGGGGGGGGGG", "topology": "circular"}),
        );
        app.run(
            "new",
            json!({"name": "pIns", "residues": "AAGCTTCCGGAATTCGGCCAAGCTT", "topology": "linear"}),
        );
        assert!(app.state.folder.is_none());
        let ids: Vec<u64> = app.session.documents().iter().map(|d| d.id).collect();
        assert_eq!(ids.len(), 2);
        app.state.new_experiment(ExperimentKind::Cloning, []);
        let Some(e) = app.state.experiments.first_mut() else {
            panic!("the experiment is open");
        };
        e.step = ExperimentStep::Insert;
        e.vector.document = Some(ids[0]);
        e.vector.enzymes = "HindIII".into();
        e.insert.document = Some(ids[1]);
        // The reaction is ready with no file named anywhere.
        let Some((command, params)) =
            crate::cloning_experiment::call(app.state.experiments.first().expect("open"), false)
        else {
            panic!("the reaction is ready to run");
        };
        assert_eq!(command, "restriction_clone");
        assert_eq!(params["vector_document"], json!(ids[0]));
        for _ in 0..3 {
            let image = app.render(900, 700, 1.0).unwrap();
            if let Ok(path) = std::env::var("SNAPPIESNAP_CLONING_TABS_PNG") {
                crate::offscreen::write_png(&path, &image).unwrap();
            }
        }
        // And it runs: the product opens as a third tab.
        app.run(command, params);
        assert_eq!(app.session.documents().len(), 3);
    }

    /// The web version's first run: no folder, no files, one click on "Try it with the
    /// sample constructs" and the experiment is ready to produce a construct. Set
    /// `SNAPPIESNAP_SAMPLE_CLONING_PNG=/path.png` to look at it.
    #[test]
    fn the_sample_cloning_experiment_runs_to_a_product() {
        use crate::state::{ExperimentKind, ExperimentStep};
        let session =
            Session::new().with_platform(std::sync::Arc::new(snappiesnap_engine::MemoryFs::new()));
        let mut app = App::headless(session);
        app.state.new_experiment(ExperimentKind::Cloning, []);
        assert!(app.session.documents().is_empty());
        app.open_sample_cloning(0);

        let Some(e) = app.state.experiments.first() else {
            panic!("the experiment is open");
        };
        assert_eq!(e.step, ExperimentStep::Vector);
        assert_eq!(e.technique, crate::state::CloningTechnique::Restriction);
        assert!(e.vector.document.is_some() && e.insert.document.is_some());
        assert!(e.vector.file.is_empty(), "no file is needed anywhere");
        assert_eq!(app.session.documents().len(), 2);
        // Neither sample counts as unsaved work, so the page can be closed without a warning.
        assert!(app.session.documents().iter().all(|d| !d.dirty));

        let Some((command, params)) = crate::cloning_experiment::call(e, false) else {
            panic!("the sample experiment is ready to run");
        };
        for step in [
            ExperimentStep::Vector,
            ExperimentStep::Insert,
            ExperimentStep::Product,
        ] {
            if let Some(e) = app.state.experiments.first_mut() {
                e.step = step;
            }
            // The construct pane offers what a document tab offers: the sequence view is
            // drawn for the insert step, the map for the others.
            app.state.cloning_sequence = step == ExperimentStep::Insert;
            app.state.orfs = std::env::var("SNAPPIESNAP_NO_ORFS").is_err();
            // A large viewport is far bigger than the default: the step's buttons have
            // to stay on screen at both sizes.
            let (w, h) = match std::env::var("SNAPPIESNAP_CLONING_SIZE").as_deref() {
                Ok("large") => (1400, 1100),
                _ => (900, 700),
            };
            for _ in 0..3 {
                let image = app.render(w, h, 1.0).unwrap();
                if let Ok(path) = std::env::var("SNAPPIESNAP_SAMPLE_CLONING_PNG") {
                    let path = path.replace(".png", &format!("_{step:?}.png"));
                    crate::offscreen::write_png(&path, &image).unwrap();
                }
            }
        }
        app.run(command, params);
        assert_eq!(app.session.documents().len(), 3);
        let product = app.session.document().expect("the product");
        assert_eq!(product.molecule.len(), 3950);
    }

    /// The Align with pane holds still while the pointer moves and the window is redrawn.
    #[test]
    fn align_with_pane_holds_still() {
        let mut x: u32 = 777;
        let c: String = (0..1500)
            .map(|_| {
                x = x.wrapping_mul(1_103_515_245).wrapping_add(12345);
                b"ACGT"[((x >> 16) % 4) as usize] as char
            })
            .collect();
        let mut session = Session::new();
        session
            .execute(
                "new",
                json!({"name": "pHold", "residues": c, "topology": "circular"}),
            )
            .unwrap();
        let construct = session.active_document_id().unwrap();
        let mut reads = Vec::new();
        for from in [100, 400, 700] {
            // The middle read differs from the construct at 431.
            let mut read = c[from..from + 500].to_string();
            if from == 400 {
                let was = read.as_bytes()[30];
                read.replace_range(30..31, if was == b'A' { "C" } else { "A" });
            }
            session
                .execute("new", json!({"name": format!("r{from}"), "residues": read}))
                .unwrap();
            reads.push(session.active_document_id().unwrap());
        }
        session.switch_document(0).unwrap();
        let mut app = App::headless(session);
        app.set_ui(crate::state::UiPatch {
            view: Some(View::Sequence),
            sequence_single_line: Some(true),
            ..Default::default()
        });
        app.state.align_panel = true;
        if let Some(list) = app.state.align_list(construct) {
            for id in reads {
                list.push(crate::state::AlignEntry {
                    source: crate::state::AlignSource::Document(id),
                    on: true,
                    collapsed: false,
                });
            }
        }
        for scale in [1.0_f32, 4.0 / 3.0] {
            let ctx = egui::Context::default();
            app.tokens.apply(&ctx);
            let mut tops = Vec::new();
            for frame in 0..60u32 {
                let mut input = egui::RawInput {
                    screen_rect: Some(egui::Rect::from_min_size(
                        egui::Pos2::ZERO,
                        egui::vec2(1400.0, 900.0),
                    )),
                    time: Some(f64::from(frame) * 0.05),
                    ..Default::default()
                };
                input
                    .viewports
                    .entry(egui::ViewportId::ROOT)
                    .or_default()
                    .native_pixels_per_point = Some(scale);
                input.events.push(egui::Event::PointerMoved(egui::pos2(
                    300.0 + frame as f32 * 13.0,
                    200.0 + (frame % 7) as f32 * 60.0,
                )));
                let mut out = ctx.run_ui(input, |ui| app.draw(ui));
                out.textures_delta.clear();
                tops.push(app.align_with.pane_rect.map(|r| (r.top(), r.bottom())));
            }
            let settled = &tops[10..];
            assert!(
                settled.iter().all(|r| r.is_some() && *r == settled[0]),
                "scale {scale}: the pane moves: {tops:?}"
            );
            // The reads take the lower half, over the navigator, on the sequence's columns.
            let (top, bottom) = settled[0].unwrap_or_default();
            let screen_h = 900.0;
            assert!(
                top > screen_h * 0.35 && top < screen_h * 0.6,
                "{top} of {screen_h}"
            );
            assert!(
                bottom < screen_h - 20.0,
                "the navigator is under the pane: {bottom}"
            );
            assert!(app.seq_view.line.is_some(), "the sequence is on one line");
        }
        // The font size preference widens the columns of the sequence and the reads alike.
        let width_at = |app: &mut App, size: f32| {
            app.state.sequence_font = size;
            let _ = app.render(1400, 900, 1.0).unwrap();
            app.seq_view.line.map(|l| (l.char_w, l.font_size))
        };
        let small = width_at(&mut app, 10.0).unwrap();
        let large = width_at(&mut app, 20.0).unwrap();
        assert_eq!((small.1, large.1), (10.0, 20.0));
        assert!((large.0 / small.0 - 2.0).abs() < 0.1, "{small:?} {large:?}");
        app.state.sequence_font = crate::state::SEQUENCE_FONT_DEFAULT;
        if let Ok(path) = std::env::var("SNAPPIESNAP_ALIGN_PANE_PNG") {
            // Each render starts afresh, scrolled to the start: ask for the reads' stretch.
            app.seq_view.scroll_to = Some(420);
            app.state.sequence_font = 18.0;
            let image = app.render(1400, 900, 4.0 / 3.0).unwrap();
            crate::offscreen::write_png(&path, &image).unwrap();
            app.set_ui(crate::state::UiPatch {
                theme: Some(crate::theme::ThemeChoice::Light),
                ..Default::default()
            });
            app.seq_view.scroll_to = Some(420);
            let image = app.render(1400, 900, 4.0 / 3.0).unwrap();
            crate::offscreen::write_png(&path.replace(".png", "_light.png"), &image).unwrap();
        }
    }

    const VECTOR_GB: &str = "LOCUS       vector            60 bp    DNA     circular SYN 01-JAN-2026\nORIGIN\n        1 aagcttccgg aattcggcca agcttacgta cgtacgtacg tacgtacgta cgtacgtacg\n//\n";
}
