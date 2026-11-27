//! The Preferences window: categories in a list on the left, the chosen category's settings
//! on the right, in titled groups.
//!
//! To add a category, add a [`PrefsPage`] variant (with its label and summary) and a page
//! function in [`page`]. To add a setting, add a [`row`] to a [`group`] of its page; the value
//! lives in `UiState` (so `ui.state`/`ui.set` can read and drive it) and is remembered in
//! `settings::Settings`.

use eframe::egui::{self, Align, CornerRadius, Layout, RichText, Sense, Stroke, Ui, vec2};
use serde::{Deserialize, Serialize};

use snappiesnap_engine::EnzymeDb;
use snappiesnap_engine::bio::primers::{BindingPrefs, ClassRules};

use crate::color_picker::{ColorChoice, color_picker};
use crate::enzyme_view::ListFrom;
use crate::state::{EnzymeSet, Translation, UiState};
use crate::theme::{ThemeChoice, Tokens};

/// A category of the Preferences window.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum PrefsPage {
    #[default]
    General,
    FilesAndTabs,
    Export,
    Dna,
    Rna,
    Gels,
    Enzymes,
    Primers,
}

impl PrefsPage {
    /// In the order of the category list.
    pub const ALL: [PrefsPage; 8] = [
        PrefsPage::General,
        PrefsPage::FilesAndTabs,
        PrefsPage::Export,
        PrefsPage::Dna,
        PrefsPage::Rna,
        PrefsPage::Gels,
        PrefsPage::Enzymes,
        PrefsPage::Primers,
    ];

    pub fn label(self) -> &'static str {
        match self {
            PrefsPage::General => "General",
            PrefsPage::FilesAndTabs => "Files and tabs",
            PrefsPage::Export => "Export",
            PrefsPage::Dna => "DNA",
            PrefsPage::Rna => "RNA",
            PrefsPage::Gels => "Gels",
            PrefsPage::Enzymes => "Enzymes",
            PrefsPage::Primers => "Primers",
        }
    }

    /// One line under the page title.
    pub fn summary(self) -> &'static str {
        match self {
            PrefsPage::General => "The window, its colours, and updates.",
            PrefsPage::FilesAndTabs => "Opening and saving files, and the document tabs.",
            PrefsPage::Export => "Images and files written for use elsewhere.",
            PrefsPage::Dna => "DNA sequences, their maps and the sequence view.",
            PrefsPage::Rna => "RNA sequences.",
            PrefsPage::Gels => "Simulated agarose gels.",
            PrefsPage::Enzymes => "Restriction enzymes and digests.",
            PrefsPage::Primers => "Primer design and binding.",
        }
    }
}

/// Smallest and largest interface size (egui's zoom factor).
pub const UI_SCALE_RANGE: std::ops::RangeInclusive<f32> = 0.75..=2.0;
/// Sizes the sequence font can take, in points.
pub const SEQUENCE_FONT_RANGE: std::ops::RangeInclusive<f32> = 9.0..=24.0;

/// Width of the category list and of the labels in front of settings.
const LIST_WIDTH: f32 = 160.0;
const LABEL_WIDTH: f32 = 170.0;

/// What Preferences shows besides `UiState`: the app's data it describes, and settings
/// kept outside `UiState`.
pub struct Extra<'a> {
    pub enzymes: &'a EnzymeDb,
    /// Which enzymes the Enzymes tab lists (`None`: the map's set).
    pub enzyme_list: &'a mut Option<ListFrom>,
    /// A REBASE download is running.
    pub downloading: bool,
    /// An update check is running, and the last one's answer.
    pub checking: bool,
    pub update: Option<&'a str>,
    /// Which primer binding sites count (the session's; the app applies a change).
    pub primer_binding: &'a mut BindingPrefs,
}

/// What the user asked for that the app carries out (files, downloads, commands).
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum PrefsAction {
    ChooseFolder,
    ImportEnzymeSet,
    RemoveEnzymeSet(String),
    DownloadRebase,
    UseBuiltinEnzymes,
    CheckForUpdates,
}

/// Draw the Preferences window if it is open; closing it clears `state.preferences`.
pub fn window(
    ctx: &egui::Context,
    t: &Tokens,
    state: &mut UiState,
    extra: Extra<'_>,
) -> Option<PrefsAction> {
    if !state.preferences {
        return None;
    }
    let mut extra = extra;
    let mut action = None;
    let screen = ctx.content_rect();
    // Fits small windows: never taller or wider than the app's window allows.
    let size = vec2(
        (screen.width() - 40.0).clamp(320.0, 760.0),
        (screen.height() - 160.0).clamp(220.0, 500.0),
    );
    let mut open = true;
    egui::Window::new("Preferences")
        .open(&mut open)
        .pivot(egui::Align2::CENTER_TOP)
        .default_pos(screen.center_top() + vec2(0.0, 100.0))
        .collapsible(false)
        .resizable(false)
        .fixed_size(size)
        .show(ctx, |ui| {
            ui.horizontal_top(|ui| {
                ui.allocate_ui_with_layout(
                    vec2(LIST_WIDTH, size.y),
                    Layout::top_down_justified(Align::Min),
                    |ui| {
                        ui.set_min_size(vec2(LIST_WIDTH, size.y));
                        category_list(ui, t, &mut state.preferences_page);
                    },
                );
                let (rule, _) = ui.allocate_exact_size(vec2(1.0, size.y), Sense::hover());
                ui.painter().line_segment(
                    [rule.center_top(), rule.center_bottom()],
                    Stroke::new(1.0, t.border),
                );
                ui.vertical(|ui| {
                    let width = (size.x - LIST_WIDTH - 2.0 * t.spacing - 1.0).max(120.0);
                    ui.set_width(width);
                    egui::ScrollArea::vertical()
                        .id_salt(("preferences-page", state.preferences_page))
                        .max_height(size.y)
                        .auto_shrink([false, false])
                        .show(ui, |ui| {
                            let shown = state.preferences_page;
                            ui.label(RichText::new(shown.label()).heading().strong());
                            ui.label(RichText::new(shown.summary()).color(t.text_muted));
                            ui.add_space(t.spacing * 2.0);
                            action = page(ui, t, state, &mut extra);
                        });
                });
            });
        });
    if !open {
        state.preferences = false;
    }
    action
}

/// The categories, the chosen one highlighted.
fn category_list(ui: &mut Ui, t: &Tokens, current: &mut PrefsPage) {
    for p in PrefsPage::ALL {
        let chosen = *current == p;
        let (rect, response) =
            ui.allocate_exact_size(vec2(LIST_WIDTH - t.spacing, 28.0), Sense::click());
        let fill = if chosen {
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
        ui.painter().text(
            rect.left_center() + vec2(10.0, 0.0),
            egui::Align2::LEFT_CENTER,
            p.label(),
            egui::FontId::proportional(14.0),
            if chosen { t.accent_text } else { t.text },
        );
        if response.clicked() {
            *current = p;
        }
    }
}

/// The settings of the chosen category.
fn page(
    ui: &mut Ui,
    t: &Tokens,
    state: &mut UiState,
    extra: &mut Extra<'_>,
) -> Option<PrefsAction> {
    match state.preferences_page {
        PrefsPage::General => general(ui, t, state, extra),
        PrefsPage::FilesAndTabs => files_and_tabs(ui, t, state),
        PrefsPage::Dna => {
            dna(ui, t, state);
            None
        }
        PrefsPage::Enzymes => enzymes(ui, t, state, extra),
        PrefsPage::Primers => {
            primers(ui, t, extra.primer_binding, &mut state.primer_tails);
            None
        }
        other @ (PrefsPage::Export | PrefsPage::Rna | PrefsPage::Gels) => {
            group(ui, t, "Nothing here yet", |ui| {
                ui.label(
                    RichText::new(format!(
                        "Settings for {} will be added here in later versions.",
                        other.label()
                    ))
                    .color(t.text_muted),
                );
            });
            None
        }
    }
}

/// General: the window, the colour theme and backgrounds, and updates.
fn general(ui: &mut Ui, t: &Tokens, state: &mut UiState, extra: &Extra<'_>) -> Option<PrefsAction> {
    let mut action = None;
    group(ui, t, "Window", |ui| {
        row(
            ui,
            t,
            "Interface size",
            Some("Text and controls everywhere (also Ctrl/Cmd + and -)."),
            |ui| {
                let mut percent = (state.ui_scale * 100.0).round();
                let range = (UI_SCALE_RANGE.start() * 100.0)..=(UI_SCALE_RANGE.end() * 100.0);
                let r = ui.add(
                    egui::Slider::new(&mut percent, range)
                        .step_by(5.0)
                        .suffix(" %"),
                );
                if r.changed() {
                    state.ui_scale = crate::state::clamp_ui_scale(percent / 100.0);
                }
                if ui
                    .add_enabled(state.ui_scale != 1.0, egui::Button::new("Reset"))
                    .clicked()
                {
                    state.ui_scale = 1.0;
                }
            },
        );
        row(ui, t, "At start", None, |ui| {
            ui.checkbox(
                &mut state.restore_window,
                "Open the window at the size it had last time",
            );
        });
    });
    group(ui, t, "Colour theme and background", |ui| {
        row(ui, t, "Theme", None, |ui| {
            for (choice, label) in [
                (ThemeChoice::System, "Like the system"),
                (ThemeChoice::Light, "Light"),
                (ThemeChoice::Dark, "Dark"),
            ] {
                ui.radio_value(&mut state.theme, choice, label);
            }
        });
        for (dark, label, theme) in [
            (false, "Background, light theme", Tokens::light()),
            (true, "Background, dark theme", Tokens::dark()),
        ] {
            row(ui, t, label, None, |ui| {
                let current = state.map_background.for_theme(dark);
                match color_picker(ui, &theme, ("map-background", dark), current, theme.canvas) {
                    Some(ColorChoice::Set(c)) => state.map_background.set(dark, Some(c)),
                    Some(ColorChoice::Default) => state.map_background.set(dark, None),
                    None => {}
                }
            });
        }
        note(
            ui,
            t,
            "The backgrounds apply to the Map and Sequence views.",
        );
    });
    group(ui, t, "Updates", |ui| {
        row(ui, t, "Version", None, |ui| {
            ui.label(format!("SnappieSnap {}", env!("CARGO_PKG_VERSION")));
        });
        row(ui, t, "New versions", None, |ui| {
            if extra.checking {
                ui.spinner();
                ui.label("Checking\u{2026}");
            } else if ui.button("Check for updates").clicked() {
                action = Some(PrefsAction::CheckForUpdates);
            }
        });
        if let Some(answer) = extra.update {
            ui.add(egui::Label::new(answer).wrap());
        }
        note(
            ui,
            t,
            "SnappieSnap only looks for a new version when you press the button.",
        );
    });
    action
}

/// Files and tabs, in the browser: a web page can't see folders, so there is no folder panel
/// to set up (its switch would do nothing).
#[cfg(target_arch = "wasm32")]
fn files_and_tabs(ui: &mut Ui, t: &Tokens, _state: &mut UiState) -> Option<PrefsAction> {
    group(ui, t, "Folder panel", |ui| {
        ui.label(
            RichText::new(
                "Not available in the browser: a web page can't see the folders on your \
                 computer. Open files with Open, or drop them on the page.",
            )
            .color(t.text_muted),
        );
    });
    None
}

/// Files and tabs: the folder panel.
#[cfg(not(target_arch = "wasm32"))]
fn files_and_tabs(ui: &mut Ui, t: &Tokens, state: &mut UiState) -> Option<PrefsAction> {
    let mut action = None;
    group(ui, t, "Folder panel", |ui| {
        row(ui, t, "Panel", None, |ui| {
            ui.checkbox(&mut state.sidebar, "Show the folder panel");
        });
        row(ui, t, "Folder", None, |ui| {
            let shown = state.folder.as_deref().unwrap_or("None chosen yet");
            ui.add(egui::Label::new(shown).truncate())
                .on_hover_text(shown);
            if ui.button("Choose\u{2026}").clicked() {
                action = Some(PrefsAction::ChooseFolder);
            }
        });
        note(
            ui,
            t,
            "Opening a file from elsewhere shows that file's folder instead.",
        );
    });
    action
}

/// DNA: what the map and the sequence view show, and translations.
fn dna(ui: &mut Ui, t: &Tokens, state: &mut UiState) {
    group(ui, t, "Map and sequence view", |ui| {
        row(ui, t, "Show", None, |ui| {
            ui.vertical(|ui| {
                ui.checkbox(&mut state.features, "Features");
                ui.checkbox(&mut state.enzymes, "Restriction sites");
                ui.checkbox(&mut state.primers, "Primer binding sites");
                ui.checkbox(&mut state.positions, "Cut positions in enzyme labels");
                ui.checkbox(&mut state.ruler, "Ruler (coordinates)");
                ui.checkbox(&mut state.orfs, "Open reading frames and translations");
            });
        });
        note(ui, t, "The same as the buttons beside the map.");
        row(
            ui,
            t,
            "Sequence font size",
            Some("The residues in the sequence view and of the reads aligned under it."),
            |ui| {
                let mut size = state.sequence_font;
                let r = ui.add(
                    egui::Slider::new(&mut size, SEQUENCE_FONT_RANGE)
                        .step_by(1.0)
                        .fixed_decimals(0)
                        .suffix(" pt"),
                );
                if r.changed() {
                    state.sequence_font = crate::state::clamp_sequence_font(size);
                }
                let default = crate::state::SEQUENCE_FONT_DEFAULT;
                if ui
                    .add_enabled(state.sequence_font != default, egui::Button::new("Reset"))
                    .clicked()
                {
                    state.sequence_font = default;
                }
            },
        );
    });
    group(ui, t, "Translation", |ui| {
        row(
            ui,
            t,
            "Reading frames",
            Some("What the ORF layer shows in the sequence view."),
            |ui| {
                egui::ComboBox::from_id_salt("prefs-translation")
                    .selected_text(state.translation.label())
                    .show_ui(ui, |ui| {
                        for (value, label) in Translation::MENU.into_iter().flatten() {
                            ui.selectable_value(&mut state.translation, value, label);
                        }
                    });
            },
        );
        row(
            ui,
            t,
            "Genetic code",
            Some("NCBI translation table."),
            |ui| {
                let name = |c: &snappiesnap_engine::bio::codes::GeneticCode| {
                    format!("{} {}", c.id, c.name)
                };
                let current = snappiesnap_engine::bio::translate::code(state.genetic_code)
                    .map_or_else(|| "1 Standard".to_string(), name);
                egui::ComboBox::from_id_salt("prefs-genetic-code")
                    .selected_text(current)
                    .width(260.0)
                    .height(360.0)
                    .show_ui(ui, |ui| {
                        for c in snappiesnap_engine::bio::codes::CODES {
                            ui.selectable_value(&mut state.genetic_code, c.id, name(c));
                        }
                    });
            },
        );
    });
}

/// Enzymes: which enzymes the map and the Enzymes tab show, the database, and the user's
/// enzyme sets.
fn enzymes(
    ui: &mut Ui,
    t: &Tokens,
    state: &mut UiState,
    extra: &mut Extra<'_>,
) -> Option<PrefsAction> {
    let db = extra.enzymes;
    let mut action = None;
    let database = match &db.version {
        Some(v) => format!("{} {v}", db.source),
        None => db.source.clone(),
    };
    let set_label = |name: &str| {
        db.set(name).map_or_else(
            || name.to_string(),
            |s| format!("{} ({})", s.name, s.enzymes.len()),
        )
    };
    group(ui, t, "Map and sequence view", |ui| {
        let set = &mut state.enzyme_set;
        if !set.only.is_empty() {
            ui.add(
                egui::Label::new(
                    RichText::new(format!(
                        "Showing only the enzymes chosen in the map's enzyme menu ({}); the \
                         choices below apply once that list is cleared there.",
                        set.only.join(", ")
                    ))
                    .color(t.text_muted),
                )
                .wrap(),
            );
        }
        ui.add_enabled_ui(set.only.is_empty(), |ui| {
            row(ui, t, "Enzymes from", None, |ui| {
                let all = format!("All enzymes ({database}, {})", db.enzymes.len());
                let shown = set.pool.as_deref().map_or_else(|| all.clone(), set_label);
                egui::ComboBox::from_id_salt("prefs-enzyme-pool")
                    .selected_text(shown)
                    .width(240.0)
                    .show_ui(ui, |ui| {
                        ui.selectable_value(&mut set.pool, None, all);
                        for user in &db.sets {
                            ui.selectable_value(
                                &mut set.pool,
                                Some(user.name.clone()),
                                set_label(&user.name),
                            );
                        }
                    });
            });
            row(ui, t, "Cutting", None, |ui| {
                for (cuts, label) in [(1, "Once"), (2, "Up to twice"), (3, "Up to 3 times")] {
                    ui.radio_value(&mut set.max_cuts, cuts, label);
                }
            });
            row(
                ui,
                t,
                "Site length",
                Some("Fewest specific bases in the site."),
                |ui| {
                    let label = |n: usize| match n {
                        0 => "Any".to_string(),
                        n => format!("{n}+ bases"),
                    };
                    egui::ComboBox::from_id_salt("prefs-site-length")
                        .selected_text(label(set.min_length))
                        .show_ui(ui, |ui| {
                            for n in [0, 4, 5, 6, 7, 8] {
                                ui.selectable_value(&mut set.min_length, n, label(n));
                            }
                        });
                },
            );
            row(ui, t, "Availability", None, |ui| {
                // A set says what is at hand, so whether an enzyme is sold doesn't matter then.
                ui.add_enabled_ui(set.pool.is_none(), |ui| {
                    ui.checkbox(&mut set.commercial_only, "Commercially available only")
                        .on_disabled_hover_text("An enzyme set lists what you have, sold or not");
                });
            });
            row(ui, t, "", None, |ui| {
                if ui
                    .button("Default")
                    .on_hover_text("Unique 6+ cutters (sold ones, unless an enzyme set is chosen)")
                    .clicked()
                {
                    *set = EnzymeSet {
                        pool: set.pool.clone(),
                        ..EnzymeSet::default()
                    };
                }
            });
        });
    });
    group(ui, t, "Enzymes tab", |ui| {
        row(ui, t, "Lists", None, |ui| {
            let list = &mut *extra.enzyme_list;
            let shown = match list {
                None => "The map's enzyme set".to_string(),
                Some(ListFrom::All) => format!("All enzymes ({database})"),
                Some(ListFrom::Set(name)) => set_label(name),
            };
            egui::ComboBox::from_id_salt("prefs-enzyme-list")
                .selected_text(shown)
                .width(240.0)
                .show_ui(ui, |ui| {
                    ui.selectable_value(list, None, "The map's enzyme set");
                    ui.selectable_value(
                        list,
                        Some(ListFrom::All),
                        format!("All enzymes ({database})"),
                    );
                    for user in &db.sets {
                        ui.selectable_value(
                            list,
                            Some(ListFrom::Set(user.name.clone())),
                            set_label(&user.name),
                        );
                    }
                });
        });
    });
    group(ui, t, "Enzyme database", |ui| {
        row(ui, t, "In use", None, |ui| {
            ui.label(format!("{database}, {} enzymes", db.enzymes.len()));
        });
        row(ui, t, "", None, |ui| {
            if extra.downloading {
                ui.spinner();
                ui.label("Downloading REBASE\u{2026}");
            } else {
                #[cfg(not(target_arch = "wasm32"))]
                {
                    let label = if db.is_builtin() {
                        "Download REBASE"
                    } else {
                        "Download REBASE again"
                    };
                    if ui
                        .button(label)
                        .on_hover_text("The full enzyme database from rebase.neb.com (about 4 MB), kept for later sessions")
                        .clicked()
                    {
                        action = Some(PrefsAction::DownloadRebase);
                    }
                }
                if !db.is_builtin() && ui.button("Use built-in").clicked() {
                    action = Some(PrefsAction::UseBuiltinEnzymes);
                }
            }
        });
    });
    group(ui, t, "Enzyme sets", |ui| {
        if db.sets.is_empty() {
            note(
                ui,
                t,
                "No enzyme sets yet. Import a text file of enzyme names, such as your lab's stock.",
            );
        }
        for user in &db.sets {
            row(ui, t, &user.name, None, |ui| {
                ui.label(
                    RichText::new(format!("{} enzymes", user.enzymes.len())).color(t.text_muted),
                );
                if ui
                    .button("Remove")
                    .on_hover_text(format!(
                        "Delete the enzyme set \u{201c}{}\u{201d}",
                        user.name
                    ))
                    .clicked()
                {
                    action = Some(PrefsAction::RemoveEnzymeSet(user.name.clone()));
                }
            });
        }
        #[cfg(not(target_arch = "wasm32"))]
        if ui
            .button("Import enzyme set\u{2026}")
            .on_hover_text(
                "A text file of enzyme names, such as a set exported from another program",
            )
            .clicked()
        {
            action = Some(PrefsAction::ImportEnzymeSet);
        }
    });
    action
}

/// A whole number with thousands separators: 1,000,000.
fn grouped(n: usize) -> String {
    let digits = n.to_string();
    let mut out = String::with_capacity(digits.len() + digits.len() / 3);
    for (i, c) in digits.chars().enumerate() {
        if i > 0 && (digits.len() - i).is_multiple_of(3) {
            out.push(',');
        }
        out.push(c);
    }
    out
}

/// Primers: which binding sites count, with separate rules on short and on long template
/// sequences.
fn primers(ui: &mut Ui, t: &Tokens, rules: &mut BindingPrefs, tails: &mut crate::state::TailStyle) {
    use crate::state::TailStyle;
    group(ui, t, "Sequence view", |ui| {
        row(
            ui,
            t,
            "5' tails",
            Some(
                "How a primer's 5' tail (the part that doesn't pair with the sequence) is shown \
                 above or below the strand it pairs with.",
            ),
            |ui| {
                ui.radio_value(tails, TailStyle::Simplified, "As a loop with its length");
                ui.radio_value(tails, TailStyle::Full, "Written out");
            },
        );
    });
    group(ui, t, "Short and long sequences", |ui| {
        row(
            ui,
            t,
            "Long sequences",
            Some(
                "Sequences this long or longer, such as chromosomes, follow the rules for long sequences.",
            ),
            |ui| {
                ui.label("From");
                ui.add(
                    egui::DragValue::new(&mut rules.long_from_bp)
                        .range(BindingPrefs::LONG_FROM)
                        .speed(1000.0)
                        .custom_formatter(|v, _| grouped(v as usize))
                        .custom_parser(|text| {
                            text.chars()
                                .filter(|c| !matches!(c, ',' | '.' | ' ' | '_'))
                                .collect::<String>()
                                .parse::<f64>()
                                .ok()
                        })
                        .suffix(" bp"),
                );
            },
        );
        row(ui, t, "", None, |ui| {
            if ui
                .add_enabled(
                    *rules != BindingPrefs::default(),
                    egui::Button::new("Default"),
                )
                .on_hover_text(
                    "Long from 1,000,000 bp. Short sequences: 10 matched bases, an isolated \
                     mismatch allowed. Long sequences: 15 matched bases, no mismatch. At least \
                     40 \u{b0}C for both. Sites matching 15 bases at the 5' end shown",
                )
                .clicked()
            {
                *rules = BindingPrefs::default();
            }
        });
        note(
            ui,
            t,
            "Exact sites over 15 bases or more (or the whole primer) always count. The rules \
             apply to the map and sequence view, the Primers tab, the lab oligo search, and \
             PCR when it also uses weaker sites.",
        );
    });
    let long_from = grouped(rules.long_from_bp);
    class_group(
        ui,
        t,
        &format!("Short sequences (under {long_from} bp)"),
        "short",
        &mut rules.short,
    );
    class_group(
        ui,
        t,
        &format!("Long sequences ({long_from} bp or more)"),
        "long",
        &mut rules.long,
    );
    group(ui, t, "Sites at the 5' end", |ui| {
        row(
            ui,
            t,
            "5' end only",
            Some(
                "Also show sites where the primer's 5' end pairs but its 3' end doesn't: it \
                 anneals there but can't be extended. For short and long sequences.",
            ),
            |ui| {
                ui.checkbox(&mut rules.five_prime_on, "At least");
                ui.add_enabled(
                    rules.five_prime_on,
                    egui::DragValue::new(&mut rules.five_prime_matched)
                        .range(ClassRules::MATCH)
                        .suffix(" bases at the 5' end"),
                );
            },
        );
    });
}

/// The binding rules on one class of template sequences.
fn class_group(ui: &mut Ui, t: &Tokens, title: &str, id: &str, rules: &mut ClassRules) {
    group(ui, t, title, |ui| {
        ui.push_id(id, |ui| {
            row(
                ui,
                t,
                "Matched bases",
                Some("Bases at the primer's 3' end that must pair with the template."),
                |ui| {
                    ui.add(
                        egui::DragValue::new(&mut rules.matched)
                            .range(ClassRules::MATCH)
                            .suffix(" bases"),
                    );
                },
            );
            row(
                ui,
                t,
                "Mismatches",
                Some(
                    "A single mismatched base with pairing bases on both sides, never in the \
                     primer's last 3 bases.",
                ),
                |ui| {
                    ui.checkbox(
                        &mut rules.isolated_mismatches,
                        "Allow one isolated mismatch",
                    );
                },
            );
            row(
                ui,
                t,
                "Minimum Tm",
                Some("Sites whose paired stretch melts lower are left out."),
                |ui| {
                    ui.checkbox(&mut rules.min_tm_on, "At least");
                    ui.add_enabled(
                        rules.min_tm_on,
                        egui::DragValue::new(&mut rules.min_tm)
                            .range(ClassRules::MIN_TM)
                            .speed(0.5)
                            .max_decimals(1)
                            .suffix(" \u{b0}C"),
                    );
                },
            );
        });
    });
}

/// A small muted line under the settings of a group.
fn note(ui: &mut Ui, t: &Tokens, text: &str) {
    ui.add(egui::Label::new(RichText::new(text).small().color(t.text_muted)).wrap());
}

/// A titled group of settings, in a frame.
pub fn group(ui: &mut Ui, t: &Tokens, title: &str, add: impl FnOnce(&mut Ui)) {
    ui.label(RichText::new(title).strong());
    ui.add_space(t.spacing * 0.5);
    egui::Frame::new()
        .stroke(Stroke::new(1.0, t.border))
        .corner_radius(CornerRadius::same(t.radius))
        .inner_margin(egui::Margin::same(10))
        .show(ui, |ui| {
            ui.set_width(ui.available_width());
            ui.spacing_mut().item_spacing.y = t.spacing;
            add(ui);
        });
    ui.add_space(t.spacing * 3.0);
}

/// One setting: a label (with an optional hint under it) and its controls beside it.
pub fn row(ui: &mut Ui, t: &Tokens, label: &str, hint: Option<&str>, add: impl FnOnce(&mut Ui)) {
    ui.horizontal_top(|ui| {
        ui.allocate_ui_with_layout(vec2(LABEL_WIDTH, 0.0), Layout::top_down(Align::Min), |ui| {
            ui.set_width(LABEL_WIDTH);
            ui.label(label);
            if let Some(hint) = hint {
                ui.add(egui::Label::new(RichText::new(hint).small().color(t.text_muted)).wrap());
            }
        });
        ui.horizontal_wrapped(add);
    });
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::App;
    use crate::state::UiPatch;
    use snappiesnap_engine::Session;

    #[test]
    fn thousands_separators() {
        for (n, text) in [
            (0, "0"),
            (999, "999"),
            (1_000, "1,000"),
            (1_000_000, "1,000,000"),
            (12_345_678, "12,345,678"),
        ] {
            assert_eq!(grouped(n), text);
        }
    }

    /// Every category draws, in both themes and in a window too small for it, and `ui.set`
    /// opens it on the chosen one.
    #[test]
    fn every_page_draws() {
        let mut app = App::headless(Session::new());
        for page in PrefsPage::ALL {
            for theme in [ThemeChoice::Light, ThemeChoice::Dark] {
                app.set_ui(UiPatch {
                    preferences: Some(true),
                    preferences_page: Some(page),
                    theme: Some(theme),
                    ..Default::default()
                });
                for (w, h) in [(1000, 760), (320, 200)] {
                    let image = app.render(w, h, 1.0).unwrap();
                    assert_eq!(image.size, [w as usize, h as usize], "{page:?}");
                    // Set `SNAPPIESNAP_PREFERENCES_DIR=/dir` to look at every page.
                    if let Ok(dir) = std::env::var("SNAPPIESNAP_PREFERENCES_DIR")
                        && w > 400
                    {
                        let name = format!("{dir}/{page:?}-{theme:?}.png");
                        crate::offscreen::write_png(&name, &image).unwrap();
                    }
                }
                assert!(app.ui_settings().preferences, "{page:?}");
                assert_eq!(app.ui_settings().preferences_page, page);
            }
        }
    }

    /// Page names are stable (`ui.set`, settings), and the interface size stays in range.
    #[test]
    fn page_names_and_scale() {
        for page in PrefsPage::ALL {
            let name = serde_json::to_value(page).unwrap();
            assert_eq!(serde_json::from_value::<PrefsPage>(name).unwrap(), page);
        }
        assert_eq!(
            serde_json::to_value(PrefsPage::FilesAndTabs).unwrap(),
            "files_and_tabs"
        );
        let mut state = UiState::default();
        for (asked, got) in [(1.5, 1.5), (9.0, 2.0), (0.1, 0.75), (f32::NAN, 1.0)] {
            state.apply(UiPatch {
                ui_scale: Some(asked),
                ..Default::default()
            });
            assert_eq!(state.ui_scale, got, "{asked}");
        }
        // The sequence font: whole points within range, the default for nonsense.
        for (asked, got) in [
            (16.0, 16.0),
            (10.4, 10.0),
            (40.0, 24.0),
            (2.0, 9.0),
            (f32::NAN, 13.0),
            (f32::INFINITY, 13.0),
        ] {
            state.apply(UiPatch {
                sequence_font: Some(asked),
                ..Default::default()
            });
            assert_eq!(state.sequence_font, got, "{asked}");
        }
    }
}
