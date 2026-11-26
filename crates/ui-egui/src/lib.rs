//! SnappieSnap `ui-egui` (layer L5).
//!
//! The thin egui shell. It reads engine state and acts only by running engine commands, so the
//! UI, the control channel and MCP behave identically. Colours and metrics come from
//! [`theme::Tokens`]; icons are drawn in code.
//!
//! Layering is enforced by `cargo xtask layers`; see `AGENTS.md`.

pub mod align_view;
pub mod align_with;
mod app;
mod browser;
mod cloning_experiment;
mod cloning_plan;
pub mod color_picker;
#[cfg(not(target_arch = "wasm32"))]
mod control;
pub mod enzyme_view;
pub mod enzymes;
pub mod experiment_view;
pub mod features_view;
pub mod find_bar;
pub mod gel_export;
pub mod gel_view;
pub mod goto_bar;
pub mod history_diagram;
pub mod history_view;
mod icons;
mod map;
pub mod msa_view;
pub mod notes_view;
pub mod offscreen;
pub mod oligo_view;
pub mod preferences;
pub mod primer_design;
pub mod reads_view;
mod seq_navigator;
mod seqview;
pub mod settings;
pub mod sidebar;
pub mod state;
pub mod svg;
/// How wide the left column is when a view splits its inputs or summary from the wide thing it
/// produces (an alignment, a grid of reads), or `None` when the page is too narrow to split.
/// Shared so the views that do it agree on where the line is.
pub(crate) fn two_columns(full: f32) -> Option<f32> {
    (full >= 1000.0).then(|| (full * 0.42).clamp(380.0, 620.0))
}

/// How many residues of an alignment fit across `room` pixels at `char_w` per character: as
/// many as there is space for, in tens so the positions above them stay readable, and never so
/// few or so many that a block is hard to follow.
pub(crate) fn alignment_columns(room: f32, char_w: f32) -> usize {
    if !room.is_finite() || !char_w.is_finite() || char_w <= 0.0 {
        return 60;
    }
    let fits = (room / char_w).max(0.0) as usize;
    (fits / 10 * 10).clamp(30, 150)
}

pub mod theme;
pub mod trace_view;
pub mod updates;
#[cfg(target_arch = "wasm32")]
mod web;
pub mod web_files;

pub use app::{App, LaunchOptions};
pub use eframe;

#[cfg(test)]
mod layout_tests {
    /// Views that put a wide result beside their inputs agree on when to split and how much
    /// each side gets, and a block of an alignment never ends up unreadably narrow or wide.
    #[test]
    fn the_shared_layout_rules_hold() {
        assert_eq!(super::two_columns(999.0), None);
        for full in [1000.0, 1400.0, 2400.0, 4000.0] {
            let left = super::two_columns(full).expect("two columns");
            assert!(
                (380.0..=620.0).contains(&left) && full - left >= 380.0,
                "{full} wide gives {left}"
            );
        }
        // Columns come in tens, within bounds, and grow with the room.
        assert_eq!(super::alignment_columns(0.0, 7.2), 30);
        assert_eq!(super::alignment_columns(432.0, 7.2), 60);
        assert_eq!(super::alignment_columns(100_000.0, 7.2), 150);
        assert!(super::alignment_columns(900.0, 7.2) > super::alignment_columns(500.0, 7.2));
        // A width that isn't a measurement (a scroll area that hasn't been laid out gives an
        // infinite one) falls back to a sensible block instead of dividing by it.
        assert_eq!(super::alignment_columns(f32::NAN, 7.2), 60);
        assert_eq!(super::alignment_columns(f32::INFINITY, 7.2), 60);
        assert_eq!(super::alignment_columns(400.0, 0.0), 60);
    }
}
