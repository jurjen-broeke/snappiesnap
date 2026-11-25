//! SnappieSnap `doc` (layer L1).
//!
//! The molecule document: a sequence with its topology, features, primers and notes. Pure data
//! with no I/O. Every edit keeps annotations consistent with the sequence (features move,
//! shrink or disappear with the residues they cover) and returns a [`DocError`] rather than
//! panicking.
//!
//! The construction history (how a molecule was made, with copies of its ancestors) is a tree
//! of [`HistoryNode`]s; undo/redo lives in the engine.
//!
//! Layering is enforced by `cargo xtask layers`; see `AGENTS.md`.

mod error;
mod feature;
mod history;
mod molecule;
mod trace;

pub use error::DocError;
pub use feature::{
    Feature, FeatureId, NewFeature, Primer, PrimerId, Qualifier, Segment, Strand, is_color,
};
pub use history::{HistoryNode, HistoryStep, MAX_DEPTH as MAX_HISTORY_DEPTH, operation_label};
pub use molecule::{Detached, EditReport, Molecule, Notes};
pub use trace::{Call, Channel, MAX_SAMPLES as MAX_TRACE_SAMPLES, Signal, Trace};

pub use snappiesnap_seq::{Alphabet, Range, Seq, SeqDelta, SeqError, Topology};
