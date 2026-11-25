//! SnappieSnap `bio` (layer L2).
//!
//! Molecular biology algorithms on the document model: restriction enzymes (site finding,
//! digests), and later Tm and primers, ORFs and translation, alignment and cloning simulation.
//! Every algorithm names its source; every function returns errors rather than panicking.
//!
//! Layering is enforced by `cargo xtask layers`; see `AGENTS.md`.

pub mod align;
pub mod assembly;
pub mod cloning;
pub mod codes;
pub mod composition;
pub mod design;
pub mod diagnostic;
pub mod enzymes;
pub mod gel;
pub mod msa;
pub mod pcr;
pub mod primers;
pub mod reads;
pub mod search;
pub mod structure;
pub mod tm;
pub mod topo;
pub mod trace;
pub mod translate;
