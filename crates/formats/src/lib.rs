//! SnappieSnap `formats` (layer L0).
//!
//! Standalone, byte-preserving readers and writers for sequence file formats. Each reader keeps
//! what it doesn't understand, and each writer reproduces unchanged input byte for byte
//! (AGENTS.md › File fidelity). Readers treat every byte as hostile: they return errors, never
//! panic, and never allocate more than the input size.
//!
//! This crate has no workspace dependencies, so it can be published on its own.
//!
//! Layering is enforced by `cargo xtask layers`; see `AGENTS.md`.

pub mod abif;
pub mod dna;
pub mod embl;
pub mod fasta;
pub mod genbank;
pub mod gff3;
pub mod sbol;
pub mod scf;
pub mod trace;
pub mod xml;
pub mod ztr;
