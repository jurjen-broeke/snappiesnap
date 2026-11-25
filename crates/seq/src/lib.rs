//! SnappieSnap `seq` (layer L0).
//!
//! Alphabets and IUPAC codes, validated residue sequences, molecule topology, and ranges that
//! may wrap past the origin of a circular molecule. Every operation checks its inputs and
//! returns a [`SeqError`] instead of panicking.
//!
//! Layering is enforced by `cargo xtask layers`; see `AGENTS.md`.

mod alphabet;
mod error;
mod range;
mod seq;

pub use alphabet::Alphabet;
pub use error::SeqError;
pub use range::{Range, Topology};
pub use seq::{MAX_LEN, Seq, SeqDelta};
