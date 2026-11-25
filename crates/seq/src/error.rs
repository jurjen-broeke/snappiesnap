use crate::Alphabet;

/// Everything that can go wrong in this crate. Messages are written for the person or agent
/// who has to fix the input.
#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum SeqError {
    #[error("invalid {alphabet} residue {found:?} at position {position}")]
    InvalidResidue {
        alphabet: Alphabet,
        position: usize,
        found: char,
    },

    #[error("position {position} is outside the sequence (length {len})")]
    OutOfBounds { position: usize, len: usize },

    #[error(
        "range {start}..{end} wraps past the origin, which is only possible on a circular molecule"
    )]
    WrapOnLinear { start: usize, end: usize },

    #[error("range {start}..{end} is empty")]
    EmptyRange { start: usize, end: usize },

    #[error("the sequence would be {requested} residues long; the maximum is {max}")]
    TooLong { requested: usize, max: usize },

    #[error("{operation} needs a DNA or RNA sequence, not {alphabet}")]
    NotNucleotide {
        operation: &'static str,
        alphabet: Alphabet,
    },

    #[error("{operation} needs a circular molecule")]
    NotCircular { operation: &'static str },
}
