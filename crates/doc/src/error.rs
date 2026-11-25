use snappiesnap_seq::SeqError;

use crate::{FeatureId, PrimerId};

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum DocError {
    #[error(transparent)]
    Seq(#[from] SeqError),

    #[error("no feature with id {0}")]
    FeatureNotFound(FeatureId),

    #[error("no primer with id {0}")]
    PrimerNotFound(PrimerId),

    #[error("invalid feature: {0}")]
    InvalidFeature(String),

    #[error("invalid primer: {0}")]
    InvalidPrimer(String),

    #[error("invalid trace: {0}")]
    InvalidTrace(String),

    #[error("invalid name: {0}")]
    InvalidName(String),

    #[error("{what} is longer than {max} characters")]
    TextTooLong { what: String, max: usize },

    #[error(
        "cannot make the molecule linear: feature {name:?} (id {id}) spans the origin; rotate the origin or edit the feature first"
    )]
    FeatureSpansOrigin { id: FeatureId, name: String },

    #[error("invalid document: {0}")]
    Inconsistent(String),
}
