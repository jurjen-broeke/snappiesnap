use std::fmt;

use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use snappiesnap_seq::{Range, Seq};

/// Longest feature name, feature type, qualifier name or value we accept (in characters).
pub(crate) const MAX_TEXT: usize = 100_000;
/// Most segments a single feature may have.
pub(crate) const MAX_SEGMENTS: usize = 10_000;
/// Most qualifiers a single feature may have.
pub(crate) const MAX_QUALIFIERS: usize = 10_000;

/// A stable identifier for a feature within one molecule.
#[derive(
    Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize, JsonSchema,
)]
#[serde(transparent)]
pub struct FeatureId(pub u64);

/// A stable identifier for a primer within one molecule.
#[derive(
    Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize, JsonSchema,
)]
#[serde(transparent)]
pub struct PrimerId(pub u64);

impl fmt::Display for FeatureId {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.0.fmt(f)
    }
}

impl fmt::Display for PrimerId {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.0.fmt(f)
    }
}

/// Which strand a feature is on.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Hash, Serialize, Deserialize, JsonSchema)]
#[serde(rename_all = "lowercase")]
pub enum Strand {
    /// Not directional (for example a binding site or a poly(A) signal).
    #[default]
    None,
    Forward,
    Reverse,
    Both,
}

impl Strand {
    /// The strand after the molecule is reverse-complemented.
    pub fn flipped(self) -> Strand {
        match self {
            Strand::Forward => Strand::Reverse,
            Strand::Reverse => Strand::Forward,
            other => other,
        }
    }
}

/// One contiguous piece of a feature (a feature made of several pieces is a join).
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Segment {
    pub range: Range,
}

/// A named list of values attached to a feature (like a GenBank `/note="..."`).
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Qualifier {
    pub name: String,
    #[serde(default)]
    pub values: Vec<String>,
}

/// An annotated region of the molecule.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Feature {
    pub id: FeatureId,
    pub name: String,
    /// The feature type, for example `CDS`, `promoter` or `rep_origin`.
    pub kind: String,
    pub strand: Strand,
    /// One or more pieces, in order along the feature (5' to 3' on its strand).
    pub segments: Vec<Segment>,
    #[serde(default)]
    pub qualifiers: Vec<Qualifier>,
    /// Display colour as `#rrggbb` (kept from `.dna` files; `None`: the default).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub color: Option<String>,
}

/// Whether `c` is a colour as `#rrggbb`.
pub fn is_color(c: &str) -> bool {
    c.len() == 7 && c.starts_with('#') && c.chars().skip(1).all(|ch| ch.is_ascii_hexdigit())
}

/// A feature to add: everything except the id, which the molecule assigns.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct NewFeature {
    pub name: String,
    /// The feature type, for example `CDS`, `promoter` or `rep_origin`.
    #[serde(default = "default_kind")]
    pub kind: String,
    #[serde(default)]
    pub strand: Strand,
    /// One or more pieces, in order along the feature.
    pub segments: Vec<Segment>,
    #[serde(default)]
    pub qualifiers: Vec<Qualifier>,
    /// Display colour as `#rrggbb`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub color: Option<String>,
}

fn default_kind() -> String {
    "misc_feature".to_string()
}

impl NewFeature {
    pub fn with_id(self, id: FeatureId) -> Feature {
        Feature {
            id,
            name: self.name,
            kind: self.kind,
            strand: self.strand,
            segments: self.segments,
            qualifiers: self.qualifiers,
            color: self.color,
        }
    }
}

/// An oligonucleotide kept with the molecule. Binding sites are computed, not stored (M5).
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Primer {
    pub id: PrimerId,
    pub name: String,
    pub seq: Seq,
    #[serde(default)]
    pub description: String,
}
