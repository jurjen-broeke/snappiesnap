use std::collections::{BTreeMap, BTreeSet};

use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use snappiesnap_seq::{MAX_LEN, Range, Seq, Topology};

use crate::feature::{MAX_QUALIFIERS, MAX_SEGMENTS, MAX_TEXT};
use crate::{
    DocError, Feature, FeatureId, HistoryNode, NewFeature, Primer, PrimerId, Segment, Trace,
};

/// Longest primer we accept.
const MAX_PRIMER_LEN: usize = 10_000;

/// Free-form information about a molecule.
#[derive(Debug, Clone, Default, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Notes {
    #[serde(default)]
    pub description: String,
    /// Other named fields (author, organism, source, ...), kept in name order.
    #[serde(default)]
    pub fields: BTreeMap<String, String>,
}

/// What an edit did to the annotations, so callers can tell the user.
#[derive(Debug, Clone, Default, PartialEq, Eq, Serialize)]
pub struct EditReport {
    /// Features removed because every residue they covered was deleted.
    pub removed_features: Vec<FeatureId>,
}

/// A molecule with its residues taken out, to keep the rest cheaply (as undo does). It becomes a
/// molecule again only with residues it is valid for ([`Detached::attach`] checks).
#[derive(Debug, Clone, PartialEq)]
pub struct Detached(Molecule);

impl Detached {
    /// The molecule with `seq` as its residues, if every invariant holds.
    pub fn attach(self, seq: Seq) -> Result<Molecule, DocError> {
        let mut m = self.0;
        m.seq = seq;
        m.validate()?;
        Ok(m)
    }
}

/// A DNA, RNA or protein molecule with its annotations.
///
/// Fields are private so every change goes through a method that keeps the document valid
/// (see [`Molecule::validate`] for the invariants).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(try_from = "MoleculeRepr", into = "MoleculeRepr")]
pub struct Molecule {
    name: String,
    seq: Seq,
    topology: Topology,
    features: Vec<Feature>,
    primers: Vec<Primer>,
    notes: Notes,
    next_feature_id: u64,
    next_primer_id: u64,
    /// How the molecule was made (see [`HistoryNode`]).
    history: Option<Box<HistoryNode>>,
    /// The sequencing trace the residues were called from, one call per residue.
    trace: Option<Trace>,
}

impl Molecule {
    pub fn new(name: &str, seq: Seq, topology: Topology) -> Result<Self, DocError> {
        check_name(name)?;
        Ok(Molecule {
            name: name.to_string(),
            seq,
            topology,
            features: Vec::new(),
            primers: Vec::new(),
            notes: Notes::default(),
            next_feature_id: 1,
            next_primer_id: 1,
            history: None,
            trace: None,
        })
    }

    /// How the molecule was made, if that is known.
    pub fn history(&self) -> Option<&HistoryNode> {
        self.history.as_deref()
    }

    /// Replace the construction history (checked against the limits in [`HistoryNode::validate`]).
    pub fn set_history(&mut self, history: Option<HistoryNode>) -> Result<(), DocError> {
        if let Some(h) = &history {
            h.validate()?;
        }
        self.history = history.map(Box::new);
        Ok(())
    }

    /// Remove and return the construction history.
    pub fn take_history(&mut self) -> Option<HistoryNode> {
        self.history.take().map(|h| *h)
    }

    pub fn name(&self) -> &str {
        &self.name
    }

    /// The sequencing trace behind the residues, if the molecule is a read.
    pub fn trace(&self) -> Option<&Trace> {
        self.trace.as_ref()
    }

    /// Attach (or remove) a trace; it needs one call per residue.
    pub fn set_trace(&mut self, trace: Option<Trace>) -> Result<(), DocError> {
        if let Some(t) = &trace {
            check_trace(t, self.len())?;
        }
        self.trace = trace;
        Ok(())
    }

    /// Take the residues out (see [`Detached`]).
    pub fn detach(mut self) -> (Detached, Seq) {
        let alphabet = self.seq.alphabet();
        let seq = std::mem::replace(&mut self.seq, Seq::empty(alphabet));
        (Detached(self), seq)
    }

    pub fn seq(&self) -> &Seq {
        &self.seq
    }

    pub fn len(&self) -> usize {
        self.seq.len()
    }

    pub fn is_empty(&self) -> bool {
        self.seq.is_empty()
    }

    pub fn topology(&self) -> Topology {
        self.topology
    }

    pub fn features(&self) -> &[Feature] {
        &self.features
    }

    pub fn feature(&self, id: FeatureId) -> Option<&Feature> {
        self.features.iter().find(|f| f.id == id)
    }

    pub fn primers(&self) -> &[Primer] {
        &self.primers
    }

    pub fn notes(&self) -> &Notes {
        &self.notes
    }

    pub fn set_notes(&mut self, notes: Notes) -> Result<(), DocError> {
        check_text("description", &notes.description)?;
        for (k, v) in &notes.fields {
            check_text("note field name", k)?;
            check_text("note field", v)?;
        }
        self.notes = notes;
        Ok(())
    }

    pub fn rename(&mut self, name: &str) -> Result<(), DocError> {
        check_name(name)?;
        self.name = name.to_string();
        Ok(())
    }

    /// Residues covered by `range` (across the origin if it wraps).
    pub fn residues(&self, range: Range) -> Result<Vec<u8>, DocError> {
        Ok(self.seq.slice(range, self.topology)?)
    }

    /// Insert residues before position `at`. Features that contain the insertion point grow;
    /// features after it move.
    pub fn insert(&mut self, at: usize, residues: &[u8]) -> Result<EditReport, DocError> {
        self.seq.insert(at, residues)?;
        let n = residues.len();
        if let Some(t) = &mut self.trace {
            t.insert(at, n);
        }
        for seg in self.features.iter_mut().flat_map(|f| f.segments.iter_mut()) {
            seg.range = seg.range.after_insert(at, n);
        }
        Ok(EditReport::default())
    }

    /// Delete the residues in `range`. Features shrink; features left with nothing are removed
    /// and listed in the report.
    pub fn delete(&mut self, range: Range) -> Result<EditReport, DocError> {
        let old_len = self.len();
        let deleted = self.seq.delete(range, self.topology)?;
        if let Some(t) = &mut self.trace {
            // Later piece first, as the residues were removed.
            let (first, second) = deleted.segments(old_len);
            t.delete(first.start, first.end);
            if let Some(second) = second {
                t.delete(second.start, second.end);
            }
        }
        Ok(self.features_after_delete(deleted, old_len))
    }

    /// Replace the residues in `range` with `residues`. All inputs are checked before anything
    /// changes, so a failed replace leaves the molecule untouched.
    pub fn replace(&mut self, range: Range, residues: &[u8]) -> Result<EditReport, DocError> {
        let range = range.validate(self.len(), self.topology)?;
        Seq::new(self.seq.alphabet(), residues)?;
        let after_delete = self.len().saturating_sub(range.len_on(self.len()));
        if after_delete.saturating_add(residues.len()) > MAX_LEN {
            return Err(snappiesnap_seq::SeqError::TooLong {
                requested: after_delete.saturating_add(residues.len()),
                max: MAX_LEN,
            }
            .into());
        }
        // Where the replaced stretch started, once the deletion is done.
        let at = if range.is_wrapping() {
            after_delete
        } else {
            range.start
        };
        let report = self.delete(range)?;
        self.insert(at, residues)?;
        Ok(report)
    }

    /// Reverse-complement the molecule. Features keep covering the same residues, now on the
    /// opposite strand, with their pieces in the new 5'→3' order.
    pub fn reverse_complement(&mut self) -> Result<(), DocError> {
        self.seq.reverse_complement()?;
        if let Some(t) = &mut self.trace {
            t.reverse_complement();
        }
        let len = self.len();
        for f in &mut self.features {
            for seg in &mut f.segments {
                seg.range = seg.range.after_reverse_complement(len);
            }
            f.segments.reverse();
            f.strand = f.strand.flipped();
        }
        Ok(())
    }

    /// Make position `origin` the new position 0 of a circular molecule.
    pub fn rotate_origin(&mut self, origin: usize) -> Result<(), DocError> {
        if self.topology != Topology::Circular {
            return Err(snappiesnap_seq::SeqError::NotCircular {
                operation: "rotating the origin",
            }
            .into());
        }
        let len = self.len();
        self.seq.rotate(origin)?;
        if let Some(t) = &mut self.trace {
            t.rotate(origin);
        }
        for seg in self.features.iter_mut().flat_map(|f| f.segments.iter_mut()) {
            seg.range = seg.range.after_rotate(origin, len);
        }
        Ok(())
    }

    /// Change the topology. A circular molecule can only become linear when no feature spans
    /// the origin.
    pub fn set_topology(&mut self, topology: Topology) -> Result<(), DocError> {
        if topology == Topology::Linear
            && let Some(f) = self
                .features
                .iter()
                .find(|f| f.segments.iter().any(|s| s.range.is_wrapping()))
        {
            return Err(DocError::FeatureSpansOrigin {
                id: f.id,
                name: f.name.clone(),
            });
        }
        self.topology = topology;
        Ok(())
    }

    pub fn add_feature(&mut self, new: NewFeature) -> Result<FeatureId, DocError> {
        let id = FeatureId(self.next_feature_id);
        let feature = self.checked_feature(new.with_id(id))?;
        self.next_feature_id = self.next_feature_id.saturating_add(1);
        self.features.push(feature);
        Ok(id)
    }

    /// Replace the feature with the same id.
    pub fn update_feature(&mut self, feature: Feature) -> Result<(), DocError> {
        let feature = self.checked_feature(feature)?;
        let slot = self
            .features
            .iter_mut()
            .find(|f| f.id == feature.id)
            .ok_or(DocError::FeatureNotFound(feature.id))?;
        *slot = feature;
        Ok(())
    }

    pub fn remove_feature(&mut self, id: FeatureId) -> Result<Feature, DocError> {
        let index = self
            .features
            .iter()
            .position(|f| f.id == id)
            .ok_or(DocError::FeatureNotFound(id))?;
        Ok(self.features.remove(index))
    }

    pub fn add_primer(
        &mut self,
        name: &str,
        seq: Seq,
        description: &str,
    ) -> Result<PrimerId, DocError> {
        let id = PrimerId(self.next_primer_id);
        let primer = checked_primer(Primer {
            id,
            name: name.to_string(),
            seq,
            description: description.to_string(),
        })?;
        self.next_primer_id = self.next_primer_id.saturating_add(1);
        self.primers.push(primer);
        Ok(id)
    }

    /// Replace a primer's name, sequence and description (its id stays).
    pub fn update_primer(&mut self, primer: Primer) -> Result<(), DocError> {
        let primer = checked_primer(primer)?;
        let slot = self
            .primers
            .iter_mut()
            .find(|p| p.id == primer.id)
            .ok_or(DocError::PrimerNotFound(primer.id))?;
        *slot = primer;
        Ok(())
    }

    pub fn remove_primer(&mut self, id: PrimerId) -> Result<Primer, DocError> {
        let index = self
            .primers
            .iter()
            .position(|p| p.id == id)
            .ok_or(DocError::PrimerNotFound(id))?;
        Ok(self.primers.remove(index))
    }

    /// Check every invariant: valid name and notes; every feature has 1..=MAX_SEGMENTS
    /// non-empty segments that are valid on this molecule; ids are unique and below the next
    /// id to hand out; primers are valid.
    pub fn validate(&self) -> Result<(), DocError> {
        check_name(&self.name)?;
        if let Some(h) = &self.history {
            h.validate()?;
        }
        if let Some(t) = &self.trace {
            check_trace(t, self.len())?;
        }
        let mut ids = BTreeSet::new();
        for f in &self.features {
            if self.checked_feature(f.clone())? != *f {
                return Err(DocError::Inconsistent(format!(
                    "feature {} has ranges that aren't normalised",
                    f.id
                )));
            }
            if !ids.insert(f.id.0) || f.id.0 >= self.next_feature_id {
                return Err(DocError::Inconsistent(format!(
                    "feature id {} is duplicated or not below the next id",
                    f.id
                )));
            }
        }
        let mut ids = BTreeSet::new();
        for p in &self.primers {
            checked_primer(p.clone())?;
            if !ids.insert(p.id.0) || p.id.0 >= self.next_primer_id {
                return Err(DocError::Inconsistent(format!(
                    "primer id {} is duplicated or not below the next id",
                    p.id
                )));
            }
        }
        Ok(())
    }

    fn features_after_delete(&mut self, deleted: Range, old_len: usize) -> EditReport {
        let mut report = EditReport::default();
        for f in &mut self.features {
            f.segments
                .retain_mut(|seg| match seg.range.after_delete(deleted, old_len) {
                    Some(r) if !r.is_empty() => {
                        seg.range = r;
                        true
                    }
                    _ => false,
                });
            if f.segments.is_empty() {
                report.removed_features.push(f.id);
            }
        }
        self.features.retain(|f| !f.segments.is_empty());
        report
    }

    /// Validate a feature against this molecule and normalise its ranges.
    fn checked_feature(&self, mut f: Feature) -> Result<Feature, DocError> {
        let bad = |msg: String| DocError::InvalidFeature(msg);
        check_text("feature name", &f.name)?;
        if f.color
            .as_deref()
            .is_some_and(|c| !crate::feature::is_color(c))
        {
            return Err(bad("colour must be written #rrggbb".into()));
        }
        if f.kind.trim().is_empty() || f.kind.chars().count() > MAX_TEXT {
            return Err(bad("type must be non-empty text".into()));
        }
        if f.segments.is_empty() || f.segments.len() > MAX_SEGMENTS {
            return Err(bad(format!(
                "needs between 1 and {MAX_SEGMENTS} segments, got {}",
                f.segments.len()
            )));
        }
        for seg in &mut f.segments {
            let r = seg.range.validate(self.len(), self.topology)?;
            if r.is_empty() {
                return Err(bad(format!(
                    "segment {}..{} is empty",
                    seg.range.start, seg.range.end
                )));
            }
            *seg = Segment { range: r };
        }
        if f.qualifiers.len() > MAX_QUALIFIERS {
            return Err(bad(format!("more than {MAX_QUALIFIERS} qualifiers")));
        }
        for q in &f.qualifiers {
            if q.name.trim().is_empty() || q.values.len() > MAX_QUALIFIERS {
                return Err(bad(
                    "qualifiers need a name and a bounded list of values".into()
                ));
            }
            check_text("qualifier name", &q.name)?;
            for v in &q.values {
                check_text("qualifier value", v)?;
            }
        }
        Ok(f)
    }
}

fn checked_primer(p: Primer) -> Result<Primer, DocError> {
    check_name(&p.name).map_err(|e| DocError::InvalidPrimer(e.to_string()))?;
    check_text("primer description", &p.description)?;
    if !p.seq.alphabet().is_nucleotide() {
        return Err(DocError::InvalidPrimer(
            "a primer must be DNA or RNA".into(),
        ));
    }
    if p.seq.is_empty() || p.seq.len() > MAX_PRIMER_LEN {
        return Err(DocError::InvalidPrimer(format!(
            "length must be 1 to {MAX_PRIMER_LEN} nt, got {}",
            p.seq.len()
        )));
    }
    Ok(p)
}

fn check_trace(t: &Trace, len: usize) -> Result<(), DocError> {
    if t.len() != len {
        return Err(DocError::InvalidTrace(format!(
            "{} calls for {len} residues",
            t.len()
        )));
    }
    check_text("sample name", t.sample_name().unwrap_or(""))
}

fn check_name(name: &str) -> Result<(), DocError> {
    if name.trim().is_empty() {
        return Err(DocError::InvalidName("a name can't be empty".into()));
    }
    if name.chars().any(char::is_control) {
        return Err(DocError::InvalidName(
            "a name can't contain control characters".into(),
        ));
    }
    check_text("name", name)
}

fn check_text(what: &str, text: &str) -> Result<(), DocError> {
    if text.chars().count() > MAX_TEXT {
        return Err(DocError::TextTooLong {
            what: what.to_string(),
            max: MAX_TEXT,
        });
    }
    Ok(())
}

/// The serialised form. `length` is written for readers' convenience and ignored on input.
#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct MoleculeRepr {
    name: String,
    topology: Topology,
    #[serde(default)]
    length: usize,
    seq: Seq,
    #[serde(default)]
    features: Vec<Feature>,
    #[serde(default)]
    primers: Vec<Primer>,
    #[serde(default)]
    notes: Notes,
    next_feature_id: u64,
    next_primer_id: u64,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    history: Option<Box<HistoryNode>>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    trace: Option<Trace>,
}

impl TryFrom<MoleculeRepr> for Molecule {
    type Error = DocError;
    fn try_from(r: MoleculeRepr) -> Result<Self, DocError> {
        let m = Molecule {
            name: r.name,
            seq: r.seq,
            topology: r.topology,
            features: r.features,
            primers: r.primers,
            notes: r.notes,
            next_feature_id: r.next_feature_id,
            next_primer_id: r.next_primer_id,
            history: r.history,
            trace: r.trace,
        };
        m.validate()?;
        Ok(m)
    }
}

impl From<Molecule> for MoleculeRepr {
    fn from(m: Molecule) -> Self {
        MoleculeRepr {
            length: m.seq.len(),
            name: m.name,
            topology: m.topology,
            seq: m.seq,
            features: m.features,
            primers: m.primers,
            notes: m.notes,
            next_feature_id: m.next_feature_id,
            next_primer_id: m.next_primer_id,
            history: m.history,
            trace: m.trace,
        }
    }
}
