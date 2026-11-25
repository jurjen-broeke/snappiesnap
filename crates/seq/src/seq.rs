use std::fmt;

use serde::{Deserialize, Serialize};

use crate::{Alphabet, Range, SeqError, Topology};

/// The longest sequence we accept (about 1.07 Gbp). Every length is checked against it, which
/// also keeps coordinate arithmetic far away from overflow on 32-bit (wasm) targets.
pub const MAX_LEN: usize = 1 << 30;

/// A validated residue sequence. Bytes are stored exactly as given (case preserved), and every
/// byte is valid in the sequence's [`Alphabet`].
#[derive(Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[serde(try_from = "SeqRepr", into = "SeqRepr")]
pub struct Seq {
    alphabet: Alphabet,
    residues: Vec<u8>,
}

impl Seq {
    pub fn new(alphabet: Alphabet, residues: impl Into<Vec<u8>>) -> Result<Self, SeqError> {
        let residues = residues.into();
        check_len(residues.len())?;
        validate(alphabet, &residues, 0)?;
        Ok(Seq { alphabet, residues })
    }

    pub fn empty(alphabet: Alphabet) -> Self {
        Seq {
            alphabet,
            residues: Vec::new(),
        }
    }

    pub fn alphabet(&self) -> Alphabet {
        self.alphabet
    }

    pub fn len(&self) -> usize {
        self.residues.len()
    }

    pub fn is_empty(&self) -> bool {
        self.residues.is_empty()
    }

    pub fn as_bytes(&self) -> &[u8] {
        &self.residues
    }

    /// The residues as text. Always succeeds because every valid residue is ASCII.
    pub fn as_str(&self) -> &str {
        std::str::from_utf8(&self.residues).unwrap_or_default()
    }

    /// The residues covered by `range`, in order along the range (across the origin if it
    /// wraps).
    pub fn slice(&self, range: Range, topology: Topology) -> Result<Vec<u8>, SeqError> {
        let range = range.validate(self.len(), topology)?;
        let (first, second) = range.segments(self.len());
        let mut out = Vec::with_capacity(range.len_on(self.len()));
        for piece in std::iter::once(first).chain(second) {
            let len = self.len();
            let bytes = self
                .residues
                .get(piece.clone())
                .ok_or(SeqError::OutOfBounds {
                    position: piece.end,
                    len,
                })?;
            out.extend_from_slice(bytes);
        }
        Ok(out)
    }

    /// Insert `residues` before position `at` (`at == len` appends).
    pub fn insert(&mut self, at: usize, residues: &[u8]) -> Result<(), SeqError> {
        if at > self.len() {
            return Err(SeqError::OutOfBounds {
                position: at,
                len: self.len(),
            });
        }
        check_len(self.len().saturating_add(residues.len()))?;
        validate(self.alphabet, residues, at)?;
        self.residues.splice(at..at, residues.iter().copied());
        Ok(())
    }

    /// Remove the residues covered by `range`. Returns the normalised range that was removed.
    pub fn delete(&mut self, range: Range, topology: Topology) -> Result<Range, SeqError> {
        let range = range.validate(self.len(), topology)?;
        let (first, second) = range.segments(self.len());
        // Remove the later piece first so the earlier piece's positions stay valid.
        self.residues.drain(first);
        if let Some(second) = second {
            self.residues.drain(second);
        }
        Ok(range)
    }

    /// Reverse-complement the whole sequence in place (DNA and RNA only).
    pub fn reverse_complement(&mut self) -> Result<(), SeqError> {
        if !self.alphabet.is_nucleotide() {
            return Err(SeqError::NotNucleotide {
                operation: "reverse complement",
                alphabet: self.alphabet,
            });
        }
        let mut out = Vec::with_capacity(self.len());
        for (position, &b) in self.residues.iter().enumerate().rev() {
            out.push(
                self.alphabet
                    .complement(b)
                    .ok_or(SeqError::InvalidResidue {
                        alphabet: self.alphabet,
                        position,
                        found: char::from(b),
                    })?,
            );
        }
        self.residues = out;
        Ok(())
    }

    /// Rotate so that position `origin` becomes position 0.
    pub fn rotate(&mut self, origin: usize) -> Result<(), SeqError> {
        if origin > self.len() {
            return Err(SeqError::OutOfBounds {
                position: origin,
                len: self.len(),
            });
        }
        if origin < self.len() {
            self.residues.rotate_left(origin);
        }
        Ok(())
    }
}

/// How one sequence turns into another: the stretch `start..start + removed` is replaced by
/// `residues` (and the alphabet becomes `alphabet`). Small when the sequences share most of their
/// residues, as after an edit: the common start and end are left out.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SeqDelta {
    alphabet: Alphabet,
    start: usize,
    removed: usize,
    residues: Vec<u8>,
}

impl SeqDelta {
    /// Bytes of residues it keeps.
    pub fn stored(&self) -> usize {
        self.residues.len()
    }
}

impl Seq {
    /// What turns `self` into `target` (see [`SeqDelta`]).
    pub fn delta_to(&self, target: &Seq) -> SeqDelta {
        let (a, b) = (self.as_bytes(), target.as_bytes());
        let prefix = a.iter().zip(b).take_while(|(x, y)| x == y).count();
        let room = a.len().min(b.len()) - prefix;
        let suffix = a
            .iter()
            .rev()
            .zip(b.iter().rev())
            .take(room)
            .take_while(|(x, y)| x == y)
            .count();
        SeqDelta {
            alphabet: target.alphabet,
            start: prefix,
            removed: a.len() - prefix - suffix,
            residues: b.get(prefix..b.len() - suffix).unwrap_or_default().to_vec(),
        }
    }

    /// `self` changed by `delta` (made by [`Seq::delta_to`] from a sequence equal to `self`).
    pub fn apply(&self, delta: &SeqDelta) -> Result<Seq, SeqError> {
        let end = delta
            .start
            .checked_add(delta.removed)
            .filter(|&e| e <= self.len())
            .ok_or(SeqError::OutOfBounds {
                position: delta.start.saturating_add(delta.removed),
                len: self.len(),
            })?;
        let mut out =
            Vec::with_capacity((self.len() - delta.removed).saturating_add(delta.residues.len()));
        out.extend_from_slice(self.residues.get(..delta.start).unwrap_or_default());
        out.extend_from_slice(&delta.residues);
        out.extend_from_slice(self.residues.get(end..).unwrap_or_default());
        Seq::new(delta.alphabet, out)
    }
}

fn check_len(requested: usize) -> Result<(), SeqError> {
    if requested > MAX_LEN {
        Err(SeqError::TooLong {
            requested,
            max: MAX_LEN,
        })
    } else {
        Ok(())
    }
}

/// `offset` is added to reported positions (where the bytes will end up in the sequence).
fn validate(alphabet: Alphabet, bytes: &[u8], offset: usize) -> Result<(), SeqError> {
    match bytes.iter().position(|&b| !alphabet.is_valid(b)) {
        None => Ok(()),
        Some(i) => Err(SeqError::InvalidResidue {
            alphabet,
            position: offset.saturating_add(i),
            found: bytes.get(i).map_or('?', |&b| char::from(b)),
        }),
    }
}

impl fmt::Debug for Seq {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        const PREVIEW: usize = 24;
        let preview = self.residues.get(..PREVIEW).unwrap_or(&self.residues);
        write!(
            f,
            "Seq({}, {} residues, {}{})",
            self.alphabet,
            self.len(),
            String::from_utf8_lossy(preview),
            if self.len() > PREVIEW { "…" } else { "" }
        )
    }
}

/// The serialised form: `{"alphabet": "dna", "residues": "ACGT"}`.
#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct SeqRepr {
    alphabet: Alphabet,
    residues: String,
}

impl TryFrom<SeqRepr> for Seq {
    type Error = SeqError;
    fn try_from(r: SeqRepr) -> Result<Self, SeqError> {
        Seq::new(r.alphabet, r.residues.into_bytes())
    }
}

impl From<Seq> for SeqRepr {
    fn from(s: Seq) -> Self {
        SeqRepr {
            alphabet: s.alphabet,
            residues: String::from_utf8_lossy(&s.residues).into_owned(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use Alphabet::*;
    use Topology::*;

    fn dna(s: &str) -> Seq {
        Seq::new(Dna, s).unwrap()
    }

    #[test]
    fn rejects_invalid_residues_with_their_position() {
        let err = Seq::new(Dna, "ACGXT").unwrap_err();
        assert_eq!(
            err,
            SeqError::InvalidResidue {
                alphabet: Dna,
                position: 3,
                found: 'X'
            }
        );
        let mut s = dna("AC");
        assert!(matches!(
            s.insert(1, b"G-"),
            Err(SeqError::InvalidResidue { position: 2, .. })
        ));
        assert_eq!(s.as_str(), "AC");
    }

    #[test]
    fn slice_across_the_origin() {
        let s = dna("AACCGGTTAC");
        assert_eq!(s.slice(Range::new(8, 2), Circular).unwrap(), b"ACAA");
        assert!(s.slice(Range::new(8, 2), Linear).is_err());
        assert_eq!(s.slice(Range::at(3), Linear).unwrap(), b"");
    }

    #[test]
    fn edits() {
        let mut s = dna("AACCGGTTAC");
        s.delete(Range::new(8, 2), Circular).unwrap();
        assert_eq!(s.as_str(), "CCGGTT");
        s.insert(6, b"aa").unwrap();
        assert_eq!(s.as_str(), "CCGGTTaa");
        s.reverse_complement().unwrap();
        assert_eq!(s.as_str(), "ttAACCGG");
        s.rotate(2).unwrap();
        assert_eq!(s.as_str(), "AACCGGtt");
        assert!(s.rotate(9).is_err());
        s.rotate(8).unwrap();
        assert_eq!(s.as_str(), "AACCGGtt");
    }

    #[test]
    fn tiny_and_empty_sequences() {
        let mut e = Seq::empty(Dna);
        assert_eq!(e.slice(Range::full(0), Linear).unwrap(), b"");
        e.reverse_complement().unwrap();
        e.rotate(0).unwrap();
        assert!(e.delete(Range::new(0, 1), Circular).is_err());
        let mut one = dna("G");
        one.reverse_complement().unwrap();
        assert_eq!(one.as_str(), "C");
        one.delete(Range::full(1), Circular).unwrap();
        assert!(one.is_empty());
    }

    #[test]
    fn protein_cannot_be_reverse_complemented() {
        let mut p = Seq::new(Protein, "MKV*").unwrap();
        assert!(matches!(
            p.reverse_complement(),
            Err(SeqError::NotNucleotide { .. })
        ));
    }

    #[test]
    fn serde_round_trip_validates() {
        let s = dna("ACGT");
        let json = serde_json::to_string(&s).unwrap();
        assert_eq!(json, r#"{"alphabet":"dna","residues":"ACGT"}"#);
        assert_eq!(serde_json::from_str::<Seq>(&json).unwrap(), s);
        assert!(serde_json::from_str::<Seq>(r#"{"alphabet":"dna","residues":"ACGU"}"#).is_err());
    }

    #[test]
    fn deltas_are_small_and_exact() {
        let dna = |s: &str| Seq::new(Alphabet::Dna, s).unwrap();
        for (a, b) in [
            ("ACGTACGT", "ACGTTACGT"),
            ("ACGTACGT", "ACGACGT"),
            ("ACGTACGT", "ACGTACGT"),
            ("ACGTACGT", ""),
            ("", "ACGT"),
            ("AAAA", "AAAAAA"),
            ("ACGT", "TGCA"),
        ] {
            let (x, y) = (dna(a), dna(b));
            let d = x.delta_to(&y);
            assert_eq!(x.apply(&d).unwrap(), y, "{a} -> {b}");
            let back = y.delta_to(&x);
            assert_eq!(y.apply(&back).unwrap(), x, "{b} -> {a}");
        }
        // One inserted base in a long sequence keeps one base.
        let long = dna(&"ACGT".repeat(10_000));
        let mut edited = long.clone();
        edited.insert(20_000, b"G").unwrap();
        assert_eq!(long.delta_to(&edited).stored(), 1);
        assert_eq!(edited.delta_to(&long).stored(), 0);
        // A delta that doesn't fit is an error, not a panic.
        assert!(dna("AC").apply(&long.delta_to(&edited)).is_err());
    }
}
