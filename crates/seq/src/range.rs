//! Ranges on linear and circular molecules.
//!
//! A [`Range`] is half-open, `[start, end)`, in 0-based residue positions. On a circular
//! molecule `start > end` means the range wraps past the origin: it covers `start..len` and
//! then `0..end`. Every edit (insert, delete, rotate, reverse-complement) has a matching
//! `after_*` method that moves a range the way an annotation on that stretch should move.

use schemars::JsonSchema;
use serde::{Deserialize, Serialize};

use crate::SeqError;

/// Whether a molecule's ends are joined.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Hash, Serialize, Deserialize, JsonSchema)]
#[serde(rename_all = "lowercase")]
pub enum Topology {
    #[default]
    Linear,
    Circular,
}

/// A half-open range `[start, end)` of 0-based positions. See the module docs for wrapping.
///
/// The fields are public so ranges can be built from parameters and files, but a range only
/// means something relative to a molecule: call [`Range::validate`] before using one that came
/// from outside.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Hash, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Range {
    /// First position (0-based, inclusive).
    pub start: usize,
    /// One past the last position (exclusive). Smaller than `start` for a range that wraps
    /// past the origin of a circular molecule.
    pub end: usize,
}

impl Range {
    pub const fn new(start: usize, end: usize) -> Self {
        Range { start, end }
    }

    /// The whole molecule.
    pub const fn full(len: usize) -> Self {
        Range { start: 0, end: len }
    }

    /// An empty range at `pos` (a cursor between two residues).
    pub const fn at(pos: usize) -> Self {
        Range {
            start: pos,
            end: pos,
        }
    }

    pub fn is_empty(self) -> bool {
        self.start == self.end
    }

    pub fn is_wrapping(self) -> bool {
        self.start > self.end
    }

    /// Check this range against a molecule of length `len` and normalise it.
    ///
    /// A wrapping range that ends at 0 becomes `start..len`, and one that starts at `len`
    /// becomes `0..end`, so every valid wrapping range has `end > 0` and `start < len`.
    pub fn validate(self, len: usize, topology: Topology) -> Result<Range, SeqError> {
        for position in [self.start, self.end] {
            if position > len {
                return Err(SeqError::OutOfBounds { position, len });
            }
        }
        if !self.is_wrapping() {
            return Ok(self);
        }
        if self.end == 0 {
            return Ok(Range::new(self.start, len));
        }
        if topology != Topology::Circular {
            return Err(SeqError::WrapOnLinear {
                start: self.start,
                end: self.end,
            });
        }
        if self.start == len {
            return Ok(Range::new(0, self.end));
        }
        Ok(self)
    }

    /// Number of residues covered, on a molecule of length `mol_len`.
    pub fn len_on(self, mol_len: usize) -> usize {
        if self.is_wrapping() {
            mol_len.saturating_sub(self.start).saturating_add(self.end)
        } else {
            self.end.saturating_sub(self.start)
        }
    }

    /// The range as one or two linear pieces, in order along the range.
    pub fn segments(
        self,
        mol_len: usize,
    ) -> (std::ops::Range<usize>, Option<std::ops::Range<usize>>) {
        if self.is_wrapping() {
            (self.start..mol_len, Some(0..self.end))
        } else {
            (self.start..self.end, None)
        }
    }

    /// Whether the residue at `pos` is inside this range.
    pub fn contains(self, pos: usize, mol_len: usize) -> bool {
        let (a, b) = self.segments(mol_len);
        a.contains(&pos) || b.is_some_and(|b| b.contains(&pos))
    }

    /// Where this range ends up after `n` residues are inserted at position `at`.
    ///
    /// An insertion strictly inside the range extends it; one at either boundary stays
    /// outside it. On a circular molecule the origin is inside a range that wraps, so
    /// inserting at `0` or at `len` extends such a range.
    pub fn after_insert(self, at: usize, n: usize) -> Range {
        if self.is_empty() {
            let p = if self.start >= at {
                self.start.saturating_add(n)
            } else {
                self.start
            };
            return Range::at(p);
        }
        let start = if self.start >= at {
            self.start.saturating_add(n)
        } else {
            self.start
        };
        let end = if self.end > at {
            self.end.saturating_add(n)
        } else {
            self.end
        };
        Range::new(start, end)
    }

    /// Where this range ends up after `deleted` is removed from a molecule of length
    /// `mol_len`. Both ranges must already be validated against that molecule.
    ///
    /// Returns `None` when every residue of a non-empty range is deleted. An empty range (a
    /// cursor) inside the deletion moves to where the deletion was.
    pub fn after_delete(self, deleted: Range, mol_len: usize) -> Option<Range> {
        if deleted.is_empty() || mol_len == 0 {
            return Some(self);
        }
        let del = Deletion::new(deleted, mol_len);
        if self.is_empty() {
            return Some(Range::at(del.map(self.start)));
        }

        let span = self.len_on(mol_len);
        // Distance from the range start, following the range (and the origin, if it wraps).
        let dist = |p: usize| {
            if p >= self.start {
                p - self.start
            } else {
                mol_len.saturating_sub(self.start).saturating_add(p)
            }
        };

        // First surviving residue, walking forwards from the start. A cyclic deletion is at
        // most two pieces, so a few jumps always suffice; the bound keeps this finite.
        let mut first = self.start;
        for _ in 0..4 {
            match del.piece_containing(first) {
                Some(piece) if piece.end >= mol_len => first = 0,
                Some(piece) => first = piece.end,
                None => break,
            }
        }
        // Last surviving residue, walking backwards from the end.
        let mut last = self.end.checked_sub(1).unwrap_or(mol_len.saturating_sub(1));
        for _ in 0..4 {
            match del.piece_containing(last) {
                Some(piece) if piece.start == 0 => last = mol_len.saturating_sub(1),
                Some(piece) => last = piece.start.saturating_sub(1),
                None => break,
            }
        }
        if del.contains(first) || del.contains(last) {
            return None;
        }
        if dist(first) >= span || dist(last) >= span || dist(last) < dist(first) {
            return None;
        }

        let start = del.map(first);
        let end = del.map(last).saturating_add(1);
        if start == end {
            // The range wrapped all the way round what is left of the molecule.
            return Some(Range::full(mol_len.saturating_sub(deleted.len_on(mol_len))));
        }
        Some(Range::new(start, end))
    }

    /// Where this range ends up when a circular molecule of length `mol_len` is rotated so
    /// that position `origin` becomes position 0. A range covering the whole molecule stays
    /// `0..len`.
    pub fn after_rotate(self, origin: usize, mol_len: usize) -> Range {
        if mol_len == 0 {
            return self;
        }
        let shift = |p: usize| {
            let p = p % mol_len;
            let o = origin % mol_len;
            if p >= o { p - o } else { mol_len - o + p }
        };
        if self.is_empty() {
            return Range::at(shift(self.start));
        }
        if self.len_on(mol_len) >= mol_len {
            return Range::full(mol_len);
        }
        let start = shift(self.start);
        let end = match shift(self.end) {
            0 => mol_len,
            e => e,
        };
        Range::new(start, end)
    }

    /// Where this range ends up when a molecule of length `mol_len` is reverse-complemented
    /// (position `p` becomes `mol_len - 1 - p`).
    pub fn after_reverse_complement(self, mol_len: usize) -> Range {
        if self.is_empty() {
            return Range::at(mol_len.saturating_sub(self.start));
        }
        Range::new(
            mol_len.saturating_sub(self.end),
            mol_len.saturating_sub(self.start),
        )
    }
}

/// A validated, non-empty deletion as one or two linear pieces.
struct Deletion {
    first: std::ops::Range<usize>,
    second: Option<std::ops::Range<usize>>,
}

impl Deletion {
    fn new(range: Range, mol_len: usize) -> Self {
        let (first, second) = range.segments(mol_len);
        Deletion { first, second }
    }

    fn pieces(&self) -> impl Iterator<Item = &std::ops::Range<usize>> {
        std::iter::once(&self.first).chain(self.second.as_ref())
    }

    fn contains(&self, p: usize) -> bool {
        self.pieces().any(|piece| piece.contains(&p))
    }

    fn piece_containing(&self, p: usize) -> Option<std::ops::Range<usize>> {
        self.pieces().find(|piece| piece.contains(&p)).cloned()
    }

    /// The new position of `p` once the deletion is removed.
    fn map(&self, p: usize) -> usize {
        let removed_before: usize = self
            .pieces()
            .map(|piece| p.min(piece.end).saturating_sub(piece.start))
            .sum();
        p.saturating_sub(removed_before)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn validate_normalises_and_rejects() {
        let c = Topology::Circular;
        let l = Topology::Linear;
        assert_eq!(Range::new(8, 2).validate(10, c), Ok(Range::new(8, 2)));
        assert_eq!(Range::new(8, 0).validate(10, l), Ok(Range::new(8, 10)));
        assert_eq!(Range::new(10, 3).validate(10, c), Ok(Range::new(0, 3)));
        assert!(matches!(
            Range::new(8, 2).validate(10, l),
            Err(SeqError::WrapOnLinear { .. })
        ));
        assert!(matches!(
            Range::new(0, 11).validate(10, l),
            Err(SeqError::OutOfBounds { .. })
        ));
        assert_eq!(Range::new(0, 0).validate(0, l), Ok(Range::at(0)));
    }

    #[test]
    fn wrapping_basics() {
        let r = Range::new(8, 2);
        assert_eq!(r.len_on(10), 4);
        assert!(r.contains(9, 10) && r.contains(0, 10) && !r.contains(5, 10));
        assert_eq!(r.segments(10), (8..10, Some(0..2)));
    }

    #[test]
    fn insert_at_the_origin_extends_a_wrapping_range() {
        assert_eq!(Range::new(8, 2).after_insert(0, 3), Range::new(11, 5));
        assert_eq!(Range::new(8, 2).after_insert(10, 3), Range::new(8, 2));
        assert_eq!(Range::new(2, 5).after_insert(2, 3), Range::new(5, 8));
        assert_eq!(Range::new(2, 5).after_insert(5, 3), Range::new(2, 5));
        assert_eq!(Range::new(2, 5).after_insert(3, 3), Range::new(2, 8));
    }

    #[test]
    fn delete_cases() {
        // Inside, overlapping each end, covering, and a deletion that wraps.
        assert_eq!(
            Range::new(2, 8).after_delete(Range::new(4, 6), 10),
            Some(Range::new(2, 6))
        );
        assert_eq!(
            Range::new(2, 8).after_delete(Range::new(0, 4), 10),
            Some(Range::new(0, 4))
        );
        assert_eq!(Range::new(2, 8).after_delete(Range::new(1, 9), 10), None);
        assert_eq!(
            Range::new(8, 2).after_delete(Range::new(9, 1), 10),
            Some(Range::new(7, 1))
        );
        assert_eq!(
            Range::full(10).after_delete(Range::new(3, 5), 10),
            Some(Range::full(8))
        );
    }

    #[test]
    fn rotate_and_reverse_complement() {
        assert_eq!(Range::new(2, 5).after_rotate(4, 10), Range::new(8, 1));
        assert_eq!(Range::new(8, 2).after_rotate(8, 10), Range::new(0, 4));
        assert_eq!(Range::new(6, 10).after_rotate(6, 10), Range::new(0, 4));
        assert_eq!(Range::full(10).after_rotate(3, 10), Range::full(10));
        assert_eq!(
            Range::new(2, 5).after_reverse_complement(10),
            Range::new(5, 8)
        );
        assert_eq!(
            Range::new(8, 2).after_reverse_complement(10),
            Range::new(8, 2)
        );
    }
}
