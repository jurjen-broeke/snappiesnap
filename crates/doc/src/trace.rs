//! The sequencing trace (chromatogram) behind a molecule read from a trace file.
//!
//! The signal never changes, so it is shared between copies (undo steps, history) behind an
//! [`Arc`]. What changes with edits is the per-residue list of calls: each residue keeps the
//! peak and quality it was called with, and residues typed in by the user have neither.
//! Reverse-complementing flips a flag instead of rewriting the signal.

use std::sync::Arc;

use serde::{Deserialize, Serialize};

use crate::DocError;

/// Most sample points per channel we accept (traces have about 10–30 thousand).
pub const MAX_SAMPLES: usize = 2_000_000;

/// The four signal channels, one value per sample point, as read (not reversed).
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Signal {
    pub a: Vec<u16>,
    pub c: Vec<u16>,
    pub g: Vec<u16>,
    pub t: Vec<u16>,
}

impl Signal {
    fn samples(&self) -> usize {
        self.a.len()
    }
}

/// What the instrument called at one residue.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct Call {
    /// The sample point of the residue's peak in the signal as read; `None` for residues that
    /// weren't called (typed in by the user).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub peak: Option<u32>,
    /// Phred quality, `-10 log10 P(error)`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub quality: Option<u8>,
}

/// A chromatogram aligned residue by residue to its molecule.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(try_from = "TraceRepr", into = "TraceRepr")]
pub struct Trace {
    signal: Arc<Signal>,
    /// One per residue of the molecule, in the molecule's current order.
    calls: Vec<Call>,
    /// The molecule was reverse-complemented: sample `i` shows the signal as read at
    /// `samples - 1 - i`, with A and T, C and G swapped.
    reversed: bool,
    /// The sample name the instrument recorded.
    sample_name: Option<String>,
}

/// A channel of the signal, by the base it reports.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Channel {
    A,
    C,
    G,
    T,
}

impl Channel {
    pub const ALL: [Channel; 4] = [Channel::A, Channel::C, Channel::G, Channel::T];

    pub fn base(self) -> u8 {
        match self {
            Channel::A => b'A',
            Channel::C => b'C',
            Channel::G => b'G',
            Channel::T => b'T',
        }
    }

    fn complement(self) -> Channel {
        match self {
            Channel::A => Channel::T,
            Channel::C => Channel::G,
            Channel::G => Channel::C,
            Channel::T => Channel::A,
        }
    }
}

impl Trace {
    /// A trace with `calls` (one per residue). Peaks outside the signal are dropped (counted
    /// in the returned number) rather than pointing nowhere.
    pub fn new(
        signal: Signal,
        calls: Vec<Call>,
        sample_name: Option<String>,
    ) -> Result<(Trace, usize), DocError> {
        let samples = signal.samples();
        if [&signal.c, &signal.g, &signal.t]
            .iter()
            .any(|c| c.len() != samples)
        {
            return Err(DocError::InvalidTrace(
                "the four channels differ in length".into(),
            ));
        }
        if samples > MAX_SAMPLES {
            return Err(DocError::InvalidTrace(format!(
                "{samples} sample points; at most {MAX_SAMPLES} are supported"
            )));
        }
        let mut dropped = 0;
        let calls = calls
            .into_iter()
            .map(|mut c| {
                if c.peak
                    .is_some_and(|p| usize::try_from(p).map_or(true, |p| p >= samples))
                {
                    c.peak = None;
                    dropped += 1;
                }
                c
            })
            .collect();
        let trace = Trace {
            signal: Arc::new(signal),
            calls,
            reversed: false,
            sample_name,
        };
        Ok((trace, dropped))
    }

    /// Sample points per channel.
    pub fn samples(&self) -> usize {
        self.signal.samples()
    }

    pub fn calls(&self) -> &[Call] {
        &self.calls
    }

    pub fn sample_name(&self) -> Option<&str> {
        self.sample_name.as_deref()
    }

    /// Whether the molecule is the reverse complement of the read.
    pub fn is_reversed(&self) -> bool {
        self.reversed
    }

    /// The signal of `channel` at sample point `i`, in the molecule's orientation.
    pub fn value(&self, channel: Channel, i: usize) -> Option<u16> {
        let (channel, i) = if self.reversed {
            (
                channel.complement(),
                self.samples().checked_sub(1)?.checked_sub(i)?,
            )
        } else {
            (channel, i)
        };
        let data = match channel {
            Channel::A => &self.signal.a,
            Channel::C => &self.signal.c,
            Channel::G => &self.signal.g,
            Channel::T => &self.signal.t,
        };
        data.get(i).copied()
    }

    /// The peak of residue `i`, as a sample point in the molecule's orientation.
    pub fn peak(&self, i: usize) -> Option<usize> {
        let p = usize::try_from(self.calls.get(i)?.peak?).ok()?;
        if self.reversed {
            self.samples().checked_sub(1)?.checked_sub(p)
        } else {
            Some(p)
        }
    }

    /// The highest value in any channel (for scaling a drawing), at least 1.
    pub fn max_value(&self) -> u16 {
        let s = &self.signal;
        [&s.a, &s.c, &s.g, &s.t]
            .iter()
            .flat_map(|c| c.iter())
            .copied()
            .max()
            .unwrap_or(0)
            .max(1)
    }

    pub(crate) fn insert(&mut self, at: usize, n: usize) {
        let at = at.min(self.calls.len());
        self.calls
            .splice(at..at, std::iter::repeat_n(Call::default(), n));
    }

    /// Remove the calls of residues `start..end` (not wrapping).
    pub(crate) fn delete(&mut self, start: usize, end: usize) {
        let end = end.min(self.calls.len());
        if start < end {
            self.calls.drain(start..end);
        }
    }

    pub(crate) fn reverse_complement(&mut self) {
        self.calls.reverse();
        self.reversed = !self.reversed;
    }

    pub(crate) fn rotate(&mut self, origin: usize) {
        if origin < self.calls.len() {
            self.calls.rotate_left(origin);
        }
    }

    pub(crate) fn len(&self) -> usize {
        self.calls.len()
    }
}

#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct TraceRepr {
    signal: Signal,
    calls: Vec<Call>,
    #[serde(default)]
    reversed: bool,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    sample_name: Option<String>,
}

impl TryFrom<TraceRepr> for Trace {
    type Error = DocError;
    fn try_from(r: TraceRepr) -> Result<Self, DocError> {
        let (mut t, dropped) = Trace::new(r.signal, r.calls, r.sample_name)?;
        if dropped > 0 {
            return Err(DocError::InvalidTrace(format!(
                "{dropped} peak(s) lie outside the signal"
            )));
        }
        t.reversed = r.reversed;
        Ok(t)
    }
}

impl From<Trace> for TraceRepr {
    fn from(t: Trace) -> Self {
        TraceRepr {
            signal: Arc::unwrap_or_clone(t.signal),
            calls: t.calls,
            reversed: t.reversed,
            sample_name: t.sample_name,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Molecule, Range, Seq, Topology};

    /// A read of `residues` whose residue `i` peaks at sample `10 i + 5` with quality `i`.
    fn read(residues: &str, topology: Topology) -> Molecule {
        let n = residues.len();
        let samples = 10 * n + 10;
        let ramp: Vec<u16> = (0..samples as u16).collect();
        let signal = Signal {
            a: ramp.clone(),
            c: vec![1; samples],
            g: vec![2; samples],
            t: vec![3; samples],
        };
        let calls = (0..n)
            .map(|i| Call {
                peak: Some(10 * i as u32 + 5),
                quality: Some(i as u8),
            })
            .collect();
        let (t, dropped) = Trace::new(signal, calls, Some("r1".into())).unwrap();
        assert_eq!(dropped, 0);
        let mut m = Molecule::new(
            "r",
            Seq::new(crate::Alphabet::Dna, residues).unwrap(),
            topology,
        )
        .unwrap();
        m.set_trace(Some(t)).unwrap();
        m
    }

    fn qualities(m: &Molecule) -> Vec<Option<u8>> {
        m.trace()
            .unwrap()
            .calls()
            .iter()
            .map(|c| c.quality)
            .collect()
    }

    #[test]
    fn calls_follow_insert_delete_and_replace() {
        let mut m = read("ACGTAC", Topology::Linear);
        m.insert(2, b"GG").unwrap();
        assert_eq!(
            qualities(&m),
            [
                Some(0),
                Some(1),
                None,
                None,
                Some(2),
                Some(3),
                Some(4),
                Some(5)
            ]
        );
        m.delete(Range::new(0, 3)).unwrap();
        assert_eq!(qualities(&m), [None, Some(2), Some(3), Some(4), Some(5)]);
        m.replace(Range::new(1, 2), b"T").unwrap();
        assert_eq!(qualities(&m), [None, None, Some(3), Some(4), Some(5)]);
        // Zero-length delete and insert change nothing.
        m.delete(Range::new(2, 2)).unwrap();
        m.insert(5, b"").unwrap();
        assert_eq!(m.trace().unwrap().calls().len(), 5);
        m.validate().unwrap();
    }

    #[test]
    fn wrapping_and_full_deletes_on_circular_reads() {
        let mut m = read("ACGTAC", Topology::Circular);
        m.delete(Range::new(4, 2)).unwrap();
        assert_eq!(qualities(&m), [Some(2), Some(3)]);
        let mut m = read("ACGTAC", Topology::Circular);
        m.delete(Range::new(0, 6)).unwrap();
        assert!(m.trace().unwrap().calls().is_empty());
        m.validate().unwrap();
        let mut m = read("ACGTAC", Topology::Circular);
        m.rotate_origin(4).unwrap();
        assert_eq!(qualities(&m)[0], Some(4));
        m.validate().unwrap();
    }

    #[test]
    fn reverse_complement_mirrors_the_signal() {
        let mut m = read("AACG", Topology::Linear);
        let samples = m.trace().unwrap().samples();
        m.reverse_complement().unwrap();
        let t = m.trace().unwrap();
        assert!(t.is_reversed());
        // Residue 0 is now the complement of the old last residue (peak 35).
        assert_eq!(t.peak(0), Some(samples - 1 - 35));
        assert_eq!(t.calls()[0].quality, Some(3));
        // The T channel shows the old A channel, mirrored.
        assert_eq!(t.value(Channel::T, 0), Some((samples - 1) as u16));
        assert_eq!(t.value(Channel::A, 0), Some(3));
        assert_eq!(t.value(Channel::A, samples), None);
        m.reverse_complement().unwrap();
        assert_eq!(m, read("AACG", Topology::Linear));
    }

    #[test]
    fn short_reads() {
        let mut m = read("A", Topology::Linear);
        m.delete(Range::new(0, 1)).unwrap();
        assert_eq!(m.len(), 0);
        m.insert(0, b"C").unwrap();
        assert_eq!(qualities(&m), [None]);
        m.reverse_complement().unwrap();
        m.validate().unwrap();
        let empty = read("", Topology::Linear);
        assert_eq!(empty.trace().unwrap().calls().len(), 0);
    }

    #[test]
    fn mismatched_and_hostile_traces_are_rejected() {
        let mut m = Molecule::new(
            "x",
            Seq::new(crate::Alphabet::Dna, "ACG").unwrap(),
            Topology::Linear,
        )
        .unwrap();
        let other = read("AC", Topology::Linear).trace().unwrap().clone();
        assert!(matches!(
            m.set_trace(Some(other)),
            Err(DocError::InvalidTrace(_))
        ));
        let uneven = Signal {
            a: vec![0; 3],
            c: vec![0; 2],
            g: vec![],
            t: vec![],
        };
        assert!(Trace::new(uneven, vec![], None).is_err());
        let signal = Signal {
            a: vec![0; 3],
            c: vec![0; 3],
            g: vec![0; 3],
            t: vec![0; 3],
        };
        let calls = vec![
            Call {
                peak: Some(3),
                quality: None,
            },
            Call {
                peak: Some(u32::MAX),
                quality: None,
            },
        ];
        let (t, dropped) = Trace::new(signal, calls, None).unwrap();
        assert_eq!(dropped, 2);
        assert!(t.calls().iter().all(|c| c.peak.is_none()));
    }

    #[test]
    fn serde_and_undo_keep_the_trace() {
        let mut m = read("ACGT", Topology::Linear);
        m.reverse_complement().unwrap();
        let json = serde_json::to_string(&m).unwrap();
        let back: Molecule = serde_json::from_str(&json).unwrap();
        assert_eq!(back, m);
        // A trace whose calls don't match the residues doesn't load.
        let bad = json.replacen("\"calls\":[{", "\"calls\":[{\"peak\":1},{", 1);
        assert!(serde_json::from_str::<Molecule>(&bad).is_err());
        // Undo keeps the molecule without residues; the signal is shared, not copied.
        let (detached, seq) = m.clone().detach();
        let again = detached.attach(seq).unwrap();
        assert!(Arc::ptr_eq(
            &again.trace().unwrap().signal,
            &m.trace().unwrap().signal
        ));
    }
}
