//! Sequencing reads: quality summaries and trimming of low-quality ends.
//!
//! Qualities are Phred values, `q = -10 log10 P(error)` (Ewing B. and Green P., "Base-calling
//! of automated sequencer traces using phred. II. Error probabilities", Genome Research
//! 8:186–194, 1998).
//!
//! Trimming follows the modified Mott algorithm that phred offers for trimming (`-trim_alt`):
//! every base scores `cutoff - P(error)`, and the read kept is the stretch with the highest
//! total score, found here as a maximum-sum contiguous subarray in one pass (Bentley J.,
//! "Programming pearls: algorithm design techniques", Communications of the ACM 27:865–873,
//! 1984). A base without a quality (typed in by the user) counts as certain.

use serde::Serialize;

/// phred's default error-probability cutoff for trimming.
pub const DEFAULT_CUTOFF: f64 = 0.05;

/// The error probability of a Phred quality.
pub fn error_probability(q: u8) -> f64 {
    10f64.powf(-f64::from(q) / 10.0)
}

/// The best stretch of the read to keep, as `start..end` residue positions; `None` when no
/// base is better than `cutoff`. Ties keep the earliest stretch.
pub fn quality_trim(quality: &[Option<u8>], cutoff: f64) -> Option<std::ops::Range<usize>> {
    let mut best = 0.0f64;
    let mut best_range = None;
    let mut sum = 0.0f64;
    let mut start = 0;
    for (i, q) in quality.iter().enumerate() {
        let score = cutoff - q.map_or(0.0, error_probability);
        if sum + score <= 0.0 {
            sum = 0.0;
            start = i.saturating_add(1);
            continue;
        }
        sum += score;
        if sum > best {
            best = sum;
            best_range = Some(start..i.saturating_add(1));
        }
    }
    best_range
}

/// Quality figures for a read.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct ReadQuality {
    /// Bases with a quality value.
    pub scored: usize,
    /// Mean Phred quality of those bases (0 when there are none).
    pub mean: f64,
    /// Bases of quality 20 or more (at most 1 error in 100).
    pub q20: usize,
    /// Bases of quality 30 or more (at most 1 error in 1,000).
    pub q30: usize,
}

pub fn read_quality(quality: &[Option<u8>]) -> ReadQuality {
    let scored: Vec<u8> = quality.iter().flatten().copied().collect();
    let total: u64 = scored.iter().map(|&q| u64::from(q)).sum();
    ReadQuality {
        scored: scored.len(),
        mean: if scored.is_empty() {
            0.0
        } else {
            total as f64 / scored.len() as f64
        },
        q20: scored.iter().filter(|&&q| q >= 20).count(),
        q30: scored.iter().filter(|&&q| q >= 30).count(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn worked_example() {
        // Scores at cutoff 0.05: q10 -0.05, q40 +0.0499, q5 -0.266. The two q40 bases (1..3)
        // add up to 0.0998; the lone q40 after the q5 base scores only 0.0499.
        let q = [10, 40, 40, 5, 40, 10].map(Some);
        assert_eq!(quality_trim(&q, 0.05), Some(1..3));
    }

    #[test]
    fn edges() {
        assert_eq!(quality_trim(&[], 0.05), None);
        assert_eq!(quality_trim(&[Some(0); 5], 0.05), None);
        assert_eq!(quality_trim(&[Some(30)], 0.05), Some(0..1));
        // Bases without a quality count as certain.
        assert_eq!(
            quality_trim(&[Some(2), None, None, Some(2)], 0.05),
            Some(1..3)
        );
        // A short bad stretch inside a good read is kept.
        let mut q = vec![Some(40); 20];
        q[10] = Some(10);
        assert_eq!(quality_trim(&q, 0.05), Some(0..20));
        // NaN or negative cutoffs keep nothing rather than misbehaving.
        assert_eq!(quality_trim(&q, f64::NAN), None);
        assert_eq!(quality_trim(&q, -1.0), None);
    }

    #[test]
    fn summary() {
        let r = read_quality(&[Some(10), Some(20), Some(30), None]);
        assert_eq!((r.scored, r.q20, r.q30), (3, 2, 1));
        assert!((r.mean - 20.0).abs() < 1e-12);
        assert_eq!(read_quality(&[]).mean, 0.0);
    }
}
