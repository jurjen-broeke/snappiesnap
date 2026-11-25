//! Secondary structure of primers: hairpins and dimers, by free energy.
//!
//! Stacks use the nearest-neighbour ΔH/ΔS of [`crate::tm`] (SantaLucia 1998) with its salt
//! correction, evaluated at 37 °C: ΔG = ΔH − T·ΔS. Duplexes between two strands (dimers) add the
//! terminal-pair initiation terms of the same paper. Hairpins add a loop penalty instead (DNA
//! hairpin loop ΔG°37 of SantaLucia J., Hicks D. (2004) "The thermodynamics of DNA structural
//! motifs", Annu. Rev. Biophys. Biomol. Struct. 33:415-440, Table 4, with its logarithmic
//! extrapolation for longer loops).
//!
//! Only ungapped, perfectly paired stretches are considered (no internal mismatches, bulges or
//! loops inside a stem), so values are those of the most stable simple structure: a fast
//! screening estimate in the spirit of the checks oligo suppliers recommend, not a full
//! folding prediction.

use serde::Serialize;

use crate::tm::TmConditions;

/// 37 °C in kelvin.
const T37: f64 = 310.15;
/// Gas constant, kcal/(K·mol).
const R_KCAL: f64 = 0.001_987;

/// Shortest stem or duplex reported.
pub const MIN_PAIRS: usize = 3;
/// Shortest hairpin loop.
pub const MIN_LOOP: usize = 3;
/// Longest primer analysed.
const MAX_LEN: usize = 200;

/// Warning thresholds (kcal/mol), as commonly recommended for PCR primers.
pub const DIMER_WARN: f64 = -9.0;
pub const DIMER_3PRIME_WARN: f64 = -5.0;
pub const HAIRPIN_WARN: f64 = -3.0;
pub const HAIRPIN_3PRIME_WARN: f64 = -2.0;

/// A duplex between two primers (or a primer and itself).
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct Dimer {
    /// ΔG°37, kcal/mol.
    pub dg: f64,
    /// Paired stretch on the first primer, `a_start..a_start + pairs` (5'->3' indices).
    pub a_start: usize,
    /// Paired stretch on the second primer, `b_start..b_start + pairs`.
    pub b_start: usize,
    pub pairs: usize,
    /// The first primer's 3' base is paired, so the dimer can be extended.
    pub a_three_prime: bool,
    /// The second primer's 3' base is paired.
    pub b_three_prime: bool,
}

/// A hairpin: `stem` bases at `start` pair with `stem` bases ending at `end` (exclusive),
/// around a loop.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct Hairpin {
    pub dg: f64,
    pub start: usize,
    pub end: usize,
    pub stem: usize,
    pub loop_len: usize,
    /// The primer's 3' end is in the stem (it can extend on itself).
    pub three_prime: bool,
}

fn norm(seq: &[u8]) -> Option<Vec<u8>> {
    let s: Vec<u8> = seq
        .iter()
        .map(|b| match b.to_ascii_uppercase() {
            b'U' => b'T',
            o => o,
        })
        .collect();
    (s.len() <= MAX_LEN && s.iter().all(|b| b"ACGT".contains(b))).then_some(s)
}

fn pairs_with(a: u8, b: u8) -> bool {
    matches!(
        (a, b),
        (b'A', b'T') | (b'T', b'A') | (b'C', b'G') | (b'G', b'C')
    )
}

/// ΔH (kcal/mol) and ΔS (cal/(K·mol), 1 M Na+) of the stacks along `s` (one strand of a
/// perfectly paired stretch), without initiation.
fn stacks(s: &[u8]) -> Option<(f64, f64)> {
    let (mut dh, mut ds) = (0.0, 0.0);
    for w in s.windows(2) {
        let (h, e) = crate::tm::nn_pair(*w.first()?, *w.get(1)?)?;
        dh += h;
        ds += e;
    }
    Some((dh, ds))
}

fn salt_ds(pairs: usize, c: &TmConditions) -> f64 {
    let na = crate::tm::sodium_equivalent(c);
    if na > 0.0 && na.is_finite() && pairs > 1 {
        0.368 * (pairs - 1) as f64 * na.ln()
    } else {
        0.0
    }
}

/// ΔG°37 of a perfectly paired duplex whose one strand is `s` (with initiation).
fn duplex_dg(s: &[u8], c: &TmConditions) -> Option<f64> {
    let (mut dh, mut ds) = stacks(s)?;
    for end in [s.first()?, s.last()?] {
        let (h, e) = crate::tm::nn_init(*end);
        dh += h;
        ds += e;
    }
    ds += salt_ds(s.len(), c);
    Some(dh - T37 * ds / 1000.0)
}

/// Hairpin loop ΔG°37 (SantaLucia & Hicks 2004, Table 4), interpolated between listed sizes and
/// extrapolated beyond 30 with 2.44·R·T·ln(n/30).
fn loop_dg(n: usize) -> f64 {
    const TABLE: [(usize, f64); 15] = [
        (3, 3.5),
        (4, 3.5),
        (5, 3.3),
        (6, 4.0),
        (7, 4.2),
        (8, 4.3),
        (9, 4.5),
        (10, 4.6),
        (12, 5.0),
        (14, 5.1),
        (16, 5.3),
        (18, 5.5),
        (20, 5.7),
        (25, 6.1),
        (30, 6.5),
    ];
    if n > 30 {
        return 6.5 + 2.44 * R_KCAL * T37 * (n as f64 / 30.0).ln();
    }
    let mut prev = (3, 3.5);
    for &(size, dg) in &TABLE {
        if n == size {
            return dg;
        }
        if n < size {
            let f = (n - prev.0) as f64 / (size - prev.0) as f64;
            return prev.1 + f * (dg - prev.1);
        }
        prev = (size, dg);
    }
    6.5
}

/// The most stable dimer between `a` and `b` (both 5'->3'), and the most stable one that pairs
/// a 3' end of either. `None` if no stretch of [`MIN_PAIRS`] pairs or the sequences aren't
/// plain DNA.
pub fn dimers(a: &[u8], b: &[u8], c: &TmConditions) -> Option<(Option<Dimer>, Option<Dimer>)> {
    let (a, b) = (norm(a)?, norm(b)?);
    let (na, nb) = (a.len(), b.len());
    if na == 0 || nb == 0 {
        return Some((None, None));
    }
    let mut best: Option<Dimer> = None;
    let mut best3: Option<Dimer> = None;
    // Antiparallel: a[i] pairs with b[j] where i + j = d (a runs 5'->3', b runs 3'->5').
    for d in 0..(na + nb).saturating_sub(1) {
        let mut i = d.saturating_sub(nb - 1);
        let last = d.min(na - 1);
        while i <= last {
            let j = d - i;
            if !a
                .get(i)
                .zip(b.get(j))
                .is_some_and(|(&x, &y)| pairs_with(x, y))
            {
                i += 1;
                continue;
            }
            // A maximal run starting at i (j decreasing as i increases).
            let mut k = 0;
            while i + k <= last
                && j >= k
                && a.get(i + k)
                    .zip(b.get(j - k))
                    .is_some_and(|(&x, &y)| pairs_with(x, y))
            {
                k += 1;
            }
            if k >= MIN_PAIRS
                && let Some(seg) = a.get(i..i + k)
                && let Some(dg) = duplex_dg(seg, c)
            {
                let dimer = Dimer {
                    dg,
                    a_start: i,
                    b_start: j + 1 - k,
                    pairs: k,
                    a_three_prime: i + k == na,
                    b_three_prime: j + 1 == nb,
                };
                let three = dimer.a_three_prime || dimer.b_three_prime;
                if best.as_ref().is_none_or(|x| dg < x.dg) {
                    best = Some(dimer.clone());
                }
                if three && best3.as_ref().is_none_or(|x| dg < x.dg) {
                    best3 = Some(dimer);
                }
            }
            i += k.max(1);
        }
    }
    Some((best, best3))
}

/// The most stable hairpin of `s` (5'->3'), if any has a stem of [`MIN_PAIRS`] or more.
pub fn hairpin(s: &[u8], c: &TmConditions) -> Option<Hairpin> {
    let s = norm(s)?;
    let n = s.len();
    let mut best: Option<Hairpin> = None;
    // Stem: s[i..i+k] pairs with s[j-k..j] reversed, loop = (j - k) - (i + k).
    for i in 0..n {
        for j in (i + 2 * MIN_PAIRS + MIN_LOOP)..=n {
            let mut k = 0;
            while i + k < j.saturating_sub(k + 1)
                && s.get(i + k)
                    .zip(s.get(j - 1 - k))
                    .is_some_and(|(&x, &y)| pairs_with(x, y))
            {
                k += 1;
            }
            // Keep a loop of at least MIN_LOOP.
            while k > 0 && (j - k).saturating_sub(i + k) < MIN_LOOP {
                k -= 1;
            }
            if k < MIN_PAIRS {
                continue;
            }
            let loop_len = (j - k) - (i + k);
            let Some(arm) = s.get(i..i + k) else { continue };
            let Some((dh, ds)) = stacks(arm) else {
                continue;
            };
            let ds = ds + salt_ds(k, c);
            let dg = dh - T37 * ds / 1000.0 + loop_dg(loop_len);
            if best.as_ref().is_none_or(|h| dg < h.dg) {
                best = Some(Hairpin {
                    dg,
                    start: i,
                    end: j,
                    stem: k,
                    loop_len,
                    three_prime: j == n,
                });
            }
        }
    }
    best
}

/// Plain-language problems with a primer's structure (empty if none).
pub fn warnings(seq: &[u8], c: &TmConditions) -> Vec<String> {
    let mut out = Vec::new();
    if let Some(h) = hairpin(seq, c) {
        let limit = if h.three_prime {
            HAIRPIN_3PRIME_WARN
        } else {
            HAIRPIN_WARN
        };
        if h.dg <= limit {
            let at = if h.three_prime { " at the 3' end" } else { "" };
            out.push(format!("hairpin{at}, \u{394}G {:.1} kcal/mol", h.dg));
        }
    }
    if let Some((any, three)) = dimers(seq, seq, c) {
        if let Some(d) = three.filter(|d| d.dg <= DIMER_3PRIME_WARN) {
            out.push(format!(
                "self-dimer at the 3' end, \u{394}G {:.1} kcal/mol",
                d.dg
            ));
        } else if let Some(d) = any.filter(|d| d.dg <= DIMER_WARN) {
            out.push(format!("self-dimer, \u{394}G {:.1} kcal/mol", d.dg));
        }
    }
    out
}

/// Cross-dimer problems between two primers (empty if none).
pub fn pair_warnings(a: &[u8], b: &[u8], c: &TmConditions) -> Vec<String> {
    let mut out = Vec::new();
    if let Some((any, three)) = dimers(a, b, c) {
        if let Some(d) = three.filter(|d| d.dg <= DIMER_3PRIME_WARN) {
            out.push(format!(
                "primers pair at a 3' end, \u{394}G {:.1} kcal/mol",
                d.dg
            ));
        } else if let Some(d) = any.filter(|d| d.dg <= DIMER_WARN) {
            out.push(format!(
                "primers pair with each other, \u{394}G {:.1} kcal/mol",
                d.dg
            ));
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn c() -> TmConditions {
        TmConditions::default()
    }

    #[test]
    fn hairpin_by_hand() {
        // GGGG/CCCC stem around a 4-nt loop: three GG/CC stacks (ΔH -8.0, ΔS -19.9 each), the
        // salt term 0.368·3·ln(0.1), and the 4-nt loop's +3.5.
        let h = hairpin(b"GGGGAAAACCCC", &c()).unwrap();
        assert_eq!(
            (h.start, h.end, h.stem, h.loop_len, h.three_prime),
            (0, 12, 4, 4, true)
        );
        let ds = 3.0 * -19.9 + 0.368 * 3.0 * (0.1f64).ln();
        let expected = 3.0 * -8.0 - T37 * ds / 1000.0 + 3.5;
        assert!((h.dg - expected).abs() < 1e-9, "{} vs {expected}", h.dg);
        assert!((h.dg - -1.20).abs() < 0.02);
        // No stem, no hairpin; too short a loop doesn't count.
        assert!(hairpin(b"AAAAAAAAAAAA", &c()).is_none());
        assert!(hairpin(b"GGGCCC", &c()).is_none());
    }

    #[test]
    fn dimers_by_hand() {
        // A palindromic 3' end pairs with itself: GAATTC/GAATTC.
        let (any, three) = dimers(b"TTTTTTGAATTC", b"TTTTTTGAATTC", &c()).unwrap();
        let three = three.unwrap();
        assert_eq!((three.a_start, three.b_start, three.pairs), (6, 6, 6));
        assert!(three.a_three_prime && three.b_three_prime);
        let mut ds = 0.0;
        let mut dh = 0.0;
        for w in b"GAATTC".windows(2) {
            let (h, s) = crate::tm::nn_pair(w[0], w[1]).unwrap();
            dh += h;
            ds += s;
        }
        // Terminal G·C at both ends.
        dh += 0.2;
        ds += -5.6 + 0.368 * 5.0 * (0.1f64).ln();
        assert!((three.dg - (dh - T37 * ds / 1000.0)).abs() < 1e-9);
        assert_eq!(any.unwrap().dg, three.dg);
        // Nothing pairs.
        assert_eq!(dimers(b"AAAAAAAA", b"CCCCCCCC", &c()), Some((None, None)));
    }

    #[test]
    fn warnings_flag_bad_primers() {
        // A strong 3' palindrome.
        let w = warnings(b"ACGTACGTGCGGCCGC", &c());
        assert!(w.iter().any(|x| x.contains("self-dimer")), "{w:?}");
        // A good, unstructured primer.
        assert!(
            warnings(b"GTAAAACGACGGCCAGTG", &c()).is_empty(),
            "{:?}",
            warnings(b"GTAAAACGACGGCCAGTG", &c())
        );
        let w = pair_warnings(b"TTTTTTTTGCGGCCGC", b"AAAAAAAAGCGGCCGC", &c());
        assert!(!w.is_empty());
    }

    #[test]
    fn loops() {
        assert_eq!(loop_dg(4), 3.5);
        assert!((loop_dg(11) - 4.8).abs() < 1e-9);
        assert!(loop_dg(60) > 6.5);
    }

    #[test]
    fn hostile_inputs() {
        for s in [
            &b""[..],
            b"A",
            b"NNNN",
            b"ACG-T",
            &[0xff; 5],
            &vec![b'G'; 500],
        ] {
            let _ = hairpin(s, &c());
            let _ = dimers(s, b"ACGT", &c());
            let _ = warnings(s, &c());
        }
        let long: Vec<u8> = (0..200).map(|i| b"ACGT"[i % 4]).collect();
        assert!(hairpin(&long, &c()).is_some());
    }
}
