//! Where primers bind a template.
//!
//! A primer is extended by polymerase from its 3' end, so what decides priming is how well its
//! 3'-terminal bases pair; a 5' tail (restriction sites, overlaps for assembly) may hang free.
//! Mismatches at the very 3' end block extension, while a single mismatch a few bases in is
//! tolerated (Kwok et al. 1990, Nucleic Acids Res. 18:999-1005). The rule used here
//! ([`BindingRules`]): the last `window` bases of the primer pair with the template with at most
//! `mismatches` mismatches, none of them in the last `clamp` bases; the annealed stretch then
//! extends 5'-ward for as long as bases keep pairing. Matching uses IUPAC codes (a primer `N`
//! pairs with anything; an ambiguous template base only counts if every base it stands for
//! would pair). Circular templates are searched across the origin.
//!
//! With `isolated` set, a mismatch is only tolerated when the bases on both sides of it pair
//! (a single mismatched base inside paired DNA, rather than a run of them), and the annealed
//! stretch extends 5'-ward across any number of such mismatches. `mismatches` is the budget
//! for the 3' window, which decides priming; behind it a repair or mutagenesis oligo carries
//! as many designed changes as its design needs, and each one still leaves the primer annealed
//! either side of it. [`BindingPrefs`] are the user's choice of rules (Preferences): how many 3'
//! bases must pair for short and for long primers, isolated mismatches, and a minimum Tm.

use serde::Serialize;
use snappiesnap_doc::{Molecule, Range, Topology};

/// Fewest bases that can count as a binding site.
pub const MIN_ANNEAL: usize = 10;

/// How strictly a primer must pair to count as binding.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, serde::Deserialize)]
#[serde(default)]
pub struct BindingRules {
    /// 3'-terminal bases that must pair (at least [`MIN_ANNEAL`]; the whole primer if shorter).
    pub window: usize,
    /// Mismatches allowed in that window.
    pub mismatches: usize,
    /// 3'-terminal bases that must pair exactly.
    pub clamp: usize,
    /// Only mismatches whose neighbours both pair are tolerated; 5' of the window, any number
    /// of them may be crossed.
    pub isolated: bool,
}

impl Default for BindingRules {
    /// A 10-base 3' window with one mismatch allowed outside the last 3 bases.
    fn default() -> Self {
        BindingRules {
            window: 10,
            mismatches: 1,
            clamp: 3,
            isolated: false,
        }
    }
}

impl BindingRules {
    /// Exact pairing of the last `window` bases.
    pub fn exact(window: usize) -> Self {
        BindingRules {
            window,
            mismatches: 0,
            clamp: window,
            isolated: false,
        }
    }
}

/// The rules for one class of templates (short or long sequences), as the user sets them in
/// Preferences: the primer's 3'-terminal `matched` bases must pair, with one isolated
/// mismatch tolerated when `isolated_mismatches` is set (never in the last 3 bases: see
/// [`BindingRules`]). With `min_tm_on`, sites whose paired stretch melts below `min_tm` °C
/// ([`BindingSite::tm`]) don't count (the value is kept while the minimum is off).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, serde::Deserialize)]
#[serde(default)]
pub struct ClassRules {
    pub matched: usize,
    pub isolated_mismatches: bool,
    pub min_tm_on: bool,
    pub min_tm: f64,
}

impl Default for ClassRules {
    /// The short-sequence defaults: 10 bases, an isolated mismatch allowed, at least 40 °C.
    fn default() -> Self {
        ClassRules {
            matched: 10,
            isolated_mismatches: true,
            min_tm_on: true,
            min_tm: 40.0,
        }
    }
}

impl ClassRules {
    /// Matched bases allowed.
    pub const MATCH: std::ops::RangeInclusive<usize> = MIN_ANNEAL..=60;
    /// Minimum Tm values allowed (°C).
    pub const MIN_TM: std::ops::RangeInclusive<f64> = 0.0..=90.0;

    /// Keep every value in range (a Tm that isn't a number becomes 40 °C).
    pub fn clamped(self) -> Self {
        ClassRules {
            matched: self.matched.clamp(*Self::MATCH.start(), *Self::MATCH.end()),
            isolated_mismatches: self.isolated_mismatches,
            min_tm_on: self.min_tm_on,
            min_tm: if self.min_tm.is_finite() {
                self.min_tm
                    .clamp(*Self::MIN_TM.start(), *Self::MIN_TM.end())
            } else {
                40.0
            },
        }
    }

    /// The search rules.
    pub fn rules(&self) -> BindingRules {
        if self.isolated_mismatches {
            BindingRules {
                window: self.matched,
                mismatches: 1,
                clamp: 3,
                isolated: true,
            }
        } else {
            BindingRules::exact(self.matched)
        }
    }
}

/// Which binding sites count, as the user sets it in Preferences: on a template of
/// `long_from_bp` bases or more (a chromosome, say) primers follow the `long` rules, on a
/// shorter one the `short` rules. Long templates have more chance matches, so they usually
/// get stricter rules. Primary sites ([`BindingSite::is_primary`]), where a primer is meant
/// to bind, always count.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, serde::Deserialize)]
#[serde(default)]
pub struct BindingPrefs {
    pub long_from_bp: usize,
    pub short: ClassRules,
    pub long: ClassRules,
    /// Also count sites where the primer's 5'-terminal `five_prime_matched` bases pair but its
    /// 3' end doesn't: it anneals there but can't be extended (for both classes).
    pub five_prime_on: bool,
    pub five_prime_matched: usize,
}

impl Default for BindingPrefs {
    /// Long from 1,000,000 bp. Short sequences: 10 matched bases, an isolated mismatch
    /// allowed. Long sequences: 15 matched bases, no mismatch. At least 40 °C for both.
    fn default() -> Self {
        BindingPrefs {
            long_from_bp: 1_000_000,
            short: ClassRules::default(),
            long: ClassRules {
                matched: 15,
                isolated_mismatches: false,
                ..ClassRules::default()
            },
            five_prime_on: true,
            five_prime_matched: 15,
        }
    }
}

impl BindingPrefs {
    /// Template lengths (bp) from which a sequence may count as long (fits 32-bit `usize`).
    pub const LONG_FROM: std::ops::RangeInclusive<usize> = 1_000..=1_000_000_000;

    /// Keep every value in range.
    pub fn clamped(self) -> Self {
        BindingPrefs {
            long_from_bp: self
                .long_from_bp
                .clamp(*Self::LONG_FROM.start(), *Self::LONG_FROM.end()),
            short: self.short.clamped(),
            long: self.long.clamped(),
            five_prime_on: self.five_prime_on,
            five_prime_matched: self
                .five_prime_matched
                .clamp(*ClassRules::MATCH.start(), *ClassRules::MATCH.end()),
        }
    }

    /// The rules for a template of `len` bases.
    pub fn class_for(&self, len: usize) -> &ClassRules {
        if len >= self.long_from_bp {
            &self.long
        } else {
            &self.short
        }
    }

    /// The search rules on a template of `len` bases.
    pub fn rules_for(&self, len: usize) -> BindingRules {
        self.class_for(len).rules()
    }
}

/// The 3' end of a site on the top strand: the position the primer is extended from.
fn three_prime(s: &BindingSite) -> (PrimerStrand, usize) {
    match s.strand {
        PrimerStrand::Forward => (s.strand, s.range.end),
        PrimerStrand::Reverse => (s.strand, s.range.start),
    }
}

/// Binding sites of `primers` (5'->3') on `m` under the user's [`BindingPrefs`]: each primer's
/// sites under the rules for `m`'s length ([`BindingPrefs::rules_for`]), plus its primary sites (which replace a
/// site found with the same 3' end), plus with `five_prime_on` its 5'-end-only sites
/// ([`five_prime_sites`]) where no site places the primer already, less those under the
/// minimum Tm in `conditions`. Ordered by primer, strand, position.
pub fn sites_with_prefs(
    m: &Molecule,
    primers: &[&[u8]],
    prefs: &BindingPrefs,
    conditions: &crate::tm::TmConditions,
) -> Vec<BindingSite> {
    let prefs = prefs.clamped();
    let class = prefs.class_for(m.len());
    let mut out = Vec::new();
    for (index, primer) in primers.iter().enumerate() {
        let n = primer.len();
        let mut mine = binding_sites(m, &[primer], class.rules());
        // A primer is meant to bind where it was designed to, so its primary site counts even
        // under rules that would miss it. Where both sets of rules find the same 3' end, keep
        // whichever anneals further: the stricter rules can only stop sooner, never show more.
        for site in binding_sites(m, &[primer], BindingRules::default()) {
            if !site.is_primary(n) {
                continue;
            }
            let end = three_prime(&site);
            if mine
                .iter()
                .any(|s| three_prime(s) == end && s.annealed >= site.annealed)
            {
                continue;
            }
            mine.retain(|s| three_prime(s) != end);
            mine.push(site);
        }
        if prefs.five_prime_on {
            for (site, end) in five_prime_sites(m, primer, prefs.five_prime_matched) {
                if !mine.iter().any(|s| three_prime(s) == end) {
                    mine.push(site);
                }
            }
        }
        mine.retain(|s| {
            s.is_primary(n)
                || !class.min_tm_on
                || s.tm(primer, conditions).is_some_and(|t| t >= class.min_tm)
        });
        mine.sort_by_key(|s| {
            (
                s.strand == PrimerStrand::Reverse,
                s.range.start,
                s.range.end,
            )
        });
        for mut s in mine {
            s.primer = index;
            out.push(s);
        }
    }
    out
}

/// Which way the primer extends: `Forward` primers have the top strand's sequence and extend
/// rightwards; `Reverse` primers pair with the top strand and extend leftwards.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[serde(rename_all = "lowercase")]
pub enum PrimerStrand {
    Forward,
    Reverse,
}

/// A binding site of a primer.
#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
pub struct BindingSite {
    /// Index of the primer in the list searched.
    pub primer: usize,
    pub strand: PrimerStrand,
    /// The annealed stretch of the template (wraps on circular molecules).
    pub range: Range,
    /// Annealed bases (the 3' part of the primer); the rest is a 5' tail.
    pub annealed: usize,
    /// Mismatches within the annealed stretch.
    pub mismatches: usize,
    /// Where the mismatches are, as distances from the primer's 3' end (0 = the 3' base).
    pub mismatch_at: Vec<usize>,
    /// Only the primer's 5' end pairs here (`annealed` 5'-terminal bases), not its 3' end: it
    /// anneals but can't be extended.
    #[serde(skip_serializing_if = "std::ops::Not::not")]
    pub five_prime_only: bool,
}

impl BindingSite {
    /// Tm (°C) of the longest stretch of the annealed part that pairs without mismatches (all
    /// of it if there are none); `primer` is the whole primer, 5'->3'. Mismatch thermodynamics
    /// aren't modelled, so this is the conservative estimate.
    pub fn tm(&self, primer: &[u8], conditions: &crate::tm::TmConditions) -> Option<f64> {
        if self.five_prime_only {
            return crate::tm::tm(primer.get(..self.annealed)?, conditions);
        }
        let n = primer.len();
        let mut cuts: Vec<usize> = self
            .mismatch_at
            .iter()
            .copied()
            .filter(|&d| d < self.annealed)
            .collect();
        cuts.sort_unstable();
        // Runs between mismatches, as distances from the 3' end: [from, to).
        let mut best = (0, 0);
        let mut from = 0;
        for &m in cuts.iter().chain(std::iter::once(&self.annealed)) {
            if m.saturating_sub(from) > best.1 - best.0 {
                best = (from, m);
            }
            from = m + 1;
        }
        // Distances from the 3' end [a, b) are primer indices n-b..n-a.
        let run = primer.get(n.checked_sub(best.1)?..n.checked_sub(best.0)?)?;
        crate::tm::tm(run, conditions)
    }

    /// A site where the primer pairs exactly over at least 15 bases (or its whole length, if
    /// shorter): where it is meant to bind, as opposed to a weaker secondary site. Deliberately
    /// strict about mismatches — [`crate::pcr`] amplifies from primary sites, so loosening it
    /// would invent products. A repair or mutagenesis oligo, which carries changes on purpose,
    /// is recognised by [`BindingSite::anneals_whole_primer`] instead.
    pub fn is_primary(&self, primer_len: usize) -> bool {
        self.mismatches == 0 && self.annealed >= primer_len.min(15)
    }

    /// Whether the primer pairs along its whole length here, mismatches and all: a repair or
    /// mutagenesis oligo at the site it was designed for. A chance site almost never anneals
    /// over a primer's entire length, so this says "drawn as a real binding site" without
    /// saying "polymerase extends from it".
    pub fn anneals_whole_primer(&self, primer_len: usize) -> bool {
        primer_len > 0 && self.annealed == primer_len
    }
}

/// Bit mask of the bases an IUPAC code stands for (A=1, C=2, G=4, T=8).
fn mask(b: u8) -> u8 {
    match b.to_ascii_uppercase() {
        b'A' => 1,
        b'C' => 2,
        b'G' => 4,
        b'T' | b'U' => 8,
        b'R' => 5,
        b'Y' => 10,
        b'S' => 6,
        b'W' => 9,
        b'K' => 12,
        b'M' => 3,
        b'B' => 14,
        b'D' => 13,
        b'H' => 11,
        b'V' => 7,
        b'N' => 15,
        _ => 0,
    }
}

fn complement_mask(m: u8) -> u8 {
    // Swap A<->T (bits 0 and 3) and C<->G (bits 1 and 2).
    ((m & 1) << 3) | ((m & 8) >> 3) | ((m & 2) << 1) | ((m & 4) >> 1)
}

/// Whether a template base can pair as primer base `p` (both given as top-strand masks).
fn pairs(p: u8, template: u8) -> bool {
    template != 0 && template & p == template
}

/// Longest primer searched.
const MAX_PRIMER: usize = 10_000;

/// Binding sites of `primers` (sequences 5'->3') on `m`, ordered by primer, strand, position.
pub fn binding_sites(m: &Molecule, primers: &[&[u8]], rules: BindingRules) -> Vec<BindingSite> {
    let seq: Vec<u8> = m.seq().as_bytes().iter().map(|&b| mask(b)).collect();
    let len = seq.len();
    let circular = m.topology() == Topology::Circular;
    let mut out = Vec::new();
    if len == 0 || !m.seq().alphabet().is_nucleotide() {
        return out;
    }
    let l = len as i64;
    let at = |i: i64| -> Option<u8> {
        let i = if circular { i.rem_euclid(l) } else { i };
        usize::try_from(i).ok().and_then(|i| seq.get(i).copied())
    };
    for (index, primer) in primers.iter().enumerate() {
        let n = primer.len();
        let k = rules.window.max(MIN_ANNEAL).min(n);
        if !(MIN_ANNEAL..=MAX_PRIMER).contains(&n) || k > len {
            continue;
        }
        let fwd: Vec<u8> = primer.iter().map(|&b| mask(b)).collect();
        if fwd.contains(&0) {
            continue;
        }
        // The primer as it lies on the top strand when it binds as a reverse primer.
        let rev: Vec<u8> = fwd.iter().rev().map(|&b| complement_mask(b)).collect();
        let anchors = if circular { len } else { len + 1 - k };
        for strand in [PrimerStrand::Forward, PrimerStrand::Reverse] {
            for a in 0..anchors {
                let a = a as i64;
                // Primer base and template base `j` positions from the primer's 3' end.
                let pair_at = |j: usize| -> bool {
                    let (p, t) = match strand {
                        // 3' end pairs with top position a + k - 1.
                        PrimerStrand::Forward => {
                            (fwd.get(n - 1 - j).copied(), at(a + k as i64 - 1 - j as i64))
                        }
                        // 3' end pairs with top position a.
                        PrimerStrand::Reverse => (rev.get(j).copied(), at(a + j as i64)),
                    };
                    matches!((p, t), (Some(p), Some(t)) if pairs(p, t))
                };
                // An isolated mismatch has a pairing base on both sides, within the primer.
                let isolated = |j: usize| j > 0 && j + 1 < n && pair_at(j - 1) && pair_at(j + 1);
                let mut mismatches = 0;
                let mut mismatch_at = Vec::new();
                let mut ok = true;
                for j in 0..k {
                    if !pair_at(j) {
                        mismatch_at.push(j);
                        mismatches += 1;
                        if j < rules.clamp
                            || mismatches > rules.mismatches
                            || (rules.isolated && !isolated(j))
                        {
                            ok = false;
                            break;
                        }
                    }
                }
                if !ok {
                    continue;
                }
                let mut annealed = k;
                while annealed < n && annealed < len {
                    if pair_at(annealed) {
                        annealed += 1;
                    } else if rules.isolated && annealed + 1 < len && isolated(annealed) {
                        mismatch_at.push(annealed);
                        mismatches += 1;
                        annealed += 1;
                    } else {
                        break;
                    }
                }
                let range = match strand {
                    PrimerStrand::Forward => {
                        let end = a + k as i64;
                        wrap(end - annealed as i64, end, l)
                    }
                    PrimerStrand::Reverse => wrap(a, a + annealed as i64, l),
                };
                out.push(BindingSite {
                    primer: index,
                    strand,
                    range,
                    annealed,
                    mismatches,
                    mismatch_at,
                    five_prime_only: false,
                });
            }
        }
    }
    out
}

/// Sites where `primer`'s 5'-terminal `min` bases or more pair exactly with `m` but the rest of
/// the primer, up to its 3' end, doesn't: the primer anneals there but polymerase can't extend
/// it. Each with the top-strand position its 3' end would have (as [`three_prime`] gives it
/// for a full site), so callers can drop places a full site already covers.
pub fn five_prime_sites(
    m: &Molecule,
    primer: &[u8],
    min: usize,
) -> Vec<(BindingSite, (PrimerStrand, usize))> {
    let seq: Vec<u8> = m.seq().as_bytes().iter().map(|&b| mask(b)).collect();
    let len = seq.len();
    let n = primer.len();
    let mut out = Vec::new();
    let min = min.max(MIN_ANNEAL);
    if len == 0 || n <= min || n > MAX_PRIMER || min > len || !m.seq().alphabet().is_nucleotide() {
        return out;
    }
    let fwd: Vec<u8> = primer.iter().map(|&b| mask(b)).collect();
    if fwd.contains(&0) {
        return out;
    }
    let circular = m.topology() == Topology::Circular;
    let l = len as i64;
    let at = |i: i64| -> Option<u8> {
        let i = if circular { i.rem_euclid(l) } else { i };
        usize::try_from(i).ok().and_then(|i| seq.get(i).copied())
    };
    // Primer base `i` (from the 5' end) pairs with top position `pos(i)`; forward primers lie
    // on the top strand 5'->3' rightwards, reverse primers pair leftwards from their 5' end.
    for strand in [PrimerStrand::Forward, PrimerStrand::Reverse] {
        let anchors: Vec<i64> = if circular {
            (0..l).collect()
        } else {
            match strand {
                PrimerStrand::Forward => (0..=l - min as i64).collect(),
                PrimerStrand::Reverse => (min as i64 - 1..l).collect(),
            }
        };
        for a in anchors {
            let pairs_at = |i: usize| -> bool {
                let p = fwd.get(i).copied();
                let t = match strand {
                    PrimerStrand::Forward => at(a + i as i64),
                    PrimerStrand::Reverse => at(a - i as i64).map(complement_mask),
                };
                matches!((p, t), (Some(p), Some(t)) if pairs(p, t))
            };
            let mut annealed = 0;
            while annealed < n && annealed < len && pairs_at(annealed) {
                annealed += 1;
            }
            // The whole primer pairs: a full site, found by the 3' search.
            if annealed < min || annealed >= n {
                continue;
            }
            let (range, end) = match strand {
                PrimerStrand::Forward => (
                    wrap(a, a + annealed as i64, l),
                    (a + n as i64).rem_euclid(l.max(1)),
                ),
                PrimerStrand::Reverse => (
                    wrap(a + 1 - annealed as i64, a + 1, l),
                    (a + 1 - n as i64).rem_euclid(l.max(1)),
                ),
            };
            // On a linear template the 3' end may lie past either end.
            let end = if circular {
                usize::try_from(end).unwrap_or(0)
            } else {
                match strand {
                    PrimerStrand::Forward => usize::try_from(a + n as i64).unwrap_or(usize::MAX),
                    PrimerStrand::Reverse => {
                        usize::try_from(a + 1 - n as i64).unwrap_or(usize::MAX)
                    }
                }
            };
            out.push((
                BindingSite {
                    primer: 0,
                    strand,
                    range,
                    annealed,
                    mismatches: 0,
                    mismatch_at: Vec::new(),
                    five_prime_only: true,
                },
                (strand, end),
            ));
        }
    }
    out
}

/// `start..end` (either may run past the origin) as a range on a molecule of length `l`.
fn wrap(start: i64, end: i64, l: i64) -> Range {
    let s = usize::try_from(start.rem_euclid(l.max(1))).unwrap_or(0);
    let e = if end > l {
        usize::try_from(end - l).unwrap_or(0)
    } else {
        usize::try_from(end).unwrap_or(0)
    };
    Range::new(s, e)
}

#[cfg(test)]
mod tests {
    use snappiesnap_doc::{Alphabet, Seq};

    use super::*;

    fn mol(s: &str, t: Topology) -> Molecule {
        Molecule::new("m", Seq::new(Alphabet::Dna, s).unwrap(), t).unwrap()
    }

    const TEMPLATE: &str = "AAAAAAAAAAACGTTGCATGCAAGCTTGGCACTGGCCGTCGTTTTACAAAAAAAAAA";

    #[test]
    fn forward_and_reverse_with_tails() {
        let m = mol(TEMPLATE, Topology::Linear);
        // Forward: a 5' tail (gaattc) on 16 matching bases.
        let fwd = b"gaattcCGTTGCATGCAAGCTT";
        // Reverse: reverse complement of TGGCACTGGCCGTCGTTTTAC with a tail.
        let rev = b"ggatccGTAAAACGACGGCCAGTGCCA";
        let sites = binding_sites(&m, &[fwd, rev], BindingRules::exact(15));
        assert_eq!(sites.len(), 2, "{sites:?}");
        assert_eq!(sites[0].strand, PrimerStrand::Forward);
        assert_eq!(sites[0].range, Range::new(11, 27));
        assert_eq!(sites[0].annealed, 16);
        assert_eq!(sites[1].strand, PrimerStrand::Reverse);
        assert_eq!(sites[1].range, Range::new(26, 47));
        assert_eq!(sites[1].annealed, 21);
    }

    #[test]
    fn a_3prime_mismatch_prevents_binding() {
        let m = mol(TEMPLATE, Topology::Linear);
        assert!(binding_sites(&m, &[b"CGTTGCATGCAAGCTA"], BindingRules::default()).is_empty());
        // A mismatch outside the 3' core only shortens the annealed part.
        let s = binding_sites(&m, &[b"CGATGCATGCAAGCTTGG"], BindingRules::exact(15));
        assert_eq!(s.len(), 1);
        assert_eq!(s[0].annealed, 15);
    }

    #[test]
    fn one_mismatch_away_from_the_3prime_end() {
        let m = mol(TEMPLATE, Topology::Linear);
        // CGTTGCATGCAAGCTT with a mismatch 3 bases in from the 3' end (A for G).
        let s = binding_sites(&m, &[b"CGTTGCATGCAAACTT"], BindingRules::default());
        assert_eq!(s.len(), 1);
        assert_eq!(
            (s[0].mismatches, s[0].annealed, &s[0].mismatch_at),
            (1, 16, &vec![3])
        );
        let primer = b"CGTTGCATGCAAACTT";
        // The longest perfectly paired run is the 12 bases 5' of the mismatch.
        let c = crate::tm::TmConditions::default();
        assert_eq!(s[0].tm(primer, &c), crate::tm::tm(b"CGTTGCATGCAA", &c));
        assert!(binding_sites(&m, &[b"CGTTGCATGCAAACTT"], BindingRules::exact(10)).is_empty());
        // In the 3' clamp it doesn't bind.
        assert!(binding_sites(&m, &[b"CGTTGCATGCAAGCAT"], BindingRules::default()).is_empty());
    }

    #[test]
    fn isolated_mismatches() {
        let m = mol(TEMPLATE, Topology::Linear);
        let isolated = BindingPrefs::default().rules_for(m.len());
        assert_eq!(
            isolated,
            BindingRules {
                window: 10,
                mismatches: 1,
                clamp: 3,
                isolated: true
            }
        );
        // Mismatches 9 and 10 bases from the 3' end: the one in the window has an unpaired
        // neighbour, so only the plain rule lets it pass.
        let p: &[u8] = b"CGTTGGTTGCAAGCTT";
        assert_eq!(binding_sites(&m, &[p], BindingRules::default()).len(), 1);
        assert!(binding_sites(&m, &[p], isolated).is_empty());
        // A single mismatch 12 bases in, outside the window: the plain rule stops the
        // annealed stretch there; isolated mismatches carry it on to the 5' end.
        let p: &[u8] = b"CGTAGCATGCAAGCTT";
        let plain = binding_sites(&m, &[p], BindingRules::default());
        assert_eq!((plain[0].annealed, plain[0].mismatches), (12, 0));
        let s = binding_sites(&m, &[p], isolated);
        assert_eq!(
            (s[0].annealed, s[0].mismatches, &s[0].mismatch_at),
            (16, 1, &vec![12])
        );
        assert_eq!(
            s[0].range,
            Range::new(plain[0].range.end - 16, plain[0].range.end)
        );
    }

    /// A repair oligo carries the changes it is meant to make, so it pairs along its whole
    /// length with several isolated mismatches in it. A synthetic pair in the shape of a real
    /// repair design: 100 bases, with mismatches 49, 74 and 79 from the 3' end.
    #[test]
    fn a_repair_oligo_anneals_across_the_changes_it_carries() {
        let template = "\
AGTCATCTAGTAGAAGCACGTGCGGAGTTATCCTTCTACGTAATCGTCTGAACAGAGGTTACCGTTAACCTGGAATCCAGAATTAAGCGGAACTTCCGTTGAAGAAGATTCTCT";
        let m = Molecule::new(
            "t",
            Seq::new(Alphabet::Dna, template).expect("a template"),
            Topology::Linear,
        )
        .expect("a molecule");
        let repair = b"GTCATCTAGTAGAAGCACGTTCGGACTTATCCTTCTACGTAATCGTCTGAGCAGAGGTTACCGTTAACCTGGAATCCAGAATTAAGCGGAACTTCCGTTG";
        assert_eq!(repair.len(), 100);
        let c = crate::tm::TmConditions::default();
        let sites = sites_with_prefs(&m, &[repair], &BindingPrefs::default(), &c);
        let site = sites.first().expect("the oligo binds its template");
        assert_eq!(
            (site.annealed, site.mismatches),
            (100, 3),
            "the whole oligo pairs, with the three changes it carries"
        );
        assert_eq!(
            site.mismatch_at,
            vec![49, 74, 79],
            "offsets from the 3' end"
        );
        assert_eq!(site.range, Range::new(1, 101));
        // Polymerase extends from it, so it is not "primary" (that would invent PCR
        // products), but it is the site the oligo was designed for.
        assert!(!site.is_primary(100));
        assert!(site.anneals_whole_primer(100));
        // Without isolated mismatches the stretch still stops at the first change.
        let strict = BindingPrefs {
            short: ClassRules {
                isolated_mismatches: false,
                ..ClassRules::default()
            },
            ..BindingPrefs::default()
        };
        let strict = sites_with_prefs(&m, &[repair], &strict, &c);
        assert_eq!(
            strict.first().map(|s| s.annealed),
            Some(49),
            "the strict rules stop at the first change"
        );
    }

    #[test]
    fn preferences() {
        let prefs = BindingPrefs::default();
        // The class follows the template's length, not the primer's.
        assert_eq!(prefs.rules_for(999_999), ClassRules::default().rules());
        assert_eq!(prefs.rules_for(1_000_000), BindingRules::exact(15));
        assert_eq!(
            ClassRules::default().rules(),
            BindingRules {
                window: 10,
                mismatches: 1,
                clamp: 3,
                isolated: true
            }
        );
        assert_eq!(
            BindingPrefs {
                long_from_bp: 0,
                short: ClassRules {
                    matched: 0,
                    isolated_mismatches: false,
                    min_tm_on: true,
                    min_tm: f64::NAN,
                },
                long: ClassRules {
                    matched: 999,
                    min_tm: 500.0,
                    ..ClassRules::default()
                },
                five_prime_on: true,
                five_prime_matched: 0,
            }
            .clamped(),
            BindingPrefs {
                long_from_bp: 1_000,
                short: ClassRules {
                    matched: 10,
                    isolated_mismatches: false,
                    min_tm_on: true,
                    min_tm: 40.0,
                },
                long: ClassRules {
                    matched: 60,
                    min_tm: 90.0,
                    ..ClassRules::default()
                },
                five_prime_on: true,
                five_prime_matched: 10,
            }
        );
        let m = mol(TEMPLATE, Topology::Linear);
        let c = crate::tm::TmConditions::default();
        // A weak site (a mismatch 3 bases in) and a primary one, on a short template.
        let weak: &[u8] = b"CGTTGCATGCAAACTT";
        let primary: &[u8] = b"CGTTGCATGCAAGCTT";
        let found = |p: &BindingPrefs| sites_with_prefs(&m, &[weak, primary], p, &c);
        let with_short = |short: ClassRules| BindingPrefs { short, ..prefs };
        let none = with_short(ClassRules {
            min_tm_on: false,
            min_tm: 80.0,
            ..prefs.short
        });
        assert_eq!(
            found(&none).iter().map(|s| s.primer).collect::<Vec<_>>(),
            vec![0, 1]
        );
        // A high minimum Tm drops the weak site but never the primary one.
        let hot = with_short(ClassRules {
            min_tm: 80.0,
            ..prefs.short
        });
        let kept = found(&hot);
        assert_eq!(kept.len(), 1);
        assert_eq!((kept[0].primer, kept[0].is_primary(16)), (1, true));
        // The long sequences' minimum doesn't apply to a short template...
        let long_hot = BindingPrefs {
            long: ClassRules {
                min_tm: 80.0,
                ..prefs.long
            },
            ..none
        };
        assert_eq!(found(&long_hot).len(), 2);
        // ...until the template counts as long.
        let as_long = BindingPrefs {
            long_from_bp: 1_000,
            ..long_hot
        };
        let big = mol(
            &format!("{TEMPLATE}{}", "C".repeat(1_000)),
            Topology::Linear,
        );
        let kept = sites_with_prefs(&big, &[weak, primary], &as_long, &c);
        assert_eq!(kept.iter().map(|s| s.primer).collect::<Vec<_>>(), vec![1]);
        // Asking for more matched bases than a primary site has still keeps it.
        let strict = with_short(ClassRules {
            matched: 60,
            min_tm_on: false,
            ..prefs.short
        });
        assert!(
            found(&strict)
                .iter()
                .any(|s| s.primer == 1 && s.is_primary(16))
        );
    }

    #[test]
    fn across_the_origin() {
        let rotated = format!("{}{}", &TEMPLATE[20..], &TEMPLATE[..20]);
        let m = mol(&rotated, Topology::Circular);
        let s = binding_sites(&m, &[b"CGTTGCATGCAAGCTT"], BindingRules::default());
        assert_eq!(s.len(), 1);
        let len = rotated.len();
        assert_eq!(s[0].range, Range::new(len - 9, 7));
        assert!(
            binding_sites(
                &mol(&rotated, Topology::Linear),
                &[b"CGTTGCATGCAAGCTT"],
                BindingRules::default()
            )
            .is_empty()
        );
    }

    #[test]
    fn hostile_inputs() {
        for s in ["", "A", "NNNNNNNNNNNNNNNNNNNN", TEMPLATE] {
            for t in [Topology::Linear, Topology::Circular] {
                let m = mol(s, t);
                let primers: [&[u8]; 5] = [
                    b"",
                    b"ACGT",
                    b"NNNNNNNNNNNNNNNNNNNNNNNNN",
                    b"ACGT-ACGTACGTACGT",
                    b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
                ];
                for k in [0, 1, 15, 1000] {
                    for isolated in [false, true] {
                        let rules = BindingRules {
                            window: k,
                            mismatches: k % 3,
                            clamp: k / 2,
                            isolated,
                        };
                        for site in binding_sites(&m, &primers, rules) {
                            assert!(site.annealed <= m.len());
                        }
                        let class = ClassRules {
                            matched: k,
                            isolated_mismatches: isolated,
                            min_tm_on: isolated,
                            min_tm: k as f64,
                        };
                        let prefs = BindingPrefs {
                            long_from_bp: k,
                            short: class,
                            long: class,
                            ..BindingPrefs::default()
                        };
                        let c = crate::tm::TmConditions::default();
                        for site in sites_with_prefs(&m, &primers, &prefs, &c) {
                            assert!(site.annealed <= m.len());
                        }
                    }
                }
            }
        }
    }

    /// Reverse complement, for building templates.
    fn rc(s: &str) -> String {
        s.chars()
            .rev()
            .map(|c| match c {
                'A' => 'T',
                'T' => 'A',
                'C' => 'G',
                _ => 'C',
            })
            .collect()
    }

    #[test]
    fn five_prime_only_sites() {
        // The primer's first 18 bases are in the template, followed by bases that don't pair
        // with its last 6.
        let primer = "ACGTTGCATGCAAGCTTGGCACTG";
        let five = &primer[..18];
        let template = format!("TTTTTTTTTT{five}AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
        let m = mol(&template, Topology::Linear);
        let found = five_prime_sites(&m, primer.as_bytes(), 15);
        assert_eq!(found.len(), 1, "{found:?}");
        let (site, end) = &found[0];
        assert!(site.five_prime_only);
        assert_eq!(site.strand, PrimerStrand::Forward);
        assert_eq!(site.range, Range::new(10, 28));
        assert_eq!(site.annealed, 18);
        // Its 3' end would be at 10 + 24.
        assert_eq!(*end, (PrimerStrand::Forward, 34));
        // Fewer matched 5' bases than asked for: nothing.
        assert!(five_prime_sites(&m, primer.as_bytes(), 19).is_empty());
        // On the other strand the same stretch is a reverse 5' site.
        let m = mol(&rc(&template), Topology::Linear);
        let found = five_prime_sites(&m, primer.as_bytes(), 15);
        assert_eq!(found.len(), 1);
        let len = template.len();
        assert_eq!(found[0].0.strand, PrimerStrand::Reverse);
        assert_eq!(found[0].0.range, Range::new(len - 28, len - 10));
        // A primer that pairs all the way is a full site, not a 5' one.
        let m = mol(&format!("TTTTTTTTTT{primer}TTTTTTTTTT"), Topology::Linear);
        assert!(five_prime_sites(&m, primer.as_bytes(), 15).is_empty());
    }

    #[test]
    fn five_prime_sites_across_the_origin() {
        let primer = "ACGTTGCATGCAAGCTTGGCACTG";
        let five = &primer[..18];
        // The 5' stretch runs across the origin: last 8 bases at the end, 10 at the start.
        let template = format!("{}AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA{}", &five[8..], &five[..8]);
        let m = mol(&template, Topology::Circular);
        let found = five_prime_sites(&m, primer.as_bytes(), 15);
        assert_eq!(found.len(), 1, "{found:?}");
        let len = template.len();
        assert_eq!(found[0].0.range, Range::new(len - 8, 10));
        assert!(five_prime_sites(&mol("", Topology::Circular), primer.as_bytes(), 15).is_empty());
        assert!(five_prime_sites(&mol("A", Topology::Circular), primer.as_bytes(), 15).is_empty());
    }

    #[test]
    fn five_prime_sites_follow_the_preferences() {
        let primer = "ACGTTGCATGCAAGCTTGGCACTG";
        let template = format!(
            "TTTTTTTTTT{}AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
            &primer[..18]
        );
        let m = mol(&template, Topology::Linear);
        let c = crate::tm::TmConditions::default();
        let on = sites_with_prefs(&m, &[primer.as_bytes()], &BindingPrefs::default(), &c);
        assert!(on.iter().any(|s| s.five_prime_only), "{on:?}");
        // Its Tm is the paired 5' stretch's.
        let site = on.iter().find(|s| s.five_prime_only).unwrap();
        assert_eq!(
            site.tm(primer.as_bytes(), &c),
            crate::tm::tm(&primer.as_bytes()[..18], &c)
        );
        let off = BindingPrefs {
            five_prime_on: false,
            ..BindingPrefs::default()
        };
        assert!(
            !sites_with_prefs(&m, &[primer.as_bytes()], &off, &c)
                .iter()
                .any(|s| s.five_prime_only)
        );
        // Where the primer pairs in full, only the full site is listed.
        let m = mol(&format!("TTTTTTTTTT{primer}TTTTTTTTTT"), Topology::Linear);
        let full = sites_with_prefs(&m, &[primer.as_bytes()], &BindingPrefs::default(), &c);
        assert_eq!(full.len(), 1);
        assert!(!full[0].five_prime_only);
    }
}
