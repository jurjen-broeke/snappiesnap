//! PCR primer design: a forward and a reverse primer that amplify a region.
//!
//! Candidates are scored with the usual design rules (Dieffenbach C. W., Lowe T. M. J.,
//! Dveksler G. S. (1993) "General concepts for PCR primer design", PCR Methods Appl.
//! 3:S30-S37): a Tm near the target and close to the partner's, 40-60 % GC, a G or C at the 3'
//! end but no more than three G/C in the last five bases, no long single-base runs or
//! dinucleotide repeats, no stable hairpins or primer-dimers (free energies from
//! [`crate::structure`]), and a single binding site on the template. Tm is the nearest-neighbour value of
//! [`crate::tm`]. Each rule adds a penalty; lower is better.

use serde::{Deserialize, Serialize};
use snappiesnap_doc::{Molecule, Range, Topology};

use crate::primers::{BindingRules, binding_sites};
use crate::structure;
use crate::tm::{TmConditions, gc_fraction, tm};

/// What to design.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(default)]
pub struct DesignOptions {
    pub target_tm: f64,
    pub min_len: usize,
    pub max_len: usize,
    /// How far outside the region a primer may start (0: the product is exactly the region).
    pub search: usize,
    /// Pairs to return.
    pub count: usize,
    pub conditions: TmConditions,
}

impl Default for DesignOptions {
    fn default() -> Self {
        DesignOptions {
            target_tm: 60.0,
            min_len: 18,
            max_len: 30,
            search: 0,
            count: 5,
            conditions: TmConditions::default(),
        }
    }
}

/// One designed primer.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct Designed {
    /// 5'->3'.
    pub seq: String,
    /// Where it anneals (top-strand coordinates; wraps on circular molecules).
    pub range: Range,
    pub tm: f64,
    pub gc_percent: f64,
    /// Rules it doesn't meet, in words.
    pub warnings: Vec<String>,
    pub penalty: f64,
}

/// A forward and reverse primer and the product they make.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct PrimerPair {
    pub forward: Designed,
    pub reverse: Designed,
    /// Product length, bp.
    pub product: usize,
    /// Pair penalty (both primers, Tm difference, primer-dimer); lower is better.
    pub penalty: f64,
    pub warnings: Vec<String>,
}

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum DesignError {
    #[error("{0}")]
    Invalid(String),
    #[error("no primers meet the length limits for this region")]
    NoCandidates,
}

const MIN_PRIMER: usize = 5;
const MAX_PRIMER: usize = 60;
const MAX_SEARCH: usize = 1000;
const MAX_COUNT: usize = 50;
/// Candidates kept per side before pairing.
const PER_SIDE: usize = 25;

fn complement(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        other => other,
    }
}

fn revcomp(s: &[u8]) -> Vec<u8> {
    s.iter().rev().map(|&b| complement(b)).collect()
}

/// Penalty for a structure ΔG beyond its warning limit (kcal/mol, both negative).
fn structure_penalty(dg: f64, limit: f64) -> f64 {
    if dg <= limit {
        2.0 * (limit - dg) + 1.0
    } else {
        0.0
    }
}

fn longest_run(s: &[u8]) -> usize {
    let mut best = 0;
    let mut run = 0;
    let mut prev = None;
    for &b in s {
        run = if Some(b) == prev { run + 1 } else { 1 };
        prev = Some(b);
        best = best.max(run);
    }
    best
}

fn longest_dinucleotide_repeat(s: &[u8]) -> usize {
    let mut best = 0;
    for start in 0..2 {
        let mut count = 0;
        let mut prev: Option<&[u8]> = None;
        let mut i = start;
        while let Some(pair) = s.get(i..i + 2) {
            count = if Some(pair) == prev && pair.first() != pair.get(1) {
                count + 1
            } else {
                1
            };
            best = best.max(count);
            prev = Some(pair);
            i += 2;
        }
    }
    best
}

/// Score one primer (5'->3') on its own. `None` if it has no Tm (ambiguous bases).
fn score(seq: &[u8], opts: &DesignOptions) -> Option<Designed> {
    let t = tm(seq, &opts.conditions)?;
    let gc = gc_fraction(seq)? * 100.0;
    let mut warnings = Vec::new();
    let mut penalty = (t - opts.target_tm).abs();
    if !(40.0..=60.0).contains(&gc) {
        penalty += (gc - 50.0).abs() / 5.0;
        warnings.push(format!("GC {gc:.0} % (aim for 40-60 %)"));
    }
    let last5 = seq.get(seq.len().saturating_sub(5)..).unwrap_or_default();
    let gc3 = last5.iter().filter(|b| matches!(b, b'G' | b'C')).count();
    if !matches!(seq.last(), Some(b'G' | b'C')) {
        penalty += 1.0;
        warnings.push("no G or C at the 3' end".into());
    }
    if gc3 > 3 {
        penalty += 2.0 * (gc3 - 3) as f64;
        warnings.push(format!(
            "{gc3} G/C in the last 5 bases (3' end may prime falsely)"
        ));
    }
    let run = longest_run(seq);
    if run >= 5 {
        penalty += 2.0 * (run - 4) as f64;
        warnings.push(format!("run of {run} identical bases"));
    }
    let di = longest_dinucleotide_repeat(seq);
    if di >= 5 {
        penalty += 2.0 * (di - 4) as f64;
        warnings.push(format!("{di} dinucleotide repeats"));
    }
    // Hairpins and self-dimers by free energy (see `crate::structure`).
    let c = &opts.conditions;
    if let Some(h) = structure::hairpin(seq, c) {
        let limit = if h.three_prime {
            structure::HAIRPIN_3PRIME_WARN
        } else {
            structure::HAIRPIN_WARN
        };
        penalty += structure_penalty(h.dg, limit);
    }
    if let Some((any, three)) = structure::dimers(seq, seq, c) {
        let p3 = three.map_or(0.0, |d| {
            structure_penalty(d.dg, structure::DIMER_3PRIME_WARN)
        });
        let p = any.map_or(0.0, |d| structure_penalty(d.dg, structure::DIMER_WARN));
        penalty += p3.max(p);
    }
    warnings.extend(structure::warnings(seq, c));
    Some(Designed {
        seq: String::from_utf8_lossy(seq).into_owned(),
        range: Range::new(0, 0),
        tm: t,
        gc_percent: gc,
        warnings,
        penalty,
    })
}

/// Design primer pairs that amplify `region` of `m`, best first.
pub fn design(
    m: &Molecule,
    region: Range,
    opts: &DesignOptions,
) -> Result<Vec<PrimerPair>, DesignError> {
    let invalid = |s: &str| DesignError::Invalid(s.to_string());
    if !m.seq().alphabet().is_nucleotide() {
        return Err(invalid("primers can only be designed on DNA"));
    }
    if !(MIN_PRIMER..=MAX_PRIMER).contains(&opts.min_len)
        || !(opts.min_len..=MAX_PRIMER).contains(&opts.max_len)
    {
        return Err(invalid("primer lengths must be within 5-60 and min <= max"));
    }
    if opts.search > MAX_SEARCH || opts.count == 0 || opts.count > MAX_COUNT {
        return Err(invalid("search must be at most 1000 and count 1-50"));
    }
    if !opts.target_tm.is_finite() || !(20.0..=90.0).contains(&opts.target_tm) {
        return Err(invalid("target Tm must be 20-90 °C"));
    }
    let len = m.len();
    let circular = m.topology() == Topology::Circular;
    let region = region
        .validate(len, m.topology())
        .map_err(|e| DesignError::Invalid(e.to_string()))?;
    let span = region.len_on(len);
    if span < 2 * opts.min_len {
        return Err(invalid("the region is shorter than two primers"));
    }
    let seq: Vec<u8> = m
        .seq()
        .as_bytes()
        .iter()
        .map(u8::to_ascii_uppercase)
        .collect();
    let l = len as i64;
    // Residues by unwrapped coordinate (relative to the region start); None outside a linear
    // molecule.
    let start = region.start as i64;
    let at = |i: i64| -> Option<u8> {
        let p = start + i;
        let p = if circular { p.rem_euclid(l) } else { p };
        usize::try_from(p).ok().and_then(|p| seq.get(p).copied())
    };
    let stretch =
        |from: i64, n: usize| -> Option<Vec<u8>> { (0..n as i64).map(|k| at(from + k)).collect() };
    let to_range = |from: i64, n: usize| -> Range {
        let s = (start + from).rem_euclid(l.max(1)) as usize;
        let e = (start + from + n as i64).rem_euclid(l.max(1)) as usize;
        Range::new(s, if e == 0 && !circular { len } else { e })
    };
    let search = opts.search as i64;
    let span_i = span as i64;
    let mut fwd = Vec::new();
    let mut rev = Vec::new();
    for shift in 0..=search {
        for n in opts.min_len..=opts.max_len {
            // Forward: 5' end `shift` bases before the region start.
            if let Some(s) = stretch(-shift, n)
                && let Some(mut d) = score(&s, opts)
            {
                d.range = to_range(-shift, n);
                d.penalty += shift as f64 * 0.01;
                fwd.push((shift, d));
            }
            // Reverse: 5' end `shift` bases after the region end.
            let from = span_i + shift - n as i64;
            if let Some(s) = stretch(from, n)
                && let Some(mut d) = score(&revcomp(&s), opts)
            {
                d.range = to_range(from, n);
                d.penalty += shift as f64 * 0.01;
                rev.push((shift, d));
            }
        }
    }
    let best = |v: &mut Vec<(i64, Designed)>| {
        v.sort_by(|a, b| a.1.penalty.total_cmp(&b.1.penalty));
        v.truncate(PER_SIDE);
    };
    best(&mut fwd);
    best(&mut rev);
    if fwd.is_empty() || rev.is_empty() {
        return Err(DesignError::NoCandidates);
    }
    let mut pairs = Vec::new();
    for (fs, f) in &fwd {
        for (rs, r) in &rev {
            let product = span + (*fs + *rs) as usize;
            if circular && product > len {
                continue;
            }
            let mut warnings = Vec::new();
            let dt = (f.tm - r.tm).abs();
            let mut penalty = f.penalty + r.penalty + dt;
            if dt > 5.0 {
                warnings.push(format!("Tm differs by {dt:.1} °C"));
            }
            let (fs, rs) = (f.seq.as_bytes(), r.seq.as_bytes());
            if let Some((any, three)) = structure::dimers(fs, rs, &opts.conditions) {
                let p3 = three.map_or(0.0, |d| {
                    structure_penalty(d.dg, structure::DIMER_3PRIME_WARN)
                });
                let p = any.map_or(0.0, |d| structure_penalty(d.dg, structure::DIMER_WARN));
                penalty += p3.max(p);
            }
            warnings.extend(structure::pair_warnings(fs, rs, &opts.conditions));
            pairs.push(PrimerPair {
                forward: f.clone(),
                reverse: r.clone(),
                product,
                penalty,
                warnings,
            });
        }
    }
    pairs.sort_by(|a, b| a.penalty.total_cmp(&b.penalty));
    pairs.truncate(opts.count);
    // Primers that also bind elsewhere can make extra products: strong sites (exact over 15+
    // bases) count heavily, weak ones a little.
    for p in &mut pairs {
        let seqs = [p.forward.seq.as_bytes(), p.reverse.seq.as_bytes()];
        let sites = binding_sites(m, &seqs, BindingRules::default());
        for (i, d) in [&mut p.forward, &mut p.reverse].into_iter().enumerate() {
            let n = d.seq.len();
            let mine: Vec<_> = sites.iter().filter(|s| s.primer == i).collect();
            let strong = mine
                .iter()
                .filter(|s| s.is_primary(n))
                .count()
                .saturating_sub(1);
            let weak = mine.iter().filter(|s| !s.is_primary(n)).count();
            let extra = 5.0 * strong as f64 + 0.5 * weak as f64;
            if strong > 0 {
                d.warnings
                    .push(format!("binds {strong} other place(s) on the template"));
            }
            if weak > 0 {
                d.warnings.push(format!("{weak} weak secondary site(s)"));
            }
            d.penalty += extra;
            p.penalty += extra;
        }
    }
    pairs.sort_by(|a, b| a.penalty.total_cmp(&b.penalty));
    Ok(pairs)
}

#[cfg(test)]
mod tests {
    use snappiesnap_doc::{Alphabet, Seq};

    use super::*;

    /// A made-up, varied 300 bp sequence.
    fn template() -> String {
        let mut x: u32 = 12345;
        (0..300)
            .map(|_| {
                x = x.wrapping_mul(1_103_515_245).wrapping_add(12345);
                b"ACGT"[((x >> 16) % 4) as usize] as char
            })
            .collect()
    }

    fn mol(s: &str, t: Topology) -> Molecule {
        Molecule::new("m", Seq::new(Alphabet::Dna, s).unwrap(), t).unwrap()
    }

    #[test]
    fn exact_region_primers() {
        let s = template();
        let m = mol(&s, Topology::Linear);
        let pairs = design(&m, Range::new(50, 250), &DesignOptions::default()).unwrap();
        assert!(!pairs.is_empty());
        let p = &pairs[0];
        assert_eq!(p.product, 200);
        // The forward primer is the top strand from the region start; the reverse primer is
        // the reverse complement ending at the region end.
        assert!(s[50..].starts_with(&p.forward.seq));
        assert_eq!(p.forward.range.start, 50);
        let rc: String = revcomp(p.reverse.seq.as_bytes())
            .iter()
            .map(|&b| b as char)
            .collect();
        assert!(s[..250].ends_with(&rc));
        assert_eq!(p.reverse.range.end, 250);
        for d in [&p.forward, &p.reverse] {
            assert!((18..=30).contains(&d.seq.len()));
            assert!((d.tm - 60.0).abs() < 8.0, "{}", d.tm);
        }
        assert!(pairs.windows(2).all(|w| w[0].penalty <= w[1].penalty));
    }

    #[test]
    fn search_margin_and_circular_wrap() {
        let s = template();
        let m = mol(&s, Topology::Circular);
        let opts = DesignOptions {
            search: 20,
            ..DesignOptions::default()
        };
        // A region across the origin.
        let pairs = design(&m, Range::new(250, 60), &opts).unwrap();
        for p in &pairs {
            assert!((110..=150).contains(&p.product), "{}", p.product);
        }
    }

    #[test]
    fn rules() {
        assert_eq!(longest_run(b"ACGTTTTTAC"), 5);
        assert_eq!(longest_dinucleotide_repeat(b"GCATATATATATGC"), 5);
        let d = score(b"ATATATATATATATATATAT", &DesignOptions::default()).unwrap();
        assert!(d.warnings.iter().any(|w| w.contains("GC")));
        assert!(d.warnings.iter().any(|w| w.contains("dinucleotide")));
    }

    #[test]
    fn bad_requests() {
        let m = mol(&template(), Topology::Linear);
        let o = DesignOptions::default();
        for (region, opts) in [
            (Range::new(0, 20), o),
            (Range::new(250, 50), o),
            (Range::new(0, 400), o),
            (
                Range::new(0, 200),
                DesignOptions {
                    min_len: 40,
                    max_len: 20,
                    ..o
                },
            ),
            (Range::new(0, 200), DesignOptions { count: 0, ..o }),
            (
                Range::new(0, 200),
                DesignOptions {
                    target_tm: f64::NAN,
                    ..o
                },
            ),
            (
                Range::new(0, 200),
                DesignOptions {
                    search: 1_000_000,
                    ..o
                },
            ),
        ] {
            assert!(design(&m, region, &opts).is_err(), "{region:?} {opts:?}");
        }
        // A primer can't start before a linear molecule: the region at the very start still
        // works without a search margin.
        assert!(design(&m, Range::new(0, 300), &o).is_ok());
    }
}
