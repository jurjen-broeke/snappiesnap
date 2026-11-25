//! Diagnostic digests: enzymes, alone or in pairs, whose fragments tell constructs apart on an
//! agarose gel (a correct clone from its empty vector, one insert orientation from the other)
//! or, for a single construct, give a pattern that confirms it.
//!
//! Fragments come from the same cut positions as [`crate::enzymes::digest`]; where each one
//! runs comes from the gel model ([`crate::gel::migration`]), so a suggestion looks on the
//! simulated gel the way it was judged. The ranking is our own heuristic, stated here so users
//! can follow it:
//! - every construct must be cut, into at most `max_bands` visible fragments (fragments under
//!   [`MIN_VISIBLE`] bp run off the gel and don't count);
//! - fragments in one lane should run apart (closer than [`MIN_GAP`] of the gel counts as one
//!   band) and inside the gel's separation range;
//! - every two constructs must differ by a band at least [`MIN_DIFFERENCE`] of the gel away
//!   from any band of the other: the larger the smallest such difference, the better;
//! - 2 to 5 bands per lane read best; one enzyme beats two (no shared buffer to find);
//! - of enzymes giving the same pattern, the most specific site wins (its information content,
//!   2 bits per A/C/G/T, less for ambiguity codes): it behaves most predictably on real DNA.
//!   Then the order of `enzymes`, so callers can put the names their users know first.

use std::collections::BTreeMap;

use serde::Serialize;
use snappiesnap_doc::{Molecule, Topology};

use crate::enzymes::{Enzyme, find_sites};
use crate::gel::{migration, separation_range};

/// Bands closer than this (fraction of the gel's length) read as one.
pub const MIN_GAP: f64 = 0.02;
/// Two constructs differ visibly when one has a band this far (fraction of the gel's length)
/// from every band of the other.
pub const MIN_DIFFERENCE: f64 = 0.04;
/// Fragments shorter than this (bp) run off the gel or are too faint to see: they don't count
/// as bands.
pub const MIN_VISIBLE: usize = 50;
/// Most constructs compared at once.
pub const MAX_CONSTRUCTS: usize = 12;
/// Most enzymes tried in pairs (the ones cutting least often are kept).
const MAX_PAIRED: usize = 150;

/// What to look for.
#[derive(Debug, Clone, PartialEq)]
pub struct Options {
    /// Agarose, % (sets the separation range).
    pub agarose: f64,
    /// 1: single enzymes only; 2: also pairs (double digests).
    pub max_enzymes: usize,
    /// Most fragments per construct.
    pub max_bands: usize,
    /// Most suggestions returned.
    pub limit: usize,
}

impl Default for Options {
    fn default() -> Self {
        Options {
            agarose: 1.0,
            max_enzymes: 2,
            max_bands: 6,
            limit: 10,
        }
    }
}

/// One suggestion.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct Suggestion {
    /// The enzymes to use together.
    pub enzymes: Vec<String>,
    /// Enzymes that cut the same way as each of `enzymes` (isoschizomers with the same cuts),
    /// in the same order.
    pub alternatives: Vec<Vec<String>>,
    /// Per construct, in the order given: fragment sizes (bp), largest first.
    pub fragments: Vec<Vec<usize>>,
    /// The smallest difference between any two constructs (fraction of the gel's length);
    /// `None` for a single construct.
    pub separation: Option<f64>,
    /// The smallest gap between two bands of one lane (fraction of the gel's length).
    pub closest_bands: Option<f64>,
    pub score: f64,
    /// What to watch out for.
    pub notes: Vec<String>,
}

/// Information content of a recognition site in bits: log2(4 / k) per position, where k is
/// how many bases the IUPAC code there stands for (A: 2 bits, R: 1 bit, N: 0).
fn specificity(site: &str) -> f64 {
    site.bytes()
        .map(|b| {
            let k: f64 = match b.to_ascii_uppercase() {
                b'A' | b'C' | b'G' | b'T' => 1.0,
                b'R' | b'Y' | b'S' | b'W' | b'K' | b'M' => 2.0,
                b'B' | b'D' | b'H' | b'V' => 3.0,
                _ => 4.0,
            };
            (4.0 / k).log2()
        })
        .sum()
}

/// Enzymes that cut the same way share one entry.
fn cut_key(e: &Enzyme) -> String {
    format!(
        "{}|{:?}|{:?}|{:?}",
        e.site, e.cut_top, e.cut_bottom, e.cut_upstream
    )
}

/// Sorted, distinct top-strand cuts of each enzyme group in `m` (inside a linear molecule).
fn cuts(m: &Molecule, enzymes: &[Enzyme]) -> Vec<Vec<usize>> {
    let len = m.len();
    let circular = m.topology() == Topology::Circular;
    let mut out = vec![Vec::new(); enzymes.len()];
    for s in find_sites(m, enzymes) {
        if let Some(list) = out.get_mut(s.enzyme) {
            for c in [s.cut, s.cut_upstream].into_iter().flatten() {
                if circular || (c > 0 && c < len) {
                    list.push(c);
                }
            }
        }
    }
    for list in &mut out {
        list.sort_unstable();
        list.dedup();
    }
    out
}

/// Fragment sizes, largest first, from sorted distinct cuts.
fn fragments(len: usize, circular: bool, cuts: &[usize]) -> Vec<usize> {
    let mut out = Vec::with_capacity(cuts.len() + 1);
    match (cuts.first(), cuts.last()) {
        (None, _) | (_, None) => out.push(len),
        (Some(&first), Some(&last)) => {
            for w in cuts.windows(2) {
                if let [a, b] = w {
                    out.push(b - a);
                }
            }
            if circular {
                out.push(len - last + first);
            } else {
                out.push(first);
                out.push(len - last);
            }
        }
    }
    out.retain(|&f| f > 0);
    out.sort_unstable_by(|a, b| b.cmp(a));
    out
}

fn merge(a: &[usize], b: &[usize]) -> Vec<usize> {
    let mut v: Vec<usize> = a.iter().chain(b).copied().collect();
    v.sort_unstable();
    v.dedup();
    v
}

/// The smallest distance from band `y` to any of `others`.
fn nearest(y: f64, others: &[f64]) -> f64 {
    others
        .iter()
        .map(|o| (o - y).abs())
        .fold(f64::INFINITY, f64::min)
}

/// How a combination did: its score, the smallest difference between constructs and between
/// bands of one lane, and notes.
struct Judgement {
    score: f64,
    separation: Option<f64>,
    closest: Option<f64>,
    notes: Vec<String>,
}

/// Judge one enzyme combination; `None` if it breaks a hard rule.
fn judge(lanes: &[Vec<usize>], opts: &Options, pair: bool) -> Option<Judgement> {
    let (lo, hi) = separation_range(opts.agarose);
    let mut notes = Vec::new();
    let mut score = 0.0;
    let tiny = lanes
        .iter()
        .flatten()
        .filter(|&&bp| bp < MIN_VISIBLE)
        .count();
    let lanes: Vec<Vec<usize>> = lanes
        .iter()
        .map(|l| l.iter().copied().filter(|&bp| bp >= MIN_VISIBLE).collect())
        .collect();
    if lanes.iter().any(Vec::is_empty) {
        return None;
    }
    if tiny > 0 {
        notes.push(format!(
            "{tiny} fragment(s) under {MIN_VISIBLE} bp run off the gel and aren't shown"
        ));
    }
    let ys: Vec<Vec<f64>> = lanes
        .iter()
        .map(|l| {
            l.iter()
                .map(|&bp| migration(bp as f64, opts.agarose))
                .collect()
        })
        .collect();
    let mut closest: Option<f64> = None;
    for (l, y) in lanes.iter().zip(&ys) {
        if l.len() > opts.max_bands {
            return None;
        }
        // Bands are sorted largest first, so migration rises along the lane.
        for w in y.windows(2) {
            if let [a, b] = w {
                let gap = b - a;
                closest = Some(closest.map_or(gap, |c: f64| c.min(gap)));
            }
        }
        score -= (l.len() as f64 - 3.5).abs().max(1.5) - 1.5;
        let outside = l
            .iter()
            .filter(|&&bp| (bp as f64) < lo || (bp as f64) > hi)
            .count();
        score -= outside as f64;
    }
    if let Some(c) = closest {
        if c < MIN_GAP {
            score -= 3.0;
            notes.push("some fragments of one lane run together".into());
        } else {
            score += (c / MIN_GAP).min(3.0);
        }
    }
    let mut separation: Option<f64> = None;
    for i in 0..ys.len() {
        for j in i + 1..ys.len() {
            let (Some(a), Some(b)) = (ys.get(i), ys.get(j)) else {
                continue;
            };
            let d = a
                .iter()
                .map(|&y| nearest(y, b))
                .chain(b.iter().map(|&y| nearest(y, a)))
                .fold(0.0, f64::max);
            separation = Some(separation.map_or(d, |s: f64| s.min(d)));
        }
    }
    if let Some(s) = separation {
        if s < MIN_DIFFERENCE {
            return None;
        }
        score += 10.0 * (s / MIN_DIFFERENCE).min(4.0);
    } else if lanes.first().is_some_and(|l| l.len() < 2) {
        // A single construct cut once confirms little.
        score -= 2.0;
    }
    if pair {
        score -= 2.0;
        notes.push("double digest: choose a buffer both enzymes work in".into());
    }
    let smallest = lanes.iter().flatten().copied().min().unwrap_or(0);
    if (smallest as f64) < lo {
        notes.push(format!(
            "a {smallest} bp fragment runs below the gel's range and may be faint"
        ));
    }
    Some(Judgement {
        score,
        separation,
        closest,
        notes,
    })
}

/// Suggestions for `molecules` (each DNA; 1 to [`MAX_CONSTRUCTS`]), best first.
pub fn design(
    molecules: &[&Molecule],
    enzymes: &[Enzyme],
    opts: &Options,
) -> Result<Vec<Suggestion>, String> {
    if molecules.is_empty() || molecules.len() > MAX_CONSTRUCTS {
        return Err(format!("give 1 to {MAX_CONSTRUCTS} constructs"));
    }
    if let Some(m) = molecules
        .iter()
        .find(|m| !m.seq().alphabet().is_nucleotide() || m.is_empty())
    {
        return Err(format!("{} isn't a DNA sequence", m.name()));
    }
    if !(1..=2).contains(&opts.max_enzymes) {
        return Err("use 1 or 2 enzymes per digest".into());
    }
    if !(1..=30).contains(&opts.max_bands) {
        return Err("allow 1 to 30 bands per lane".into());
    }
    // One representative per way of cutting; the rest are alternatives.
    let mut groups: BTreeMap<String, usize> = BTreeMap::new();
    let mut reps: Vec<Enzyme> = Vec::new();
    let mut alts: Vec<Vec<String>> = Vec::new();
    for e in enzymes.iter().filter(|e| e.cut_top.is_some()) {
        match groups.get(&cut_key(e)) {
            Some(&i) => {
                if let Some(a) = alts.get_mut(i) {
                    a.push(e.name.clone());
                }
            }
            None => {
                groups.insert(cut_key(e), reps.len());
                reps.push(e.clone());
                alts.push(Vec::new());
            }
        }
    }
    // cuts[m][e]
    let all: Vec<Vec<Vec<usize>>> = molecules.iter().map(|m| cuts(m, &reps)).collect();
    let lanes_for = |set: &[usize]| -> Option<Vec<Vec<usize>>> {
        molecules
            .iter()
            .zip(&all)
            .map(|(m, per)| {
                let mut c: Vec<usize> = Vec::new();
                for &e in set {
                    c = merge(&c, per.get(e)?);
                }
                (!c.is_empty()).then(|| fragments(m.len(), m.topology() == Topology::Circular, &c))
            })
            .collect()
    };
    // With each suggestion: its sites' specificity and its enzymes' places in `enzymes`.
    let mut found: Vec<(Suggestion, f64, Vec<usize>)> = Vec::new();
    let mut consider = |set: &[usize]| {
        let Some(lanes) = lanes_for(set) else {
            return;
        };
        if let Some(Judgement {
            score,
            separation,
            closest,
            notes,
        }) = judge(&lanes, opts, set.len() > 1)
        {
            let spec = set
                .iter()
                .filter_map(|&i| reps.get(i).map(|e| specificity(&e.site)))
                .sum::<f64>();
            let found_one = Suggestion {
                enzymes: set
                    .iter()
                    .filter_map(|&i| reps.get(i).map(|e| e.name.clone()))
                    .collect(),
                alternatives: set
                    .iter()
                    .map(|&i| alts.get(i).cloned().unwrap_or_default())
                    .collect(),
                fragments: lanes,
                separation,
                closest_bands: closest,
                score,
                notes,
            };
            found.push((found_one, spec, set.to_vec()));
        }
    };
    let cut_count = |e: usize| -> usize {
        all.iter()
            .map(|per| per.get(e).map_or(0, Vec::len))
            .max()
            .unwrap_or(0)
    };
    for e in 0..reps.len() {
        consider(&[e]);
    }
    if opts.max_enzymes >= 2 {
        // Enzymes that cut no construct add nothing to a pair.
        let mut paired: Vec<usize> = (0..reps.len())
            .filter(|&e| (1..opts.max_bands).contains(&cut_count(e)))
            .collect();
        paired.sort_by_key(|&e| cut_count(e));
        paired.truncate(MAX_PAIRED);
        paired.sort_unstable();
        for (k, &a) in paired.iter().enumerate() {
            for &b in paired.iter().skip(k + 1) {
                consider(&[a, b]);
            }
        }
    }
    // Best first; ties by fewer enzymes, more specific sites, the caller's order, so results
    // are the same on every run.
    found.sort_by(|(x, xs, xi), (y, ys, yi)| {
        y.score
            .total_cmp(&x.score)
            .then(x.enzymes.len().cmp(&y.enzymes.len()))
            .then(ys.total_cmp(xs))
            .then(xi.cmp(yi))
    });
    // Different enzymes giving the same pattern add nothing: keep the best of each.
    let mut seen: Vec<Vec<Vec<usize>>> = Vec::new();
    let mut out = Vec::new();
    for (s, _, _) in found {
        if !seen.contains(&s.fragments) {
            seen.push(s.fragments.clone());
            out.push(s);
        }
        if out.len() >= opts.limit {
            break;
        }
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use snappiesnap_doc::{Alphabet, Seq};

    fn mol(name: &str, residues: &str, topology: Topology) -> Molecule {
        Molecule::new(name, Seq::new(Alphabet::Dna, residues).unwrap(), topology).unwrap()
    }

    fn enzymes() -> Vec<Enzyme> {
        vec![
            Enzyme::parse("EcoRI", "G^AATTC", true).unwrap(),
            Enzyme::parse("BamHI", "G^GATCC", true).unwrap(),
            Enzyme::parse("HindIII", "A^AGCTT", true).unwrap(),
            // Same site and cut as EcoRI: an alternative, not a second suggestion.
            Enzyme::parse("EcoRI-twin", "G^AATTC", true).unwrap(),
        ]
    }

    /// A filler without any of the test sites.
    fn filler(n: usize) -> String {
        "ACGTTGCA".repeat(n / 8 + 1)[..n].to_string()
    }

    #[test]
    fn site_information() {
        assert_eq!(specificity("GAATTC"), 12.0);
        assert_eq!(specificity("RAATTY"), 10.0);
        assert_eq!(specificity("AGCT"), 8.0);
        assert_eq!(specificity("GCNNNNNNNGC"), 8.0);
    }

    #[test]
    fn fragment_sizes_by_hand() {
        assert_eq!(fragments(100, false, &[]), vec![100]);
        assert_eq!(fragments(100, false, &[30, 80]), vec![50, 30, 20]);
        assert_eq!(fragments(100, true, &[30]), vec![100]);
        assert_eq!(fragments(100, true, &[30, 80]), vec![50, 50]);
        assert_eq!(fragments(100, true, &[10, 30, 80]), vec![50, 30, 20]);
    }

    #[test]
    fn tells_a_clone_from_its_empty_vector() {
        // The vector has EcoRI and HindIII sites 3,000 bp apart around the cloning site; the
        // clone has a 1,500 bp insert with a BamHI site between them.
        let vector = format!("GAATTC{}AAGCTT{}", filler(2994), filler(2988));
        let clone = format!(
            "GAATTC{}GGATCC{}AAGCTT{}",
            filler(2994),
            filler(1494),
            filler(2988)
        );
        let v = mol("vector", &vector, Topology::Circular);
        let c = mol("clone", &clone, Topology::Circular);
        let found = design(&[&v, &c], &enzymes(), &Options::default()).unwrap();
        assert!(!found.is_empty());
        let best = &found[0];
        assert!(best.separation.unwrap() >= MIN_DIFFERENCE, "{best:?}");
        // EcoRI + HindIII cut both (cuts after the first base of each site): the 5,994 bp vector
        // into 3,000 + 2,994 bp, the 7,494 bp clone into 4,500 + 2,994 bp.
        let pair = found
            .iter()
            .find(|s| s.enzymes == ["EcoRI", "HindIII"])
            .expect("EcoRI + HindIII");
        assert_eq!(pair.fragments, vec![vec![3000, 2994], vec![4500, 2994]]);
        assert_eq!(
            pair.alternatives,
            vec![vec!["EcoRI-twin".to_string()], vec![]]
        );
        assert!(pair.notes.iter().any(|n| n.contains("buffer")));
        // BamHI alone leaves the vector uncut: not suggested.
        assert!(!found.iter().any(|s| s.enzymes == ["BamHI"]));
        // The same pattern twice is listed once.
        let mut patterns: Vec<_> = found.iter().map(|s| &s.fragments).collect();
        patterns.dedup();
        assert_eq!(patterns.len(), found.len());
    }

    #[test]
    fn pairs_are_found_among_many_enzymes_that_dont_cut() {
        // 300 enzymes with sites absent from the test sequences must not crowd out the pair.
        let mut many = enzymes();
        for i in 0..300usize {
            let site: String = (0..8)
                .map(|k| ['C', 'G'][(i >> k) & 1])
                .chain("TTAACC".chars())
                .collect();
            many.push(Enzyme::parse(&format!("Fake{i}"), &format!("{site}^"), true).unwrap());
        }
        let vector = format!("GAATTC{}AAGCTT{}", filler(2994), filler(2988));
        let clone = format!(
            "GAATTC{}GGATCC{}AAGCTT{}",
            filler(2994),
            filler(1494),
            filler(2988)
        );
        let v = mol("vector", &vector, Topology::Circular);
        let c = mol("clone", &clone, Topology::Circular);
        let opts = Options {
            limit: 50,
            ..Options::default()
        };
        let found = design(&[&v, &c], &many, &opts).unwrap();
        assert!(
            found.iter().any(|s| s.enzymes == ["EcoRI", "HindIII"]),
            "{found:?}"
        );
    }

    #[test]
    fn tiny_fragments_are_not_bands() {
        // Two EcoRI sites 20 bp apart: the 20 bp piece runs off the gel.
        let s = format!("GAATTC{}GAATTC{}", filler(14), filler(3000));
        let a = mol("a", &s, Topology::Circular);
        let found = design(&[&a], &enzymes()[..1], &Options::default()).unwrap();
        let ecori = found.first().unwrap();
        assert_eq!(ecori.fragments, vec![vec![3006, 20]]);
        assert!(
            ecori.notes.iter().any(|n| n.contains("under 50 bp")),
            "{ecori:?}"
        );
    }

    #[test]
    fn identical_constructs_cannot_be_told_apart() {
        let s = format!("GAATTC{}AAGCTT{}", filler(2000), filler(2000));
        let a = mol("a", &s, Topology::Circular);
        let b = mol("b", &s, Topology::Circular);
        assert!(
            design(&[&a, &b], &enzymes(), &Options::default())
                .unwrap()
                .is_empty()
        );
    }

    #[test]
    fn single_constructs_and_limits() {
        let s = format!(
            "GAATTC{}AAGCTT{}GGATCC{}",
            filler(1000),
            filler(2000),
            filler(500)
        );
        let a = mol("a", &s, Topology::Linear);
        let found = design(&[&a], &enzymes(), &Options::default()).unwrap();
        assert!(found.iter().all(|f| f.separation.is_none()));
        let singles = design(
            &[&a],
            &enzymes(),
            &Options {
                max_enzymes: 1,
                ..Options::default()
            },
        )
        .unwrap();
        assert!(singles.iter().all(|f| f.enzymes.len() == 1));
        assert!(design(&[], &enzymes(), &Options::default()).is_err());
        let p = Molecule::new(
            "p",
            Seq::new(Alphabet::Protein, "MKV").unwrap(),
            Topology::Linear,
        )
        .unwrap();
        assert!(design(&[&p], &enzymes(), &Options::default()).is_err());
        let e = mol("e", "", Topology::Linear);
        assert!(design(&[&e], &enzymes(), &Options::default()).is_err());
        assert!(design(&[&a], &[], &Options::default()).unwrap().is_empty());
    }
}
