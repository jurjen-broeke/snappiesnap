//! Restriction enzymes: recognition sites, cut positions, site finding and digests.
//!
//! Sites use the conventional notation of the enzyme literature and supplier catalogues:
//! IUPAC letters 5'→3' on one strand, with the top-strand cut marked by `^` (`G^AATTC`), or for
//! enzymes that cut outside their site, `(n/m)` after it: the top strand is cut `n` bases and
//! the bottom strand `m` bases past the site's 3' end, in top-strand coordinates
//! (`GGTCTC(1/5)`). Negative values cut inside or before the site. Enzymes that cut on both
//! sides of their site also have a pair before it, counted upstream of the site's 5' end:
//! `(8/13)GACNNNNNNTGG(12/7)`.

mod builtin;
mod rebase;

pub use builtin::{builtin, builtin_suppliers};
pub use rebase::{Rebase, RebaseError, parse_withrefm};

use serde::Serialize;
use snappiesnap_doc::{Molecule, Range, Topology};

/// Why an enzyme definition couldn't be read.
#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
#[error("enzyme {name}: {reason}")]
pub struct EnzymeError {
    pub name: String,
    pub reason: String,
}

/// Longest recognition site we accept.
const MAX_SITE: usize = 64;
/// Furthest a cut may be from its site.
const MAX_OFFSET: i32 = 1000;

/// A restriction enzyme.
#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
pub struct Enzyme {
    pub name: String,
    /// Recognition site, upper-case IUPAC, 5'→3' on the top strand.
    pub site: String,
    /// Where the top strand is cut, as a count of positions from the site's 5' end (may be
    /// negative or beyond the site). `None` if the cut position is unknown.
    pub cut_top: Option<i32>,
    /// The same for the bottom strand, in top-strand coordinates.
    pub cut_bottom: Option<i32>,
    /// For enzymes that cut on both sides: the upstream (top, bottom) cut, in the same
    /// coordinates (negative).
    pub cut_upstream: Option<(i32, i32)>,
    /// Sold by at least one supplier.
    pub commercial: bool,
}

impl Enzyme {
    /// Parse an enzyme from its name and site notation (`G^AATTC`, `GGTCTC(1/5)`, `GCANNNN^NTGC`).
    pub fn parse(name: &str, notation: &str, commercial: bool) -> Result<Enzyme, EnzymeError> {
        let err = |reason: &str| EnzymeError {
            name: name.to_string(),
            reason: reason.to_string(),
        };
        let pair = |inner: &str| -> Result<(i32, i32), EnzymeError> {
            let (n, m) = inner.split_once('/').ok_or_else(|| err("expected (n/m)"))?;
            let n: i32 = n.trim().parse().map_err(|_| err("bad cut offset"))?;
            let m: i32 = m.trim().parse().map_err(|_| err("bad cut offset"))?;
            if n.abs() > MAX_OFFSET || m.abs() > MAX_OFFSET {
                return Err(err("cut offset too large"));
            }
            Ok((n, m))
        };
        let mut body = notation.trim();
        let mut upstream = None;
        if let Some(rest) = body.strip_prefix('(') {
            let (inner, rest) = rest.split_once(')').ok_or_else(|| err("unclosed '('"))?;
            upstream = Some(pair(inner)?);
            body = rest;
        }
        let offsets = match body.split_once('(') {
            Some((b, rest)) => {
                let inner = rest.strip_suffix(')').ok_or_else(|| err("unclosed '('"))?;
                body = b;
                Some(pair(inner)?)
            }
            None => None,
        };
        if upstream.is_some() && (offsets.is_none() || body.contains('^')) {
            return Err(err("an upstream cut needs a downstream (n/m) cut"));
        }
        let caret = body.find('^');
        let site: String = body
            .chars()
            .filter(|c| *c != '^')
            .collect::<String>()
            .to_ascii_uppercase();
        if site.is_empty() || site.len() > MAX_SITE {
            return Err(err("site must be 1 to 64 bases"));
        }
        if !site.bytes().all(is_iupac) {
            return Err(err("site has letters that aren't IUPAC nucleotide codes"));
        }
        if body.matches('^').count() > 1 {
            return Err(err("more than one '^'"));
        }
        let len = site.len() as i32;
        let (cut_top, cut_bottom) = match (caret, offsets) {
            (_, Some((n, m))) => (Some(len + n), Some(len + m)),
            (Some(c), None) => {
                // The bottom strand is cut symmetrically.
                let c = c as i32;
                (Some(c), Some(len - c))
            }
            (None, None) => (None, None),
        };
        Ok(Enzyme {
            name: name.to_string(),
            site,
            cut_top,
            cut_bottom,
            cut_upstream: upstream.map(|(n, m)| (-n, -m)),
            commercial,
        })
    }

    /// The site in standard notation (`G^AATTC`, `GGTCTC(1/5)`); [`Enzyme::parse`] reads it back.
    pub fn notation(&self) -> String {
        let len = self.site.len() as i32;
        let (Some(top), Some(bottom)) = (self.cut_top, self.cut_bottom) else {
            return self.site.clone();
        };
        match self.cut_upstream {
            Some((t, b)) => format!(
                "({}/{}){}({}/{})",
                -t,
                -b,
                self.site,
                top - len,
                bottom - len
            ),
            None if (0..=len).contains(&top) && bottom == len - top => {
                let (left, right) = self.site.split_at(top as usize);
                format!("{left}^{right}")
            }
            None => format!("{}({}/{})", self.site, top - len, bottom - len),
        }
    }

    /// Whether the cut positions are known.
    pub fn cuts(&self) -> bool {
        self.cut_top.is_some() && self.cut_bottom.is_some()
    }

    /// Bases in the site that aren't `N` (how specific the enzyme is).
    pub fn specific_length(&self) -> usize {
        self.site.bytes().filter(|&b| b != b'N').count()
    }

    /// Overhang: positive for a 5' overhang, negative for 3', zero for blunt.
    pub fn overhang(&self) -> Option<i32> {
        Some(self.cut_bottom? - self.cut_top?)
    }

    /// Whether the site reads the same on both strands.
    pub fn palindromic(&self) -> bool {
        reverse_complement(self.site.as_bytes()) == self.site.as_bytes()
    }
}

fn is_iupac(b: u8) -> bool {
    matches!(
        b,
        b'A' | b'C'
            | b'G'
            | b'T'
            | b'R'
            | b'Y'
            | b'S'
            | b'W'
            | b'K'
            | b'M'
            | b'B'
            | b'D'
            | b'H'
            | b'V'
            | b'N'
    )
}

/// The set of bases an IUPAC code stands for, as a bit mask (A=1, C=2, G=4, T=8).
fn mask(b: u8) -> u8 {
    match b.to_ascii_uppercase() {
        b'A' => 1,
        b'C' => 2,
        b'G' => 4,
        b'T' | b'U' => 8,
        b'R' => 1 | 4,
        b'Y' => 2 | 8,
        b'S' => 2 | 4,
        b'W' => 1 | 8,
        b'K' => 4 | 8,
        b'M' => 1 | 2,
        b'B' => 2 | 4 | 8,
        b'D' => 1 | 4 | 8,
        b'H' => 1 | 2 | 8,
        b'V' => 1 | 2 | 4,
        b'N' => 15,
        _ => 0,
    }
}

fn reverse_complement(s: &[u8]) -> Vec<u8> {
    s.iter()
        .rev()
        .map(|&b| match b {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            b'R' => b'Y',
            b'Y' => b'R',
            b'K' => b'M',
            b'M' => b'K',
            b'B' => b'V',
            b'V' => b'B',
            b'D' => b'H',
            b'H' => b'D',
            other => other,
        })
        .collect()
}

/// Whether the sequence base can be the pattern base. Ambiguous sequence bases (N, R, ...)
/// only match where every base they stand for would match, so sites are never invented.
fn matches(pattern: u8, base: u8) -> bool {
    let b = mask(base);
    b != 0 && b & mask(pattern) == b
}

/// Which strand a site was found on.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[serde(rename_all = "lowercase")]
pub enum SiteStrand {
    Forward,
    Reverse,
}

/// A recognition site found in a molecule.
#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
pub struct Site {
    /// Index of the enzyme in the list searched.
    pub enzyme: usize,
    /// The site's residues (wraps past the origin on circular molecules).
    pub range: Range,
    pub strand: SiteStrand,
    /// Top-strand cut as the number of residues before it (`0..=len`; `0..len` on circular
    /// molecules). `None` if the cut position is unknown or falls outside a linear molecule.
    pub cut: Option<usize>,
    /// Bottom-strand cut in the same coordinates.
    pub cut_bottom: Option<usize>,
    /// Top-strand position of the second cut of an enzyme that cuts on both sides.
    pub cut_upstream: Option<usize>,
}

/// Find every site of every enzyme in `m`, in enzyme order, then strand, then position.
/// Palindromic sites are reported once. On circular molecules sites may span the origin.
/// Ambiguous residues in the molecule never create a site.
pub fn find_sites(m: &Molecule, enzymes: &[Enzyme]) -> Vec<Site> {
    let seq = m.seq().as_bytes();
    let len = seq.len();
    let circular = m.topology() == Topology::Circular;
    let mut out = Vec::new();
    if len == 0 || !m.seq().alphabet().is_nucleotide() {
        return out;
    }
    for (index, e) in enzymes.iter().enumerate() {
        let fwd = e.site.as_bytes();
        let k = fwd.len();
        if k == 0 || k > len {
            continue;
        }
        let rev = reverse_complement(fwd);
        let mut strands = vec![(SiteStrand::Forward, fwd)];
        if rev != fwd {
            strands.push((SiteStrand::Reverse, rev.as_slice()));
        }
        let starts = if circular { len } else { len - k + 1 };
        let k = k as i64;
        for (strand, pattern) in strands {
            for start in 0..starts {
                let hit = pattern
                    .iter()
                    .enumerate()
                    .all(|(j, &p)| seq.get((start + j) % len).is_some_and(|&b| matches(p, b)));
                if !hit {
                    continue;
                }
                // On the reverse strand the enzyme's top strand is our bottom strand, read from
                // the other end of the site.
                let orient = |t: Option<i32>, b: Option<i32>| match strand {
                    SiteStrand::Forward => (t.map(i64::from), b.map(i64::from)),
                    SiteStrand::Reverse => {
                        (b.map(|c| k - i64::from(c)), t.map(|c| k - i64::from(c)))
                    }
                };
                let (top, bottom) = orient(e.cut_top, e.cut_bottom);
                let place = |offset: Option<i64>| -> Option<usize> {
                    let pos = start as i64 + offset?;
                    if circular {
                        usize::try_from(pos.rem_euclid(len as i64)).ok()
                    } else {
                        usize::try_from(pos).ok().filter(|&p| p <= len)
                    }
                };
                // A cut only happens if both strands are cut within the molecule.
                let (cut, cut_bottom) = match (place(top), place(bottom)) {
                    (Some(t), Some(b)) => (Some(t), Some(b)),
                    _ => (None, None),
                };
                let cut_upstream = e.cut_upstream.and_then(|(t, b)| {
                    let (t, b) = orient(Some(t), Some(b));
                    place(b)?;
                    place(t)
                });
                let end = start + k as usize;
                out.push(Site {
                    enzyme: index,
                    range: Range::new(start, if end > len { end - len } else { end }),
                    strand,
                    cut,
                    cut_bottom,
                    cut_upstream,
                });
            }
        }
    }
    out
}

/// A fragment of a digest.
#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
pub struct Fragment {
    /// From one top-strand cut to the next (wraps on circular molecules). A circular molecule
    /// cut once gives one fragment with `start == end` and `length` the whole molecule.
    pub range: Range,
    pub length: usize,
    /// Enzyme at each end (`None` for an end of a linear molecule).
    pub left: Option<String>,
    pub right: Option<String>,
}

/// Cut `m` at every site of `enzymes`; fragments sorted longest first. A molecule without
/// cuts is one fragment. Fragment lengths always add up to the molecule's length.
pub fn digest(m: &Molecule, enzymes: &[Enzyme]) -> Vec<Fragment> {
    let len = m.len();
    let circular = m.topology() == Topology::Circular;
    let mut cuts: Vec<(usize, String)> = find_sites(m, enzymes)
        .into_iter()
        .flat_map(|s| {
            let name = enzymes.get(s.enzyme).map(|e| e.name.clone());
            [s.cut, s.cut_upstream]
                .into_iter()
                .filter_map(move |c| Some((c?, name.clone()?)))
        })
        .filter(|(c, _)| circular || (*c > 0 && *c < len))
        .collect();
    cuts.sort();
    cuts.dedup_by(|a, b| a.0 == b.0);
    let mut out = Vec::new();
    if cuts.is_empty() {
        out.push(Fragment {
            range: Range::full(len),
            length: len,
            left: None,
            right: None,
        });
    } else if circular {
        let next = cuts.iter().cycle().skip(1);
        for ((a, an), (b, bn)) in cuts.iter().zip(next) {
            let length = match (b + len - a) % len {
                0 => len,
                l => l,
            };
            out.push(Fragment {
                range: Range::new(*a, *b),
                length,
                left: Some(an.clone()),
                right: Some(bn.clone()),
            });
        }
    } else {
        let mut prev: (usize, Option<String>) = (0, None);
        for (c, name) in cuts {
            out.push(Fragment {
                range: Range::new(prev.0, c),
                length: c - prev.0,
                left: prev.1.take(),
                right: Some(name.clone()),
            });
            prev = (c, Some(name));
        }
        out.push(Fragment {
            range: Range::new(prev.0, len),
            length: len - prev.0,
            left: prev.1,
            right: None,
        });
    }
    out.sort_by(|a, b| {
        b.length
            .cmp(&a.length)
            .then(a.range.start.cmp(&b.range.start))
    });
    out
}

/// Group sites by enzyme: `(enzyme index, its sites)` for every enzyme with at least one site.
pub fn sites_by_enzyme(sites: Vec<Site>) -> Vec<(usize, Vec<Site>)> {
    let mut out: Vec<(usize, Vec<Site>)> = Vec::new();
    for s in sites {
        match out.last_mut() {
            Some((i, v)) if *i == s.enzyme => v.push(s),
            _ => out.push((s.enzyme, vec![s])),
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use snappiesnap_doc::{Alphabet, Seq};

    use super::*;

    fn mol(s: &str, topology: Topology) -> Molecule {
        Molecule::new("m", Seq::new(Alphabet::Dna, s).unwrap(), topology).unwrap()
    }

    fn enz(name: &str, site: &str) -> Enzyme {
        Enzyme::parse(name, site, true).unwrap()
    }

    #[test]
    fn notation() {
        let e = enz("EcoRI", "G^AATTC");
        assert_eq!(
            (e.cut_top, e.cut_bottom, e.overhang()),
            (Some(1), Some(5), Some(4))
        );
        assert!(e.palindromic());
        let e = enz("BsaI", "ggtctc(1/5)");
        assert_eq!(
            (e.site.as_str(), e.cut_top, e.cut_bottom),
            ("GGTCTC", Some(7), Some(11))
        );
        assert!(!e.palindromic());
        assert_eq!(enz("PstI", "CTGCA^G").overhang(), Some(-4));
        assert_eq!(enz("SmaI", "CCC^GGG").overhang(), Some(0));
        assert_eq!(enz("BglI", "GCCNNNN^NGGC").specific_length(), 6);
        assert_eq!(enz("BaeI", "ACNNNNGTAYC(-10/-15)").cut_top, Some(1));
        assert!(!Enzyme::parse("X", "GGATCC", false).unwrap().cuts());
        let e = enz("Bsp24I", "(8/13)GACNNNNNNTGG(12/7)");
        assert_eq!(
            (e.cut_top, e.cut_bottom, e.cut_upstream),
            (Some(24), Some(19), Some((-8, -13)))
        );
        for bad in [
            "",
            "^",
            "?",
            "G^AA^TTC",
            "GAXTTC",
            "GGTCTC(1/5",
            "GGTCTC(a/5)",
            "GGTCTC(1)",
            "(1/5)GGTCTC",
            "(1/5)G^GTCTC(1/5)",
            "GGTCTC(99999/5)",
        ] {
            assert!(Enzyme::parse("Bad", bad, true).is_err(), "{bad:?}");
        }
    }

    #[test]
    fn notation_round_trips() {
        for n in [
            "G^AATTC",
            "^GATC",
            "CATG^",
            "GGTCTC(1/5)",
            "GAATGC(1/-1)",
            "GCANNNN^NTGC",
            "(8/13)GACNNNNNNTGG(12/7)",
            "GGATCC",
            "CAGGTACCCTTTAAACCTACTAACCC(-12/-16)",
        ] {
            let e = enz("E", n);
            assert_eq!(e.notation(), n);
            assert_eq!(enz("E", &e.notation()), e);
        }
        for e in builtin() {
            assert_eq!(Enzyme::parse(&e.name, &e.notation(), true).as_ref(), Ok(&e));
        }
    }

    #[test]
    fn sites_on_both_strands_and_across_the_origin() {
        let bsai = enz("BsaI", "GGTCTC(1/5)");
        let ecori = enz("EcoRI", "G^AATTC");
        // EcoRI across the origin (...G|AATTC...), BsaI forward at 10 and reverse (GAGACC) at 22.
        let s = "AATTCAAAAAGGTCTCAAAAAAGAGACCAAAG";
        let m = mol(s, Topology::Circular);
        let sites = find_sites(&m, &[bsai, ecori.clone()]);
        let eco: Vec<&Site> = sites.iter().filter(|s| s.enzyme == 1).collect();
        assert_eq!(eco.len(), 1);
        assert_eq!(eco[0].range, Range::new(31, 5));
        assert_eq!(eco[0].cut, Some(0), "cut right at the origin");
        let bsa: Vec<&Site> = sites.iter().filter(|s| s.enzyme == 0).collect();
        assert_eq!(bsa.len(), 2);
        assert_eq!(
            (bsa[0].strand, bsa[0].range),
            (SiteStrand::Forward, Range::new(10, 16))
        );
        assert_eq!((bsa[0].cut, bsa[0].cut_bottom), (Some(17), Some(21)));
        // GAGACC at 22..28: the top strand is cut 5 bases before it, the bottom strand 1.
        assert_eq!(
            (bsa[1].strand, bsa[1].range),
            (SiteStrand::Reverse, Range::new(22, 28))
        );
        assert_eq!((bsa[1].cut, bsa[1].cut_bottom), (Some(17), Some(21)));
        // On a linear molecule the origin-spanning site doesn't exist.
        assert!(find_sites(&mol(s, Topology::Linear), &[ecori]).is_empty());
    }

    #[test]
    fn cuts_outside_a_linear_molecule_are_dropped() {
        let bsai = enz("BsaI", "GGTCTC(1/5)");
        let sites = find_sites(&mol("AAGGTCTCAA", Topology::Linear), &[bsai]);
        assert_eq!(sites.len(), 1);
        assert_eq!((sites[0].cut, sites[0].cut_bottom), (None, None));
    }

    #[test]
    fn cuts_on_both_sides() {
        let e = enz("Bsp24I", "(8/13)GACNNNNNNTGG(12/7)");
        let s = format!("{}GACAAAAAATGG{}", "C".repeat(30), "C".repeat(30));
        let m = mol(&s, Topology::Linear);
        let sites = find_sites(&m, std::slice::from_ref(&e));
        assert_eq!(sites.len(), 1);
        assert_eq!(
            (sites[0].cut, sites[0].cut_upstream),
            (Some(30 + 24), Some(30 - 8))
        );
        let lens: Vec<usize> = digest(&m, &[e]).iter().map(|f| f.length).collect();
        assert_eq!(lens, vec![32, 22, 18]);
    }

    #[test]
    fn ambiguity_codes() {
        let acci = enz("AccI", "GT^MKAC");
        let m = mol("AAGTAGACAAGTCTACAAGTNNACAA", Topology::Linear);
        assert_eq!(
            find_sites(&m, &[acci]).len(),
            2,
            "GTAGAC and GTCTAC, not GTNNAC"
        );
    }

    #[test]
    fn digests() {
        let ecori = enz("EcoRI", "G^AATTC");
        let bamhi = enz("BamHI", "G^GATCC");
        let s = format!(
            "{}GAATTC{}GGATCC{}",
            "A".repeat(10),
            "A".repeat(20),
            "A".repeat(30)
        );
        let lens = |t, e: &[Enzyme]| -> Vec<usize> {
            digest(&mol(&s, t), e).iter().map(|f| f.length).collect()
        };
        let both = [ecori.clone(), bamhi];
        assert_eq!(lens(Topology::Circular, &both), vec![46, 26]);
        assert_eq!(lens(Topology::Linear, &both), vec![35, 26, 11]);
        assert_eq!(
            lens(Topology::Circular, std::slice::from_ref(&ecori)),
            vec![72]
        );
        let frags = digest(&mol(&s, Topology::Linear), &both);
        assert_eq!(frags[0].left.as_deref(), Some("BamHI"));
        assert_eq!(frags[0].right, None);
        assert_eq!(
            lens(Topology::Circular, &[enz("NotI", "GC^GGCCGC")]),
            vec![72]
        );
    }

    #[test]
    fn grouping() {
        let e = [enz("EcoRI", "G^AATTC"), enz("BamHI", "G^GATCC")];
        let m = mol("GAATTCGAATTCGGATCC", Topology::Linear);
        let groups = sites_by_enzyme(find_sites(&m, &e));
        let counts: Vec<(usize, usize)> = groups.iter().map(|(i, v)| (*i, v.len())).collect();
        assert_eq!(counts, vec![(0, 2), (1, 1)]);
    }

    #[test]
    fn hostile_molecules() {
        let all = builtin();
        for s in ["", "A", "N", "NNNNNNNNNNNN", "GAATTC", "AATTCG", "GCGGCCGC"] {
            for t in [Topology::Linear, Topology::Circular] {
                let m = mol(s, t);
                for site in find_sites(&m, &all) {
                    assert!(site.cut.is_none_or(|c| c <= m.len()));
                }
                let total: usize = digest(&m, &all).iter().map(|f| f.length).sum();
                assert_eq!(total, m.len(), "{s:?} {t:?}");
            }
        }
    }
}
