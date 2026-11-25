//! Assembly methods built on [`crate::cloning`]: homology-based assembly (Gibson/HiFi,
//! In-Fusion), Golden Gate, and site-directed mutagenesis.
//!
//! - Homology-based assembly (Gibson D. G. et al. (2009) Nat. Methods 6:343-345; In-Fusion uses
//!   15 bp end homology): consecutive pieces are joined where the end of one is identical to the
//!   start of the next over at least `min_overlap` bases; the shared stretch appears once. A piece
//!   is turned around when only its reverse complement fits.
//! - Golden Gate (Engler C., Kandzia R., Marillonnet S. (2008) PLoS ONE 3:e3647): every part is
//!   cut with a Type IIS enzyme, and fragments are chained by matching overhangs until the
//!   chain closes into a circle, each part used once.
//! - Site-directed mutagenesis by back-to-back primers: the whole plasmid is copied by PCR, the
//!   primers' changes included, and the linear product is ligated closed (blunt ends).

use serde::Serialize;
use snappiesnap_doc::{Molecule, Topology};

use crate::cloning::{DsFragment, Part, assemble, digest_fragments, ligate};
use crate::enzymes::Enzyme;
use crate::pcr::{DEFAULT_MAX_PRODUCT, pcr, product_molecule};
use crate::tm::{TmConditions, tm};

/// One junction of an assembly.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct Junction {
    /// The pieces it joins (index of the left piece; the right one follows, the last wraps to
    /// the first).
    pub after: usize,
    /// Overlap length (bp) for homology assembly, or overhang length for Golden Gate.
    pub length: usize,
    pub bases: String,
    /// Tm of the overlap (homology assembly only).
    pub tm: Option<f64>,
}

/// An assembled construct.
#[derive(Debug, Clone, PartialEq)]
pub struct Assembly {
    pub part: Part,
    pub junctions: Vec<Junction>,
    /// Pieces that were used turned around.
    pub reversed: Vec<usize>,
}

/// Longest suffix of `a` that equals a prefix of `b` (at least `min`, at most `max` bases).
fn overlap(a: &[u8], b: &[u8], min: usize, max: usize) -> Option<usize> {
    let limit = max.min(a.len()).min(b.len());
    (min..=limit).rev().find(|&k| {
        let tail = a.get(a.len() - k..);
        let head = b.get(..k);
        matches!((tail, head), (Some(t), Some(h)) if t.eq_ignore_ascii_case(h))
    })
}

/// Longest overlap considered.
const MAX_OVERLAP: usize = 200;

/// Join `pieces` (linear, in order) by end homology, closing the circle if `circular`.
pub fn homology_assembly(
    pieces: &[Part],
    min_overlap: usize,
    circular: bool,
) -> Result<Assembly, String> {
    if pieces.is_empty() {
        return Err("no pieces".to_string());
    }
    if pieces.len() > 50 {
        return Err("at most 50 pieces".to_string());
    }
    let min = min_overlap.max(5);
    let conditions = TmConditions::default();
    let mut oriented: Vec<Part> = Vec::new();
    let mut reversed = Vec::new();
    let mut junctions = Vec::new();
    for (i, p) in pieces.iter().enumerate() {
        let Some(prev) = oriented.last() else {
            oriented.push(p.clone());
            continue;
        };
        let fwd = overlap(&prev.seq, &p.seq, min, MAX_OVERLAP);
        let rev_part = p.reversed();
        let rev = overlap(&prev.seq, &rev_part.seq, min, MAX_OVERLAP);
        let (part, k) = match (fwd, rev) {
            (Some(k), _) => (p.clone(), k),
            (None, Some(k)) => {
                reversed.push(i);
                (rev_part, k)
            }
            (None, None) => {
                return Err(format!(
                    "pieces {} and {} don't share at least {min} bp at their ends",
                    i,
                    i + 1
                ));
            }
        };
        let bases = part.seq.get(..k).unwrap_or_default().to_vec();
        junctions.push(Junction {
            after: i - 1,
            length: k,
            bases: String::from_utf8_lossy(&bases).into_owned(),
            tm: tm(&bases, &conditions),
        });
        oriented.push(part);
    }
    // Join, dropping each shared stretch from the end of the left piece.
    let mut joined = oriented.first().cloned().unwrap_or(Part {
        seq: Vec::new(),
        features: Vec::new(),
    });
    for (j, p) in junctions.iter().zip(oriented.iter().skip(1)) {
        joined = crate::cloning::join_trimmed(joined, j.length, p);
    }
    if circular {
        let first = oriented.first().map(|p| p.seq.clone()).unwrap_or_default();
        let k = overlap(&joined.seq, &first, min, MAX_OVERLAP).ok_or_else(|| {
            format!("the last piece and the first don't share at least {min} bp, so the circle can't close")
        })?;
        if k >= joined.seq.len() {
            return Err("the overlap is as long as the construct".to_string());
        }
        let bases = first.get(..k).unwrap_or_default().to_vec();
        junctions.push(Junction {
            after: pieces.len() - 1,
            length: k,
            bases: String::from_utf8_lossy(&bases).into_owned(),
            tm: tm(&bases, &conditions),
        });
        joined = crate::cloning::trim_end(joined, k);
    }
    Ok(Assembly {
        part: joined,
        junctions,
        reversed,
    })
}

/// Golden Gate: cut every source with `enzyme` and chain fragments by their overhangs, starting
/// from the largest fragment of the first source (the destination vector), each source used
/// once, until the chain closes into a circle.
pub fn golden_gate(sources: &[Molecule], enzyme: &Enzyme) -> Result<Assembly, String> {
    if sources.len() < 2 {
        return Err("give a vector and at least one part".to_string());
    }
    if sources.len() > 30 {
        return Err("at most 30 parts".to_string());
    }
    let cut: Vec<Vec<DsFragment>> = sources
        .iter()
        .map(|m| digest_fragments(m, std::slice::from_ref(enzyme)))
        .collect();
    for (i, f) in cut.iter().enumerate() {
        if f.is_empty() {
            return Err(format!(
                "{} isn't cut by {}",
                sources.get(i).map_or("?", |m| m.name()),
                enzyme.name
            ));
        }
    }
    let start = cut
        .first()
        .and_then(|f| f.iter().max_by_key(|x| x.len()))
        .cloned()
        .ok_or("the vector isn't cut")?;
    // Depth-first search over parts not yet used.
    fn extend(
        chain: &mut Vec<(usize, DsFragment)>,
        used: &mut Vec<bool>,
        cut: &[Vec<DsFragment>],
        steps: &mut usize,
    ) -> bool {
        *steps += 1;
        if *steps > 100_000 {
            return false;
        }
        let (Some(first_left), Some(last_right)) = (
            chain.first().map(|f| f.1.left.clone()),
            chain.last().map(|f| f.1.right.clone()),
        ) else {
            return false;
        };
        if chain.len() > 1 && last_right.ligates(&first_left) && used.iter().all(|&u| u) {
            return true;
        }
        for (src, frags) in cut.iter().enumerate() {
            if used.get(src).copied().unwrap_or(true) {
                continue;
            }
            for f in frags {
                for cand in [f.clone(), f.reversed()] {
                    if last_right.ligates(&cand.left) {
                        if let Some(u) = used.get_mut(src) {
                            *u = true;
                        }
                        chain.push((src, cand));
                        if extend(chain, used, cut, steps) {
                            return true;
                        }
                        chain.pop();
                        if let Some(u) = used.get_mut(src) {
                            *u = false;
                        }
                    }
                }
            }
        }
        false
    }
    let mut used = vec![false; sources.len()];
    if let Some(u) = used.first_mut() {
        *u = true;
    }
    let mut chain = vec![(0, start)];
    let mut steps = 0;
    if !extend(&mut chain, &mut used, &cut, &mut steps) {
        return Err("the parts' overhangs don't chain into a circle with the vector".to_string());
    }
    let fragments: Vec<DsFragment> = chain.iter().map(|(_, f)| f.clone()).collect();
    let part = ligate(&fragments, true)?;
    let junctions = fragments
        .iter()
        .enumerate()
        .map(|(i, f)| {
            let bases = match &f.right {
                crate::cloning::End::Overhang5(s) | crate::cloning::End::Overhang3(s) => s.clone(),
                crate::cloning::End::Blunt => String::new(),
            };
            Junction {
                after: i,
                length: bases.len(),
                bases,
                tm: None,
            }
        })
        .collect();
    Ok(Assembly {
        part,
        junctions,
        reversed: Vec::new(),
    })
}

/// Site-directed mutagenesis with back-to-back primers: the whole-plasmid PCR product, ligated
/// closed, with the template's origin kept where it still exists.
pub fn mutagenesis(
    template: &Molecule,
    forward: &[u8],
    reverse: &[u8],
    names: [&str; 2],
    name: &str,
) -> Result<Molecule, String> {
    if template.topology() != Topology::Circular {
        return Err("mutagenesis by back-to-back primers needs a circular template".to_string());
    }
    let products = pcr(template, [forward, reverse], None, DEFAULT_MAX_PRODUCT);
    // The product that copies (nearly) the whole plasmid.
    let len = template.len();
    let p = products
        .iter()
        .filter(|p| p.span + 200 >= len)
        .max_by_key(|p| p.span)
        .ok_or("the primers don't copy the whole plasmid back to back")?;
    let mut m = product_molecule(template, p, names, name).map_err(|e| e.to_string())?;
    m.set_topology(Topology::Circular)
        .map_err(|e| e.to_string())?;
    // Put the template's origin back at 0: it sits after the forward primer's tail.
    let origin_in_template = (len - p.template.start) % len;
    let origin = p.forward_tail + origin_in_template;
    if origin < m.len() {
        m.rotate_origin(origin).map_err(|e| e.to_string())?;
    }
    Ok(m)
}

/// Build the molecule for an assembly.
pub fn to_molecule(a: &Assembly, name: &str, circular: bool) -> Result<Molecule, String> {
    let topology = if circular {
        Topology::Circular
    } else {
        Topology::Linear
    };
    assemble(std::slice::from_ref(&a.part), topology, name).map_err(|e| e.to_string())
}

#[cfg(test)]
mod tests {
    use snappiesnap_doc::{Alphabet, NewFeature, Range, Segment, Seq, Strand};

    use super::*;
    use crate::cloning::{extract, revcomp};

    fn rnd(seed: u32, n: usize) -> String {
        let mut x = seed;
        (0..n)
            .map(|_| {
                x = x.wrapping_mul(1_103_515_245).wrapping_add(12345);
                b"ACGT"[((x >> 16) % 4) as usize] as char
            })
            .collect()
    }

    fn part(s: &str) -> Part {
        Part {
            seq: s.as_bytes().to_vec(),
            features: Vec::new(),
        }
    }

    #[test]
    fn gibson_three_pieces_into_a_circle() {
        let (a, b, c) = (rnd(1, 300), rnd(2, 200), rnd(3, 250));
        let ov = |s: &str, from_end: bool| -> String {
            if from_end {
                s[s.len() - 20..].to_string()
            } else {
                s[..20].to_string()
            }
        };
        // Each piece carries 20 bp of its neighbour at its start.
        let p1 = format!("{}{}", ov(&c, true), a);
        let p2 = format!("{}{}", ov(&a, true), b);
        let p3 = format!("{}{}", ov(&b, true), c);
        let asm = homology_assembly(&[part(&p1), part(&p2), part(&p3)], 15, true).unwrap();
        assert_eq!(asm.part.seq.len(), 750);
        assert_eq!(asm.junctions.len(), 3);
        assert!(
            asm.junctions
                .iter()
                .all(|j| j.length == 20 && j.tm.is_some())
        );
        let s = String::from_utf8(asm.part.seq).unwrap();
        let doubled = format!("{a}{b}{c}{a}{b}{c}");
        assert!(doubled.contains(&s));
        // A piece given the wrong way round is turned.
        let rev2 = String::from_utf8(revcomp(p2.as_bytes())).unwrap();
        let asm = homology_assembly(&[part(&p1), part(&rev2), part(&p3)], 15, true).unwrap();
        assert_eq!(asm.reversed, vec![1]);
        assert_eq!(asm.part.seq.len(), 750);
        // Too little homology.
        assert!(homology_assembly(&[part(&p1), part(&p2), part(&p3)], 25, true).is_err());
    }

    #[test]
    fn features_survive_assembly() {
        let (a, b) = (rnd(4, 100), rnd(5, 100));
        let mut m = Molecule::new(
            "a",
            Seq::new(Alphabet::Dna, format!("{}{a}", &b[80..])).unwrap(),
            Topology::Linear,
        )
        .unwrap();
        m.add_feature(NewFeature {
            color: None,
            name: "gene".into(),
            kind: "CDS".into(),
            strand: Strand::Forward,
            segments: vec![Segment {
                range: Range::new(40, 90),
            }],
            qualifiers: vec![],
        })
        .unwrap();
        let p1 = extract(&m, 0, m.len());
        let p2 = part(&format!("{}{b}", &a[80..]));
        let asm = homology_assembly(&[p1, p2], 15, true).unwrap();
        assert_eq!(asm.part.seq.len(), 200);
        assert!(asm.part.features.iter().any(|f| f.name == "gene"));
    }

    #[test]
    fn golden_gate_three_parts() {
        let bsai = Enzyme::parse("BsaI", "GGTCTC(1/5)", true).unwrap();
        // Vector: backbone, GGAG, then a reverse BsaI site (GAGACC) pointing into the dropout,
        // the dropout, a forward site (GGTCTC) pointing into it, CGCT, backbone. Cutting leaves
        // the backbone with a CGCT end on the left and a GGAG end on the right.
        let vec = format!(
            "{}GGAGAGAGACC{}GGTCTCACGCT{}",
            rnd(6, 60),
            rnd(7, 30),
            rnd(8, 60)
        );
        let p1 = format!("AAGGTCTCAGGAG{}AATGTGAGACCAA", rnd(9, 40));
        let p2 = format!("AAGGTCTCAAATG{}CGCTTGAGACCAA", rnd(10, 50));
        let ms: Vec<Molecule> = [vec.clone(), p1, p2]
            .iter()
            .enumerate()
            .map(|(i, s)| {
                Molecule::new(
                    &format!("m{i}"),
                    Seq::new(Alphabet::Dna, s.as_str()).unwrap(),
                    if i == 0 {
                        Topology::Circular
                    } else {
                        Topology::Linear
                    },
                )
                .unwrap()
            })
            .collect();
        let asm = golden_gate(&ms, &bsai).unwrap();
        let s = String::from_utf8(asm.part.seq.clone()).unwrap();
        assert!(!s.contains("GGTCTC") && !s.contains("GAGACC"), "{s}");
        let overhangs: Vec<&str> = asm.junctions.iter().map(|j| j.bases.as_str()).collect();
        assert_eq!(overhangs.len(), 3);
        for o in ["GGAG", "AATG", "CGCT"] {
            assert!(overhangs.contains(&o), "{overhangs:?}");
        }
        // A part that doesn't fit breaks the chain.
        let bad = Molecule::new(
            "bad",
            Seq::new(
                Alphabet::Dna,
                format!("AAGGTCTCATTTT{}CCCCTGAGACCAA", rnd(11, 30)),
            )
            .unwrap(),
            Topology::Linear,
        )
        .unwrap();
        assert!(golden_gate(&[ms[0].clone(), bad], &bsai).is_err());
    }

    #[test]
    fn mutagenesis_changes_a_base_and_keeps_the_origin() {
        let t = rnd(12, 400);
        let m = Molecule::new(
            "p",
            Seq::new(Alphabet::Dna, t.as_str()).unwrap(),
            Topology::Circular,
        )
        .unwrap();
        // Change base 210 (in the forward primer's 5' part) from whatever it is to something else.
        let mut fwd: Vec<u8> = t[205..230].bytes().collect();
        let old = fwd[5];
        fwd[5] = if old == b'A' { b'C' } else { b'A' };
        let rev = revcomp(&t.as_bytes()[180..205]);
        let mutant = mutagenesis(&m, &fwd, &rev, ["F", "R"], "mut").unwrap();
        let s = mutant.seq().as_bytes();
        assert_eq!(s.len(), 400);
        assert_eq!(s[210], fwd[5]);
        let mut expected = t.as_bytes().to_vec();
        expected[210] = fwd[5];
        assert_eq!(s, &expected[..]);
        assert_eq!(mutant.topology(), Topology::Circular);
        // A deletion: the forward primer starts 10 bases further on.
        let fwd_del = &t.as_bytes()[215..240];
        let del = mutagenesis(&m, fwd_del, &rev, ["F", "R"], "del").unwrap();
        assert_eq!(del.len(), 390);
    }

    #[test]
    fn hostile() {
        assert!(homology_assembly(&[], 15, true).is_err());
        let _ = homology_assembly(&[part(""), part("")], 0, true);
        let _ = homology_assembly(
            &[
                part("AAAAAAAAAAAAAAAAAAAAAA"),
                part("AAAAAAAAAAAAAAAAAAAAAA"),
            ],
            15,
            true,
        );
        let e = Enzyme::parse("BsaI", "GGTCTC(1/5)", true).unwrap();
        let m = Molecule::new(
            "x",
            Seq::new(Alphabet::Dna, "ACGT").unwrap(),
            Topology::Circular,
        )
        .unwrap();
        assert!(golden_gate(&[m.clone(), m.clone()], &e).is_err());
        assert!(mutagenesis(&m, b"ACGT", b"ACGT", ["a", "b"], "x").is_err());
    }
}
