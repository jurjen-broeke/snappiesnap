//! PCR simulation: the products two primers make from a template.
//!
//! A product forms between a primer bound on the top strand (extending rightwards) and a primer
//! bound on the bottom strand (extending leftwards) further along; either primer may take either
//! role, and a single primer can pair with itself. Primer binding follows [`crate::primers`]
//! (primary sites by default). After the first cycles the product's ends are the primers
//! themselves, so 5' tails are added and any primer/template mismatches become part of the
//! product: product = forward primer + template between the annealed stretches + reverse
//! complement of the reverse primer.

use serde::Serialize;
use snappiesnap_doc::{
    Alphabet, DocError, Molecule, NewFeature, Qualifier, Range, Segment, Seq, Strand, Topology,
};

use crate::primers::{
    BindingPrefs, BindingRules, BindingSite, PrimerStrand, binding_sites, sites_with_prefs,
};

/// Longest product reported by default.
pub const DEFAULT_MAX_PRODUCT: usize = 50_000;
/// Most products reported (primers that bind everywhere would otherwise make millions).
pub const MAX_PRODUCTS: usize = 200;

/// A PCR product.
#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
pub struct Product {
    pub length: usize,
    /// The product, 5'->3' on the template's top strand.
    #[serde(skip)]
    pub seq: Vec<u8>,
    /// Which primer (0 or 1) primes on the top strand, and which on the bottom.
    pub forward: usize,
    pub reverse: usize,
    /// Where the primers anneal on the template.
    pub forward_site: Range,
    pub reverse_site: Range,
    /// The template stretch copied, from the forward site's start to the reverse site's end
    /// (wraps on circular templates; `start == end` with `span == len` for the whole circle).
    pub template: Range,
    /// Length of that stretch.
    pub span: usize,
    /// 5' tails (primer bases outside the annealed stretch) at the forward and reverse ends.
    pub forward_tail: usize,
    pub reverse_tail: usize,
}

fn complement(b: u8) -> u8 {
    match b.to_ascii_uppercase() {
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
        o => o,
    }
}

fn revcomp(s: &[u8]) -> Vec<u8> {
    s.iter().rev().map(|&b| complement(b)).collect()
}

/// The products `primers` make on `m` (at most [`MAX_PRODUCTS`]), shortest first. Only primary
/// binding sites are used, unless `secondary` gives the rules for weaker ones as well
/// ([`sites_with_prefs`], Tm at the default conditions). Products longer than `max_len` are
/// left out.
pub fn pcr(
    m: &Molecule,
    primers: [&[u8]; 2],
    secondary: Option<&BindingPrefs>,
    max_len: usize,
) -> Vec<Product> {
    let len = m.len();
    let circular = m.topology() == Topology::Circular;
    let mut out = Vec::new();
    if len == 0 {
        return out;
    }
    let sites: Vec<BindingSite> = match secondary {
        Some(prefs) => sites_with_prefs(m, &primers, prefs, &crate::tm::TmConditions::default()),
        None => binding_sites(m, &primers, BindingRules::default())
            .into_iter()
            .filter(|s| primers.get(s.primer).is_some_and(|p| s.is_primary(p.len())))
            .collect(),
    };
    let seq = m.seq().as_bytes();
    let fwd_sites = sites.iter().filter(|s| s.strand == PrimerStrand::Forward);
    'outer: for f in fwd_sites {
        for r in sites.iter().filter(|s| s.strand == PrimerStrand::Reverse) {
            let (Some(fp), Some(rp)) = (primers.get(f.primer), primers.get(r.primer)) else {
                continue;
            };
            // Offset of the reverse site from the forward site's start, along the top strand.
            let fs = f.range.start;
            let rs = r.range.start;
            let d = if circular {
                (rs + len - fs) % len
            } else if rs >= fs {
                rs - fs
            } else {
                continue;
            };
            // The annealed stretches mustn't overlap, and on a circle the product can't lap the
            // template (back-to-back primers copy exactly all of it).
            if d < f.annealed || (circular && d + r.annealed > len) {
                continue;
            }
            if out.len() >= MAX_PRODUCTS {
                break 'outer;
            }
            let between = d - f.annealed;
            let length = fp.len() + between + rp.len();
            if length > max_len {
                continue;
            }
            let mut product = Vec::with_capacity(length);
            product.extend(fp.iter().map(u8::to_ascii_uppercase));
            for k in 0..between {
                let p = (fs + f.annealed + k) % len;
                product.push(seq.get(p).copied().unwrap_or(b'N').to_ascii_uppercase());
            }
            product.extend(revcomp(rp));
            let end = (fs + d + r.annealed) % len;
            out.push(Product {
                length,
                seq: product,
                forward: f.primer,
                reverse: r.primer,
                forward_site: f.range,
                reverse_site: r.range,
                template: Range::new(fs, if !circular && end == 0 { len } else { end }),
                span: d + r.annealed,
                forward_tail: fp.len().saturating_sub(f.annealed),
                reverse_tail: rp.len().saturating_sub(r.annealed),
            });
        }
    }
    out.sort_by_key(|p| (p.length, p.forward_site.start));
    out.dedup_by(|a, b| a.seq == b.seq);
    out
}

/// The product as a molecule: linear, with the template's features inside the copied stretch
/// (clipped and noted as partial where they extend beyond it), and the two primers as
/// `primer_bind` features.
pub fn product_molecule(
    m: &Molecule,
    p: &Product,
    primer_names: [&str; 2],
    name: &str,
) -> Result<Molecule, DocError> {
    let mut out = Molecule::new(
        name,
        Seq::new(Alphabet::Dna, p.seq.clone())?,
        Topology::Linear,
    )?;
    let len = m.len().max(1);
    // Template offset 0 (the forward site's start) sits after the forward primer's tail.
    let offset = p.forward_tail;
    // Features overlapping the copied stretch, clipped to it (a clipped feature is noted as
    // partial).
    let span = p.span as i64;
    let l = len as i64;
    for f in m.features() {
        let mut segments = Vec::new();
        let mut clipped = false;
        for seg in &f.segments {
            let s_len = seg.range.len_on(m.len()) as i64;
            // Signed offset of the segment from the copied stretch's start.
            let mut o = (seg.range.start as i64 - p.template.start as i64).rem_euclid(l);
            if o >= span && o + s_len > l {
                o -= l;
            }
            let (a, b) = (o.max(0), (o + s_len).min(span));
            if s_len == 0 || a >= b {
                clipped = true;
                continue;
            }
            clipped |= a != o || b != o + s_len;
            let start = usize::try_from(a).unwrap_or(0) + offset;
            let end = usize::try_from(b).unwrap_or(0) + offset;
            segments.push(Segment {
                range: Range::new(start, end),
            });
        }
        if segments.is_empty() {
            continue;
        }
        let mut qualifiers = f.qualifiers.clone();
        if clipped {
            qualifiers.push(Qualifier {
                name: "note".to_string(),
                values: vec![
                    "partial: only part of this feature is in the PCR product".to_string(),
                ],
            });
        }
        out.add_feature(NewFeature {
            color: f.color.clone(),
            name: f.name.clone(),
            kind: f.kind.clone(),
            strand: f.strand,
            segments,
            qualifiers,
        })?;
    }
    let total = p.seq.len();
    let fwd_len = p.forward_tail + p.forward_site.len_on(m.len());
    let rev_len = p.reverse_tail + p.reverse_site.len_on(m.len());
    for (which, range, strand) in [
        (
            p.forward,
            Range::new(0, fwd_len.min(total)),
            Strand::Forward,
        ),
        (
            p.reverse,
            Range::new(total.saturating_sub(rev_len), total),
            Strand::Reverse,
        ),
    ] {
        out.add_feature(NewFeature {
            color: None,
            name: primer_names
                .get(which)
                .copied()
                .unwrap_or("primer")
                .to_string(),
            kind: "primer_bind".to_string(),
            strand,
            segments: vec![Segment { range }],
            qualifiers: vec![Qualifier {
                name: "note".to_string(),
                values: vec!["PCR primer".to_string()],
            }],
        })?;
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use snappiesnap_doc::Segment;

    use super::*;

    /// A made-up 400 bp sequence.
    fn template() -> String {
        let mut x: u32 = 4242;
        (0..400)
            .map(|_| {
                x = x.wrapping_mul(1_103_515_245).wrapping_add(12345);
                b"ACGT"[((x >> 16) % 4) as usize] as char
            })
            .collect()
    }

    fn rc(s: &str) -> String {
        String::from_utf8(revcomp(s.as_bytes())).unwrap()
    }

    #[test]
    fn product_with_tails() {
        let t = template();
        let mut m = Molecule::new(
            "t",
            Seq::new(Alphabet::Dna, t.as_str()).unwrap(),
            Topology::Linear,
        )
        .unwrap();
        m.add_feature(NewFeature {
            color: None,
            name: "inside".into(),
            kind: "misc_feature".into(),
            strand: Strand::Forward,
            segments: vec![Segment {
                range: Range::new(150, 200),
            }],
            qualifiers: vec![],
        })
        .unwrap();
        m.add_feature(NewFeature {
            color: None,
            name: "outside".into(),
            kind: "misc_feature".into(),
            strand: Strand::None,
            segments: vec![Segment {
                range: Range::new(10, 60),
            }],
            qualifiers: vec![],
        })
        .unwrap();
        let fwd = format!("gaattc{}", &t[100..120]);
        let rev = format!("ggatcc{}", rc(&t[280..300]));
        let products = pcr(
            &m,
            [fwd.as_bytes(), rev.as_bytes()],
            None,
            DEFAULT_MAX_PRODUCT,
        );
        assert_eq!(products.len(), 1, "{products:?}");
        let p = &products[0];
        assert_eq!(p.length, 6 + 200 + 6);
        assert_eq!(
            String::from_utf8(p.seq.clone()).unwrap(),
            format!("GAATTC{}{}", &t[100..300], rc("GGATCC"))
        );
        // (The tail's last base may happen to pair too, so check the bookkeeping, not 6/6.)
        assert_eq!(p.forward_tail + p.span + p.reverse_tail, p.length);
        assert_eq!((p.forward, p.reverse), (0, 1));
        assert_eq!(p.template.end, 300);
        let pm = product_molecule(&m, p, ["F", "R"], "product").unwrap();
        let names: Vec<(&str, Range)> = pm
            .features()
            .iter()
            .map(|f| (f.name.as_str(), f.segments[0].range))
            .collect();
        assert_eq!(
            names,
            vec![
                ("inside", Range::new(56, 106)),
                ("F", Range::new(0, 26)),
                ("R", Range::new(186, 212))
            ]
        );
        // A feature hanging over the forward end is clipped and noted.
        let mut m2 = m.clone();
        m2.add_feature(NewFeature {
            color: None,
            name: "over".into(),
            kind: "misc_feature".into(),
            strand: Strand::None,
            segments: vec![Segment {
                range: Range::new(80, 120),
            }],
            qualifiers: vec![],
        })
        .unwrap();
        let pm = product_molecule(&m2, p, ["F", "R"], "product").unwrap();
        let over = pm.features().iter().find(|f| f.name == "over").unwrap();
        assert_eq!(over.segments[0].range.start, p.forward_tail);
        assert!(
            over.qualifiers
                .iter()
                .any(|q| q.values.iter().any(|v| v.starts_with("partial")))
        );
    }

    #[test]
    fn across_the_origin_and_either_orientation() {
        let t = template();
        let m = Molecule::new(
            "t",
            Seq::new(Alphabet::Dna, t.as_str()).unwrap(),
            Topology::Circular,
        )
        .unwrap();
        // "Reverse" given first: roles are found from binding.
        let a = rc(&t[30..50]);
        let b = t[350..370].to_string();
        let products = pcr(&m, [a.as_bytes(), b.as_bytes()], None, DEFAULT_MAX_PRODUCT);
        assert_eq!(products.len(), 1);
        let p = &products[0];
        assert_eq!((p.forward, p.reverse), (1, 0));
        assert_eq!(p.length, 400 - 350 + 50);
        assert_eq!(p.template, Range::new(350, 50));
        let expected = format!("{}{}", &t[350..], &t[..50]);
        assert_eq!(String::from_utf8(p.seq.clone()).unwrap(), expected);
        // On a linear template there's no product across the end.
        let lin = Molecule::new(
            "t",
            Seq::new(Alphabet::Dna, t.as_str()).unwrap(),
            Topology::Linear,
        )
        .unwrap();
        assert!(
            pcr(
                &lin,
                [a.as_bytes(), b.as_bytes()],
                None,
                DEFAULT_MAX_PRODUCT
            )
            .is_empty()
        );
        // Too long for the limit.
        assert!(pcr(&m, [a.as_bytes(), b.as_bytes()], None, 50).is_empty());
    }

    #[test]
    fn back_to_back_primers_copy_the_whole_plasmid() {
        let t = template();
        let m = Molecule::new(
            "t",
            Seq::new(Alphabet::Dna, t.as_str()).unwrap(),
            Topology::Circular,
        )
        .unwrap();
        // Forward starts at 200; reverse ends at 200 (inverse PCR around the plasmid).
        let fwd = &t[200..220];
        let rev = rc(&t[180..200]);
        let p = pcr(
            &m,
            [fwd.as_bytes(), rev.as_bytes()],
            None,
            DEFAULT_MAX_PRODUCT,
        );
        assert_eq!(p.len(), 1);
        assert_eq!((p[0].length, p[0].span), (400, 400));
        assert_eq!(
            String::from_utf8(p[0].seq.clone()).unwrap(),
            format!("{}{}", &t[200..], &t[..200])
        );
    }

    #[test]
    fn mutagenic_primer_changes_the_product() {
        let t = template();
        let m = Molecule::new(
            "t",
            Seq::new(Alphabet::Dna, t.as_str()).unwrap(),
            Topology::Linear,
        )
        .unwrap();
        // A mismatch 8 bases from the 3' end: the site is secondary, so ask for those too.
        let mut fwd: Vec<u8> = t[100..125].bytes().collect();
        fwd[16] = if fwd[16] == b'A' { b'C' } else { b'A' };
        let rev = rc(&t[280..300]);
        assert!(pcr(&m, [&fwd, rev.as_bytes()], None, DEFAULT_MAX_PRODUCT).is_empty());
        let d = BindingPrefs::default();
        let weak = BindingPrefs {
            long: crate::primers::ClassRules {
                min_tm_on: false,
                ..d.long
            },
            ..d
        };
        let p = pcr(&m, [&fwd, rev.as_bytes()], Some(&weak), DEFAULT_MAX_PRODUCT);
        assert_eq!(p.len(), 1);
        assert_eq!(p[0].seq.get(16), fwd.get(16));
        assert_eq!(p[0].length, 200);
    }

    #[test]
    fn hostile() {
        let m = Molecule::new(
            "t",
            Seq::new(Alphabet::Dna, "").unwrap(),
            Topology::Circular,
        )
        .unwrap();
        assert!(pcr(&m, [b"ACGT", b""], Some(&BindingPrefs::default()), 10).is_empty());
        let t = template();
        let m = Molecule::new(
            "t",
            Seq::new(Alphabet::Dna, t.as_str()).unwrap(),
            Topology::Circular,
        )
        .unwrap();
        // The same primer twice, a palindrome, nonsense.
        for p in [
            &b"NNNNNNNNNNNNNNNNNNNN"[..],
            b"",
            b"ACGTACGTACGTACGTACGT",
            &t.as_bytes()[0..20],
        ] {
            for q in pcr(&m, [p, p], Some(&BindingPrefs::default()), usize::MAX) {
                assert!(q.length <= 3 * t.len());
                let _ = product_molecule(&m, &q, ["a", "b"], "p");
            }
        }
    }
}
