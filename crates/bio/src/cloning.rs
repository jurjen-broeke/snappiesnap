//! Cloning simulation: cutting molecules into double-stranded fragments with their ends,
//! ligating compatible ends, and Gateway recombination.
//!
//! A fragment is kept as the stretch of sequence that either strand covers, with the kind of end
//! at each side: blunt, or a 5' or 3' overhang whose bases are written as they read on the top
//! strand. Two ends ligate when they are the same kind with the same bases (DNA ligase joins
//! annealed complementary overhangs, or two blunt ends). Features travel with the sequence;
//! a feature cut by a fragment end is clipped and noted as partial.
//!
//! Gateway recombination (Hartley J. L., Temple G. F., Brasch M. A. (2000) Genome Res.
//! 10:1788-1795) exchanges the stretch between two att sites: the reaction crosses over inside
//! the core shared by the recombining sites, `TTTGTACAAAAAAGC` for att1 and `CTTTCTTGTACAAAG`
//! for att2 (found inside attB1 `ACAAGTTTGTACAAAAAAGCAGGCT` and attB2
//! `ACCCAGCTTTCTTGTACAAAGTGGT`). LR (attL × attR → attB) and BP (attB × attP → attL) both put
//! the donor's stretch between its cores in place of the acceptor's.
//!
//! Because all four kinds of site share those cores, which kind a site is can only be read
//! from the arms around it, and a pair that cannot recombine is refused rather than swapped:
//! see [`AttKind`] and [`att_kinds`].

use serde::Serialize;
use snappiesnap_doc::{
    Alphabet, DocError, Molecule, NewFeature, Qualifier, Range, Segment, Seq, Strand, Topology,
};

use crate::enzymes::{Enzyme, find_sites};

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

pub(crate) fn revcomp(s: &[u8]) -> Vec<u8> {
    s.iter().rev().map(|&b| complement(b)).collect()
}

/// A stretch of sequence with the features on it (in its own coordinates).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Part {
    pub seq: Vec<u8>,
    pub features: Vec<NewFeature>,
}

fn partial_note() -> Qualifier {
    Qualifier {
        name: "note".to_string(),
        values: vec!["partial: only part of this feature is in the construct".to_string()],
    }
}

/// The stretch `start..start + length` of `m` (wrapping on circular molecules), with the
/// features on it clipped to it.
pub fn extract(m: &Molecule, start: usize, length: usize) -> Part {
    let len = m.len();
    if len == 0 {
        return Part {
            seq: Vec::new(),
            features: Vec::new(),
        };
    }
    let circular = m.topology() == Topology::Circular;
    let bases = m.seq().as_bytes();
    let length = if circular {
        length.min(len)
    } else {
        length.min(len.saturating_sub(start))
    };
    let seq: Vec<u8> = (0..length)
        .filter_map(|k| bases.get((start + k) % len).copied())
        .collect();
    let (l, span, s0) = (len as i64, length as i64, start as i64);
    let mut features = Vec::new();
    for f in m.features() {
        let mut segments = Vec::new();
        let mut clipped = false;
        for seg in &f.segments {
            let s_len = seg.range.len_on(len) as i64;
            let mut o = (seg.range.start as i64 - s0).rem_euclid(l);
            if circular && o >= span && o + s_len > l {
                o -= l;
            }
            if !circular {
                o = seg.range.start as i64 - s0;
            }
            let (a, b) = (o.max(0), (o + s_len).min(span));
            if s_len == 0 || a >= b {
                clipped = true;
                continue;
            }
            clipped |= a != o || b != o + s_len;
            segments.push(Segment {
                range: Range::new(a as usize, b as usize),
            });
        }
        if segments.is_empty() {
            continue;
        }
        let mut qualifiers = f.qualifiers.clone();
        if clipped {
            qualifiers.push(partial_note());
        }
        features.push(NewFeature {
            color: f.color.clone(),
            name: f.name.clone(),
            kind: f.kind.clone(),
            strand: f.strand,
            segments,
            qualifiers,
        });
    }
    Part { seq, features }
}

impl Part {
    /// The part turned around (reverse complement), its features with it.
    pub fn reversed(&self) -> Part {
        let n = self.seq.len();
        let features = self
            .features
            .iter()
            .map(|f| NewFeature {
                segments: f
                    .segments
                    .iter()
                    .rev()
                    .map(|s| Segment {
                        range: Range::new(n - s.range.end, n - s.range.start),
                    })
                    .collect(),
                strand: f.strand.flipped(),
                ..f.clone()
            })
            .collect();
        Part {
            seq: revcomp(&self.seq),
            features,
        }
    }

    /// The last `n` bases are dropped.
    fn trim_end(&self, n: usize) -> Part {
        self.window(0, self.seq.len().saturating_sub(n))
    }

    fn window(&self, start: usize, length: usize) -> Part {
        let end = (start + length).min(self.seq.len());
        let seq = self.seq.get(start..end).unwrap_or_default().to_vec();
        let mut features = Vec::new();
        for f in &self.features {
            let mut clipped = false;
            let segments: Vec<Segment> = f
                .segments
                .iter()
                .filter_map(|s| {
                    let (a, b) = (s.range.start.max(start), s.range.end.min(end));
                    clipped |= a != s.range.start || b != s.range.end;
                    (a < b).then(|| Segment {
                        range: Range::new(a - start, b - start),
                    })
                })
                .collect();
            if segments.is_empty() {
                continue;
            }
            let mut qualifiers = f.qualifiers.clone();
            if clipped && !qualifiers.contains(&partial_note()) {
                qualifiers.push(partial_note());
            }
            features.push(NewFeature {
                segments,
                qualifiers,
                ..f.clone()
            });
        }
        Part { seq, features }
    }

    /// This part followed by `other`.
    fn join(mut self, other: &Part) -> Part {
        let offset = self.seq.len();
        self.seq.extend_from_slice(&other.seq);
        self.features.extend(other.features.iter().map(|f| {
            NewFeature {
                segments: f
                    .segments
                    .iter()
                    .map(|s| Segment {
                        range: Range::new(s.range.start + offset, s.range.end + offset),
                    })
                    .collect(),
                ..f.clone()
            }
        }));
        self
    }
}

/// `left` followed by `right`, with the last `shared` bases of `left` dropped (they are the
/// first bases of `right`).
pub(crate) fn join_trimmed(left: Part, shared: usize, right: &Part) -> Part {
    left.trim_end(shared).join(right)
}

/// `p` without its last `n` bases.
pub(crate) fn trim_end(p: Part, n: usize) -> Part {
    p.trim_end(n)
}

/// Build a molecule from `parts` in order, circular or linear.
pub fn assemble(parts: &[Part], topology: Topology, name: &str) -> Result<Molecule, DocError> {
    let mut all = Part {
        seq: Vec::new(),
        features: Vec::new(),
    };
    for p in parts {
        all = all.join(p);
    }
    let mut m = Molecule::new(name, Seq::new(Alphabet::Dna, all.seq)?, topology)?;
    for f in all.features {
        m.add_feature(f)?;
    }
    Ok(m)
}

/// One end of a double-stranded fragment.
#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
#[serde(tag = "kind", content = "bases", rename_all = "snake_case")]
pub enum End {
    Blunt,
    /// A 5' overhang, its bases as they read on the top strand.
    Overhang5(String),
    /// A 3' overhang, likewise.
    Overhang3(String),
}

impl End {
    fn len(&self) -> usize {
        match self {
            End::Blunt => 0,
            End::Overhang5(s) | End::Overhang3(s) => s.len(),
        }
    }

    /// The same end seen from the other strand (after turning the fragment around).
    pub fn flipped(&self) -> End {
        let rc = |s: &str| String::from_utf8_lossy(&revcomp(s.as_bytes())).into_owned();
        match self {
            End::Blunt => End::Blunt,
            End::Overhang5(s) => End::Overhang5(rc(s)),
            End::Overhang3(s) => End::Overhang3(rc(s)),
        }
    }

    /// Whether a right end and a left end anneal and ligate.
    pub fn ligates(&self, left_of_next: &End) -> bool {
        self == left_of_next
    }

    pub fn describe(&self) -> String {
        match self {
            End::Blunt => "blunt".to_string(),
            End::Overhang5(s) => format!("5' {s}"),
            End::Overhang3(s) => format!("3' {s}"),
        }
    }
}

/// A double-stranded fragment: the bases either strand covers, and its two ends. The overhang
/// bases are part of `part.seq`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DsFragment {
    pub part: Part,
    pub left: End,
    pub right: End,
    /// Where it came from: enzyme at each end (`None` for a molecule end).
    pub left_enzyme: Option<String>,
    pub right_enzyme: Option<String>,
}

impl DsFragment {
    /// Length of the double-stranded part plus overhangs.
    pub fn len(&self) -> usize {
        self.part.seq.len()
    }

    pub fn is_empty(&self) -> bool {
        self.part.seq.is_empty()
    }

    pub fn reversed(&self) -> DsFragment {
        DsFragment {
            part: self.part.reversed(),
            left: self.right.flipped(),
            right: self.left.flipped(),
            left_enzyme: self.right_enzyme.clone(),
            right_enzyme: self.left_enzyme.clone(),
        }
    }
}

/// A cut: the top- and bottom-strand cut positions and the enzyme.
#[derive(Debug, Clone)]
struct Cut {
    top: usize,
    bottom: usize,
    enzyme: String,
}

/// Cut `m` with `enzymes` into double-stranded fragments, in order along the molecule (a
/// circular molecule's last fragment runs across the origin). An uncut molecule gives no
/// fragments.
pub fn digest_fragments(m: &Molecule, enzymes: &[Enzyme]) -> Vec<DsFragment> {
    let len = m.len();
    let circular = m.topology() == Topology::Circular;
    let mut cuts: Vec<Cut> = find_sites(m, enzymes)
        .into_iter()
        .filter_map(|s| {
            Some(Cut {
                top: s.cut?,
                bottom: s.cut_bottom?,
                enzyme: enzymes.get(s.enzyme)?.name.clone(),
            })
        })
        .filter(|c| circular || (c.top.min(c.bottom) > 0 && c.top.max(c.bottom) < len))
        .collect();
    cuts.sort_by_key(|c| c.top);
    cuts.dedup_by(|a, b| a.top == b.top && a.bottom == b.bottom);
    if cuts.is_empty() || len == 0 {
        return Vec::new();
    }
    let seq = m.seq().as_bytes();
    let bases = |from: usize, to: usize| -> String {
        // `from..to` along the top strand, wrapping.
        let n = (to + len - from) % len;
        (0..n)
            .filter_map(|k| seq.get((from + k) % len).map(|&b| b as char))
            .collect()
    };
    // The two strands' cuts of one site, in order along the top strand: a 5' overhang when
    // the top strand is cut first, a 3' overhang when the bottom strand is (measured round a
    // circle, so a site across the origin is read correctly).
    let ends_of = |c: &Cut| -> (usize, usize, End) {
        let d = if circular {
            (c.bottom + len - c.top) % len
        } else {
            c.bottom.wrapping_sub(c.top)
        };
        if c.top == c.bottom {
            (c.top, c.top, End::Blunt)
        } else if (circular && d < len / 2) || (!circular && c.top < c.bottom) {
            (c.top, c.bottom, End::Overhang5(bases(c.top, c.bottom)))
        } else {
            (c.bottom, c.top, End::Overhang3(bases(c.bottom, c.top)))
        }
    };
    let mut out = Vec::new();
    let n = cuts.len();
    let pairs: Vec<(Option<&Cut>, Option<&Cut>)> = if circular {
        (0..n)
            .map(|i| (cuts.get(i), cuts.get((i + 1) % n)))
            .collect()
    } else {
        std::iter::once((None, cuts.first()))
            .chain((0..n.saturating_sub(1)).map(|i| (cuts.get(i), cuts.get(i + 1))))
            .chain(std::iter::once((cuts.last(), None)))
            .collect()
    };
    for (a, b) in pairs {
        let (start, left) = a.map_or((0, End::Blunt), |c| {
            let (first, _, e) = ends_of(c);
            (first, e)
        });
        let (end, right) = b.map_or((len, End::Blunt), |c| {
            let (_, last, e) = ends_of(c);
            (last, e)
        });
        let length = if circular {
            match (end + len - start) % len {
                0 => len,
                l => l,
            }
        } else {
            end.saturating_sub(start)
        };
        // A single cut in a circle opens it: the fragment runs from the cut round to itself,
        // overhang included at both ends.
        let length = if circular && n == 1 {
            len + left.len()
        } else {
            length
        };
        let part = extract_long(m, start, length);
        out.push(DsFragment {
            part,
            left,
            right,
            left_enzyme: a.map(|c| c.enzyme.clone()),
            right_enzyme: b.map(|c| c.enzyme.clone()),
        });
    }
    out
}

/// Like [`extract`], but may run past the full length of a circle (an opened plasmid's
/// overhang repeats at both ends).
fn extract_long(m: &Molecule, start: usize, length: usize) -> Part {
    let len = m.len();
    if length <= len || m.topology() != Topology::Circular {
        return extract(m, start, length);
    }
    let first = extract(m, start, len);
    let extra = extract(m, start, length - len);
    first.join(&Part {
        seq: extra.seq,
        features: Vec::new(),
    })
}

/// Ligate fragments in order; circularise if `circular` (the last fragment's right end must
/// ligate to the first's left end). Returns the joined part, or which junction doesn't fit.
pub fn ligate(fragments: &[DsFragment], circular: bool) -> Result<Part, String> {
    let Some(first) = fragments.first() else {
        return Err("nothing to ligate".to_string());
    };
    let mut part = first.part.clone();
    let mut right = first.right.clone();
    for (i, f) in fragments.iter().enumerate().skip(1) {
        if !right.ligates(&f.left) {
            return Err(format!(
                "fragment {} ends {} but fragment {} starts {}",
                i,
                right.describe(),
                i + 1,
                f.left.describe()
            ));
        }
        // The shared overhang is kept once: drop it from the left fragment.
        part = part.trim_end(right.len()).join(&f.part);
        right = f.right.clone();
    }
    if circular {
        if !right.ligates(&first.left) {
            return Err(format!(
                "the last fragment ends {} but the first starts {}",
                right.describe(),
                first.left.describe()
            ));
        }
        part = part.trim_end(right.len());
    }
    Ok(part)
}

/// How a cut piece can go into a cut backbone.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Fit {
    /// Its ends don't match the backbone's: it won't go in as cut.
    No,
    /// One way round only: the two ends differ, which is what makes a clone directional.
    OneWay,
    /// Either way round, so half the colonies carry it backwards. This is what one enzyme at
    /// both ends, or two blunt ends, gives you.
    EitherWay,
}

/// Whether a piece fits a backbone, and whether its direction is settled. Both orientations
/// are tried, as [`restriction_clone`] does: a piece that ligates both ways has no defined
/// direction in the product.
pub fn how_it_fits(backbone: (&End, &End), piece: (&End, &End)) -> Fit {
    let (bl, br) = backbone;
    let (pl, pr) = piece;
    // The backbone's right end meets the piece's left end, and the piece's right end comes
    // back round to the backbone's left end.
    let forward = br.ligates(pl) && pr.ligates(bl);
    // Turned around, the piece's ends swap and each is read from the other strand.
    let (rl, rr) = (pr.flipped(), pl.flipped());
    let reverse = br.ligates(&rl) && rr.ligates(bl);
    match (forward, reverse) {
        (false, false) => Fit::No,
        (true, true) => Fit::EitherWay,
        _ => Fit::OneWay,
    }
}

/// A construct from restriction cloning.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Construct {
    pub part: Part,
    /// How it was made, in words.
    pub description: String,
}

/// Cut `vector` and `insert` with their enzymes and ligate the vector's largest fragment with
/// each insert fragment whose ends fit, in each orientation that fits. Directional
/// (non-symmetric) clones come first.
pub fn restriction_clone(
    vector: &Molecule,
    vector_enzymes: &[Enzyme],
    insert: &Molecule,
    insert_enzymes: &[Enzyme],
) -> Result<Vec<Construct>, String> {
    let backbone = digest_fragments(vector, vector_enzymes)
        .into_iter()
        .max_by_key(DsFragment::len)
        .ok_or("the enzymes don't cut the vector")?;
    let inserts = if insert_enzymes.is_empty() {
        // An uncut linear insert (blunt PCR product, for example).
        if insert.topology() == Topology::Circular {
            return Err("name enzymes to cut the insert out of a circular molecule".to_string());
        }
        vec![DsFragment {
            part: extract(insert, 0, insert.len()),
            left: End::Blunt,
            right: End::Blunt,
            left_enzyme: None,
            right_enzyme: None,
        }]
    } else {
        let all = digest_fragments(insert, insert_enzymes);
        if all.is_empty() {
            return Err("the enzymes don't cut the insert".to_string());
        }
        all
    };
    let mut out = Vec::new();
    for ins in &inserts {
        for (flip, candidate) in [(false, ins.clone()), (true, ins.reversed())] {
            if let Ok(part) = ligate(&[backbone.clone(), candidate.clone()], true) {
                let ends = |f: &DsFragment| {
                    format!(
                        "{}\u{2013}{}",
                        f.left_enzyme.as_deref().unwrap_or("end"),
                        f.right_enzyme.as_deref().unwrap_or("end")
                    )
                };
                out.push(Construct {
                    part,
                    description: format!(
                        "{} bp insert ({}){} into the {} bp vector backbone ({})",
                        candidate.len(),
                        ends(ins),
                        if flip { " reversed" } else { "" },
                        backbone.len(),
                        ends(&backbone)
                    ),
                });
            }
        }
    }
    if out.is_empty() {
        return Err(format!(
            "no insert fragment fits the vector's ends ({} and {})",
            backbone.left.describe(),
            backbone.right.describe()
        ));
    }
    // An insert that only fits one way is the intended (directional) clone.
    let count = out.len();
    out.sort_by_key(|c| {
        (
            count > 1 && c.description.contains("reversed"),
            std::cmp::Reverse(c.part.seq.len()),
        )
    });
    Ok(out)
}

/// The Gateway att cores (see the module docs).
pub const ATT1_CORE: &[u8] = b"TTTGTACAAAAAAGC";
pub const ATT2_CORE: &[u8] = b"CTTTCTTGTACAAAG";

/// Which kind of att site: the arms flanking the core say so. attB and attP are the parents,
/// attL and attR the children of a BP reaction (Hartley, Temple & Brasch, Genome Res. 10:1788
/// (2000)). Each site is one arm from the phage (P) and one from the bacterium (B), or two of
/// a kind:
///
/// | site | arm facing the cassette | arm facing the backbone |
/// |---|---|---|
/// | attB | B | B |
/// | attL | B | P |
/// | attR | P | B |
/// | attP | P | P |
///
/// An entry clone therefore carries its insert between two B arms, a destination vector its
/// ccdB cassette between two P arms, and only attL and attR recombine (LR), as only attB and
/// attP do (BP).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AttKind {
    B,
    L,
    R,
    P,
}

impl AttKind {
    pub fn label(self) -> &'static str {
        match self {
            AttKind::B => "attB",
            AttKind::L => "attL",
            AttKind::R => "attR",
            AttKind::P => "attP",
        }
    }
}

/// Which arm: from the bacterial attB site, or from the phage attP site.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Arm {
    B,
    P,
}

/// The bases immediately flanking each core in the att sites SnappieSnap ships, read outwards
/// from the core. A B arm is short, so it is compared over its own length.
const ARM1_BACKBONE_B: &[u8] = b"GAACA"; // attB1/attR1, 5' of core 1, read outwards
const ARM1_BACKBONE_P: &[u8] = b"CAACCGTAAT"; // attL1/attP1, likewise
const ARM1_CASSETTE_B: &[u8] = b"AGGCT"; // attB1/attL1, 3' of core 1
const ARM1_CASSETTE_P: &[u8] = b"GAACGAGAA"; // attR1/attP1, likewise (one variable base first)
const ARM2_CASSETTE_B: &[u8] = b"GACCCA"; // attB2/attL2, 5' of core 2, read outwards
const ARM2_CASSETTE_P: &[u8] = b"GACTTGCTCT"; // attR2/attP2, likewise
const ARM2_BACKBONE_B: &[u8] = b"TGGT"; // attB2/attR2, 3' of core 2
const ARM2_BACKBONE_P: &[u8] = b"GGCATTATAA"; // attL2/attP2, likewise (two variable bases first)

/// How many of `pattern` differ from `text` (counting a short text as all different).
fn differences(pattern: &[u8], text: &[u8]) -> usize {
    if text.len() < pattern.len() {
        return pattern.len();
    }
    pattern
        .iter()
        .zip(text)
        .filter(|(a, b)| !a.eq_ignore_ascii_case(b))
        .count()
}

/// Which arm `text` is, or `None` when neither fits clearly. A B arm is only 4-6 bases, so it
/// must match exactly; the longer P arm may differ in one base, which leaves the two apart.
/// `p_skip` passes over the bases where attL and attP (or attR and attP) differ, next to the
/// core, so one constant serves both.
fn arm(text: &[u8], b: &[u8], p: &[u8], p_skip: usize) -> Option<Arm> {
    let (db, dp) = (
        differences(b, text),
        differences(p, text.get(p_skip..).unwrap_or_default()),
    );
    match (db == 0, dp <= 1) {
        (true, false) => Some(Arm::B),
        (false, true) => Some(Arm::P),
        _ => None,
    }
}

/// `n` bases of `s` before `at`, read outwards from `at` (so the base just before `at` comes
/// first); wrapping when the molecule is circular.
fn before(s: &[u8], at: usize, n: usize, circular: bool) -> Vec<u8> {
    let len = s.len();
    if len == 0 {
        return Vec::new();
    }
    let mut out = Vec::with_capacity(n);
    for k in 1..=n {
        let i = match at.checked_sub(k) {
            Some(i) => i,
            None if circular => (at + len - k % len) % len,
            None => break,
        };
        match s.get(i) {
            Some(b) => out.push(*b),
            None => break,
        }
    }
    out
}

/// `n` bases of `s` from `at` onwards, wrapping when the molecule is circular.
fn after(s: &[u8], at: usize, n: usize, circular: bool) -> Vec<u8> {
    let len = s.len();
    if len == 0 {
        return Vec::new();
    }
    let mut out = Vec::with_capacity(n);
    for k in 0..n {
        let i = at.saturating_add(k);
        let i = if i < len {
            i
        } else if circular {
            i % len
        } else {
            break;
        };
        match s.get(i) {
            Some(b) => out.push(*b),
            None => break,
        }
    }
    out
}

/// What kind of att site sits at each core of `m`: site 1 first, then site 2. `None` where the
/// arms are too unusual to tell, which is reported rather than guessed.
pub fn att_kinds(m: &Molecule) -> Result<[Option<AttKind>; 2], String> {
    let (mol, a1, a2) = oriented(m)?;
    let circular = mol.topology() == Topology::Circular;
    let s = mol.seq().as_bytes();
    // Site 1's 3' arm and site 2's 5' arm face the cassette between them.
    let site1 = {
        let backbone = arm(
            &before(s, a1, ARM1_BACKBONE_P.len(), circular),
            ARM1_BACKBONE_B,
            ARM1_BACKBONE_P,
            0,
        );
        let cassette = arm(
            &after(s, a1 + ATT1_CORE.len(), ARM1_CASSETTE_P.len() + 1, circular),
            ARM1_CASSETTE_B,
            ARM1_CASSETTE_P,
            1,
        );
        kind_of(backbone, cassette)
    };
    let site2 = {
        let cassette = arm(
            &before(s, a2, ARM2_CASSETTE_P.len(), circular),
            ARM2_CASSETTE_B,
            ARM2_CASSETTE_P,
            0,
        );
        let backbone = arm(
            &after(s, a2 + ATT2_CORE.len(), ARM2_BACKBONE_P.len() + 2, circular),
            ARM2_BACKBONE_B,
            ARM2_BACKBONE_P,
            2,
        );
        kind_of(backbone, cassette)
    };
    Ok([site1, site2])
}

/// Whether these two molecules can react the way `reaction` says.
///
/// LR joins an entry clone (attL) to a destination vector (attR); BP joins an attB insert to a
/// donor vector (attP). The core is the same in all four kinds of site, so a cassette swap
/// would otherwise "work" between molecules that never recombine at the bench. Sites whose
/// arms are too unusual to name are let through: see [`att_kinds`].
fn check_pair(
    acceptor: &Molecule,
    donor: &Molecule,
    reaction: GatewayReaction,
) -> Result<(), String> {
    let (want_acceptor, want_donor) = match reaction {
        GatewayReaction::Lr => (AttKind::R, AttKind::L),
        GatewayReaction::Bp => (AttKind::P, AttKind::B),
    };
    // The kind of a molecule's sites, when both say the same thing.
    let kind = |m: &Molecule| -> Option<AttKind> {
        match att_kinds(m).ok()? {
            [Some(a), Some(b)] if a == b => Some(a),
            _ => None,
        }
    };
    let (have_acceptor, have_donor) = (kind(acceptor), kind(donor));
    let wrong = |have: Option<AttKind>, want: AttKind| have.is_some_and(|h| h != want);
    if !wrong(have_acceptor, want_acceptor) && !wrong(have_donor, want_donor) {
        return Ok(());
    }
    // Swapping them over would be a real reaction: say so, it is the likeliest slip.
    if have_acceptor == Some(want_donor) && have_donor == Some(want_acceptor) {
        return Err(format!(
            "these are the other way round: {} carries {} sites and {} carries {}, so {} is the {} and {} the {}",
            acceptor.name(),
            want_donor.label(),
            donor.name(),
            want_acceptor.label(),
            donor.name(),
            acceptor_word(reaction),
            acceptor.name(),
            donor_word(reaction),
        ));
    }
    let describe = |m: &Molecule, have: Option<AttKind>, want: AttKind, role: &str| match have {
        Some(h) if h != want => format!(
            "{} carries {} sites, but the {role} of this reaction needs {}",
            m.name(),
            h.label(),
            want.label()
        ),
        _ => String::new(),
    };
    let mut why: Vec<String> = Vec::new();
    for s in [
        describe(
            acceptor,
            have_acceptor,
            want_acceptor,
            acceptor_word(reaction),
        ),
        describe(donor, have_donor, want_donor, donor_word(reaction)),
    ] {
        if !s.is_empty() {
            why.push(s);
        }
    }
    Err(why.join("; "))
}

/// What each side of the reaction is called.
fn acceptor_word(reaction: GatewayReaction) -> &'static str {
    match reaction {
        GatewayReaction::Lr => "destination vector",
        GatewayReaction::Bp => "donor vector",
    }
}

fn donor_word(reaction: GatewayReaction) -> &'static str {
    match reaction {
        GatewayReaction::Lr => "entry clone",
        GatewayReaction::Bp => "attB insert",
    }
}

/// The site's kind from its two arms (see [`AttKind`]).
fn kind_of(backbone: Option<Arm>, cassette: Option<Arm>) -> Option<AttKind> {
    match (backbone?, cassette?) {
        (Arm::B, Arm::B) => Some(AttKind::B),
        (Arm::P, Arm::B) => Some(AttKind::L),
        (Arm::B, Arm::P) => Some(AttKind::R),
        (Arm::P, Arm::P) => Some(AttKind::P),
    }
}

/// Which Gateway reaction: decides how the new sites are named.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, serde::Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum GatewayReaction {
    /// Entry clone (attL) × destination vector (attR) → expression clone (attB).
    Lr,
    /// attB-flanked insert × donor vector (attP) → entry clone (attL).
    Bp,
}

/// Whether two ranges on a molecule of length `len` share a position.
fn overlaps(a: Range, b: Range, len: usize) -> bool {
    let (a1, a2) = a.segments(len);
    let (b1, b2) = b.segments(len);
    std::iter::once(a1).chain(a2).any(|x| {
        std::iter::once(b1.clone())
            .chain(b2.clone())
            .any(|y| x.start < y.end && y.start < x.end)
    })
}

/// Where a core sits: start on the top strand of a (possibly reverse-complemented) molecule.
fn find_unique(seq: &[u8], core: &[u8], circular: bool) -> Result<usize, String> {
    let n = seq.len();
    let k = core.len();
    let starts = if circular { n } else { n.saturating_sub(k) + 1 };
    let mut hits = (0..starts).filter(|&i| {
        (0..k).all(|j| {
            seq.get((i + j) % n.max(1)).map(u8::to_ascii_uppercase) == core.get(j).copied()
        })
    });
    let first = hits.next().ok_or("not found")?;
    if hits.next().is_some() {
        return Err("found more than once".to_string());
    }
    Ok(first)
}

/// A molecule oriented so both cores read on the top strand, with their positions.
fn oriented(m: &Molecule) -> Result<(Molecule, usize, usize), String> {
    let circular = m.topology() == Topology::Circular;
    for flip in [false, true] {
        let mol = if flip {
            let mut r = m.clone();
            r.reverse_complement().map_err(|e| e.to_string())?;
            r
        } else {
            m.clone()
        };
        let s = mol.seq().as_bytes();
        if let (Ok(a), Ok(b)) = (
            find_unique(s, ATT1_CORE, circular),
            find_unique(s, ATT2_CORE, circular),
        ) {
            return Ok((mol, a, b));
        }
    }
    Err(format!(
        "{}: needs exactly one att1 and one att2 site",
        m.name()
    ))
}

/// Gateway recombination: the acceptor (destination or donor vector) keeps its backbone, and
/// the stretch between its att1 and att2 cores is replaced by the donor's (entry clone or attB
/// insert). The product keeps the acceptor's origin and orientation; the new sites are
/// annotated (attB1/attB2 after LR, attL1/attL2 after BP).
pub fn gateway(
    acceptor: &Molecule,
    donor: &Molecule,
    reaction: GatewayReaction,
    name: &str,
) -> Result<Molecule, String> {
    let (acc, a1, a2) = oriented(acceptor)?;
    let (don, d1, d2) = oriented(donor)?;
    check_pair(acceptor, donor, reaction)?;
    let (alen, dlen) = (acc.len(), don.len());
    let k1 = ATT1_CORE.len();
    // Acceptor: from the end of its att1 core round to the start of its att2 core is replaced.
    // The kept backbone runs from the att2 core start to the att1 core end (wrapping).
    let backbone_len = (a1 + k1 + alen - a2) % alen.max(1);
    let insert_start = (d1 + k1) % dlen.max(1);
    let insert_len = (d2 + dlen - insert_start) % dlen.max(1);
    if acc.topology() == Topology::Linear && a2 < a1 {
        return Err("the acceptor's att2 site comes before its att1 site".to_string());
    }
    if insert_len == 0 || insert_len > dlen {
        return Err("the donor's att sites are in the wrong order".to_string());
    }
    let backbone = extract(&acc, a2, backbone_len);
    let insert = extract(&don, insert_start, insert_len);
    // Keep the acceptor's origin: the product reads acceptor[0..a1+k1] + insert + acceptor
    // [a2..]. Build it as backbone (from a2) + insert, then rotate.
    let joined = backbone.join(&insert);
    let total = joined.seq.len();
    let mut m = assemble(&[joined], Topology::Circular, name).map_err(|e| e.to_string())?;
    // Backbone starts at a2; the acceptor's origin is alen - a2 further on.
    let origin = (alen - a2) % total.max(1);
    m.rotate_origin(origin).map_err(|e| e.to_string())?;
    // The new att sites: the 25-bp attB around each core (5 bp before core 1, 6 before core 2),
    // or the attL sites named over the same cores.
    let s = m.seq().as_bytes().to_vec();
    let c1 = find_unique(&s, ATT1_CORE, true)?;
    let c2 = find_unique(&s, ATT2_CORE, true)?;
    let (n1, n2, kind) = match reaction {
        GatewayReaction::Lr => ("attB1", "attB2", "protein_bind"),
        GatewayReaction::Bp => ("attL1", "attL2", "protein_bind"),
    };
    let len = m.len();
    // The old sites, cut through by the crossover, became the new ones: drop their pieces.
    let near_core = |r: Range| {
        [c1, c2].iter().any(|&c| {
            let site = Range::new((c + len - 6) % len, (c + 21) % len);
            overlaps(r, site, len)
        })
    };
    let replaced: Vec<_> = m
        .features()
        .iter()
        .filter(|f| f.qualifiers.contains(&partial_note()))
        .filter(|f| f.segments.iter().any(|s| near_core(s.range)))
        .map(|f| f.id)
        .collect();
    for id in replaced {
        m.remove_feature(id).map_err(|e| e.to_string())?;
    }
    for (name, core, before, strand) in [(n1, c1, 5, Strand::Forward), (n2, c2, 6, Strand::Reverse)]
    {
        let start = (core + len - before) % len;
        let end = (start + 25) % len;
        m.add_feature(NewFeature {
            color: None,
            name: name.to_string(),
            kind: kind.to_string(),
            strand,
            segments: vec![Segment {
                range: Range::new(start, if end == 0 { len } else { end }),
            }],
            qualifiers: vec![Qualifier {
                name: "note".to_string(),
                values: vec![format!(
                    "made by Gateway {} recombination",
                    match reaction {
                        GatewayReaction::Lr => "LR",
                        GatewayReaction::Bp => "BP",
                    }
                )],
            }],
        })
        .map_err(|e| e.to_string())?;
    }
    Ok(m)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn mol(s: &str, t: Topology) -> Molecule {
        Molecule::new("m", Seq::new(Alphabet::Dna, s).unwrap(), t).unwrap()
    }

    fn enz(n: &str, site: &str) -> Enzyme {
        Enzyme::parse(n, site, true).unwrap()
    }

    #[test]
    fn digest_ends() {
        let ecori = enz("EcoRI", "G^AATTC");
        let pst = enz("PstI", "CTGCA^G");
        let s = format!(
            "{}GAATTC{}CTGCAG{}",
            "A".repeat(10),
            "C".repeat(20),
            "T".repeat(10)
        );
        let f = digest_fragments(&mol(&s, Topology::Linear), &[ecori, pst]);
        assert_eq!(f.len(), 3);
        // Middle fragment: AATTC…CTGCA plus overhangs; left 5' AATT, right 3' TGCA.
        assert_eq!(f[1].left, End::Overhang5("AATT".into()));
        assert_eq!(f[1].right, End::Overhang3("TGCA".into()));
        assert_eq!(
            String::from_utf8(f[1].part.seq.clone()).unwrap(),
            format!("AATTC{}CTGCA", "C".repeat(20))
        );
        assert_eq!(f[0].right, End::Overhang5("AATT".into()));
        assert_eq!(f[2].left, End::Overhang3("TGCA".into()));
        // Turning the middle fragment around: ends swap and read on the other strand.
        let r = f[1].reversed();
        assert_eq!(r.left, End::Overhang3("TGCA".into()));
        assert_eq!(r.right, End::Overhang5("AATT".into()));
    }

    #[test]
    fn religating_a_digest_restores_the_molecule() {
        let ecori = enz("EcoRI", "G^AATTC");
        let bam = enz("BamHI", "G^GATCC");
        let s = format!(
            "{}GAATTC{}GGATCC{}",
            "A".repeat(15),
            "C".repeat(30),
            "T".repeat(25)
        );
        let m = mol(&s, Topology::Circular);
        let f = digest_fragments(&m, &[ecori.clone(), bam]);
        assert_eq!(f.len(), 2);
        let back = ligate(&f, true).unwrap();
        assert_eq!(back.seq.len(), s.len());
        // Same circle (it starts at the first cut).
        let doubled = format!("{s}{s}");
        assert!(doubled.contains(&String::from_utf8(back.seq.clone()).unwrap()));
        // A single cut opens the plasmid with the overhang at both ends; closing it again works.
        let one = digest_fragments(&m, &[ecori]);
        assert_eq!(one.len(), 1);
        assert_eq!(one[0].len(), s.len() + 4);
        assert_eq!(ligate(&one, true).unwrap().seq.len(), s.len());
        // Incompatible ends don't ligate.
        let a = DsFragment {
            right: End::Overhang5("AATT".into()),
            ..f[0].clone()
        };
        let b = DsFragment {
            left: End::Overhang5("GATC".into()),
            ..f[1].clone()
        };
        assert!(ligate(&[a, b], false).is_err());
    }

    /// Two different enzymes give a directional clone; one enzyme at both ends, or blunt
    /// ends, let the piece go in either way round, which is why half the colonies carry it
    /// backwards. Ends that don't match don't go in at all.
    #[test]
    fn how_a_piece_fits_its_backbone() {
        let eco = End::Overhang5("AATT".into());
        let bam = End::Overhang5("GATC".into());
        // Backbone cut EcoRI...BamHI takes an EcoRI-BamHI piece one way round only.
        assert_eq!(
            how_it_fits((&eco, &bam), (&bam, &eco)),
            Fit::OneWay,
            "different ends settle the direction"
        );
        // One enzyme at both ends: it goes in either way.
        assert_eq!(
            how_it_fits((&eco, &eco), (&eco, &eco)),
            Fit::EitherWay,
            "the same cut at both ends leaves the direction open"
        );
        // Blunt both ends, likewise.
        assert_eq!(
            how_it_fits((&End::Blunt, &End::Blunt), (&End::Blunt, &End::Blunt)),
            Fit::EitherWay
        );
        // An EcoRI piece into a BamHI backbone: no.
        assert_eq!(how_it_fits((&bam, &bam), (&eco, &eco)), Fit::No);
        // A 5' overhang doesn't anneal to a 3' one of the same bases.
        assert_eq!(
            how_it_fits(
                (
                    &End::Overhang3("AATT".into()),
                    &End::Overhang3("AATT".into())
                ),
                (&eco, &eco)
            ),
            Fit::No
        );
        // The answer agrees with the reaction itself: pUC-style EcoRI/BamHI cloning gives
        // one construct, and EcoRI at both ends gives two (the second one reversed).
        let (eco_e, bam_e) = (enz("EcoRI", "G^AATTC"), enz("BamHI", "G^GATCC"));
        let v = mol(
            &format!("{0}GAATTC{0}GGATCC{0}", "ACGTACGTAC".repeat(12)),
            Topology::Circular,
        );
        let i = mol(
            &format!("TTTT{0}GAATTC{0}GGATCC{0}TTTT", "GTCAGTCAGT".repeat(6)),
            Topology::Linear,
        );
        let enzymes = [eco_e.clone(), bam_e];
        let both = restriction_clone(&v, &enzymes, &i, &enzymes).unwrap();
        assert_eq!(both.len(), 1, "directional: one construct");
        let one = [eco_e];
        // The same enzyme at both ends of the piece, so it has no preferred direction.
        let flanked = mol(
            &format!(
                "TTTT{0}GAATTC{1}GAATTC{0}TTTT",
                "ACGT".repeat(8),
                "GTCAGTCAGT".repeat(6)
            ),
            Topology::Linear,
        );
        let same = restriction_clone(&v, &one, &flanked, &one).unwrap();
        assert_eq!(same.len(), 2, "either way round: two constructs");
        assert!(same.iter().any(|c| c.description.contains("reversed")));
    }

    #[test]
    fn a_site_across_the_origin() {
        // EcoRI's site split by the origin: G|AATTC with the G last. Top cut at 0 (after the
        // G), bottom cut at 4.
        let ecori = enz("EcoRI", "G^AATTC");
        let s = format!("AATTC{}G", "T".repeat(40));
        let f = digest_fragments(&mol(&s, Topology::Circular), std::slice::from_ref(&ecori));
        assert_eq!(f.len(), 1);
        assert_eq!(
            (f[0].left.clone(), f[0].right.clone()),
            (End::Overhang5("AATT".into()), End::Overhang5("AATT".into()))
        );
        assert!(
            String::from_utf8(f[0].part.seq.clone())
                .unwrap()
                .starts_with("AATTC")
        );
        // PstI (3' overhang) split by the origin: CTGCA|G with the cut site across it.
        let pst = enz("PstI", "CTGCA^G");
        let s = format!("GCAG{}CT", "A".repeat(40));
        let f = digest_fragments(&mol(&s, Topology::Circular), &[pst]);
        assert_eq!(f[0].left, End::Overhang3("TGCA".into()));
        assert_eq!(ligate(&f, true).unwrap().seq.len(), s.len());
    }

    #[test]
    fn directional_restriction_cloning() {
        let ecori = enz("EcoRI", "G^AATTC");
        let bam = enz("BamHI", "G^GATCC");
        let vector = mol(
            &format!(
                "{}GAATTC{}GGATCC{}",
                "A".repeat(40),
                "C".repeat(10),
                "G".repeat(50)
            ),
            Topology::Circular,
        );
        let mut insert = mol(
            &format!("TTGAATTC{}GGATCCTT", "ATGCATGCAT".repeat(3)),
            Topology::Linear,
        );
        insert
            .add_feature(NewFeature {
                color: None,
                name: "gene".into(),
                kind: "CDS".into(),
                strand: Strand::Forward,
                segments: vec![Segment {
                    range: Range::new(8, 38),
                }],
                qualifiers: vec![],
            })
            .unwrap();
        let both = [ecori, bam];
        let c = restriction_clone(&vector, &both, &insert, &both).unwrap();
        assert_eq!(c.len(), 1, "{c:?}");
        let product = &c[0].part;
        // Vector minus the small EcoRI-BamHI piece plus the insert's EcoRI-BamHI piece.
        assert_eq!(product.seq.len(), 40 + 6 + 30 + 6 + 50);
        assert!(product.features.iter().any(|f| f.name == "gene"));
        // One enzyme: both orientations fit.
        let one = [enz("EcoRI", "G^AATTC")];
        let ins2 = mol(
            &format!("TTGAATTC{}GAATTCTT", "ATGCATGCAT".repeat(3)),
            Topology::Linear,
        );
        let c = restriction_clone(&vector, &one, &ins2, &one).unwrap();
        assert_eq!(c.len(), 2);
        assert!(c[1].description.contains("reversed"));
        // Nothing fits.
        assert!(restriction_clone(&vector, &both, &ins2, &one).is_err());
    }

    #[test]
    fn gateway_lr_exchanges_the_cassette() {
        // A made-up destination: backbone, attR1-like core, a cassette, att2 core, backbone.
        let dest = format!(
            "{}ACAAGTTTGTACAAAAAAGCTGAAC{}ACCCAGCTTTCTTGTACAAAGTGGT{}",
            "G".repeat(30),
            "CCDB".replace('D', "G").repeat(10),
            "C".repeat(30)
        );
        let entry = format!(
            "{}CAACTTTGTACAAAAAAGCAGGCT{}ACCCAGCTTTCTTGTACAAAGTTGG{}",
            "T".repeat(20),
            "ATGAAAGGGCCC".repeat(5),
            "A".repeat(20)
        );
        let d = mol(&dest, Topology::Circular);
        let e = mol(&entry, Topology::Circular);
        let p = gateway(&d, &e, GatewayReaction::Lr, "expr").unwrap();
        let s = String::from_utf8(p.seq().as_bytes().to_vec()).unwrap();
        // Acceptor up to the end of core 1, then the entry's stretch, then the acceptor from
        // core 2: the cassette is gone and the gene is in.
        assert!(s.starts_with(&"G".repeat(30)));
        assert!(s.contains("ATGAAAGGGCCC") && !s.contains("CCGB"));
        let expected = format!(
            "{}ACAAGTTTGTACAAAAAAGCAGGCT{}ACCCAGCTTTCTTGTACAAAGTGGT{}",
            "G".repeat(30),
            "ATGAAAGGGCCC".repeat(5),
            "C".repeat(30)
        );
        assert_eq!(s, expected);
        let names: Vec<&str> = p.features().iter().map(|f| f.name.as_str()).collect();
        assert_eq!(names, vec!["attB1", "attB2"]);
        assert_eq!(p.features()[0].segments[0].range, Range::new(30, 55));
        // Works whichever way round the donor is stored.
        let mut e_rc = e.clone();
        e_rc.reverse_complement().unwrap();
        assert_eq!(
            gateway(&d, &e_rc, GatewayReaction::Lr, "x").unwrap().seq(),
            p.seq()
        );
        // Missing sites are an error.
        assert!(
            gateway(
                &d,
                &mol("ACGTACGT", Topology::Circular),
                GatewayReaction::Lr,
                "x"
            )
            .is_err()
        );
    }

    #[test]
    fn hostile() {
        let e = [enz("EcoRI", "G^AATTC")];
        for s in ["", "A", "GAATTC", "GAATTCGAATTC"] {
            for t in [Topology::Linear, Topology::Circular] {
                let m = mol(s, t);
                let f = digest_fragments(&m, &e);
                let _ = ligate(&f, true);
                let _ = ligate(&f, false);
                let _ = restriction_clone(&m, &e, &m, &e);
                let _ = gateway(&m, &m, GatewayReaction::Bp, "x");
            }
        }
    }
}

#[cfg(test)]
mod att_kind_tests {
    use super::*;
    use snappiesnap_doc::{Alphabet, Seq};

    /// The eight att sites SnappieSnap ships, pasted into a plasmid, must each be recognised
    /// for what they are: an entry clone reads attL, a destination vector attR, and so on.
    fn plasmid(site1: &str, cassette: &str, site2: &str) -> Molecule {
        let seq = format!(
            "{}{site1}{cassette}{site2}{}",
            "G".repeat(60),
            "C".repeat(60)
        );
        Molecule::new(
            "p",
            Seq::new(Alphabet::Dna, seq.into_bytes()).unwrap(),
            Topology::Circular,
        )
        .unwrap()
    }

    const ATTB1: &str = "ACAAGTTTGTACAAAAAAGCAGGCT";
    const ATTB2: &str = "ACCCAGCTTTCTTGTACAAAGTGGT";
    const ATTL1: &str = "CAAATAATGATTTTATTTTGACTGATAGTGACCTGTTCGTTGCAACAAATTGATGAGCAATGCTTTTTTATAATGCCAACTTTGTACAAAAAAGCAGGCT";
    const ATTL2: &str = "ACCCAGCTTTCTTGTACAAAGTTGGCATTATAAAAAAGCATTGCTTATCAATTTGTTGCAACGAACAGGTCACTATCAGTCAAAATAAAATCATTATTTG";

    const ATTR1: &str = "ACAAGTTTGTACAAAAAAGCTGAACGAGAAACGTAAAATGATATAAATATCAATATATTAAATTAGATTTTGCATAAAAAACAGACTACATAATACTGTAAAACACAACATATCCAGTCACTATG";
    const ATTR2: &str = "CATAGTGACTGGATATGTTGTGTTTTACAGTATTATGTAGTCTGTTTTTTATGCAAAATCTAATTTAATATATTGATATTTATATCATTTTACGTTTCTCGTTCAGCTTTCTTGTACAAAGTGGT";
    const ATTP1: &str = "AAATAATGATTTTATTTTGACTGATAGTGACCTGTTCGTTGCAACAAATTGATGAGCAATGCTTTTTTATAATGCCAAGTTTGTACAAAAAAGCAGAACGAGAAACGTAAAATGATATAAATATCAATATATTAAATTAGATTTTGCATAAAAAACAGACTACATAATACTGTAAAACACAACATATCCAGTCACTATGAATCAACTACTTAGATGGTATTAGTGACCTGTA";
    const ATTP2: &str = "TACAGGTCACTAATACCATCTAAGTAGTTGATTCATAGTGACTGGATATGTTGTGTTTTACAGTATTATGTAGTCTGTTTTTTATGCAAAATCTAATTTAATATATTGATATTTATATCATTTTACGTTTCTCGTTCAGCTTTCTTGTACAAAGTGGGCATTATAAGAAAGCATTGCTTATCAATTTGTTGCAACGAACAGGTCACTATCAGTCAAAATAAAATCATTATTT";

    #[test]
    fn att_sites_are_told_apart_by_their_arms() {
        let expression = plasmid(ATTB1, "ATGAAACCCGGGTTTAAACCCGGG", ATTB2);
        assert_eq!(
            att_kinds(&expression).unwrap(),
            [Some(AttKind::B), Some(AttKind::B)],
            "an expression clone reads attB"
        );
        let entry = plasmid(ATTL1, "ATGAAACCCGGGTTTAAACCCGGG", ATTL2);
        assert_eq!(
            att_kinds(&entry).unwrap(),
            [Some(AttKind::L), Some(AttKind::L)],
            "an entry clone reads attL"
        );
        let destination = plasmid(ATTR1, "ATGAAACCCGGGTTTAAACCCGGG", ATTR2);
        assert_eq!(
            att_kinds(&destination).unwrap(),
            [Some(AttKind::R), Some(AttKind::R)],
            "a destination vector reads attR"
        );
        let donor = plasmid(ATTP1, "ATGAAACCCGGGTTTAAACCCGGG", ATTP2);
        assert_eq!(
            att_kinds(&donor).unwrap(),
            [Some(AttKind::P), Some(AttKind::P)],
            "a donor vector reads attP"
        );
    }

    /// Only attL with attR recombines in an LR reaction, and only attB with attP in a BP one.
    /// The cores are the same in all four, so a pair that cannot react must be refused rather
    /// than quietly swapped.
    #[test]
    fn only_the_pairs_that_recombine_are_allowed() {
        let cassette = "ATGAAACCCGGGTTTAAACCCGGG";
        let entry = plasmid(ATTL1, cassette, ATTL2);
        let destination = plasmid(ATTR1, cassette, ATTR2);
        let expression = plasmid(ATTB1, cassette, ATTB2);
        let donor = plasmid(ATTP1, cassette, ATTP2);

        // The reactions that happen at the bench.
        assert!(gateway(&destination, &entry, GatewayReaction::Lr, "p").is_ok());
        assert!(gateway(&donor, &expression, GatewayReaction::Bp, "p").is_ok());

        // An entry clone and an expression clone share the att cores but cannot recombine.
        let e = gateway(&entry, &expression, GatewayReaction::Lr, "p").unwrap_err();
        assert!(e.contains("attL") && e.contains("attR"), "{e}");
        // The right pair, the wrong way round, says so.
        let e = gateway(&entry, &destination, GatewayReaction::Lr, "p").unwrap_err();
        assert!(e.contains("the other way round"), "{e}");
        // BP needs attP and attB, not the children of an earlier BP.
        let e = gateway(&destination, &entry, GatewayReaction::Bp, "p").unwrap_err();
        assert!(e.contains("attP"), "{e}");
    }
}
