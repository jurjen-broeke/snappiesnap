//! Pairwise alignment of nucleotide sequences with affine gap costs.
//!
//! - Global alignment: Needleman S. B., Wunsch C. D. (1970) J. Mol. Biol. 48:443-453.
//! - Local alignment: Smith T. F., Waterman M. S. (1981) J. Mol. Biol. 147:195-197.
//! - Affine gaps in O(nm) time with three matrices: Gotoh O. (1982) J. Mol. Biol. 162:705-708.
//!
//! A gap of `k` bases costs `gap_open + k * gap_extend`. The default scores are NCBI blastn's
//! (reward 2, penalty -3, gap existence 5, gap extension 2; BLAST+ user manual). IUPAC codes
//! score 0 where they could stand for the other base, and as a mismatch where they can't.
//!
//! Besides global and local alignment, [`Mode::Fit`] places the whole query inside the
//! reference (end gaps in the reference are free): the right mode for a read or a construct
//! against a longer reference.
//!
//! Scores are kept for two rows only; the traceback takes one byte per cell, and alignments
//! larger than [`MAX_CELLS`] are refused rather than run out of memory.

use serde::Serialize;

/// Largest `(reference + 1) * (query + 1)` we align (one byte each).
pub const MAX_CELLS: usize = 200_000_000;

/// How the ends are treated.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum Mode {
    /// Both sequences end to end.
    Global,
    /// The best-scoring stretch of each.
    Local,
    /// The whole query, anywhere in the reference.
    Fit,
    /// Both end to end, but gaps at either end of either sequence are free (they overlap,
    /// or one lies inside the other): for merging profiles of sequences of different lengths.
    Overlap,
}

/// Scores (gap costs are positive numbers that are subtracted).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
pub struct Scoring {
    pub matched: i32,
    pub mismatch: i32,
    pub gap_open: i32,
    pub gap_extend: i32,
}

impl Default for Scoring {
    fn default() -> Self {
        Scoring {
            matched: 2,
            mismatch: -3,
            gap_open: 5,
            gap_extend: 2,
        }
    }
}

/// One column kind of an alignment.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum Op {
    /// The same base.
    Match,
    /// Different bases (or an ambiguity code).
    Mismatch,
    /// A query base with no reference base.
    Insertion,
    /// A reference base with no query base.
    Deletion,
}

/// An alignment of a query to a reference. Ranges are 0-based, half-open.
#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
pub struct Alignment {
    pub score: i32,
    pub reference_start: usize,
    pub reference_end: usize,
    pub query_start: usize,
    pub query_end: usize,
    /// Columns, run-length encoded.
    pub ops: Vec<(Op, usize)>,
}

impl Alignment {
    fn count(&self, op: Op) -> usize {
        self.ops
            .iter()
            .filter(|(o, _)| *o == op)
            .map(|(_, n)| n)
            .sum()
    }
    pub fn matches(&self) -> usize {
        self.count(Op::Match)
    }
    pub fn mismatches(&self) -> usize {
        self.count(Op::Mismatch)
    }
    pub fn insertions(&self) -> usize {
        self.count(Op::Insertion)
    }
    pub fn deletions(&self) -> usize {
        self.count(Op::Deletion)
    }
    /// Alignment columns.
    pub fn columns(&self) -> usize {
        self.ops.iter().map(|(_, n)| n).sum()
    }
    /// Identical columns over all columns (0 for an empty alignment).
    pub fn identity(&self) -> f64 {
        let c = self.columns();
        if c == 0 {
            0.0
        } else {
            self.matches() as f64 / c as f64
        }
    }
    /// The two aligned rows, with `-` for gaps.
    pub fn rows(&self, reference: &[u8], query: &[u8]) -> (String, String) {
        let (mut r, mut q) = (String::new(), String::new());
        let (mut i, mut j) = (self.reference_start, self.query_start);
        let at = |s: &[u8], k: usize| s.get(k).map_or('?', |&b| b as char);
        for &(op, n) in &self.ops {
            for _ in 0..n {
                match op {
                    Op::Match | Op::Mismatch => {
                        r.push(at(reference, i));
                        q.push(at(query, j));
                        i += 1;
                        j += 1;
                    }
                    Op::Insertion => {
                        r.push('-');
                        q.push(at(query, j));
                        j += 1;
                    }
                    Op::Deletion => {
                        r.push(at(reference, i));
                        q.push('-');
                        i += 1;
                    }
                }
            }
        }
        (r, q)
    }
}

/// IUPAC nucleotide code as a set of bases (A=1, C=2, G=4, T/U=8); 0 if not a code.
fn bases(b: u8) -> u8 {
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

/// Score and column kind of two residues.
pub(crate) fn pair(a: u8, b: u8, s: &Scoring) -> (i32, Op) {
    let (x, y) = (bases(a), bases(b));
    if x != 0 && x == y && x.count_ones() == 1 {
        (s.matched, Op::Match)
    } else if x & y != 0 {
        // An ambiguity that could be the same base: neither rewarded nor penalised.
        (0, Op::Mismatch)
    } else {
        (s.mismatch, Op::Mismatch)
    }
}

const NEG: i32 = i32::MIN / 4;

// Traceback byte: bits 0-1 where H came from, bit 2 E extended, bit 3 F extended.
const FROM_DIAG: u8 = 0;
const FROM_E: u8 = 1;
const FROM_F: u8 = 2;
const FROM_START: u8 = 3;
const E_EXTEND: u8 = 4;
const F_EXTEND: u8 = 8;

/// Align nucleotide `query` to `reference`.
pub fn align(
    reference: &[u8],
    query: &[u8],
    mode: Mode,
    scoring: &Scoring,
) -> Result<Alignment, String> {
    align_by(
        reference,
        query,
        mode,
        scoring.gap_open,
        scoring.gap_extend,
        |a, b| pair(a, b, scoring),
    )
}

/// The alignment itself, with any residue scores (`pair` gives the score and column kind).
fn align_by(
    reference: &[u8],
    query: &[u8],
    mode: Mode,
    gap_open: i32,
    gap_extend: i32,
    pair: impl Fn(u8, u8) -> (i32, Op),
) -> Result<Alignment, String> {
    let at = |s: &[u8], k: usize| s.get(k).copied().unwrap_or(b'N');
    align_core(
        reference.len(),
        query.len(),
        mode,
        gap_open,
        gap_extend,
        |i, j| pair(at(reference, i), at(query, j)),
    )
}

/// The alignment of `n` reference positions with `m` query positions, where `pair(i, j)` scores
/// putting reference position `i` against query position `j` (0-based): residues, or the
/// columns of two profiles (see [`crate::msa`]).
pub(crate) fn align_core(
    n: usize,
    m: usize,
    mode: Mode,
    gap_open: i32,
    gap_extend: i32,
    pair: impl Fn(usize, usize) -> (i32, Op),
) -> Result<Alignment, String> {
    if gap_open < 0 || gap_extend < 0 {
        return Err("gap costs must not be negative".into());
    }
    let cells = (n + 1)
        .checked_mul(m + 1)
        .filter(|c| *c <= MAX_CELLS)
        .ok_or_else(|| {
            format!(
                "{n} x {m} residues is too large to align at once (at most {MAX_CELLS} cells); align a shorter stretch"
            )
        })?;
    let width = m + 1;
    let mut trace = vec![0u8; cells];
    let gap = |k: usize| -> i32 {
        let k = i32::try_from(k).unwrap_or(i32::MAX / 4);
        gap_open.saturating_add(k.saturating_mul(gap_extend))
    };
    let open = gap_open.saturating_add(gap_extend);
    let ext = gap_extend;

    // Row 0.
    let mut h_prev = vec![0i32; width];
    let mut e_row = vec![NEG; width];
    for j in 1..=m {
        let (score, t) = match mode {
            Mode::Local | Mode::Overlap => (0, FROM_START),
            Mode::Global | Mode::Fit => (-gap(j), FROM_F | if j > 1 { F_EXTEND } else { 0 }),
        };
        if let Some(h) = h_prev.get_mut(j) {
            *h = score;
        }
        if let Some(c) = trace.get_mut(j) {
            *c = t;
        }
    }
    if let Some(c) = trace.get_mut(0) {
        *c = FROM_START;
    }

    let mut best = (NEG, 0usize, 0usize);
    if mode == Mode::Local {
        best = (0, 0, 0);
    }
    if matches!(mode, Mode::Fit | Mode::Overlap) && m == 0 {
        best = (0, 0, 0);
    }
    let mut h_cur = vec![0i32; width];
    for i in 1..=n {
        let row = i * width;
        // Column 0.
        let (score, t) = match mode {
            Mode::Global => (-gap(i), FROM_E | if i > 1 { E_EXTEND } else { 0 }),
            Mode::Local | Mode::Fit | Mode::Overlap => (0, FROM_START),
        };
        if let Some(h) = h_cur.get_mut(0) {
            *h = score;
        }
        if let Some(e) = e_row.get_mut(0) {
            *e = if mode == Mode::Global { score } else { NEG };
        }
        if let Some(c) = trace.get_mut(row) {
            *c = t;
        }
        let mut f = NEG;
        for j in 1..=m {
            let mut t = 0u8;
            // E: a gap in the query (reference base i-1 unpaired), from the row above.
            let up_h = h_prev.get(j).copied().unwrap_or(NEG);
            let up_e = e_row.get(j).copied().unwrap_or(NEG);
            let e = if up_e.saturating_sub(ext) > up_h.saturating_sub(open) {
                t |= E_EXTEND;
                up_e.saturating_sub(ext)
            } else {
                up_h.saturating_sub(open)
            };
            // F: a gap in the reference (query base j-1 unpaired), from the left.
            let left_h = h_cur.get(j - 1).copied().unwrap_or(NEG);
            f = if f.saturating_sub(ext) > left_h.saturating_sub(open) {
                t |= F_EXTEND;
                f.saturating_sub(ext)
            } else {
                left_h.saturating_sub(open)
            };
            let diag = h_prev
                .get(j - 1)
                .copied()
                .unwrap_or(NEG)
                .saturating_add(pair(i - 1, j - 1).0);
            let (mut h, mut from) = (diag, FROM_DIAG);
            if e > h {
                (h, from) = (e, FROM_E);
            }
            if f > h {
                (h, from) = (f, FROM_F);
            }
            if mode == Mode::Local && h <= 0 {
                (h, from) = (0, FROM_START);
            }
            t |= from;
            if let Some(c) = trace.get_mut(row + j) {
                *c = t;
            }
            if let Some(x) = h_cur.get_mut(j) {
                *x = h;
            }
            if let Some(x) = e_row.get_mut(j) {
                *x = e;
            }
            if mode == Mode::Local && h > best.0 {
                best = (h, i, j);
            }
        }
        if matches!(mode, Mode::Fit | Mode::Overlap) {
            let h = h_cur.get(m).copied().unwrap_or(NEG);
            if h > best.0 {
                best = (h, i, m);
            }
        }
        std::mem::swap(&mut h_prev, &mut h_cur);
    }
    if mode == Mode::Overlap {
        // The last row too: the query may run on past the reference's end.
        for j in 0..=m {
            let h = h_prev.get(j).copied().unwrap_or(NEG);
            if h > best.0 {
                best = (h, n, j);
            }
        }
        if n == 0 && best.0 == NEG {
            best = (0, 0, 0);
        }
    }
    let (score, end_i, end_j) = match mode {
        Mode::Global => (h_prev.get(m).copied().unwrap_or(0), n, m),
        Mode::Fit if n == 0 => (h_prev.get(m).copied().unwrap_or(0), 0, m),
        Mode::Overlap if n == 0 || m == 0 => (0, n, 0),
        _ => best,
    };

    // Trace back from the end cell.
    #[derive(Clone, Copy, PartialEq)]
    enum State {
        H,
        E,
        F,
    }
    let (mut i, mut j) = (end_i, end_j);
    let mut state = State::H;
    let mut ops: Vec<(Op, usize)> = Vec::new();
    let mut push = |op: Op| match ops.last_mut() {
        Some((o, k)) if *o == op => *k += 1,
        _ => ops.push((op, 1)),
    };
    // Each step moves at least one cell up or left, or ends: at most n + m + 1 steps.
    for _ in 0..=(n + m + 1) {
        let t = trace.get(i * width + j).copied().unwrap_or(FROM_START);
        match state {
            State::H => match t & 3 {
                FROM_DIAG if i > 0 && j > 0 => {
                    let (_, op) = pair(i - 1, j - 1);
                    push(op);
                    i -= 1;
                    j -= 1;
                }
                FROM_E if i > 0 => state = State::E,
                FROM_F if j > 0 => state = State::F,
                _ => break,
            },
            State::E => {
                if i == 0 {
                    break;
                }
                push(Op::Deletion);
                if t & E_EXTEND == 0 {
                    state = State::H;
                }
                i -= 1;
            }
            State::F => {
                if j == 0 {
                    break;
                }
                push(Op::Insertion);
                if t & F_EXTEND == 0 {
                    state = State::H;
                }
                j -= 1;
            }
        }
    }
    ops.reverse();
    Ok(Alignment {
        score,
        reference_start: i,
        reference_end: end_i,
        query_start: j,
        query_end: end_j,
        ops,
    })
}

// --- Proteins ----------------------------------------------------------------------------

/// Symbols of [`BLOSUM62`], in its row and column order.
const BLOSUM_ORDER: &[u8; 24] = b"ARNDCQEGHILKMFPSTWYVBZX*";

/// BLOSUM62 (Henikoff S., Henikoff J. G. (1992) Proc. Natl. Acad. Sci. USA 89:10915-10919), in
/// the 24-symbol layout NCBI distributes: the 20 amino acids, B (D or N), Z (E or Q), X (any)
/// and `*` (stop).
#[rustfmt::skip]
const BLOSUM62: [[i8; 24]; 24] = [
    //A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V   B   Z   X   *
    [ 4, -1, -2, -2,  0, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -3, -2,  0, -2, -1,  0, -4], // A
    [-1,  5,  0, -2, -3,  1,  0, -2,  0, -3, -2,  2, -1, -3, -2, -1, -1, -3, -2, -3, -1,  0, -1, -4], // R
    [-2,  0,  6,  1, -3,  0,  0,  0,  1, -3, -3,  0, -2, -3, -2,  1,  0, -4, -2, -3,  3,  0, -1, -4], // N
    [-2, -2,  1,  6, -3,  0,  2, -1, -1, -3, -4, -1, -3, -3, -1,  0, -1, -4, -3, -3,  4,  1, -1, -4], // D
    [ 0, -3, -3, -3,  9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1, -3, -3, -2, -4], // C
    [-1,  1,  0,  0, -3,  5,  2, -2,  0, -3, -2,  1,  0, -3, -1,  0, -1, -2, -1, -2,  0,  3, -1, -4], // Q
    [-1,  0,  0,  2, -4,  2,  5, -2,  0, -3, -3,  1, -2, -3, -1,  0, -1, -3, -2, -2,  1,  4, -1, -4], // E
    [ 0, -2,  0, -1, -3, -2, -2,  6, -2, -4, -4, -2, -3, -3, -2,  0, -2, -2, -3, -3, -1, -2, -1, -4], // G
    [-2,  0,  1, -1, -3,  0,  0, -2,  8, -3, -3, -1, -2, -1, -2, -1, -2, -2,  2, -3,  0,  0, -1, -4], // H
    [-1, -3, -3, -3, -1, -3, -3, -4, -3,  4,  2, -3,  1,  0, -3, -2, -1, -3, -1,  3, -3, -3, -1, -4], // I
    [-1, -2, -3, -4, -1, -2, -3, -4, -3,  2,  4, -2,  2,  0, -3, -2, -1, -2, -1,  1, -4, -3, -1, -4], // L
    [-1,  2,  0, -1, -3,  1,  1, -2, -1, -3, -2,  5, -1, -3, -1,  0, -1, -3, -2, -2,  0,  1, -1, -4], // K
    [-1, -1, -2, -3, -1,  0, -2, -3, -2,  1,  2, -1,  5,  0, -2, -1, -1, -1, -1,  1, -3, -1, -1, -4], // M
    [-2, -3, -3, -3, -2, -3, -3, -3, -1,  0,  0, -3,  0,  6, -4, -2, -2,  1,  3, -1, -3, -3, -1, -4], // F
    [-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4,  7, -1, -1, -4, -3, -2, -2, -1, -2, -4], // P
    [ 1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -2,  0, -1, -2, -1,  4,  1, -3, -2, -2,  0,  0,  0, -4], // S
    [ 0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  1,  5, -2, -2,  0, -1, -1,  0, -4], // T
    [-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1,  1, -4, -3, -2, 11,  2, -3, -4, -3, -2, -4], // W
    [-2, -2, -2, -3, -2, -1, -2, -3,  2, -1, -1, -2, -1,  3, -3, -2, -2,  2,  7, -1, -3, -2, -1, -4], // Y
    [ 0, -3, -3, -3, -1, -2, -2, -3, -3,  3,  1, -2,  1, -1, -2, -2,  0, -3, -1,  4, -3, -2, -1, -4], // V
    [-2, -1,  3,  4, -3,  0,  1, -1,  0, -3, -4,  0, -3, -3, -2,  0, -1, -4, -3, -3,  4,  1, -1, -4], // B
    [-1,  0,  0,  1, -3,  3,  4, -2,  0, -3, -3,  1, -1, -3, -1,  0, -1, -3, -2, -2,  1,  4, -1, -4], // Z
    [ 0, -1, -1, -1, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -2,  0,  0, -2, -1, -1, -1, -1, -1, -4], // X
    [-4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4,  1], // *
];

/// The BLOSUM62 row of a residue (any other letter counts as X).
fn blosum_index(b: u8) -> usize {
    let b = b.to_ascii_uppercase();
    BLOSUM_ORDER.iter().position(|&c| c == b).unwrap_or(22)
}

/// The BLOSUM62 score of two residues.
pub fn blosum62(a: u8, b: u8) -> i32 {
    BLOSUM62
        .get(blosum_index(a))
        .and_then(|row| row.get(blosum_index(b)))
        .map_or(-1, |&v| i32::from(v))
}

/// Gap costs for proteins: NCBI BLASTP's defaults (existence 11, extension 1).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
pub struct ProteinScoring {
    pub gap_open: i32,
    pub gap_extend: i32,
}

impl Default for ProteinScoring {
    fn default() -> Self {
        ProteinScoring {
            gap_open: 11,
            gap_extend: 1,
        }
    }
}

/// Align protein `query` to `reference` with BLOSUM62. Identical residues are matches; other
/// pairs are mismatches (see [`positives`] for similar ones).
pub fn align_protein(
    reference: &[u8],
    query: &[u8],
    mode: Mode,
    scoring: &ProteinScoring,
) -> Result<Alignment, String> {
    align_by(
        reference,
        query,
        mode,
        scoring.gap_open,
        scoring.gap_extend,
        |a, b| {
            let op = if a.eq_ignore_ascii_case(&b) {
                Op::Match
            } else {
                Op::Mismatch
            };
            (blosum62(a, b), op)
        },
    )
}

/// Aligned residue pairs that score above zero in BLOSUM62 (identical or similar).
pub fn positives(a: &Alignment, reference: &[u8], query: &[u8]) -> usize {
    let (mut i, mut j, mut n) = (a.reference_start, a.query_start, 0usize);
    for &(op, k) in &a.ops {
        match op {
            Op::Match | Op::Mismatch => {
                for d in 0..k {
                    if let (Some(&x), Some(&y)) = (reference.get(i + d), query.get(j + d))
                        && blosum62(x, y) > 0
                    {
                        n += 1;
                    }
                }
                i += k;
                j += k;
            }
            Op::Insertion => j += k,
            Op::Deletion => i += k,
        }
    }
    n
}

// --- Aligning to a molecule ----------------------------------------------------------------

/// Most differences listed for one alignment.
pub const MAX_DIFFERENCES: usize = 10_000;

/// Which strand of the reference the query matched.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum Strand {
    Forward,
    Reverse,
}

/// One difference between the query and the reference, in reference coordinates (0-based;
/// bases as they read on the reference's top strand).
#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub enum Difference {
    /// The query has `query` where the reference has `reference`, at `position`.
    Substitution {
        position: usize,
        reference: String,
        query: String,
    },
    /// The query has extra `bases` just before reference `position`.
    Insertion { position: usize, bases: String },
    /// The query lacks reference `start..start + length` (`bases`).
    Deletion {
        start: usize,
        length: usize,
        bases: String,
    },
}

/// A query aligned to a molecule.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct MoleculeAlignment {
    pub strand: Strand,
    pub score: i32,
    /// Where the query lies on the molecule: `start..end` (0-based, half-open; on a circular
    /// molecule `end` may be smaller than `start` when it crosses the origin, and `end ==
    /// start` with `length == len` is the whole circle).
    pub start: usize,
    pub end: usize,
    /// Reference bases covered.
    pub length: usize,
    /// The aligned part of the query (in the query's own orientation).
    pub query_start: usize,
    pub query_end: usize,
    pub query_length: usize,
    pub matches: usize,
    pub mismatches: usize,
    pub insertions: usize,
    pub deletions: usize,
    pub identity: f64,
    pub differences: Vec<Difference>,
    /// Whether the list was cut at [`MAX_DIFFERENCES`].
    pub differences_truncated: bool,
    /// The aligned rows (reference top strand, and the query as it pairs with it).
    pub reference_row: String,
    pub query_row: String,
    /// Proteins: aligned pairs that score above zero in BLOSUM62 (identical or similar).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub positives: Option<usize>,
}

fn revcomp(s: &[u8]) -> Vec<u8> {
    s.iter()
        .rev()
        .map(|&b| {
            let c = match b.to_ascii_uppercase() {
                b'A' => b'T',
                b'T' | b'U' => b'A',
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
            };
            if b.is_ascii_lowercase() {
                c.to_ascii_lowercase()
            } else {
                c
            }
        })
        .collect()
}

/// Align `query` to `molecule`, on both strands (the better one is kept, forward on a tie).
/// On a circular molecule the query may cross the origin.
pub fn align_to_molecule(
    molecule: &snappiesnap_doc::Molecule,
    query: &[u8],
    mode: Mode,
    scoring: &Scoring,
) -> Result<MoleculeAlignment, String> {
    let seq = molecule.seq().as_bytes();
    let len = seq.len();
    let circular = molecule.topology() == snappiesnap_doc::Topology::Circular;
    // A circular reference continues past its origin, far enough for the query to cross it.
    let mut reference = seq.to_vec();
    if circular && mode != Mode::Global {
        let extra = query
            .len()
            .saturating_add(query.len() / 10)
            .saturating_add(50)
            .min(len);
        reference.extend_from_slice(seq.get(..extra).unwrap_or_default());
    }
    let forward = align(&reference, query, mode, scoring)?;
    let rc = revcomp(query);
    let reverse = align(&reference, &rc, mode, scoring)?;
    let (strand, a, q) = if reverse.score > forward.score {
        (Strand::Reverse, reverse, rc.as_slice())
    } else {
        (Strand::Forward, forward, query)
    };
    Ok(report(
        &a,
        &reference,
        q,
        query.len(),
        strand,
        len,
        circular,
        None,
    ))
}

/// Align protein `query` to protein `reference` (a protein molecule, or a CDS translation).
pub fn align_protein_report(
    reference: &[u8],
    query: &[u8],
    mode: Mode,
    scoring: &ProteinScoring,
) -> Result<MoleculeAlignment, String> {
    let a = align_protein(reference, query, mode, scoring)?;
    let positives = positives(&a, reference, query);
    Ok(report(
        &a,
        reference,
        query,
        query.len(),
        Strand::Forward,
        reference.len(),
        false,
        Some(positives),
    ))
}

/// The alignment `a` of `q` (as aligned: reverse-complemented for the reverse strand) to
/// `reference` (which may run past the origin of a circular molecule of length `len`), in the
/// molecule's coordinates.
#[allow(clippy::too_many_arguments)]
fn report(
    a: &Alignment,
    reference: &[u8],
    q: &[u8],
    query_len: usize,
    strand: Strand,
    len: usize,
    circular: bool,
    positives: Option<usize>,
) -> MoleculeAlignment {
    let wrap = |p: usize| if len == 0 { 0 } else { p % len };
    let mut differences = Vec::new();
    let mut truncated = false;
    let (mut i, mut j) = (a.reference_start, a.query_start);
    let text = |s: &[u8], from: usize, n: usize| {
        String::from_utf8_lossy(s.get(from..from + n).unwrap_or_default()).into_owned()
    };
    for &(op, n) in &a.ops {
        let room = differences.len() < MAX_DIFFERENCES;
        match op {
            Op::Match => {
                i += n;
                j += n;
            }
            Op::Mismatch => {
                for k in 0..n {
                    if differences.len() >= MAX_DIFFERENCES {
                        truncated = true;
                        break;
                    }
                    differences.push(Difference::Substitution {
                        position: wrap(i + k),
                        reference: text(reference, i + k, 1),
                        query: text(q, j + k, 1),
                    });
                }
                i += n;
                j += n;
            }
            Op::Insertion => {
                if room {
                    differences.push(Difference::Insertion {
                        position: wrap(i),
                        bases: text(q, j, n),
                    });
                } else {
                    truncated = true;
                }
                j += n;
            }
            Op::Deletion => {
                if room {
                    differences.push(Difference::Deletion {
                        start: wrap(i),
                        length: n,
                        bases: text(reference, i, n),
                    });
                } else {
                    truncated = true;
                }
                i += n;
            }
        }
    }
    let covered = a.reference_end - a.reference_start;
    let (reference_row, query_row) = a.rows(reference, q);
    // Query coordinates in the query's own orientation.
    let (qs, qe) = match strand {
        Strand::Forward => (a.query_start, a.query_end),
        Strand::Reverse => (
            query_len.saturating_sub(a.query_end),
            query_len.saturating_sub(a.query_start),
        ),
    };
    MoleculeAlignment {
        strand,
        score: a.score,
        start: wrap(a.reference_start),
        end: if circular {
            wrap(a.reference_end)
        } else {
            a.reference_end
        },
        length: covered,
        query_start: qs,
        query_end: qe,
        query_length: query_len,
        matches: a.matches(),
        mismatches: a.mismatches(),
        insertions: a.insertions(),
        deletions: a.deletions(),
        identity: a.identity(),
        differences,
        differences_truncated: truncated,
        reference_row,
        query_row,
        positives,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rows(r: &str, q: &str, mode: Mode) -> (String, String, Alignment) {
        let a = align(r.as_bytes(), q.as_bytes(), mode, &Scoring::default()).unwrap();
        let (x, y) = a.rows(r.as_bytes(), q.as_bytes());
        (x, y, a)
    }

    #[test]
    fn identical_sequences() {
        let (r, q, a) = rows("ACGTACGT", "ACGTACGT", Mode::Global);
        assert_eq!((r.as_str(), q.as_str()), ("ACGTACGT", "ACGTACGT"));
        assert_eq!(a.score, 16);
        assert_eq!(a.ops, vec![(Op::Match, 8)]);
        assert!((a.identity() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn substitution_insertion_and_deletion_by_hand() {
        // One substitution: 7 matches (14) and one mismatch (-3).
        let (_, _, a) = rows("ACGTACGT", "ACGAACGT", Mode::Global);
        assert_eq!(a.score, 11);
        assert_eq!(a.mismatches(), 1);
        // One base deleted from the query: 11 matches (22) and a 1-base gap (5 + 2).
        let (r, q, a) = rows("AAAACCCCGGGG", "AAAACCCGGGG", Mode::Global);
        assert_eq!(a.score, 22 - 7);
        assert_eq!((a.deletions(), a.insertions()), (1, 0));
        assert_eq!(r.len(), q.len());
        assert!(q.contains('-'));
        // Two extra query bases: one gap of 2 (5 + 4), 12 matches.
        let (r, _, a) = rows("AAAACCCCGGGG", "AAAACCTTCCGGGG", Mode::Global);
        assert_eq!(a.score, 24 - 9);
        assert_eq!(a.insertions(), 2);
        assert_eq!(r, "AAAACC--CCGGGG");
    }

    #[test]
    fn affine_gaps_prefer_one_long_gap() {
        // Deleting 4 bases as one gap (5 + 8 = 13) beats two gaps (2 x (5 + 4) = 18).
        let (_, q, a) = rows("ACGTTTTTTGCA", "ACGTTGCA", Mode::Global);
        assert_eq!(a.deletions(), 4);
        assert_eq!(q.matches("-").count(), 4);
        assert_eq!(
            a.ops.iter().filter(|(o, _)| *o == Op::Deletion).count(),
            1,
            "one gap: {q}"
        );
    }

    #[test]
    fn fit_places_the_query_inside_the_reference() {
        let reference = "TTTTTTTTTTACGGATCCAGTTTTTTTTTT";
        let (_, _, a) = rows(reference, "ACGGATCCAG", Mode::Fit);
        assert_eq!((a.reference_start, a.reference_end), (10, 20));
        assert_eq!((a.query_start, a.query_end), (0, 10));
        assert_eq!(a.score, 20);
        // Global would pay for the reference's ends instead.
        let (_, _, g) = rows(reference, "ACGGATCCAG", Mode::Global);
        assert!(g.score < 0);
    }

    #[test]
    fn local_finds_the_shared_stretch() {
        let (_, _, a) = rows("GGGGGGACGTACGTGGGGG", "CCCCACGTACGTCCCC", Mode::Local);
        assert_eq!((a.reference_start, a.reference_end), (6, 14));
        assert_eq!((a.query_start, a.query_end), (4, 12));
        assert_eq!(a.score, 16);
        // Nothing in common: an empty alignment with score 0.
        let (_, _, a) = rows("AAAA", "CCCC", Mode::Local);
        assert_eq!((a.score, a.columns()), (0, 0));
    }

    #[test]
    fn ambiguity_codes_are_neutral_where_they_fit() {
        let (_, _, a) = rows("ACGTACGT", "ACGNACGT", Mode::Global);
        assert_eq!(a.score, 14);
        assert_eq!(a.mismatches(), 1);
        let (_, _, a) = rows("ACGTACGT", "ACGRACGT", Mode::Global);
        // R (A or G) can't be T.
        assert_eq!(a.score, 11);
    }

    #[test]
    fn empty_and_one_base_inputs() {
        for mode in [Mode::Global, Mode::Local, Mode::Fit] {
            for (r, q) in [("", ""), ("A", ""), ("", "A"), ("A", "A"), ("A", "C")] {
                let a = align(r.as_bytes(), q.as_bytes(), mode, &Scoring::default()).unwrap();
                let (x, y) = a.rows(r.as_bytes(), q.as_bytes());
                assert_eq!(x.len(), y.len(), "{mode:?} {r:?} {q:?}");
                assert_eq!(
                    a.query_end - a.query_start,
                    y.chars().filter(|c| *c != '-').count()
                );
                if mode != Mode::Local {
                    assert_eq!(
                        (a.query_start, a.query_end),
                        (0, q.len()),
                        "{mode:?} {r:?} {q:?}"
                    );
                }
            }
        }
        let a = align(b"", b"ACG", Mode::Global, &Scoring::default()).unwrap();
        assert_eq!((a.score, a.insertions()), (-11, 3));
    }

    #[test]
    fn too_large_or_bad_scores_are_errors() {
        let big = vec![b'A'; 20_000];
        assert!(align(&big, &big, Mode::Global, &Scoring::default()).is_err());
        let bad = Scoring {
            gap_open: -1,
            ..Scoring::default()
        };
        assert!(align(b"A", b"A", Mode::Global, &bad).is_err());
    }

    /// The score of an alignment, recomputed from its columns.
    fn rescore(a: &Alignment, r: &[u8], q: &[u8], s: &Scoring) -> i32 {
        let (mut i, mut j, mut total) = (a.reference_start, a.query_start, 0i32);
        for &(op, n) in &a.ops {
            match op {
                Op::Match | Op::Mismatch => {
                    for _ in 0..n {
                        total += pair(r[i], q[j], s).0;
                        i += 1;
                        j += 1;
                    }
                }
                Op::Insertion => {
                    total -= s.gap_open + n as i32 * s.gap_extend;
                    j += n;
                }
                Op::Deletion => {
                    total -= s.gap_open + n as i32 * s.gap_extend;
                    i += n;
                }
            }
        }
        assert_eq!((i, j), (a.reference_end, a.query_end));
        total
    }

    /// The best global score by exhaustive recursion (affine gaps): an independent oracle for
    /// short sequences. `gap` is the kind of gap open at this point (0 none, 1 E, 2 F).
    fn brute(r: &[u8], q: &[u8], gap: u8, s: &Scoring) -> i32 {
        let mut best = i32::MIN / 2;
        if r.is_empty() && q.is_empty() {
            return 0;
        }
        if let (Some(&a), Some(&b)) = (r.first(), q.first()) {
            best = best.max(pair(a, b, s).0 + brute(&r[1..], &q[1..], 0, s));
        }
        if !r.is_empty() {
            let cost = if gap == 1 {
                s.gap_extend
            } else {
                s.gap_open + s.gap_extend
            };
            best = best.max(brute(&r[1..], q, 1, s) - cost);
        }
        if !q.is_empty() {
            let cost = if gap == 2 {
                s.gap_extend
            } else {
                s.gap_open + s.gap_extend
            };
            best = best.max(brute(r, &q[1..], 2, s) - cost);
        }
        best
    }

    #[test]
    fn random_pairs_agree_with_the_oracle_and_their_own_columns() {
        let mut seed = 0x5eed_1234_u64;
        let mut next = |k: u64| {
            seed = seed
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (seed >> 33) % k
        };
        let s = Scoring::default();
        for _ in 0..400 {
            let (lr, lq) = (next(7) as usize, next(7) as usize);
            let r: Vec<u8> = (0..lr).map(|_| b"ACGTN"[next(5) as usize]).collect();
            let q: Vec<u8> = (0..lq).map(|_| b"ACGT"[next(4) as usize]).collect();
            let g = align(&r, &q, Mode::Global, &s).unwrap();
            assert_eq!(g.score, brute(&r, &q, 0, &s), "{r:?} {q:?}");
            for mode in [Mode::Global, Mode::Local, Mode::Fit, Mode::Overlap] {
                let a = align(&r, &q, mode, &s).unwrap();
                assert_eq!(rescore(&a, &r, &q, &s), a.score, "{mode:?} {r:?} {q:?}");
            }
        }
        // Longer pairs: the columns still add up to the score.
        for _ in 0..50 {
            let r: Vec<u8> = (0..300).map(|_| b"ACGT"[next(4) as usize]).collect();
            let mut q = r[50..250].to_vec();
            for _ in 0..10 {
                let k = next(q.len() as u64) as usize;
                match next(3) {
                    0 => q[k] = b"ACGT"[next(4) as usize],
                    1 => {
                        q.remove(k);
                    }
                    _ => q.insert(k, b'A'),
                }
            }
            for mode in [Mode::Global, Mode::Local, Mode::Fit, Mode::Overlap] {
                let a = align(&r, &q, mode, &s).unwrap();
                assert_eq!(rescore(&a, &r, &q, &s), a.score);
            }
            let fit = align(&r, &q, Mode::Fit, &s).unwrap();
            assert!(
                fit.reference_start >= 30 && fit.reference_end <= 270,
                "{fit:?}"
            );
        }
    }

    #[test]
    fn reads_align_to_either_strand_and_across_the_origin() {
        use snappiesnap_doc::{Alphabet, Molecule, Seq, Topology};
        let plasmid = "ATGCGTACGTTAGCCGATCGGATCCTTAGGCATGCAAGCTTGACCTGAAGTCCGATAC";
        let m = Molecule::new(
            "p",
            Seq::new(Alphabet::Dna, plasmid).unwrap(),
            Topology::Circular,
        )
        .unwrap();
        let s = Scoring::default();
        // A read across the origin with one substitution, one insertion, one deletion.
        let mut read = format!("{}{}", &plasmid[40..], &plasmid[..20]).into_bytes();
        // position 5 of the read = plasmid 45: substitute
        read[5] = if read[5] == b'A' { b'C' } else { b'A' };
        read.insert(12, b'G'); // before read base 12 = plasmid 52
        read.remove(26); // read 26 (after the insert) = plasmid 7 (the plasmid is 58 bp)
        let a = align_to_molecule(&m, &read, Mode::Fit, &s).unwrap();
        assert_eq!(a.strand, Strand::Forward);
        assert_eq!((a.start, a.end), (40, 20));
        assert_eq!((a.mismatches, a.insertions, a.deletions), (1, 1, 1));
        assert!(matches!(
            a.differences[0],
            Difference::Substitution { position: 45, .. }
        ));
        assert!(matches!(
            a.differences[1],
            Difference::Insertion { position: 52, .. }
        ));
        assert!(matches!(
            a.differences[2],
            Difference::Deletion {
                start: 7,
                length: 1,
                ..
            }
        ));
        // The same read from the other strand.
        let rc = revcomp(&read);
        let b = align_to_molecule(&m, &rc, Mode::Fit, &s).unwrap();
        assert_eq!(b.strand, Strand::Reverse);
        assert_eq!((b.start, b.end, b.score), (a.start, a.end, a.score));
        assert_eq!((b.query_start, b.query_end), (0, rc.len()));
        assert_eq!(b.differences, a.differences);
        // A linear molecule doesn't wrap.
        let lin = Molecule::new(
            "l",
            Seq::new(Alphabet::Dna, plasmid).unwrap(),
            Topology::Linear,
        )
        .unwrap();
        let c = align_to_molecule(&lin, &plasmid.as_bytes()[10..30], Mode::Fit, &s).unwrap();
        assert_eq!((c.start, c.end, c.identity), (10, 30, 1.0));
    }

    #[test]
    fn a_whole_rotated_plasmid_aligns_end_to_end() {
        use snappiesnap_doc::{Alphabet, Molecule, Seq, Topology};
        let plasmid = "ATGCGTACGTTAGCCGATCGGATCCTTAGGCATGCAAGCTTGACCTGAAGTCCGATAC";
        let m = Molecule::new(
            "p",
            Seq::new(Alphabet::Dna, plasmid).unwrap(),
            Topology::Circular,
        )
        .unwrap();
        let rotated = format!("{}{}", &plasmid[17..], &plasmid[..17]);
        let a = align_to_molecule(&m, rotated.as_bytes(), Mode::Fit, &Scoring::default()).unwrap();
        assert_eq!((a.start, a.end, a.length), (17, 17, plasmid.len()));
        assert!(a.differences.is_empty());
    }

    #[test]
    fn blosum62_is_the_published_table() {
        // Symmetric, with the published diagonal.
        for (i, row) in BLOSUM62.iter().enumerate() {
            for (j, v) in row.iter().enumerate() {
                assert_eq!(*v, BLOSUM62[j][i], "{i} {j}");
            }
        }
        let diagonal: Vec<i8> = (0..24).map(|i| BLOSUM62[i][i]).collect();
        assert_eq!(
            diagonal,
            vec![
                4, 5, 6, 6, 9, 5, 5, 6, 8, 4, 4, 5, 5, 6, 7, 4, 5, 11, 7, 4, 4, 4, -1, 1
            ]
        );
        // A few well-known pairs.
        assert_eq!(blosum62(b'W', b'C'), -2);
        assert_eq!(blosum62(b'I', b'V'), 3);
        assert_eq!(blosum62(b'D', b'E'), 2);
        assert_eq!(blosum62(b'F', b'Y'), 3);
        assert_eq!(blosum62(b'a', b'S'), 1);
        // Unknown letters count as X.
        assert_eq!(blosum62(b'U', b'A'), 0);
        assert_eq!(blosum62(b'*', b'*'), 1);
    }

    #[test]
    fn proteins_align_with_blosum62() {
        let s = ProteinScoring::default();
        // Identical: the sum of the diagonal.
        let p = b"MKTAYIAKQR";
        let a = align_protein(p, p, Mode::Global, &s).unwrap();
        let diag: i32 = p.iter().map(|&c| blosum62(c, c)).sum();
        assert_eq!(a.score, diag);
        // One residue missing: one gap of 1 (11 + 1).
        let a = align_protein(p, b"MKTAYAKQR", Mode::Global, &s).unwrap();
        assert_eq!(a.score, diag - blosum62(b'I', b'I') - 12);
        assert_eq!(a.deletions(), 1);
        // A conservative change (I -> V) is a positive mismatch.
        let a = align_protein(p, b"MKTAYVAKQR", Mode::Global, &s).unwrap();
        assert_eq!(a.mismatches(), 1);
        assert_eq!(positives(&a, p, b"MKTAYVAKQR"), 10);
        // Local: the shared domain inside unrelated ends.
        let a = align_protein(b"GGGGGWWCHHKKGGGGG", b"PPPWWCHHKKPPP", Mode::Local, &s).unwrap();
        assert_eq!((a.reference_start, a.reference_end), (5, 12));
    }

    #[test]
    fn overlap_ends_are_free_on_both_sides() {
        // The query's start overlaps the reference's end.
        let (_, _, a) = rows("TTTTTTACGTACGT", "ACGTACGTGGGGGG", Mode::Overlap);
        assert_eq!((a.reference_start, a.reference_end), (6, 14));
        assert_eq!((a.query_start, a.query_end), (0, 8));
        assert_eq!(a.score, 16);
        // One inside the other, either way round.
        let (_, _, a) = rows("ACGT", "GGACGTGG", Mode::Overlap);
        assert_eq!((a.query_start, a.query_end, a.score), (2, 6, 8));
        let (_, _, a) = rows("GGACGTGG", "ACGT", Mode::Overlap);
        assert_eq!((a.reference_start, a.reference_end, a.score), (2, 6, 8));
        for (r, q) in [("", ""), ("A", ""), ("", "A")] {
            let a = align(
                r.as_bytes(),
                q.as_bytes(),
                Mode::Overlap,
                &Scoring::default(),
            )
            .unwrap();
            assert_eq!(a.score, 0);
        }
    }
}
