//! Multiple sequence alignment by progressive profile alignment.
//!
//! 1. Distances from shared k-mers (Edgar R. C. (2004) BMC Bioinformatics 5:113, the k-mer
//!    distance of MUSCLE's first stage): `1 - F`, with `F` the shared k-mers (counted with
//!    multiplicity) over the most the shorter sequence could share.
//! 2. A guide tree by UPGMA (Sokal R. R., Michener C. D. (1958) Univ. Kansas Sci. Bull.
//!    38:1409-1438).
//! 3. Profiles merged along the tree, closest first (Feng D.-F., Doolittle R. F. (1987) J. Mol.
//!    Evol. 25:351-360; Thompson J. D., Higgins D. G., Gibson T. J. (1994) Nucleic Acids Res.
//!    22:4673-4680), each merge a global affine-gap alignment of two profiles ([`crate::align`])
//!    scored by the average sum of pairs of their columns: nucleotides with the blastn scores,
//!    proteins with BLOSUM62. Gaps in a column score zero against anything.
//!
//! Every row of the result, without its gaps, is exactly the sequence it came from.

use serde::Serialize;

use crate::align::{Mode, Op, ProteinScoring, Scoring, align_core, blosum62, pair};

/// Most sequences aligned at once.
pub const MAX_SEQUENCES: usize = 500;

/// Gap character in the rows.
pub const GAP: u8 = b'-';

/// Scores are scaled to integers with this factor (column scores are averages).
const SCALE: f64 = 100.0;

/// Nucleotides or amino acids.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum Kind {
    Dna,
    Protein,
}

impl Kind {
    /// Nucleotides if every residue is a nucleotide code and at least 90% are A, C, G, T, U or
    /// N (many amino-acid letters are also IUPAC nucleotide codes), otherwise protein.
    pub fn guess(seqs: &[&[u8]]) -> Kind {
        let code = |b: &u8| b"ACGTURYSWKMBDHVN".contains(&b.to_ascii_uppercase());
        let plain = |b: &&u8| b"ACGTUN".contains(&b.to_ascii_uppercase());
        let all: usize = seqs.iter().map(|s| s.len()).sum();
        let plain_n: usize = seqs.iter().map(|s| s.iter().filter(plain).count()).sum();
        if seqs.iter().all(|s| s.iter().all(code)) && plain_n * 10 >= all * 9 {
            Kind::Dna
        } else {
            Kind::Protein
        }
    }

    /// Symbols scored (an index per residue class).
    fn symbols(self) -> &'static [u8] {
        match self {
            Kind::Dna => b"ACGTN",
            Kind::Protein => b"ARNDCQEGHILKMFPSTWYVBZX*",
        }
    }

    fn index(self, b: u8) -> usize {
        let b = b.to_ascii_uppercase();
        let b = if b == b'U' { b'T' } else { b };
        let syms = self.symbols();
        syms.iter().position(|&c| c == b).unwrap_or(
            syms.len()
                .saturating_sub(if self == Kind::Dna { 1 } else { 2 }),
        )
    }

    fn gaps(self) -> (i32, i32) {
        match self {
            Kind::Dna => {
                let s = Scoring::default();
                (s.gap_open, s.gap_extend)
            }
            Kind::Protein => {
                let s = ProteinScoring::default();
                (s.gap_open, s.gap_extend)
            }
        }
    }

    fn k(self) -> usize {
        match self {
            Kind::Dna => 6,
            Kind::Protein => 3,
        }
    }
}

/// An alignment of several sequences.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct Msa {
    pub kind: Kind,
    /// One row per input sequence, in input order, all the same length, `-` for gaps.
    pub rows: Vec<String>,
    /// Per column: the residue most sequences have there if more than half do, `-` if more
    /// than half have a gap, otherwise `N` (nucleotides) or `X` (proteins).
    pub consensus: String,
    /// Per row: residues equal to the consensus, over the row's residues.
    pub identity: Vec<f64>,
    /// Columns where every sequence has the same residue (no gaps).
    pub identical_columns: usize,
}

/// A group of aligned rows and its column profile.
struct Group {
    /// Input indices of the rows.
    members: Vec<usize>,
    rows: Vec<Vec<u8>>,
}

impl Group {
    fn columns(&self) -> usize {
        self.rows.first().map_or(0, Vec::len)
    }

    /// Per column: (symbol, count) of the residues in it.
    fn profile(&self, kind: Kind) -> Vec<Vec<(usize, f64)>> {
        let k = kind.symbols().len();
        (0..self.columns())
            .map(|c| {
                let mut counts = vec![0.0f64; k];
                for r in &self.rows {
                    if let Some(&b) = r.get(c)
                        && b != GAP
                        && let Some(x) = counts.get_mut(kind.index(b))
                    {
                        *x += 1.0;
                    }
                }
                counts
                    .into_iter()
                    .enumerate()
                    .filter(|(_, n)| *n > 0.0)
                    .collect()
            })
            .collect()
    }
}

/// Score of symbol `x` against symbol `y`.
fn matrix(kind: Kind) -> Vec<Vec<f64>> {
    let syms = kind.symbols();
    let s = Scoring::default();
    syms.iter()
        .map(|&a| {
            syms.iter()
                .map(|&b| match kind {
                    Kind::Dna => f64::from(pair(a, b, &s).0),
                    Kind::Protein => f64::from(blosum62(a, b)),
                })
                .collect()
        })
        .collect()
}

/// Shared k-mer distance of two sequences (0 identical .. 1 nothing shared).
fn kmer_distance(a: &[u8], b: &[u8], k: usize) -> f64 {
    use std::collections::HashMap;
    let count = |s: &[u8]| {
        let mut m: HashMap<Vec<u8>, u32> = HashMap::new();
        for w in s.windows(k) {
            *m.entry(w.iter().map(u8::to_ascii_uppercase).collect())
                .or_insert(0) += 1;
        }
        m
    };
    let most = a.len().min(b.len()).saturating_sub(k) + 1;
    if a.len() < k || b.len() < k || most == 0 {
        return 1.0;
    }
    let (ca, cb) = (count(a), count(b));
    let shared: u32 = ca
        .iter()
        .map(|(w, n)| (*n).min(cb.get(w).copied().unwrap_or(0)))
        .sum();
    1.0 - (f64::from(shared) / most as f64).min(1.0)
}

/// Align `seqs` (at least two).
pub fn align_many(seqs: &[&[u8]], kind: Kind) -> Result<Msa, String> {
    let n = seqs.len();
    if n < 2 {
        return Err("give at least two sequences".into());
    }
    if n > MAX_SEQUENCES {
        return Err(format!("at most {MAX_SEQUENCES} sequences at once"));
    }
    if let Some(i) = seqs.iter().position(|s| s.is_empty()) {
        return Err(format!("sequence {} is empty", i + 1));
    }
    let k = kind.k();
    // Distances, then UPGMA: repeatedly merge the closest two clusters.
    let mut dist = vec![vec![0.0f64; n]; n];
    for i in 0..n {
        for j in (i + 1)..n {
            let d = match (seqs.get(i), seqs.get(j)) {
                (Some(a), Some(b)) => kmer_distance(a, b, k),
                _ => 1.0,
            };
            if let Some(row) = dist.get_mut(i)
                && let Some(x) = row.get_mut(j)
            {
                *x = d;
            }
            if let Some(row) = dist.get_mut(j)
                && let Some(x) = row.get_mut(i)
            {
                *x = d;
            }
        }
    }
    let mut groups: Vec<Option<Group>> = seqs
        .iter()
        .enumerate()
        .map(|(i, s)| {
            Some(Group {
                members: vec![i],
                rows: vec![s.to_vec()],
            })
        })
        .collect();
    let mut sizes = vec![1usize; n];
    let scores = matrix(kind);
    let (gap_open, gap_extend) = kind.gaps();
    let syms = kind.symbols().len();
    for _ in 1..n {
        // The closest pair of live clusters (lowest indices on ties).
        let mut best: Option<(f64, usize, usize)> = None;
        for i in 0..n {
            if groups.get(i).is_none_or(Option::is_none) {
                continue;
            }
            for j in (i + 1)..n {
                if groups.get(j).is_none_or(Option::is_none) {
                    continue;
                }
                let d = dist.get(i).and_then(|r| r.get(j)).copied().unwrap_or(1.0);
                if best.is_none_or(|(b, _, _)| d < b) {
                    best = Some((d, i, j));
                }
            }
        }
        let Some((_, i, j)) = best else { break };
        let (Some(a), Some(b)) = (
            groups.get_mut(i).and_then(Option::take),
            groups.get_mut(j).and_then(Option::take),
        ) else {
            return Err("internal: a cluster went missing".into());
        };
        let merged = merge(a, b, kind, &scores, syms, gap_open, gap_extend)?;
        // UPGMA: the new cluster's distance is the size-weighted average.
        let (si, sj) = (
            sizes.get(i).copied().unwrap_or(1),
            sizes.get(j).copied().unwrap_or(1),
        );
        for x in 0..n {
            if x == i || groups.get(x).is_none_or(Option::is_none) {
                continue;
            }
            let di = dist.get(i).and_then(|r| r.get(x)).copied().unwrap_or(1.0);
            let dj = dist.get(j).and_then(|r| r.get(x)).copied().unwrap_or(1.0);
            let d = (di * si as f64 + dj * sj as f64) / (si + sj) as f64;
            if let Some(r) = dist.get_mut(i)
                && let Some(v) = r.get_mut(x)
            {
                *v = d;
            }
            if let Some(r) = dist.get_mut(x)
                && let Some(v) = r.get_mut(i)
            {
                *v = d;
            }
        }
        if let Some(s) = sizes.get_mut(i) {
            *s = si + sj;
        }
        if let Some(g) = groups.get_mut(i) {
            *g = Some(merged);
        }
    }
    let all = groups
        .into_iter()
        .flatten()
        .next()
        .ok_or_else(|| "internal: nothing aligned".to_string())?;
    // Rows back in input order.
    let mut rows = vec![Vec::new(); n];
    for (member, row) in all.members.iter().zip(all.rows) {
        if let Some(slot) = rows.get_mut(*member) {
            *slot = row;
        }
    }
    Ok(summarise(kind, rows))
}

/// Align two groups' profiles and join their rows.
fn merge(
    a: Group,
    b: Group,
    kind: Kind,
    scores: &[Vec<f64>],
    syms: usize,
    gap_open: i32,
    gap_extend: i32,
) -> Result<Group, String> {
    let (pa, pb) = (a.profile(kind), b.profile(kind));
    let (na, nb) = (a.rows.len() as f64, b.rows.len() as f64);
    // Each column of `b` mixed with the scores: weight[x] = sum_y count_b[y] * S[x][y] / nb.
    let weights: Vec<Vec<f64>> = pb
        .iter()
        .map(|col| {
            (0..syms)
                .map(|x| {
                    col.iter()
                        .map(|&(y, c)| {
                            c * scores.get(x).and_then(|r| r.get(y)).copied().unwrap_or(0.0)
                        })
                        .sum::<f64>()
                        / nb
                })
                .collect()
        })
        .collect();
    let scale = |v: i32| v.saturating_mul(SCALE as i32);
    let alignment = align_core(
        pa.len(),
        pb.len(),
        Mode::Overlap,
        scale(gap_open),
        scale(gap_extend),
        |i, j| {
            let s: f64 = match (pa.get(i), weights.get(j)) {
                (Some(col), Some(w)) => {
                    col.iter()
                        .map(|&(x, c)| c * w.get(x).copied().unwrap_or(0.0))
                        .sum::<f64>()
                        / na
                }
                _ => 0.0,
            };
            ((s * SCALE).round() as i32, Op::Match)
        },
    )?;
    let mut rows_a: Vec<Vec<u8>> = vec![Vec::new(); a.rows.len()];
    let mut rows_b: Vec<Vec<u8>> = vec![Vec::new(); b.rows.len()];
    let column = |rows: &[Vec<u8>], c: usize, out: &mut [Vec<u8>]| {
        for (r, o) in rows.iter().zip(out.iter_mut()) {
            o.push(r.get(c).copied().unwrap_or(GAP));
        }
    };
    let gaps = |out: &mut [Vec<u8>]| {
        for o in out.iter_mut() {
            o.push(GAP);
        }
    };
    // Leading ends (free end gaps): each against gaps.
    for c in 0..alignment.reference_start {
        column(&a.rows, c, &mut rows_a);
        gaps(&mut rows_b);
    }
    for c in 0..alignment.query_start {
        gaps(&mut rows_a);
        column(&b.rows, c, &mut rows_b);
    }
    let (mut i, mut j) = (alignment.reference_start, alignment.query_start);
    for &(op, k) in &alignment.ops {
        for _ in 0..k {
            match op {
                Op::Match | Op::Mismatch => {
                    column(&a.rows, i, &mut rows_a);
                    column(&b.rows, j, &mut rows_b);
                    i += 1;
                    j += 1;
                }
                Op::Deletion => {
                    column(&a.rows, i, &mut rows_a);
                    gaps(&mut rows_b);
                    i += 1;
                }
                Op::Insertion => {
                    gaps(&mut rows_a);
                    column(&b.rows, j, &mut rows_b);
                    j += 1;
                }
            }
        }
    }
    // Trailing ends.
    for c in i..a.columns() {
        column(&a.rows, c, &mut rows_a);
        gaps(&mut rows_b);
    }
    for c in j..b.columns() {
        gaps(&mut rows_a);
        column(&b.rows, c, &mut rows_b);
    }
    let mut members = a.members;
    members.extend(b.members);
    rows_a.extend(rows_b);
    Ok(Group {
        members,
        rows: rows_a,
    })
}

fn summarise(kind: Kind, rows: Vec<Vec<u8>>) -> Msa {
    let width = rows.first().map_or(0, Vec::len);
    let n = rows.len();
    let unknown = if kind == Kind::Dna { b'N' } else { b'X' };
    let mut consensus = Vec::with_capacity(width);
    let mut identical_columns = 0usize;
    for c in 0..width {
        let mut counts: Vec<(u8, usize)> = Vec::new();
        let mut gaps = 0usize;
        for r in &rows {
            match r.get(c).copied() {
                Some(GAP) | None => gaps += 1,
                Some(b) => {
                    let b = b.to_ascii_uppercase();
                    match counts.iter_mut().find(|(x, _)| *x == b) {
                        Some((_, k)) => *k += 1,
                        None => counts.push((b, 1)),
                    }
                }
            }
        }
        let top = counts.iter().max_by_key(|(_, k)| *k).copied();
        if gaps == 0 && counts.len() == 1 {
            identical_columns += 1;
        }
        consensus.push(match top {
            _ if gaps * 2 > n => GAP,
            Some((b, k)) if k * 2 > n => b,
            _ => unknown,
        });
    }
    let identity = rows
        .iter()
        .map(|r| {
            let (mut same, mut residues) = (0usize, 0usize);
            for (&b, &c) in r.iter().zip(&consensus) {
                if b != GAP {
                    residues += 1;
                    if b.eq_ignore_ascii_case(&c) {
                        same += 1;
                    }
                }
            }
            if residues == 0 {
                0.0
            } else {
                same as f64 / residues as f64
            }
        })
        .collect();
    Msa {
        kind,
        rows: rows
            .into_iter()
            .map(|r| String::from_utf8_lossy(&r).into_owned())
            .collect(),
        consensus: String::from_utf8_lossy(&consensus).into_owned(),
        identity,
        identical_columns,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn ungapped(row: &str) -> String {
        row.chars().filter(|c| *c != '-').collect()
    }

    fn check(seqs: &[&str], kind: Kind) -> Msa {
        let input: Vec<&[u8]> = seqs.iter().map(|s| s.as_bytes()).collect();
        let m = align_many(&input, kind).unwrap();
        assert_eq!(m.rows.len(), seqs.len());
        let width = m.rows[0].len();
        for (row, s) in m.rows.iter().zip(seqs) {
            assert_eq!(row.len(), width);
            assert_eq!(ungapped(row), *s, "a row must be its sequence");
        }
        assert_eq!(m.consensus.len(), width);
        // No column is all gaps.
        for c in 0..width {
            assert!(m.rows.iter().any(|r| r.as_bytes()[c] != GAP), "column {c}");
        }
        m
    }

    #[test]
    fn identical_sequences_need_no_gaps() {
        let s = "ACGTTGCAAGGCTTACGATCGGATCC";
        let m = check(&[s, s, s], Kind::Dna);
        assert!(m.rows.iter().all(|r| r == s));
        assert_eq!(m.consensus, s);
        assert_eq!(m.identical_columns, s.len());
        assert!(m.identity.iter().all(|i| (*i - 1.0).abs() < 1e-12));
    }

    #[test]
    fn known_indels_line_up() {
        // A base deleted in one, a base inserted in another, a substitution in a third.
        let base = "ATGGCTAGCAAAGGAGAAGAACTTTTCACTGGAGTTGTCCCAATTCTTGTTGAATTAGATGGTGATGTTAATGG";
        let del = base.replacen("CTTTTCACT", "CTTTCACT", 1);
        let ins = base.replacen("GGAGTTGTC", "GGAGTTTGTC", 1);
        let sub = base.replacen("AATTCTTG", "AATACTTG", 1);
        let m = check(&[base, &del, &ins, &sub], Kind::Dna);
        // One extra column for the insertion; the deletion is one gap in its row.
        assert_eq!(m.rows[0].len(), base.len() + 1);
        assert_eq!(m.rows[1].matches('-').count(), 2);
        assert_eq!(m.rows[2].matches('-').count(), 0);
        assert_eq!(m.rows[3].matches('-').count(), 1);
        // The consensus is the base sequence plus a gap at the insertion.
        assert_eq!(ungapped(&m.consensus), base);
    }

    #[test]
    fn different_lengths_and_order_are_kept() {
        let long = "GGGGGATGCGTACGTTAGCCGATCGGATCCTTAGGCATGCAAGCTTGGGGG";
        let short = "ACGTTAGCCGATCGGATCCTTAGG";
        let mid = "ATGCGTACGTTAGCCGATCGGATCCTTAGGCATGCAAG";
        let m = check(&[short, long, mid], Kind::Dna);
        assert_eq!(m.rows[1], long);
        // The short one sits in the middle, flanked by gaps.
        assert!(m.rows[0].starts_with("-----") && m.rows[0].ends_with("-----"));
    }

    #[test]
    fn proteins_align_with_blosum62() {
        let seqs = [
            "MKTAYIAKQRQISFVKSHFSRQ",
            "MKTAYVAKQRQISFVKSHFSRQ",
            "MKTAYIAKQRISFVKSHFSRQ",
        ];
        let m = check(&seqs, Kind::Protein);
        assert_eq!(m.rows[2].matches('-').count(), 1);
        assert_eq!(Kind::guess(&[b"ACGT", b"acgu"]), Kind::Dna);
        assert_eq!(Kind::guess(&[b"ACGT", b"MKV"]), Kind::Protein);
    }

    #[test]
    fn random_families_keep_every_sequence() {
        let mut seed = 0x0a11_9ce5_u64;
        let mut next = |k: u64| {
            seed = seed
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (seed >> 33) % k
        };
        for _ in 0..20 {
            let base: Vec<u8> = (0..(20 + next(60) as usize))
                .map(|_| b"ACGT"[next(4) as usize])
                .collect();
            let family: Vec<String> = (0..(2 + next(6) as usize))
                .map(|_| {
                    let mut s = base.clone();
                    for _ in 0..next(6) {
                        let k = next(s.len() as u64) as usize;
                        match next(3) {
                            0 => s[k] = b"ACGT"[next(4) as usize],
                            1 if s.len() > 1 => {
                                s.remove(k);
                            }
                            _ => s.insert(k, b"ACGT"[next(4) as usize]),
                        }
                    }
                    String::from_utf8(s).unwrap()
                })
                .collect();
            let refs: Vec<&str> = family.iter().map(String::as_str).collect();
            check(&refs, Kind::Dna);
        }
    }

    #[test]
    fn bad_input_is_an_error() {
        assert!(align_many(&[b"ACGT"], Kind::Dna).is_err());
        assert!(align_many(&[b"ACGT", b""], Kind::Dna).is_err());
        // Very short sequences still align (no shared k-mers).
        let m = align_many(&[b"A", b"C", b"AC"], Kind::Dna).unwrap();
        assert_eq!(m.rows.len(), 3);
    }
}
