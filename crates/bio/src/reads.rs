//! Checking a construct with sequencing reads: each read is trimmed of its low-quality ends
//! ([`crate::trace::quality_trim`]), aligned to the construct on both strands
//! ([`crate::align::align_to_molecule`], whole read inside the construct), and its differences
//! get the quality of the read's own call. The reads are then pooled: which stretches of the
//! construct they cover, and for each difference how many covering reads show it.
//!
//! A difference is **confirmed** when every read covering it shows it and there are at least
//! two; **one read** when only one read covers it; **disputed** when some covering reads show
//! the construct's base instead (often a poor call in one read). Judging the calls is left to
//! the user, with the traces; nothing here changes the construct.

use std::collections::BTreeMap;

use serde::Serialize;
use snappiesnap_doc::{Molecule, Topology};

use crate::align::{Difference, Mode, MoleculeAlignment, Scoring, Strand, align_to_molecule};
use crate::trace::{DEFAULT_CUTOFF, quality_trim};

/// Most reads checked at once.
pub const MAX_READS: usize = 200;
/// Longest construct checked (coverage is counted per base).
pub const MAX_CONSTRUCT: usize = 20_000_000;

/// A read to check.
#[derive(Debug, Clone, PartialEq)]
pub struct Read {
    pub name: String,
    /// Bases as called, in the read's own orientation.
    pub residues: Vec<u8>,
    /// Phred quality per base, if the read came from a trace.
    pub quality: Option<Vec<Option<u8>>>,
    /// The stretch to keep (0-based, end exclusive), when it was already chosen, as for reads
    /// matched in a `.dna` file; it replaces the quality trimming.
    pub keep: Option<(usize, usize)>,
}

/// One difference as one read shows it.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct ReadDifference {
    #[serde(flatten)]
    pub difference: Difference,
    /// Where it is in the read (0-based, the read's own orientation, untrimmed): the base
    /// itself, the first inserted base, or for a deletion the base before the gap.
    pub read_position: Option<usize>,
    /// Quality of that call (for a deletion, the lower of the bases around the gap).
    pub quality: Option<u8>,
}

/// How one read lies on the construct.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct ReadResult {
    pub name: String,
    pub read_length: usize,
    /// The stretch kept after trimming (`start..end` of the read); `None` if not trimmed.
    pub kept: Option<(usize, usize)>,
    /// Why the read couldn't be placed, if it couldn't.
    pub error: Option<String>,
    pub strand: Option<Strand>,
    /// Covered stretch of the construct (0-based, half-open; wraps on circles).
    pub start: usize,
    pub end: usize,
    pub length: usize,
    pub identity: f64,
    pub differences: Vec<ReadDifference>,
    /// The read's bases as they pair with the construct's top strand, one per construct
    /// position from `start` (`-` where the read has a deletion), so the rows of several reads
    /// line up in one grid. Only built when the caller asks for rows.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub row: Option<String>,
    /// The quality of each base of `row`, the same length (`None` for a gap, or a read that
    /// carries no qualities).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub row_quality: Option<Vec<Option<u8>>>,
    /// Offsets in `row` that the read has extra bases before. Insertions get no column of
    /// their own, so one grid serves every read; the bases themselves are in `differences`.
    #[serde(skip_serializing_if = "Vec::is_empty")]
    #[serde(default)]
    pub row_insertions: Vec<usize>,
    /// The base of the untrimmed read (its own orientation) behind each column of `row`, the
    /// same length (`None` for a gap), so a trace can be drawn under the columns it belongs to.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub row_positions: Option<Vec<Option<usize>>>,
    /// Construct position → base of the untrimmed read (its own orientation) aligned to it.
    #[serde(skip)]
    pub positions: BTreeMap<usize, usize>,
}

/// One difference pooled over the reads.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct Finding {
    #[serde(flatten)]
    pub difference: Difference,
    /// Reads showing it: (read index, position in the read, quality).
    pub supported_by: Vec<(usize, Option<usize>, Option<u8>)>,
    /// Reads covering its position, with the base of each read aligned there (for a read
    /// showing an insertion, its first inserted base; for a deletion, the base before it).
    pub covered_by: Vec<(usize, Option<usize>)>,
    /// "confirmed", "one_read" or "disputed".
    pub verdict: &'static str,
}

/// The whole check.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct Check {
    pub reads: Vec<ReadResult>,
    /// Stretches no read covers (0-based, half-open; one may wrap the origin).
    pub uncovered: Vec<(usize, usize)>,
    /// Bases covered by at least one read, and by at least two.
    pub covered: usize,
    pub covered_twice: usize,
    pub findings: Vec<Finding>,
}

/// The position of a difference in the construct.
fn position(d: &Difference) -> usize {
    match d {
        Difference::Substitution { position, .. } | Difference::Insertion { position, .. } => {
            *position
        }
        Difference::Deletion { start, .. } => *start,
    }
}

/// Read positions under the alignment: construct position → aligned read index for paired
/// columns, and → index of the first inserted base for insertions before it.
fn column_maps(
    a: &MoleculeAlignment,
    len: usize,
) -> (BTreeMap<usize, usize>, BTreeMap<usize, usize>) {
    let mut paired = BTreeMap::new();
    let mut inserted = BTreeMap::new();
    let mut r = a.start;
    let mut q = match a.strand {
        Strand::Forward => a.query_start,
        Strand::Reverse => a.query_length.saturating_sub(a.query_end),
    };
    let wrap = |p: usize| if len == 0 { 0 } else { p % len };
    for (rc, qc) in a.reference_row.chars().zip(a.query_row.chars()) {
        match (rc == '-', qc == '-') {
            (false, false) => {
                paired.insert(wrap(r), q);
                r += 1;
                q += 1;
            }
            (true, _) => {
                inserted.entry(wrap(r)).or_insert(q);
                q += 1;
            }
            (false, true) => r += 1,
        }
    }
    (paired, inserted)
}

/// Whether construct position `p` lies in `start..end` (wrapping on circles; `start == end`
/// with `length == len` is the whole circle).
fn covers(start: usize, end: usize, length: usize, len: usize, p: usize) -> bool {
    if length >= len {
        true
    } else if start <= end {
        (start..end).contains(&p)
    } else {
        p >= start || p < end
    }
}

/// Check `molecule` with `reads`. With `rows`, every placed read also gets its bases in the
/// construct's orientation ([`ReadResult::row`]), so they can be shown one under the other.
pub fn check(molecule: &Molecule, reads: &[Read], trim: bool, rows: bool) -> Result<Check, String> {
    let len = molecule.len();
    if len == 0 || !molecule.seq().alphabet().is_nucleotide() {
        return Err("the construct must be a DNA sequence".into());
    }
    if len > MAX_CONSTRUCT {
        return Err(format!(
            "constructs up to {MAX_CONSTRUCT} bp can be checked"
        ));
    }
    if reads.is_empty() || reads.len() > MAX_READS {
        return Err(format!("give 1 to {MAX_READS} reads"));
    }
    let circular = molecule.topology() == Topology::Circular;
    let mut results = Vec::with_capacity(reads.len());
    for read in reads {
        results.push(place(molecule, read, trim, len, rows));
    }
    // Coverage per base.
    let mut depth = vec![0u16; len];
    for r in results.iter().filter(|r| r.error.is_none()) {
        let mut p = r.start;
        for _ in 0..r.length.min(len) {
            if let Some(d) = depth.get_mut(p) {
                *d = d.saturating_add(1);
            }
            p = if p + 1 == len { 0 } else { p + 1 };
        }
    }
    let covered = depth.iter().filter(|&&d| d > 0).count();
    let covered_twice = depth.iter().filter(|&&d| d > 1).count();
    let mut uncovered: Vec<(usize, usize)> = Vec::new();
    let mut i = 0;
    while i < len {
        if depth.get(i) == Some(&0) {
            let s = i;
            while depth.get(i) == Some(&0) {
                i += 1;
            }
            uncovered.push((s, i));
        } else {
            i += 1;
        }
    }
    // On a circle, gaps at both ends are one gap across the origin.
    if circular && uncovered.len() > 1 {
        let first = uncovered.first().copied();
        let last = uncovered.last().copied();
        if let (Some((0, fe)), Some((ls, le))) = (first, last)
            && le == len
        {
            uncovered.remove(0);
            uncovered.pop();
            uncovered.push((ls, fe));
        }
    }
    // Pool the differences.
    let mut pooled: BTreeMap<(usize, String), Finding> = BTreeMap::new();
    for (k, r) in results.iter().enumerate() {
        for d in &r.differences {
            let key = (position(&d.difference), format!("{:?}", d.difference));
            pooled
                .entry(key)
                .or_insert_with(|| Finding {
                    difference: d.difference.clone(),
                    supported_by: Vec::new(),
                    covered_by: Vec::new(),
                    verdict: "",
                })
                .supported_by
                .push((k, d.read_position, d.quality));
        }
    }
    let mut findings: Vec<Finding> = pooled.into_values().collect();
    for f in &mut findings {
        let p = position(&f.difference);
        let mut covered: Vec<(usize, Option<usize>)> = results
            .iter()
            .enumerate()
            .filter(|(_, r)| r.error.is_none() && covers(r.start, r.end, r.length, len, p))
            .map(|(k, r)| (k, r.positions.get(&p).copied()))
            .collect();
        // Reads showing it: their own position for it (an insertion just past a read's end
        // still counts as covered).
        for (k, pos, _) in &f.supported_by {
            match covered.iter_mut().find(|(c, _)| c == k) {
                Some(entry) => entry.1 = *pos,
                None => covered.push((*k, *pos)),
            }
        }
        covered.sort_unstable();
        f.covered_by = covered;
        f.verdict = if f.supported_by.len() < f.covered_by.len() {
            "disputed"
        } else if f.supported_by.len() >= 2 {
            "confirmed"
        } else {
            "one_read"
        };
    }
    Ok(Check {
        reads: results,
        uncovered,
        covered,
        covered_twice,
        findings,
    })
}

/// Trim and align one read.
fn place(molecule: &Molecule, read: &Read, trim: bool, len: usize, rows: bool) -> ReadResult {
    let mut out = ReadResult {
        name: read.name.clone(),
        read_length: read.residues.len(),
        kept: None,
        error: None,
        strand: None,
        start: 0,
        end: 0,
        length: 0,
        identity: 0.0,
        differences: Vec::new(),
        row: None,
        row_quality: None,
        row_insertions: Vec::new(),
        row_positions: None,
        positions: BTreeMap::new(),
    };
    let keep = match (read.keep, &read.quality, trim) {
        (Some((a, b)), _, _) => a.min(read.residues.len())..b.min(read.residues.len()),
        (None, quality, trim) => match (quality, trim) {
            (Some(q), true) => match quality_trim(q, DEFAULT_CUTOFF) {
                Some(k) => k,
                None => {
                    out.error = Some("no stretch of the read is of good enough quality".into());
                    return out;
                }
            },
            _ => 0..read.residues.len(),
        },
    };
    if read.keep.is_some() || (trim && read.quality.is_some()) {
        out.kept = Some((keep.start, keep.end));
    }
    let query = read.residues.get(keep.clone()).unwrap_or(&[]);
    if query.len() < 20 {
        out.error = Some("fewer than 20 bases to align".into());
        return out;
    }
    let a = match align_to_molecule(molecule, query, Mode::Fit, &Scoring::default()) {
        Ok(a) => a,
        Err(e) => {
            out.error = Some(e);
            return out;
        }
    };
    // A read that doesn't belong to this construct aligns, but badly.
    if a.identity < 0.8 {
        out.error = Some(format!(
            "doesn't match the construct ({:.0}% identity)",
            a.identity * 100.0
        ));
        return out;
    }
    let (paired, inserted) = column_maps(&a, len);
    // Aligned read index → index in the untrimmed read, in its own orientation.
    let qlen = query.len();
    let original = |qa: usize| -> Option<usize> {
        let t = match a.strand {
            Strand::Forward => qa,
            Strand::Reverse => qlen.checked_sub(1)?.checked_sub(qa)?,
        };
        keep.start.checked_add(t)
    };
    let quality_at =
        |i: Option<usize>| -> Option<u8> { read.quality.as_ref()?.get(i?).copied().flatten() };
    for d in &a.differences {
        let (read_position, quality) = match d {
            Difference::Substitution { position, .. } => {
                let i = paired.get(position).and_then(|&q| original(q));
                (i, quality_at(i))
            }
            Difference::Insertion { position, .. } => {
                let i = inserted.get(position).and_then(|&q| original(q));
                (i, quality_at(i))
            }
            Difference::Deletion { start, length, .. } => {
                let before = start
                    .checked_sub(1)
                    .or_else(|| len.checked_sub(1))
                    .and_then(|p| paired.get(&p))
                    .and_then(|&q| original(q));
                let after = paired
                    .get(&((start + length) % len.max(1)))
                    .and_then(|&q| original(q));
                let q = match (quality_at(before), quality_at(after)) {
                    (Some(a), Some(b)) => Some(a.min(b)),
                    (a, b) => a.or(b),
                };
                (before.or(after), q)
            }
        };
        out.differences.push(ReadDifference {
            difference: d.clone(),
            read_position,
            quality,
        });
    }
    if rows {
        // One column per construct position the read covers, in the construct's orientation.
        // A column the read has no base for is a gap; extra bases of the read get no column,
        // only a mark on the column that follows them.
        let mut row = String::with_capacity(a.length);
        let mut row_quality = Vec::with_capacity(a.length);
        let mut insertions = Vec::new();
        let mut row_positions = Vec::with_capacity(a.length);
        let mut q = match a.strand {
            Strand::Forward => a.query_start,
            Strand::Reverse => a.query_length.saturating_sub(a.query_end),
        };
        for (rc, qc) in a.reference_row.chars().zip(a.query_row.chars()) {
            if rc == '-' {
                // An extra base of the read, before the column that comes next.
                if insertions.last() != Some(&row.len()) {
                    insertions.push(row.len());
                }
                q += 1;
                continue;
            }
            if qc == '-' {
                row.push('-');
                row_quality.push(None);
                row_positions.push(None);
            } else {
                row.push(qc.to_ascii_uppercase());
                row_quality.push(quality_at(original(q)));
                row_positions.push(original(q));
                q += 1;
            }
        }
        out.row_quality = read.quality.as_ref().map(|_| row_quality);
        out.row = Some(row);
        out.row_insertions = insertions;
        out.row_positions = Some(row_positions);
    }
    out.positions = paired
        .iter()
        .filter_map(|(&p, &q)| Some((p, original(q)?)))
        .collect();
    out.strand = Some(a.strand);
    out.start = a.start;
    out.end = a.end;
    out.length = a.length;
    out.identity = a.identity;
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use snappiesnap_doc::{Alphabet, Seq};

    /// A deterministic, low-repeat test sequence.
    fn construct(n: usize) -> String {
        let mut x: u32 = 12345;
        (0..n)
            .map(|_| {
                x = x.wrapping_mul(1_103_515_245).wrapping_add(12345);
                ['A', 'C', 'G', 'T'][((x >> 16) & 3) as usize]
            })
            .collect()
    }

    fn revcomp(s: &str) -> String {
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

    fn read(name: &str, residues: &str, q: u8) -> Read {
        Read {
            name: name.into(),
            residues: residues.as_bytes().to_vec(),
            quality: Some(vec![Some(q); residues.len()]),
            keep: None,
        }
    }

    #[test]
    fn pooled_check() {
        let c = construct(1000);
        let m = Molecule::new(
            "c",
            Seq::new(Alphabet::Dna, c.as_str()).unwrap(),
            Topology::Circular,
        )
        .unwrap();
        // Read 1: forward 100..500 with a substitution at 300 (reported in two reads).
        let mut r1: Vec<u8> = c.as_bytes()[100..500].to_vec();
        let sub = if r1[200] == b'A' { b'C' } else { b'A' };
        r1[200] = sub;
        // Read 2: reverse complement of 250..700, also with the substitution at 300, and with
        // a poor call at 600 (quality 8) that read 3 doesn't share.
        let mut fwd2: Vec<u8> = c.as_bytes()[250..700].to_vec();
        fwd2[50] = sub;
        let odd = if fwd2[350] == b'G' { b'T' } else { b'G' };
        fwd2[350] = odd;
        let r2 = revcomp(std::str::from_utf8(&fwd2).unwrap());
        let mut q2 = vec![Some(40u8); r2.len()];
        q2[r2.len() - 1 - 350] = Some(8);
        // Read 3: forward 550..980, exact; plus 30 poor bases at each end to trim.
        let r3 = format!("{}{}{}", "N".repeat(30), &c[550..980], "N".repeat(30));
        let mut q3 = vec![Some(40u8); r3.len()];
        for q in q3.iter_mut().take(30) {
            *q = Some(2);
        }
        for q in q3.iter_mut().rev().take(30) {
            *q = Some(2);
        }
        let reads = vec![
            read("r1", std::str::from_utf8(&r1).unwrap(), 35),
            Read {
                name: "r2".into(),
                residues: r2.into_bytes(),
                quality: Some(q2),
                keep: None,
            },
            Read {
                name: "r3".into(),
                residues: r3.into_bytes(),
                quality: Some(q3),
                keep: None,
            },
        ];
        let check = check(&m, &reads, true, false).unwrap();
        let [a, b, c3] = &check.reads[..] else {
            panic!()
        };
        assert_eq!(
            (a.strand, a.start, a.end),
            (Some(Strand::Forward), 100, 500)
        );
        assert_eq!(
            (b.strand, b.start, b.end),
            (Some(Strand::Reverse), 250, 700)
        );
        assert_eq!((c3.start, c3.end, c3.kept), (550, 980, Some((30, 460))));
        assert!(c3.differences.is_empty());
        // Covered 100..980; the rest wraps the origin.
        assert_eq!(check.covered, 880);
        assert_eq!(check.uncovered, vec![(980, 100)]);
        assert_eq!(check.covered_twice, (500 - 250) + (700 - 550));
        let at = |p: usize| {
            check
                .findings
                .iter()
                .find(|f| position(&f.difference) == p)
                .unwrap()
        };
        let s300 = at(300);
        assert_eq!(s300.verdict, "confirmed");
        assert_eq!(s300.supported_by.len(), 2);
        // In read 1 it is base 200, quality 35; in read 2 (reverse) base 449 - 50 = 399.
        assert_eq!(s300.supported_by[0], (0, Some(200), Some(35)));
        assert_eq!(s300.supported_by[1], (1, Some(449 - 50), Some(40)));
        let s600 = at(600);
        assert_eq!(s600.verdict, "disputed");
        assert_eq!(s600.supported_by, vec![(1, Some(449 - 350), Some(8))]);
        // Read 3 (forward from 550, trimmed by 30) has construct base 600 at read base 80.
        assert_eq!(
            s600.covered_by,
            vec![(1, Some(449 - 350)), (2, Some(30 + 50))]
        );
    }

    #[test]
    fn insertions_deletions_and_bad_reads() {
        let c = construct(800);
        let m = Molecule::new(
            "c",
            Seq::new(Alphabet::Dna, c.as_str()).unwrap(),
            Topology::Linear,
        )
        .unwrap();
        // 3 extra bases after construct position 400, and 2 missing at 600..602.
        let r = format!("{}GGG{}{}", &c[200..400], &c[400..600], &c[602..750]);
        let check = check(&m, &[read("indels", &r, 30)], true, false).unwrap();
        let r0 = &check.reads[0];
        assert_eq!(r0.error, None);
        assert_eq!(r0.row, None, "rows are built only when asked for");
        // The same read with its row: one column per construct position it covers, the extra
        // bases marked rather than given columns of their own, the missing ones as gaps.
        let rows = super::check(&m, &[read("indels", &r, 30)], true, true).unwrap();
        let r1 = &rows.reads[0];
        let row = r1.row.as_deref().unwrap_or_default();
        assert_eq!(row.len(), r1.length, "a column for every construct base");
        assert_eq!(row.matches('-').count(), 2, "the two missing bases: {row}");
        assert_eq!(r1.row_insertions, vec![400 - r1.start]);
        assert_eq!(r1.row_quality.as_ref().map(Vec::len), Some(row.len()));
        // Every column knows the read base behind it: one step per column, three more across
        // the extra bases, none for the missing ones.
        let at = r1.row_positions.clone().unwrap_or_default();
        assert_eq!(at.len(), row.len());
        let col = |p: usize| at.get(p - r1.start).copied().flatten();
        assert_eq!(col(r1.start), Some(r1.start - 200));
        assert_eq!(col(399), Some(199));
        assert_eq!(col(400), Some(203), "after the three extra bases");
        let gaps: Vec<usize> = (0..at.len()).filter(|&i| at[i].is_none()).collect();
        let dashes: Vec<usize> = row
            .bytes()
            .enumerate()
            .filter(|(_, b)| *b == b'-')
            .map(|(i, _)| i)
            .collect();
        assert_eq!(gaps, dashes, "no read base behind a gap");
        // A read of the other strand counts down along the construct.
        let rc: String = r
            .bytes()
            .rev()
            .map(|b| match b {
                b'A' => 'T',
                b'C' => 'G',
                b'G' => 'C',
                _ => 'A',
            })
            .collect();
        let back = super::check(&m, &[read("back", &rc, 30)], false, true).unwrap();
        let b0 = &back.reads[0];
        assert_eq!(b0.strand, Some(Strand::Reverse));
        let at = b0.row_positions.clone().unwrap_or_default();
        assert_eq!(at.first().copied().flatten(), Some(rc.len() - 1));
        assert_eq!(at.last().copied().flatten(), Some(0));
        // Column by column, the row is the construct except where the read differs.
        let bases: Vec<u8> = c.bytes().collect();
        let same = row
            .bytes()
            .enumerate()
            .filter(|(i, b)| bases.get(r1.start + i) == Some(b))
            .count();
        assert_eq!(same, row.len() - 2, "only the gaps differ: {row}");
        let ins = r0
            .differences
            .iter()
            .find(|d| matches!(d.difference, Difference::Insertion { .. }))
            .unwrap();
        assert_eq!(position(&ins.difference), 400);
        assert_eq!(ins.read_position, Some(200));
        let del = r0
            .differences
            .iter()
            .find(|d| matches!(d.difference, Difference::Deletion { .. }))
            .unwrap();
        // Deleting 600..602 can be placed elsewhere in a repeat (here CCTCC: 599..601 gives the
        // same sequence); any placement that gives the same sequence is right.
        let Difference::Deletion {
            start, length: 2, ..
        } = del.difference
        else {
            panic!("{del:?}");
        };
        let without = |s: usize| format!("{}{}", &c[..s], &c[s + 2..]);
        assert_eq!(without(start), without(600), "{del:?}");
        // The read base before the gap: construct start - 1, 200 bases into the read plus 3.
        assert_eq!(del.read_position, Some(start - 1 - 200 + 3));
        assert_eq!(
            check
                .findings
                .iter()
                .filter(|f| f.verdict == "one_read")
                .count(),
            2
        );
        // Unrelated or empty reads are reported, not placed.
        let junk = construct(3000)[2000..2400].to_string();
        let check2 = super::check(
            &m,
            &[read("junk", &junk, 30), read("short", "ACGT", 30)],
            true,
            false,
        )
        .unwrap();
        assert!(
            check2.reads[0]
                .error
                .as_deref()
                .unwrap()
                .contains("identity")
        );
        assert!(
            check2.reads[1]
                .error
                .as_deref()
                .unwrap()
                .contains("20 bases")
        );
        assert_eq!(check2.covered, 0);
        assert_eq!(check2.uncovered, vec![(0, 800)]);
        let all_poor = Read {
            name: "poor".into(),
            residues: c.as_bytes()[0..100].to_vec(),
            quality: Some(vec![Some(1); 100]),
            keep: None,
        };
        let check3 = super::check(&m, &[all_poor], true, false).unwrap();
        assert!(
            check3.reads[0]
                .error
                .as_deref()
                .unwrap()
                .contains("quality")
        );
        assert!(super::check(&m, &[], true, false).is_err());
        // A stretch chosen beforehand (as stored with matched reads) replaces the trimming.
        let chosen = Read {
            name: "chosen".into(),
            residues: c.as_bytes()[0..100].to_vec(),
            quality: Some(vec![Some(1); 100]),
            keep: Some((10, 90)),
        };
        let check4 = super::check(&m, &[chosen], true, false).unwrap();
        assert_eq!(check4.reads[0].error, None);
        assert_eq!(check4.reads[0].kept, Some((10, 90)));
        assert_eq!((check4.reads[0].start, check4.reads[0].end), (10, 90));
    }
}
