//! SCF, the Standard Chromatogram Format (`.scf`), as published by Dear S. and Staden R., "A
//! standard file format for data from DNA sequencing instruments", DNA Sequence 3:107–110
//! (1992), and the Staden package's public description of versions 2 and 3.
//!
//! A 128-byte header of big-endian 32-bit fields (magic `.scf`, sample count and offset, base
//! count, clip points and offset, comment size and offset, version text, sample size in bytes,
//! code set, private data size and offset) is followed by three sections:
//! - **samples**: version 3 stores each channel in turn (A, C, G, T), every channel encoded as
//!   second differences; earlier versions store A, C, G, T interleaved per sample point, raw;
//! - **bases**: version 3 stores all peak indices (32-bit), then the A, C, G and T
//!   probabilities (one byte each), then the base letters, then three spare bytes, each as a
//!   column; earlier versions store one 12-byte record per base in that order;
//! - **comments**: `name=value` lines.
//!
//! The quality of a base is the probability value of its own channel (version 3 stores Phred
//! values there); `N` and other ambiguity codes take the highest of the four.

use crate::trace::Chromatogram;

/// Why an SCF file couldn't be read.
#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum ScfError {
    #[error("not an SCF file: it must start with \".scf\"")]
    NotScf,

    #[error("the SCF file is damaged: {0}")]
    Malformed(String),

    #[error("the file has no base calls")]
    NoBaseCalls,
}

fn be_u32(b: &[u8], at: usize) -> Option<u32> {
    Some(u32::from_be_bytes(
        b.get(at..at.checked_add(4)?)?.try_into().ok()?,
    ))
}

fn usize_at(b: &[u8], at: usize) -> Result<usize, ScfError> {
    be_u32(b, at)
        .and_then(|v| usize::try_from(v).ok())
        .ok_or_else(|| ScfError::Malformed("the header is cut short".into()))
}

/// The bytes `len` long at `offset`, or an error naming `what`.
fn section<'a>(
    b: &'a [u8],
    offset: usize,
    len: Option<usize>,
    what: &str,
) -> Result<&'a [u8], ScfError> {
    len.and_then(|n| offset.checked_add(n))
        .and_then(|end| b.get(offset..end))
        .ok_or_else(|| ScfError::Malformed(format!("the {what} run past the end of the file")))
}

/// Read the chromatogram of an SCF file.
pub fn read_chromatogram(bytes: &[u8]) -> Result<Chromatogram, ScfError> {
    if bytes.get(..4) != Some(b".scf") {
        return Err(ScfError::NotScf);
    }
    let samples = usize_at(bytes, 4)?;
    let samples_offset = usize_at(bytes, 8)?;
    let bases = usize_at(bytes, 12)?;
    let bases_offset = usize_at(bytes, 24)?;
    let comments_size = usize_at(bytes, 28)?;
    let comments_offset = usize_at(bytes, 32)?;
    let version = bytes
        .get(36..40)
        .ok_or_else(|| ScfError::Malformed("the header is cut short".into()))?;
    let v3 = version.first().is_some_and(|&b| b >= b'3');
    let sample_size = usize_at(bytes, 40)?;
    if !matches!(sample_size, 1 | 2) {
        return Err(ScfError::Malformed(format!(
            "sample size must be 1 or 2 bytes, not {sample_size}"
        )));
    }
    if bases == 0 {
        return Err(ScfError::NoBaseCalls);
    }
    // Lengths come from the header: `section` checks them against the file before any copy,
    // so allocations are bounded by the input size.
    let raw = section(
        bytes,
        samples_offset,
        samples
            .checked_mul(4)
            .and_then(|n| n.checked_mul(sample_size)),
        "samples",
    )?;
    let value = |c: &[u8]| match c {
        [a] => u16::from(*a),
        [a, b] => u16::from_be_bytes([*a, *b]),
        _ => 0,
    };
    let mut channels: [Vec<u16>; 4] = Default::default();
    if v3 {
        let per_channel = samples.saturating_mul(sample_size);
        for (i, c) in channels.iter_mut().enumerate() {
            let start = i.saturating_mul(per_channel);
            let data = raw
                .get(start..start.saturating_add(per_channel))
                .unwrap_or(&[]);
            *c = undelta(
                data.chunks_exact(sample_size).map(value).collect(),
                sample_size,
            );
        }
    } else {
        for point in raw.chunks_exact(sample_size.saturating_mul(4)) {
            for (c, v) in channels.iter_mut().zip(point.chunks_exact(sample_size)) {
                c.push(value(v));
            }
        }
    }
    let raw = section(bytes, bases_offset, bases.checked_mul(12), "base calls")?;
    let mut letters = Vec::with_capacity(bases);
    let mut peaks = Vec::with_capacity(bases);
    let mut quality = Vec::with_capacity(bases);
    for i in 0..bases {
        let (peak, probs, letter) = if v3 {
            let peak = be_u32(raw, i.saturating_mul(4));
            let prob = |c: usize| {
                raw.get(bases.saturating_mul(4 + c).saturating_add(i))
                    .copied()
                    .unwrap_or(0)
            };
            let letter = raw.get(bases.saturating_mul(8).saturating_add(i)).copied();
            (peak, [prob(0), prob(1), prob(2), prob(3)], letter)
        } else {
            let at = i.saturating_mul(12);
            let rec = raw.get(at..at.saturating_add(12)).unwrap_or(&[]);
            let p = |k: usize| rec.get(k).copied().unwrap_or(0);
            (
                be_u32(rec, 0),
                [p(4), p(5), p(6), p(7)],
                rec.get(8).copied(),
            )
        };
        let (Some(peak), Some(letter)) = (peak, letter) else {
            return Err(ScfError::Malformed("a base call is cut short".into()));
        };
        let q = match letter.to_ascii_uppercase() {
            b'A' => probs[0],
            b'C' => probs[1],
            b'G' => probs[2],
            b'T' => probs[3],
            _ => probs.iter().copied().max().unwrap_or(0),
        };
        letters.push(letter);
        peaks.push(peak);
        quality.push(q);
    }
    let comments = section(bytes, comments_offset, Some(comments_size), "comments")
        .map(|c| String::from_utf8_lossy(c).into_owned())
        .unwrap_or_default();
    let sample_name = comments.lines().find_map(|l| {
        let (k, v) = l.split_once('=')?;
        let v = v.trim_end_matches('\0').trim();
        (k.trim().eq_ignore_ascii_case("NAME") && !v.is_empty()).then(|| v.to_string())
    });
    let c = Chromatogram {
        sample_name,
        bases: letters,
        peaks,
        quality: Some(quality),
        channels,
    };
    c.check().map_err(ScfError::Malformed)?;
    Ok(c)
}

/// Undo the version 3 second-difference encoding: two running sums, wrapping at the sample
/// width as the encoder's arithmetic does.
fn undelta(mut v: Vec<u16>, sample_size: usize) -> Vec<u16> {
    let mask = if sample_size == 1 { 0xFF } else { 0xFFFF };
    for _ in 0..2 {
        let mut sum: u16 = 0;
        for x in &mut v {
            sum = sum.wrapping_add(*x) & mask;
            *x = sum;
        }
    }
    v
}

/// The version 3 second-difference encoding of one channel (16-bit samples).
fn delta(mut v: Vec<u16>) -> Vec<u16> {
    for _ in 0..2 {
        let mut prev: u16 = 0;
        for x in &mut v {
            let cur = *x;
            *x = cur.wrapping_sub(prev);
            prev = cur;
        }
    }
    v
}

/// Write a chromatogram as an SCF version 3.00 file with 16-bit samples. Used for synthetic
/// test traces; files read from disk are saved from their original bytes.
pub fn write_chromatogram(c: &Chromatogram) -> Result<Vec<u8>, String> {
    c.check()?;
    let too_big = || "the trace is too large for an SCF file".to_string();
    let u32_of = |n: usize| u32::try_from(n).map_err(|_| too_big());
    let samples = c.samples();
    let n = c.bases.len();
    let comments = match &c.sample_name {
        Some(name) => format!("NAME={name}\n").into_bytes(),
        None => Vec::new(),
    };
    let samples_offset = 128usize;
    let bases_offset = samples
        .checked_mul(8)
        .and_then(|s| s.checked_add(samples_offset))
        .ok_or_else(too_big)?;
    let comments_offset = n
        .checked_mul(12)
        .and_then(|b| b.checked_add(bases_offset))
        .ok_or_else(too_big)?;
    let mut out = Vec::new();
    out.extend_from_slice(b".scf");
    for v in [
        u32_of(samples)?,
        u32_of(samples_offset)?,
        u32_of(n)?,
        0,
        0,
        u32_of(bases_offset)?,
        u32_of(comments.len())?,
        u32_of(comments_offset)?,
    ] {
        out.extend_from_slice(&v.to_be_bytes());
    }
    out.extend_from_slice(b"3.00");
    for v in [2u32, 0, 0, 0] {
        out.extend_from_slice(&v.to_be_bytes());
    }
    out.resize(128, 0);
    for ch in &c.channels {
        for v in delta(ch.clone()) {
            out.extend_from_slice(&v.to_be_bytes());
        }
    }
    for p in &c.peaks {
        out.extend_from_slice(&p.to_be_bytes());
    }
    let quality = c.quality.clone().unwrap_or_else(|| vec![0; n]);
    for channel in *b"ACGT" {
        for (b, q) in c.bases.iter().zip(&quality) {
            let own = b.to_ascii_uppercase() == channel
                || !matches!(b.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T');
            out.push(if own { *q } else { 0 });
        }
    }
    out.extend_from_slice(&c.bases);
    out.resize(out.len().saturating_add(n.saturating_mul(3)), 0);
    out.extend_from_slice(&comments);
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample() -> Chromatogram {
        let mut channels: [Vec<u16>; 4] = Default::default();
        for (i, c) in channels.iter_mut().enumerate() {
            // Ramps that cross 0 and 65535 when differenced, to exercise wrapping.
            *c = (0..50u16)
                .map(|x| (x * 1300).wrapping_add(i as u16 * 7))
                .collect();
        }
        Chromatogram {
            sample_name: Some("clone 7 fwd".into()),
            bases: b"ACGTN".to_vec(),
            peaks: vec![3, 13, 23, 33, 43],
            quality: Some(vec![10, 20, 30, 40, 5]),
            channels,
        }
    }

    #[test]
    fn version3_round_trip() {
        let c = sample();
        let bytes = write_chromatogram(&c).unwrap();
        assert_eq!(read_chromatogram(&bytes).unwrap(), c);
    }

    #[test]
    fn second_differences_by_hand() {
        // 10, 12, 15: first differences 10, 2, 3; second differences 10, -8, 1.
        assert_eq!(delta(vec![10, 12, 15]), vec![10, 0u16.wrapping_sub(8), 1]);
        assert_eq!(
            undelta(vec![10, 0u16.wrapping_sub(8), 1], 2),
            vec![10, 12, 15]
        );
        // 8-bit samples wrap at 256: 250, 250 + 10 - 256 = 4; then 250, 250 + 4.
        assert_eq!(undelta(vec![250, 10], 1), vec![250, 254]);
    }

    #[test]
    fn version2_layout() {
        // Two samples interleaved A C G T, one base record of 12 bytes.
        let mut b = Vec::new();
        b.extend_from_slice(b".scf");
        for v in [2u32, 128, 1, 0, 0, 128 + 8, 0, 0] {
            b.extend_from_slice(&v.to_be_bytes());
        }
        b.extend_from_slice(b"2.00");
        b.extend_from_slice(&1u32.to_be_bytes());
        b.resize(128, 0);
        b.extend_from_slice(&[1, 2, 3, 4, 5, 6, 7, 8]);
        b.extend_from_slice(&1u32.to_be_bytes());
        b.extend_from_slice(&[0, 0, 33, 0, b'G', 0, 0, 0]);
        let c = read_chromatogram(&b).unwrap();
        assert_eq!(c.channels, [vec![1, 5], vec![2, 6], vec![3, 7], vec![4, 8]]);
        assert_eq!(c.bases, b"G");
        assert_eq!(c.peaks, vec![1]);
        assert_eq!(c.quality, Some(vec![33]));
        assert_eq!(c.sample_name, None);
    }

    #[test]
    fn hostile_files_are_errors() {
        assert_eq!(read_chromatogram(b"ABIF"), Err(ScfError::NotScf));
        let good = write_chromatogram(&sample()).unwrap();
        for cut in [4, 40, 127, 128, 300, good.len() - 20] {
            assert!(read_chromatogram(&good[..cut]).is_err(), "cut at {cut}");
        }
        // Absurd sample count.
        let mut big = good.clone();
        big[4..8].copy_from_slice(&u32::MAX.to_be_bytes());
        assert!(read_chromatogram(&big).is_err());
        // Sample size 3.
        let mut odd = good;
        odd[40..44].copy_from_slice(&3u32.to_be_bytes());
        assert!(read_chromatogram(&odd).is_err());
    }
}
