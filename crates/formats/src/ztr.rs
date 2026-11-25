//! ZTR, the compressed sequencing trace format of the Staden package, as published in "ZTR
//! SPEC v1.2" (Bonfield J.K. and Staden R., "ZTR: a new format for DNA sequence trace data",
//! Bioinformatics 18:3–10, 2002; staden.sourceforge.net/ztr.html). SnapGene `.dna` files keep
//! the reads matched to a construct in this format (see [`crate::dna`]).
//!
//! A file is an 8-byte magic (`AE 5A 54 52 0D 0A 1A 0A`), a major and a minor version byte,
//! then chunks: a 4-byte type, a 4-byte big-endian metadata length and the metadata, a 4-byte
//! big-endian data length and the data. A chunk's data starts with a format byte; anything
//! other than 0 (raw) is a filter whose output is again data starting with a format byte, so
//! filters stack. The spec leaves some details open; files written by the reference application
//! settle them as follows (checked on six traces: every Q30 base's own channel is the highest
//! at its peak):
//! - the 4-byte uncompressed lengths of RLE (1) and zlib (2) are **little-endian**;
//! - the delta filters (64, 65, 66) and the word-size reductions (70, 71) cover the whole inner
//!   buffer, including its own format byte and padding;
//! - the follow filter (72) stores `follow[previous] - value`.
//!
//! Chunks read: `SMP4` (or four `SAMP`) signal, `BASE` calls, `BPOS` peaks, `CNF4` confidences,
//! `TEXT` (the `NAME` entry is the sample name). Others are skipped. The integer Chebyshev
//! filter (74) has no published formula and is refused.

use crate::trace::Chromatogram;

/// Why a ZTR trace couldn't be read.
#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum ZtrError {
    #[error("not a ZTR trace")]
    NotZtr,

    #[error("ZTR version {0} isn't supported (version 1 is)")]
    Version(u8),

    #[error("the ZTR trace is damaged: {0}")]
    Malformed(String),

    #[error("the ZTR trace uses a compression SnappieSnap doesn't know (format {0})")]
    Unsupported(u8),
}

/// The ZTR magic.
pub const MAGIC: [u8; 8] = [0xAE, 0x5A, 0x54, 0x52, 0x0D, 0x0A, 0x1A, 0x0A];
/// Largest decompressed chunk accepted.
const MAX_CHUNK: usize = 64 << 20;
/// Most filters stacked on one chunk.
const MAX_LAYERS: usize = 16;

fn bad(m: impl Into<String>) -> ZtrError {
    ZtrError::Malformed(m.into())
}

fn le_len(d: &[u8]) -> Result<usize, ZtrError> {
    let b: [u8; 4] = d
        .get(1..5)
        .and_then(|b| b.try_into().ok())
        .ok_or_else(|| bad("a length is cut short"))?;
    let n = u32::from_le_bytes(b) as usize;
    if n > MAX_CHUNK {
        return Err(bad(format!("a chunk claims {n} bytes")));
    }
    Ok(n)
}

/// Run-length decoding: `guard, count, value` is `count` copies of `value`; `guard, 0` is the
/// guard byte itself.
fn unrle(d: &[u8]) -> Result<Vec<u8>, ZtrError> {
    let n = le_len(d)?;
    let guard = *d.get(5).ok_or_else(|| bad("RLE header cut short"))?;
    let mut out = Vec::with_capacity(n);
    let mut i = 6;
    while let Some(&c) = d.get(i) {
        if c == guard {
            let count = *d.get(i + 1).ok_or_else(|| bad("RLE run cut short"))?;
            if count == 0 {
                out.push(guard);
                i += 2;
            } else {
                let v = *d.get(i + 2).ok_or_else(|| bad("RLE run cut short"))?;
                out.extend(std::iter::repeat_n(v, usize::from(count)));
                i += 3;
            }
        } else {
            out.push(c);
            i += 1;
        }
        if out.len() > n {
            return Err(bad("RLE data longer than declared"));
        }
    }
    if out.len() != n {
        return Err(bad("RLE data shorter than declared"));
    }
    Ok(out)
}

fn unzlib(d: &[u8]) -> Result<Vec<u8>, ZtrError> {
    let n = le_len(d)?;
    let out = miniz_oxide::inflate::decompress_to_vec_zlib_with_limit(d.get(5..).unwrap_or(&[]), n)
        .map_err(|_| bad("zlib data can't be decompressed"))?;
    if out.len() != n {
        return Err(bad("zlib data has the wrong length"));
    }
    Ok(out)
}

/// Undo `level` rounds of differencing on `width`-byte big-endian words after a `header`-byte
/// filter header.
fn undelta(d: &[u8], width: usize, header: usize) -> Result<Vec<u8>, ZtrError> {
    let level = *d.get(1).ok_or_else(|| bad("delta header cut short"))?;
    if level > 8 {
        return Err(bad(format!("delta level {level}")));
    }
    let body = d
        .get(header..)
        .ok_or_else(|| bad("delta header cut short"))?;
    if !body.len().is_multiple_of(width) {
        return Err(bad("delta data isn't whole words"));
    }
    let mask: u64 = if width >= 8 {
        u64::MAX
    } else {
        (1u64 << (8 * width)) - 1
    };
    let mut vals: Vec<u64> = body
        .chunks_exact(width)
        .map(|w| w.iter().fold(0u64, |acc, &b| (acc << 8) | u64::from(b)))
        .collect();
    for _ in 0..level {
        let mut acc = 0u64;
        for v in &mut vals {
            acc = acc.wrapping_add(*v) & mask;
            *v = acc;
        }
    }
    let mut out = Vec::with_capacity(body.len());
    for v in vals {
        for k in (0..width).rev() {
            out.push(((v >> (8 * k)) & 0xFF) as u8);
        }
    }
    Ok(out)
}

/// Undo a reduction of `width`-byte signed words to single bytes: a byte from -127 to 127 is
/// the value, -128 is followed by the word itself.
fn unshrink(d: &[u8], width: usize) -> Result<Vec<u8>, ZtrError> {
    let mut out = Vec::new();
    let mut i = 1;
    while let Some(&b) = d.get(i) {
        if b == 0x80 {
            let w = d
                .get(i + 1..i + 1 + width)
                .ok_or_else(|| bad("escaped word cut short"))?;
            out.extend_from_slice(w);
            i += 1 + width;
        } else {
            let v = i64::from(b as i8);
            let bytes = v.to_be_bytes();
            out.extend_from_slice(bytes.get(8 - width..).unwrap_or(&[]));
            i += 1;
        }
        if out.len() > MAX_CHUNK {
            return Err(bad("expanded data too large"));
        }
    }
    Ok(out)
}

/// Undo the follow filter: a 256-byte table of predicted successors, then the first byte as
/// is and each later byte as `follow[previous] - value`.
fn unfollow(d: &[u8]) -> Result<Vec<u8>, ZtrError> {
    let table = d.get(1..257).ok_or_else(|| bad("follow table cut short"))?;
    let data = d.get(257..).unwrap_or(&[]);
    let mut out: Vec<u8> = Vec::with_capacity(data.len());
    for &b in data {
        let v = match out.last() {
            None => b,
            Some(&prev) => table
                .get(usize::from(prev))
                .copied()
                .unwrap_or(0)
                .wrapping_sub(b),
        };
        out.push(v);
    }
    Ok(out)
}

/// Undo every filter of a chunk's data, down to raw data (its format byte 0 included).
pub fn unfilter(data: &[u8]) -> Result<Vec<u8>, ZtrError> {
    let mut d = data.to_vec();
    for _ in 0..MAX_LAYERS {
        let format = *d.first().ok_or_else(|| bad("empty chunk data"))?;
        d = match format {
            0 => return Ok(d),
            1 => unrle(&d)?,
            2 => unzlib(&d)?,
            64 => undelta(&d, 1, 2)?,
            65 => undelta(&d, 2, 2)?,
            66 => undelta(&d, 4, 4)?,
            70 => unshrink(&d, 2)?,
            71 => unshrink(&d, 4)?,
            72 => unfollow(&d)?,
            other => return Err(ZtrError::Unsupported(other)),
        };
    }
    Err(bad("too many filters on one chunk"))
}

fn be_u32(d: &[u8], at: usize) -> Result<usize, ZtrError> {
    d.get(at..at + 4)
        .and_then(|b| b.try_into().ok())
        .map(|b: [u8; 4]| u32::from_be_bytes(b) as usize)
        .ok_or_else(|| bad("a chunk header is cut short"))
}

/// Big-endian `width`-byte words after `skip` bytes (the raw format byte and padding).
fn words(d: &[u8], skip: usize, width: usize) -> Result<Vec<u32>, ZtrError> {
    let body = d.get(skip..).ok_or_else(|| bad("raw chunk cut short"))?;
    Ok(body
        .chunks_exact(width)
        .map(|w| w.iter().fold(0u32, |acc, &b| (acc << 8) | u32::from(b)))
        .collect())
}

/// A chunk: type, metadata, data.
type Chunk<'a> = ([u8; 4], &'a [u8], &'a [u8]);

/// The chunks of a ZTR file, in file order.
fn chunks(bytes: &[u8]) -> Result<Vec<Chunk<'_>>, ZtrError> {
    if bytes.get(..8) != Some(&MAGIC[..]) {
        return Err(ZtrError::NotZtr);
    }
    let major = *bytes.get(8).ok_or_else(|| bad("the header is cut short"))?;
    if major != 1 {
        return Err(ZtrError::Version(major));
    }
    let mut out = Vec::new();
    let mut at = 10;
    while at < bytes.len() {
        let kind: [u8; 4] = bytes
            .get(at..at + 4)
            .and_then(|b| b.try_into().ok())
            .ok_or_else(|| bad("a chunk type is cut short"))?;
        let meta_len = be_u32(bytes, at + 4)?;
        let meta_start = at + 8;
        let meta = bytes
            .get(meta_start..meta_start.saturating_add(meta_len))
            .ok_or_else(|| bad("chunk metadata runs past the end"))?;
        let data_at = meta_start + meta_len;
        let data_len = be_u32(bytes, data_at)?;
        let data = bytes
            .get(data_at + 4..(data_at + 4).saturating_add(data_len))
            .ok_or_else(|| bad("chunk data runs past the end"))?;
        at = data_at + 4 + data_len;
        out.push((kind, meta, data));
    }
    Ok(out)
}

/// The base calls of a ZTR file (its `BASE` chunk), also when it holds no signal, as for a
/// plain sequence stored in the format.
pub fn read_bases(bytes: &[u8]) -> Result<Vec<u8>, ZtrError> {
    let base = chunks(bytes)?
        .into_iter()
        .find(|(kind, _, _)| kind == b"BASE")
        .ok_or_else(|| bad("there are no base calls"))?;
    Ok(unfilter(base.2)?.get(1..).unwrap_or(&[]).to_vec())
}

/// Read a ZTR trace.
pub fn read_chromatogram(bytes: &[u8]) -> Result<Chromatogram, ZtrError> {
    let mut smp4: Option<Vec<u8>> = None;
    let mut samp: [Option<Vec<u8>>; 4] = Default::default();
    let (mut base, mut bpos, mut cnf4, mut text) = (None, None, None, None);
    for (kind, meta, data) in chunks(bytes)? {
        match &kind {
            b"SMP4" => smp4 = Some(unfilter(data)?),
            b"SAMP" => {
                let lane = match meta.first() {
                    Some(b'A') => 0,
                    Some(b'C') => 1,
                    Some(b'G') => 2,
                    Some(b'T') => 3,
                    _ => continue,
                };
                if let Some(slot) = samp.get_mut(lane) {
                    *slot = Some(unfilter(data)?);
                }
            }
            b"BASE" => base = Some(unfilter(data)?),
            b"BPOS" => bpos = Some(unfilter(data)?),
            b"CNF4" => cnf4 = Some(unfilter(data)?),
            b"TEXT" => text = Some(unfilter(data)?),
            _ => {}
        }
    }
    // Signal: SMP4 is A, C, G, T one after the other, 16-bit, after the format byte and one
    // byte of padding; SAMP is one channel each, laid out the same.
    let channels: [Vec<u16>; 4] = match (smp4, samp) {
        (Some(d), _) => {
            let all = words(&d, 2, 2)?;
            let n = all.len() / 4;
            let lane = |k: usize| -> Vec<u16> {
                all.get(k * n..(k + 1) * n)
                    .unwrap_or(&[])
                    .iter()
                    .map(|&v| u16::try_from(v).unwrap_or(u16::MAX))
                    .collect()
            };
            [lane(0), lane(1), lane(2), lane(3)]
        }
        (None, [Some(a), Some(c), Some(g), Some(t)]) => {
            let lane = |d: &[u8]| -> Result<Vec<u16>, ZtrError> {
                Ok(words(d, 2, 2)?
                    .into_iter()
                    .map(|v| u16::try_from(v).unwrap_or(u16::MAX))
                    .collect())
            };
            [lane(&a)?, lane(&c)?, lane(&g)?, lane(&t)?]
        }
        _ => return Err(bad("the trace has no signal")),
    };
    let bases = base
        .ok_or_else(|| bad("the trace has no base calls"))?
        .get(1..)
        .unwrap_or(&[])
        .to_vec();
    let peaks: Vec<u32> = match bpos {
        Some(d) => words(&d, 4, 4)?,
        None => return Err(bad("the trace has no peak positions")),
    };
    // CNF4: the called bases' confidences first (one byte each), then the other three per base.
    let quality = cnf4.map(|d| {
        d.get(1..1 + bases.len())
            .unwrap_or(&[])
            .iter()
            .map(|&q| q as i8)
            .map(|q| u8::try_from(q.max(0)).unwrap_or(0))
            .collect::<Vec<u8>>()
    });
    let quality = quality.filter(|q| q.len() == bases.len());
    let sample_name = text.and_then(|d| {
        let body = d.get(1..)?;
        let mut parts = body.split(|b| *b == 0);
        while let Some(key) = parts.next() {
            let value = parts.next()?;
            if key == b"NAME" {
                return Some(String::from_utf8_lossy(value).into_owned());
            }
        }
        None
    });
    let c = Chromatogram {
        sample_name,
        bases,
        peaks,
        quality,
        channels,
    };
    c.check().map_err(bad)?;
    Ok(c)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn spec_examples() {
        // RLE: guard 8; 5 nines and a literal 8 (written 8 0). Our files store the length LE.
        let rle = [1, 10, 0, 0, 0, 8, 20, 8, 5, 9, 10, 9, 8, 0, 7];
        assert_eq!(unrle(&rle).unwrap(), vec![20, 9, 9, 9, 9, 9, 10, 9, 8, 7]);
        // 8-bit delta, level 1 and 2 (the spec's examples).
        assert_eq!(
            undelta(&[64, 1, 10, 10, 246, 190, 246, 71], 1, 2).unwrap(),
            vec![10, 20, 10, 200, 190, 5]
        );
        assert_eq!(
            undelta(&[64, 2, 10, 0, 236, 200, 56, 81], 1, 2).unwrap(),
            vec![10, 20, 10, 200, 190, 5]
        );
        // 16-bit delta: 0x1020 0x1f0f... the spec's 0x10 0x20 0x30 0x10 → 0x10 0x20 0x1f 0xf0.
        assert_eq!(
            undelta(&[65, 1, 0x10, 0x20, 0x1f, 0xf0], 2, 2).unwrap(),
            vec![0x10, 0x20, 0x30, 0x10]
        );
        // 16-to-8: 10, 5, -5, 200 (escaped), -800 (escaped).
        let shrunk = [70, 10, 5, 0xfb, 0x80, 0, 200, 0x80, 0xfc, 0xe0];
        assert_eq!(
            unshrink(&shrunk, 2).unwrap(),
            vec![0, 10, 0, 5, 0xff, 0xfb, 0, 200, 0xfc, 0xe0]
        );
        // 32-to-8 keeps whole words: a small value and an escaped one.
        assert_eq!(
            unshrink(&[71, 2, 0x80, 0x42, 1, 0, 0], 4).unwrap(),
            vec![0, 0, 0, 2, 0x42, 1, 0, 0]
        );
    }

    #[test]
    fn follow_and_stacks() {
        // A table that predicts each byte to follow itself: stored = prev - value.
        let mut table = vec![72];
        table.extend(0..=255u8);
        table.extend([5, 0, 255, 2]);
        assert_eq!(unfollow(&table).unwrap(), vec![5, 5, 6, 4]);
        // zlib around raw data, through `unfilter`.
        let raw = [0u8, 1, 2, 3];
        let z = miniz_oxide::deflate::compress_to_vec_zlib(&raw, 6);
        let mut chunk = vec![2];
        chunk.extend(4u32.to_le_bytes());
        chunk.extend(z);
        assert_eq!(unfilter(&chunk).unwrap(), raw);
        assert_eq!(unfilter(&[74, 0, 1, 2]), Err(ZtrError::Unsupported(74)));
    }

    /// A tiny raw trace: two bases, eight samples.
    fn tiny() -> Vec<u8> {
        let mut f = MAGIC.to_vec();
        f.extend([1, 2]);
        let mut chunk = |kind: &[u8; 4], data: Vec<u8>| {
            f.extend(kind);
            f.extend(0u32.to_be_bytes());
            f.extend(u32::try_from(data.len()).unwrap().to_be_bytes());
            f.extend(data);
        };
        let mut smp = vec![0, 0];
        for ch in 0..4u16 {
            for k in 0..2u16 {
                smp.extend((ch * 100 + k).to_be_bytes());
            }
        }
        chunk(b"SMP4", smp);
        chunk(b"BASE", vec![0, b'A', b'T']);
        chunk(b"BPOS", vec![0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1]);
        chunk(b"CNF4", vec![0, 30, 20, 1, 2, 3, 4, 5, 6]);
        chunk(b"TEXT", b"\0NAME\0read 7\0\0".to_vec());
        chunk(b"COMM", vec![0, b'x']);
        f
    }

    #[test]
    fn reads_a_raw_trace() {
        let c = read_chromatogram(&tiny()).unwrap();
        assert_eq!(c.bases, b"AT");
        assert_eq!(c.peaks, vec![0, 1]);
        assert_eq!(c.quality, Some(vec![30, 20]));
        assert_eq!(c.channels[2], vec![200, 201]);
        assert_eq!(c.sample_name.as_deref(), Some("read 7"));
        assert_eq!(read_bases(&tiny()).unwrap(), b"AT");
    }

    #[test]
    fn hostile_traces() {
        assert_eq!(read_chromatogram(b"ZTR"), Err(ZtrError::NotZtr));
        let good = tiny();
        for cut in 0..good.len() {
            let _ = read_chromatogram(&good[..cut]);
        }
        let mut v2 = good.clone();
        v2[8] = 2;
        assert_eq!(read_chromatogram(&v2), Err(ZtrError::Version(2)));
        // A zlib chunk claiming a huge size.
        let mut huge = vec![2];
        huge.extend(u32::MAX.to_le_bytes());
        assert!(unfilter(&huge).is_err());
        // RLE runs longer than declared.
        assert!(unrle(&[1, 1, 0, 0, 0, 8, 8, 200, 1]).is_err());
        assert!(unfilter(&[]).is_err());
    }
}
