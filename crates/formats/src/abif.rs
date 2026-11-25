//! ABIF, the Applied Biosystems sequencing trace format (`.ab1`, `.abi`, `.fsa`), as described
//! in Applied Biosystems' public "ABIF File Format Specification and Sample File Schema" (2009).
//!
//! A file starts with the four bytes `ABIF`, a version number and the directory entry of the
//! root ("tdir"), which points at a directory of 28-byte entries. Each entry names a tag (four
//! characters plus a number, such as `PBAS` 2), the type and count of its elements and where
//! its data is. Data of four bytes or less sits in the entry's offset field itself. Every
//! number is big-endian.
//!
//! The tags a chromatogram needs:
//! - `DATA` 9–12: the analysed signal, one channel each, in the base order given by `FWO_` 1;
//! - `PBAS` 1 (edited by the user) or 2 (as called): the bases;
//! - `PLOC` 1 or 2: each base's peak position in the signal;
//! - `PCON` 1 or 2: each base's quality value;
//! - `SMPL` 1: the sample name.
//!
//! The file is kept as bytes by the caller; this module only reads it, so nothing is lost.

use crate::trace::{CHANNEL_BASES, Chromatogram};

/// Why an ABIF file couldn't be read.
#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum AbifError {
    #[error("not an ABIF file: it must start with \"ABIF\"")]
    NotAbif,

    #[error("the ABIF file is damaged: {0}")]
    Malformed(String),

    #[error(
        "the file has no base calls: it holds no sequencing read (a fragment analysis file, \
         perhaps)"
    )]
    NoBaseCalls,
}

/// One directory entry.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Entry {
    /// The tag name, four bytes such as `PBAS`.
    pub name: [u8; 4],
    /// The tag number (the 2 in `PBAS` 2).
    pub number: i32,
    /// The element type code (2 = char, 4 = short, 5 = long, 18 = pString, ...).
    pub element_type: i16,
    pub element_size: i16,
    pub elements: i32,
    pub data_size: i32,
    pub data_offset: i32,
    /// Where the entry is in the file (its offset field holds data of four bytes or less).
    position: usize,
}

/// A parsed ABIF directory over the file's bytes.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AbifFile<'a> {
    bytes: &'a [u8],
    pub version: i16,
    pub entries: Vec<Entry>,
}

const ENTRY_SIZE: usize = 28;

fn be_i16(b: &[u8], at: usize) -> Option<i16> {
    Some(i16::from_be_bytes(
        b.get(at..at.checked_add(2)?)?.try_into().ok()?,
    ))
}

fn be_i32(b: &[u8], at: usize) -> Option<i32> {
    Some(i32::from_be_bytes(
        b.get(at..at.checked_add(4)?)?.try_into().ok()?,
    ))
}

fn read_entry(b: &[u8], at: usize) -> Option<Entry> {
    let field = |i: usize| at.checked_add(i);
    Some(Entry {
        name: b.get(at..field(4)?)?.try_into().ok()?,
        number: be_i32(b, field(4)?)?,
        element_type: be_i16(b, field(8)?)?,
        element_size: be_i16(b, field(10)?)?,
        elements: be_i32(b, field(12)?)?,
        data_size: be_i32(b, field(16)?)?,
        data_offset: be_i32(b, field(20)?)?,
        position: at,
    })
}

impl<'a> AbifFile<'a> {
    /// Read the header and directory. Entries are checked only when their data is read.
    pub fn parse(bytes: &'a [u8]) -> Result<Self, AbifError> {
        if bytes.get(..4) != Some(b"ABIF") {
            return Err(AbifError::NotAbif);
        }
        let bad = |m: &str| AbifError::Malformed(m.to_string());
        let version = be_i16(bytes, 4).ok_or_else(|| bad("the header is cut short"))?;
        let root = read_entry(bytes, 6).ok_or_else(|| bad("the header is cut short"))?;
        let count = usize::try_from(root.elements)
            .map_err(|_| bad("the directory has a negative number of entries"))?;
        let start = usize::try_from(root.data_offset)
            .map_err(|_| bad("the directory has a negative offset"))?;
        // Each entry needs 28 bytes of the file, which bounds `count` by the input size.
        let end = count
            .checked_mul(ENTRY_SIZE)
            .and_then(|n| n.checked_add(start))
            .filter(|&end| end <= bytes.len())
            .ok_or_else(|| bad("the directory runs past the end of the file"))?;
        let mut entries = Vec::with_capacity(count);
        let mut at = start;
        while at < end {
            entries
                .push(read_entry(bytes, at).ok_or_else(|| bad("a directory entry is cut short"))?);
            at = at.saturating_add(ENTRY_SIZE);
        }
        Ok(AbifFile {
            bytes,
            version,
            entries,
        })
    }

    /// The entry for tag `name` number `number`.
    pub fn entry(&self, name: &[u8; 4], number: i32) -> Option<&Entry> {
        self.entries
            .iter()
            .find(|e| &e.name == name && e.number == number)
    }

    /// The data bytes of an entry.
    pub fn data(&self, e: &Entry) -> Result<&'a [u8], AbifError> {
        let tag = tag_name(e);
        let bad = |m: String| AbifError::Malformed(format!("tag {tag}: {m}"));
        let size = usize::try_from(e.data_size).map_err(|_| bad("negative data size".into()))?;
        let elements =
            usize::try_from(e.elements).map_err(|_| bad("negative element count".into()))?;
        let element_size =
            usize::try_from(e.element_size).map_err(|_| bad("negative element size".into()))?;
        if elements.checked_mul(element_size).is_none_or(|n| n > size) {
            return Err(bad("its elements don't fit in its data".into()));
        }
        if size <= 4 {
            let field = e.position.saturating_add(20);
            return self
                .bytes
                .get(field..field.saturating_add(size))
                .ok_or_else(|| bad("the entry is cut short".into()));
        }
        let start =
            usize::try_from(e.data_offset).map_err(|_| bad("negative data offset".into()))?;
        start
            .checked_add(size)
            .and_then(|end| self.bytes.get(start..end))
            .ok_or_else(|| bad("its data runs past the end of the file".into()))
    }

    /// A tag's elements as integers (byte, char, word, short or long elements).
    pub fn integers(&self, e: &Entry) -> Result<Vec<i64>, AbifError> {
        let data = self.data(e)?;
        let width = match e.element_type {
            1 | 2 => 1,
            3 | 4 => 2,
            5 => 4,
            t => {
                return Err(AbifError::Malformed(format!(
                    "tag {} has element type {t}, not numbers",
                    tag_name(e)
                )));
            }
        };
        let count = usize::try_from(e.elements).unwrap_or(0);
        Ok(data
            .chunks_exact(width)
            .take(count)
            .map(|c| match (e.element_type, c) {
                (1 | 2, [a]) => i64::from(*a),
                (3, [a, b]) => i64::from(u16::from_be_bytes([*a, *b])),
                (4, [a, b]) => i64::from(i16::from_be_bytes([*a, *b])),
                (5, [a, b, c, d]) => i64::from(i32::from_be_bytes([*a, *b, *c, *d])),
                _ => 0,
            })
            .collect())
    }

    /// A tag's text: a pString (length byte first), a cString (NUL-terminated) or a char
    /// array. Bytes that aren't UTF-8 are read as Latin-1.
    pub fn text(&self, e: &Entry) -> Result<String, AbifError> {
        let data = self.data(e)?;
        let raw = match e.element_type {
            18 => {
                let n = usize::from(data.first().copied().unwrap_or(0));
                data.get(1..)
                    .unwrap_or(&[])
                    .get(..n)
                    .unwrap_or(data.get(1..).unwrap_or(&[]))
            }
            19 => data.split(|b| *b == 0).next().unwrap_or(&[]),
            _ => data,
        };
        Ok(match std::str::from_utf8(raw) {
            Ok(s) => s.to_string(),
            Err(_) => raw.iter().map(|&b| char::from(b)).collect(),
        })
    }
}

fn tag_name(e: &Entry) -> String {
    let name: String = e
        .name
        .iter()
        .map(|&b| {
            if b.is_ascii_graphic() {
                char::from(b)
            } else {
                '?'
            }
        })
        .collect();
    format!("{name} {}", e.number)
}

/// Read the chromatogram of a sequencing trace. The user-edited calls (`PBAS` 1 with `PLOC` 1
/// and `PCON` 1) are preferred when present; otherwise the basecaller's (number 2).
pub fn read_chromatogram(bytes: &[u8]) -> Result<Chromatogram, AbifError> {
    let file = AbifFile::parse(bytes)?;
    let bad = |m: String| AbifError::Malformed(m);
    let number = if file.entry(b"PBAS", 1).is_some() {
        1
    } else {
        2
    };
    let bases_entry = file.entry(b"PBAS", number).ok_or(AbifError::NoBaseCalls)?;
    let bases = file.data(bases_entry)?.to_vec();
    if bases.is_empty() {
        return Err(AbifError::NoBaseCalls);
    }
    let peaks = file
        .entry(b"PLOC", number)
        .ok_or_else(|| bad("the bases have no peak positions (PLOC)".into()))?;
    let peaks = file
        .integers(peaks)?
        .into_iter()
        .map(|p| u32::try_from(p).map_err(|_| bad(format!("negative peak position {p}"))))
        .collect::<Result<Vec<u32>, _>>()?;
    let quality = match file.entry(b"PCON", number) {
        Some(e) => Some(
            file.integers(e)?
                .into_iter()
                .map(|q| u8::try_from(q.clamp(0, 255)).unwrap_or(0))
                .collect(),
        ),
        None => None,
    };
    let order = match file.entry(b"FWO_", 1) {
        Some(e) => file.data(e)?.to_vec(),
        None => b"GATC".to_vec(),
    };
    let mut channels: [Vec<u16>; 4] = Default::default();
    let mut seen = [false; 4];
    for (i, base) in order.iter().take(4).enumerate() {
        let slot = CHANNEL_BASES
            .iter()
            .position(|b| b.eq_ignore_ascii_case(base))
            .ok_or_else(|| {
                bad(format!(
                    "unknown base {:?} in the channel order (FWO_)",
                    char::from(*base)
                ))
            })?;
        if std::mem::replace(
            seen.get_mut(slot).ok_or_else(|| bad("channel".into()))?,
            true,
        ) {
            return Err(bad(
                "a base appears twice in the channel order (FWO_)".into()
            ));
        }
        let tag = i32::try_from(i).unwrap_or(0).saturating_add(9);
        let e = file
            .entry(b"DATA", tag)
            .ok_or_else(|| bad(format!("the signal channel DATA {tag} is missing")))?;
        // Analysed signal can dip just below zero after baseline correction; it is drawn from 0.
        let values = file
            .integers(e)?
            .into_iter()
            .map(|v| u16::try_from(v.clamp(0, i64::from(u16::MAX))).unwrap_or(0))
            .collect();
        if let Some(c) = channels.get_mut(slot) {
            *c = values;
        }
    }
    if seen.contains(&false) {
        return Err(bad(
            "the channel order (FWO_) doesn't name all four bases".into()
        ));
    }
    let sample_name = match file.entry(b"SMPL", 1) {
        Some(e) => Some(file.text(e)?).filter(|s| !s.trim().is_empty()),
        None => None,
    };
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

/// Write a chromatogram as a minimal ABIF file: the basecaller tags (number 2), the four
/// signal channels, the channel order and the sample name. Used for synthetic test traces;
/// files read from disk are saved from their original bytes.
pub fn write_chromatogram(c: &Chromatogram) -> Result<Vec<u8>, String> {
    c.check()?;
    struct Tag {
        name: [u8; 4],
        number: i32,
        element_type: i16,
        element_size: i16,
        elements: usize,
        data: Vec<u8>,
    }
    let shorts = |v: &[u16]| -> Result<Vec<u8>, String> {
        let mut out = Vec::with_capacity(v.len().saturating_mul(2));
        for &x in v {
            let x = i16::try_from(x).map_err(|_| format!("signal value {x} is above 32767"))?;
            out.extend_from_slice(&x.to_be_bytes());
        }
        Ok(out)
    };
    let mut tags = Vec::new();
    // Channels in the file's usual order, G A T C.
    for (i, base) in b"GATC".iter().enumerate() {
        let slot = CHANNEL_BASES.iter().position(|b| b == base).unwrap_or(0);
        let values = c.channels.get(slot).map(Vec::as_slice).unwrap_or(&[]);
        tags.push(Tag {
            name: *b"DATA",
            number: 9 + i32::try_from(i).unwrap_or(0),
            element_type: 4,
            element_size: 2,
            elements: values.len(),
            data: shorts(values)?,
        });
    }
    tags.push(Tag {
        name: *b"FWO_",
        number: 1,
        element_type: 2,
        element_size: 1,
        elements: 4,
        data: b"GATC".to_vec(),
    });
    tags.push(Tag {
        name: *b"PBAS",
        number: 2,
        element_type: 2,
        element_size: 1,
        elements: c.bases.len(),
        data: c.bases.clone(),
    });
    let mut ploc = Vec::new();
    for &p in &c.peaks {
        let p = i16::try_from(p).map_err(|_| format!("peak position {p} is above 32767"))?;
        ploc.extend_from_slice(&p.to_be_bytes());
    }
    tags.push(Tag {
        name: *b"PLOC",
        number: 2,
        element_type: 4,
        element_size: 2,
        elements: c.peaks.len(),
        data: ploc,
    });
    if let Some(q) = &c.quality {
        tags.push(Tag {
            name: *b"PCON",
            number: 2,
            element_type: 2,
            element_size: 1,
            elements: q.len(),
            data: q.clone(),
        });
    }
    if let Some(name) = &c.sample_name {
        let bytes = name.as_bytes();
        let n = u8::try_from(bytes.len()).map_err(|_| "the sample name is too long".to_string())?;
        let mut data = vec![n];
        data.extend_from_slice(bytes);
        tags.push(Tag {
            name: *b"SMPL",
            number: 1,
            element_type: 18,
            element_size: 1,
            elements: data.len(),
            data,
        });
    }
    let too_big = || "the trace is too large for an ABIF file".to_string();
    let i32_of = |n: usize| i32::try_from(n).map_err(|_| too_big());
    let mut body = Vec::new();
    let mut directory = Vec::new();
    for t in &tags {
        let mut offset = [0u8; 4];
        if t.data.len() <= 4 {
            for (o, d) in offset.iter_mut().zip(&t.data) {
                *o = *d;
            }
        } else {
            offset = i32_of(128usize.saturating_add(body.len()))?.to_be_bytes();
            body.extend_from_slice(&t.data);
        }
        directory.extend_from_slice(&t.name);
        directory.extend_from_slice(&t.number.to_be_bytes());
        directory.extend_from_slice(&t.element_type.to_be_bytes());
        directory.extend_from_slice(&t.element_size.to_be_bytes());
        directory.extend_from_slice(&i32_of(t.elements)?.to_be_bytes());
        directory.extend_from_slice(&i32_of(t.data.len())?.to_be_bytes());
        directory.extend_from_slice(&offset);
        directory.extend_from_slice(&0i32.to_be_bytes());
    }
    let mut out = Vec::with_capacity(128 + body.len() + directory.len());
    out.extend_from_slice(b"ABIF");
    out.extend_from_slice(&101i16.to_be_bytes());
    out.extend_from_slice(b"tdir");
    out.extend_from_slice(&1i32.to_be_bytes());
    out.extend_from_slice(&1023i16.to_be_bytes());
    out.extend_from_slice(&28i16.to_be_bytes());
    out.extend_from_slice(&i32_of(tags.len())?.to_be_bytes());
    out.extend_from_slice(&i32_of(directory.len())?.to_be_bytes());
    out.extend_from_slice(&i32_of(128usize.saturating_add(body.len()))?.to_be_bytes());
    out.extend_from_slice(&0i32.to_be_bytes());
    out.resize(128, 0);
    out.extend_from_slice(&body);
    out.extend_from_slice(&directory);
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    pub(crate) fn sample() -> Chromatogram {
        // Four bases, peaks every 10 samples, each base's channel peaking at its peak.
        let bases = b"GATC".to_vec();
        let peaks = vec![5, 15, 25, 35];
        let mut channels: [Vec<u16>; 4] = Default::default();
        for c in &mut channels {
            *c = vec![0; 40];
        }
        for (b, p) in bases.iter().zip(&peaks) {
            let slot = CHANNEL_BASES.iter().position(|x| x == b).unwrap();
            for d in 0..5u16 {
                let v = 1000 - d * 200;
                channels[slot][*p as usize + d as usize] = v;
                channels[slot][*p as usize - d as usize] = v;
            }
        }
        Chromatogram {
            sample_name: Some("read 1".into()),
            bases,
            peaks,
            quality: Some(vec![12, 40, 3, 60]),
            channels,
        }
    }

    #[test]
    fn written_file_reads_back() {
        let c = sample();
        let bytes = write_chromatogram(&c).unwrap();
        assert_eq!(&bytes[..4], b"ABIF");
        assert_eq!(read_chromatogram(&bytes).unwrap(), c);
    }

    #[test]
    fn short_data_sits_in_the_offset_field() {
        // FWO_ (4 bytes) is stored inline: its offset field holds "GATC".
        let bytes = write_chromatogram(&sample()).unwrap();
        let file = AbifFile::parse(&bytes).unwrap();
        let e = file.entry(b"FWO_", 1).unwrap();
        assert_eq!(&e.data_offset.to_be_bytes(), b"GATC");
        assert_eq!(file.data(e).unwrap(), b"GATC");
    }

    #[test]
    fn channel_order_follows_fwo() {
        // Rewrite FWO_ to "ACGT": DATA 9 is then read as A, so A gets the G signal.
        let c = sample();
        let mut bytes = write_chromatogram(&c).unwrap();
        // Found through the directory: the four bases (PBAS) are stored inline as "GATC" too.
        let at = AbifFile::parse(&bytes)
            .unwrap()
            .entry(b"FWO_", 1)
            .unwrap()
            .position
            + 20;
        bytes[at..at + 4].copy_from_slice(b"ACGT");
        let r = read_chromatogram(&bytes).unwrap();
        assert_eq!(r.channels[0], c.channels[2]);
    }

    #[test]
    fn edited_calls_win() {
        let mut bytes = write_chromatogram(&sample()).unwrap();
        // Rename PBAS/PLOC/PCON 2 to number 1; reading still works and uses them.
        let mut file_entries = Vec::new();
        {
            let f = AbifFile::parse(&bytes).unwrap();
            let dir = 128
                + f.entries
                    .iter()
                    .filter(|e| e.data_size > 4)
                    .map(|e| e.data_size as usize)
                    .sum::<usize>();
            for (i, e) in f.entries.iter().enumerate() {
                if [*b"PBAS", *b"PLOC", *b"PCON"].contains(&e.name) {
                    file_entries.push(dir + i * 28 + 4);
                }
            }
        }
        for at in file_entries {
            bytes[at..at + 4].copy_from_slice(&1i32.to_be_bytes());
        }
        assert_eq!(read_chromatogram(&bytes).unwrap().bases, b"GATC");
    }

    #[test]
    fn hostile_files_are_errors() {
        assert_eq!(read_chromatogram(b""), Err(AbifError::NotAbif));
        assert!(matches!(
            read_chromatogram(b"ABIF"),
            Err(AbifError::Malformed(_))
        ));
        let good = write_chromatogram(&sample()).unwrap();
        for cut in [10, 127, 128, 200, good.len() - 1] {
            assert!(read_chromatogram(&good[..cut]).is_err(), "cut at {cut}");
        }
        // A directory claiming two billion entries.
        let mut big = good.clone();
        big[18..22].copy_from_slice(&i32::MAX.to_be_bytes());
        assert!(read_chromatogram(&big).is_err());
    }

    #[test]
    fn no_base_calls() {
        let mut c = sample();
        c.bases.clear();
        c.peaks.clear();
        c.quality = Some(vec![]);
        let bytes = write_chromatogram(&c).unwrap();
        assert_eq!(read_chromatogram(&bytes), Err(AbifError::NoBaseCalls));
    }
}
