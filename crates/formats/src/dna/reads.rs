//! Sequencing reads matched to a construct, as the reference application stores them in
//! `.dna` files: packet 0x11 lists them (`<AlignableSequences>`, one `<Sequence>` each with its
//! `name`, `ID`, `isTrace` and `trimmedRange`, 0-based and inclusive), and each read is a 0x10
//! packet: u32 ID, u8 kind (0x12: ZTR, [`crate::ztr`]), u32 length, then the data. A trace has
//! its signal; a plain sequence (an assembly, say) is a ZTR with base calls only. Learned from
//! `.dna` files we authored ourselves; the packets are kept as read and written back unchanged.

use super::{DnaFile, tag, views::utf8};
use crate::trace::Chromatogram;
use crate::xml::{self, Token, Tokenizer};

/// Most matched reads read from one file.
const MAX_READS: usize = 1000;

/// One matched read.
#[derive(Debug, Clone, PartialEq)]
pub struct MatchedRead {
    pub id: u32,
    pub name: String,
    /// The stretch kept after trimming: 0-based, end exclusive (stored inclusive).
    pub trimmed: Option<(usize, usize)>,
    /// The read, or why it can't be read (a kind or compression we don't know).
    pub read: Result<ReadData, String>,
}

/// A matched read's bases, and its trace if it is one.
#[derive(Debug, Clone, PartialEq)]
pub struct ReadData {
    pub bases: Vec<u8>,
    pub trace: Option<Chromatogram>,
}

/// `a..b` (0-based, inclusive) as `a..b+1`.
fn range(text: &str) -> Option<(usize, usize)> {
    let (a, b) = text.split_once("..")?;
    let a: usize = a.trim().parse().ok()?;
    let b: usize = b.trim().parse().ok()?;
    (a <= b).then(|| (a, b.saturating_add(1)))
}

impl DnaFile {
    /// The reads matched to this construct, in list order (empty if there are none).
    pub fn matched_reads(&self) -> Result<Vec<MatchedRead>, super::DnaError> {
        let Some(list) = self.packet(tag::ALIGNABLE) else {
            return Ok(Vec::new());
        };
        let text = utf8(&list.payload, tag::ALIGNABLE)?;
        let bad = |reason: String| super::DnaError::InvalidPacket {
            tag: tag::ALIGNABLE,
            reason,
        };
        let mut out = Vec::new();
        let mut tokens = Tokenizer::new(text);
        while let Some(tok) = tokens.next_token() {
            let (token, _) = tok.map_err(|e| bad(e.to_string()))?;
            let Token::Start {
                name: "Sequence",
                attrs,
                ..
            } = token
            else {
                continue;
            };
            let get = |n: &str| {
                attrs
                    .iter()
                    .find(|a| a.name == n)
                    .map(|a| xml::decode(a.raw))
            };
            let Some(id) = get("ID").and_then(|v| v.parse::<u32>().ok()) else {
                continue;
            };
            if out.len() >= MAX_READS {
                break;
            }
            let is_trace = get("isTrace").as_deref() == Some("1");
            out.push(MatchedRead {
                id,
                name: get("name").unwrap_or_else(|| format!("read {id}")),
                trimmed: get("trimmedRange").as_deref().and_then(range),
                read: self.read_payload(id, is_trace),
            });
        }
        Ok(out)
    }

    /// The data stored for read `id` (packet 0x10).
    fn read_payload(&self, id: u32, is_trace: bool) -> Result<ReadData, String> {
        for p in self.packets().iter().filter(|p| p.tag == tag::TRACE) {
            let d = &p.payload;
            let Some(pid) = d
                .get(0..4)
                .and_then(|b| b.try_into().ok())
                .map(u32::from_be_bytes)
            else {
                continue;
            };
            if pid != id {
                continue;
            }
            let kind = d.get(4).copied().unwrap_or(0);
            let len = d
                .get(5..9)
                .and_then(|b| b.try_into().ok())
                .map(u32::from_be_bytes)
                .unwrap_or(0) as usize;
            let body = d
                .get(9..9usize.saturating_add(len))
                .ok_or_else(|| "the trace is cut short".to_string())?;
            if kind != 0x12 {
                return Err(format!("a read of kind {kind:#04x}, which isn't known yet"));
            }
            return if is_trace {
                let c = crate::ztr::read_chromatogram(body).map_err(|e| e.to_string())?;
                Ok(ReadData {
                    bases: c.bases.clone(),
                    trace: Some(c),
                })
            } else {
                Ok(ReadData {
                    bases: crate::ztr::read_bases(body).map_err(|e| e.to_string())?,
                    trace: None,
                })
            };
        }
        Err("its data isn't in the file".to_string())
    }
}

#[cfg(test)]
mod tests {
    use super::super::{MAGIC, Packet};
    use super::*;

    #[test]
    fn ranges() {
        assert_eq!(range("17..868"), Some((17, 869)));
        assert_eq!(range("0..4194"), Some((0, 4195)));
        assert_eq!(range("0..0"), Some((0, 1)));
        assert_eq!(range("9..5"), None);
        assert_eq!(range("x"), None);
    }

    #[test]
    fn listed_reads_without_traces() {
        let mut cookie = MAGIC.to_vec();
        cookie.extend_from_slice(&[0, 1, 0, 15, 0, 19]);
        let list = r#"<AlignableSequences trimStringency="Medium"><Sequence name="r1" ID="0" isTrace="1" trimmedRange="5..10"/><Sequence name="s2" ID="1" isTrace="0"/><Sequence name="bad"/></AlignableSequences>"#;
        let f = DnaFile::from_packets(vec![
            Packet {
                tag: tag::COOKIE,
                payload: cookie,
            },
            Packet {
                tag: tag::SEQUENCE,
                payload: b"\x00ACGT".to_vec(),
            },
            Packet {
                tag: tag::ALIGNABLE,
                payload: list.as_bytes().to_vec(),
            },
            Packet {
                tag: tag::TRACE,
                payload: vec![0, 0, 0, 0, 0x99, 0, 0, 0, 1, 7],
            },
        ])
        .unwrap();
        let reads = f.matched_reads().unwrap();
        assert_eq!(reads.len(), 2);
        assert_eq!(reads[0].name, "r1");
        assert_eq!(reads[0].trimmed, Some((5, 11)));
        assert!(reads[0].read.as_ref().unwrap_err().contains("0x99"));
        assert!(
            reads[1]
                .read
                .as_ref()
                .unwrap_err()
                .contains("isn't in the file")
        );
    }
}
