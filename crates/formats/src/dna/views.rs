//! Typed views of the cookie, sequence and notes packets.

use super::{DnaError, DnaFile, MAGIC, tag};
use crate::xml::{self, Token, Tokenizer};

/// The cookie: what kind of sequence the file holds and which format versions wrote it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Cookie {
    /// 1 = DNA (seen in every file so far).
    pub sequence_type: u16,
    /// *Hypothesis:* the format version of the writer (15 in all files so far).
    pub export_version: u16,
    /// *Hypothesis:* the minimum reader version (19 or 20 so far).
    pub import_version: u16,
}

/// The sequence packet: one flags byte, then the residues.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SequencePacket {
    /// Bit 0: circular (confirmed on linear and circular files). Other bits are kept as read.
    pub flags: u8,
    pub residues: Vec<u8>,
}

impl SequencePacket {
    pub const CIRCULAR: u8 = 0x01;

    pub fn is_circular(&self) -> bool {
        self.flags & Self::CIRCULAR != 0
    }

    pub fn set_circular(&mut self, circular: bool) {
        if circular {
            self.flags |= Self::CIRCULAR;
        } else {
            self.flags &= !Self::CIRCULAR;
        }
    }

    pub fn to_payload(&self) -> Vec<u8> {
        let mut out = Vec::with_capacity(self.residues.len().saturating_add(1));
        out.push(self.flags);
        out.extend_from_slice(&self.residues);
        out
    }
}

/// One element of the notes packet, such as `<Description>…</Description>`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct NoteField {
    pub name: String,
    /// Attributes, decoded, in file order (for example `UTC` on dates).
    pub attrs: Vec<(String, String)>,
    /// The text content, decoded.
    pub text: String,
}

impl DnaFile {
    pub fn cookie(&self) -> Result<Cookie, DnaError> {
        let p = self.payload(tag::COOKIE)?;
        let invalid = |reason: &str| DnaError::InvalidPacket {
            tag: tag::COOKIE,
            reason: reason.into(),
        };
        let fields = p
            .strip_prefix(MAGIC.as_slice())
            .ok_or_else(|| invalid("missing magic"))?;
        let u16_at = |i: usize| {
            fields
                .get(i..i.saturating_add(2))
                .and_then(|b| <[u8; 2]>::try_from(b).ok())
                .map(u16::from_be_bytes)
                .ok_or_else(|| invalid("shorter than 14 bytes"))
        };
        Ok(Cookie {
            sequence_type: u16_at(0)?,
            export_version: u16_at(2)?,
            import_version: u16_at(4)?,
        })
    }

    pub fn sequence(&self) -> Result<SequencePacket, DnaError> {
        let p = self.payload(tag::SEQUENCE)?;
        let (&flags, residues) = p.split_first().ok_or(DnaError::InvalidPacket {
            tag: tag::SEQUENCE,
            reason: "empty (no flags byte)".into(),
        })?;
        Ok(SequencePacket {
            flags,
            residues: residues.to_vec(),
        })
    }

    /// The notes as a flat list of fields, in file order. Read-only: the packet itself is
    /// always kept as read.
    pub fn notes(&self) -> Result<Vec<NoteField>, DnaError> {
        let text = utf8(self.payload(tag::NOTES)?, tag::NOTES)?;
        parse_notes(text).map_err(|e| DnaError::InvalidPacket {
            tag: tag::NOTES,
            reason: e.to_string(),
        })
    }
}

pub(crate) fn utf8(payload: &[u8], tag: u8) -> Result<&str, DnaError> {
    std::str::from_utf8(payload).map_err(|e| DnaError::InvalidPacket {
        tag,
        reason: format!("not UTF-8: {e}"),
    })
}

/// Fields are the elements directly inside `<Notes>`; nested markup is flattened into text.
fn parse_notes(text: &str) -> Result<Vec<NoteField>, xml::XmlError> {
    let mut tokens = Tokenizer::new(text);
    let mut fields = Vec::new();
    let mut depth = 0usize;
    let mut current: Option<NoteField> = None;
    while let Some(tok) = tokens.next_token() {
        match tok?.0 {
            Token::Start {
                name,
                attrs,
                self_closing,
            } => {
                if depth == 1 {
                    let field = NoteField {
                        name: name.to_string(),
                        attrs: attrs
                            .iter()
                            .map(|a| (a.name.to_string(), xml::decode(a.raw)))
                            .collect(),
                        text: String::new(),
                    };
                    if self_closing {
                        fields.push(field);
                    } else {
                        current = Some(field);
                    }
                }
                if !self_closing {
                    depth = depth.saturating_add(1);
                }
            }
            Token::End { .. } => {
                depth = depth.saturating_sub(1);
                if depth == 1
                    && let Some(field) = current.take()
                {
                    fields.push(field);
                }
            }
            Token::Text(t) => {
                if let Some(field) = current.as_mut() {
                    field.text.push_str(&xml::decode(t));
                }
            }
            Token::Other(_) => {}
        }
    }
    Ok(fields)
}
