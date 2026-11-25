use super::{MAGIC, tag};

/// Why a `.dna` file couldn't be read or written.
#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum DnaError {
    #[error("not a SnapGene .dna file (no \"SnapGene\" cookie at the start)")]
    NotDna,

    #[error(
        "file is truncated: packet 0x{tag:02X} at byte {offset} declares {declared} bytes, but only {available} remain"
    )]
    Truncated {
        offset: usize,
        tag: u8,
        declared: usize,
        available: usize,
    },

    #[error("file ends inside a packet header at byte {offset}")]
    TruncatedHeader { offset: usize },

    #[error("packet 0x{tag:02X} is {len} bytes, more than the format's 4 GiB limit")]
    TooLarge { tag: u8, len: usize },

    #[error("packet 0x{tag:02X} is invalid: {reason}")]
    InvalidPacket { tag: u8, reason: String },

    #[error("the file has no packet 0x{tag:02X}")]
    MissingPacket { tag: u8 },
}

/// One packet: its tag and its payload exactly as stored.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Packet {
    pub tag: u8,
    pub payload: Vec<u8>,
}

/// A `.dna` file as an ordered list of packets.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DnaFile {
    packets: Vec<Packet>,
}

/// Bytes in a packet header: tag + big-endian u32 length.
const HEADER: usize = 5;

impl DnaFile {
    /// Read a file. Every declared length is checked against the input before anything is
    /// copied, so hostile input can't cause a panic or an oversized allocation.
    pub fn parse(bytes: &[u8]) -> Result<Self, DnaError> {
        let mut packets = Vec::new();
        let mut offset = 0usize;
        while offset < bytes.len() {
            let header = bytes
                .get(offset..offset.saturating_add(HEADER))
                .ok_or(DnaError::TruncatedHeader { offset })?;
            let (&tag, len_bytes) = header
                .split_first()
                .ok_or(DnaError::TruncatedHeader { offset })?;
            let len_array: [u8; 4] = len_bytes
                .try_into()
                .map_err(|_| DnaError::TruncatedHeader { offset })?;
            let declared = usize::try_from(u32::from_be_bytes(len_array)).unwrap_or(usize::MAX);
            let start = offset.saturating_add(HEADER);
            let available = bytes.len().saturating_sub(start);
            let payload = start
                .checked_add(declared)
                .and_then(|end| bytes.get(start..end))
                .ok_or(DnaError::Truncated {
                    offset,
                    tag,
                    declared,
                    available,
                })?;
            packets.push(Packet {
                tag,
                payload: payload.to_vec(),
            });
            offset = start.saturating_add(declared);
        }
        let file = DnaFile { packets };
        if !file.has_cookie() {
            return Err(DnaError::NotDna);
        }
        Ok(file)
    }

    /// Build a file from packets. The first must be a cookie.
    pub fn from_packets(packets: Vec<Packet>) -> Result<Self, DnaError> {
        let file = DnaFile { packets };
        if file.has_cookie() {
            Ok(file)
        } else {
            Err(DnaError::NotDna)
        }
    }

    fn has_cookie(&self) -> bool {
        self.packets
            .first()
            .is_some_and(|p| p.tag == tag::COOKIE && p.payload.starts_with(MAGIC))
    }

    /// The file's bytes. For a file read with [`DnaFile::parse`] and not modified, this is
    /// exactly the input.
    pub fn to_bytes(&self) -> Result<Vec<u8>, DnaError> {
        let total = self
            .packets
            .iter()
            .map(|p| p.payload.len().saturating_add(HEADER))
            .fold(0usize, usize::saturating_add);
        let mut out = Vec::with_capacity(total);
        for p in &self.packets {
            let len = u32::try_from(p.payload.len()).map_err(|_| DnaError::TooLarge {
                tag: p.tag,
                len: p.payload.len(),
            })?;
            out.push(p.tag);
            out.extend_from_slice(&len.to_be_bytes());
            out.extend_from_slice(&p.payload);
        }
        Ok(out)
    }

    pub fn packets(&self) -> &[Packet] {
        &self.packets
    }

    /// The first packet with `tag`.
    pub fn packet(&self, tag: u8) -> Option<&Packet> {
        self.packets.iter().find(|p| p.tag == tag)
    }

    /// The payload of the first packet with `tag`, or an error naming the missing packet.
    pub fn payload(&self, tag: u8) -> Result<&[u8], DnaError> {
        self.packet(tag)
            .map(|p| p.payload.as_slice())
            .ok_or(DnaError::MissingPacket { tag })
    }

    /// Replace the payload of the first packet with `tag`, keeping its position, or append a
    /// new packet if there is none. The cookie can't be replaced by something without the
    /// magic.
    pub fn set_payload(&mut self, tag: u8, payload: Vec<u8>) -> Result<(), DnaError> {
        if tag == tag::COOKIE && !payload.starts_with(MAGIC) {
            return Err(DnaError::InvalidPacket {
                tag,
                reason: "the cookie must start with \"SnapGene\"".into(),
            });
        }
        match self.packets.iter_mut().find(|p| p.tag == tag) {
            Some(p) => p.payload = payload,
            None => self.packets.push(Packet { tag, payload }),
        }
        Ok(())
    }

    /// Remove every packet with `tag` (never the cookie). Returns how many were removed.
    pub fn remove(&mut self, tag: u8) -> usize {
        if tag == tag::COOKIE {
            return 0;
        }
        let before = self.packets.len();
        self.packets.retain(|p| p.tag != tag);
        before.saturating_sub(self.packets.len())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn packet(tag: u8, payload: &[u8]) -> Vec<u8> {
        let mut v = vec![tag];
        v.extend_from_slice(&(payload.len() as u32).to_be_bytes());
        v.extend_from_slice(payload);
        v
    }

    fn cookie() -> Vec<u8> {
        packet(tag::COOKIE, b"SnapGene\x00\x01\x00\x0f\x00\x13")
    }

    #[test]
    fn round_trips_known_unknown_empty_and_repeated_packets() {
        let mut bytes = cookie();
        bytes.extend(packet(tag::SEQUENCE, b"\x02ACGT"));
        bytes.extend(packet(0x7F, b""));
        bytes.extend(packet(0x42, b"\x00\xff opaque"));
        bytes.extend(packet(0x42, b"second"));
        let file = DnaFile::parse(&bytes).unwrap();
        assert_eq!(file.packets().len(), 5);
        assert_eq!(file.payload(tag::SEQUENCE).unwrap(), b"\x02ACGT");
        assert_eq!(file.packet(0x42).unwrap().payload, b"\x00\xff opaque");
        assert_eq!(file.to_bytes().unwrap(), bytes);
    }

    #[test]
    fn rejects_files_that_are_not_dna() {
        assert_eq!(DnaFile::parse(b""), Err(DnaError::NotDna));
        assert_eq!(
            DnaFile::parse(&packet(tag::SEQUENCE, b"\x02ACGT")),
            Err(DnaError::NotDna)
        );
        assert_eq!(
            DnaFile::parse(&packet(tag::COOKIE, b"NotSnapG")),
            Err(DnaError::NotDna)
        );
    }

    #[test]
    fn reports_truncation_precisely() {
        let mut bytes = cookie();
        bytes.extend(packet(tag::SEQUENCE, b"\x02ACGT"));
        let cut = bytes.len() - 2;
        assert_eq!(
            DnaFile::parse(&bytes[..cut]),
            Err(DnaError::Truncated {
                offset: 19,
                tag: 0,
                declared: 5,
                available: 3
            })
        );
        assert_eq!(
            DnaFile::parse(&bytes[..21]),
            Err(DnaError::TruncatedHeader { offset: 19 })
        );
        // A length field claiming 4 GiB must not allocate anything.
        let mut huge = cookie();
        huge.extend([0x00, 0xff, 0xff, 0xff, 0xff, b'A']);
        assert!(matches!(
            DnaFile::parse(&huge),
            Err(DnaError::Truncated { .. })
        ));
    }

    #[test]
    fn edits_keep_packet_positions() {
        let mut bytes = cookie();
        bytes.extend(packet(tag::SEQUENCE, b"\x02ACGT"));
        bytes.extend(packet(tag::NOTES, b"<Notes/>"));
        let mut file = DnaFile::parse(&bytes).unwrap();
        file.set_payload(tag::SEQUENCE, b"\x03GG".to_vec()).unwrap();
        file.set_payload(tag::FEATURES, b"<Features/>".to_vec())
            .unwrap();
        let tags: Vec<u8> = file.packets().iter().map(|p| p.tag).collect();
        assert_eq!(
            tags,
            vec![tag::COOKIE, tag::SEQUENCE, tag::NOTES, tag::FEATURES]
        );
        assert!(file.set_payload(tag::COOKIE, b"nope".to_vec()).is_err());
        assert_eq!(file.remove(tag::COOKIE), 0);
        assert_eq!(file.remove(tag::NOTES), 1);
        assert_eq!(
            file.payload(tag::NOTES),
            Err(DnaError::MissingPacket { tag: tag::NOTES })
        );
    }
}
