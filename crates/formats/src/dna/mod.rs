//! SnapGene `.dna` files.
//!
//! Learned from public third-party descriptions of the format and from files we inspected
//! ourselves (never from the reference application's code). A file is a sequence
//! of packets: a 1-byte tag, a 4-byte big-endian payload length, then the payload. The first
//! packet is the "cookie", which starts with the ASCII magic `SnapGene`.
//!
//! [`DnaFile`] keeps every packet's payload exactly as read, in file order, including tags we
//! don't understand. Writing an unmodified file gives back the identical bytes; an edit replaces
//! only the payloads it touches.

mod container;
mod features;
pub mod history;
mod reads;
mod views;

pub use container::{DnaError, DnaFile, Packet};
pub use features::{Attrs, DEFAULT_COLOR, FeatureChild, FeatureXml, FeaturesXml, SegmentInfo};
pub use reads::{MatchedRead, ReadData};
pub use views::{Cookie, NoteField, SequencePacket};

/// Tags of the packets we know something about.
pub mod tag {
    /// Residues, preceded by one flags byte (bit 0: circular).
    pub const SEQUENCE: u8 = 0x00;
    /// Primers (XML).
    pub const PRIMERS: u8 = 0x05;
    /// Notes: description, dates, author, ... (XML).
    pub const NOTES: u8 = 0x06;
    /// Additional sequence properties: end stickiness and modifications (XML).
    pub const PROPERTIES: u8 = 0x08;
    /// The cookie: magic, sequence type and format versions. Always first.
    pub const COOKIE: u8 = 0x09;
    /// Features (XML).
    pub const FEATURES: u8 = 0x0A;
    /// A sequencing trace matched to the construct (see [`super::MatchedRead`]).
    pub const TRACE: u8 = 0x10;
    /// The list of matched reads (XML, `<AlignableSequences>`).
    pub const ALIGNABLE: u8 = 0x11;
}

/// The magic at the start of the cookie payload.
pub const MAGIC: &[u8; 8] = b"SnapGene";
