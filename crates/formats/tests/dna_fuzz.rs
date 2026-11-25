//! Hostile input for the `.dna` reader: arbitrary bytes never panic, and anything that parses
//! writes back byte for byte.

use proptest::prelude::*;
use snappiesnap_formats::dna::{DnaFile, FeaturesXml, MAGIC, Packet, tag};

fn with_cookie(rest: &[u8]) -> Vec<u8> {
    let mut v = vec![tag::COOKIE, 0, 0, 0, 14];
    v.extend_from_slice(MAGIC);
    v.extend_from_slice(&[0, 1, 0, 15, 0, 19]);
    v.extend_from_slice(rest);
    v
}

proptest! {
    #[test]
    fn arbitrary_bytes_never_panic(bytes in proptest::collection::vec(any::<u8>(), 0..512)) {
        for input in [bytes.clone(), with_cookie(&bytes)] {
            if let Ok(file) = DnaFile::parse(&input) {
                prop_assert_eq!(file.to_bytes().ok(), Some(input));
            }
        }
    }

    #[test]
    fn any_packet_list_round_trips(
        packets in proptest::collection::vec(
            (any::<u8>(), proptest::collection::vec(any::<u8>(), 0..64)),
            0..12,
        )
    ) {
        let mut list = vec![Packet { tag: tag::COOKIE, payload: with_cookie(&[])[5..].to_vec() }];
        list.extend(packets.into_iter().map(|(tag, payload)| Packet { tag, payload }));
        let Ok(file) = DnaFile::from_packets(list) else {
            return Ok(());
        };
        let bytes = file.to_bytes().unwrap_or_default();
        prop_assert_eq!(DnaFile::parse(&bytes).ok(), Some(file));
    }

    /// Flipping or truncating bytes of a valid file never panics.
    #[test]
    fn corrupted_files_never_panic(
        cut in 0usize..80,
        flips in proptest::collection::vec((0usize..80, any::<u8>()), 0..6),
    ) {
        let mut bytes = with_cookie(&[tag::SEQUENCE, 0, 0, 0, 5, 2, b'A', b'C', b'G', b'T']);
        bytes.extend_from_slice(&[tag::NOTES, 0, 0, 0, 8]);
        bytes.extend_from_slice(b"<Notes/>");
        for (i, b) in flips {
            if let Some(slot) = bytes.get_mut(i) {
                *slot = b;
            }
        }
        bytes.truncate(cut.min(bytes.len()));
        let _ = DnaFile::parse(&bytes);
    }

    /// Arbitrary feature XML never panics, and whatever parses is written back exactly.
    #[test]
    fn arbitrary_feature_xml_never_panics(
        parts in proptest::collection::vec(
            prop_oneof![
                Just("<Features nextValidID=\"1\">".to_string()),
                Just("</Features>".to_string()),
                Just("<Feature name=\"a&amp;b\" type=\"CDS\">".to_string()),
                Just("</Feature>".to_string()),
                Just("<Feature/>".to_string()),
                Just("<Segment range=\"1-3\" color=\"#fff\"/>".to_string()),
                Just("<Q name=\"note\"><V text=\"x\"/></Q>".to_string()),
                Just("<Q name=\"n\">".to_string()),
                Just("</Q>".to_string()),
                Just("<V int=\"1\"/>".to_string()),
                Just("<?xml version=\"1.0\"?>".to_string()),
                ".{0,8}",
            ],
            0..16,
        )
    ) {
        let text = parts.concat();
        if let Ok(parsed) = FeaturesXml::parse(&text) {
            prop_assert_eq!(parsed.to_text(), text.clone());
            for f in parsed.features() {
                let _ = f.segments();
                let _ = f.qualifiers();
                let _ = f.to_xml();
            }
        }
    }
}
