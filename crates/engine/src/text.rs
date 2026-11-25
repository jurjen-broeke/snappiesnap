//! Text files from other programs: lists of enzymes or oligos exported as plain text.

/// Text from bytes: UTF-16 LE or BE when a byte-order mark says so (as some Windows programs
/// export), otherwise UTF-8 (a BOM is dropped). Invalid bytes become U+FFFD, which the parsers
/// then reject where it matters (names, sequences).
pub fn decode_text(bytes: &[u8]) -> String {
    let utf16 = |rest: &[u8], le: bool| {
        let units = rest.as_chunks::<2>().0.iter().map(|&pair| {
            if le {
                u16::from_le_bytes(pair)
            } else {
                u16::from_be_bytes(pair)
            }
        });
        char::decode_utf16(units)
            .map(|c| c.unwrap_or(char::REPLACEMENT_CHARACTER))
            .collect()
    };
    match bytes {
        [0xFF, 0xFE, rest @ ..] => utf16(rest, true),
        [0xFE, 0xFF, rest @ ..] => utf16(rest, false),
        [0xEF, 0xBB, 0xBF, rest @ ..] => String::from_utf8_lossy(rest).into_owned(),
        _ => String::from_utf8_lossy(bytes).into_owned(),
    }
}

#[cfg(test)]
mod tests {
    use super::decode_text;

    #[test]
    fn encodings() {
        let utf16 = |text: &str, le: bool| {
            let mut out = if le {
                vec![0xFF, 0xFE]
            } else {
                vec![0xFE, 0xFF]
            };
            for u in text.encode_utf16() {
                out.extend(if le { u.to_le_bytes() } else { u.to_be_bytes() });
            }
            out
        };
        assert_eq!(decode_text(&utf16("a\tb\r\n", true)), "a\tb\r\n");
        assert_eq!(decode_text(&utf16("Ωx", false)), "Ωx");
        assert_eq!(decode_text(b"\xEF\xBB\xBFabc"), "abc");
        assert_eq!(decode_text(b"abc"), "abc");
        // An odd trailing byte and a lone surrogate never panic.
        assert_eq!(decode_text(&[0xFF, 0xFE, 0x41, 0x00, 0x42]), "A");
        assert_eq!(decode_text(&[0xFF, 0xFE, 0x00, 0xD8]), "\u{FFFD}");
    }
}
