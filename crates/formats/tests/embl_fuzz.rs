//! Hostile input for the EMBL reader: nothing panics, and whatever parses is written back
//! exactly.

use proptest::prelude::*;
use snappiesnap_formats::embl::{EmblFile, ft_text, sequence_text, sq_line};

fn line() -> impl Strategy<Value = String> {
    prop_oneof![
        Just("ID   P; SV 1; circular; other DNA; STD; SYN; 10 BP.".to_string()),
        Just("ID   OLD standard; linear DNA; SYN; 5 BP.".to_string()),
        Just("XX".to_string()),
        Just("DE   something".to_string()),
        Just("FH   Key             Location/Qualifiers".to_string()),
        Just("FT   CDS             complement(join(1..3,5..9))".to_string()),
        Just("FT   gene            <1..>10".to_string()),
        Just("FT                   /label=x y".to_string()),
        Just("FT                   /note=\"open quote".to_string()),
        Just("FT                   more\"\"text\"".to_string()),
        Just("SQ   Sequence 10 BP; 3 A; 3 C; 2 G; 2 T; 0 other;".to_string()),
        Just(
            "     acgtacgtac                                                               10"
                .to_string()
        ),
        Just("//".to_string()),
        "[ -~]{0,30}",
    ]
}

proptest! {
    #[test]
    fn arbitrary_embl_never_panics(lines in proptest::collection::vec(line(), 0..20)) {
        let text = lines.join("\n") + "\n";
        if let Ok(file) = EmblFile::parse(text.as_bytes()) {
            prop_assert_eq!(file.to_text(), text);
            for r in &file.records {
                let _ = r.id.to_line();
                let _ = sq_line(&r.residues);
                let _ = sequence_text(&r.residues);
                for f in &r.features {
                    let _ = ft_text(&f.entry);
                }
            }
        }
    }

    #[test]
    fn arbitrary_bytes_never_panic(bytes in proptest::collection::vec(any::<u8>(), 0..300)) {
        let _ = EmblFile::parse(&bytes);
    }
}
