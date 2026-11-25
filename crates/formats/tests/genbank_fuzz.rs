//! Hostile input for the GenBank and FASTA readers: nothing panics, and whatever parses is
//! written back exactly.

use proptest::prelude::*;
use snappiesnap_formats::fasta::FastaFile;
use snappiesnap_formats::genbank::{GenBankFile, Location, format_location};

fn line() -> impl Strategy<Value = String> {
    prop_oneof![
        Just(
            "LOCUS       p                         10 bp    DNA     circular SYN 01-JAN-2020"
                .to_string()
        ),
        Just("FEATURES             Location/Qualifiers".to_string()),
        Just("     CDS             complement(join(1..3,5..9))".to_string()),
        Just("     gene            <1..>10".to_string()),
        Just("                     /label=x y".to_string()),
        Just("                     /note=\"open quote".to_string()),
        Just("                     more\"\"text\"".to_string()),
        Just("                     /pseudo".to_string()),
        Just("ORIGIN".to_string()),
        Just("        1 acgtacgtac".to_string()),
        Just("//".to_string()),
        "[ -~]{0,30}",
    ]
}

proptest! {
    #[test]
    fn arbitrary_genbank_never_panics(lines in proptest::collection::vec(line(), 0..20)) {
        let text = lines.join("\n") + "\n";
        if let Ok(file) = GenBankFile::parse(text.as_bytes()) {
            prop_assert_eq!(file.to_text(), text);
            for r in &file.records {
                let _ = r.locus.to_line();
                for f in &r.features {
                    let _ = f.entry.to_text();
                    if let Ok(loc) = Location::parse(&f.entry.location) {
                        let _ = loc.segments();
                    }
                }
            }
        }
    }

    #[test]
    fn arbitrary_bytes_never_panic(bytes in proptest::collection::vec(any::<u8>(), 0..300)) {
        let _ = GenBankFile::parse(&bytes);
        let _ = FastaFile::parse(&bytes);
    }

    #[test]
    fn arbitrary_locations_never_panic(text in "[0-9.,<>^():a-z]{0,40}") {
        if let Ok(loc) = Location::parse(&text)
            && let Some(segs) = loc.segments()
        {
            // What we write parses back to the same segments.
            let again = Location::parse(&format_location(&segs)).ok().and_then(|l| l.segments());
            prop_assert_eq!(again, Some(segs));
        }
    }
}
