//! Hostile input for the GFF3 reader: nothing panics, and whatever parses is written back
//! exactly; rows written from their parsed form read back the same.

use proptest::prelude::*;
use snappiesnap_formats::gff3::{Gff3File, Line};

fn line() -> impl Strategy<Value = String> {
    prop_oneof![
        Just("##sequence-region c 1 40".to_string()),
        Just("# comment".to_string()),
        Just(String::new()),
        Just("c\t.\tregion\t1\t40\t.\t+\t.\tID=c;Is_circular=true".to_string()),
        Just("c\tsrc\tCDS\t35\t50\t.\t-\t0\tID=x;Name=a%3Bb;Note=1,2".to_string()),
        Just("c\t.\tgene\t2\t9\t.\t?\t.\t.".to_string()),
        Just("c\t.\tgene\t2\t9\t.\t+\t.\tID=%zz;=;a=".to_string()),
        Just("##FASTA".to_string()),
        Just(">c desc".to_string()),
        Just("ACGTACGTNN".to_string()),
        "[ -~\t]{0,40}",
    ]
}

proptest! {
    #[test]
    fn arbitrary_gff3_never_panics(lines in proptest::collection::vec(line(), 0..20)) {
        let text = format!("##gff-version 3\n{}\n", lines.join("\n"));
        if let Ok(file) = Gff3File::parse(text.as_bytes()) {
            prop_assert_eq!(file.to_text(), text);
            let _ = file.seqids();
            for l in &file.lines {
                if let Line::Feature { row, .. } = l {
                    let again = Gff3File::parse(format!("##gff-version 3\n{}", row.to_line()).as_bytes());
                    match again.ok().and_then(|f| f.lines.into_iter().nth(1)) {
                        Some(Line::Feature { row: back, .. }) => prop_assert_eq!(&back, row),
                        other => prop_assert!(false, "row didn't read back: {other:?}"),
                    }
                }
            }
        }
    }

    #[test]
    fn arbitrary_bytes_never_panic(bytes in proptest::collection::vec(any::<u8>(), 0..300)) {
        let _ = Gff3File::parse(&bytes);
        let mut with_header = b"##gff-version 3\n".to_vec();
        with_header.extend(bytes);
        let _ = Gff3File::parse(&with_header);
    }
}
