//! Hostile input for the SBOL reader: arbitrary text and XML-shaped input never panic.

use proptest::prelude::*;
use snappiesnap_formats::sbol::{SbolDocument, parse_tree};

fn piece() -> impl Strategy<Value = String> {
    prop_oneof![
        Just("<rdf:RDF xmlns:rdf=\"http://www.w3.org/1999/02/22-rdf-syntax-ns#\" xmlns:sbol=\"http://sbols.org/v2#\">".to_string()),
        Just("</rdf:RDF>".to_string()),
        Just("<sbol:ComponentDefinition rdf:about=\"http://x/cd\">".to_string()),
        Just("</sbol:ComponentDefinition>".to_string()),
        Just("<sbol:sequenceAnnotation><sbol:SequenceAnnotation rdf:about=\"http://x/a\">".to_string()),
        Just("</sbol:SequenceAnnotation></sbol:sequenceAnnotation>".to_string()),
        Just("<sbol:location><sbol:Range rdf:about=\"http://x/r\"><sbol:start>3</sbol:start><sbol:end>99999999999999999999</sbol:end></sbol:Range></sbol:location>".to_string()),
        Just("<sbol:location><sbol:Cut><sbol:at>-1</sbol:at></sbol:Cut></sbol:location>".to_string()),
        Just("<sbol:sequence rdf:resource=\"http://x/s\"/>".to_string()),
        Just("<sbol:Sequence rdf:about=\"http://x/s\"><sbol:elements>acgt</sbol:elements></sbol:Sequence>".to_string()),
        Just("<![CDATA[x]]>".to_string()),
        Just("<!-- c -->".to_string()),
        Just("&amp;&lt;&#x41;&bogus;".to_string()),
        "[ -~]{0,20}",
    ]
}

proptest! {
    #[test]
    fn xml_shaped_input_never_panics(pieces in proptest::collection::vec(piece(), 0..24)) {
        let text = pieces.concat();
        let _ = parse_tree(&text);
        if let Ok(d) = SbolDocument::parse(text.as_bytes()) {
            let _ = d.main_definition();
            for cd in &d.component_definitions {
                for s in &cd.sequences {
                    let _ = d.sequence(s);
                }
            }
        }
    }

    #[test]
    fn arbitrary_bytes_never_panic(bytes in proptest::collection::vec(any::<u8>(), 0..300)) {
        let _ = SbolDocument::parse(&bytes);
    }
}
