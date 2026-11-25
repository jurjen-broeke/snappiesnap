//! Saving `.dna` documents: unchanged content is byte-exact, edits touch only what changed.

// Test code: the never-crash lints apply to production code only (see clippy.toml).
#![allow(clippy::unwrap_used, clippy::indexing_slicing)]

use snappiesnap_doc::{
    Alphabet, Molecule, NewFeature, Qualifier, Range, Segment, Seq, Strand, Topology,
};
use snappiesnap_formats::dna::{DnaFile, FeatureXml, FeaturesXml, MAGIC, Packet, tag};
use snappiesnap_io::{Format, IoError, Source, export, import};

/// 16 bp circular with an unmodelled packet, two features (one with an unknown child element
/// and attributes we don't model) and notes.
fn sample() -> Vec<u8> {
    let mut cookie = MAGIC.to_vec();
    cookie.extend_from_slice(&[0, 1, 0, 15, 0, 19]);
    let features = concat!(
        r##"<?xml version="1.0"?><Features nextValidID="5">"##,
        r##"<Feature recentID="3" name="kept" directionality="1" translationMW="1.23" type="CDS" allowSegmentOverlaps="0" consecutiveTranslationNumbering="1"><Segment range="2-7" color="#ff0000" type="standard" translated="1"/><Q name="note"><V text="&lt;b&gt;x&lt;/b&gt;"/></Q><Odd/></Feature>"##,
        r##"<Feature recentID="4" name="span" directionality="2" type="misc_feature" allowSegmentOverlaps="0" consecutiveTranslationNumbering="1"><Segment range="15-2" color="#00ff00" type="standard"/></Feature>"##,
        "</Features>\n"
    );
    DnaFile::from_packets(vec![
        Packet { tag: tag::COOKIE, payload: cookie },
        Packet { tag: tag::SEQUENCE, payload: b"\x1fACGTACGTAAAACCCC".to_vec() },
        Packet { tag: 0x03, payload: vec![1, 2, 3, 0xff] },
        Packet { tag: tag::FEATURES, payload: features.as_bytes().to_vec() },
        Packet {
            tag: tag::NOTES,
            payload: b"<Notes>\n<Type>Synthetic</Type>\n<Description>&lt;html>d&lt;/html></Description>\n</Notes>".to_vec(),
        },
    ])
    .unwrap()
    .to_bytes()
    .unwrap()
}

fn open(bytes: &[u8]) -> (Molecule, Source) {
    let i = import("sample.dna", bytes).unwrap();
    assert!(i.warnings.is_empty(), "{:?}", i.warnings);
    (i.molecule, i.source)
}

fn save(m: &Molecule, src: &Source) -> Vec<u8> {
    export(m, Some(src), Format::Dna).unwrap().bytes
}

fn packet(bytes: &[u8], t: u8) -> Vec<u8> {
    DnaFile::parse(bytes).unwrap().payload(t).unwrap().to_vec()
}

fn features_text(bytes: &[u8]) -> String {
    String::from_utf8(packet(bytes, tag::FEATURES)).unwrap()
}

#[test]
fn unchanged_documents_are_byte_exact() {
    let bytes = sample();
    let (m, src) = open(&bytes);
    assert_eq!(save(&m, &src), bytes);
    // Saving again from the new source is still exact.
    let e = export(&m, Some(&src), Format::Dna).unwrap();
    let (out, src2) = (e.bytes, e.source);
    assert_eq!(save(&m, &src2), out);
}

#[test]
fn sequence_edits_rewrite_only_the_sequence_and_moved_segments() {
    let bytes = sample();
    let (mut m, src) = open(&bytes);
    m.insert(10, b"GG").unwrap();
    let out = save(&m, &src);
    let file = DnaFile::parse(&out).unwrap();
    assert_eq!(
        file.payload(tag::SEQUENCE).unwrap(),
        b"\x1fACGTACGTAAGGAACCCC"
    );
    assert_eq!(file.payload(0x03).unwrap(), &[1, 2, 3, 0xff]);
    assert_eq!(
        file.payload(tag::NOTES).unwrap(),
        packet(&bytes, tag::NOTES).as_slice()
    );
    let text = features_text(&out);
    // "kept" didn't move: its text is unchanged, including what we don't model.
    assert!(text.contains(r##"<Feature recentID="3" name="kept" directionality="1" translationMW="1.23" type="CDS" allowSegmentOverlaps="0" consecutiveTranslationNumbering="1"><Segment range="2-7" color="#ff0000" type="standard" translated="1"/><Q name="note"><V text="&lt;b&gt;x&lt;/b&gt;"/></Q><Odd/></Feature>"##));
    // "span" moved by two: only its range changed.
    assert!(
        text.contains(r##"<Segment range="17-2" color="#00ff00" type="standard"/>"##),
        "{text}"
    );
    let (back, _) = open(&out);
    assert_eq!(back.seq(), m.seq());
    assert_eq!(back.features()[1].segments, m.features()[1].segments);
}

#[test]
fn topology_changes_only_the_circular_bit() {
    let bytes = sample();
    let (mut m, src) = open(&bytes);
    m.rotate_origin(14).unwrap();
    m.set_topology(Topology::Linear).unwrap();
    let out = save(&m, &src);
    let seq = packet(&out, tag::SEQUENCE);
    assert_eq!(seq[0], 0x1e);
    let (back, _) = open(&out);
    assert_eq!(back.topology(), Topology::Linear);
    assert_eq!(back.seq(), m.seq());
}

#[test]
fn feature_add_edit_and_remove() {
    let bytes = sample();
    let (mut m, src) = open(&bytes);
    let kept = m.features()[0].clone();
    let span_id = m.features()[1].id;
    m.add_feature(NewFeature {
        color: None,
        name: "start-codon".into(),
        kind: "misc_feature".into(),
        strand: Strand::None,
        segments: vec![Segment {
            range: Range::new(0, 3),
        }],
        qualifiers: vec![],
    })
    .unwrap();
    let mut renamed = kept.clone();
    renamed.name = "renamed & <b>".into();
    renamed.qualifiers = vec![Qualifier {
        name: "note".into(),
        values: vec!["plain".into()],
    }];
    m.update_feature(renamed).unwrap();
    m.remove_feature(span_id).unwrap();
    let out = save(&m, &src);
    let text = features_text(&out);
    assert!(
        text.starts_with(r##"<?xml version="1.0"?><Features nextValidID="6">"##),
        "{text}"
    );
    // The edited feature keeps what we don't model (translationMW, <Odd/>, segment colour).
    assert!(
        text.contains(
            r##"name="renamed &amp; &lt;b&gt;" directionality="1" translationMW="1.23""##
        ),
        "{text}"
    );
    assert!(
        text.contains(r##"<Q name="note"><V text="plain"/></Q><Odd/>"##),
        "{text}"
    );
    assert!(!text.contains(r##"name="span""##));
    assert!(text.ends_with(r##"<Feature recentID="5" name="start-codon" type="misc_feature" allowSegmentOverlaps="0" consecutiveTranslationNumbering="1"><Segment range="1-3" color="#a6acb3" type="standard"/></Feature></Features>
"##), "{text}");
    let (back, back_src) = open(&out);
    assert_eq!(back.features().len(), 2);
    // Saving the re-opened file again is exact.
    assert_eq!(save(&back, &back_src), out);
}

#[test]
fn reverse_multi_segment_features_round_trip() {
    let bytes = sample();
    let (mut m, src) = open(&bytes);
    m.add_feature(NewFeature {
        color: None,
        name: "join".into(),
        kind: "promoter".into(),
        strand: Strand::Reverse,
        segments: vec![
            Segment {
                range: Range::new(9, 12),
            },
            Segment {
                range: Range::new(4, 6),
            },
        ],
        qualifiers: vec![],
    })
    .unwrap();
    let out = save(&m, &src);
    assert!(features_text(&out).contains(r##"directionality="2" type="promoter" allowSegmentOverlaps="0" consecutiveTranslationNumbering="1"><Segment range="5-6" color="#a6acb3" type="standard"/><Segment range="10-12""##));
    let (back, _) = open(&out);
    let join = back.features().iter().find(|f| f.name == "join").unwrap();
    assert_eq!(
        join.segments,
        m.features()
            .iter()
            .find(|f| f.name == "join")
            .unwrap()
            .segments
    );
}

#[test]
fn primers_and_notes_are_saved() {
    let (mut m, src) = open(&sample());
    let before = sample();
    m.add_primer("p", Seq::new(Alphabet::Dna, "ACGTACGTAA").unwrap(), "")
        .unwrap();
    let mut notes = m.notes().clone();
    notes.description = "changed".into();
    notes.fields.insert("CreatedBy".into(), "Lab".into());
    m.set_notes(notes).unwrap();
    let bytes = export(&m, Some(&src), Format::Dna).unwrap().bytes;
    let back = import("sample.dna", &bytes).unwrap();
    assert!(back.warnings.is_empty(), "{:?}", back.warnings);
    let p = &back.molecule.primers()[0];
    assert_eq!(
        (p.name.as_str(), p.seq.as_bytes()),
        ("p", &b"ACGTACGTAA"[..])
    );
    assert_eq!(
        back.molecule.notes().description,
        "<html><body>changed</body></html>"
    );
    assert_eq!(
        back.molecule
            .notes()
            .fields
            .get("CreatedBy")
            .map(String::as_str),
        Some("Lab")
    );
    // Only those packets changed; the unknown 0x03 and the features are as they were.
    for t in [0x03, tag::FEATURES, tag::SEQUENCE] {
        assert_eq!(packet(&bytes, t), packet(&before, t));
    }
    // Saved again unchanged: the same bytes.
    assert_eq!(save(&back.molecule, &back.source), bytes);
    let _ = IoError::UnknownFormat;
}

#[test]
fn new_documents_can_be_saved_and_reopened() {
    let mut m = Molecule::new(
        "pNew",
        Seq::new(Alphabet::Dna, "ACGTACGTAC").unwrap(),
        Topology::Circular,
    )
    .unwrap();
    m.add_feature(NewFeature {
        color: None,
        name: "f".into(),
        kind: "misc_feature".into(),
        strand: Strand::Forward,
        segments: vec![Segment {
            range: Range::new(8, 2),
        }],
        qualifiers: vec![],
    })
    .unwrap();
    let e = export(&m, None, Format::Dna).unwrap();
    let (bytes, src) = (e.bytes, e.source);
    let i = import("pNew.dna", &bytes).unwrap();
    assert_eq!(i.molecule.seq(), m.seq());
    assert_eq!(i.molecule.topology(), Topology::Circular);
    assert_eq!(i.molecule.features()[0].segments, m.features()[0].segments);
    assert_eq!(packet(&bytes, tag::SEQUENCE)[0], 0x03);
    // The returned source makes the next save exact.
    assert_eq!(export(&m, Some(&src), Format::Dna).unwrap().bytes, bytes);
    let empty = FeaturesXml::parse(&features_text(&bytes)).unwrap();
    assert_eq!(empty.next_valid_id(), Some(1));
    let _ = FeatureXml::new(0, None, "", "", None, &[], &[]);
}

#[test]
fn primers_are_shown_and_kept() {
    let primers = concat!(
        r##"<?xml version="1.0"?><Primers nextValidID="2"><HybridizationParams minContinuousMatchLen="10"/>"##,
        r##"<Primer recentID="0" name="fwd &amp; tail" sequence="ggACGTACGTAA" description="&lt;html&gt;&lt;body&gt;&lt;/body&gt;&lt;/html&gt;"><BindingSite location="0-9" boundStrand="0"/></Primer>"##,
        r##"<Primer recentID="1" name="bad" sequence="AC!GT"/>"##,
        "</Primers>\n"
    );
    let mut packets = DnaFile::parse(&sample()).unwrap().packets().to_vec();
    packets.push(Packet {
        tag: tag::PRIMERS,
        payload: primers.as_bytes().to_vec(),
    });
    let bytes = DnaFile::from_packets(packets).unwrap().to_bytes().unwrap();
    let i = import("sample.dna", &bytes).unwrap();
    assert_eq!(i.molecule.primers().len(), 1);
    let p = &i.molecule.primers()[0];
    assert_eq!(
        (p.name.as_str(), p.seq.as_bytes()),
        ("fwd & tail", &b"ggACGTACGTAA"[..])
    );
    assert_eq!(i.warnings.len(), 1, "{:?}", i.warnings);
    assert!(i.warnings[0].contains("bad"));
    // Unchanged, the file is written back exactly, primers packet included.
    assert_eq!(save(&i.molecule, &i.source), bytes);
    // A sequence edit moves the shown primer's binding site; the one we couldn't read is kept.
    let mut m = i.molecule.clone();
    m.insert(0, b"T").unwrap();
    let after = String::from_utf8(packet(&save(&m, &i.source), tag::PRIMERS).to_vec()).unwrap();
    assert!(
        after.contains(r##"<Primer recentID="0" name="fwd &amp; tail" sequence="ggACGTACGTAA" description="&lt;html&gt;&lt;body&gt;&lt;/body&gt;&lt;/html&gt;"><BindingSite location="1-10" boundStrand="0" annealedBases="ACGTACGTAA""##),
        "{after}"
    );
    assert!(after.contains(r##"<Primer recentID="1" name="bad" sequence="AC!GT"/>"##));
}

#[test]
fn feature_colours_carry_into_new_files() {
    // Read from a file, then written into a new one (as when a construct is made from it).
    let (m, _) = open(&sample());
    let colours: Vec<Option<&str>> = m.features().iter().map(|f| f.color.as_deref()).collect();
    assert_eq!(colours, vec![Some("#ff0000"), Some("#00ff00")]);
    let fresh = Molecule::new("fresh", m.seq().clone(), m.topology()).unwrap();
    let mut copy = fresh;
    for f in m.features() {
        copy.add_feature(NewFeature {
            color: f.color.clone(),
            name: f.name.clone(),
            kind: f.kind.clone(),
            strand: f.strand,
            segments: f.segments.clone(),
            qualifiers: Vec::new(),
        })
        .unwrap();
    }
    let out = export(&copy, None, Format::Dna).unwrap();
    let back = import("fresh.dna", &out.bytes).unwrap().molecule;
    let colours: Vec<Option<&str>> = back.features().iter().map(|f| f.color.as_deref()).collect();
    assert_eq!(colours, vec![Some("#ff0000"), Some("#00ff00")]);
    // Not a colour: refused, not written.
    let bad = NewFeature {
        color: Some("red".into()),
        name: "x".into(),
        kind: "misc_feature".into(),
        strand: Strand::None,
        segments: vec![Segment {
            range: Range::new(0, 2),
        }],
        qualifiers: Vec::new(),
    };
    assert!(copy.add_feature(bad).is_err());
}
