//! Shared test fixtures.

#![allow(dead_code, clippy::unwrap_used)]

use std::sync::Arc;

use snappiesnap_engine::{MemoryFs, Session};
use snappiesnap_formats::dna::{DnaFile, FeatureXml, FeaturesXml, MAGIC, Packet, tag};

/// A small synthetic circular `.dna` file (16 bp) with one feature across the origin.
pub fn example_dna() -> Vec<u8> {
    let mut cookie = MAGIC.to_vec();
    cookie.extend_from_slice(&[0, 1, 0, 15, 0, 19]);
    let mut features = FeaturesXml::empty();
    features.push(FeatureXml::new(
        0,
        None,
        "span",
        "misc_feature",
        Some(1),
        &[(15, 2)],
        &[("note".into(), vec!["synthetic".into()])],
    ));
    features.set_next_valid_id(1);
    DnaFile::from_packets(vec![
        Packet {
            tag: tag::COOKIE,
            payload: cookie,
        },
        Packet {
            tag: tag::SEQUENCE,
            payload: b"\x03ACGTACGTAAAACCCC".to_vec(),
        },
        Packet {
            tag: tag::FEATURES,
            payload: features.to_text().into_bytes(),
        },
    ])
    .unwrap()
    .to_bytes()
    .unwrap()
}

/// A synthetic sequencing read (`.ab1`): 4 poor bases, the first 40 bases of the "template"
/// state's sequence, then 4 poor bases. Base `i` peaks at sample `12 i + 6`.
pub fn read_ab1() -> Vec<u8> {
    use snappiesnap_formats::trace::{CHANNEL_BASES, Chromatogram};
    let core = &b"GATTACAGATTACACCGTTGCATGCAAGCTTGGCACTGGCCGTCG"[..40];
    let bases: Vec<u8> = [&b"NTGA"[..], core, b"CNAN"].concat();
    let n = bases.len();
    let peaks: Vec<u32> = (0..n as u32).map(|i| 12 * i + 6).collect();
    let quality: Vec<u8> = (0..n)
        .map(|i| if i < 4 || i >= n - 4 { 4 } else { 38 })
        .collect();
    let mut channels: [Vec<u16>; 4] = Default::default();
    for c in &mut channels {
        *c = vec![20; 12 * n];
    }
    for (b, p) in bases.iter().zip(&peaks) {
        if let Some(slot) = CHANNEL_BASES.iter().position(|x| x == b) {
            for d in 0..6u32 {
                let v = 900 - 140 * d as u16;
                channels[slot][(p + d) as usize] = v;
                channels[slot][(p - d) as usize] = v;
            }
        }
    }
    snappiesnap_formats::abif::write_chromatogram(&Chromatogram {
        sample_name: Some("read 1 fwd".into()),
        bases,
        peaks,
        quality: Some(quality),
        channels,
    })
    .unwrap()
}

/// A gel experiment (`.sslab`) with two lanes of `example.dna` (found by its relative path),
/// and one of a construct that only exists as the experiment's copy.
pub fn experiment_file() -> Vec<u8> {
    use snappiesnap_engine::experiment::fingerprint;
    let genbank = |name: &str, residues: &str| {
        format!(
            "LOCUS       {name} {} bp    DNA     circular UNK 01-JAN-2026\nORIGIN\n        1 {}\n//\n",
            residues.len(),
            residues.to_lowercase()
        )
    };
    serde_json::to_vec_pretty(&serde_json::json!({
        "format": "snappiesnap-experiment",
        "version": 1,
        "kind": "gel",
        "name": "Fixture run",
        "gel": {"agarose": 1.2, "ladder": "1 kb"},
        "constructs": [
            {"name": "example", "file": "../example.dna", "length": 16, "topology": "circular",
             "fingerprint": fingerprint(b"ACGTACGTAAAACCCC"),
             "genbank": genbank("example", "ACGTACGTAAAACCCC")},
            {"name": "pCopy", "length": 12, "topology": "circular",
             "fingerprint": fingerprint(b"GAATTCAAAAAA"),
             "genbank": genbank("pCopy", "GAATTCAAAAAA")}
        ],
        "lanes": [
            {"construct": 0, "method": {"kind": "uncut"}},
            {"construct": 1, "method": {"kind": "digest", "enzymes": ["EcoRI"]}},
            {"construct": 0, "method": null}
        ]
    }))
    .unwrap()
}

/// A session whose file system is in memory and holds `example.dna`, so tests (including
/// fuzzed `open`/`save` calls) never touch real files.
pub fn session() -> Session {
    let mut s = Session::new().with_platform(Arc::new(
        MemoryFs::new()
            .with_file("example.dna", example_dna())
            .with_file("read.ab1", read_ab1())
            .with_file("runs/experiment.sslab", experiment_file())
            .with_file("experiment.sslab", experiment_file())
            .with_file("withrefm.txt", WITHREFM.as_bytes().to_vec())
            .with_file(
                "insert.fa",
                b">insert\nAAGCTTCCGGAATTCGGCCAAGCTT\n".to_vec(),
            )
            .with_file("dest.gb", genbank("dest", DEST).into_bytes())
            .with_file("entry.gb", genbank("entry", ENTRY).into_bytes())
            .with_file(
                "gibson_a.fa",
                format!(">a\n{}{}\n", &GIBSON_B[GIBSON_B.len() - 20..], GIBSON_A).into_bytes(),
            )
            .with_file(
                "gibson_b.fa",
                format!(">b\n{}{}\n", &GIBSON_A[GIBSON_A.len() - 20..], GIBSON_B).into_bytes(),
            )
            .with_file("gg_vector.gb", genbank("ggvec", GG_VECTOR).into_bytes())
            .with_file("gg_part.fa", format!(">part\n{GG_PART}\n").into_bytes())
            .with_file(
                "mut_template.gb",
                genbank("muttemplate", MUT_TEMPLATE).into_bytes(),
            )
            .with_file("history.dna", history_dna())
            .with_file("matched.dna", matched_dna())
            // GFF3 annotations without a sequence, for any document of 10 bp or more.
            .with_file(
                "annotations.gff3",
                b"##gff-version 3\ndemo\tlab\tpromoter\t2\t6\t.\t+\t.\tID=p1;Name=demo promoter\ndemo\tlab\tCDS\t3\t9\t.\t-\t0\tID=c1;Name=demo cds;gene=x\n".to_vec(),
            )
            // A folder of annotated files for the feature library, and a library already built
            // (its one feature is in history.dna's sequence).
            .with_file("lab/marker.gb", annotated("marker", "TTTTACGTACGTAAAACCCCTTTT", "start block", 5, 20)
                .into_bytes())
            // An oligo list (tab-separated, with a header), and a library already built whose
            // oligos M13F/M13R bind example sequences of the tests.
            .with_file(
                "oligos.txt",
                b"Name\tSequence\nM13F\tGTAAAACGACGGCCAGT\t\nT7\tTAATACGACTCACTATAGGG\tpromoter\n".to_vec(),
            )
            // A folder of cloning primers the way a lab keeps them: one bare sequence per
            // file, named only by the file name. The .dna there is not an oligo list.
            .with_file("primers/SP_201.txt", b"ggatccagatctacgtacgtacgttgcatgcatgcaag".to_vec())
            .with_file(
                "primers/SP_202.txt",
                b"gtagcttgagctcgatacgtacgtacgtacgtacgtacgagag\r\n".to_vec(),
            )
            .with_file("primers/notes.dna", b"not an oligo list".to_vec())
            .with_file(
                "app-data/oligo_library.json",
                br#"{"oligos":[{"name":"M13F","residues":"GTAAAACGACGGCCAGT"},{"name":"tail","residues":"gaattcACGTTTTTACGTCAGT"}]}"#.to_vec(),
            )
            // An exported enzyme set, and one set the user already has.
            .with_file(
                "enzymes.txt",
                b"LAB ENZYMES\r\n\r\nEcoRI, BamHI,\r\nNotI, BlnI\r\n".to_vec(),
            )
            .with_file(
                "app-data/enzyme_sets.json",
                br#"{"sets":[{"name":"Lab freezer","enzymes":["EcoRI","SmaI"]}]}"#.to_vec(),
            )
            .with_file(
                "app-data/feature_library.json",
                br#"{"folder":"lab","files":1,"features":[{"name":"start block","kind":"misc_feature","residues":"ACGTACGTAAAACCCC","directional":true,"source":"marker.gb"}]}"#.to_vec(),
            )
            .with_url(
                snappiesnap_engine::enzymes::REBASE_URL,
                WITHREFM.as_bytes().to_vec(),
            )
            // Two made-up vector records, as feature_library_download receives them.
            .with_url(
                VECTORS_URL,
                format!(
                    "{}\n{}",
                    annotated("vecA", "GGGGCATCATCATCATCATCACGGGG", "His tag", 5, 22),
                    annotated("vecB", "TTTTGACTACAAGGACGACGATGACAAGTTTT", "FLAG tag", 5, 28)
                )
                .into_bytes(),
            ),
    ));
    s.load_enzyme_sets();
    s
}

/// A made-up file in REBASE's withrefm layout (not REBASE data).
pub const WITHREFM: &str = "
REBASE version 999                                              withrefm.999

<1>AaaI
<2>
<3>G^AATTC
<4>
<5>
<6>
<7>XY
<8>

<1>BbbI
<2>
<3>GGTCTC(1/5)
<4>
<5>
<6>
<7>
<8>
";

/// A made-up Gateway destination vector: att1 and att2 cores around a cassette.
pub const DEST: &str =
    "GGGGGGGGGGACAAGTTTGTACAAAAAAGCTGAACCCGGCCGGCCGGACCCAGCTTTCTTGTACAAAGTGGTCCCCCCCCCC";
/// A made-up entry clone with an insert between its cores.
pub const ENTRY: &str =
    "TTTTTTCAACTTTGTACAAAAAAGCAGGCTATGAAAGGGCCCATGAAAACCCAGCTTTCTTGTACAAAGTTGGAAAAAA";

/// A minimal circular GenBank record.
/// A GenBank record with one forward misc_feature at `start..end` (1-based, inclusive).
/// Where the test download of vector records comes from.
pub const VECTORS_URL: &str = "http://fixtures.test/vectors.gb";

pub fn annotated(name: &str, seq: &str, feature: &str, start: usize, end: usize) -> String {
    let plain = genbank(name, seq);
    plain.replacen(
        "ORIGIN\n",
        &format!(
            "FEATURES             Location/Qualifiers\n     misc_feature    {start}..{end}\n                     /label=\"{feature}\"\nORIGIN\n"
        ),
        1,
    )
}

pub fn genbank(name: &str, seq: &str) -> String {
    let mut out = format!(
        "LOCUS       {name:<16} {:>11} bp    DNA     circular SYN 01-JAN-2026\nORIGIN\n",
        seq.len()
    );
    for (i, line) in seq.as_bytes().chunks(60).enumerate() {
        out.push_str(&format!("{:>9}", i * 60 + 1));
        for group in line.chunks(10) {
            out.push(' ');
            out.push_str(&String::from_utf8_lossy(group).to_lowercase());
        }
        out.push('\n');
    }
    out.push_str("//\n");
    out
}

/// Two made-up pieces for homology assembly (each file adds 20 bp of the other at its start).
pub const GIBSON_A: &str = "ATGCGTACGTTAGCCGATCGATCGGCTAGCTAGGCTTACCGATGCAAGT";
pub const GIBSON_B: &str = "TTGACCGGTAACGTTGCAGCTAGCCATGGCATCGATGCATGGCCAATCG";
/// A made-up Golden Gate vector (BsaI sites pointing into a dropout, ends GGAG and CGCT) and a
/// part with matching overhangs.
pub const GG_VECTOR: &str =
    "ACGTTGCAACGTAACCGGTTGGAGAGAGACCTTTTTTTTTTGGTCTCACGCTAACCGGTTAACGTACGTT";
pub const GG_PART: &str = "AAGGTCTCAGGAGATGCCCGGGAAATTTCGCTTGAGACCAA";

/// A small circular plasmid used as a mutagenesis template, from a file.
pub const MUT_TEMPLATE: &str =
    "ACGTTGCAACGTAACCGGTTGGCATGCATTTCGAGCTAGCATCGGATCGATCGTAGCTAGGCATCGATCGACGACTAGCAGCTAG";

fn xz(data: &[u8]) -> Vec<u8> {
    let mut out = Vec::new();
    lzma_rs::xz_compress(&mut std::io::BufReader::new(data), &mut out).unwrap();
    out
}

fn packet(tag: u8, payload: &[u8]) -> Vec<u8> {
    let mut v = vec![tag];
    v.extend((payload.len() as u32).to_be_bytes());
    v.extend_from_slice(payload);
    v
}

/// The construct of `matched.dna`.
pub const MATCHED_CONSTRUCT: &str = "GATTACAGATTACACCGTTGCATGCAAGCTTGGCACTGGCCGTCGTTTTACAACGTCGTGA";

/// A made-up .dna file with one sequencing read matched to its 60 bp linear construct, stored
/// as the reference app does (packet 0x11 lists it, 0x10 holds a raw ZTR trace): bases 10..40
/// of the construct, with base 20 read as A, peaks every 10 samples.
pub fn matched_dna() -> Vec<u8> {
    let mut bases: Vec<u8> = MATCHED_CONSTRUCT.as_bytes()[10..40].to_vec();
    bases[10] = b'A';
    let n = bases.len();
    let mut ztr = vec![0xAE, 0x5A, 0x54, 0x52, 0x0D, 0x0A, 0x1A, 0x0A, 1, 2];
    let mut chunk = |kind: &[u8; 4], data: Vec<u8>| {
        ztr.extend(kind);
        ztr.extend(0u32.to_be_bytes());
        ztr.extend((data.len() as u32).to_be_bytes());
        ztr.extend(data);
    };
    let samples = 10 * n;
    let mut smp = vec![0, 0];
    for channel in b"ACGT" {
        for k in 0..samples {
            let peak = k % 10 == 5 && bases[k / 10] == *channel;
            smp.extend((if peak { 800u16 } else { 20 }).to_be_bytes());
        }
    }
    chunk(b"SMP4", smp);
    let mut base = vec![0];
    base.extend(&bases);
    chunk(b"BASE", base);
    let mut bpos = vec![0, 0, 0, 0];
    for i in 0..n as u32 {
        bpos.extend((10 * i + 5).to_be_bytes());
    }
    chunk(b"BPOS", bpos);
    let mut cnf = vec![0];
    cnf.extend(std::iter::repeat_n(40u8, n));
    cnf.extend(std::iter::repeat_n(0u8, 3 * n));
    chunk(b"CNF4", cnf);
    let mut trace = 0u32.to_be_bytes().to_vec();
    trace.push(0x12);
    trace.extend((ztr.len() as u32).to_be_bytes());
    trace.extend(ztr);
    let list = r#"<AlignableSequences trimStringency="Medium"><Sequence name="read 1" ID="0" isTrace="1" sortOrder="0" trimmedRange="0..29"/></AlignableSequences>"#;
    let mut cookie = MAGIC.to_vec();
    cookie.extend_from_slice(&[0, 1, 0, 15, 0, 20]);
    let mut seq = vec![0u8];
    seq.extend(MATCHED_CONSTRUCT.as_bytes());
    DnaFile::from_packets(vec![
        Packet {
            tag: tag::COOKIE,
            payload: cookie,
        },
        Packet {
            tag: tag::SEQUENCE,
            payload: seq,
        },
        Packet {
            tag: 0x10,
            payload: trace,
        },
        Packet {
            tag: 0x11,
            payload: list.as_bytes().to_vec(),
        },
    ])
    .unwrap()
    .to_bytes()
    .unwrap()
}

/// A made-up .dna file with a construction history: an 18 bp construct ("history") made by
/// inserting TT into a stored 16 bp starting plasmid ("start", with one feature).
pub fn history_dna() -> Vec<u8> {
    let start = b"ACGTACGTAAAACCCC";
    // The stored step: its 2-bit sequence and its packets (a features packet).
    let mut seq = (start.len() as u32).to_be_bytes().to_vec();
    seq.push(0x03);
    seq.extend(1u32.to_be_bytes());
    seq.extend(0u32.to_be_bytes());
    seq.push(1);
    seq.extend((start.len() as u32).to_be_bytes());
    for c in start.chunks(4) {
        let mut byte = 0u8;
        for &b in c {
            byte = (byte << 2)
                | match b {
                    b'G' => 0,
                    b'A' => 1,
                    b'T' => 2,
                    _ => 3,
                };
        }
        seq.push(byte);
    }
    let features = r##"<?xml version="1.0"?><Features nextValidID="1"><Feature recentID="0" name="start feature" directionality="1" type="misc_feature"><Segment range="2-9" color="#a6acb3" type="standard"/></Feature></Features>"##;
    let mut record = 1u32.to_be_bytes().to_vec();
    record.extend(packet(0x01, &seq));
    record.extend(packet(0x1e, &xz(&packet(0x0a, features.as_bytes()))));
    let tree = r#"<?xml version="1.0" encoding="UTF-8"?><HistoryTree><Node name="history.dna" type="DNA" seqLen="18" strandedness="double" ID="3" circular="1" operation="insert"><Node name="start.dna" type="DNA" seqLen="16" strandedness="double" ID="1" circular="1" resurrectable="1" operation="invalid"/></Node></HistoryTree>"#;
    let mut cookie = MAGIC.to_vec();
    cookie.extend_from_slice(&[0, 1, 0, 15, 0, 19]);
    DnaFile::from_packets(vec![
        Packet {
            tag: tag::COOKIE,
            payload: cookie,
        },
        Packet {
            tag: tag::SEQUENCE,
            payload: b"\x03ACGTACGTAAAACCCCTT".to_vec(),
        },
        Packet {
            tag: 0x07,
            payload: xz(tree.as_bytes()),
        },
        Packet {
            tag: 0x0b,
            payload: record,
        },
    ])
    .unwrap()
    .to_bytes()
    .unwrap()
}
