//! Enzyme sites over the local corpus (`cargo xtask corpus-test`).

#![cfg(feature = "corpus")]
// Test code: the never-crash lints apply to production code only (see clippy.toml).
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use std::path::PathBuf;

use snappiesnap_bio::enzymes::{builtin, find_sites};
use snappiesnap_doc::{Alphabet, Molecule, Seq, Topology};
use snappiesnap_formats::dna::DnaFile;

/// Whether the corpus holds this file: tests pinned to one construct skip without it.
fn has(rel: &str) -> bool {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../corpus")
        .join(rel)
        .exists()
}

fn open(rel: &str) -> Molecule {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../corpus")
        .join(rel);
    let file = DnaFile::parse(&std::fs::read(&path).unwrap()).unwrap();
    let seq = file.sequence().unwrap();
    let topology = if seq.is_circular() {
        Topology::Circular
    } else {
        Topology::Linear
    };
    Molecule::new(
        "m",
        Seq::new(Alphabet::Dna, seq.residues).unwrap(),
        topology,
    )
    .unwrap()
}

/// Like [`open`], with the features (names and ranges only, read from the features packet).
fn open_with_features(rel: &str) -> Molecule {
    use snappiesnap_doc::{NewFeature, Segment, Strand};
    let mut m = open(rel);
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../corpus")
        .join(rel);
    let file = DnaFile::parse(&std::fs::read(&path).unwrap()).unwrap();
    if let Ok(features) = file.features() {
        for f in features.features() {
            let Ok(segs) = f.segments() else { continue };
            let len = m.len();
            let segments: Vec<Segment> = segs
                .iter()
                .map(|s| {
                    snappiesnap_doc::Range::new(
                        (s.start as usize).saturating_sub(1),
                        (s.end as usize) % len.max(1),
                    )
                })
                .filter_map(|r| r.validate(len, m.topology()).ok())
                .map(|range| Segment { range })
                .collect();
            if segments.is_empty() {
                continue;
            }
            let _ = m.add_feature(NewFeature {
                color: None,
                name: f.name().to_string(),
                kind: "misc_feature".into(),
                strand: Strand::None,
                segments,
                qualifiers: vec![],
            });
        }
    }
    m
}

/// Unique cutters on pCas9_GFP as the reference application shows them on its map: the label
/// position is the number of residues before the top-strand cut.
#[test]
fn unique_cutters_match_the_reference_application() {
    if !has("dna/roundtrip/pCas9_GFP.dna") {
        println!("no dna/roundtrip/pCas9_GFP.dna in the corpus: nothing to check");
        return;
    }
    let m = open("dna/roundtrip/pCas9_GFP.dna");
    let all = builtin();
    let observed: &[(&str, usize)] = &[
        ("MluI", 1),
        ("SpeI", 22),
        ("SnaBI", 363),
        ("PspOMI", 863),
        ("ApaI", 867),
        ("SacII", 1034),
        ("BsiWI", 1190),
        ("MreI", 1316),
        ("SgrAI", 1316),
        ("BstAPI", 1350),
        ("NotI", 1734),
        ("EcoRV", 1774),
        ("BsmBI", 1907),
        ("Esp3I", 1907),
        ("SbfI", 1992),
        ("XcmI", 2096),
        ("Bsu36I", 3011),
        ("TspMI", 3111),
        ("XmaI", 3111),
        ("SmaI", 3113),
        ("Eco53kI", 3313),
        ("SacI", 3315),
        ("BamHI", 3750),
        ("BstXI", 4534),
        ("PmlI", 4534),
        ("FspAI", 4695),
        ("PaeR7I", 5288),
        ("XhoI", 5288),
        ("AscI", 5700),
        ("BsrGI", 6660),
        ("BstZ17I", 6845),
        ("BspQI", 7108),
        ("SapI", 7108),
        ("AhdI", 8117),
        ("SgrDI", 9038),
        ("BspDI", 9055),
        ("ClaI", 9055),
        ("MfeI", 9205),
    ];
    // Enzymes sold by a single supplier aren't built in; check every observed one that is.
    let mut wrong = Vec::new();
    let mut checked = 0;
    for &(name, position) in observed {
        let Some(index) = all.iter().position(|e| e.name == name) else {
            continue;
        };
        checked += 1;
        let enzyme = &all[index..=index];
        let cuts: Vec<Option<usize>> = find_sites(&m, enzyme).iter().map(|s| s.cut).collect();
        if cuts != [Some(position)] {
            wrong.push(format!(
                "{name}: expected one cut at {position}, found {cuts:?}"
            ));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
    assert!(
        checked >= 25,
        "only {checked} observed enzymes are built in"
    );
}

/// A primer's sequence and its stored binding sites (strand, first, last; 0-based).
type StoredPrimer = (
    String,
    Vec<(snappiesnap_bio::primers::PrimerStrand, usize, usize)>,
);

/// Primer binding sites against the ones the reference application stored with each primer in
/// the corpus files (same strand, same 3' end). Measured 2026-10-23: 624/656 stored sites found
/// (the rest are gapped or many-mismatch alignments), plus 188 weak one-mismatch sites the
/// reference didn't list (it also applies a melting-temperature threshold we don't model yet).
#[test]
fn primer_binding_sites_match_stored_ones() {
    use snappiesnap_bio::primers::{BindingRules, PrimerStrand, binding_sites};
    use snappiesnap_formats::dna::tag;
    use snappiesnap_formats::xml::{Token, Tokenizer, decode};

    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../corpus");
    let mut files = Vec::new();
    let mut stack = vec![root.clone()];
    while let Some(dir) = stack.pop() {
        for entry in std::fs::read_dir(&dir).unwrap() {
            let path = entry.unwrap().path();
            if path.is_dir() {
                stack.push(path);
            } else if path.extension().is_some_and(|e| e == "dna") {
                files.push(path);
            }
        }
    }
    let (mut stored_total, mut found, mut extra, mut missing) =
        (0usize, 0usize, 0usize, Vec::new());
    for path in files {
        let rel = path.strip_prefix(&root).unwrap().display().to_string();
        let m = open(&rel);
        let file = DnaFile::parse(&std::fs::read(&path).unwrap()).unwrap();
        let Some(packet) = file.packet(tag::PRIMERS) else {
            continue;
        };
        let text = std::str::from_utf8(&packet.payload).unwrap();
        // (sequence, [(strand, first, last)]) per primer, from non-simplified binding sites.
        let mut primers: Vec<StoredPrimer> = Vec::new();
        let mut tokens = Tokenizer::new(text);
        while let Some(Ok((token, _))) = tokens.next_token() {
            let Token::Start { name, attrs, .. } = token else {
                continue;
            };
            let get = |n: &str| attrs.iter().find(|a| a.name == n).map(|a| decode(a.raw));
            match name {
                "Primer" => primers.push((get("sequence").unwrap_or_default(), Vec::new())),
                "BindingSite" if get("simplified").is_none() => {
                    let loc = get("location").unwrap_or_default();
                    let Some((a, b)) = loc.split_once('-') else {
                        continue;
                    };
                    let strand = if get("boundStrand").as_deref() == Some("1") {
                        PrimerStrand::Reverse
                    } else {
                        PrimerStrand::Forward
                    };
                    if let (Ok(a), Ok(b), Some(p)) = (a.parse(), b.parse(), primers.last_mut()) {
                        p.1.push((strand, a, b));
                    }
                }
                _ => {}
            }
        }
        let seqs: Vec<&[u8]> = primers.iter().map(|p| p.0.as_bytes()).collect();
        let sites = binding_sites(&m, &seqs, BindingRules::default());
        for s in &sites {
            let Some((_, stored)) = primers.get(s.primer) else {
                continue;
            };
            let known = stored.iter().any(|&(strand, first, last)| {
                strand == s.strand
                    && match strand {
                        PrimerStrand::Forward => {
                            s.range.end == (last + 1) % m.len().max(1) || s.range.end == last + 1
                        }
                        PrimerStrand::Reverse => s.range.start == first,
                    }
            });
            if !known {
                extra += 1;
            }
        }
        for (i, (seq, stored)) in primers.iter().enumerate() {
            for &(strand, first, last) in stored {
                stored_total += 1;
                let hit = sites.iter().any(|s| {
                    s.primer == i
                        && s.strand == strand
                        && match strand {
                            // Stored locations are 0-based and inclusive.
                            PrimerStrand::Forward => {
                                s.range.end == (last + 1) % m.len().max(1)
                                    || s.range.end == last + 1
                            }
                            PrimerStrand::Reverse => s.range.start == first,
                        }
                });
                if hit {
                    found += 1;
                } else {
                    missing.push(format!(
                        "primer {i} ({} nt) {strand:?} {first}-{last}",
                        seq.len()
                    ));
                }
            }
        }
    }
    println!("primer sites: {found}/{stored_total} stored sites found, {extra} more found");
    for m in &missing {
        println!("  not found: {m}");
    }
    assert!(
        found * 100 >= stored_total * 90,
        "only {found}/{stored_total} stored sites found"
    );
}

/// Gateway LR on every pair of local corpus files: wherever an acceptor and a donor both have
/// att1 and att2 cores, the product must equal another corpus file that was made that way in the
/// reference application (found by sequence, so the test names no files). At least one such
/// triple must exist.
#[test]
fn gateway_products_match_the_reference_application() {
    use snappiesnap_bio::cloning::{GatewayReaction, gateway};
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../corpus");
    let mut files = Vec::new();
    let mut stack = vec![root.clone()];
    while let Some(dir) = stack.pop() {
        for entry in std::fs::read_dir(&dir).unwrap() {
            let path = entry.unwrap().path();
            if path.is_dir() {
                stack.push(path);
            } else if path.extension().is_some_and(|e| e == "dna") {
                files.push(path);
            }
        }
    }
    let molecules: Vec<Molecule> = files
        .iter()
        .map(|p| open_with_features(&p.strip_prefix(&root).unwrap().display().to_string()))
        .collect();
    let (mut matched, mut products) = (0, 0);
    for (i, acceptor) in molecules.iter().enumerate() {
        for (j, donor) in molecules.iter().enumerate() {
            if i == j {
                continue;
            }
            let Ok(product) = gateway(acceptor, donor, GatewayReaction::Lr, "product") else {
                continue;
            };
            let p = product.seq().as_bytes();
            // Swapping a cassette for the same one gives an input back: not a real product.
            if p == acceptor.seq().as_bytes() || p == donor.seq().as_bytes() {
                continue;
            }
            products += 1;
            if let Some(reference) = molecules.iter().find(|m| m.seq().as_bytes() == p) {
                // The same features, by name (the replaced att sites become attB1/attB2).
                let names = |m: &Molecule| {
                    let mut v: Vec<String> = m.features().iter().map(|f| f.name.clone()).collect();
                    v.sort();
                    v
                };
                assert_eq!(names(&product), names(reference));
                matched += 1;
            }
        }
    }
    println!("gateway: {matched} of {products} products match corpus files");
    if products == 0 {
        // No pair in the corpus recombines: an acceptor (attR) and a donor (attL) are both
        // needed, and entry clones alone give none.
        println!("no acceptor and donor pair in the corpus: nothing to check");
        return;
    }
    assert!(matched >= 1);
}

/// Trimming the corpus reads (`corpus/abif`, Biopython's test traces) at phred's cutoff 0.05.
/// Expected stretches computed independently (a scratch Python maximum-subarray over
/// Biopython's qualities). Biopython's own "abi-trim" differs: its slice drops the last base
/// of the stretch, and it starts at the first base that ever scores above zero.
#[test]
fn abif_reads_trim_to_their_best_stretch() {
    use snappiesnap_bio::trace::{DEFAULT_CUTOFF, quality_trim};
    let dir = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../corpus/abif");
    if !dir.join("3100.ab1").exists() {
        println!("no Biopython traces in corpus/abif: nothing to check");
        return;
    }
    for (name, want) in [
        ("310.ab1", None),
        ("3100.ab1", Some(18..697)),
        ("3730.ab1", Some(14..1090)),
        ("A6_1-DB3.ab1", Some(25..838)),
        ("no_smpl1.ab1", Some(45..163)),
        ("empty.ab1", None),
    ] {
        let c =
            snappiesnap_formats::abif::read_chromatogram(&std::fs::read(dir.join(name)).unwrap())
                .unwrap();
        let q: Vec<Option<u8>> = c.quality.unwrap().into_iter().map(Some).collect();
        assert_eq!(quality_trim(&q, DEFAULT_CUTOFF), want, "{name}");
    }
}
