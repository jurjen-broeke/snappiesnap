//! Hostile input for the trace readers (ABIF and SCF): arbitrary bytes and damaged copies of
//! valid files never panic, and whatever reads is internally consistent.
// Test code: the never-crash lints apply to production code only (see clippy.toml).
#![allow(clippy::unwrap_used, clippy::indexing_slicing)]

use proptest::prelude::*;
use snappiesnap_formats::trace::Chromatogram;
use snappiesnap_formats::{abif, scf};

fn chromatogram() -> impl Strategy<Value = Chromatogram> {
    (1usize..12, 1usize..60, any::<bool>()).prop_flat_map(|(bases, samples, quality)| {
        (
            proptest::collection::vec(
                prop_oneof![Just(b'A'), Just(b'C'), Just(b'G'), Just(b'T'), Just(b'N')],
                bases,
            ),
            proptest::collection::vec(0u32..samples as u32, bases),
            proptest::collection::vec(0u8..70, bases),
            proptest::collection::vec(proptest::collection::vec(0u16..32_000, samples), 4),
        )
            .prop_map(move |(b, p, q, ch)| Chromatogram {
                sample_name: Some("s".into()),
                bases: b,
                peaks: p,
                quality: quality.then_some(q),
                channels: [ch[0].clone(), ch[1].clone(), ch[2].clone(), ch[3].clone()],
            })
    })
}

fn check(c: &Chromatogram) {
    assert!(c.check().is_ok());
}

proptest! {
    #[test]
    fn arbitrary_bytes_never_panic(bytes in proptest::collection::vec(any::<u8>(), 0..400)) {
        for prefix in [&b"ABIF"[..], b".scf", b""] {
            let mut b = prefix.to_vec();
            b.extend_from_slice(&bytes);
            if let Ok(c) = abif::read_chromatogram(&b) { check(&c) }
            if let Ok(c) = scf::read_chromatogram(&b) { check(&c) }
        }
    }

    #[test]
    fn written_traces_read_back(c in chromatogram()) {
        let a = abif::read_chromatogram(&abif::write_chromatogram(&c).unwrap()).unwrap();
        prop_assert_eq!(&a, &c);
        let s = scf::read_chromatogram(&scf::write_chromatogram(&c).unwrap()).unwrap();
        prop_assert_eq!(&s.bases, &c.bases);
        prop_assert_eq!(&s.channels, &c.channels);
        prop_assert_eq!(&s.peaks, &c.peaks);
    }

    #[test]
    fn damaged_files_never_panic(
        c in chromatogram(),
        hits in proptest::collection::vec((any::<prop::sample::Index>(), any::<u8>()), 1..8),
        cut in any::<prop::sample::Index>(),
    ) {
        for mut bytes in [abif::write_chromatogram(&c).unwrap(), scf::write_chromatogram(&c).unwrap()] {
            for (at, v) in &hits {
                let i = at.index(bytes.len());
                bytes[i] = *v;
            }
            if let Ok(r) = abif::read_chromatogram(&bytes) { check(&r) }
            if let Ok(r) = scf::read_chromatogram(&bytes) { check(&r) }
            let short = &bytes[..cut.index(bytes.len())];
            let _ = abif::read_chromatogram(short);
            let _ = scf::read_chromatogram(short);
        }
    }
}
