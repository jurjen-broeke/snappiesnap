//! Property tests: every `Range::after_*` rule agrees with a naive model.
//!
//! The model labels each residue with its original position, performs the edit on the list of
//! labels, and computes which new positions hold residues the range originally covered. The
//! range returned by the `after_*` method must cover exactly those positions.

use proptest::prelude::*;
use snappiesnap_seq::{Range, Topology};

/// A molecule length, a topology, and a valid non-empty range on it.
fn molecule_and_range() -> impl Strategy<Value = (usize, Topology, Range)> {
    (1usize..40, any::<bool>())
        .prop_flat_map(|(len, circular)| {
            let topo = if circular {
                Topology::Circular
            } else {
                Topology::Linear
            };
            (Just(len), Just(topo), 0..len, 1..=len)
        })
        .prop_map(|(len, topo, start, span)| {
            let range = if topo == Topology::Circular {
                Range::new(start, (start + span) % len)
            } else {
                let span = span.min(len - start);
                Range::new(start, start + span)
            };
            let range = range.validate(len, topo).unwrap_or(Range::full(len));
            (len, topo, range)
        })
}

fn covered(range: Range, len: usize) -> Vec<usize> {
    (0..len).filter(|&p| range.contains(p, len)).collect()
}

/// New positions whose original label is in `orig`.
fn expected(labels: &[Option<usize>], orig: &[usize]) -> Vec<usize> {
    labels
        .iter()
        .enumerate()
        .filter(|(_, l)| l.is_some_and(|l| orig.contains(&l)))
        .map(|(i, _)| i)
        .collect()
}

proptest! {
    #[test]
    fn delete_matches_model(
        (len, topo, feature) in molecule_and_range(),
        (start, end) in (0usize..40, 0usize..40),
    ) {
        let Ok(del) = Range::new(start % (len + 1), end % (len + 1)).validate(len, topo) else {
            return Ok(());
        };
        prop_assume!(!del.is_empty());
        let orig = covered(feature, len);
        let labels: Vec<Option<usize>> =
            (0..len).filter(|p| !del.contains(*p, len)).map(Some).collect();
        let new_len = labels.len();
        let want = expected(&labels, &orig);
        match feature.after_delete(del, len) {
            None => prop_assert!(want.is_empty(), "lost {want:?}"),
            Some(r) => {
                prop_assert!(r.validate(new_len, topo).is_ok());
                prop_assert_eq!(covered(r, new_len), want);
            }
        }
    }

    #[test]
    fn insert_matches_model(
        (len, topo, feature) in molecule_and_range(),
        at in 0usize..41,
        n in 1usize..5,
    ) {
        let at = at % (len + 1);
        let orig = covered(feature, len);
        let inside = if feature.is_wrapping() {
            at > feature.start || at < feature.end
        } else {
            feature.start < at && at < feature.end
        };
        let mut labels: Vec<Option<usize>> = (0..len).map(Some).collect();
        labels.splice(at..at, std::iter::repeat_n(None, n));
        let mut want = expected(&labels, &orig);
        if inside {
            want.extend(at..at + n);
            want.sort_unstable();
        }
        let r = feature.after_insert(at, n);
        prop_assert!(r.validate(len + n, topo).is_ok());
        prop_assert_eq!(covered(r, len + n), want);
    }

    #[test]
    fn rotate_matches_model((len, _, feature) in molecule_and_range(), origin in 0usize..40) {
        let origin = origin % len;
        let orig = covered(feature, len);
        let mut labels: Vec<Option<usize>> = (0..len).map(Some).collect();
        labels.rotate_left(origin);
        let r = feature.after_rotate(origin, len);
        prop_assert!(r.validate(len, Topology::Circular).is_ok());
        prop_assert_eq!(covered(r, len), expected(&labels, &orig));
    }

    #[test]
    fn reverse_complement_matches_model((len, topo, feature) in molecule_and_range()) {
        let orig = covered(feature, len);
        let labels: Vec<Option<usize>> = (0..len).rev().map(Some).collect();
        let r = feature.after_reverse_complement(len);
        prop_assert!(r.validate(len, topo).is_ok());
        prop_assert_eq!(covered(r, len), expected(&labels, &orig));
    }

    #[test]
    fn validate_never_panics(start in any::<usize>(), end in any::<usize>(), len in 0usize..100) {
        for topo in [Topology::Linear, Topology::Circular] {
            if let Ok(r) = Range::new(start, end).validate(len, topo) {
                prop_assert!(r.start <= len && r.end <= len);
                let _ = r.len_on(len);
            }
        }
    }
}
