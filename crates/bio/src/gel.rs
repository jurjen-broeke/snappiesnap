//! Agarose gel simulation: where fragments run and how bright they are.
//!
//! Migration is log-linear in fragment length across the gel's separation range, which depends
//! on the agarose concentration (Sambrook J., Russell D. W. (2001) Molecular Cloning, 3rd ed.,
//! Table 5-1: 0.3 % separates 5-60 kb, 0.6 % 1-20 kb, 0.7 % 0.8-10 kb, 0.9 % 0.5-7 kb,
//! 1.2 % 0.4-6 kb, 1.5 % 0.2-3 kb, 2 % 0.1-2 kb; interpolated between those). Fragments
//! outside the range bunch up near the wells or the front. Brightness follows DNA mass (length ×
//! copies), as with intercalating stains. Fragments too close to resolve merge into one band.
//!
//! Supercoiled circular DNA runs faster than linear DNA of the same length by an amount that
//! depends on the gel and the conditions; it is shown at an approximate apparent size of 70 % of
//! its length.

use serde::Serialize;

/// A size ladder: name and band sizes (bp). Generic round-number ladders, not any vendor's.
pub const LADDERS: &[(&str, &[usize])] = &[
    (
        "1 kb",
        &[
            10000, 8000, 6000, 5000, 4000, 3000, 2500, 2000, 1500, 1000, 750, 500, 250,
        ],
    ),
    (
        "100 bp",
        &[
            3000, 2000, 1500, 1200, 1000, 900, 800, 700, 600, 500, 400, 300, 200, 100,
        ],
    ),
    (
        "High range",
        &[
            48000, 30000, 20000, 15000, 10000, 8000, 6000, 5000, 4000, 3000,
        ],
    ),
];

/// Apparent size of supercoiled circular DNA relative to linear DNA of the same length.
pub const SUPERCOILED_FACTOR: f64 = 0.7;

/// Bands closer than this (fraction of the gel's length) look like one.
const RESOLUTION: f64 = 0.008;

/// Agarose concentrations (%) with their separation ranges (bp), from the table above.
const RANGES: [(f64, f64, f64); 7] = [
    (0.3, 5000.0, 60000.0),
    (0.6, 1000.0, 20000.0),
    (0.7, 800.0, 10000.0),
    (0.9, 500.0, 7000.0),
    (1.2, 400.0, 6000.0),
    (1.5, 200.0, 3000.0),
    (2.0, 100.0, 2000.0),
];

/// The separation range (smallest, largest bp) at `agarose` % (clamped to 0.3-2 %).
pub fn separation_range(agarose: f64) -> (f64, f64) {
    let a = if agarose.is_finite() {
        agarose.clamp(0.3, 2.0)
    } else {
        1.0
    };
    let mut prev = RANGES[0];
    for r in RANGES {
        if a == r.0 {
            return (r.1, r.2);
        }
        if a < r.0 {
            let f = if r.0 > prev.0 {
                (a - prev.0) / (r.0 - prev.0)
            } else {
                0.0
            };
            // Interpolate the logarithms of the limits.
            let lo = (prev.1.ln() + f * (r.1.ln() - prev.1.ln())).exp();
            let hi = (prev.2.ln() + f * (r.2.ln() - prev.2.ln())).exp();
            return (lo, hi);
        }
        prev = r;
    }
    (prev.1, prev.2)
}

/// Where a fragment of `bp` ends up: 0 at the wells, 1 at the dye front.
pub fn migration(bp: f64, agarose: f64) -> f64 {
    let (lo, hi) = separation_range(agarose);
    if !(bp.is_finite() && bp > 0.0) {
        return 1.0;
    }
    let (top, bottom) = (0.08, 0.92);
    let x = bp.ln();
    let (llo, lhi) = (lo.ln(), hi.ln());
    if x >= lhi {
        // Larger than the range: crowd towards the wells.
        top * (-(x - lhi)).exp()
    } else if x <= llo {
        // Smaller than the range: crowd towards the front.
        bottom + (1.0 - bottom) * (1.0 - (-(llo - x)).exp())
    } else {
        top + (bottom - top) * (lhi - x) / (lhi - llo)
    }
}

/// A fragment to run: its length, copies, and whether it is supercoiled circular DNA.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct Fragment {
    pub bp: usize,
    pub copies: f64,
    pub supercoiled: bool,
    /// What it is, for the user ("EcoRI – BamHI", "PCR product").
    pub label: String,
}

/// A band: one or more fragments that run together.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct Band {
    /// Position, 0 at the wells to 1 at the front.
    pub y: f64,
    /// Length (bp) of the (largest) fragment in the band.
    pub bp: usize,
    /// DNA mass relative to the brightest band in the lane (0..=1).
    pub intensity: f64,
    /// Fragments in the band.
    pub fragments: Vec<Fragment>,
}

/// Run `fragments` in one lane at `agarose` %.
pub fn lane(fragments: &[Fragment], agarose: f64) -> Vec<Band> {
    let mut placed: Vec<(f64, &Fragment)> = fragments
        .iter()
        .filter(|f| f.bp > 0 && f.copies > 0.0 && f.copies.is_finite())
        .map(|f| {
            let apparent = if f.supercoiled {
                f.bp as f64 * SUPERCOILED_FACTOR
            } else {
                f.bp as f64
            };
            (migration(apparent, agarose), f)
        })
        .collect();
    placed.sort_by(|a, b| a.0.total_cmp(&b.0));
    let mut bands: Vec<Band> = Vec::new();
    for (y, f) in placed {
        let mass = f.bp as f64 * f.copies;
        match bands.last_mut() {
            Some(b) if (y - b.y).abs() < RESOLUTION => {
                b.intensity += mass;
                b.bp = b.bp.max(f.bp);
                b.fragments.push(f.clone());
            }
            _ => bands.push(Band {
                y,
                bp: f.bp,
                intensity: mass,
                fragments: vec![f.clone()],
            }),
        }
    }
    let max = bands.iter().map(|b| b.intensity).fold(0.0, f64::max);
    if max > 0.0 {
        for b in &mut bands {
            b.intensity /= max;
        }
    }
    bands
}

/// A ladder's bands (equal mass per band, as ladders are usually made).
pub fn ladder(name: &str, agarose: f64) -> Option<Vec<Band>> {
    let (_, sizes) = LADDERS.iter().find(|(n, _)| n.eq_ignore_ascii_case(name))?;
    let fragments: Vec<Fragment> = sizes
        .iter()
        .map(|&bp| Fragment {
            bp,
            copies: 1.0 / bp as f64,
            supercoiled: false,
            label: format!("{bp} bp"),
        })
        .collect();
    Some(lane(&fragments, agarose))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn linear(bp: usize) -> Fragment {
        Fragment {
            bp,
            copies: 1.0,
            supercoiled: false,
            label: String::new(),
        }
    }

    #[test]
    fn ranges_follow_the_table() {
        assert_eq!(separation_range(0.7), (800.0, 10000.0));
        assert_eq!(separation_range(2.0), (100.0, 2000.0));
        let (lo, hi) = separation_range(1.0);
        assert!(
            lo > 400.0 && lo < 500.0 && hi > 6000.0 && hi < 7000.0,
            "{lo} {hi}"
        );
        assert_eq!(separation_range(f64::NAN), separation_range(1.0));
        assert_eq!(separation_range(9.0), (100.0, 2000.0));
    }

    #[test]
    fn bigger_runs_slower() {
        for agarose in [0.7, 1.0, 2.0] {
            let ys: Vec<f64> = [50.0, 100.0, 500.0, 1000.0, 5000.0, 20000.0, 100000.0]
                .iter()
                .map(|&bp| migration(bp, agarose))
                .collect();
            assert!(ys.windows(2).all(|w| w[0] > w[1]), "{ys:?}");
            assert!(ys.iter().all(|y| (0.0..=1.0).contains(y)));
        }
        // Log-linear inside the range: 1 kb is halfway (in log) between 500 bp and 2 kb.
        let (a, b, c) = (
            migration(500.0, 1.2),
            migration(1000.0, 1.2),
            migration(2000.0, 1.2),
        );
        assert!(((a - b) - (b - c)).abs() < 1e-9);
    }

    #[test]
    fn bands_merge_and_shine_by_mass() {
        let bands = lane(&[linear(3000), linear(3010), linear(500)], 1.0);
        assert_eq!(bands.len(), 2);
        assert_eq!(bands[0].fragments.len(), 2);
        assert_eq!(bands[0].intensity, 1.0);
        assert!((bands[1].intensity - 500.0 / 6010.0).abs() < 1e-9);
        let sc = Fragment {
            supercoiled: true,
            ..linear(5000)
        };
        let b = lane(&[sc], 1.0);
        assert!((b[0].y - migration(3500.0, 1.0)).abs() < 1e-12);
    }

    #[test]
    fn ladders_and_hostile_input() {
        assert_eq!(ladder("1 KB", 1.0).unwrap().len(), 13);
        assert!(ladder("nope", 1.0).is_none());
        let weird = [
            linear(0),
            Fragment {
                copies: f64::NAN,
                ..linear(10)
            },
            Fragment {
                copies: -1.0,
                ..linear(10)
            },
            linear(usize::MAX),
        ];
        let bands = lane(&weird, f64::INFINITY);
        assert_eq!(bands.len(), 1);
        assert!(bands[0].y.is_finite());
    }
}
