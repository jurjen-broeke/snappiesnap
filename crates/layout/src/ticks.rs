/// A ruler tick at a 0-based position, labelled with its 1-based coordinate.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Tick {
    pub pos: usize,
    pub label: usize,
}

/// A round tick spacing (1, 2 or 5 × 10ⁿ) giving at most about `target` ticks over `len`.
pub fn tick_step(len: usize, target: usize) -> usize {
    let target = target.max(1);
    let raw = len.div_ceil(target).max(1);
    let mut magnitude = 1usize;
    while magnitude.saturating_mul(10) <= raw {
        magnitude = magnitude.saturating_mul(10);
    }
    for m in [1, 2, 5, 10] {
        let step = magnitude.saturating_mul(m);
        if step >= raw {
            return step;
        }
    }
    magnitude.saturating_mul(10)
}

/// Ticks at every multiple of the step, labelled 1-based (position `p` is base `p`), skipping
/// position 0 and any too close to the end.
pub fn ticks(len: usize, target: usize) -> Vec<Tick> {
    if len == 0 {
        return Vec::new();
    }
    let step = tick_step(len, target);
    (1..)
        .map(|k: usize| k.saturating_mul(step))
        .take_while(|&p| p < len)
        .filter(|&p| len - p >= step / 5)
        .map(|p| Tick { pos: p, label: p })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn steps_are_round() {
        assert_eq!(tick_step(9271, 10), 1000);
        assert_eq!(tick_step(711, 10), 100);
        assert_eq!(tick_step(3200, 10), 500);
        assert_eq!(tick_step(1, 10), 1);
        assert_eq!(tick_step(0, 0), 1);
        assert_eq!(tick_step(usize::MAX, 10), tick_step(usize::MAX, 10));
    }

    #[test]
    fn ticks_avoid_the_ends() {
        let t = ticks(9271, 10);
        assert_eq!(t.first().map(|t| t.pos), Some(1000));
        assert_eq!(t.last().map(|t| t.pos), Some(9000));
        assert!(ticks(0, 10).is_empty());
        assert_eq!(ticks(1050, 10).last().map(|t| t.pos), Some(1000));
    }
}
