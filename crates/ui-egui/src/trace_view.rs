//! The Trace tab: the chromatogram of a sequencing read. The four channels are drawn as curves,
//! each base call is written above its peak with a bar for its quality, and the ends that
//! trimming would remove are shaded. Clicking or dragging over the bases selects them, so the
//! selection is shared with the sequence view and edits.

use eframe::egui::{
    self, Align2, Color32, FontId, Pos2, Rect, RichText, Sense, Shape, Stroke, pos2, vec2,
};
use snappiesnap_engine::bio::trace::{DEFAULT_CUTOFF, quality_trim, read_quality};
use snappiesnap_engine::doc::{Channel, Molecule, Range, Trace};

use crate::state::{clamp_trace_gain, clamp_trace_zoom};
use crate::theme::Tokens;

/// Space left and right of the signal, in points.
const MARGIN: f32 = 24.0;
const RULER: f32 = 16.0;
const LETTERS: f32 = 18.0;
const QUALITY: f32 = 26.0;
/// Quality drawn as a full bar (Phred 60: one error in a million).
const QUALITY_FULL: f32 = 60.0;
/// Quality below this is drawn muted (Phred 20: one error in 100).
const GOOD_QUALITY: u8 = 20;
/// Most sample points per channel drawn per frame; denser views draw each column's highest.
const MAX_POINTS: usize = 8_000;

/// What the user asked for.
#[derive(Debug, Clone, PartialEq)]
pub enum TraceAction {
    Select(Range),
    Trim,
}

/// Drag state between frames.
#[derive(Debug, Default)]
pub struct TraceViewState {
    anchor: Option<usize>,
}

/// Where each residue is drawn, in sample points: at its peak, or spread evenly between the
/// called neighbours for residues typed in by the user.
pub(crate) fn positions(t: &Trace, n: usize) -> Vec<f32> {
    let peaks: Vec<Option<f32>> = (0..n).map(|i| t.peak(i).map(|p| p as f32)).collect();
    let spacing = {
        let known: Vec<f32> = peaks.iter().flatten().copied().collect();
        match (known.first(), known.last()) {
            (Some(a), Some(b)) if known.len() > 1 => ((b - a) / (known.len() - 1) as f32).max(1.0),
            _ => 12.0,
        }
    };
    let mut out = vec![0.0; n];
    let mut i = 0;
    while i < n {
        if let Some(p) = peaks.get(i).copied().flatten() {
            if let Some(o) = out.get_mut(i) {
                *o = p;
            }
            i += 1;
            continue;
        }
        // A run of residues without peaks, between the called residues around it.
        let end = (i..n)
            .find(|&j| peaks.get(j).copied().flatten().is_some())
            .unwrap_or(n);
        let before = i.checked_sub(1).and_then(|j| out.get(j).copied());
        let after = peaks.get(end).copied().flatten();
        let run = (end - i) as f32;
        for (k, j) in (i..end).enumerate() {
            let k = k as f32 + 1.0;
            let x = match (before, after) {
                (Some(a), Some(b)) => a + (b - a) * k / (run + 1.0),
                (Some(a), None) => a + spacing * k,
                (None, Some(b)) => b - spacing * (run + 1.0 - k),
                (None, None) => spacing * (k - 0.5),
            };
            if let Some(o) = out.get_mut(j) {
                *o = x;
            }
        }
        i = end;
    }
    out
}

/// The highest signal over the residues in `keep` (the whole read when nothing is kept), at
/// least 1.
fn scale_peak(t: &Trace, keep: Option<std::ops::Range<usize>>, n: usize) -> u16 {
    let keep = keep.unwrap_or(0..n);
    let peaks: Vec<usize> = keep.filter_map(|i| t.peak(i)).collect();
    let (Some(&lo), Some(&hi)) = (peaks.iter().min(), peaks.iter().max()) else {
        return t.max_value();
    };
    (lo..=hi.min(t.samples().saturating_sub(1)))
        .flat_map(|i| Channel::ALL.into_iter().filter_map(move |c| t.value(c, i)))
        .max()
        .unwrap_or(1)
        .max(1)
}

/// The toolbar, then the chromatogram. Returns what the user asked for.
pub fn show(
    ui: &mut egui::Ui,
    t: &Tokens,
    m: &Molecule,
    selection: Option<Range>,
    zoom: &mut f32,
    gain: &mut f32,
    st: &mut TraceViewState,
) -> Vec<TraceAction> {
    let mut actions = Vec::new();
    let Some(trace) = m.trace() else {
        ui.label(
            RichText::new("This document has no chromatogram. Open a read from an .ab1 or .scf file to see one.")
                .color(t.text_muted),
        );
        return actions;
    };
    let n = m.len();
    let quality: Vec<Option<u8>> = trace.calls().iter().map(|c| c.quality).collect();
    let keep = quality_trim(&quality, DEFAULT_CUTOFF);
    let q = read_quality(&quality);
    ui.horizontal(|ui| {
        let name = trace.sample_name().unwrap_or(m.name());
        let mut parts = vec![name.to_string(), format!("{n} bases")];
        if q.scored > 0 {
            parts.push(format!("mean quality {:.0}", q.mean));
            parts.push(format!("{} bases \u{2265} Q20", q.q20));
        }
        if trace.is_reversed() {
            parts.push("reverse complement".into());
        }
        ui.label(RichText::new(parts.join("  \u{b7}  ")).color(t.text_muted));
        ui.separator();
        let trimmable = keep.as_ref().is_some_and(|k| k.start > 0 || k.end < n);
        let trim = ui
            .add_enabled(trimmable, egui::Button::new("Trim ends").small())
            .on_hover_text(match &keep {
                Some(k) => format!(
                    "Delete the shaded low-quality ends: {} bases at the start, {} at the end (undoable)",
                    k.start,
                    n.saturating_sub(k.end)
                ),
                None => "No stretch of this read is of good quality".into(),
            });
        if trim.clicked() {
            actions.push(TraceAction::Trim);
        }
        ui.separator();
        ui.label(RichText::new("Width").color(t.text_muted));
        if ui.small_button("\u{2212}").on_hover_text("Narrower").clicked() {
            *zoom = clamp_trace_zoom(*zoom / 1.5);
        }
        if ui.small_button("+").on_hover_text("Wider").clicked() {
            *zoom = clamp_trace_zoom(*zoom * 1.5);
        }
        ui.label(RichText::new("Height").color(t.text_muted));
        if ui.small_button("\u{2212}").on_hover_text("Lower peaks").clicked() {
            *gain = clamp_trace_gain(*gain / 1.5);
        }
        if ui.small_button("+").on_hover_text("Higher peaks").clicked() {
            *gain = clamp_trace_gain(*gain * 1.5);
        }
        ui.separator();
        for (c, color) in Channel::ALL.iter().zip(t.trace_channels) {
            ui.label(RichText::new(char::from(c.base()).to_string()).color(color).strong());
        }
    });

    let xs = positions(trace, n);
    let samples = trace.samples();
    let last_x = xs.last().copied().unwrap_or(0.0).max(samples as f32);
    let width = MARGIN * 2.0 + last_x * *zoom;
    let height = ui.available_height().max(160.0);
    // Scaled to the highest peak of the stretch trimming keeps: the poor start of a read often
    // has dye blobs far above the real peaks.
    let max = f32::from(scale_peak(trace, keep.clone(), n)) / *gain;
    egui::ScrollArea::horizontal()
        .id_salt("trace-view")
        .auto_shrink([false, false])
        .show_viewport(ui, |ui, viewport| {
            let (rect, response) =
                ui.allocate_exact_size(vec2(width, height), Sense::click_and_drag());
            let painter = ui.painter_at(rect);
            painter.rect_filled(rect, 0.0, t.canvas);
            let x_of = |s: f32| rect.left() + MARGIN + s * *zoom;
            let visible = Rect::from_min_max(
                pos2(rect.left() + viewport.left(), rect.top()),
                pos2(rect.left() + viewport.right(), rect.bottom()),
            );
            // Residue i covers the half-way points to its neighbours.
            let edge = |i: usize| -> f32 {
                let here = xs.get(i).copied();
                let before = i.checked_sub(1).and_then(|j| xs.get(j).copied());
                match (before, here) {
                    (Some(a), Some(b)) => x_of((a + b) / 2.0),
                    (None, Some(b)) => x_of(b) - 6.0 * *zoom,
                    (Some(a), None) => x_of(a) + 6.0 * *zoom,
                    (None, None) => x_of(0.0),
                }
            };
            let top = rect.top();
            let signal_top = top + RULER + LETTERS + QUALITY + 6.0;
            let signal_bottom = rect.bottom() - 8.0;
            let span = (signal_bottom - signal_top).max(10.0);

            // Ends that trimming would remove.
            let shade = Color32::from_rgba_unmultiplied(
                t.text_muted.r(),
                t.text_muted.g(),
                t.text_muted.b(),
                0x26,
            );
            let (cut_start, cut_end) = match &keep {
                Some(k) => (k.start, k.end),
                None => (n, n),
            };
            if cut_start > 0 {
                let r = Rect::from_min_max(
                    pos2(rect.left(), top),
                    pos2(edge(cut_start), rect.bottom()),
                );
                painter.rect_filled(r, 0.0, shade);
            }
            if cut_end < n {
                let r =
                    Rect::from_min_max(pos2(edge(cut_end), top), pos2(rect.right(), rect.bottom()));
                painter.rect_filled(r, 0.0, shade);
            }

            // The selection.
            if let Some(sel) = selection {
                let (a, b) = sel.segments(n);
                for piece in std::iter::once(a).chain(b) {
                    if piece.is_empty() {
                        let x = edge(piece.start);
                        painter.line_segment(
                            [pos2(x, top), pos2(x, rect.bottom())],
                            Stroke::new(1.5, t.accent),
                        );
                    } else {
                        let r = Rect::from_min_max(
                            pos2(edge(piece.start), top),
                            pos2(edge(piece.end), rect.bottom()),
                        );
                        painter.rect_filled(r, 0.0, t.selection);
                    }
                }
            }

            // Signal: only the visible samples, at most MAX_POINTS per channel.
            let first =
                (((visible.left() - x_of(0.0)) / *zoom).floor().max(0.0) as usize).min(samples);
            let last = ((((visible.right() - x_of(0.0)) / *zoom).ceil().max(0.0) as usize) + 1)
                .min(samples);
            let step = (last.saturating_sub(first) / MAX_POINTS).max(1);
            for (c, color) in Channel::ALL.iter().zip(t.trace_channels) {
                let mut points = Vec::with_capacity(last.saturating_sub(first) / step + 1);
                let mut i = first;
                while i < last {
                    let end = i.saturating_add(step).min(last);
                    let v = (i..end)
                        .filter_map(|k| trace.value(*c, k))
                        .max()
                        .unwrap_or(0);
                    let y = signal_bottom - (f32::from(v) / max).min(1.0) * span;
                    points.push(pos2(x_of(i as f32), y));
                    i = end;
                }
                if points.len() > 1 {
                    painter.add(Shape::line(points, Stroke::new(1.2, color)));
                }
            }

            // Ruler, letters and quality bars for the visible residues.
            let font = FontId::monospace(13.0);
            let small = FontId::proportional(10.0);
            let residues = m.seq().as_bytes();
            let label_every = if *zoom * 12.0 < 9.0 { 50 } else { 10 };
            for (i, &x) in xs.iter().enumerate() {
                let px = x_of(x);
                if px < visible.left() - 20.0 || px > visible.right() + 20.0 {
                    continue;
                }
                let base = residues.get(i).copied().unwrap_or(b'N');
                let call = trace.calls().get(i).copied().unwrap_or_default();
                if (i + 1) % label_every == 0 {
                    painter.text(
                        pos2(px, top + 2.0),
                        Align2::CENTER_TOP,
                        (i + 1).to_string(),
                        small.clone(),
                        t.text_muted,
                    );
                }
                if *zoom * 12.0 >= 7.0 {
                    let color = match base.to_ascii_uppercase() {
                        b'A' => t.trace_channels[0],
                        b'C' => t.trace_channels[1],
                        b'G' => t.trace_channels[2],
                        b'T' => t.trace_channels[3],
                        _ => t.text_muted,
                    };
                    let color = if call.quality.is_some_and(|q| q < GOOD_QUALITY) {
                        color.gamma_multiply(0.55)
                    } else {
                        color
                    };
                    painter.text(
                        pos2(px, top + RULER + LETTERS / 2.0),
                        Align2::CENTER_CENTER,
                        char::from(base).to_string(),
                        font.clone(),
                        color,
                    );
                }
                if let Some(q) = call.quality {
                    let h = (f32::from(q) / QUALITY_FULL).min(1.0) * (QUALITY - 4.0);
                    let w = (*zoom * 12.0 * 0.7).clamp(1.0, 10.0);
                    let bottom = top + RULER + LETTERS + QUALITY - 2.0;
                    let fill = if q >= GOOD_QUALITY {
                        t.accent.gamma_multiply(0.6)
                    } else {
                        t.error.gamma_multiply(0.6)
                    };
                    painter.rect_filled(
                        Rect::from_min_max(
                            pos2(px - w / 2.0, bottom - h),
                            pos2(px + w / 2.0, bottom),
                        ),
                        0.0,
                        fill,
                    );
                }
                if call.peak.is_none() {
                    // Typed in, not called: marked with a dot under the letter.
                    painter.circle_filled(pos2(px, top + RULER + LETTERS + 1.0), 1.5, t.accent);
                }
            }
            painter.line_segment(
                [
                    pos2(rect.left(), signal_bottom),
                    pos2(rect.right(), signal_bottom),
                ],
                Stroke::new(1.0, t.border),
            );

            // Click or drag over the bases to select them.
            // The residue whose position is nearest the pointer.
            let residue_at = |p: Pos2| -> usize {
                let s = (p.x - x_of(0.0)) / *zoom;
                let j = xs.partition_point(|&x| x < s);
                let nearer_before = j
                    .checked_sub(1)
                    .and_then(|i| Some((xs.get(i)?, xs.get(j))))
                    .is_some_and(|(a, b)| b.is_none_or(|b| s - a < b - s));
                if nearer_before {
                    j - 1
                } else {
                    j.min(n.saturating_sub(1))
                }
            };
            if n > 0
                && let Some(p) = response.interact_pointer_pos()
            {
                let i = residue_at(p);
                if response.drag_started() || response.clicked() {
                    st.anchor = Some(i);
                }
                if let Some(a) = st.anchor {
                    let (lo, hi) = (a.min(i), a.max(i));
                    let range = Range::new(lo, hi.saturating_add(1).min(n));
                    if selection != Some(range) {
                        actions.push(TraceAction::Select(range));
                    }
                }
            }
            if response.drag_stopped() || response.clicked() {
                st.anchor = None;
            }
        });
    actions
}

/// A snippet of a read's chromatogram: `half` bases either side of read base `center`, in a
/// `size` box. With `reverse` (the read lies on the construct's other strand) it is drawn as
/// the construct's top strand reads it: mirrored, with complementary bases and channels. The
/// centre base is highlighted.
pub fn snippet(
    ui: &mut egui::Ui,
    t: &Tokens,
    m: &Molecule,
    center: usize,
    half: usize,
    reverse: bool,
    size: egui::Vec2,
) {
    let (rect, _) = ui.allocate_exact_size(size, Sense::hover());
    let painter = ui.painter_at(rect);
    painter.rect_filled(rect, 2.0, t.canvas);
    let Some(trace) = m.trace() else {
        painter.text(
            rect.center(),
            Align2::CENTER_CENTER,
            "no trace (sequence only)",
            FontId::proportional(11.0),
            t.text_muted,
        );
        return;
    };
    let n = m.len();
    if n == 0 || center >= n {
        return;
    }
    let xs = positions(trace, n);
    let lo = center.saturating_sub(half);
    let hi = center.saturating_add(half).saturating_add(1).min(n);
    let (Some(&first), Some(&last)) = (xs.get(lo), xs.get(hi.saturating_sub(1))) else {
        return;
    };
    let pad = 6.0;
    let (s0, s1) = (first - pad, last + pad);
    let span = (s1 - s0).max(1.0);
    let x_of = |sample: f32| -> f32 {
        let f = (sample - s0) / span;
        if reverse {
            rect.right() - f * rect.width()
        } else {
            rect.left() + f * rect.width()
        }
    };
    let letters_h = 15.0;
    let bottom = rect.bottom() - 2.0;
    let top = rect.top() + letters_h + 2.0;
    // The centre base.
    if let Some(&c) = xs.get(center) {
        let w = (rect.width() / (hi - lo).max(1) as f32).max(4.0);
        painter.rect_filled(
            Rect::from_center_size(pos2(x_of(c), rect.center().y), vec2(w, rect.height())),
            0.0,
            t.selection,
        );
    }
    let first_sample = s0.floor().max(0.0) as usize;
    let last_sample = (s1.ceil().max(0.0) as usize).min(trace.samples());
    let max = (first_sample..last_sample)
        .flat_map(|i| {
            Channel::ALL
                .into_iter()
                .filter_map(move |c| trace.value(c, i))
        })
        .max()
        .unwrap_or(1)
        .max(1);
    let complement = |c: Channel| match c {
        Channel::A => Channel::T,
        Channel::C => Channel::G,
        Channel::G => Channel::C,
        Channel::T => Channel::A,
    };
    let color_of = |c: Channel| {
        let shown = if reverse { complement(c) } else { c };
        let i = Channel::ALL.iter().position(|x| *x == shown).unwrap_or(0);
        t.trace_channels.get(i).copied().unwrap_or(t.text)
    };
    for c in Channel::ALL {
        let points: Vec<Pos2> = (first_sample..last_sample)
            .filter_map(|i| {
                let v = f32::from(trace.value(c, i)?) / f32::from(max);
                Some(pos2(x_of(i as f32), bottom - v * (bottom - top)))
            })
            .collect();
        if points.len() > 1 {
            painter.add(Shape::line(points, Stroke::new(1.0, color_of(c))));
        }
    }
    let residues = m.seq().as_bytes();
    for i in lo..hi {
        let (Some(&x), Some(&b)) = (xs.get(i), residues.get(i)) else {
            continue;
        };
        let base = if reverse {
            snappiesnap_engine::doc::Alphabet::Dna
                .complement(b)
                .unwrap_or(b'N')
        } else {
            b
        };
        let channel = match base.to_ascii_uppercase() {
            b'A' => Some(0),
            b'C' => Some(1),
            b'G' => Some(2),
            b'T' => Some(3),
            _ => None,
        };
        let color = channel
            .and_then(|k| t.trace_channels.get(k).copied())
            .unwrap_or(t.text_muted);
        let quality = trace.calls().get(i).and_then(|c| c.quality);
        let color = if quality.is_some_and(|q| q < GOOD_QUALITY) {
            color.gamma_multiply(0.55)
        } else {
            color
        };
        painter.text(
            pos2(x_of(x), rect.top() + 2.0),
            Align2::CENTER_TOP,
            char::from(base).to_string(),
            FontId::monospace(12.0),
            color,
        );
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use snappiesnap_engine::doc::{Alphabet, Call, Seq, Signal, Topology};

    fn read(peaks: &[Option<u32>]) -> Molecule {
        let n = peaks.len();
        let signal = Signal {
            a: vec![1; 200],
            c: vec![1; 200],
            g: vec![1; 200],
            t: vec![1; 200],
        };
        let calls = peaks
            .iter()
            .map(|&peak| Call {
                peak,
                quality: Some(30),
            })
            .collect();
        let (t, _) = Trace::new(signal, calls, None).unwrap();
        let mut m = Molecule::new(
            "r",
            Seq::new(Alphabet::Dna, "A".repeat(n)).unwrap(),
            Topology::Linear,
        )
        .unwrap();
        m.set_trace(Some(t)).unwrap();
        m
    }

    #[test]
    fn typed_in_residues_sit_between_their_neighbours() {
        let m = read(&[Some(10), None, None, Some(40), None]);
        let xs = positions(m.trace().unwrap(), m.len());
        assert_eq!(xs, vec![10.0, 20.0, 30.0, 40.0, 70.0]);
        let m = read(&[None, Some(50)]);
        assert_eq!(positions(m.trace().unwrap(), 2), vec![38.0, 50.0]);
        let m = read(&[None, None]);
        assert_eq!(positions(m.trace().unwrap(), 2), vec![6.0, 18.0]);
        assert!(positions(read(&[]).trace().unwrap(), 0).is_empty());
    }
}
