//! The Gel view: simulated agarose gels, one for restriction digests (uncut molecule and
//! digests) and one for PCR products, each with its own ladder and lanes. Bands come from
//! `snappiesnap_engine::gel`, the same code the `gel` command uses; the gels are UI state
//! ([`GelSetup`]) so agents can set them with `ui.set`. A gel across several constructs is an
//! experiment ([`crate::experiment_view`], from the toolbar's Experiments menu).

use eframe::egui::{self, Align2, FontId, Pos2, Rect, RichText, Sense, Stroke, Ui, pos2, vec2};
use snappiesnap_engine::Session;
use snappiesnap_engine::bio::gel::{Band, LADDERS};
use snappiesnap_engine::gel::{LaneSpec, bands};

use crate::state::{GelLanes, GelMode, GelSetup};
use crate::theme::Tokens;

/// Inputs of the "add lane" controls.
#[derive(Default)]
pub struct GelInputs {
    pub enzymes: String,
    pub forward: String,
    pub reverse: String,
}

/// A lane as drawn: title and bands, or why it can't be run.
pub(crate) type Lanes = Vec<snappiesnap_engine::gel::RunLane>;

/// Bands per (document revision, setup), so nothing is recomputed per frame.
#[derive(Default)]
pub struct GelCache {
    key: Option<(u64, String)>,
    lanes: Lanes,
}

impl GelCache {
    pub(crate) fn lanes(&mut self, session: &Session, setup: &GelLanes) -> &Lanes {
        let key = (
            session.revision(),
            serde_json::to_string(setup).unwrap_or_default(),
        );
        if self.key.as_ref() != Some(&key) {
            self.lanes = match session.document() {
                Some(doc) => {
                    let ladder = LaneSpec::Ladder {
                        name: setup.ladder.clone(),
                    };
                    std::iter::once(&ladder)
                        .chain(setup.lanes.iter())
                        .map(|spec| {
                            let b = bands(&doc.molecule, session.enzymes(), spec, setup.agarose);
                            (spec.title(), b)
                        })
                        .collect()
                }
                None => Vec::new(),
            };
            self.key = Some(key);
        }
        &self.lanes
    }
}

/// Draw the Gel tab. Returns true when the user asked to save the gel as an image.
pub fn show(
    ui: &mut Ui,
    t: &Tokens,
    session: &Session,
    setup: &mut GelSetup,
    inputs: &mut GelInputs,
    cache: &mut GelCache,
) -> bool {
    let Some(doc) = session.document() else {
        return false;
    };
    let mut save = false;
    let primers: Vec<String> = doc
        .molecule
        .primers()
        .iter()
        .map(|p| p.name.clone())
        .collect();
    // Two separate gels, one per method.
    ui.horizontal(|ui| {
        for (mode, label) in [
            (GelMode::Digest, "Restriction digest"),
            (GelMode::Pcr, "PCR"),
        ] {
            let on = setup.mode == mode;
            let text = RichText::new(label).color(if on { t.accent_text } else { t.text_muted });
            if ui.add(egui::Button::selectable(on, text)).clicked() {
                setup.mode = mode;
            }
        }
    });
    ui.add_space(t.spacing);
    let mode = setup.mode;
    let gel = setup.current();
    ui.horizontal_top(|ui| {
        // Controls.
        ui.vertical(|ui| {
            ui.set_width(270.0);
            ui.horizontal(|ui| {
                ui.label("Agarose");
                ui.add(
                    egui::DragValue::new(&mut gel.agarose)
                        .range(0.3..=2.0)
                        .speed(0.05)
                        .fixed_decimals(1)
                        .suffix(" %"),
                );
                ui.label("Ladder");
                egui::ComboBox::from_id_salt(("ladder", mode == GelMode::Pcr))
                    .selected_text(&gel.ladder)
                    .show_ui(ui, |ui| {
                        for (name, _) in LADDERS {
                            ui.selectable_value(&mut gel.ladder, name.to_string(), *name);
                        }
                    });
            });
            ui.add_space(t.spacing);
            ui.label(RichText::new("Lanes").strong());
            let mut remove = None;
            for (i, lane) in gel.lanes.iter().enumerate() {
                ui.horizontal(|ui| {
                    ui.label(RichText::new(format!("{}", i + 2)).color(t.text_muted));
                    ui.add(egui::Label::new(lane.title()).truncate());
                    ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                        if ui.small_button("Remove").clicked() {
                            remove = Some(i);
                        }
                    });
                });
            }
            if gel.lanes.is_empty() {
                ui.label(RichText::new("Only the ladder so far.").color(t.text_muted));
            }
            if let Some(i) = remove {
                gel.lanes.remove(i);
            }
            ui.add_space(t.spacing);
            ui.label(RichText::new("Add a lane").strong());
            match mode {
                GelMode::Digest => digest_controls(ui, gel, inputs),
                GelMode::Pcr => pcr_controls(ui, t, gel, inputs, &primers),
            }
            ui.add_space(t.spacing * 2.0);
            save = ui
                .button("Save image\u{2026}")
                .on_hover_text("Save this gel as a PNG or JPEG image, as shown or inverted")
                .clicked();
            ui.add_space(t.spacing);
            let note = match mode {
                GelMode::Digest => "Migration: log-linear over the agarose's separation range (Sambrook & Russell 2001). Brightness follows DNA mass. Uncut plasmid is supercoiled and is placed at about 70 % of its length, which varies in practice.",
                GelMode::Pcr => "Migration: log-linear over the agarose's separation range (Sambrook & Russell 2001). Brightness follows DNA mass. Each lane holds the products of one primer pair.",
            };
            ui.label(RichText::new(note).small().color(t.text_muted));
        });
        ui.separator();
        let lanes = cache.lanes(session, gel);
        draw_gel(ui, t, lanes, None);
    });
    save
}

fn digest_controls(ui: &mut Ui, gel: &mut GelLanes, inputs: &mut GelInputs) {
    if ui.button("Uncut").clicked() {
        gel.lanes.push(LaneSpec::Uncut);
    }
    ui.horizontal(|ui| {
        ui.add(
            egui::TextEdit::singleline(&mut inputs.enzymes)
                .hint_text("EcoRI, BamHI")
                .desired_width(170.0),
        );
        let enzymes: Vec<String> = inputs
            .enzymes
            .split([',', ' ', '+'])
            .map(str::trim)
            .filter(|s| !s.is_empty())
            .map(str::to_string)
            .collect();
        if ui
            .add_enabled(!enzymes.is_empty(), egui::Button::new("Digest"))
            .clicked()
        {
            gel.lanes.push(LaneSpec::Digest { enzymes });
        }
    });
}

fn pcr_controls(
    ui: &mut Ui,
    t: &Tokens,
    gel: &mut GelLanes,
    inputs: &mut GelInputs,
    primers: &[String],
) {
    if primers.is_empty() {
        ui.label(
            RichText::new("Add primers to simulate a PCR.")
                .small()
                .color(t.text_muted),
        );
        return;
    }
    ui.horizontal(|ui| {
        for (k, (id, value)) in [("fwd", &mut inputs.forward), ("rev", &mut inputs.reverse)]
            .into_iter()
            .enumerate()
        {
            if !primers.contains(value) {
                // Forward defaults to the first primer, reverse to the second.
                *value = primers
                    .get(k)
                    .or(primers.first())
                    .cloned()
                    .unwrap_or_default();
            }
            egui::ComboBox::from_id_salt(id)
                .width(80.0)
                .selected_text(value.as_str())
                .show_ui(ui, |ui| {
                    for p in primers {
                        ui.selectable_value(value, p.clone(), p);
                    }
                });
        }
        if ui.button("Add PCR").clicked() {
            gel.lanes.push(LaneSpec::Pcr {
                forward: inputs.forward.clone(),
                reverse: inputs.reverse.clone(),
            });
        }
    });
}

/// Draw a gel. `selected` (a lane index, 0 = the ladder) is outlined; returns the lane the
/// user clicked, if any.
pub(crate) fn draw_gel(
    ui: &mut Ui,
    t: &Tokens,
    lanes: &[snappiesnap_engine::gel::RunLane],
    selected: Option<usize>,
) -> Option<usize> {
    let avail = ui.available_size();
    let lane_w = 54.0;
    let width = (lanes.len() as f32 * lane_w + 70.0).min(avail.x.max(200.0));
    let height = avail.y.clamp(200.0, 900.0);
    let sense = if selected.is_some() {
        Sense::click()
    } else {
        Sense::hover()
    };
    let (rect, response) = ui.allocate_exact_size(vec2(width, height), sense);
    let p = ui.painter_at(rect);
    let label_w = 52.0;
    let gel = Rect::from_min_max(rect.min + vec2(label_w, 30.0), rect.max - vec2(4.0, 8.0));
    p.rect_filled(gel, 6.0, t.gel);
    let lane_w = (gel.width() / lanes.len().max(1) as f32).min(lane_w * 1.4);
    let top = gel.top() + 14.0;
    let run = gel.height() - 24.0;
    let y_at = |y: f64| top + run * y as f32;
    let hover = response.hover_pos();
    let mut hovered: Option<(&Band, &str)> = None;
    let clicked = response
        .clicked()
        .then(|| response.interact_pointer_pos())
        .flatten()
        .filter(|p| gel.expand2(vec2(0.0, 30.0)).contains(*p))
        .map(|p| ((p.x - gel.left()) / lane_w).floor())
        .filter(|i| *i >= 0.0 && (*i as usize) < lanes.len())
        .map(|i| i as usize);
    for (i, (title, result)) in lanes.iter().enumerate() {
        let cx = gel.left() + lane_w * (i as f32 + 0.5);
        let band_w = lane_w * 0.72;
        if selected == Some(i) {
            p.rect_stroke(
                Rect::from_center_size(
                    pos2(cx, gel.center().y),
                    vec2(lane_w - 4.0, gel.height() - 4.0),
                ),
                4.0,
                Stroke::new(1.5, t.accent),
                egui::StrokeKind::Inside,
            );
        }
        // Title above the well.
        let galley = p.layout_no_wrap(title.clone(), FontId::proportional(10.0), t.text_muted);
        let title_rect =
            Align2::CENTER_BOTTOM.anchor_size(pos2(cx, gel.top() - 4.0), galley.size());
        if galley.size().x <= lane_w + 6.0 {
            p.galley(title_rect.min, galley, t.text_muted);
        } else {
            p.text(
                pos2(cx, gel.top() - 4.0),
                Align2::CENTER_BOTTOM,
                format!("{}", i + 1),
                FontId::proportional(10.0),
                t.text_muted,
            );
        }
        if let Some(h) = hover
            && (h.x - cx).abs() < lane_w / 2.0
            && h.y < gel.top()
        {
            response.clone().on_hover_text_at_pointer(title.as_str());
        }
        // Well.
        p.rect_filled(
            Rect::from_center_size(pos2(cx, gel.top() + 7.0), vec2(band_w, 5.0)),
            1.0,
            t.gel_well,
        );
        let bands = match result {
            Ok(b) => b,
            Err(e) => {
                p.text(
                    pos2(cx, gel.top() + 30.0),
                    Align2::CENTER_TOP,
                    "!",
                    FontId::proportional(14.0),
                    t.error,
                );
                if let Some(h) = hover
                    && (h.x - cx).abs() < lane_w / 2.0
                    && h.y > gel.top()
                {
                    response.clone().on_hover_text_at_pointer(e.as_str());
                }
                continue;
            }
        };
        let mut last_label: Option<f32> = None;
        for band in bands {
            let y = y_at(band.y);
            // A faint band stays visible: brightness grows with the square root of the mass.
            let alpha = (0.18 + 0.82 * band.intensity.clamp(0.0, 1.0).sqrt()) as f32;
            let r = Rect::from_center_size(pos2(cx, y), vec2(band_w, 4.0));
            p.rect_filled(r, 1.5, t.gel_band.gamma_multiply(alpha));
            // Ladder sizes on the left, skipping any that would overlap the one above.
            if i == 0 && last_label.is_none_or(|prev: f32| y - prev >= 11.0) {
                last_label = Some(y);
                p.text(
                    pos2(gel.left() - 6.0, y),
                    Align2::RIGHT_CENTER,
                    size_text(band.bp),
                    FontId::proportional(10.0),
                    t.text_muted,
                );
                p.line_segment(
                    [pos2(gel.left() - 4.0, y), pos2(gel.left(), y)],
                    Stroke::new(1.0, t.text_muted),
                );
            }
            if hover.is_some_and(|h: Pos2| r.expand2(vec2(0.0, 3.0)).contains(h)) {
                hovered = Some((band, title.as_str()));
                p.rect_stroke(
                    r.expand(2.0),
                    2.0,
                    Stroke::new(1.0, t.accent),
                    egui::StrokeKind::Outside,
                );
            }
        }
    }
    if let Some((band, title)) = hovered {
        let mut text = format!("{title}\n");
        for f in &band.fragments {
            text.push_str(&format!("{} bp  {}\n", group(f.bp), f.label));
        }
        response.on_hover_text_at_pointer(text.trim_end().to_string());
    }
    clicked
}

/// "10 kb", "1.5 kb", "750 bp".
fn size_text(bp: usize) -> String {
    if bp >= 1000 && bp.is_multiple_of(100) {
        let kb = bp as f64 / 1000.0;
        if bp.is_multiple_of(1000) {
            format!("{kb:.0} kb")
        } else {
            format!("{kb:.1} kb")
        }
    } else {
        format!("{bp}")
    }
}

/// Thousands separated: 12,345.
fn group(n: usize) -> String {
    let s = n.to_string();
    let mut out = String::new();
    for (i, c) in s.chars().enumerate() {
        if i > 0 && (s.len() - i).is_multiple_of(3) {
            out.push(',');
        }
        out.push(c);
    }
    out
}

#[cfg(test)]
mod tests {
    use snappiesnap_engine::Document;
    use snappiesnap_engine::doc::{Alphabet, Molecule, Seq, Topology};

    use super::*;

    #[test]
    fn labels() {
        assert_eq!(size_text(10000), "10 kb");
        assert_eq!(size_text(1500), "1.5 kb");
        assert_eq!(size_text(750), "750");
        assert_eq!(group(1234567), "1,234,567");
        assert_eq!(group(12), "12");
    }

    /// The view draws a gel with a ladder, uncut, digest, PCR and a broken lane. Set
    /// `SNAPPIESNAP_GEL_PNG=/path.png` to look at it.
    #[test]
    fn draws() {
        let template = format!(
            "{}GAATTC{}GGATCC{}CACCGTTGCATGCAAGCTTG{}TCGTGACTGGGAAAACCCTG{}",
            "A".repeat(300),
            "C".repeat(1200),
            "T".repeat(500),
            "G".repeat(900),
            "A".repeat(2000)
        );
        let mut m = Molecule::new(
            "m",
            Seq::new(Alphabet::Dna, template).unwrap(),
            Topology::Circular,
        )
        .unwrap();
        m.add_primer(
            "f",
            Seq::new(Alphabet::Dna, "CACCGTTGCATGCAAGCTTG").unwrap(),
            "",
        )
        .unwrap();
        m.add_primer(
            "r",
            Seq::new(Alphabet::Dna, "CAGGGTTTTCCCAGTCACGA").unwrap(),
            "",
        )
        .unwrap();
        let mut session = Session::new();
        session.open_document(Document::new(m));
        let digest = GelLanes {
            lanes: vec![
                LaneSpec::Uncut,
                LaneSpec::Digest {
                    enzymes: vec!["EcoRI".into()],
                },
                LaneSpec::Digest {
                    enzymes: vec!["EcoRI".into(), "BamHI".into()],
                },
                LaneSpec::Digest {
                    enzymes: vec!["NoSuchI".into()],
                },
            ],
            ..GelLanes::default()
        };
        let pcr = GelLanes {
            lanes: vec![LaneSpec::Pcr {
                forward: "f".into(),
                reverse: "r".into(),
            }],
            ..GelSetup::default().pcr
        };
        let mut setup = GelSetup {
            digest: digest.clone(),
            pcr: pcr.clone(),
            ..GelSetup::default()
        };
        let mut cache = GelCache::default();
        let lanes = cache.lanes(&session, &digest).clone();
        assert_eq!(lanes.len(), 5);
        assert!(lanes[4].1.is_err());
        let lanes = cache.lanes(&session, &pcr).clone();
        let product = lanes[1].1.as_ref().unwrap();
        assert_eq!(product.len(), 1);
        assert_eq!(product[0].bp, 940);
        assert_eq!(lanes[0].0, "100 bp");
        let t = Tokens::dark();
        let mut inputs = GelInputs::default();
        for mode in [GelMode::Digest, GelMode::Pcr] {
            setup.mode = mode;
            let image = crate::offscreen::render(
                900,
                520,
                1.0,
                |ctx| t.apply(ctx),
                |ui| {
                    egui::Frame::new().fill(t.panel).show(ui, |ui| {
                        show(ui, &t, &session, &mut setup, &mut inputs, &mut cache);
                    });
                },
            )
            .unwrap();
            if let Ok(path) = std::env::var("SNAPPIESNAP_GEL_PNG") {
                let path = path.replace(".png", &format!("_{mode:?}.png"));
                crate::offscreen::write_png(&path, &image).unwrap();
            }
        }
    }
}
