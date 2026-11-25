//! The colour picker used across the app (the map background in Preferences, feature colours):
//! a swatch that opens a picker with a saturation/value field and a hue bar, a `#rrggbb` box,
//! the theme's feature palette as quick choices, and "Default" to go back to the automatic
//! colour.

use eframe::egui::{
    self, Color32, CornerRadius, Mesh, Pos2, Rect, RichText, Sense, Shape, Stroke, Ui,
    ecolor::Hsva, pos2, vec2,
};

use crate::theme::{Tokens, hex, parse_hex};

/// What the user chose in a [`color_picker`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ColorChoice {
    Set(Color32),
    /// Go back to the default (the theme's or the automatic colour).
    Default,
}

/// A swatch showing `current` (or `default`, marked as automatic, when there is none) that
/// opens the picker when clicked. Returns the user's choice in the frame they make it.
pub fn color_picker(
    ui: &mut Ui,
    t: &Tokens,
    id_salt: impl std::hash::Hash + std::fmt::Debug,
    current: Option<Color32>,
    default: Color32,
) -> Option<ColorChoice> {
    let shown = current.unwrap_or(default);
    let (rect, response) = ui.allocate_exact_size(vec2(44.0, 20.0), Sense::click());
    let p = ui.painter();
    p.rect_filled(rect, CornerRadius::same(4), shown);
    let border = if response.hovered() {
        t.accent
    } else {
        t.border
    };
    p.rect_stroke(
        rect,
        CornerRadius::same(4),
        Stroke::new(1.0, border),
        egui::StrokeKind::Inside,
    );
    if current.is_none() {
        // A small "A" for automatic, in whichever of black or white reads on the swatch.
        p.text(
            rect.right_center() - vec2(5.0, 0.0),
            egui::Align2::RIGHT_CENTER,
            "A",
            egui::FontId::proportional(10.0),
            readable_on(shown),
        );
    }
    let response = response.on_hover_text(match current {
        Some(c) => format!("{} (click to change)", hex(c)),
        None => format!("Default, {} (click to change)", hex(default)),
    });
    let id = ui.make_persistent_id(("color-picker", id_salt));
    let mut choice = None;
    egui::Popup::menu(&response)
        .id(id)
        .close_behavior(egui::PopupCloseBehavior::CloseOnClickOutside)
        .show(|ui| {
            ui.set_min_width(220.0);
            // Hue, saturation and value are kept between frames, so the hue survives greys.
            let hsva_id = id.with("hsva");
            let mut hsva = ui
                .data_mut(|d| d.get_temp::<Hsva>(hsva_id))
                .filter(|h| Color32::from(*h) == shown)
                .unwrap_or_else(|| Hsva::from(shown));
            if hsv_field(ui, t, &mut hsva) {
                choice = Some(ColorChoice::Set(Color32::from(hsva)));
            }
            ui.data_mut(|d| d.insert_temp(hsva_id, hsva));
            ui.add_space(t.spacing);
            // The hex box keeps what is being typed until it is a whole colour.
            let text_id = id.with("hex");
            let mut text = ui
                .data_mut(|d| d.get_temp::<String>(text_id))
                .unwrap_or_else(|| hex(shown));
            if let Some(ColorChoice::Set(c)) = choice {
                // Picked in the field above: the box follows.
                text = hex(c);
            }
            ui.horizontal(|ui| {
                ui.label(RichText::new("Hex").color(t.text_muted));
                let r = ui.add(egui::TextEdit::singleline(&mut text).desired_width(80.0));
                if r.changed()
                    && let Some(c) = parse_hex(&text)
                {
                    choice = Some(ColorChoice::Set(c));
                }
            });
            ui.add_space(t.spacing);
            ui.horizontal_wrapped(|ui| {
                ui.spacing_mut().item_spacing = vec2(4.0, 4.0);
                for c in t.feature_palette {
                    let (r, swatch) = ui.allocate_exact_size(vec2(20.0, 20.0), Sense::click());
                    ui.painter().rect_filled(r, CornerRadius::same(3), c);
                    if swatch.hovered() {
                        ui.painter().rect_stroke(
                            r,
                            CornerRadius::same(3),
                            Stroke::new(1.5, t.text),
                            egui::StrokeKind::Outside,
                        );
                    }
                    if swatch.on_hover_text(hex(c)).clicked() {
                        choice = Some(ColorChoice::Set(c));
                        text = hex(c);
                    }
                }
            });
            ui.data_mut(|d| d.insert_temp(text_id, text));
            ui.add_space(t.spacing);
            if ui
                .add_enabled(current.is_some(), egui::Button::new("Default"))
                .on_hover_text(format!("Back to the default colour ({})", hex(default)))
                .clicked()
            {
                choice = Some(ColorChoice::Default);
                ui.data_mut(|d| d.remove::<String>(text_id));
                ui.close();
            }
        });
    choice
}

/// A saturation (across) / value (down) field for the current hue, and a hue bar below it.
/// Returns whether the user changed the colour.
fn hsv_field(ui: &mut Ui, t: &Tokens, hsva: &mut Hsva) -> bool {
    const STEPS: u16 = 16;
    let mut changed = false;
    let width = 200.0;
    // Saturation and value.
    let (rect, r) = ui.allocate_exact_size(vec2(width, 130.0), Sense::click_and_drag());
    let hue = hsva.h;
    paint_gradient(ui, rect, STEPS, STEPS, |x, y| {
        Color32::from(Hsva::new(hue, x, 1.0 - y, 1.0))
    });
    if let Some(p) = r.interact_pointer_pos()
        && (r.dragged() || r.clicked())
    {
        let (x, y) = fraction(rect, p);
        hsva.s = x;
        hsva.v = 1.0 - y;
        changed = true;
    }
    let marker = pos2(
        rect.left() + hsva.s * rect.width(),
        rect.top() + (1.0 - hsva.v) * rect.height(),
    );
    ui.painter().circle_stroke(
        marker,
        5.0,
        Stroke::new(1.5, readable_on(Color32::from(*hsva))),
    );
    ui.add_space(4.0);
    // Hue.
    let (bar, r) = ui.allocate_exact_size(vec2(width, 14.0), Sense::click_and_drag());
    paint_gradient(ui, bar, 36, 1, |x, _| {
        Color32::from(Hsva::new(x, 1.0, 1.0, 1.0))
    });
    if let Some(p) = r.interact_pointer_pos()
        && (r.dragged() || r.clicked())
    {
        hsva.h = fraction(bar, p).0;
        changed = true;
    }
    let x = bar.left() + hsva.h * bar.width();
    ui.painter().rect_stroke(
        Rect::from_center_size(pos2(x, bar.center().y), vec2(5.0, bar.height() + 4.0)),
        CornerRadius::same(2),
        Stroke::new(1.5, t.text),
        egui::StrokeKind::Middle,
    );
    changed
}

/// Where `p` lies in `rect`, as fractions 0..=1 across and down.
fn fraction(rect: Rect, p: Pos2) -> (f32, f32) {
    let f = |v: f32, lo: f32, size: f32| {
        if size > 0.0 {
            ((v - lo) / size).clamp(0.0, 1.0)
        } else {
            0.0
        }
    };
    (
        f(p.x, rect.left(), rect.width()),
        f(p.y, rect.top(), rect.height()),
    )
}

/// Fill `rect` with a grid of `nx` by `ny` cells whose corners take `color_at(x, y)` (fractions
/// across and down), blended between corners.
fn paint_gradient(ui: &Ui, rect: Rect, nx: u16, ny: u16, color_at: impl Fn(f32, f32) -> Color32) {
    let mut mesh = Mesh::default();
    for j in 0..=ny {
        for i in 0..=nx {
            let (x, y) = (f32::from(i) / f32::from(nx), f32::from(j) / f32::from(ny));
            mesh.colored_vertex(
                pos2(
                    rect.left() + x * rect.width(),
                    rect.top() + y * rect.height(),
                ),
                color_at(x, y),
            );
        }
    }
    let row = u32::from(nx) + 1;
    for j in 0..u32::from(ny) {
        for i in 0..u32::from(nx) {
            let a = j * row + i;
            mesh.add_triangle(a, a + 1, a + row);
            mesh.add_triangle(a + 1, a + row + 1, a + row);
        }
    }
    ui.painter().add(Shape::mesh(mesh));
}

/// Black or white, whichever reads better on `background` (by relative luminance, ITU-R
/// BT.709 weights on sRGB values; good enough for a small label).
pub fn readable_on(background: Color32) -> Color32 {
    let l = 0.2126 * f32::from(background.r())
        + 0.7152 * f32::from(background.g())
        + 0.0722 * f32::from(background.b());
    if l > 140.0 {
        Color32::BLACK
    } else {
        Color32::WHITE
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn readable_labels() {
        assert_eq!(readable_on(Color32::WHITE), Color32::BLACK);
        assert_eq!(
            readable_on(Color32::from_rgb(0x15, 0x17, 0x1b)),
            Color32::WHITE
        );
    }

    /// The swatch draws with and without a colour of its own, and its popup opens; set
    /// `SNAPPIESNAP_COLOR_PICKER_PNG=/path.png` to look at it.
    #[test]
    fn draws() {
        let t = Tokens::dark();
        let ctx = egui::Context::default();
        t.apply(&ctx);
        let mut opened = false;
        for frame in 0..3 {
            let mut out = ctx.run_ui(egui::RawInput::default(), |ui| {
                ui.push_id("color-picker-test", |ui| {
                    color_picker(ui, &t, "a", None, t.feature_color(0));
                    let r = color_picker(ui, &t, "b", Some(Color32::RED), t.canvas);
                    assert_eq!(r, None);
                    if frame == 1 {
                        let popup = ui.make_persistent_id(("color-picker", "b"));
                        egui::Popup::open_id(ui.ctx(), popup);
                        opened = true;
                    }
                });
            });
            out.textures_delta.clear();
        }
        assert!(opened);
        let image = crate::offscreen::render(
            300,
            380,
            1.0,
            |ctx| t.apply(ctx),
            |ui| {
                let popup = ui.make_persistent_id(("color-picker", "c"));
                egui::Popup::open_id(ui.ctx(), popup);
                color_picker(
                    ui,
                    &t,
                    "c",
                    Some(Color32::from_rgb(0x2b, 0x7b, 0xc4)),
                    t.canvas,
                );
            },
        )
        .unwrap();
        if let Ok(path) = std::env::var("SNAPPIESNAP_COLOR_PICKER_PNG") {
            crate::offscreen::write_png(&path, &image).unwrap();
        }
    }
}
