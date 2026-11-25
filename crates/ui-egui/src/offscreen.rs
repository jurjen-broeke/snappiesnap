//! Offscreen rendering: run UI code in a windowless egui context and rasterise egui's output in
//! software. No GPU or display is needed, so `ui.render` works with a locked screen, in CI and
//! in tests, and the result is deterministic.
//!
//! The rasteriser fills egui's textured, premultiplied-alpha triangles and blends them in gamma
//! space, as egui's own renderers do.

use std::collections::HashMap;

use eframe::egui::{
    self, Color32, ColorImage, Pos2, RawInput, Rect, TextureId, Ui, ViewportId, epaint, pos2, vec2,
};

/// Largest image side we render, in pixels.
pub const MAX_SIDE: u32 = 8192;

/// Render `draw` into a `width`×`height` pixel image at `scale` pixels per point.
pub fn render(
    width: u32,
    height: u32,
    scale: f32,
    setup: impl FnOnce(&egui::Context),
    mut draw: impl FnMut(&mut Ui),
) -> Result<ColorImage, String> {
    if !(1..=MAX_SIDE).contains(&width) || !(1..=MAX_SIDE).contains(&height) {
        return Err(format!("image size must be 1..={MAX_SIDE} pixels per side"));
    }
    if !(0.25..=8.0).contains(&scale) {
        return Err("scale must be between 0.25 and 8".into());
    }
    let ctx = egui::Context::default();
    setup(&ctx);
    let mut textures: HashMap<TextureId, ColorImage> = HashMap::new();
    let points = vec2(width as f32 / scale, height as f32 / scale);
    let mut output = None;
    // A few passes so layout that depends on the previous frame (text sizes) settles.
    for pass in 0..3u32 {
        let mut input = RawInput {
            screen_rect: Some(Rect::from_min_size(Pos2::ZERO, points)),
            time: Some(f64::from(pass) * 0.1),
            ..Default::default()
        };
        input
            .viewports
            .entry(ViewportId::ROOT)
            .or_default()
            .native_pixels_per_point = Some(scale);
        let mut out = ctx.run_ui(input, &mut draw);
        apply_textures(&mut textures, &out.textures_delta);
        // Handled: tell egui so (it checks that texture updates aren't silently dropped).
        out.textures_delta.clear();
        output = Some(out);
    }
    let out = output.ok_or("nothing was drawn")?;
    let primitives = ctx.tessellate(out.shapes, out.pixels_per_point);
    Ok(rasterize(
        &primitives,
        &textures,
        width,
        height,
        out.pixels_per_point,
    ))
}

/// Run `draw` in a windowless egui context of `width`×`height` points and return what it
/// painted, as egui shapes (for vector output, see [`crate::svg`]).
pub fn shapes(
    width: u32,
    height: u32,
    setup: impl FnOnce(&egui::Context),
    mut draw: impl FnMut(&mut Ui),
) -> Result<Vec<epaint::ClippedShape>, String> {
    if !(1..=MAX_SIDE).contains(&width) || !(1..=MAX_SIDE).contains(&height) {
        return Err(format!("image size must be 1..={MAX_SIDE} points per side"));
    }
    let ctx = egui::Context::default();
    setup(&ctx);
    let points = vec2(width as f32, height as f32);
    let mut shapes = Vec::new();
    // A few passes so layout that depends on the previous frame (text sizes) settles.
    for pass in 0..3u32 {
        let mut input = RawInput {
            screen_rect: Some(Rect::from_min_size(Pos2::ZERO, points)),
            time: Some(f64::from(pass) * 0.1),
            ..Default::default()
        };
        input
            .viewports
            .entry(ViewportId::ROOT)
            .or_default()
            .native_pixels_per_point = Some(1.0);
        let mut out = ctx.run_ui(input, &mut draw);
        // Handled: textures aren't needed for vector output.
        out.textures_delta.clear();
        shapes = out.shapes;
    }
    Ok(shapes)
}

/// Write an image as an RGBA PNG file.
#[cfg(not(target_arch = "wasm32"))]
pub fn write_png(path: &str, image: &ColorImage) -> Result<(), String> {
    std::fs::write(path, png_bytes(image)?).map_err(|e| e.to_string())
}

/// An image as RGBA PNG bytes.
pub fn png_bytes(image: &ColorImage) -> Result<Vec<u8>, String> {
    let [w, h] = image.size;
    let width = u32::try_from(w).map_err(|e| e.to_string())?;
    let height = u32::try_from(h).map_err(|e| e.to_string())?;
    let mut out = Vec::new();
    let mut encoder = png::Encoder::new(&mut out, width, height);
    encoder.set_color(png::ColorType::Rgba);
    encoder.set_depth(png::BitDepth::Eight);
    let mut writer = encoder.write_header().map_err(|e| e.to_string())?;
    let mut rgba = Vec::with_capacity(image.pixels.len().saturating_mul(4));
    for p in &image.pixels {
        rgba.extend_from_slice(&p.to_srgba_unmultiplied());
    }
    writer.write_image_data(&rgba).map_err(|e| e.to_string())?;
    writer.finish().map_err(|e| e.to_string())?;
    Ok(out)
}

fn apply_textures(
    textures: &mut HashMap<TextureId, ColorImage>,
    delta: &epaint::textures::TexturesDelta,
) {
    for (id, deltas) in &delta.set {
        for d in deltas {
            let epaint::ImageData::Color(image) = &d.image;
            match d.pos {
                None => {
                    textures.insert(*id, (**image).clone());
                }
                Some([x0, y0]) => {
                    let Some(target) = textures.get_mut(id) else {
                        continue;
                    };
                    let [tw, _] = target.size;
                    let [pw, ph] = image.size;
                    for y in 0..ph {
                        for x in 0..pw {
                            let src = image.pixels.get(y * pw + x).copied();
                            let dst = (y0 + y)
                                .checked_mul(tw)
                                .and_then(|row| row.checked_add(x0 + x));
                            if let (Some(src), Some(dst)) = (src, dst)
                                && let Some(slot) = target.pixels.get_mut(dst)
                            {
                                *slot = src;
                            }
                        }
                    }
                }
            }
        }
    }
    for id in &delta.free {
        textures.remove(id);
    }
}

type Rgba = [f32; 4];

fn to_f(c: Color32) -> Rgba {
    let [r, g, b, a] = c.to_array();
    [
        f32::from(r) / 255.0,
        f32::from(g) / 255.0,
        f32::from(b) / 255.0,
        f32::from(a) / 255.0,
    ]
}

/// Bilinear sample of a premultiplied texture at normalised `uv`.
fn sample(tex: Option<&ColorImage>, u: f32, v: f32) -> Rgba {
    let Some(tex) = tex else {
        return [1.0; 4];
    };
    let [w, h] = tex.size;
    if w == 0 || h == 0 {
        return [1.0; 4];
    }
    let fx = (u * w as f32 - 0.5).clamp(0.0, (w - 1) as f32);
    let fy = (v * h as f32 - 0.5).clamp(0.0, (h - 1) as f32);
    let (x0, y0) = (fx.floor() as usize, fy.floor() as usize);
    let (x1, y1) = ((x0 + 1).min(w - 1), (y0 + 1).min(h - 1));
    let (tx, ty) = (fx - x0 as f32, fy - y0 as f32);
    let px = |x: usize, y: usize| tex.pixels.get(y * w + x).copied().map_or([0.0; 4], to_f);
    let lerp = |t: f32| move |a: f32, b: f32| a + (b - a) * t;
    let top = zip4(px(x0, y0), px(x1, y0), lerp(tx));
    let bottom = zip4(px(x0, y1), px(x1, y1), lerp(tx));
    zip4(top, bottom, lerp(ty))
}

/// Combine two colours channel by channel.
fn zip4(a: Rgba, b: Rgba, f: impl Fn(f32, f32) -> f32) -> Rgba {
    let [a0, a1, a2, a3] = a;
    let [b0, b1, b2, b3] = b;
    [f(a0, b0), f(a1, b1), f(a2, b2), f(a3, b3)]
}

fn rasterize(
    primitives: &[epaint::ClippedPrimitive],
    textures: &HashMap<TextureId, ColorImage>,
    width: u32,
    height: u32,
    ppp: f32,
) -> ColorImage {
    let (w, h) = (width as usize, height as usize);
    let mut buf: Vec<Rgba> = vec![[0.0, 0.0, 0.0, 1.0]; w * h];
    for prim in primitives {
        let epaint::Primitive::Mesh(mesh) = &prim.primitive else {
            continue;
        };
        let clip = Rect::from_min_max(
            pos2(prim.clip_rect.min.x * ppp, prim.clip_rect.min.y * ppp),
            pos2(prim.clip_rect.max.x * ppp, prim.clip_rect.max.y * ppp),
        )
        .intersect(Rect::from_min_size(Pos2::ZERO, vec2(w as f32, h as f32)));
        if !clip.is_positive() {
            continue;
        }
        let tex = textures.get(&mesh.texture_id);
        let (triangles, _) = mesh.indices.as_chunks::<3>();
        for &[i0, i1, i2] in triangles {
            let vertex = |i: u32| mesh.vertices.get(i as usize);
            let (Some(a), Some(b), Some(c)) = (vertex(i0), vertex(i1), vertex(i2)) else {
                continue;
            };
            fill_triangle(&mut buf, w, clip, ppp, [a, b, c], tex);
        }
    }
    let pixels = buf
        .iter()
        .map(|&[r, g, b, _]| {
            let q = |x: f32| (x.clamp(0.0, 1.0) * 255.0).round() as u8;
            Color32::from_rgba_premultiplied(q(r), q(g), q(b), 255)
        })
        .collect();
    ColorImage::new([w, h], pixels)
}

fn fill_triangle(
    buf: &mut [Rgba],
    w: usize,
    clip: Rect,
    ppp: f32,
    v: [&epaint::Vertex; 3],
    tex: Option<&ColorImage>,
) {
    let p = v.map(|v| pos2(v.pos.x * ppp, v.pos.y * ppp));
    let [p0, p1, p2] = p;
    let [v0, v1, v2] = v;
    let area = (p1.x - p0.x) * (p2.y - p0.y) - (p2.x - p0.x) * (p1.y - p0.y);
    if area.abs() < 1e-6 {
        return;
    }
    let min_x = p
        .iter()
        .map(|q| q.x)
        .fold(f32::MAX, f32::min)
        .max(clip.min.x)
        .floor();
    let max_x = p
        .iter()
        .map(|q| q.x)
        .fold(f32::MIN, f32::max)
        .min(clip.max.x)
        .ceil();
    let min_y = p
        .iter()
        .map(|q| q.y)
        .fold(f32::MAX, f32::min)
        .max(clip.min.y)
        .floor();
    let max_y = p
        .iter()
        .map(|q| q.y)
        .fold(f32::MIN, f32::max)
        .min(clip.max.y)
        .ceil();
    if min_x >= max_x || min_y >= max_y {
        return;
    }
    let [c0, c1, c2] = [v0, v1, v2].map(|v| to_f(v.color));
    let edge = |a: Pos2, b: Pos2, x: f32, y: f32| (b.x - a.x) * (y - a.y) - (b.y - a.y) * (x - a.x);
    for y in (min_y as usize)..(max_y as usize) {
        let py = y as f32 + 0.5;
        for x in (min_x as usize)..(max_x as usize) {
            let px = x as f32 + 0.5;
            let w0 = edge(p1, p2, px, py) / area;
            let w1 = edge(p2, p0, px, py) / area;
            let w2 = 1.0 - w0 - w1;
            if w0 < -1e-4 || w1 < -1e-4 || w2 < -1e-4 {
                continue;
            }
            let u = w0 * v0.uv.x + w1 * v1.uv.x + w2 * v2.uv.x;
            let vv = w0 * v0.uv.y + w1 * v1.uv.y + w2 * v2.uv.y;
            let color = zip4(zip4(c0, c1, |a, b| w0 * a + w1 * b), c2, |ab, c| {
                ab + w2 * c
            });
            let src = zip4(color, sample(tex, u, vv), |c, t| c * t);
            if let Some(dst) = buf.get_mut(y * w + x) {
                let [_, _, _, alpha] = src;
                let keep = 1.0 - alpha;
                *dst = zip4(src, *dst, |s, d| s + d * keep);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn renders_shapes_and_text_deterministically() {
        let draw = |ui: &mut Ui| {
            let painter = ui.painter();
            painter.rect_filled(
                Rect::from_min_size(pos2(0.0, 0.0), vec2(100.0, 60.0)),
                0.0,
                Color32::from_rgb(200, 40, 40),
            );
            ui.label("Hello");
        };
        let a = render(120, 80, 1.0, |_| {}, draw).unwrap();
        let b = render(120, 80, 1.0, |_| {}, draw).unwrap();
        assert_eq!(a.size, [120, 80]);
        assert_eq!(a.pixels, b.pixels);
        // Inside the red rectangle, away from the label.
        let red = a.pixels[50 * 120 + 50];
        assert!(red.r() > 150 && red.g() < 80, "{red:?}");
        // Text pixels exist: something differs from the red fill near the label.
        let label_area = &a.pixels[2 * 120..20 * 120];
        assert!(label_area.iter().any(|p| p.g() > 100));
    }

    #[test]
    fn rejects_bad_sizes() {
        assert!(render(0, 10, 1.0, |_| {}, |_| {}).is_err());
        assert!(render(10, MAX_SIDE + 1, 1.0, |_| {}, |_| {}).is_err());
        assert!(render(10, 10, 0.0, |_| {}, |_| {}).is_err());
    }
}
