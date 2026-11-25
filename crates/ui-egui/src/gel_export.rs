//! Gels as image files (PNG or JPEG): the gel drawn exactly as in the Gel tab or an experiment
//! tab ([`crate::gel_view::draw_gel`], rendered offscreen), with a legend naming each lane
//! underneath, since lane titles are often too long to fit over their lanes. As shown, bands
//! are light on a dark gel, like a photo on a transilluminator; inverted, they are dark on
//! white, for print and lab notebooks.

use eframe::egui::{self, Color32, ColorImage, RichText};
use snappiesnap_engine::gel::RunLane;

use crate::theme::Tokens;

/// Width of one lane in the image, in points (the drawing takes at most this much).
const LANE: f32 = 64.0;
/// Pixels per point: images are drawn at twice screen resolution, so they print crisply.
const SCALE: f32 = 2.0;
/// Height of the gel itself, in points.
const GEL_HEIGHT: f32 = 620.0;

/// The gel as an image: `title` above, the gel, and a legend of the lanes. `invert` gives dark
/// bands on white.
pub fn gel_image(title: &str, lanes: &[RunLane], invert: bool) -> Result<ColorImage, String> {
    if lanes.is_empty() {
        return Err("the gel has no lanes".to_string());
    }
    // Always drawn in the dark theme: a gel photo is dark, and inverting it gives a clean white.
    let t = Tokens {
        canvas: Tokens::dark().gel,
        panel: Tokens::dark().gel,
        ..Tokens::dark()
    };
    let width = (lanes.len() as f32 * LANE + 90.0).max(420.0);
    let legend = lanes.len() as f32 * 17.0 + 30.0;
    let height = GEL_HEIGHT + legend + 50.0;
    let px = |points: f32| (points * SCALE).ceil() as u32;
    let mut image = crate::offscreen::render(
        px(width),
        px(height),
        SCALE,
        |ctx| t.apply(ctx),
        |ui| {
            egui::Frame::new()
                .fill(t.gel)
                .inner_margin(12)
                .show(ui, |ui| {
                    ui.set_min_size(ui.available_size());
                    ui.label(RichText::new(title).strong().color(t.text));
                    ui.allocate_ui(egui::vec2(width - 24.0, GEL_HEIGHT), |ui| {
                        crate::gel_view::draw_gel(ui, &t, lanes, None);
                    });
                    ui.add_space(8.0);
                    for (i, (lane, result)) in lanes.iter().enumerate() {
                        let mut text = format!("Lane {}: {lane}", i + 1);
                        if let Err(e) = result {
                            text.push_str(&format!(" (not run: {e})"));
                        }
                        ui.label(RichText::new(text).size(11.0).color(t.text));
                    }
                });
        },
    )?;
    if invert {
        for p in &mut image.pixels {
            let [r, g, b, a] = p.to_srgba_unmultiplied();
            *p = Color32::from_rgba_unmultiplied(255 - r, 255 - g, 255 - b, a);
        }
    }
    Ok(image)
}

/// Save an image as PNG or JPEG, by the file's extension (`.png`, `.jpg`, `.jpeg`).
#[cfg(not(target_arch = "wasm32"))]
pub fn save(path: &str, image: &ColorImage) -> Result<(), String> {
    let ext = std::path::Path::new(path)
        .extension()
        .map(|e| e.to_string_lossy().to_ascii_lowercase())
        .unwrap_or_default();
    match ext.as_str() {
        "png" => crate::offscreen::write_png(path, image),
        "jpg" | "jpeg" => write_jpeg(path, image),
        _ => Err(format!(
            "{path}: save the gel as .png or .jpg (the file name says which)"
        )),
    }
}

/// JPEG at quality 92: small files, no visible loss on bands and text.
#[cfg(not(target_arch = "wasm32"))]
fn write_jpeg(path: &str, image: &ColorImage) -> Result<(), String> {
    let [w, h] = image.size;
    let width = u32::try_from(w).map_err(|e| e.to_string())?;
    let height = u32::try_from(h).map_err(|e| e.to_string())?;
    let mut rgb = Vec::with_capacity(image.pixels.len().saturating_mul(3));
    for p in &image.pixels {
        let [r, g, b, _] = p.to_srgba_unmultiplied();
        rgb.extend_from_slice(&[r, g, b]);
    }
    let file = std::fs::File::create(path).map_err(|e| format!("cannot write {path}: {e}"))?;
    let mut out = std::io::BufWriter::new(file);
    image::codecs::jpeg::JpegEncoder::new_with_quality(&mut out, 92)
        .encode(&rgb, width, height, image::ExtendedColorType::Rgb8)
        .map_err(|e| format!("cannot write {path}: {e}"))
}

#[cfg(test)]
mod tests {
    use snappiesnap_engine::bio::gel::ladder;

    use super::*;

    fn lanes() -> Vec<RunLane> {
        vec![
            ("1 kb".to_string(), Ok(ladder("1 kb", 1.0).unwrap())),
            (
                "pSample: EcoRI + BamHI".to_string(),
                Ok(ladder("100 bp", 1.0).unwrap()),
            ),
            ("broken".to_string(), Err("no such enzyme".to_string())),
        ]
    }

    /// As shown the image is mostly dark, inverted mostly light; both save as PNG and JPEG.
    /// Set `SNAPPIESNAP_GEL_EXPORT_DIR=/dir` to look at them.
    #[test]
    fn images_and_files() {
        let brightness = |img: &ColorImage| {
            let sum: u64 = img
                .pixels
                .iter()
                .map(|p| u64::from(p.r()) + u64::from(p.g()) + u64::from(p.b()))
                .sum();
            sum / (img.pixels.len() as u64 * 3)
        };
        let shown = gel_image("Digest of pSample", &lanes(), false).unwrap();
        let inverted = gel_image("Digest of pSample", &lanes(), true).unwrap();
        assert_eq!(shown.size, inverted.size);
        assert!(brightness(&shown) < 80, "{}", brightness(&shown));
        assert!(brightness(&inverted) > 175, "{}", brightness(&inverted));
        assert!(gel_image("empty", &[], false).is_err());

        let dir = std::env::var("SNAPPIESNAP_GEL_EXPORT_DIR").map_or_else(
            |_| std::env::temp_dir().join(format!("snappiesnap-gel-{}", std::process::id())),
            std::path::PathBuf::from,
        );
        std::fs::create_dir_all(&dir).unwrap();
        for (name, img) in [("shown.png", &shown), ("inverted.jpg", &inverted)] {
            let path = dir.join(name);
            save(path.to_str().unwrap(), img).unwrap();
            let bytes = std::fs::read(&path).unwrap();
            let magic: &[u8] = if name.ends_with("png") {
                b"\x89PNG"
            } else {
                b"\xFF\xD8\xFF"
            };
            assert!(bytes.starts_with(magic), "{name}");
        }
        assert!(save(dir.join("x.gif").to_str().unwrap(), &shown).is_err());
        if std::env::var("SNAPPIESNAP_GEL_EXPORT_DIR").is_err() {
            let _ = std::fs::remove_dir_all(&dir);
        }
    }
}
