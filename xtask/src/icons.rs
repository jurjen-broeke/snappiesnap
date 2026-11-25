//! `cargo xtask icons`: render the application icon pack from the brand drawing.
//!
//! The master is `assets/brand/snappiesnap.svg`; everything under `assets/brand/icons/` is
//! generated from it and is rewritten on every run (commit the result):
//!
//! - `linux/hicolor/<n>x<n>/apps/snappiesnap.png` and `scalable/apps/snappiesnap.svg`: the
//!   freedesktop icon theme layout, installed under `/usr/share/icons/hicolor/`.
//! - `windows/snappiesnap.ico`: embedded in `snappiesnap.exe` and used by both installers.
//! - `macos/SnappieSnap.icns`: the bundle icon (`CFBundleIconFile`). macOS icons sit on a
//!   rounded-square tile, so these get one; `macos/snappiesnap-512.png` is the same tile for the
//!   Dock when the app runs outside a bundle (`cargo run`).
//! - `web/`: favicons, PWA icons (plain and maskable), the Apple touch icon and a web manifest.
//!
//! The drawing is used as is, with no tile, everywhere but macOS and the opaque web icons.
//! Rendering is pure Rust (resvg) so the pack can be rebuilt on any platform with the same result.
//! `.ico` (Microsoft ICONDIR with PNG entries) and `.icns` (Apple icon family with PNG entries)
//! are written here; both are simple containers.

use std::fs;
use std::path::Path;

use resvg::tiny_skia::{
    self, Color, FillRule, GradientStop, LinearGradient, Paint, Path as SkPath, PathBuilder,
    Pixmap, Point, SpreadMode, Transform,
};
use resvg::usvg;
use serde_json::json;

use crate::util::{self, Result};

const MASTER: &str = "assets/brand/snappiesnap.svg";
const OUT: &str = "assets/brand/icons";

/// freedesktop hicolor sizes.
const LINUX: &[u32] = &[16, 22, 24, 32, 48, 64, 128, 256, 512];
/// Windows shell sizes at 100–200 % display scaling, up to the 256 px Explorer tile.
const WINDOWS: &[u32] = &[16, 20, 24, 32, 40, 48, 64, 96, 128, 256];
/// `.icns` entries as (OSType, pixel size): 16–512 pt at 1x and 2x.
const MACOS: &[([u8; 4], u32)] = &[
    (*b"icp4", 16),
    (*b"icp5", 32),
    (*b"ic11", 32),
    (*b"ic12", 64),
    (*b"ic07", 128),
    (*b"ic13", 256),
    (*b"ic08", 256),
    (*b"ic14", 512),
    (*b"ic09", 512),
    (*b"ic10", 1024),
];
const FAVICON: &[u32] = &[16, 32, 48];

/// The tile behind the drawing on macOS and on opaque web icons: a light water blue, top to bottom.
const TILE_TOP: (u8, u8, u8) = (0xEC, 0xF8, 0xFF);
const TILE_BOTTOM: (u8, u8, u8) = (0xBF, 0xE2, 0xF4);
/// The drawing's head green, for the web manifest's theme colour.
const THEME: &str = "#18bc18";

pub fn run() -> Result {
    let root = util::root()?;
    let master = root.join(MASTER);
    let svg = fs::read(&master).map_err(|e| format!("cannot read {}: {e}", master.display()))?;
    let tree = usvg::Tree::from_data(&svg, &usvg::Options::default())
        .map_err(|e| format!("{MASTER}: {e}"))?;
    let out = root.join(OUT);
    if out.exists() {
        // Start clean, so a size dropped from the lists above doesn't linger.
        fs::remove_dir_all(&out).map_err(|e| format!("cannot clear {}: {e}", out.display()))?;
    }

    for &size in LINUX {
        let path = format!("linux/hicolor/{size}x{size}/apps/snappiesnap.png");
        save(&out.join(path), &png(&plain(&tree, size)?)?)?;
    }
    save(
        &out.join("linux/hicolor/scalable/apps/snappiesnap.svg"),
        &svg,
    )?;

    save(
        &out.join("windows/snappiesnap.ico"),
        &ico(&rendered(WINDOWS, |s| plain(&tree, s))?)?,
    )?;

    let mut entries = Vec::new();
    for &(kind, size) in MACOS {
        entries.push((kind, png(&tile(&tree, size)?)?));
    }
    save(&out.join("macos/SnappieSnap.icns"), &icns(&entries)?)?;
    save(
        &out.join("macos/snappiesnap-512.png"),
        &png(&tile(&tree, 512)?)?,
    )?;

    let web = out.join("web");
    save(
        &web.join("favicon.ico"),
        &ico(&rendered(FAVICON, |s| plain(&tree, s))?)?,
    )?;
    save(&web.join("favicon.svg"), &svg)?;
    save(&web.join("icon-192.png"), &png(&plain(&tree, 192)?)?)?;
    save(&web.join("icon-512.png"), &png(&plain(&tree, 512)?)?)?;
    // iOS rounds the corners itself and shows transparency as black: full-bleed tile.
    save(
        &web.join("apple-touch-icon.png"),
        &png(&full_bleed(&tree, 180, 0.80)?)?,
    )?;
    // Maskable icons may be cropped to a circle of 80 % of the width (W3C manifest spec).
    save(
        &web.join("icon-maskable-512.png"),
        &png(&full_bleed(&tree, 512, 0.72)?)?,
    )?;
    save(&web.join("manifest.webmanifest"), manifest()?.as_bytes())?;

    println!("icons: wrote the icon pack to {OUT}/ from {MASTER}");
    Ok(())
}

/// The drawing alone on a transparent square, with a margin of 1/32 of the size.
fn plain(tree: &usvg::Tree, size: u32) -> Result<Pixmap> {
    let mut pixmap = canvas(size)?;
    let s = size as f32;
    draw(tree, &mut pixmap, s / 2.0, s / 2.0, s * 30.0 / 32.0);
    Ok(pixmap)
}

/// The drawing on a macOS-style tile: Apple's icon grid puts an 824 px rounded square at
/// (100, 100) on a 1024 px canvas, with a soft shadow below it.
fn tile(tree: &usvg::Tree, size: u32) -> Result<Pixmap> {
    let mut pixmap = canvas(size)?;
    let s = size as f32;
    let body = s * 824.0 / 1024.0;
    let centre = s / 2.0;
    // Shadow: stacked, slightly larger copies of the tile at low opacity stand in for a blur.
    const LAYERS: u8 = 12;
    let blur = s * 24.0 / 1024.0;
    let drop = s * 10.0 / 1024.0;
    for i in 1..=LAYERS {
        let grow = blur * f32::from(i) / f32::from(LAYERS);
        if let Some(path) = squircle(centre, centre + drop, body / 2.0 + grow) {
            let mut paint = Paint::default();
            paint.set_color_rgba8(0, 0, 0, 7);
            paint.anti_alias = true;
            pixmap.fill_path(
                &path,
                &paint,
                FillRule::Winding,
                Transform::identity(),
                None,
            );
        }
    }
    let shape = squircle(centre, centre, body / 2.0).ok_or("cannot build the tile shape")?;
    let top = centre - body / 2.0;
    pixmap.fill_path(
        &shape,
        &gradient(top, top + body)?,
        FillRule::Winding,
        Transform::identity(),
        None,
    );
    draw(tree, &mut pixmap, centre, centre, body * 0.84);
    Ok(pixmap)
}

/// The drawing on a square tile that fills the canvas, scaled to `extent` of its width.
fn full_bleed(tree: &usvg::Tree, size: u32, extent: f32) -> Result<Pixmap> {
    let mut pixmap = canvas(size)?;
    let s = size as f32;
    let rect = tiny_skia::Rect::from_xywh(0.0, 0.0, s, s).ok_or("empty icon")?;
    pixmap.fill_rect(rect, &gradient(0.0, s)?, Transform::identity(), None);
    draw(tree, &mut pixmap, s / 2.0, s / 2.0, s * extent);
    Ok(pixmap)
}

fn canvas(size: u32) -> Result<Pixmap> {
    Pixmap::new(size, size).ok_or_else(|| format!("cannot allocate a {size} px icon"))
}

/// Draw the whole drawing centred on (`cx`, `cy`), its longer side `extent` pixels.
fn draw(tree: &usvg::Tree, pixmap: &mut Pixmap, cx: f32, cy: f32, extent: f32) {
    let size = tree.size();
    let scale = extent / size.width().max(size.height());
    let transform = Transform::from_translate(
        cx - size.width() * scale / 2.0,
        cy - size.height() * scale / 2.0,
    )
    .pre_scale(scale, scale);
    resvg::render(tree, transform, &mut pixmap.as_mut());
}

fn gradient(top: f32, bottom: f32) -> Result<Paint<'static>> {
    let rgb = |(r, g, b): (u8, u8, u8)| Color::from_rgba8(r, g, b, 255);
    let shader = LinearGradient::new(
        Point::from_xy(0.0, top),
        Point::from_xy(0.0, bottom),
        vec![
            GradientStop::new(0.0, rgb(TILE_TOP)),
            GradientStop::new(1.0, rgb(TILE_BOTTOM)),
        ],
        SpreadMode::Pad,
        Transform::identity(),
    )
    .ok_or("cannot build the tile gradient")?;
    Ok(Paint {
        shader,
        anti_alias: true,
        ..Default::default()
    })
}

/// A superellipse |x|^5 + |y|^5 = r^5: close to the continuous-corner rounded square of
/// macOS app icons.
fn squircle(cx: f32, cy: f32, r: f32) -> Option<SkPath> {
    const STEPS: u16 = 256;
    const EXPONENT: f32 = 2.0 / 5.0;
    let mut pb = PathBuilder::new();
    for i in 0..STEPS {
        let t = f32::from(i) / f32::from(STEPS) * std::f32::consts::TAU;
        let (sin, cos) = t.sin_cos();
        let x = cx + r * cos.signum() * cos.abs().powf(EXPONENT);
        let y = cy + r * sin.signum() * sin.abs().powf(EXPONENT);
        if i == 0 {
            pb.move_to(x, y);
        } else {
            pb.line_to(x, y);
        }
    }
    pb.close();
    pb.finish()
}

fn rendered(sizes: &[u32], render: impl Fn(u32) -> Result<Pixmap>) -> Result<Vec<(u32, Vec<u8>)>> {
    sizes.iter().map(|&s| Ok((s, png(&render(s)?)?))).collect()
}

fn png(pixmap: &Pixmap) -> Result<Vec<u8>> {
    pixmap
        .encode_png()
        .map_err(|e| format!("cannot encode PNG: {e}"))
}

/// A Windows `.ico`: ICONDIR header, one 16-byte ICONDIRENTRY per image, then the PNG data
/// (PNG entries are supported since Windows Vista). A width or height of 256 is stored as 0.
fn ico(images: &[(u32, Vec<u8>)]) -> Result<Vec<u8>> {
    let count = u16::try_from(images.len()).map_err(|_| "too many icon images")?;
    let mut out = Vec::new();
    out.extend_from_slice(&0u16.to_le_bytes()); // reserved
    out.extend_from_slice(&1u16.to_le_bytes()); // type: icon
    out.extend_from_slice(&count.to_le_bytes());
    let mut offset = 6 + 16 * images.len();
    for (size, data) in images {
        if !(1..=256).contains(size) {
            return Err(format!("an .ico image can't be {size} px"));
        }
        let side = u8::try_from(*size).unwrap_or(0); // 256 is written as 0
        out.extend_from_slice(&[side, side, 0, 0]); // width, height, palette, reserved
        out.extend_from_slice(&1u16.to_le_bytes()); // colour planes
        out.extend_from_slice(&32u16.to_le_bytes()); // bits per pixel
        out.extend_from_slice(&len32(data.len())?.to_le_bytes());
        out.extend_from_slice(&len32(offset)?.to_le_bytes());
        offset = offset.saturating_add(data.len());
    }
    for (_, data) in images {
        out.extend_from_slice(data);
    }
    Ok(out)
}

/// An Apple `.icns`: the `icns` magic and total length, then (OSType, length including the
/// 8-byte header, PNG data) per entry. All lengths are big-endian.
fn icns(entries: &[([u8; 4], Vec<u8>)]) -> Result<Vec<u8>> {
    let mut body = Vec::new();
    for (kind, data) in entries {
        body.extend_from_slice(kind);
        body.extend_from_slice(&len32(data.len().saturating_add(8))?.to_be_bytes());
        body.extend_from_slice(data);
    }
    let mut out = Vec::with_capacity(body.len().saturating_add(8));
    out.extend_from_slice(b"icns");
    out.extend_from_slice(&len32(body.len().saturating_add(8))?.to_be_bytes());
    out.extend_from_slice(&body);
    Ok(out)
}

fn len32(n: usize) -> Result<u32> {
    u32::try_from(n).map_err(|_| "icon file too large".to_string())
}

fn manifest() -> Result<String> {
    let hex = |(r, g, b): (u8, u8, u8)| format!("#{r:02x}{g:02x}{b:02x}");
    let manifest = json!({
        "name": "SnappieSnap",
        "short_name": "SnappieSnap",
        "start_url": ".",
        "display": "standalone",
        "theme_color": THEME,
        "background_color": hex(TILE_TOP),
        "icons": [
            { "src": "favicon.svg", "sizes": "any", "type": "image/svg+xml" },
            { "src": "icon-192.png", "sizes": "192x192", "type": "image/png" },
            { "src": "icon-512.png", "sizes": "512x512", "type": "image/png" },
            { "src": "icon-maskable-512.png", "sizes": "512x512", "type": "image/png", "purpose": "maskable" },
        ],
    });
    serde_json::to_string_pretty(&manifest)
        .map(|s| s + "\n")
        .map_err(|e| format!("web manifest: {e}"))
}

fn save(path: &Path, bytes: &[u8]) -> Result {
    if let Some(dir) = path.parent() {
        fs::create_dir_all(dir).map_err(|e| format!("cannot create {}: {e}", dir.display()))?;
    }
    fs::write(path, bytes).map_err(|e| format!("cannot write {}: {e}", path.display()))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ico_directory_points_at_each_image() {
        let images = vec![(16, vec![1, 2, 3]), (256, vec![4, 5])];
        let ico = ico(&images).unwrap();
        assert_eq!(&ico[..6], &[0, 0, 1, 0, 2, 0]);
        // First entry: 16 px, 3 bytes at offset 6 + 2 * 16.
        assert_eq!(ico[6], 16);
        assert_eq!(&ico[14..18], &3u32.to_le_bytes());
        assert_eq!(&ico[18..22], &38u32.to_le_bytes());
        // Second entry: 256 px is stored as 0, its data follows the first image.
        assert_eq!(ico[22], 0);
        assert_eq!(&ico[34..38], &41u32.to_le_bytes());
        assert_eq!(&ico[38..], &[1, 2, 3, 4, 5]);
        assert!(super::ico(&[(512, vec![])]).is_err());
    }

    #[test]
    fn icns_lengths_include_headers() {
        let icns = icns(&[(*b"ic07", vec![9; 4])]).unwrap();
        assert_eq!(&icns[..4], b"icns");
        assert_eq!(&icns[4..8], &20u32.to_be_bytes());
        assert_eq!(&icns[8..12], b"ic07");
        assert_eq!(&icns[12..16], &12u32.to_be_bytes());
        assert_eq!(icns.len(), 20);
    }

    #[test]
    fn squircle_fits_its_box() {
        let path = squircle(50.0, 50.0, 40.0).unwrap();
        let b = path.bounds();
        assert!((b.left() - 10.0).abs() < 0.01 && (b.right() - 90.0).abs() < 0.01);
        assert!((b.top() - 10.0).abs() < 0.01 && (b.bottom() - 90.0).abs() < 0.01);
    }
}
