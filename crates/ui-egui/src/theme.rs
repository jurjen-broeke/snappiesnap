//! Design tokens. Every colour, radius and spacing in the UI comes from [`Tokens`]; nothing is
//! hard-coded in widgets (AGENTS.md › Architecture rules). The palette is SnappieSnap's own.

use eframe::egui::{self, Color32, CornerRadius, Stroke, Visuals};
use serde::{Deserialize, Serialize};

/// Which theme to use.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum ThemeChoice {
    /// Follow the operating system.
    #[default]
    System,
    Light,
    Dark,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Tokens {
    pub dark: bool,
    /// Window and panel background.
    pub panel: Color32,
    /// Toolbar and status bar background.
    pub chrome: Color32,
    /// Background of the map and sequence canvases.
    pub canvas: Color32,
    pub text: Color32,
    pub text_muted: Color32,
    pub border: Color32,
    pub accent: Color32,
    pub accent_text: Color32,
    /// Selection highlight on canvases.
    pub selection: Color32,
    /// Backbone lines of maps.
    pub backbone: Color32,
    pub error: Color32,
    /// Something to weigh, not something wrong: a suggestion that needs an enzyme the lab
    /// hasn't got, a datum we couldn't use. Distinct from [`Tokens::error`], which marks a
    /// failure.
    pub warning: Color32,
    /// The insert in the cloning schematic, drawn against [`Tokens::backbone`] for the
    /// vector: the lab convention of a black vector and a red insert, so it reads the same
    /// way in both themes (the "black" line is the canvas foreground).
    pub insert: Color32,
    /// A mark the pointer is over, drawn to stand out from the accent-coloured marks around
    /// it: the cut of the enzyme whose name is hovered in the sequence view.
    pub highlight: Color32,
    /// A simulated agarose gel: its background, wells and stained bands (dark in both themes,
    /// like a gel on a transilluminator).
    pub gel: Color32,
    pub gel_well: Color32,
    pub gel_band: Color32,
    /// Feature colours, cycled by feature index when a feature has no colour of its own.
    pub feature_palette: [Color32; 8],
    /// Chromatogram channels A, C, G and T (the field's usual green, blue, black, red; G is
    /// light on the dark canvas).
    pub trace_channels: [Color32; 4],
    pub radius: u8,
    pub spacing: f32,
    pub toolbar_icon: f32,
}

impl Tokens {
    pub fn dark() -> Self {
        Tokens {
            dark: true,
            panel: Color32::from_rgb(0x22, 0x25, 0x2b),
            chrome: Color32::from_rgb(0x1b, 0x1d, 0x22),
            canvas: Color32::from_rgb(0x15, 0x17, 0x1b),
            text: Color32::from_rgb(0xe6, 0xe8, 0xeb),
            text_muted: Color32::from_rgb(0x9a, 0xa1, 0xab),
            border: Color32::from_rgb(0x34, 0x38, 0x42),
            accent: Color32::from_rgb(0x4f, 0x8c, 0xff),
            accent_text: Color32::WHITE,
            selection: Color32::from_rgba_unmultiplied(0x4f, 0x8c, 0xff, 0x55),
            backbone: Color32::from_rgb(0xd8, 0xdb, 0xe0),
            error: Color32::from_rgb(0xff, 0x6b, 0x6b),
            warning: Color32::from_rgb(0xe8, 0xb3, 0x39),
            insert: Color32::from_rgb(0xf0, 0x55, 0x55),
            highlight: Color32::from_rgb(0xff, 0x4f, 0xc4),
            gel: Color32::from_rgb(0x10, 0x12, 0x16),
            gel_well: Color32::from_rgb(0x2a, 0x2e, 0x36),
            gel_band: Color32::from_rgb(0xf3, 0xf6, 0xff),
            feature_palette: [
                Color32::from_rgb(0x5d, 0xa9, 0xe9),
                Color32::from_rgb(0x8f, 0xd1, 0x7f),
                Color32::from_rgb(0xf2, 0xb1, 0x4f),
                Color32::from_rgb(0xd9, 0x7b, 0xc9),
                Color32::from_rgb(0x6f, 0xd3, 0xc8),
                Color32::from_rgb(0xef, 0x86, 0x6b),
                Color32::from_rgb(0xb3, 0x9d, 0xf2),
                Color32::from_rgb(0xd6, 0xd0, 0x5e),
            ],
            trace_channels: [
                Color32::from_rgb(0x4c, 0xc2, 0x5a),
                Color32::from_rgb(0x5a, 0x9b, 0xff),
                Color32::from_rgb(0xd8, 0xdb, 0xe0),
                Color32::from_rgb(0xff, 0x5f, 0x5f),
            ],
            radius: 5,
            spacing: 6.0,
            toolbar_icon: 18.0,
        }
    }

    pub fn light() -> Self {
        Tokens {
            dark: false,
            panel: Color32::from_rgb(0xf5, 0xf6, 0xf8),
            chrome: Color32::from_rgb(0xea, 0xec, 0xf0),
            canvas: Color32::WHITE,
            text: Color32::from_rgb(0x1d, 0x21, 0x28),
            text_muted: Color32::from_rgb(0x5d, 0x66, 0x73),
            border: Color32::from_rgb(0xd0, 0xd5, 0xdc),
            accent: Color32::from_rgb(0x2f, 0x6f, 0xe4),
            accent_text: Color32::WHITE,
            selection: Color32::from_rgba_unmultiplied(0x2f, 0x6f, 0xe4, 0x40),
            backbone: Color32::from_rgb(0x3a, 0x3f, 0x47),
            error: Color32::from_rgb(0xc6, 0x28, 0x28),
            warning: Color32::from_rgb(0x8a, 0x5a, 0x00),
            insert: Color32::from_rgb(0xcc, 0x22, 0x22),
            highlight: Color32::from_rgb(0xc2, 0x18, 0x8a),
            gel: Color32::from_rgb(0x1a, 0x1c, 0x22),
            gel_well: Color32::from_rgb(0x3a, 0x3f, 0x4a),
            gel_band: Color32::from_rgb(0xf3, 0xf6, 0xff),
            feature_palette: [
                Color32::from_rgb(0x2b, 0x7b, 0xc4),
                Color32::from_rgb(0x3f, 0x9a, 0x3a),
                Color32::from_rgb(0xd0, 0x85, 0x12),
                Color32::from_rgb(0xae, 0x3f, 0x9d),
                Color32::from_rgb(0x1f, 0x9a, 0x8e),
                Color32::from_rgb(0xc9, 0x4f, 0x2e),
                Color32::from_rgb(0x6b, 0x55, 0xc4),
                Color32::from_rgb(0x9a, 0x92, 0x10),
            ],
            trace_channels: [
                Color32::from_rgb(0x1e, 0x9a, 0x2c),
                Color32::from_rgb(0x1f, 0x5f, 0xd6),
                Color32::from_rgb(0x20, 0x22, 0x26),
                Color32::from_rgb(0xd6, 0x26, 0x26),
            ],
            radius: 5,
            spacing: 6.0,
            toolbar_icon: 18.0,
        }
    }

    pub fn for_choice(choice: ThemeChoice, system_dark: bool) -> Self {
        match choice {
            ThemeChoice::Dark => Tokens::dark(),
            ThemeChoice::Light => Tokens::light(),
            ThemeChoice::System if system_dark => Tokens::dark(),
            ThemeChoice::System => Tokens::light(),
        }
    }

    /// The colour of a feature: its own stored colour (set with the colour picker, or read
    /// from a file), else one chosen from its name and type. Either way the same feature has
    /// the same colour in every molecule it is carried into (clones, PCR products, history
    /// steps). The choice uses FNV-1a, so it is the same on every platform and run.
    pub fn color_of(&self, f: &snappiesnap_engine::doc::Feature) -> Color32 {
        if let Some(c) = f.color.as_deref().and_then(parse_hex) {
            return c;
        }
        self.automatic_color_of(f)
    }

    /// The colour a feature gets when it has none of its own (see [`Tokens::color_of`]).
    pub fn automatic_color_of(&self, f: &snappiesnap_engine::doc::Feature) -> Color32 {
        self.automatic_color(&f.name, &f.kind)
    }

    /// The automatic colour of a feature with this name and type.
    pub fn automatic_color(&self, name: &str, kind: &str) -> Color32 {
        let mut h: u64 = 0xcbf2_9ce4_8422_2325;
        for b in name.bytes().chain([0]).chain(kind.bytes()) {
            h ^= u64::from(b);
            h = h.wrapping_mul(0x0000_0100_0000_01b3);
        }
        self.feature_color((h % 1_000_003) as usize)
    }

    /// The colour of feature `index` of `m` (see [`Tokens::color_of`]).
    pub fn color_at(&self, m: &snappiesnap_engine::doc::Molecule, index: usize) -> Color32 {
        m.features()
            .get(index)
            .map_or(self.accent, |f| self.color_of(f))
    }

    pub fn feature_color(&self, index: usize) -> Color32 {
        self.feature_palette
            .get(index % self.feature_palette.len())
            .copied()
            .unwrap_or(self.accent)
    }

    /// Apply the tokens to egui's widget visuals.
    pub fn apply(&self, ctx: &egui::Context) {
        let mut v = if self.dark {
            Visuals::dark()
        } else {
            Visuals::light()
        };
        let r = CornerRadius::same(self.radius);
        v.panel_fill = self.panel;
        v.window_fill = self.panel;
        v.extreme_bg_color = self.canvas;
        v.faint_bg_color = self.chrome;
        v.override_text_color = Some(self.text);
        v.selection.bg_fill = self.accent;
        v.selection.stroke = Stroke::new(1.0, self.accent_text);
        v.hyperlink_color = self.accent;
        v.error_fg_color = self.error;
        v.window_corner_radius = r;
        v.menu_corner_radius = r;
        for w in [
            &mut v.widgets.noninteractive,
            &mut v.widgets.inactive,
            &mut v.widgets.hovered,
            &mut v.widgets.active,
            &mut v.widgets.open,
        ] {
            w.corner_radius = r;
        }
        v.widgets.noninteractive.bg_stroke = Stroke::new(1.0, self.border);
        ctx.set_visuals(v);
    }
}

/// A colour written `#rrggbb` (case-insensitive), as feature colours and preferences store it.
pub fn parse_hex(text: &str) -> Option<Color32> {
    let digits = text.trim().strip_prefix('#')?;
    if digits.len() != 6 || !digits.chars().all(|c| c.is_ascii_hexdigit()) {
        return None;
    }
    let byte = |i: usize| {
        digits
            .get(i..i + 2)
            .and_then(|h| u8::from_str_radix(h, 16).ok())
    };
    Some(Color32::from_rgb(byte(0)?, byte(2)?, byte(4)?))
}

/// A colour as `#rrggbb` (lower case; alpha is dropped).
pub fn hex(c: Color32) -> String {
    format!("#{:02x}{:02x}{:02x}", c.r(), c.g(), c.b())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hex_colours_round_trip() {
        assert_eq!(parse_hex("#FF8000"), Some(Color32::from_rgb(255, 128, 0)));
        assert_eq!(parse_hex(" #0a0B0c "), Some(Color32::from_rgb(10, 11, 12)));
        assert_eq!(hex(Color32::from_rgb(10, 11, 12)), "#0a0b0c");
        for bad in ["", "#", "#12345", "#1234567", "123456", "#12345g", "#ééé"] {
            assert_eq!(parse_hex(bad), None, "{bad}");
        }
    }

    #[test]
    fn a_stored_colour_wins() {
        use snappiesnap_engine::doc::{FeatureId, NewFeature, Range, Segment, Strand};
        let mut f = NewFeature {
            color: Some("#123456".into()),
            name: "AmpR".into(),
            kind: "CDS".into(),
            strand: Strand::Forward,
            segments: vec![Segment {
                range: Range::new(0, 4),
            }],
            qualifiers: Vec::new(),
        }
        .with_id(FeatureId(1));
        let t = Tokens::dark();
        assert_eq!(t.color_of(&f), Color32::from_rgb(0x12, 0x34, 0x56));
        f.color = None;
        assert_eq!(t.color_of(&f), t.automatic_color_of(&f));
    }

    #[test]
    fn a_feature_keeps_its_colour_wherever_it_goes() {
        use snappiesnap_engine::doc::{
            Alphabet, FeatureId, Molecule, NewFeature, Range, Segment, Seq, Strand, Topology,
        };
        let feature = |name: &str, kind: &str, at: usize| NewFeature {
            color: None,
            name: name.into(),
            kind: kind.into(),
            strand: Strand::Forward,
            segments: vec![Segment {
                range: Range::new(at, at + 4),
            }],
            qualifiers: Vec::new(),
        };
        let seq = Seq::new(Alphabet::Dna, "ACGT".repeat(10)).unwrap();
        let mut parent = Molecule::new("parent", seq.clone(), Topology::Circular).unwrap();
        parent.add_feature(feature("AmpR", "CDS", 0)).unwrap();
        // In the product it is the third feature, not the first.
        let mut product = Molecule::new("product", seq, Topology::Linear).unwrap();
        product
            .add_feature(feature("ori", "rep_origin", 0))
            .unwrap();
        product.add_feature(feature("lacZ", "CDS", 8)).unwrap();
        product.add_feature(feature("AmpR", "CDS", 20)).unwrap();
        for t in [Tokens::light(), Tokens::dark()] {
            assert_eq!(t.color_at(&parent, 0), t.color_at(&product, 2));
            // Fixed values: the same on every platform and in every run.
            assert_eq!(t.color_at(&product, 2), t.color_of(&product.features()[2]));
        }
        let t = Tokens::light();
        let colours: std::collections::HashSet<_> =
            ["ori", "lacZ", "AmpR", "CMV", "SV40", "f1 ori"]
                .iter()
                .map(|n| t.color_of(&feature(n, "misc_feature", 0).with_id(FeatureId(1))))
                .collect();
        assert!(colours.len() >= 4, "names spread over the palette");
    }

    #[test]
    fn themes_differ_and_choose_correctly() {
        assert!(Tokens::dark().dark && !Tokens::light().dark);
        assert_eq!(
            Tokens::for_choice(ThemeChoice::System, true),
            Tokens::dark()
        );
        assert_eq!(
            Tokens::for_choice(ThemeChoice::Light, true),
            Tokens::light()
        );
        let t = Tokens::dark();
        assert_eq!(t.feature_color(0), t.feature_color(8));
    }
}
