//! Settings remembered between sessions: theme, view, map style, the folder panel's folder,
//! visibility and expanded folders. Stored as JSON in a folder of the user's own
//! (`snappiesnap_engine::paths`; override with an absolute `SNAPPIESNAP_SETTINGS=/path/settings.json`). A missing file means defaults; an
//! unreadable or corrupt one is reported and ignored, never a crash. Writes are atomic.

use std::path::PathBuf;

use serde::{Deserialize, Serialize};

use crate::enzyme_view::ListFrom;
use crate::preferences::PrefsPage;
use crate::state::{CanvasColors, EnzymeSet, MapStyle, Translation, View};
use crate::theme::ThemeChoice;
use snappiesnap_engine::bio::primers::BindingPrefs;

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(default)]
pub struct Settings {
    pub theme: ThemeChoice,
    pub view: View,
    pub map_style: MapStyle,
    pub folder: Option<String>,
    pub sidebar: bool,
    pub expanded: Vec<String>,
    /// Map layers (see `UiState`).
    pub features: bool,
    pub enzymes: bool,
    pub primers: bool,
    pub positions: bool,
    pub ruler: bool,
    pub orfs: bool,
    /// Which enzymes the map and the sequence view show (including the user's enzyme set). A
    /// broken value falls back to the default without losing the other settings.
    #[serde(deserialize_with = "lenient")]
    pub enzyme_set: EnzymeSet,
    /// What the Enzymes view lists: `"all"`, `{"set": "name"}`, or nothing yet (the map's
    /// set). A broken value means nothing yet.
    #[serde(deserialize_with = "lenient_list")]
    pub enzyme_list: Option<ListFrom>,
    /// The map and sequence background per theme (Preferences). A broken value means the
    /// themes' own.
    #[serde(deserialize_with = "lenient_background")]
    pub map_background: CanvasColors,
    /// The category Preferences last showed. One this version doesn't know means General.
    #[serde(deserialize_with = "lenient_default")]
    pub preferences_page: PrefsPage,
    /// Interface size (1 = normal), kept within Preferences' range.
    #[serde(deserialize_with = "lenient_scale")]
    pub ui_scale: f32,
    /// Size of the residues in the sequence view and the Align with pane, in points, kept
    /// within Preferences' range.
    #[serde(deserialize_with = "lenient_sequence_font")]
    pub sequence_font: f32,
    /// Reopen the window at its last size.
    #[serde(deserialize_with = "lenient_restore")]
    pub restore_window: bool,
    /// Look for a new version once at startup (installed builds only). A broken value means
    /// yes, like a missing one.
    #[serde(deserialize_with = "lenient_restore")]
    pub check_updates: bool,
    /// What the ORF layer shows in the sequence view, and the genetic code of translations
    /// (an NCBI table number; one that doesn't exist means the standard code).
    #[serde(deserialize_with = "lenient_default")]
    pub translation: Translation,
    #[serde(deserialize_with = "lenient_code")]
    pub genetic_code: u8,
    /// How primer tails show in the sequence view. One this version doesn't know means the
    /// default.
    #[serde(deserialize_with = "lenient_default")]
    pub primer_tails: crate::state::TailStyle,
    /// Which primer binding sites count (Preferences > Primers). A broken value means the
    /// defaults; values out of range are brought into range when used.
    #[serde(deserialize_with = "lenient_default")]
    pub primer_binding: BindingPrefs,
    /// The window's size when the app last ran. A broken value means none.
    #[serde(deserialize_with = "lenient_default")]
    pub window: Option<WindowSize>,
}

/// The size of the main window, in points (before the interface size is applied).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct WindowSize {
    pub width: f32,
    pub height: f32,
    #[serde(default)]
    pub maximized: bool,
}

impl WindowSize {
    /// Smallest and largest size we reopen at; the window's minimum is 640 x 400.
    const MIN: [f32; 2] = [640.0, 400.0];
    const MAX: f32 = 16_384.0;

    /// The size to open at, or `None` when the stored one isn't usable.
    pub fn usable(self) -> Option<[f32; 2]> {
        let ok = |v: f32| v.is_finite() && v > 0.0;
        (ok(self.width) && ok(self.height)).then(|| {
            [
                self.width.clamp(Self::MIN[0], Self::MAX),
                self.height.clamp(Self::MIN[1], Self::MAX),
            ]
        })
    }
}

/// A value of `T`, or `T`'s default when the stored one doesn't parse, without losing the
/// other settings.
fn lenient_default<'de, D: serde::Deserializer<'de>, T: serde::de::DeserializeOwned + Default>(
    d: D,
) -> Result<T, D::Error> {
    let value = serde_json::Value::deserialize(d)?;
    Ok(serde_json::from_value(value).unwrap_or_default())
}

fn lenient_scale<'de, D: serde::Deserializer<'de>>(d: D) -> Result<f32, D::Error> {
    let value = serde_json::Value::deserialize(d)?;
    Ok(value
        .as_f64()
        .map_or(1.0, |z| crate::state::clamp_ui_scale(z as f32)))
}

fn lenient_sequence_font<'de, D: serde::Deserializer<'de>>(d: D) -> Result<f32, D::Error> {
    let value = serde_json::Value::deserialize(d)?;
    Ok(value
        .as_f64()
        .map_or(crate::state::SEQUENCE_FONT_DEFAULT, |z| {
            crate::state::clamp_sequence_font(z as f32)
        }))
}

fn lenient_code<'de, D: serde::Deserializer<'de>>(d: D) -> Result<u8, D::Error> {
    let value = serde_json::Value::deserialize(d)?;
    Ok(value
        .as_u64()
        .and_then(|c| u8::try_from(c).ok())
        .filter(|&c| snappiesnap_engine::bio::translate::code(c).is_some())
        .unwrap_or(1))
}

fn lenient_restore<'de, D: serde::Deserializer<'de>>(d: D) -> Result<bool, D::Error> {
    let value = serde_json::Value::deserialize(d)?;
    Ok(value.as_bool().unwrap_or(true))
}

fn lenient_background<'de, D: serde::Deserializer<'de>>(d: D) -> Result<CanvasColors, D::Error> {
    let value = serde_json::Value::deserialize(d)?;
    Ok(serde_json::from_value::<CanvasColors>(value)
        .map(CanvasColors::clamped)
        .unwrap_or_default())
}

fn lenient_list<'de, D: serde::Deserializer<'de>>(d: D) -> Result<Option<ListFrom>, D::Error> {
    let value = serde_json::Value::deserialize(d)?;
    Ok(serde_json::from_value::<ListFrom>(value)
        .ok()
        .map(|l| match l {
            ListFrom::Set(name) => ListFrom::Set(
                name.trim()
                    .chars()
                    .take(snappiesnap_engine::enzymes::MAX_SET_NAME)
                    .collect(),
            ),
            all => all,
        }))
}

fn lenient<'de, D: serde::Deserializer<'de>>(d: D) -> Result<EnzymeSet, D::Error> {
    let value = serde_json::Value::deserialize(d)?;
    Ok(serde_json::from_value::<EnzymeSet>(value)
        .map(EnzymeSet::clamped)
        .unwrap_or_default())
}

impl Default for Settings {
    fn default() -> Self {
        Settings {
            theme: ThemeChoice::default(),
            view: View::default(),
            map_style: MapStyle::default(),
            folder: None,
            sidebar: true,
            expanded: Vec::new(),
            features: true,
            enzymes: true,
            primers: true,
            positions: true,
            ruler: true,
            orfs: false,
            enzyme_set: EnzymeSet::default(),
            enzyme_list: None,
            map_background: CanvasColors::default(),
            preferences_page: PrefsPage::default(),
            ui_scale: 1.0,
            sequence_font: crate::state::SEQUENCE_FONT_DEFAULT,
            restore_window: true,
            check_updates: true,
            translation: Translation::default(),
            genetic_code: 1,
            primer_tails: crate::state::TailStyle::default(),
            primer_binding: BindingPrefs::default(),
            window: None,
        }
    }
}

/// Largest settings file we read.
const MAX_BYTES: u64 = 1 << 20;

/// Where settings live on this platform: a folder of the user's own (see
/// `snappiesnap_engine::paths`), never the application's folder.
pub fn default_path() -> Option<PathBuf> {
    snappiesnap_engine::paths::settings_file()
}

/// Load settings. Returns defaults and a message if the file can't be used.
pub fn load(path: &std::path::Path) -> (Settings, Option<String>) {
    use std::io::Read;
    let file = match std::fs::File::open(path) {
        Ok(f) => f,
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => return (Settings::default(), None),
        Err(e) => {
            return (
                Settings::default(),
                Some(format!("settings not loaded: {e}")),
            );
        }
    };
    let mut text = String::new();
    if let Err(e) = file.take(MAX_BYTES).read_to_string(&mut text) {
        return (
            Settings::default(),
            Some(format!("settings not loaded: {e}")),
        );
    }
    match serde_json::from_str(&text) {
        Ok(s) => (s, None),
        Err(e) => (
            Settings::default(),
            Some(format!("settings file ignored ({e}); defaults are used")),
        ),
    }
}

/// Save settings atomically (temporary file, then rename), creating the folder if needed.
pub fn save(path: &std::path::Path, settings: &Settings) -> std::io::Result<()> {
    if let Some(dir) = path.parent() {
        std::fs::create_dir_all(dir)?;
    }
    let text = serde_json::to_string_pretty(settings).map_err(std::io::Error::other)?;
    let temp = path.with_extension("json.tmp");
    std::fs::write(&temp, text)?;
    std::fs::rename(&temp, path)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn round_trip_and_bad_files() {
        let dir = std::env::temp_dir().join(format!("snappiesnap-settings-{}", std::process::id()));
        let path = dir.join("nested/settings.json");
        assert_eq!(load(&path), (Settings::default(), None));
        let s = Settings {
            theme: ThemeChoice::Light,
            view: View::Sequence,
            folder: Some("/x/y".into()),
            sidebar: false,
            expanded: vec!["/x/y/z".into()],
            enzyme_set: EnzymeSet {
                pool: Some("Lab stock".into()),
                max_cuts: 2,
                ..EnzymeSet::default()
            },
            enzyme_list: Some(ListFrom::All),
            map_background: CanvasColors {
                light: Some("#fffff0".into()),
                dark: None,
            },
            ..Settings::default()
        };
        save(&path, &s).unwrap();
        assert_eq!(load(&path), (s, None));
        std::fs::write(&path, b"{not json").unwrap();
        let (loaded, msg) = load(&path);
        assert_eq!(loaded, Settings::default());
        assert!(msg.unwrap().contains("ignored"));
        // Unknown fields from newer versions are fine.
        std::fs::write(&path, br#"{"theme": "dark", "future": 1}"#).unwrap();
        assert_eq!(load(&path).0.theme, ThemeChoice::Dark);
        // A broken enzyme choice is reset on its own.
        std::fs::write(
            &path,
            br#"{"theme": "dark", "enzyme_set": {"max_cuts": "x"}}"#,
        )
        .unwrap();
        let (loaded, msg) = load(&path);
        assert_eq!((loaded.theme, msg), (ThemeChoice::Dark, None));
        assert_eq!(loaded.enzyme_set, EnzymeSet::default());
        // The Enzymes view's list: a set by name, or a broken value meaning no choice yet.
        std::fs::write(&path, br#"{"enzyme_list": {"set": " Lab "}}"#).unwrap();
        assert_eq!(load(&path).0.enzyme_list, Some(ListFrom::Set("Lab".into())));
        std::fs::write(&path, br#"{"theme": "dark", "enzyme_list": 7}"#).unwrap();
        let (loaded, msg) = load(&path);
        assert_eq!(
            (loaded.theme, loaded.enzyme_list, msg),
            (ThemeChoice::Dark, None, None)
        );
        // Preferences: remembered; a page or window from a newer version, or an absurd
        // size, falls back on its own.
        let s = Settings {
            preferences_page: PrefsPage::Enzymes,
            ui_scale: 1.25,
            sequence_font: 17.0,
            restore_window: false,
            check_updates: false,
            translation: Translation::Six,
            genetic_code: 11,
            primer_binding: BindingPrefs {
                long_from_bp: 50_000,
                long: snappiesnap_engine::bio::primers::ClassRules {
                    matched: 20,
                    isolated_mismatches: false,
                    min_tm_on: false,
                    min_tm: 52.5,
                },
                ..BindingPrefs::default()
            },
            window: Some(WindowSize {
                width: 1400.0,
                height: 900.0,
                maximized: true,
            }),
            ..Settings::default()
        };
        save(&path, &s).unwrap();
        assert_eq!(load(&path), (s, None));
        std::fs::write(
            &path,
            br#"{"theme": "dark", "preferences_page": "plugins", "ui_scale": 40, "sequence_font": "big",
                "restore_window": "yes", "check_updates": 0, "window": {"width": "wide"},
                "translation": "nine", "genetic_code": 7, "primer_binding": {"short": {"min_tm": "hot"}}}"#,
        )
        .unwrap();
        let (loaded, msg) = load(&path);
        assert_eq!((loaded.theme, msg), (ThemeChoice::Dark, None));
        assert_eq!(
            (
                loaded.preferences_page,
                loaded.ui_scale,
                loaded.restore_window,
                loaded.check_updates,
                loaded.window
            ),
            (PrefsPage::General, 2.0, true, true, None)
        );
        assert_eq!(loaded.sequence_font, crate::state::SEQUENCE_FONT_DEFAULT);
        assert_eq!(
            (
                loaded.translation,
                loaded.genetic_code,
                loaded.primer_binding
            ),
            (Translation::default(), 1, BindingPrefs::default())
        );
        std::fs::remove_dir_all(&dir).unwrap();
    }

    #[test]
    fn window_sizes_we_reopen_at() {
        let size = |width, height| WindowSize {
            width,
            height,
            maximized: false,
        };
        assert_eq!(size(1400.0, 900.0).usable(), Some([1400.0, 900.0]));
        assert_eq!(size(10.0, 1e9).usable(), Some([640.0, 16_384.0]));
        for bad in [size(f32::NAN, 900.0), size(1400.0, -1.0), size(0.0, 0.0)] {
            assert_eq!(bad.usable(), None, "{bad:?}");
        }
    }
}
