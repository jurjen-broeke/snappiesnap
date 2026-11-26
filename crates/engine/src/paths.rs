//! Where SnappieSnap keeps its own files: always a folder of the user's own, never next to the
//! application or in the working folder (which may be read-only, shared, or watched by
//! security software).
//!
//! | | data (libraries, downloads, caches) | settings |
//! |---|---|---|
//! | macOS | `~/Library/Application Support/SnappieSnap` | the same folder |
//! | Windows | `%APPDATA%\SnappieSnap` (else `%USERPROFILE%\AppData\Roaming\SnappieSnap`) | the same folder |
//! | Linux and others | `$XDG_DATA_HOME/snappiesnap` (else `~/.local/share/snappiesnap`) | `$XDG_CONFIG_HOME/snappiesnap` (else `~/.config/snappiesnap`) |
//!
//! `SNAPPIESNAP_DATA` and `SNAPPIESNAP_SETTINGS` override these. Like relative XDG variables
//! (which the XDG Base Directory specification says to ignore), relative overrides are
//! ignored, so nothing lands in whatever folder the app was started from. With no home folder
//! at all there is no location, and callers report that instead of guessing.

use std::ffi::OsString;
use std::path::PathBuf;

/// The operating system, for choosing locations.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Os {
    MacOs,
    Windows,
    Other,
}

impl Os {
    pub fn current() -> Os {
        if cfg!(target_os = "macos") {
            Os::MacOs
        } else if cfg!(target_os = "windows") {
            Os::Windows
        } else {
            Os::Other
        }
    }
}

/// An environment variable as an absolute path; None when unset, empty or relative.
fn absolute(env: &dyn Fn(&str) -> Option<OsString>, name: &str) -> Option<PathBuf> {
    let value = env(name)?;
    let path = PathBuf::from(value);
    (!path.as_os_str().is_empty() && is_absolute(&path)).then_some(path)
}

/// Absolute on any OS's terms: `/…`, or a Windows drive or UNC path.
///
/// Decided from the path's own text, never from the host's rules, so a location computed for
/// one [`Os`] is the same on every platform. `Path::is_absolute` cannot do this: on Windows it
/// rejects `/Users/ann`, which would make a macOS location depend on the machine working it
/// out, and would silently ignore a rooted `SNAPPIESNAP_DATA` on Windows.
fn is_absolute(path: &std::path::Path) -> bool {
    let s = path.to_string_lossy();
    let b = s.as_bytes();
    // Rooted: POSIX `/…`, and on Windows a UNC path, a verbatim prefix or a drive-rooted `\…`.
    if matches!(b.first(), Some(b'/' | b'\\')) {
        return true;
    }
    // A Windows drive: `C:\…` or `C:/…`. `C:x` is relative to that drive's own folder.
    if b.len() >= 3
        && b.first().is_some_and(u8::is_ascii_alphabetic)
        && b.get(1) == Some(&b':')
        && matches!(b.get(2), Some(b'\\' | b'/'))
    {
        return true;
    }
    // Anything else the host considers absolute.
    path.is_absolute()
}

/// The user's own application folder (shared by data and, on macOS and Windows, settings).
fn app_folder(os: Os, env: &dyn Fn(&str) -> Option<OsString>) -> Option<PathBuf> {
    match os {
        Os::MacOs => Some(absolute(env, "HOME")?.join("Library/Application Support/SnappieSnap")),
        Os::Windows => absolute(env, "APPDATA")
            .or_else(|| {
                Some(
                    absolute(env, "USERPROFILE")?
                        .join("AppData")
                        .join("Roaming"),
                )
            })
            .map(|d| d.join("SnappieSnap")),
        Os::Other => absolute(env, "XDG_DATA_HOME")
            .or_else(|| Some(absolute(env, "HOME")?.join(".local/share")))
            .map(|d| d.join("snappiesnap")),
    }
}

/// The folder for SnappieSnap's data: libraries, downloads and caches.
pub fn data_dir_with(os: Os, env: &dyn Fn(&str) -> Option<OsString>) -> Option<PathBuf> {
    absolute(env, "SNAPPIESNAP_DATA").or_else(|| app_folder(os, env))
}

/// The settings file.
pub fn settings_file_with(os: Os, env: &dyn Fn(&str) -> Option<OsString>) -> Option<PathBuf> {
    if let Some(p) = absolute(env, "SNAPPIESNAP_SETTINGS") {
        return Some(p);
    }
    let dir = match os {
        Os::MacOs | Os::Windows => app_folder(os, env)?,
        Os::Other => absolute(env, "XDG_CONFIG_HOME")
            .or_else(|| Some(absolute(env, "HOME")?.join(".config")))?
            .join("snappiesnap"),
    };
    Some(dir.join("settings.json"))
}

fn process_env(name: &str) -> Option<OsString> {
    std::env::var_os(name)
}

/// [`data_dir_with`] for this computer.
pub fn data_dir() -> Option<PathBuf> {
    data_dir_with(Os::current(), &process_env)
}

/// [`settings_file_with`] for this computer.
pub fn settings_file() -> Option<PathBuf> {
    settings_file_with(Os::current(), &process_env)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn env(pairs: &'static [(&'static str, &'static str)]) -> impl Fn(&str) -> Option<OsString> {
        move |name| {
            pairs
                .iter()
                .find(|(k, _)| *k == name)
                .map(|(_, v)| OsString::from(*v))
        }
    }

    /// Absoluteness is read from the path's text, not the host's rules: the same answer on
    /// Linux, macOS and Windows. Without this, a location worked out for one OS would depend on
    /// the machine working it out, and a rooted override would be dropped on Windows.
    #[test]
    fn absolute_paths_are_judged_the_same_on_every_host() {
        for yes in [
            "/Users/ann",
            "/home/ann",
            "/",
            r"C:\Users\ann",
            "C:/Users/ann",
            r"c:\x",
            r"\\server\share\x",
            r"\\?\C:\x",
            r"\Users\ann",
        ] {
            assert!(is_absolute(std::path::Path::new(yes)), "{yes}");
        }
        for no in ["", "data", "rel/x", r"rel\x", "C:x", "./x", "../x", "~/x"] {
            assert!(!is_absolute(std::path::Path::new(no)), "{no}");
        }
    }

    #[test]
    fn locations_are_in_the_users_home() {
        let mac = env(&[("HOME", "/Users/ann")]);
        assert_eq!(
            data_dir_with(Os::MacOs, &mac),
            Some(PathBuf::from(
                "/Users/ann/Library/Application Support/SnappieSnap"
            ))
        );
        assert_eq!(
            settings_file_with(Os::MacOs, &mac),
            Some(PathBuf::from(
                "/Users/ann/Library/Application Support/SnappieSnap/settings.json"
            ))
        );
        let linux = env(&[("HOME", "/home/ann")]);
        assert_eq!(
            data_dir_with(Os::Other, &linux),
            Some(PathBuf::from("/home/ann/.local/share/snappiesnap"))
        );
        assert_eq!(
            settings_file_with(Os::Other, &linux),
            Some(PathBuf::from("/home/ann/.config/snappiesnap/settings.json"))
        );
        let xdg = env(&[("HOME", "/home/ann"), ("XDG_DATA_HOME", "/data/ann")]);
        assert_eq!(
            data_dir_with(Os::Other, &xdg),
            Some(PathBuf::from("/data/ann/snappiesnap"))
        );
    }

    #[test]
    fn windows_uses_appdata_or_the_profile() {
        let appdata = env(&[("APPDATA", r"C:\Users\ann\AppData\Roaming")]);
        let dir = data_dir_with(Os::Windows, &appdata).unwrap();
        assert!(
            dir.to_string_lossy()
                .starts_with(r"C:\Users\ann\AppData\Roaming")
        );
        assert!(dir.ends_with("SnappieSnap"));
        let profile = env(&[("USERPROFILE", r"C:\Users\ann")]);
        let dir = data_dir_with(Os::Windows, &profile).unwrap();
        assert!(dir.to_string_lossy().starts_with(r"C:\Users\ann"));
        assert!(dir.ends_with("SnappieSnap"));
        assert!(dir.to_string_lossy().contains("Roaming"));
    }

    #[test]
    fn relative_or_missing_locations_are_never_used() {
        // Relative overrides and XDG variables fall back to the home folder…
        let relative = env(&[
            ("HOME", "/home/ann"),
            ("SNAPPIESNAP_DATA", "data"),
            ("SNAPPIESNAP_SETTINGS", "s.json"),
            ("XDG_DATA_HOME", "rel"),
        ]);
        assert_eq!(
            data_dir_with(Os::Other, &relative),
            Some(PathBuf::from("/home/ann/.local/share/snappiesnap"))
        );
        assert_eq!(
            settings_file_with(Os::Other, &relative),
            Some(PathBuf::from("/home/ann/.config/snappiesnap/settings.json"))
        );
        // …and with no home at all there is no location (not the working folder).
        let none = env(&[("HOME", "relative/home")]);
        assert_eq!(data_dir_with(Os::Other, &none), None);
        assert_eq!(data_dir_with(Os::MacOs, &none), None);
        assert_eq!(data_dir_with(Os::Windows, &none), None);
        // Absolute overrides win.
        let set = env(&[("HOME", "/home/ann"), ("SNAPPIESNAP_DATA", "/srv/snap")]);
        assert_eq!(
            data_dir_with(Os::Other, &set),
            Some(PathBuf::from("/srv/snap"))
        );
    }
}
