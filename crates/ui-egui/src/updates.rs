//! Checking for a newer SnappieSnap. The installed app (built with the `update-check` feature
//! of the `snappiesnap` app, as the release packages are) checks once at startup, unless
//! Preferences turns that off; anyone can check from Preferences with "Check for updates".
//! The check downloads a small JSON file announcing the latest release,
//! `{"version": "0.3.0", "url": "https://..."}`, from the project's GitHub repository and
//! nowhere else, and compares versions. It sends nothing about the user (a generic
//! `SnappieSnap-dev` User-Agent). Nothing is installed: the user downloads the new version in
//! their browser and installs it themselves.

use serde::Deserialize;

/// The project's GitHub repository: the only place the check reads from or links to.
pub const REPOSITORY: &str = "https://github.com/jurjen-broeke/snappiesnap";

/// Where the latest release is announced (an `https://` address of the JSON file above).
/// The release workflow uploads `latest.json` as an asset of every release, and GitHub serves
/// the newest one at this address. `None` disables the check, and the button then says so.
pub const RELEASE_URL: Option<&str> =
    Some("https://github.com/jurjen-broeke/snappiesnap/releases/latest/download/latest.json");

/// A newer release: its version and the pages to offer for it, all in the repository.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct NewVersion {
    /// `major.minor.patch`, rebuilt from the parsed numbers (never the announcement's text).
    pub version: String,
    /// The release's page on GitHub (its notes and every package).
    pub page: String,
    /// The package for this computer, or the release page where the user has a choice
    /// (Linux: .deb or .rpm).
    pub download: String,
}

impl NewVersion {
    fn new((major, minor, patch): (u64, u64, u64)) -> Self {
        let version = format!("{major}.{minor}.{patch}");
        let page = format!("{REPOSITORY}/releases/tag/v{version}");
        let file = if cfg!(target_os = "macos") {
            Some(format!("SnappieSnap-{version}-macos-universal.dmg"))
        } else if cfg!(target_os = "windows") {
            Some(format!("SnappieSnap-{version}-windows-x64-setup.exe"))
        } else {
            None
        };
        let download = match file {
            Some(file) => format!("{REPOSITORY}/releases/download/v{version}/{file}"),
            None => page.clone(),
        };
        NewVersion {
            version,
            page,
            download,
        }
    }
}

/// Largest announcement we read.
const MAX_BYTES: usize = 64 * 1024;

/// What a check found, in words for the Preferences window.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CheckResult {
    UpToDate,
    Newer(NewVersion),
    /// Nothing to check against, or the check failed: why.
    Unavailable(String),
}

impl CheckResult {
    pub fn message(&self) -> String {
        match self {
            CheckResult::UpToDate => format!(
                "SnappieSnap {} is the latest version.",
                env!("CARGO_PKG_VERSION")
            ),
            CheckResult::Newer(new) => format!(
                "SnappieSnap {} is available (this is {}). Download it with the button in the toolbar, or from {}",
                new.version,
                env!("CARGO_PKG_VERSION"),
                new.page
            ),
            CheckResult::Unavailable(why) => why.clone(),
        }
    }
}

/// The announcement. Its `url` is not used: every link is built from the version, so a
/// changed announcement can never send the user outside the repository.
#[derive(Deserialize)]
struct Release {
    version: String,
}

/// A `major.minor.patch` version (an optional leading `v`; anything after `-` or `+` is
/// ignored). `None` for anything else.
pub fn parse_version(text: &str) -> Option<(u64, u64, u64)> {
    let text = text.trim();
    let text = text.strip_prefix('v').unwrap_or(text);
    let core = text.split(['-', '+']).next()?;
    let mut parts = core.split('.').map(|p| p.parse::<u64>().ok());
    let version = (parts.next()??, parts.next()??, parts.next()??);
    parts.next().is_none().then_some(version)
}

/// Compare the running version with an announcement's bytes.
pub fn compare(current: &str, announcement: &[u8]) -> CheckResult {
    let unavailable = |why: &str| CheckResult::Unavailable(format!("Update check failed: {why}."));
    if announcement.len() > MAX_BYTES {
        return unavailable("the release announcement is too large");
    }
    let Ok(release) = serde_json::from_slice::<Release>(announcement) else {
        return unavailable("the release announcement could not be read");
    };
    let (Some(running), Some(latest)) = (parse_version(current), parse_version(&release.version))
    else {
        return unavailable("the release announcement has no usable version");
    };
    if latest > running {
        CheckResult::Newer(NewVersion::new(latest))
    } else {
        CheckResult::UpToDate
    }
}

/// Start a check against [`RELEASE_URL`]. Without a release address the answer comes at once;
/// otherwise the file is downloaded on a background thread.
pub fn start(
    platform: std::sync::Arc<dyn snappiesnap_engine::Platform>,
) -> std::sync::mpsc::Receiver<CheckResult> {
    start_from(platform, RELEASE_URL)
}

/// The same check against a given address, so both answers can be tested without a network.
fn start_from(
    platform: std::sync::Arc<dyn snappiesnap_engine::Platform>,
    release_url: Option<&'static str>,
) -> std::sync::mpsc::Receiver<CheckResult> {
    let (tx, rx) = std::sync::mpsc::channel();
    let Some(url) = release_url else {
        let _ = tx.send(CheckResult::Unavailable(
            "Releases aren't published yet, so there is nothing to check against.".into(),
        ));
        return rx;
    };
    #[cfg(not(target_arch = "wasm32"))]
    std::thread::spawn(move || {
        let result = match platform.fetch(url) {
            Ok(bytes) => compare(env!("CARGO_PKG_VERSION"), &bytes),
            Err(e) => CheckResult::Unavailable(format!("Update check failed: {e}.")),
        };
        let _ = tx.send(result);
    });
    #[cfg(target_arch = "wasm32")]
    {
        let _ = (platform, url);
        let _ = tx.send(CheckResult::Unavailable(
            "The web version is always the latest one.".into(),
        ));
    }
    rx
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn versions() {
        assert_eq!(parse_version("0.2.0"), Some((0, 2, 0)));
        assert_eq!(parse_version(" v1.10.3-beta+5 "), Some((1, 10, 3)));
        for bad in [
            "",
            "1",
            "1.2",
            "1.2.3.4",
            "a.b.c",
            "1..3",
            "-1.0.0",
            "99999999999999999999.0.0",
        ] {
            assert_eq!(parse_version(bad), None, "{bad}");
        }
    }

    #[test]
    fn announcements() {
        let newer = compare("0.2.0", br#"{"version": "0.10.0", "url": "https://x/y"}"#);
        let CheckResult::Newer(new) = newer else {
            panic!("{newer:?}");
        };
        assert_eq!(new.version, "0.10.0");
        assert_eq!(
            new.page,
            "https://github.com/jurjen-broeke/snappiesnap/releases/tag/v0.10.0"
        );
        assert_eq!(
            compare("0.2.0", br#"{"version": "0.2.0"}"#),
            CheckResult::UpToDate
        );
        assert_eq!(
            compare("0.2.0", br#"{"version": "0.1.9"}"#),
            CheckResult::UpToDate
        );
        // The announcement's link and the text around the version are never passed on:
        // every address is built from the parsed version, inside the repository.
        let CheckResult::Newer(new) = compare(
            "0.2.0",
            br#"{"version": " v1.0.0-rc1 ", "url": "https://elsewhere.example/x"}"#,
        ) else {
            panic!("newer");
        };
        assert_eq!(new.version, "1.0.0");
        for link in [&new.page, &new.download] {
            assert!(
                link.starts_with(&format!("{REPOSITORY}/releases/")),
                "{link}"
            );
            assert!(link.contains("/v1.0.0"), "{link}");
        }
        for bad in [
            &b"not json"[..],
            br#"{"version": 3}"#,
            br#"{"version": "soon"}"#,
            &vec![b' '; MAX_BYTES + 1],
        ] {
            assert!(
                matches!(compare("0.2.0", bad), CheckResult::Unavailable(_)),
                "{bad:?}"
            );
        }
    }

    fn platform() -> std::sync::Arc<dyn snappiesnap_engine::Platform> {
        std::sync::Arc::new(snappiesnap_engine::MemoryFs::default())
    }

    #[test]
    fn without_an_address_says_so() {
        let result = start_from(platform(), None).recv().unwrap();
        assert!(result.message().contains("aren't published"), "{result:?}");
    }

    /// A download that doesn't arrive is reported, not treated as "up to date": a check that
    /// failed must never tell the user they have the newest version.
    #[test]
    fn a_download_that_fails_says_why() {
        let result = start_from(platform(), Some("https://example.invalid/latest.json"))
            .recv()
            .unwrap();
        assert!(matches!(result, CheckResult::Unavailable(_)), "{result:?}");
        assert!(
            result.message().contains("Update check failed"),
            "{result:?}"
        );
    }

    /// The address shipped in the build is an `https://` one in the repository, so the check
    /// can never be downgraded to plain HTTP or pointed elsewhere by a careless edit.
    #[test]
    fn the_built_in_address_is_the_repository() {
        assert!(REPOSITORY.starts_with("https://github.com/"));
        if let Some(url) = RELEASE_URL {
            assert!(url.starts_with(&format!("{REPOSITORY}/releases/")), "{url}");
        }
    }

    /// The download is this platform's package, named as the release workflow names it.
    #[test]
    fn downloads_name_this_platforms_package() {
        let new = NewVersion::new((0, 45, 0));
        let base = format!("{REPOSITORY}/releases/download/v0.45.0/");
        if cfg!(target_os = "macos") {
            assert_eq!(
                new.download,
                format!("{base}SnappieSnap-0.45.0-macos-universal.dmg")
            );
        } else if cfg!(target_os = "windows") {
            assert_eq!(
                new.download,
                format!("{base}SnappieSnap-0.45.0-windows-x64-setup.exe")
            );
        } else {
            assert_eq!(new.download, new.page);
        }
    }
}
