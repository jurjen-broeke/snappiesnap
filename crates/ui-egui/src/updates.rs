//! Checking for a newer SnappieSnap, only when the user presses "Check for updates" in
//! Preferences; never automatically. The check downloads a small JSON file announcing the
//! latest release, `{"version": "0.3.0", "url": "https://..."}`, and compares versions. It
//! sends nothing about the user (a generic `SnappieSnap-dev` User-Agent).

use serde::Deserialize;

/// Where the latest release is announced (an `https://` address of the JSON file above).
/// The release workflow uploads `latest.json` as an asset of every release, and GitHub serves
/// the newest one at this address. `None` disables the check, and the button then says so.
pub const RELEASE_URL: Option<&str> =
    Some("https://github.com/jurjen-broeke/snappiesnap/releases/latest/download/latest.json");

/// Largest announcement we read.
const MAX_BYTES: usize = 64 * 1024;

/// What a check found, in words for the Preferences window.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CheckResult {
    UpToDate,
    Newer {
        version: String,
        url: Option<String>,
    },
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
            CheckResult::Newer { version, url } => match url {
                Some(url) => format!("SnappieSnap {version} is available: {url}"),
                None => format!("SnappieSnap {version} is available."),
            },
            CheckResult::Unavailable(why) => why.clone(),
        }
    }
}

#[derive(Deserialize)]
struct Release {
    version: String,
    #[serde(default)]
    url: Option<String>,
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
        CheckResult::Newer {
            version: release.version.trim().chars().take(40).collect(),
            // Only web addresses are offered to the user.
            url: release
                .url
                .filter(|u| u.starts_with("https://") && u.len() <= 500),
        }
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
        assert_eq!(
            newer,
            CheckResult::Newer {
                version: "0.10.0".into(),
                url: Some("https://x/y".into())
            }
        );
        assert_eq!(
            compare("0.2.0", br#"{"version": "0.2.0"}"#),
            CheckResult::UpToDate
        );
        assert_eq!(
            compare("0.2.0", br#"{"version": "0.1.9"}"#),
            CheckResult::UpToDate
        );
        // Only https links are passed on.
        assert_eq!(
            compare("0.2.0", br#"{"version": "1.0.0", "url": "javascript:x"}"#),
            CheckResult::Newer {
                version: "1.0.0".into(),
                url: None
            }
        );
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

    /// The address shipped in the build is an `https://` one, so the check can never be
    /// downgraded to plain HTTP by a careless edit.
    #[test]
    fn the_built_in_address_is_https() {
        if let Some(url) = RELEASE_URL {
            assert!(url.starts_with("https://"), "{url}");
        }
    }
}
