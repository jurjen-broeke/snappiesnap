//! File access for commands, behind a trait (AGENTS.md › Platform services).
//!
//! The engine never touches `std::fs` directly: the desktop app and CLI use [`NativeFs`], the
//! web build and tests use [`MemoryFs`]. Writes are atomic on native file systems (write a
//! temporary file, then rename), so a crash or full disk can't leave a half-written document.

use std::collections::BTreeMap;
use std::io;
use std::sync::{Mutex, PoisonError};

/// Largest file a command will read (1 GiB).
pub const MAX_FILE: u64 = 1 << 30;

pub trait Platform: Send + Sync {
    fn read(&self, path: &str) -> io::Result<Vec<u8>>;
    fn write(&self, path: &str, bytes: &[u8]) -> io::Result<()>;

    /// Today's date (year, month, day) for file metadata, if the platform has a clock.
    fn today(&self) -> Option<(i32, u32, u32)> {
        None
    }

    /// Whether a file or folder exists at `path`.
    fn exists(&self, path: &str) -> bool {
        self.read(path).is_ok()
    }

    /// Rename or move a file or folder. Never replaces an existing entry.
    fn rename(&self, from: &str, to: &str) -> io::Result<()> {
        let _ = (from, to);
        Err(io::Error::new(
            io::ErrorKind::Unsupported,
            "this platform can't rename files",
        ))
    }

    /// Create a folder (its parent must exist).
    fn create_dir(&self, path: &str) -> io::Result<()> {
        let _ = path;
        Err(io::Error::new(
            io::ErrorKind::Unsupported,
            "this platform can't create folders",
        ))
    }

    /// Download a public data file. Only plain `http://` is supported on native platforms.
    fn fetch(&self, url: &str) -> io::Result<Vec<u8>> {
        let _ = url;
        Err(io::Error::new(
            io::ErrorKind::Unsupported,
            "this platform can't download files",
        ))
    }

    /// A folder for the app's own data, such as downloaded databases (created if needed).
    fn data_dir(&self) -> Option<String> {
        None
    }

    /// The entries of a folder, folders first, then by name. Hidden entries (starting with `.`)
    /// are left out, and at most [`MAX_ENTRIES`] are returned.
    fn list_dir(&self, path: &str) -> io::Result<Vec<DirEntry>> {
        let _ = path;
        Err(io::Error::new(
            io::ErrorKind::Unsupported,
            "this platform can't list folders",
        ))
    }
}

/// Most folder entries returned by [`Platform::list_dir`].
pub const MAX_ENTRIES: usize = 10_000;

/// One entry of a folder listing.
#[derive(Debug, Clone, PartialEq, Eq, serde::Serialize)]
pub struct DirEntry {
    pub name: String,
    /// The full path (folder path joined with the name).
    pub path: String,
    pub is_dir: bool,
}

fn sort_entries(entries: &mut Vec<DirEntry>) {
    entries.sort_by(|a, b| {
        b.is_dir
            .cmp(&a.is_dir)
            .then_with(|| a.name.to_lowercase().cmp(&b.name.to_lowercase()))
    });
    entries.truncate(MAX_ENTRIES);
}

#[cfg_attr(target_arch = "wasm32", allow(dead_code))]
/// The civil date (year, month, day) of a day count since 1970-01-01 (Howard Hinnant's
/// `civil_from_days`).
pub fn civil_from_days(days: i64) -> (i32, u32, u32) {
    let z = days + 719_468;
    let era = z.div_euclid(146_097);
    let doe = z.rem_euclid(146_097);
    let yoe = (doe - doe / 1460 + doe / 36_524 - doe / 146_096) / 365;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let day = doy - (153 * mp + 2) / 5 + 1;
    let month = if mp < 10 { mp + 3 } else { mp - 9 };
    let year = yoe + era * 400 + i64::from(month <= 2);
    (
        i32::try_from(year).unwrap_or(1970),
        u32::try_from(month).unwrap_or(1),
        u32::try_from(day).unwrap_or(1),
    )
}

/// The real file system.
#[cfg(not(target_arch = "wasm32"))]
#[derive(Debug, Default, Clone, Copy)]
pub struct NativeFs;

#[cfg(not(target_arch = "wasm32"))]
impl Platform for NativeFs {
    fn read(&self, path: &str) -> io::Result<Vec<u8>> {
        use std::io::Read;
        let file = std::fs::File::open(path)?;
        if file.metadata()?.len() > MAX_FILE {
            return Err(io::Error::other(format!(
                "file is larger than {} MiB",
                MAX_FILE >> 20
            )));
        }
        let mut bytes = Vec::new();
        file.take(MAX_FILE).read_to_end(&mut bytes)?;
        Ok(bytes)
    }

    fn write(&self, path: &str, bytes: &[u8]) -> io::Result<()> {
        use std::io::Write;
        if path.trim().is_empty() {
            return Err(io::Error::new(io::ErrorKind::InvalidInput, "empty path"));
        }
        let target = std::path::Path::new(path);
        let file_name = target
            .file_name()
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "path has no file name"))?;
        let mut temp_name = std::ffi::OsString::from(".");
        temp_name.push(file_name);
        temp_name.push(".snappiesnap-tmp");
        let temp = target.with_file_name(temp_name);
        let result = (|| {
            let mut f = std::fs::File::create(&temp)?;
            f.write_all(bytes)?;
            f.sync_all()?;
            std::fs::rename(&temp, target)
        })();
        if result.is_err() {
            let _ = std::fs::remove_file(&temp);
        }
        result
    }

    fn exists(&self, path: &str) -> bool {
        std::path::Path::new(path).exists()
    }

    fn rename(&self, from: &str, to: &str) -> io::Result<()> {
        if std::path::Path::new(to).exists() {
            return Err(io::Error::new(
                io::ErrorKind::AlreadyExists,
                format!("{to} already exists"),
            ));
        }
        std::fs::rename(from, to)
    }

    fn create_dir(&self, path: &str) -> io::Result<()> {
        std::fs::create_dir(path)
    }

    fn list_dir(&self, path: &str) -> io::Result<Vec<DirEntry>> {
        let mut entries = Vec::new();
        for entry in std::fs::read_dir(path)? {
            let Ok(entry) = entry else { continue };
            let name = entry.file_name().to_string_lossy().into_owned();
            if name.starts_with('.') {
                continue;
            }
            // Follow symlinks so linked folders behave like folders.
            let is_dir = std::fs::metadata(entry.path()).is_ok_and(|m| m.is_dir());
            entries.push(DirEntry {
                path: entry.path().to_string_lossy().into_owned(),
                name,
                is_dir,
            });
            if entries.len() > MAX_ENTRIES * 2 {
                break;
            }
        }
        sort_entries(&mut entries);
        Ok(entries)
    }

    /// `http://` with the built-in client; `https://` with the operating system's `curl`
    /// (macOS, Linux and Windows 10+ ship it), which keeps TLS and the system's trusted
    /// certificates up to date. The URL is one argument, never a shell string.
    fn fetch(&self, url: &str) -> io::Result<Vec<u8>> {
        if url.trim_start().starts_with("https://") {
            return curl(url.trim());
        }
        crate::http::get(url)
    }

    /// The user's own data folder (see [`crate::paths`]); never the application's folder.
    fn data_dir(&self) -> Option<String> {
        let dir = crate::paths::data_dir()?;
        std::fs::create_dir_all(&dir).ok()?;
        Some(dir.to_string_lossy().into_owned())
    }

    fn today(&self) -> Option<(i32, u32, u32)> {
        let secs = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .ok()?
            .as_secs();
        Some(civil_from_days(i64::try_from(secs / 86_400).ok()?))
    }
}

/// An in-memory file system (tests, and the web build until it gets browser storage).
#[derive(Debug, Default)]
pub struct MemoryFs {
    files: Mutex<BTreeMap<String, Vec<u8>>>,
    /// What [`Platform::fetch`] returns per URL.
    urls: BTreeMap<String, Vec<u8>>,
}

impl MemoryFs {
    pub fn new() -> Self {
        MemoryFs::default()
    }

    pub fn with_file(self, path: &str, bytes: Vec<u8>) -> Self {
        self.files
            .lock()
            .unwrap_or_else(PoisonError::into_inner)
            .insert(path.to_string(), bytes);
        self
    }

    /// Serve `bytes` for downloads of `url`.
    pub fn with_url(mut self, url: &str, bytes: Vec<u8>) -> Self {
        self.urls.insert(url.to_string(), bytes);
        self
    }

    pub fn get(&self, path: &str) -> Option<Vec<u8>> {
        self.files
            .lock()
            .unwrap_or_else(PoisonError::into_inner)
            .get(path)
            .cloned()
    }
}

impl Platform for MemoryFs {
    fn fetch(&self, url: &str) -> io::Result<Vec<u8>> {
        self.urls
            .get(url)
            .cloned()
            .ok_or_else(|| io::Error::new(io::ErrorKind::NotFound, format!("{url}: not found")))
    }

    fn data_dir(&self) -> Option<String> {
        Some("app-data".to_string())
    }

    fn read(&self, path: &str) -> io::Result<Vec<u8>> {
        self.get(path)
            .ok_or_else(|| io::Error::new(io::ErrorKind::NotFound, format!("{path}: not found")))
    }

    fn write(&self, path: &str, bytes: &[u8]) -> io::Result<()> {
        if path.trim().is_empty() {
            return Err(io::Error::new(io::ErrorKind::InvalidInput, "empty path"));
        }
        self.files
            .lock()
            .unwrap_or_else(PoisonError::into_inner)
            .insert(path.to_string(), bytes.to_vec());
        Ok(())
    }

    fn exists(&self, path: &str) -> bool {
        let prefix = format!("{}/", path.trim_end_matches('/'));
        let files = self.files.lock().unwrap_or_else(PoisonError::into_inner);
        files.contains_key(path) || files.keys().any(|k| k.starts_with(&prefix))
    }

    /// Renames a file, or every file under a folder.
    fn rename(&self, from: &str, to: &str) -> io::Result<()> {
        if self.exists(to) {
            return Err(io::Error::new(
                io::ErrorKind::AlreadyExists,
                format!("{to} already exists"),
            ));
        }
        let mut files = self.files.lock().unwrap_or_else(PoisonError::into_inner);
        let from_dir = format!("{}/", from.trim_end_matches('/'));
        let moved: Vec<String> = files
            .keys()
            .filter(|k| k.as_str() == from || k.starts_with(&from_dir))
            .cloned()
            .collect();
        if moved.is_empty() {
            return Err(io::Error::new(
                io::ErrorKind::NotFound,
                format!("{from}: not found"),
            ));
        }
        for key in moved {
            if let Some(bytes) = files.remove(&key) {
                let new_key = format!("{to}{}", key.get(from.len()..).unwrap_or(""));
                files.insert(new_key, bytes);
            }
        }
        Ok(())
    }

    /// Folders exist through their contents, so an empty folder holds a hidden marker.
    fn create_dir(&self, path: &str) -> io::Result<()> {
        if self.exists(path) {
            return Err(io::Error::new(
                io::ErrorKind::AlreadyExists,
                format!("{path} already exists"),
            ));
        }
        self.files
            .lock()
            .unwrap_or_else(PoisonError::into_inner)
            .insert(
                format!("{}/.folder", path.trim_end_matches('/')),
                Vec::new(),
            );
        Ok(())
    }

    /// Paths use `/`; folders are implied by the paths of the files in them. `.` is the root.
    fn list_dir(&self, path: &str) -> io::Result<Vec<DirEntry>> {
        let path = if path == "." { "" } else { path };
        let prefix = if path.is_empty() || path.ends_with('/') {
            path.to_string()
        } else {
            format!("{path}/")
        };
        let files = self.files.lock().unwrap_or_else(PoisonError::into_inner);
        let mut entries: Vec<DirEntry> = Vec::new();
        for key in files.keys() {
            let Some(rest) = key.strip_prefix(&prefix) else {
                continue;
            };
            let (name, is_dir) = match rest.split_once('/') {
                Some((dir, _)) => (dir, true),
                None => (rest, false),
            };
            if name.is_empty() || name.starts_with('.') || entries.iter().any(|e| e.name == name) {
                continue;
            }
            entries.push(DirEntry {
                name: name.to_string(),
                path: format!("{prefix}{name}"),
                is_dir,
            });
        }
        if entries.is_empty() && !prefix.is_empty() {
            return Err(io::Error::new(
                io::ErrorKind::NotFound,
                format!("{path}: no such folder"),
            ));
        }
        sort_entries(&mut entries);
        Ok(entries)
    }
}

/// The platform a new session uses: the real file system natively, memory on the web.
pub fn default_platform() -> std::sync::Arc<dyn Platform> {
    #[cfg(not(target_arch = "wasm32"))]
    {
        std::sync::Arc::new(NativeFs)
    }
    #[cfg(target_arch = "wasm32")]
    {
        std::sync::Arc::new(MemoryFs::new())
    }
}

/// Download an `https://` URL with the system's curl: HTTPS only (also on redirects), the
/// generic user agent, and time and size limits.
#[cfg(not(target_arch = "wasm32"))]
fn curl(url: &str) -> io::Result<Vec<u8>> {
    if url.len() > 8192 || url.chars().any(|c| c.is_whitespace() || c.is_control()) {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "not a valid URL",
        ));
    }
    let max = crate::http::MAX_DOWNLOAD.to_string();
    let out = std::process::Command::new("curl")
        .args([
            "--silent",
            "--show-error",
            "--fail",
            "--location",
            "--proto",
            "=https",
            "--proto-redir",
            "=https",
            "--max-time",
            "300",
            "--max-filesize",
            &max,
            "--user-agent",
            "SnappieSnap-dev",
            "--",
            url,
        ])
        .stdin(std::process::Stdio::null())
        .output()
        .map_err(|e| {
            io::Error::new(
                e.kind(),
                format!("cannot run curl, which HTTPS downloads use ({e})"),
            )
        })?;
    if !out.status.success() {
        let why = String::from_utf8_lossy(&out.stderr);
        return Err(io::Error::other(format!(
            "download failed: {}",
            why.trim().chars().take(300).collect::<String>()
        )));
    }
    if out.stdout.len() > crate::http::MAX_DOWNLOAD {
        return Err(io::Error::other("the download is too large"));
    }
    Ok(out.stdout)
}

#[cfg(all(test, not(target_arch = "wasm32")))]
mod tests {
    use super::*;

    #[test]
    fn native_write_is_atomic_and_readable() {
        let dir = std::env::temp_dir().join(format!("snappiesnap-test-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("x.dna");
        let p = path.to_string_lossy().into_owned();
        NativeFs.write(&p, b"one").unwrap();
        NativeFs.write(&p, b"two").unwrap();
        assert_eq!(NativeFs.read(&p).unwrap(), b"two");
        // No temporary file is left behind.
        let leftovers = std::fs::read_dir(&dir).unwrap().count();
        assert_eq!(leftovers, 1);
        assert!(NativeFs.write("", b"x").is_err());
        assert!(
            NativeFs
                .read(&dir.join("missing").to_string_lossy())
                .is_err()
        );
        std::fs::remove_dir_all(&dir).unwrap();
    }

    #[test]
    fn civil_dates() {
        assert_eq!(civil_from_days(0), (1970, 1, 1));
        assert_eq!(civil_from_days(20_383), (2025, 10, 22));
        assert_eq!(civil_from_days(-1), (1969, 12, 31));
        assert_eq!(civil_from_days(11_016), (2000, 2, 29));
        assert!(NativeFs.today().is_some_and(|(y, _, _)| y >= 2024));
    }

    #[test]
    fn listing_folders() {
        let fs = MemoryFs::new()
            .with_file("lab/b.dna", vec![])
            .with_file("lab/A.gb", vec![])
            .with_file("lab/sub/c.fa", vec![])
            .with_file("lab/.hidden", vec![]);
        let names: Vec<(String, bool)> = fs
            .list_dir("lab")
            .unwrap()
            .into_iter()
            .map(|e| (e.name, e.is_dir))
            .collect();
        assert_eq!(
            names,
            vec![
                ("sub".into(), true),
                ("A.gb".into(), false),
                ("b.dna".into(), false)
            ]
        );
        assert_eq!(fs.list_dir("lab/sub").unwrap()[0].path, "lab/sub/c.fa");
        assert!(fs.list_dir("nope").is_err());

        // Rename, create, never overwrite.
        fs.rename("lab/sub", "lab/moved").unwrap();
        assert!(fs.exists("lab/moved/c.fa") && !fs.exists("lab/sub"));
        assert!(fs.rename("lab/b.dna", "lab/A.gb").is_err());
        fs.create_dir("lab/empty").unwrap();
        assert!(
            fs.list_dir("lab")
                .unwrap()
                .iter()
                .any(|e| e.name == "empty" && e.is_dir)
        );
        assert!(fs.create_dir("lab/empty").is_err());

        let dir = std::env::temp_dir().join(format!("snappiesnap-ls-{}", std::process::id()));
        std::fs::create_dir_all(dir.join("inner")).unwrap();
        std::fs::write(dir.join("x.dna"), b"").unwrap();
        std::fs::write(dir.join(".secret"), b"").unwrap();
        let listed = NativeFs.list_dir(&dir.to_string_lossy()).unwrap();
        let names: Vec<&str> = listed.iter().map(|e| e.name.as_str()).collect();
        assert_eq!(names, vec!["inner", "x.dna"]);
        assert!(
            NativeFs
                .list_dir(&dir.join("missing").to_string_lossy())
                .is_err()
        );
        let x = dir.join("x.dna").to_string_lossy().into_owned();
        let y = dir.join("inner/y.dna").to_string_lossy().into_owned();
        NativeFs.rename(&x, &y).unwrap();
        assert!(NativeFs.exists(&y) && !NativeFs.exists(&x));
        std::fs::write(&x, b"").unwrap();
        assert!(NativeFs.rename(&x, &y).is_err(), "never overwrites");
        let nested = dir.join("made").to_string_lossy().into_owned();
        NativeFs.create_dir(&nested).unwrap();
        assert!(NativeFs.create_dir(&nested).is_err());
        std::fs::remove_dir_all(&dir).unwrap();
    }

    #[test]
    fn memory_fs() {
        let fs = MemoryFs::new().with_file("a", b"1".to_vec());
        assert_eq!(fs.read("a").unwrap(), b"1");
        assert!(fs.read("b").is_err());
        fs.write("b", b"2").unwrap();
        assert_eq!(fs.get("b").unwrap(), b"2");
    }
}
