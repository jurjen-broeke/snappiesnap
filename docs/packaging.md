# Packaging and releases

**Status:** 2026-11-25. Packages are built by GitHub Actions on GitHub-hosted runners, and a
version tag publishes a GitHub release.
- **Tested:** the macOS `.dmg` was built and run on a Mac.
- **Not yet tested on GitHub Actions:** the Linux, Windows, web and publish jobs. They run for
  the first time on the first version tag pushed to GitHub; use the workflow's manual run
  (**Actions → release → Run workflow**, "dry run" left on) to build every package without
  publishing anything.
- **Signing is left out for now:** without certificates the Windows and macOS packages are
  unsigned, so other computers warn about them (see *Signing* below).

Everything is driven by two workflows, and no server of our own is involved:

| Workflow | Runs on | Does |
|---|---|---|
| [`.github/workflows/ci.yml`](../.github/workflows/ci.yml) | every push to `main`, every pull request | `cargo xtask ci` on Linux; `cargo test --workspace` on macOS and Windows |
| [`.github/workflows/release.yml`](../.github/workflows/release.yml) | a `v<version>` tag, or a manual run | the gate, then every package, then the GitHub release |

**Targets:** x86-64 for Linux (Debian and Red Hat families) and Windows, universal for macOS.
ARM Linux and Windows on ARM aren't packaged: those users run the web version.

**Why `cross` for Linux:** the GitHub runner is a current Ubuntu, so binaries built on it
directly would need its glibc. The Linux job builds in cross's CentOS 7 image (glibc 2.17)
instead, and then refuses any binary needing glibc newer than 2.28 — so a changed image cannot
silently raise what the packages require. The images are pinned per target in
[`Cross.toml`](../Cross.toml).

**Why MinGW for Windows:** the `x86_64-pc-windows-gnu` target needs no Microsoft SDK licence, so
the Windows binaries and both installers are cross-built on Linux (`makensis` for the `.exe`,
`wixl` for the `.msi`). On Debian and Ubuntu `wixl` is its own package: the `msitools` package
ships `msiinfo` and `msibuild` but not `wixl`.

**Why macOS needs a Mac:** building for macOS needs Apple's SDK, whose licence allows use on
Apple hardware only. The job runs on GitHub's `macos-latest` runner.

## What the packages contain

| Package | Installs | Rights needed |
|---|---|---|
| `snappiesnap_<v>_amd64.deb`, `snappiesnap-<v>-1.x86_64.rpm` | `/usr/bin/snappiesnap`, `/usr/bin/snappiesnap-cli`, a desktop entry and its icons | the system's package manager (sudo) |
| `SnappieSnap-<v>-windows-x64-setup.exe` | `%LOCALAPPDATA%\Programs\SnappieSnap` and a Start menu entry for the current user | none: no administrator rights, no UAC prompt |
| `SnappieSnap-<v>-windows-x64.msi` | `C:\Program Files\SnappieSnap`, for all users; meant for deployment by IT | administrator |
| `SnappieSnap-<v>-macos-universal.dmg` | `SnappieSnap.app`, dragged to Applications (or `~/Applications`) | none |
| `SnappieSnap-<v>-web.tar.gz` | nothing: static files for a web server ([web.md](web.md)) | none for users; the server's own |

Every release also carries `latest.json`, which the app's **Check for updates** button reads
from the latest release (`crates/ui-egui/src/updates.rs`). Nothing is checked automatically and
nothing about the user is ever sent.

SnappieSnap's own files (settings, libraries, downloads) are always kept in the user's own
folder, never next to the program. See the README.

**Linux dependencies:** the display, Vulkan and portal libraries are loaded at run time and are
only recommended dependencies, so the packages install on a minimal system.

## Releasing a version

1. Bump the version and its `CHANGELOG.md` entry as usual (AGENTS.md › Versioning), commit, and
   check that `main` is green.
2. Tag the commit with its version and push the tag. Pushing the tag is the decision to release:

   ```sh
   git tag v0.39.0 && git push origin v0.39.0
   ```

   The tag must match `version` in `Cargo.toml`, or the workflow stops before building anything.
3. The workflow runs `cargo xtask ci`, builds the Linux, Windows, macOS and web packages in
   parallel, and publishes the release with this version's changelog entry as its notes and
   every package attached.
4. If a job fails partway (a network error, say), re-run it. Publishing again is safe: an
   existing release keeps its notes and its packages are replaced rather than duplicated.

To build the packages without releasing anything, use **Actions → release → Run workflow** with
"dry run" left on. Artifacts are attached to that workflow run instead.

`packaging/macos/build-dmg.sh` also works by hand on any Mac: it writes
`dist/SnappieSnap-<v>-macos-universal.dmg`.

## Signing: avoiding security warnings on other computers

Unsigned packages work, but other computers warn about them:
- macOS: "cannot be opened because the developer cannot be verified".
- Windows: a SmartScreen warning.

Signing needs certificates that cost money and belong to a person or organisation. The release
workflow signs automatically once these **repository secrets** exist
(**Settings → Secrets and variables → Actions**), and warns in the job log while they don't:

| Platform | What is needed | Secrets |
|---|---|---|
| macOS | Apple Developer Program (about USD 99 a year): a "Developer ID Application" certificate exported as a `.p12`, and an App Store Connect API key for notarisation | `MACOS_CERT_P12_BASE64`, `MACOS_CERT_PASSWORD`, `MACOS_SIGN_IDENTITY`, `MACOS_NOTARY_API_KEY_BASE64`, `MACOS_NOTARY_API_KEY_ID`, `MACOS_NOTARY_API_ISSUER_ID`; optionally `MACOS_BUNDLE_ID` (default `org.snappiesnap.SnappieSnap`) |
| Windows | an Authenticode code-signing certificate as a `.p12`. Since 2023, new certificates are issued on hardware tokens or cloud HSMs; a cloud service such as Azure Trusted Signing needs `jsign` in place of `osslsigncode`. | `WINDOWS_CERT_P12_BASE64`, `WINDOWS_CERT_PASSWORD` |
| Linux | optional: GPG-signed packages or a repository | not set up |

Store a certificate with `base64 -i cert.p12 | pbcopy` (macOS) or `base64 -w0 cert.p12` (Linux)
and paste the result as the secret's value. The workflow decodes it to a temporary file, and on
macOS imports it into a keychain created for that run only.

## Open decisions before the first signed release

- **Bundle identifier and publisher name.** A reverse-DNS name on a domain you control, plus the
  name shown as publisher. Both are hard to change after the first release, because Windows
  recognises upgrades by them and macOS ties the bundle id to the signing identity. The default
  is `org.snappiesnap.SnappieSnap`.
- **The MSI's `UpgradeCode`** in `packaging/windows/snappiesnap.wxs` must never change: it is how
  Windows recognises newer versions.

## Files

| File | Purpose |
|---|---|
| `.github/workflows/ci.yml` | the quality gate on every push and pull request |
| `.github/workflows/release.yml` | the packages and the release, from a version tag |
| `Cross.toml` | the `cross` images per target (CentOS 7 for Linux, MinGW-w64 for Windows) |
| `apps/snappiesnap/Cargo.toml` (`[package.metadata.deb]`, `[package.metadata.generate-rpm]`) | Linux package contents |
| `packaging/linux/snappiesnap.desktop` | Linux menu entry |
| `packaging/windows/installer.nsi`, `packaging/windows/snappiesnap.wxs` | Windows installers |
| `packaging/macos/Info.plist`, `packaging/macos/build-dmg.sh` | macOS app bundle and `.dmg` |
| `apps/snappiesnap-web/index.html`, `apps/snappiesnap-web/Trunk.toml` | the web version ([web.md](web.md)) |
