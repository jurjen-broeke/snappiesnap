# Brand

`snappiesnap.svg` is the master drawing of the SnappieSnap logo. Every file in `icons/` is
generated from it by `cargo xtask icons`; edit the SVG, rerun the task and commit both.

| Platform | Files | Used by |
|---|---|---|
| Windows | `icons/windows/snappiesnap.ico` (16–256 px) | embedded in `snappiesnap.exe` (`apps/snappiesnap/build.rs`); the NSIS and MSI installers |
| macOS | `icons/macos/SnappieSnap.icns` (16–1024 px, on a rounded tile) | `SnappieSnap.app` (`packaging/macos/`) |
| macOS | `icons/macos/snappiesnap-512.png` | the Dock icon when the app runs outside a bundle |
| Linux | `icons/linux/hicolor/<size>/apps/snappiesnap.{png,svg}` | `/usr/share/icons/hicolor/` in the `.deb` and `.rpm`; `Icon=snappiesnap` in the desktop entry; the window icon (256 px) |
| Web | `icons/web/*` | the WASM build's page (below) |

For the web build, copy `icons/web/` next to `index.html` and add to its `<head>`:

```html
<link rel="icon" href="favicon.ico" sizes="48x48">
<link rel="icon" href="favicon.svg" type="image/svg+xml">
<link rel="apple-touch-icon" href="apple-touch-icon.png">
<link rel="manifest" href="manifest.webmanifest">
<meta name="theme-color" content="#18bc18">
```

With trunk, `<link data-trunk rel="copy-dir" href="../../assets/brand/icons/web">` copies the
folder (adjust the `href` paths above to `web/…`).
