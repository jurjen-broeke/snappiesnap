//! Opening web pages in the user's default browser: on desktop through the operating system's
//! own launcher (`open` on macOS, `xdg-open` on Linux, the URL protocol handler on Windows),
//! on the web in a new tab. Only `https://` addresses are opened.

use eframe::egui;

/// Open `url` in the default browser.
pub fn open(ctx: &egui::Context, url: &str) -> Result<(), String> {
    if !url.starts_with("https://") || url.chars().any(|c| c.is_whitespace() || c.is_control()) {
        return Err("only https:// addresses are opened".into());
    }
    #[cfg(target_arch = "wasm32")]
    {
        ctx.open_url(egui::OpenUrl::new_tab(url));
        Ok(())
    }
    #[cfg(not(target_arch = "wasm32"))]
    {
        let _ = ctx;
        let mut command = if cfg!(target_os = "macos") {
            std::process::Command::new("open")
        } else if cfg!(target_os = "windows") {
            let mut c = std::process::Command::new("rundll32");
            c.arg("url.dll,FileProtocolHandler");
            c
        } else {
            std::process::Command::new("xdg-open")
        };
        let mut child = command
            .arg(url)
            .stdin(std::process::Stdio::null())
            .stdout(std::process::Stdio::null())
            .stderr(std::process::Stdio::null())
            .spawn()
            .map_err(|e| format!("couldn't open the browser ({e})"))?;
        // Reap the launcher when it exits, without waiting for it here.
        std::thread::spawn(move || {
            let _ = child.wait();
        });
        Ok(())
    }
}
