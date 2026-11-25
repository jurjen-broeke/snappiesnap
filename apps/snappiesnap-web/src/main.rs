//! SnappieSnap in the browser: the desktop app's UI and engine compiled to WebAssembly.
//!
//! Build it with trunk (`trunk build --release` in this folder; see `docs/web.md`). The result
//! is a folder of static files that any web server can serve. Everything runs in the user's
//! browser: files are opened with the browser's file picker or by dropping them on the page,
//! kept in memory, and saved as downloads. Nothing is sent to the server.
//!
//! A native build of this crate only explains how to build it, so `cargo build --workspace`
//! keeps working.

#[cfg(not(target_arch = "wasm32"))]
fn main() -> std::process::ExitCode {
    eprintln!(
        "snappiesnap-web runs in a browser. Build it with `trunk build --release` in \
         apps/snappiesnap-web (see docs/web.md)."
    );
    std::process::ExitCode::from(2)
}

#[cfg(target_arch = "wasm32")]
fn main() {
    wasm_bindgen_futures::spawn_local(async {
        if let Err(e) = web::start().await {
            web::show_message(&format!("SnappieSnap couldn't start: {e}"));
        }
    });
}

#[cfg(target_arch = "wasm32")]
mod web {
    use snappiesnap_engine::Session;
    use snappiesnap_ui_egui::eframe;
    use snappiesnap_ui_egui::{App, LaunchOptions};
    use wasm_bindgen::JsCast;

    /// The page's canvas and message box (see `index.html`).
    const CANVAS: &str = "snappiesnap";
    const MESSAGE: &str = "message";

    fn document() -> Result<web_sys::Document, String> {
        web_sys::window()
            .and_then(|w| w.document())
            .ok_or_else(|| "this page has no document".to_string())
    }

    /// Show `text` in the page's message box (startup problems, an internal error).
    pub fn show_message(text: &str) {
        if let Some(el) = document().ok().and_then(|d| d.get_element_by_id(MESSAGE)) {
            el.set_text_content(Some(text));
            let _ = el.set_attribute("class", "shown");
        }
    }

    fn hide_message() {
        if let Some(el) = document().ok().and_then(|d| d.get_element_by_id(MESSAGE)) {
            let _ = el.set_attribute("class", "");
        }
    }

    pub async fn start() -> Result<(), String> {
        let canvas: web_sys::HtmlCanvasElement = document()?
            .get_element_by_id(CANVAS)
            .ok_or("the page has no canvas")?
            .dyn_into()
            .map_err(|_| "the page's canvas isn't a canvas".to_string())?;
        // The runner logs a panic to the console; also tell the user, who would otherwise see
        // a frozen page. (Unwinding isn't available in the browser, so the app can't recover.)
        let runner = eframe::WebRunner::new();
        let previous = std::panic::take_hook();
        std::panic::set_hook(Box::new(move |info| {
            previous(info);
            show_message(
                "SnappieSnap stopped because of an internal error (please report it). \
                 Reload the page to start again; changes that weren't downloaded are lost.",
            );
        }));
        let mut session = Session::new();
        match snappiesnap_engine::sample::load(&mut session, true) {
            // A demo, not the user's work: closing the page needn't ask about it.
            Ok(()) => {
                if let Some(doc) = session.document_mut() {
                    doc.dirty = false;
                }
            }
            Err(e) => session.log(
                snappiesnap_engine::Level::Warn,
                "web",
                format!("the sample plasmid failed to load: {e}"),
            ),
        }
        let options = eframe::WebOptions::default();
        runner
            .start(
                canvas,
                options,
                Box::new(move |cc| {
                    Ok(Box::new(App::new(
                        &cc.egui_ctx,
                        session,
                        &LaunchOptions::default(),
                    )))
                }),
            )
            .await
            .map_err(|e| {
                e.as_string()
                    .unwrap_or_else(|| "the browser couldn't start the graphics".to_string())
            })?;
        hide_message();
        Ok(())
    }
}
