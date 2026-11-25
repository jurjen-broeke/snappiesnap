//! The browser side of the web build (wasm only): the file picker, dropped files, downloads and
//! the page title. Files arrive asynchronously in an [`Inbox`], which the app empties at the
//! start of each frame (see [`crate::web_files`]). Everything is done through `web-sys`; there
//! is no handwritten JavaScript.

use std::sync::Arc;
use std::sync::atomic::{AtomicBool, Ordering};

use eframe::egui;
use wasm_bindgen::JsCast;
use wasm_bindgen::closure::Closure;

use crate::web_files::{Inbox, Incoming, size_allowed};

/// A readable message for a JavaScript error value.
fn describe(e: &wasm_bindgen::JsValue) -> String {
    e.as_string()
        .or_else(|| {
            e.dyn_ref::<js_sys::Error>()
                .map(|e| String::from(e.message()))
        })
        .unwrap_or_else(|| "unknown browser error".to_string())
}

fn document() -> Result<web_sys::Document, String> {
    web_sys::window()
        .and_then(|w| w.document())
        .ok_or_else(|| "this page has no document".to_string())
}

/// Ask the browser for files (`accept` is a comma-separated list of extensions). Chosen files
/// are read in the background and land in `inbox`.
pub fn pick_files(
    inbox: &Inbox,
    accept: &str,
    multiple: bool,
    ctx: &egui::Context,
) -> Result<(), String> {
    let input: web_sys::HtmlInputElement = document()?
        .create_element("input")
        .map_err(|e| describe(&e))?
        .dyn_into()
        .map_err(|_| "the browser didn't make a file field".to_string())?;
    input.set_type("file");
    input.set_accept(accept);
    input.set_multiple(multiple);
    let (inbox, ctx, field) = (inbox.clone(), ctx.clone(), input.clone());
    // Called once, when files were chosen (a cancelled picker never calls it).
    let on_change = Closure::once_into_js(move || {
        let Some(files) = field.files() else {
            return;
        };
        for i in 0..files.length() {
            if let Some(file) = files.get(i) {
                read_file(inbox.clone(), file, ctx.clone());
            }
        }
    });
    input.set_onchange(Some(on_change.unchecked_ref()));
    input.click();
    Ok(())
}

/// Read a picked file in the background.
fn read_file(inbox: Inbox, file: web_sys::File, ctx: egui::Context) {
    let name = file.name();
    if !size_allowed(file.size()) {
        inbox.push(Incoming {
            bytes: Err(format!("{name} is too large to open")),
            name,
        });
        ctx.request_repaint();
        return;
    }
    wasm_bindgen_futures::spawn_local(async move {
        let bytes = wasm_bindgen_futures::JsFuture::from(file.array_buffer())
            .await
            .map(|buffer| js_sys::Uint8Array::new(&buffer).to_vec())
            .map_err(|e| format!("the browser couldn't read {name}: {}", describe(&e)));
        inbox.push(Incoming { name, bytes });
        ctx.request_repaint();
    });
}

/// Read files dropped on the page in the background.
pub fn read_dropped(inbox: &Inbox, files: Vec<egui::DroppedFileHandle>, ctx: &egui::Context) {
    for file in files {
        let (inbox, ctx) = (inbox.clone(), ctx.clone());
        wasm_bindgen_futures::spawn_local(async move {
            let name = file.path().to_string_lossy().into_owned();
            let bytes = file
                .bytes_async()
                .await
                .and_then(|b| {
                    if size_allowed(b.len() as f64) {
                        Ok(b)
                    } else {
                        Err("too large to open".to_string())
                    }
                })
                .map_err(|e| format!("{name}: {e}"));
            inbox.push(Incoming { name, bytes });
            ctx.request_repaint();
        });
    }
}

/// Hand `bytes` to the browser as a download named `name`.
pub fn download(name: &str, bytes: &[u8]) -> Result<(), String> {
    let window = web_sys::window().ok_or("this page has no window")?;
    let document = document()?;
    let parts = js_sys::Array::of1(&js_sys::Uint8Array::from(bytes));
    let options = web_sys::BlobPropertyBag::new();
    options.set_type("application/octet-stream");
    let blob = web_sys::Blob::new_with_u8_array_sequence_and_options(&parts, &options)
        .map_err(|e| describe(&e))?;
    let url = web_sys::Url::create_object_url_with_blob(&blob).map_err(|e| describe(&e))?;
    let link: web_sys::HtmlAnchorElement = document
        .create_element("a")
        .map_err(|e| describe(&e))?
        .dyn_into()
        .map_err(|_| "the browser didn't make a link".to_string())?;
    link.set_href(&url);
    link.set_download(name);
    // Some browsers only follow links that are part of the page.
    let body = document.body().ok_or("this page has no body")?;
    body.append_child(&link).map_err(|e| describe(&e))?;
    link.click();
    link.remove();
    // Free the file's memory once the browser has surely taken it.
    let revoke = Closure::once_into_js(move || {
        let _ = web_sys::Url::revoke_object_url(&url);
    });
    window
        .set_timeout_with_callback_and_timeout_and_arguments_0(revoke.unchecked_ref(), 60_000)
        .map_err(|e| describe(&e))?;
    Ok(())
}

/// Have the browser ask before the page is closed or reloaded while `unsaved` is set: the web
/// build keeps documents in memory only, so leaving would lose them.
pub fn warn_before_leaving(unsaved: Arc<AtomicBool>) -> Result<(), String> {
    let window = web_sys::window().ok_or("this page has no window")?;
    let on_leave = Closure::<dyn Fn(web_sys::BeforeUnloadEvent)>::new(
        move |event: web_sys::BeforeUnloadEvent| {
            if unsaved.load(Ordering::Relaxed) {
                event.prevent_default();
                // Older browsers ask only when a return value is set.
                event.set_return_value("unsaved changes");
            }
        },
    );
    window
        .add_event_listener_with_callback("beforeunload", on_leave.as_ref().unchecked_ref())
        .map_err(|e| describe(&e))?;
    // The listener lives as long as the page.
    on_leave.forget();
    Ok(())
}

/// Set the browser tab's title.
pub fn set_title(title: &str) {
    if let Ok(d) = document() {
        d.set_title(title);
    }
}
