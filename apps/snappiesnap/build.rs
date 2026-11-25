//! Embeds the application icon in `snappiesnap.exe`. Does nothing for other targets.

use embed_resource::CompilationResult;

fn main() {
    println!("cargo:rerun-if-changed=snappiesnap.rc");
    println!("cargo:rerun-if-changed=../../assets/brand/icons/windows/snappiesnap.ico");
    match embed_resource::compile("snappiesnap.rc", embed_resource::NONE) {
        CompilationResult::NotWindows | CompilationResult::Ok => {}
        // No resource compiler (rc.exe or windres): still build, without an icon.
        CompilationResult::NotAttempted(why) => {
            println!("cargo:warning=snappiesnap.exe gets no icon: {why}");
        }
        CompilationResult::Failed(why) => {
            eprintln!("error: cannot embed the application icon: {why}");
            std::process::exit(1);
        }
    }
}
