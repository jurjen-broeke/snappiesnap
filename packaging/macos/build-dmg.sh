#!/usr/bin/env bash
# Build SnappieSnap.app as a universal binary (Apple silicon + Intel) and wrap it in a .dmg.
# Runs on any Mac with Xcode's command-line tools and rustup: the release workflow's macos job,
# or by hand.
#
# Optional signing and notarisation (environment, set from repository secrets in CI):
#   MACOS_SIGN_IDENTITY   "Developer ID Application: <name> (<team id>)", in the keychain
#   MACOS_NOTARY_PROFILE  a notarytool keychain profile (xcrun notarytool store-credentials)
# Without them the app is signed ad hoc: it runs, but Gatekeeper warns on other Macs.
set -euo pipefail

cd "$(dirname "$0")/../.."
VERSION="$(sed -n 's/^version = "\(.*\)"/\1/p' Cargo.toml | head -n1)"
BUILD="${GITHUB_RUN_NUMBER:-0}"
BUNDLE_ID="${MACOS_BUNDLE_ID:-org.snappiesnap.SnappieSnap}"
OUT="dist"
APP="${OUT}/SnappieSnap.app"
DMG="${OUT}/SnappieSnap-${VERSION}-macos-universal.dmg"
export MACOSX_DEPLOYMENT_TARGET=11.0

rustup toolchain install
rustup target add aarch64-apple-darwin x86_64-apple-darwin
for target in aarch64-apple-darwin x86_64-apple-darwin; do
  cargo build --release --locked --target "$target" -p snappiesnap -p snappiesnap-cli --features snappiesnap/update-check
done

rm -rf "$OUT" && mkdir -p "$APP/Contents/MacOS" "$APP/Contents/Resources"
for bin in snappiesnap snappiesnap-cli; do
  lipo -create -output "$APP/Contents/MacOS/$bin" \
    "target/aarch64-apple-darwin/release/$bin" "target/x86_64-apple-darwin/release/$bin"
done
sed -e "s/__VERSION__/${VERSION}/" -e "s/__BUILD__/${BUILD}/" -e "s/__BUNDLE_ID__/${BUNDLE_ID}/" \
  packaging/macos/Info.plist > "$APP/Contents/Info.plist"
cp README.md ATTRIBUTION.md LICENSE-MIT LICENSE-APACHE "$APP/Contents/Resources/"
cp assets/brand/icons/macos/SnappieSnap.icns "$APP/Contents/Resources/"
lipo -info "$APP/Contents/MacOS/snappiesnap"

if [[ -n "${MACOS_SIGN_IDENTITY:-}" ]]; then
  codesign --force --options runtime --timestamp --sign "$MACOS_SIGN_IDENTITY" \
    "$APP/Contents/MacOS/snappiesnap-cli"
  codesign --force --options runtime --timestamp --sign "$MACOS_SIGN_IDENTITY" "$APP"
else
  echo "warning: MACOS_SIGN_IDENTITY not set; signing ad hoc (Gatekeeper will warn)" >&2
  codesign --force --sign - "$APP/Contents/MacOS/snappiesnap-cli"
  codesign --force --sign - "$APP"
fi
codesign --verify --deep --strict "$APP"

STAGE="$(mktemp -d)"
cp -R "$APP" "$STAGE/"
ln -s /Applications "$STAGE/Applications"
hdiutil create -volname "SnappieSnap ${VERSION}" -srcfolder "$STAGE" -fs HFS+ -format UDZO -ov "$DMG"
rm -rf "$STAGE"

if [[ -n "${MACOS_SIGN_IDENTITY:-}" ]]; then
  codesign --force --timestamp --sign "$MACOS_SIGN_IDENTITY" "$DMG"
  if [[ -n "${MACOS_NOTARY_PROFILE:-}" ]]; then
    xcrun notarytool submit "$DMG" --keychain-profile "$MACOS_NOTARY_PROFILE" --wait
    xcrun stapler staple "$DMG"
  fi
fi
echo "built $DMG"
