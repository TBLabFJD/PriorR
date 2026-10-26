#!/bin/bash
# Builds the macOS installer of PriorR (PriorR-<version>.pkg) on a Mac (or in GitHub Actions):
#   bash installers/macos/build_pkg.sh
# The installer is a wizard (welcome, license, install) that puts PriorR.app in the Applications
# folder of the user (no administrator needed). PriorR.app opens a Terminal window running PriorR
# (closing it stops PriorR). The data of the lab are not inside the app: the first start asks
# where the shared folder is. If R is missing, the installer says so and opens the download page.
set -euo pipefail
VERSION=4.0
HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=$(cd "$HERE/../.." && pwd)
BUILD="$HERE/build"
rm -rf "$BUILD"
APP="$BUILD/root/PriorR.app"
RES="$APP/Contents/Resources"
mkdir -p "$APP/Contents/MacOS" "$RES/PriorR"

# PriorR itself (without the data of this computer)
rsync -a --exclude 'cache/' --exclude 'shiny_bookmarks/' "$ROOT/app" "$RES/PriorR/"
cp "$ROOT/launch.R" "$ROOT/PriorR.command" "$ROOT/priorr.cfg.example" "$ROOT/README.md" "$RES/PriorR/"
chmod +x "$RES/PriorR/PriorR.command"

# Icon (.icns) from the PNG
ICONSET="$BUILD/priorr.iconset"
mkdir -p "$ICONSET"
for size in 16 32 128 256 512; do
  sips -z $size $size "$ROOT/installers/priorr_icon.png" --out "$ICONSET/icon_${size}x${size}.png" >/dev/null
  sips -z $((size * 2)) $((size * 2)) "$ROOT/installers/priorr_icon.png" --out "$ICONSET/icon_${size}x${size}@2x.png" >/dev/null
done
iconutil -c icns "$ICONSET" -o "$RES/priorr.icns"

cat > "$APP/Contents/Info.plist" <<EOF
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE plist PUBLIC "-//Apple//DTD PLIST 1.0//EN" "http://www.apple.com/DTDs/PropertyList-1.0.dtd">
<plist version="1.0">
<dict>
  <key>CFBundleName</key><string>PriorR</string>
  <key>CFBundleDisplayName</key><string>PriorR</string>
  <key>CFBundleIdentifier</key><string>es.iis-fjd.priorr</string>
  <key>CFBundleVersion</key><string>$VERSION</string>
  <key>CFBundleShortVersionString</key><string>$VERSION</string>
  <key>CFBundlePackageType</key><string>APPL</string>
  <key>CFBundleExecutable</key><string>PriorR</string>
  <key>CFBundleIconFile</key><string>priorr</string>
  <key>LSMinimumSystemVersion</key><string>11.0</string>
</dict>
</plist>
EOF

# The app opens PriorR.command in Terminal (the window of PriorR, like on Windows)
cat > "$APP/Contents/MacOS/PriorR" <<'EOF'
#!/bin/bash
HERE=$(cd "$(dirname "$0")/../Resources/PriorR" && pwd)
open -a Terminal "$HERE/PriorR.command"
EOF
chmod +x "$APP/Contents/MacOS/PriorR"

# After installing: is R there?
mkdir -p "$BUILD/scripts"
cat > "$BUILD/scripts/postinstall" <<'EOF'
#!/bin/bash
for r in "$(command -v Rscript)" /usr/local/bin/Rscript /opt/homebrew/bin/Rscript /Library/Frameworks/R.framework/Resources/bin/Rscript; do
  [ -n "$r" ] && [ -x "$r" ] && exit 0
done
answer=$(osascript -e 'display dialog "PriorR needs R, which is not installed on this Mac.\n\nInstall R from CRAN and then open PriorR." buttons {"Later", "Open the download page"} default button 2 with title "PriorR" with icon caution' 2>/dev/null || true)
case "$answer" in *"Open the download page"*) open "https://cran.r-project.org/bin/macosx/" ;; esac
exit 0
EOF
chmod +x "$BUILD/scripts/postinstall"

# Component package, then the installer with its wizard pages, for the user's own Applications
pkgbuild --root "$BUILD/root" --identifier es.iis-fjd.priorr --version "$VERSION" \
  --install-location /Applications --scripts "$BUILD/scripts" "$BUILD/PriorR-component.pkg"
mkdir -p "$BUILD/resources"
cp "$ROOT/installers/LICENSE.txt" "$BUILD/resources/license.txt"
cat > "$BUILD/resources/welcome.html" <<EOF
<html><body style="font-family: -apple-system, Helvetica; font-size: 13px;">
<h2>PriorR v$VERSION</h2>
<p>Filtering and prioritization of disease-linked genetic variants.</p>
<p>PriorR will be installed in the Applications folder of your user (no administrator needed).
It needs <b>R</b> (<a href="https://cran.r-project.org/bin/macosx/">cran.r-project.org</a>); the first
start installs the R packages it uses (internet needed) and asks for the shared folder of the lab.</p>
<p style="color: #666;">Bioinformatics Unit - Genetics and Genomics Department, IIS-Fundación Jiménez Díaz</p>
</body></html>
EOF
cat > "$BUILD/distribution.xml" <<EOF
<?xml version="1.0" encoding="utf-8"?>
<installer-gui-script minSpecVersion="2">
  <title>PriorR v$VERSION</title>
  <welcome file="welcome.html" mime-type="text/html"/>
  <license file="license.txt" mime-type="text/plain"/>
  <domains enable_anywhere="false" enable_currentUserHome="true" enable_localSystem="false"/>
  <options customize="never" require-scripts="false" hostArchitectures="x86_64,arm64"/>
  <choices-outline><line choice="priorr"/></choices-outline>
  <choice id="priorr" title="PriorR"><pkg-ref id="es.iis-fjd.priorr"/></choice>
  <pkg-ref id="es.iis-fjd.priorr" version="$VERSION" onConclusion="none">PriorR-component.pkg</pkg-ref>
</installer-gui-script>
EOF
mkdir -p "$HERE/Output"
productbuild --distribution "$BUILD/distribution.xml" --resources "$BUILD/resources" --package-path "$BUILD" \
  "$HERE/Output/PriorR-$VERSION.pkg"
echo "Installer: $HERE/Output/PriorR-$VERSION.pkg"
