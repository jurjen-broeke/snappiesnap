; SnappieSnap installer (.exe), built on Linux with NSIS (makensis).
; Installs for the current user only, into %LOCALAPPDATA%\Programs\SnappieSnap: no
; administrator rights, no UAC prompt. SnappieSnap's own files live in %APPDATA%\SnappieSnap
; and are kept on uninstall (the user's library and settings).
;
; makensis -DVERSION=0.0.1 -DSRC=target/x86_64-pc-windows-gnu/release -DOUT=dist/x.exe installer.nsi

Unicode true
!include "MUI2.nsh"

!ifndef VERSION
  !error "pass -DVERSION=x.y.z"
!endif
!ifndef SRC
  !error "pass -DSRC=<folder with snappiesnap.exe>"
!endif
!ifndef OUT
  !define OUT "SnappieSnap-${VERSION}-windows-x64-setup.exe"
!endif

Name "SnappieSnap"
OutFile "${OUT}"
RequestExecutionLevel user
InstallDir "$LOCALAPPDATA\Programs\SnappieSnap"
InstallDirRegKey HKCU "Software\SnappieSnap" "InstallDir"
SetCompressor /SOLID lzma

VIProductVersion "${VERSION}.0"
VIAddVersionKey "ProductName" "SnappieSnap"
VIAddVersionKey "FileDescription" "SnappieSnap installer"
VIAddVersionKey "FileVersion" "${VERSION}"
VIAddVersionKey "ProductVersion" "${VERSION}"
VIAddVersionKey "LegalCopyright" "SnappieSnap contributors"

!define UNINSTALL_KEY "Software\Microsoft\Windows\CurrentVersion\Uninstall\SnappieSnap"
!define MUI_ICON "assets\brand\icons\windows\snappiesnap.ico"
!define MUI_UNICON "assets\brand\icons\windows\snappiesnap.ico"

!insertmacro MUI_PAGE_DIRECTORY
!insertmacro MUI_PAGE_INSTFILES
!insertmacro MUI_PAGE_FINISH
!insertmacro MUI_UNPAGE_CONFIRM
!insertmacro MUI_UNPAGE_INSTFILES
!insertmacro MUI_LANGUAGE "English"

Section "SnappieSnap"
  SetShellVarContext current
  SetOutPath "$INSTDIR"
  File "${SRC}\snappiesnap.exe"
  File "${SRC}\snappiesnap-cli.exe"
  File "/oname=README.txt" "README.md"
  File "/oname=ATTRIBUTION.txt" "ATTRIBUTION.md"
  File "/oname=LICENSE-MIT.txt" "LICENSE-MIT"
  File "/oname=LICENSE-APACHE.txt" "LICENSE-APACHE"
  WriteUninstaller "$INSTDIR\uninstall.exe"
  CreateShortcut "$SMPROGRAMS\SnappieSnap.lnk" "$INSTDIR\snappiesnap.exe"
  WriteRegStr HKCU "Software\SnappieSnap" "InstallDir" "$INSTDIR"
  WriteRegStr HKCU "${UNINSTALL_KEY}" "DisplayName" "SnappieSnap"
  WriteRegStr HKCU "${UNINSTALL_KEY}" "DisplayVersion" "${VERSION}"
  WriteRegStr HKCU "${UNINSTALL_KEY}" "Publisher" "SnappieSnap contributors"
  WriteRegStr HKCU "${UNINSTALL_KEY}" "InstallLocation" "$INSTDIR"
  WriteRegStr HKCU "${UNINSTALL_KEY}" "DisplayIcon" "$INSTDIR\snappiesnap.exe"
  WriteRegStr HKCU "${UNINSTALL_KEY}" "UninstallString" '"$INSTDIR\uninstall.exe"'
  WriteRegDWORD HKCU "${UNINSTALL_KEY}" "NoModify" 1
  WriteRegDWORD HKCU "${UNINSTALL_KEY}" "NoRepair" 1
SectionEnd

Section "Uninstall"
  SetShellVarContext current
  Delete "$SMPROGRAMS\SnappieSnap.lnk"
  Delete "$INSTDIR\snappiesnap.exe"
  Delete "$INSTDIR\snappiesnap-cli.exe"
  Delete "$INSTDIR\README.txt"
  Delete "$INSTDIR\ATTRIBUTION.txt"
  Delete "$INSTDIR\LICENSE-MIT.txt"
  Delete "$INSTDIR\LICENSE-APACHE.txt"
  Delete "$INSTDIR\uninstall.exe"
  RMDir "$INSTDIR"
  DeleteRegKey HKCU "${UNINSTALL_KEY}"
  DeleteRegKey HKCU "Software\SnappieSnap"
SectionEnd
