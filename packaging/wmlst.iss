; SPDX-License-Identifier: GPL-2.0-only
; Copyright (C) 2025-2026 IOWA-Tech - Giovanni Lorenzin
; Copyright (C) Torsten Seemann (upstream `mlst`, from which WMLST is ported)
;
; Inno Setup 6 installer for WMLST -- docs/ARCHITECTURE.md section 13.1.
;
; Build order:
;   python packaging\make_icon.py
;   pyinstaller --noconfirm --clean packaging\wmlst.spec      -> dist\WMLST\
;   iscc packaging\wmlst.iss                                  -> packaging\Output\
;
; ---------------------------------------------------------------------------
; THIS INSTALLER IS NOT CODE-SIGNED. READ THIS BEFORE YOU SHIP IT.
; ---------------------------------------------------------------------------
; An OV/EV Authenticode certificate costs several hundred euro a year and is not
; something an independent tool can assume. The honest consequence, which the
; README, docs\WINDOWS.md and every GitHub Release body must state in the same
; words:
;
;   Windows SmartScreen will show "Windows protected your PC" the first time
;   this installer is run. That is a REPUTATION signal, not a virus detection:
;   SmartScreen has simply never seen this file before. The way through it is
;   More info -> Run anyway.
;
; What we do instead of signing:
;   * every release is built by a public GitHub Actions run, and the build log
;     is linked from the release body;
;   * SHA256SUMS.txt is published next to the binaries so anyone can verify the
;     download matches what CI produced.
;
; What we must NOT do, and what the documentation must never suggest:
;   * claim or imply that an unsigned binary is "safe" or "verified";
;   * tell users to disable SmartScreen or Defender, or to edit the registry;
;   * ship a self-signed certificate. SmartScreen ignores untrusted roots, so it
;     changes nothing except making the build look like it is hiding something.
; ---------------------------------------------------------------------------

#define MyAppName        "WMLST"
#define MyAppVersion     "1.0.2"
#define MyAppPublisher   "IOWA-Tech"
#define MyAppAuthor      "Giovanni Lorenzin"
#define MyAppURL         "https://github.com/iowa69/WMLST"
#define MyAppExeName     "WMLST.exe"
#define MyAppCliName     "wmlst-cli.exe"
#define MyAppId          "{{6C4E9E2B-1D57-4B9A-9F3C-7A21D5B0E8C4}"
#define SourceDir        "..\dist\WMLST"
#define RepoRoot         ".."

[Setup]
AppId={#MyAppId}
AppName={#MyAppName}
AppVersion={#MyAppVersion}
AppVerName={#MyAppName} {#MyAppVersion}
AppPublisher={#MyAppPublisher}
AppCopyright=Copyright (C) 2025-2026 {#MyAppPublisher} - {#MyAppAuthor}. GPL-2.0-only.
AppPublisherURL={#MyAppURL}
AppSupportURL={#MyAppURL}/issues
AppUpdatesURL={#MyAppURL}/releases
AppComments=MLST typing for Windows. A port of mlst by Torsten Seemann. Allele data from PubMLST.
AppContact={#MyAppAuthor}
VersionInfoVersion={#MyAppVersion}
VersionInfoCompany={#MyAppPublisher}
VersionInfoDescription={#MyAppName} Setup
VersionInfoCopyright=Copyright (C) 2025-2026 {#MyAppPublisher} - {#MyAppAuthor}

; PrivilegesRequired=lowest means a per-user install under %LOCALAPPDATA%, with
; no UAC prompt at all. It also matches where the BLAST+ bootstrap puts its
; binaries, so a user who can install WMLST can always finish the setup.
PrivilegesRequired=lowest
PrivilegesRequiredOverridesAllowed=dialog
DefaultDirName={autopf}\{#MyAppPublisher}\{#MyAppName}
DefaultGroupName={#MyAppPublisher}
DisableProgramGroupPage=yes
DisableDirPage=no
AllowNoIcons=yes
UsePreviousAppDir=yes

; Windows 10 1809 is the floor: it is the first release with a built-in UTF-8
; console and long-path support, both of which section 8.5 depends on.
MinVersion=10.0.17763
ArchitecturesAllowed=x64compatible
ArchitecturesInstallIn64BitMode=x64compatible

; The database compresses about 20:1, so lzma2/max plus solid compression turns
; a 195 MB payload into a ~28 MB installer. SolidCompression is what makes the
; 1,108 small allele files cheap.
Compression=lzma2/max
SolidCompression=yes
LZMAUseSeparateProcess=yes
LZMANumBlockThreads=4
InternalCompressLevel=max

OutputDir=Output
OutputBaseFilename=WMLST-{#MyAppVersion}-win64-setup
SetupIconFile=wmlst.ico
UninstallDisplayIcon={app}\{#MyAppExeName}
UninstallDisplayName={#MyAppName} {#MyAppVersion}
WizardStyle=modern
WizardSizePercent=110
ShowLanguageDialog=no
CloseApplications=yes
RestartApplications=no

; GPLv2 section 1 requires the licence to travel with every copy. Showing it in
; the wizard AND installing it into {app} is the part installer authors forget.
LicenseFile={#RepoRoot}\LICENSE
InfoBeforeFile=SMARTSCREEN.txt

[Languages]
Name: "english"; MessagesFile: "compiler:Default.isl"

[Tasks]
; Checked by default: the desktop shortcut is how most users will start WMLST.
Name: "desktopicon"; Description: "{cm:CreateDesktopIcon}"; GroupDescription: "{cm:AdditionalIcons}"
Name: "addtopath"; Description: "Add the command-line tool (wmlst-cli.exe) to my PATH"; GroupDescription: "Command line"; Flags: unchecked

[Files]
; The whole PyInstaller one-folder tree, including _internal\db (162 schemes).
Source: "{#SourceDir}\*"; DestDir: "{app}"; Flags: ignoreversion recursesubdirs createallsubdirs

; Licence and attribution, installed as .txt so Notepad opens them on a
; double-click. Required by GPLv2 sections 1 and 3(a).
Source: "{#RepoRoot}\LICENSE";      DestDir: "{app}"; DestName: "LICENSE.txt";  Flags: ignoreversion
Source: "{#RepoRoot}\NOTICE";       DestDir: "{app}"; DestName: "NOTICE.txt";   Flags: ignoreversion
Source: "{#RepoRoot}\README.md";    DestDir: "{app}"; DestName: "README.md";    Flags: ignoreversion
Source: "{#RepoRoot}\CITATION.cff"; DestDir: "{app}"; DestName: "CITATION.cff"; Flags: ignoreversion
Source: "SMARTSCREEN.txt";          DestDir: "{app}"; DestName: "SMARTSCREEN.txt"; Flags: ignoreversion

[Icons]
Name: "{group}\{#MyAppName}"; Filename: "{app}\{#MyAppExeName}"; Comment: "MLST typing for Windows - {#MyAppPublisher}"
Name: "{group}\{#MyAppName} command line"; Filename: "{sys}\cmd.exe"; Parameters: "/K ""{app}\{#MyAppCliName}"" --help"; IconFilename: "{app}\{#MyAppCliName}"; Comment: "Open a prompt with wmlst-cli"
Name: "{group}\Licence and attribution"; Filename: "{app}\NOTICE.txt"
Name: "{group}\{cm:UninstallProgram,{#MyAppName}}"; Filename: "{uninstallexe}"
Name: "{autodesktop}\{#MyAppName}"; Filename: "{app}\{#MyAppExeName}"; Tasks: desktopicon

[Registry]
Root: HKCU; Subkey: "Environment"; ValueType: expandsz; ValueName: "Path"; \
    ValueData: "{olddata};{app}"; Tasks: addtopath; Check: NeedsAddPath(ExpandConstant('{app}'))

[Run]
Filename: "{app}\{#MyAppExeName}"; Description: "{cm:LaunchProgram,{#MyAppName}}"; Flags: nowait postinstall skipifsilent
Filename: "{app}\NOTICE.txt"; Description: "Show the licence and citation notice"; Flags: shellexec nowait postinstall skipifsilent unchecked

[UninstallDelete]
; ONLY derived data. The BLAST index and the bootstrapped BLAST+ binaries are
; rebuildable; a user's result files, exported reports and any database they
; updated in place are NOT, and are deliberately left alone. If you ever feel
; tempted to add {userdocs} or {localappdata}\IOWA-Tech\WMLST\results here:
; don't. Deleting a scientist's output on uninstall is unforgivable.
Type: filesandordirs; Name: "{app}\_internal\db\blast"
Type: filesandordirs; Name: "{app}\_internal\__pycache__"
Type: files;          Name: "{app}\*.log"
Type: dirifempty;     Name: "{app}"

[Code]
function NeedsAddPath(Param: string): boolean;
var
  OrigPath: string;
begin
  if not RegQueryStringValue(HKEY_CURRENT_USER, 'Environment', 'Path', OrigPath) then
  begin
    Result := True;
    exit;
  end;
  Result := Pos(';' + Uppercase(Param) + ';', ';' + Uppercase(OrigPath) + ';') = 0;
end;

procedure CurStepChanged(CurStep: TSetupStep);
begin
  if CurStep = ssPostInstall then
  begin
    { The BLAST+ index is derived, not shipped. Building it here would add 12 s
      to every install and would have to be redone after every database update,
      so WMLST builds it lazily on first use instead and says so in the GUI. }
  end;
end;
