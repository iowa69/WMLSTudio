# Portable WMLSTudio for Windows 11

WMLSTudio uses native Qt widgets. The portable folder includes the Python
interpreter, Qt libraries, the compiled exact-match engine, and a static scheme
snapshot. The end user does not install Python, WSL, Docker, Conda, or a browser
server. Windows binaries require a Windows Python build runtime. The supported
CI recipe runs on Windows; the same Windows runtime can be exercised under Wine
for an experimental build. Linux or Wine success does not establish Windows 11
compatibility.

## Running a Windows build

Extract the entire ZIP into a writable folder, such as
`Documents\WMLSTudio`, then double-click `WMLSTudio.exe`. Keep `_internal` next
to the executable. Do not run the executable inside the ZIP or move it alone.
The frozen application stores its local workspace in the adjacent `Data`
directory. Back up that directory and the original sequence files. Moving the
application does not move external FASTA/FASTQ files referenced by a project.

Import assembled FASTA files, select the correct organism's scheme, then run
typing. FASTQ files receive a clearly labelled sampled quality summary and must
be assembled before this engine can type them. HYDRA JSON import displays
previously computed AMR/virulence evidence; it does not run HYDRA on Windows.

`WMLSTudio-CLI.exe --help` opens the optional command-line companion. Neither
executable downloads references at runtime. The bundle's manifest records its
reference snapshot so results can be tied to the data actually used.

## Build on Windows

Developers need Git, Python 3.12, and [uv](https://docs.astral.sh/uv/). From the
checked-out repository in PowerShell:

```powershell
.\studio_packaging\build_windows.ps1
```

This installs dependencies from `uv.lock`, explicitly downloads the pinned WMLST
reference snapshot, preserves version-matched dependency license texts, runs the
test suite, freezes both entry points with
PyInstaller, checks the frozen desktop and CLI, and creates
`dist\WMLSTudio-Windows-x64.zip` with a SHA-256 displayed in the terminal.
The data source is pinned to commit
`6cad46ffd9dfddfaa55f7993cf391f80a70556d3`; no moving branch is used for downloads.

For an offline reference cache, while Python dependencies are already available:

```powershell
.\studio_packaging\build_windows.ps1 -SchemeSource C:\reference-data\pubmlst
```

The cache must contain one directory per scheme, with allele FASTA files and
optional ST profile tables. This mode records the local source revision when
available and hashes every copied file. An existing staging directory containing
stale unrelated data is refused; use a fresh checkout or staging destination.

The CI workflow `.github/workflows/studio.yml` runs on Ubuntu and Windows and
attaches development packages and test reports as workflow artifacts. It does
not create releases or publish a download site. An unsigned development package
is not a code-signed production release.

## Validation and distribution limits

The frozen smoke test covers launch, Qt resource discovery, screenshot generation,
and clean exit; the test suite covers scientific edge cases and UI workflows.
Release validation still requires a clean Windows 11 machine without Python,
offline startup, paths containing spaces and non-ASCII characters, high-DPI
displays, cancellation during large jobs, and recovery after interrupted jobs.
Track performed checks in `artifacts`; do not substitute planned CI for observed
Windows results.

Read [the capability audit](STUDIO_CAPABILITIES.md) and
[third-party notices](../studio_packaging/THIRD_PARTY_NOTICES.md). Reference-data
redistribution rights must be established for a public release; preserved
download metadata is provenance, not an independent license determination.
