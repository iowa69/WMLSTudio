# WMLST on Windows

Everything on this page is about the Windows build specifically. For what WMLST *is*,
start at the [README](../README.md).

---

## 1. SmartScreen, honestly

**WMLST is not code-signed.**

The first time you run `WMLST-1.0.0-win64-setup.exe`, Windows will almost certainly show:

> **Windows protected your PC**
> Microsoft Defender SmartScreen prevented an unrecognised app from starting.

To continue: **More info → Run anyway**.

### What that message actually means

SmartScreen scores files on *reputation*, not on content. A file it has not seen before,
from a publisher without an Authenticode certificate, gets this warning by default. It
is not a malware detection, and clearing it is not a matter of "fixing" anything in the
build: it goes away once enough people have downloaded the file, or the day someone pays
for a certificate (several hundred euro a year, renewed annually — which an independent
tool cannot simply assume).

### What we do instead

- Every release is built by a **public GitHub Actions run**, and the release notes link
  the build log for that exact commit.
- **`SHA256SUMS.txt`** is published next to the binaries.

```powershell
Get-FileHash .\WMLST-1.0.0-win64-setup.exe -Algorithm SHA256
```

Compare the result with the line in `SHA256SUMS.txt`. If they match, you have exactly
what CI produced from exactly that source.

### What we will never tell you to do

Disable SmartScreen. Disable Defender. Edit the registry. Install a self-signed
certificate (SmartScreen ignores untrusted roots, so it changes nothing except making
the build look evasive). If a page claiming to be about WMLST suggests any of these,
it is not us.

### Antivirus false positives

An unsigned PyInstaller executable is a well-known false-positive trigger. WMLST is
deliberately **not** UPX-compressed, because packing makes the heuristics worse and
corrupts the MSVC runtime DLLs CPython ships. If your endpoint protection quarantines
`WMLST.exe`, add an exclusion for the install folder or use the Python package instead
(`pip install wmlst`).

---

## 2. What the installer does

| | |
| :-- | :-- |
| Install location | `%LOCALAPPDATA%\Programs\IOWA-Tech\WMLST` (per-user) |
| Administrator rights | **none required** |
| Registry | only a `PATH` entry, and only if you tick that optional task |
| Start Menu | `IOWA-Tech → WMLST` |
| Desktop shortcut | optional, off by default |
| Size on disk | ~195 MB, of which ~116 MB is the allele database |

Two executables are installed from one shared payload:

- **`WMLST.exe`** — the graphical application, built windowed so it never flashes a
  console.
- **`wmlst-cli.exe`** — the command line, output-compatible with `mlst` 2.35.0.

`LICENSE.txt`, `NOTICE.txt` and `SMARTSCREEN.txt` are installed alongside them. That is
a GPLv2 §1 requirement, not decoration.

Uninstalling removes the program and the derived BLAST index. It **does not** remove
your results, your exported reports, or a database you updated yourself.

---

## 3. NCBI BLAST+ on first run

WMLST drives `blastn`. The repository and the installer do not ship NCBI's binaries;
the application fetches them once:

```powershell
wmlst-cli --bootstrap-blast
```

- Source: the official NCBI HTTPS endpoint, pinned to **2.17.0** (never `LATEST/`, so an
  upstream bump cannot silently change your results).
- Archive: 143,400,333 bytes, MD5-verified before anything is unpacked.
- Only ~35 MB is extracted: `blastn.exe`, `makeblastdb.exe`, `nghttp2.dll` and the NCBI
  licence notices.
- Destination: `%LOCALAPPDATA%\IOWA-Tech\WMLST\blast\` — per-user, no elevation, and
  deliberately outside any OneDrive-synced folder.

Then the search index is built once, in about twelve seconds:

```powershell
wmlst-cli --make-blast-db
wmlst-cli --check
```

### If BLAST fails immediately

A `0xC0000135` (`STATUS_DLL_NOT_FOUND`) or `0xC0000142` exit code with no output means a
missing **Microsoft Visual C++ runtime**, not a WMLST fault. `blastn.exe` imports
`MSVCP140.dll`, `VCRUNTIME140.dll` and `VCRUNTIME140_1.dll`, and `MSVCP140.dll` is not
guaranteed present on a clean Windows 11. Install the
[VC++ 2015-2022 x64 redistributable](https://aka.ms/vs/17/release/vc_redist.x64.exe).

---

## 4. Windows-specific behaviour worth knowing

- **Output is always LF-terminated**, never CRLF, so a file written on Windows is
  byte-identical to one written on Linux. This matters more than it sounds: `blastn`
  writes CRLF on Windows, and the last field of its tabular output is the strand, so an
  unstripped `\r` would make `"plus\r"` never compare equal to `"minus"` and would emit
  un-reverse-complemented novel alleles with wrong MD5s.
- **Non-ASCII and long paths work.** Filenames outside the active code page, paths over
  260 characters, spaces, and shell metacharacters such as `'&^%!()` are all handled;
  nothing is passed through a shell.
- **`%TEMP%` is not used for scratch files.** Working files live beside the output, are
  flushed, `fsync`-ed and closed before `blastn` is launched, and are cleaned up on
  cancellation — no orphaned `blastn.exe` is left behind.
- **No console windows.** Child processes are created with the no-window flag, so a
  twenty-genome scan from the GUI shows exactly zero flashing black boxes.
- **`MLST_DBDIR` is honoured** alongside `WMLST_DBDIR`, so a script written for upstream
  `mlst` keeps working.

---

## 5. Portable use, and locked-down machines

If you cannot install software, download `WMLST-1.0.0-win64-portable.zip`, unblock it
(right-click → Properties → **Unblock**), and extract it anywhere you can write —
including a USB stick. Run `WMLST.exe` from the extracted folder.

For a fully air-gapped machine, bootstrap BLAST+ and update the database on a connected
machine, then move a verified snapshot across with `--export-bundle` / `--import-bundle`.

---

*WMLST — IOWA-Tech · Giovanni Lorenzin · GPL-2.0-only.
Port of [`mlst`](https://github.com/tseemann/mlst) by Torsten Seemann.
Allele data © PubMLST — cite Jolley, Bray & Maiden (2018), PMID 30345391.
Research use only. Not a medical device.*
