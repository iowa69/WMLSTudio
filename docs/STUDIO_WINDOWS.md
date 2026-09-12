# Portable WMLSTudio for Windows 11

The portable Windows build contains native Qt, Python, the typing engine,
Pyrodigal, a dedicated HYDRA worker, official Windows BLAST+, a verified native
SKESA tool package, classical reference schemes and an NCBI AMR starter snapshot.
No WSL, Docker, Conda, browser server or separately installed Python is required
by the end user.

This is a non-commercial research workbench. Windows 11 is the target platform;
hosted Windows tests and Windows binaries under Wine are distinct from a clean
Windows 11 desktop acceptance test. See the release's observed validation results.

## Open and retest

Extract the entire ZIP to a writable folder such as `Documents\WMLSTudio`,
then double-click `WMLSTudio.exe`. Keep all files and `_internal` together.
The other executables are the optional command-line utility and internal HYDRA
worker; normal users do not need to launch them.

The local workspace is stored under adjacent `Data`. Back up that directory,
your project-specific managed-input folders and original files. Project data
save automatically. Open a copy of an old project for your first retest because
the new schema cannot be opened by the 0.1 application.

Input paths are not automatically rewritten when files or an application/project
folder move, including moves between drive letters. Frozen profiles remain usable
if inputs are missing; rerunning requires the input bytes to be available.

To reconnect a moved input, select one sample and choose **Samples → Relink input
(same bytes)…**, then select the relocated original or a byte-identical copy.
WMLSTudio checks the entire file against its recorded SHA-256 before changing the
path. Existing profiles and evidence stay unchanged; relinking does not retype
the sample. A different hash or missing recorded fingerprint is rejected—import
that file as a new sample instead. Recompression or edited headers can change the
file hash even when sequence letters appear unchanged. A copy outside its known
managed location is treated as an external input, not automatically made eligible
for managed-file cleanup.

Portable profile bundles intentionally carry results and metadata without
sequence files; they do not locate or restore missing FASTA/FASTQ files.

Use **Help → Practice project** for the synthetic demonstration. For real data,
import and assign workflows, review the launch plan, and select optional paired
read assembly/HYDRA. The SKESA build targets x64/SSE4.2. Its memory setting must
exceed the engine's fixed reserve; the default is 8 GB. Do not allocate more
memory/threads than the computer can reasonably supply.

## Reproducible Windows build

Developers need Git, Python 3.12, uv and a successfully validated native SKESA
artifact. The dedicated **Native Windows SKESA** workflow builds the pinned
source, applies the reviewed Windows compatibility patch, packages its DLL
closure and corresponding source, then performs genuine assembly/adapter tests.

The main `.github/workflows/studio.yml` invokes that reusable build before its
native Windows test/freeze job. It downloads the verified tool artifact into
the resource staging area; missing SKESA or starter references are build errors,
not silently disabled features.

For a local PowerShell build:

```powershell
.\studio_packaging\build_windows.ps1 -SkesaSource C:\validated-tools\skesa
```

See `Get-Help .\studio_packaging\build_windows.ps1` or the script parameters
for local scheme and HYDRA reference cache options. Dependencies are installed
from `uv.lock`. The reference staging command uses the pinned WMLST commit
`6cad46ffd9dfddfaa55f7993cf391f80a70556d3`, preserves per-file hashes and
refuses inconsistent staging directories.

BLAST+ archives are SHA-256 pinned. An explicit build-time NCBI download records
the provider version before/after download and all companion file hashes.
Already-staged verified references can be reused. Runtime analysis never
downloads databases; desktop reference updates require explicit user action.

The build preserves dependency notices, runs tests, freezes the desktop/CLI/HYDRA
entry points, and checks actual frozen typing, HYDRA and SKESA execution before
creating `dist/WMLSTudio-Windows-x64.zip`. Source/material provenance accompanies
the native tools. A build recipe is not a claim that any particular run passed.

## Validation gates

Required package checks include offline startup, native menu repainting,
reference discovery, genuine positive-control typing/AMR, real native assembly,
spaces/Unicode/comma paths, cancellation, complete-pair validation, interrupted
jobs, input immutability, stale-evidence handling, sample selection and report scope.

Local real-data checks and synthetic controls are kept distinct. A small number
of successful assemblies or exact cross-platform calls do not establish clinical
sensitivity, specificity, contamination detection or organism-module equivalence.

Clean Windows 11 acceptance still needs to cover high-DPI displays, ordinary
restricted user accounts, drive relocation, storage failures and realistic larger
cohorts without a development environment installed. The package is unsigned.
Do not disable operating-system protections to run it.

Read the [capability audit](STUDIO_CAPABILITIES.md), [microbiology workflow
contract](MICROBIOLOGY_WORKFLOWS.md), [HYDRA integration](HYDRA_INTEGRATION.md)
and [third-party notices](../studio_packaging/THIRD_PARTY_NOTICES.md).
cgMLST.org database contents are not redistributed in this package. Software
licenses do not transfer database rights.
