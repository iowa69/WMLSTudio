# WMLSTudio

A local microbial-genomics workbench for non-commercial research on Windows:
import, assemble paired reads, type, investigate AMR evidence, compare saved
isolates and report a defined cohort.

Native dark Qt desktop. No browser server, WSL, Docker or separate Python
installation for the portable package.

*IOWA-BioTech — Giovanni Lorenzin*

## Download and retest

### [⬇ Download WMLSTudio — Windows portable ZIP](https://github.com/iowa69/WMLSTudio/releases/download/v0.2.0-workbench.1/WMLSTudio-Windows-x64.zip)

**v0.2.0-workbench.1 · Windows x64 · 147.6 MB ZIP · Non-commercial research**

[Release notes](https://github.com/iowa69/WMLSTudio/releases/tag/v0.2.0-workbench.1)
· [SHA-256 checksum](https://github.com/iowa69/WMLSTudio/releases/download/v0.2.0-workbench.1/WMLSTudio-Windows-x64.zip.sha256)
· [Build/test results](https://github.com/iowa69/WMLSTudio/actions/workflows/studio.yml)
· [Previous 0.1 release](https://github.com/iowa69/WMLSTudio/releases/tag/v0.1.0-preview.1)

This repository is private: sign in with an account that has access to download.
The release notes distinguish observed Linux, native Windows CI and Wine checks
from the still-required clean Windows 11 desktop acceptance test.
The published binary was built from `e5d5f72`: **466 tests passed on each hosted
platform**, followed by frozen-executable checks in
[this successful build](https://github.com/iowa69/WMLSTudio/actions/runs/34713168070).

> Research software, not a validated diagnostic device or an established
> SeqSphere+ equivalent. Review the [capability audit](docs/STUDIO_CAPABILITIES.md).
> The application is unsigned. Check the source and checksum before deciding to
> run it; do not disable Defender or SmartScreen. A checksum establishes file
> integrity, not safety or clinical validity.

## Three steps to the desktop

### 1. Extract the entire folder

Choose **Extract All**, use a writable folder such as Documents, then open
`WMLSTudio/WMLSTudio.exe`. Keep every executable and the `_internal` folder
together. Do not launch from inside the ZIP or move only the main executable.
The extracted application occupies about 462 MB, before your projects and inputs.

Keep your existing data backed up. Open a **copy** of an older project for the
first retest: the workbench upgrades its project schema, and the old 0.1
application cannot read the upgraded project.

### 2. Explore the practice project

Use **Help → Practice project**. Seven explicitly synthetic assemblies let you
test navigation, sample selection, comparisons, graph options and reports.

![Native WMLSTudio workbench with the synthetic practice cohort](docs/images/wmlstudio-workbench.png)

*Native Windows application; synthetic demonstration data, not a clinical cohort.*

### 3. Import and review your own samples

Drop FASTA/FASTQ files or a folder into Overview. The import dialog supports:

- Automatic classical-MLST discovery, manual organism/scheme, or unknown.
- Per-file assignments and apply-to-selected/apply-to-all controls.
- Optional local input copies, organism/ST organization and ST filename suffixes.

On **Samples**, review the selection and choose **Analyse pending…** or
**Analyse selected…**. The launch dialog explains the workflow and optional
paired-short-read assembly and HYDRA analyses. Original input files stay unchanged.

Accepted formats: `.fasta`, `.fa`, `.fna`, `.fastq`, `.fq`, including
gzip and bzip2. Unassembled reads receive labelled sampled QC. The paired-read
workflow validates every pair before native SKESA; mate records and read
provenance remain linked to the resulting assembly.

## One connected workflow

| Workspace | What you can do |
| --- | --- |
| Overview | Browse the library by organism/ST, open projects and research saved isolates |
| Samples | Assign workflows, select cohorts, inspect QC/typing, edit annotations, create collections and inspect history |
| Compare | Choose a cohort and scheme snapshot, call additional cgMLST/wgMLST profiles, reuse saved profiles, style/export the minimum spanning forest |
| Scheme library | Import local schemes, browse online catalogs, install versioned snapshots and create a local ad-hoc cohort scheme |
| HYDRA insights | Run native assembly AMR searches, map external reports explicitly, inspect all features and primary-aware AMR matrices |
| Reports | Include/exclude isolates, highlight investigation groups, export PDF/HTML/CSV/TSV/JSON and portable profile bundles |
| Settings and help | Reduced motion, local data location, workflow guidance and scientific limitations |

Menus provide the canonical actions; **Ctrl+K** opens command search.
**Ctrl+R** analyses selected samples. **Alt+1…7** switches workspace pages.

## Typing and comparison

Classical STs require a complete, unambiguous registered profile. Automatic
scheme matching is **provisional organism evidence**, not independent species
confirmation or a contamination assessment.

Additional cgMLST/wgMLST results do not overwrite classical MLST. Large schemes
use exact-first matching with complete-CDS guarded novel calling. Novel sequences
receive full SHA-256 identifiers; they are not silently assigned public allele
numbers. Missing, mixed, duplicated and ambiguous loci remain explicit.

Compare only the intended cohort against one identical scheme snapshot. Shared
loci and excluded pairs remain visible. Graph options include metadata/manual
colors, labels, draggable layout, identical-genotype pies and cluster halos.
Export PNG, SVG, GraphML, MST-topology Newick, pairwise JSON or a distance-matrix TSV.
A close edge is not proof of transmission.

## HYDRA is an analysis engine, not just an imported page

The portable assembly runtime uses pinned HYDRA 1.4.0 and native BLAST+ 2.17.0.
An NCBI nucleotide/protein/mutation starter snapshot is supplied; additional
databases require explicit installation and review of provider terms.

Run-plan controls expose nucleotide identity/coverage, translated protein search,
protein thresholds, CPU allocation and organism-specific mutation evidence.
Protein “complete coverage” distinguishes complete from partial evidence; it
does not discard every shorter hit. Unknown organisms do not silently receive a
species-specific mutation catalog.

Results join typing, QC and metadata through stable sample IDs. Ambiguous names
in an imported report require explicit mapping. Primary/secondary detections
are distinguished; “no report” is unknown, not a negative result.
Gene/mutation evidence is not measured susceptibility.

See the [engine integration and real-data checks](docs/HYDRA_INTEGRATION.md).

## Reuse, storage and reports

Projects save locally and automatically. Managed copies can be organized by
organism and ST; optional ST suffixes apply to those copies, never originals.
Read pairs retain both source records and the derived assembly's provenance.

After moving an input, select its sample and use **Samples → Relink input
(same bytes)…**. The complete file SHA-256 must match recorded evidence; profiles
remain unchanged. Paths are not rewritten automatically when folders or drive
letters change.

**Research saved library** searches indexed projects without re-reading genomes.
Bring selected frozen profiles into the active project, or import a collaborator's
profile table/bundle. No FASTA is required to compare existing profiles.

TSV preserves full-length `SHA256_`/`NOVEL_` identifiers, but imported hash IDs
remain unverified: the table does not establish sequence or CDS validation.
Use a full portable profile bundle for lossless transfer of validated novel-call
evidence.

Reports use an explicit cohort and optional highlighted groups. JSON and portable
bundles retain full additional profiles and linked evidence. PDF summarizes
additional schemes and AMR evidence. Back up the adjacent `Data` directory,
project-specific managed-input folders and original sequences. Local storage is
not an encrypted clinical-record system.

## References and rights

The classical reference cache contains 162 staged schemes; it is a snapshot, not
a promise that every scheme is current. Online providers can require credentials
or restrict access by submission date. Updates create new snapshots and retain
old ones; analysis never silently downloads references or uploads genomes.

**cgMLST.org database contents are not included in the ZIP.** Individual research
downloads are subject to its [server policy](https://www.cgmlst.org/serverpolicy.html);
database-driven product/service use requires permission. PubMLST has separate
[terms](https://pubmlst.org/terms-conditions). Confirm rights for your actual use
and sharing. Software licensing does not grant database redistribution rights.

The complete Kleborate/Kaptive, AMRFinderPlus, agr/SCCmec/spa, MOB-recon and
abricate execution stack, direct-read AMR/pileup, fastp/SPAdes, independent species
confirmation and clinical validation remain unfinished. The
[scientific workflow contract](docs/MICROBIOLOGY_WORKFLOWS.md) makes those gaps
explicit.

## Verify the download

```powershell
Get-FileHash .\WMLSTudio-Windows-x64.zip -Algorithm SHA256
```

Compare it with the checksum asset linked above.

## Develop and test

```bash
uv sync --locked --group dev
uv run python studio_scripts/stage_schemes.py --download
uv run python -m wmlstudio
QT_QPA_PLATFORM=offscreen uv run pytest
uv run ruff check src/wmlstudio studio_tests studio_scripts studio_packaging
```

Native engine staging and Windows packaging are documented in
[STUDIO_WINDOWS.md](docs/STUDIO_WINDOWS.md). On Windows use PowerShell's
`$env:QT_QPA_PLATFORM = "offscreen"` for headless tests; unset it for normal use.

The source lives in `src/wmlstudio`; tests in `studio_tests`; reproducible build
and validation drivers in `studio_packaging` and `studio_scripts`. Original
sibling projects and local sequence-bearing artifacts are not committed.

## Credits and licenses

Built by Giovanni Lorenzin / IOWA-BioTech. Scientific baseline:
[MLSTudio](https://github.com/iowa69/mlstudio); reference snapshot staging:
[WMLST](https://github.com/iowa69/WMLST); AMR engine:
[HYDRA](https://github.com/iowa69/hydra). Cite the original schemes and scientific
tools used in your analyses.

WMLSTudio's source retains its stated license. The combined portable distribution
includes components under other terms, notably GPL Pyrodigal and AGPL-containing
SKESA code. Corresponding source/build materials and notices accompany those
components. See [third-party notices](studio_packaging/THIRD_PARTY_NOTICES.md);
the complete binary package must not be described as MIT-only.
