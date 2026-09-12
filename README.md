# WMLSTudio

A native desktop workspace for bacterial sequence typing, designed for people using
Windows 11 without a command line. Built with Qt for Python; no browser or web server.

## What works in this preview

- A native, animated Qt workspace with a compact-screen layout and reduced-motion setting.
- Drag-and-drop sequence import; uncompressed/gzip/bzip2 FASTA and FASTQ.
- Assembly QC and exact known-allele MLST, with 162 existing PubMLST schemes.
- Local MLST/cgMLST scheme import with content fingerprints and background validation.
- Bounded FASTQ QC, explicitly labelled as a sample of the file when appropriate.
- Saved SQLite projects, sample metadata, cancellable jobs and interrupted-job recovery.
- Interactive minimum spanning forests, shared-locus evidence and adjustable visual groups.
- HYDRA JSON import: native AMR, virulence, mutation, plasmid and lineage evidence views.
- PDF/HTML reports, CSV/TSV/JSON exports and full HYDRA evidence export.

This is a research preview. It does **not** yet assemble raw reads, infer novel
alleles, run HYDRA's engines, perform phenotype prediction, or establish SeqSphere+
parity. Imported cgMLST uses exact known alleles; a missing match is not a novel
allele call. See [the capability audit](docs/STUDIO_CAPABILITIES.md) and
[HYDRA integration assessment](docs/HYDRA_INTEGRATION.md).

## Using the application

Extract the portable Windows ZIP into a writable folder, then double-click
`WMLSTudio.exe` inside `WMLSTudio`. Python, Qt and the compiled matching engine are
included. Keep the `_internal` folder beside the executable. Workspaces and local
scheme imports are saved under `Data` beside the executable. No administrator
privileges, browser, server or WSL are required by the application.

Start with **Practice project** to analyse seven clearly labelled synthetic
assemblies. For your own data, import files, select an organism-appropriate scheme
on **Samples**, and click **Analyse pending**. FASTQ files receive quality checks;
use an assembly FASTA for exact typing. Open **Compare** after typing assemblies.

Projects save automatically. **Reports → Save project copy** carries results and
metadata to another computer; original sequence files are referenced in place
and must be copied separately to rerun them. Closing during a job requests safe
cancellation. Failed/interrupted samples can be retried. Two application windows
cannot open the same project concurrently.

Every typed result retains its input SHA-256, scheme fingerprint, parameters,
software version and individual allele evidence. Original inputs are never edited;
exports refuse to replace sequence inputs or the current project.

## Develop and test

```bash
uv sync --locked
uv run python studio_scripts/stage_schemes.py --source ../wmlst/db/pubmlst
uv run python -m wmlstudio
```

Alternatively stage the pinned source snapshot with `--download`. Without staged
schemes, QC, custom scheme import and the practice project still work. This build
does not download databases at application startup.

```bash
QT_QPA_PLATFORM=offscreen uv run pytest
uv run ruff check src/wmlstudio studio_tests studio_scripts studio_packaging
QT_QPA_PLATFORM=offscreen uv run python studio_scripts/capture_desktop.py
```

On Windows, set `$env:QT_QPA_PLATFORM = "offscreen"` in PowerShell for headless
tests, then `uv run pytest`. Unset that variable before using the GUI normally.
Linux desktop development needs the standard Qt platform libraries, including
`libxcb-cursor0` for X11. Headless tests use the offscreen platform.

The command-line companion uses the same engine:

```bash
uv run wmlstudio-cli type isolate.fasta --scheme path/to/scheme -o results.json
uv run wmlstudio-cli qc reads.fastq.gz --max-reads 10000 -o quality.json
uv run wmlstudio-cli check-pair isolate_R1.fastq.gz isolate_R2.fastq.gz
uv run wmlstudio-cli compare results.json -o distances.json
```

See [Windows build instructions](docs/STUDIO_WINDOWS.md) for packaging and
acceptance checks. Windows CI is provided; configuration alone does not establish
a successful Windows 11 test.

This repository contains the new native implementation in `src/wmlstudio`.
Legacy WMLST source/history and reference databases are not included. Existing
sibling projects are reference inputs only; schemes are staged explicitly using
the commands above. Local sequence data and generated builds stay out of Git.

## Environment

Python 3.12, uv lockfile, PySide6 native widgets and pyahocorasick. Development and
real-data validation run on Linux; Windows binaries are built with Windows Python.
Database sources already present: `../wmlst/db/pubmlst` and MLSTudio local caches.

## Layout

- `src/wmlstudio`: application, analysis, project files and bundled demonstration.
- `studio_tests`: scientific, persistence and GUI tests.
- `studio_scripts`: build and reproducible validation drivers.
- `tasks`: methods, decisions, outstanding work and test hypotheses.
- `results`: generated local validation evidence, excluded from version control.

Measured validation and exact commands are recorded in `tasks/METHODS.md`,
`SUMMARY.md`, and the generated validation reports. Reference data and dependency
terms are documented in `studio_packaging/THIRD_PARTY_NOTICES.md`.
