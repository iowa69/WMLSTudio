# Changelog

All notable changes to WMLST are recorded here.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and WMLST
uses [Semantic Versioning](https://semver.org/spec/v2.0.0.html) with one project-specific
rule: **any change that can alter a typing result is a MAJOR change**, however small the
diff looks. Reproducibility is the product.

The bundled allele database has its own version (`db/VERSION.txt`), which moves
independently of the software version. Both appear in every JSON and HTML report.

## [Unreleased]

Nothing yet.

## [1.0.2] - 2026-09-11

### Fixed

- **A restored or copied database no longer looks out of date.** Index staleness
  was decided by comparing modification times, and mtimes do not survive a git
  checkout, a CI cache restore, a backup restore, or a copy between filesystems.
  Since 1.0.1 rebuilds a stale index automatically, a false positive cost the
  user a needless one-minute rebuild on start-up. The index now records a
  content fingerprint of the allele files it was built from, and staleness is
  decided by comparing that. Modification times are still used as a fallback for
  an index built by an older WMLST, and a genuine content change is still
  detected.
- A logger reference in the new stamp writer was misspelled and would have
  raised `NameError` on a read-only installation, i.e. exactly where the
  fallback it guards is needed.

## [1.0.1] - 2026-09-11

### Fixed

- **First launch after an install now works with no intervention.** The installer
  ships the allele files but not the 180 MB derived BLAST index, so the first
  start always found it missing. WMLST asked the user to go and rebuild it from
  the Database tab; it now builds the index itself, with progress, and says so.
- **A missing index is no longer misreported as a missing BLAST+ installation.**
  `blastn` runs with the database directory as its working directory, so when
  that directory did not exist the launch failed with a bare "No such file or
  directory" attributed to the `blastn` path. On first run that surfaced as
  *"WMLST needs its search engine"* and offered a 137 MB download that would not
  have fixed anything. `run_blastn` now checks the index first and raises a
  `DatabaseMissingError` that names the real problem.
- **Files passed on the command line no longer race the first-run index build.**
  `wmlst-gui sample.fna` began analysing before the index finished building and
  failed with the misleading error above. They are now held until it is ready.

## [1.0.0] - 2026-09-11

First public release. A complete Windows-native port of
[`mlst` 2.35.0](https://github.com/tseemann/mlst) by Torsten Seemann.

### Added

- **Byte-compatible CLI.** TSV, `--full`, `--csv`, `--legacy`, `--json`, `--novel`,
  `--list`, `--longlist` and `--info` reproduce real Perl `mlst` 2.35.0 output
  byte for byte, verified against a golden corpus generated from this database.
- **Bundled PubMLST database.** 162 schemes, 1,108 loci, 225,777 alleles, snapshot
  `2025-12-29`, shipped inside the installer and the wheel. Nothing to download
  before the first run.
- **Graphical interface** (`WMLST.exe`, `wmlst-gui`): drag and drop assemblies,
  sortable results, and export to TSV, CSV, JSON or a self-contained HTML report.
  Drag-and-drop uses `tkinterdnd2` when it is installed and degrades to a Browse
  button when it is not.
- **Self-contained HTML report** (`--html`): one file, no CDN, no web fonts, no
  network access, safe to email.
- **NCBI BLAST+ bootstrap** (`--bootstrap-blast`): downloads the official 2.17.0
  archive over HTTPS, verifies its MD5, and extracts only the ~35 MB actually
  needed into `%LOCALAPPDATA%\IOWA-Tech\WMLST\blast\`. No administrator rights.
- **Database updater** (`--update-db`, `wmlst-update-db`, and the GUI Database tab):
  content-hash change detection rather than date comparison, per-scheme atomic
  commits, resumable, with rollback and offline import/export bundles.
- **Pure-Python `any2fasta`** replacement, reproducing the same conversions and the
  same failure strings without a Perl dependency.
- **Windows correctness work** that has no upstream equivalent: CRLF stripped from
  every `blastn` output row before `sstrand` is compared, UTF-8 forced on all
  streams, no `NamedTemporaryFile`, long-path and non-ACP filename support, no
  console window under `pythonw`, and a stale-index guard.
- **Extras that never change a compatible output**: `--jobs`, `--blast-timeout`,
  `--evidence-tsv`, `--html-evidence`, and `--repair-locus-ids`
  (off by default, and its help text says results will not match upstream).
- `MLST_DBDIR` is honoured alongside `WMLST_DBDIR`, so scripts written for upstream
  `mlst` keep working.
- Packaging: PyInstaller one-folder build producing `WMLST.exe` (windowed) and
  `wmlst-cli.exe` (console) from one shared payload, an Inno Setup per-user
  installer, a portable zip, a wheel and an sdist, all with `SHA256SUMS.txt`.

### Deliberate differences from upstream

Full register in `docs/ARCHITECTURE.md` section 15. None of them changes a number
except `--repair-locus-ids`, which is off by default.

- Ordering is deterministic everywhere. Upstream iterates randomised Perl hashes,
  so exact score ties are decided by a coin flip and cannot be reproduced
  run-to-run; WMLST breaks ties by scheme name.
- `err()` prints `ERROR: `, not upstream's `ERRPR: ` typo.
- A failing input file no longer aborts a whole batch; it is reported and skipped.
- Relative `--fofn` entries resolve against the FOFN's own directory, and blank
  lines in a FOFN are skipped rather than raising `Unable to read from ''`.
- JSON key order is fixed rather than randomised.

### Known limitations

- **The Windows binaries are not code-signed.** SmartScreen will show
  *"Windows protected your PC"* on first run; **More info → Run anyway**. This is a
  reputation signal, not a detection. `SHA256SUMS.txt` and the public CI build log
  are the provenance record.
- MLST needs each locus on a single contig. Highly fragmented assemblies produce
  partial (`?`) calls or none at all; this is inherent to the method.
- `--info` reports `DATE` as `Unknown` for all 162 shipped schemes, exactly as
  upstream does with this database.

[Unreleased]: https://github.com/iowa69/WMLST/compare/v1.0.0...HEAD
[1.0.0]: https://github.com/iowa69/WMLST/releases/tag/v1.0.0
