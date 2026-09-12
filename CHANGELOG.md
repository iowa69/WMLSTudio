# Changelog

All notable changes to WMLST are recorded here.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and WMLST
uses [Semantic Versioning](https://semver.org/spec/v2.0.0.html) with one project-specific
rule: **any change that can alter a typing result is a MAJOR change**, however small the
diff looks. Reproducibility is the product.

The bundled allele database has its own version (`db/VERSION.txt`), which moves
independently of the software version. Both appear in every JSON and HTML report.

## [Unreleased]

### Changed

- **A tie between two schemes is now broken by allele registry depth, not by
  alphabetical order.** *(This can alter a typing result, so it is a MAJOR
  change under the rule above.)* WMLST's scheme score is coarse — a 7-locus
  scheme has eight reachable values — so two schemes reaching 100 on the same
  assembly is routine. The Achtman *E. coli* scheme is built from housekeeping
  genes conserved across the Enterobacteriaceae and carries off-species alleles,
  so a *Klebsiella pneumoniae* draft scores 100 in both `klebsiella` and
  `ecoli_achtman_4`. Upstream `mlst` breaks that tie with a coin flip; WMLST used
  to take the alphabetically first scheme name — reproducible, but no more
  meaningful, and it reported five *K. pneumoniae* as `ecoli_achtman_4` ST 14464
  instead of the carbapenem-resistant ST 258 / ST 11 clones.

  Where, and only where, the scores are exactly equal, WMLST now prefers the
  scheme whose called alleles sit lowest in each locus' allele registry.
  PubMLST issues allele numbers in order of first observation, so a scheme's own
  species keeps matching the long-established alleles it has been depositing
  since the scheme opened, while an off-species coincidence can only match rare,
  late-registered variants. The scheme name remains the final key, so the result
  stays fully deterministic. Measured over a 210-genome labelled RefSeq corpus
  (30 assemblies each of *S. aureus*, *A. baumannii*, *E. coli*, *E. faecium*,
  *K. pneumoniae*, *P. aeruginosa*, *E. cloacae*): 17 genomes tie, 17/17 now
  resolve to the correct species where 0/17 did before, and the scheme call is
  correct on 210/210 rather than 193/210. Only the 17 tied rows changed; no
  other call moved, and the 17 byte-identical goldens are untouched.

### Added

- **A tie is now visible instead of silently resolved.** The upstream
  `WARNING: a(st)==b(st) score=N` line on stderr is unchanged, and
  `SampleResult.tied` now carries the tied alternatives, so the HTML report names
  both schemes with their sequence types ("2 schemes fit this assembly equally
  well at score 100 — klebsiella ST 258 (reported) and ecoli_achtman_4
  ST 14464"), flags the tied rows in the runner-up table, and the GUI shows the
  same on the summary card, in an expansion row and in the status line. The
  compat TSV/CSV/JSON row formats are byte-identity surfaces and are unchanged.

## [1.2.1] - 2026-09-12

### Changed

- **The download link works.** It pointed at a GitHub Pages site that had never
  been enabled, so every Pages run had failed and the URL returned 404. Pages is
  enabled and deployed, and the README now leads with the releases page, which
  cannot break, offering the landing page as the alternative.
- **The Tree page is documented** in the README, with a screenshot of each view,
  what the numbers on the edges mean and what a clonal group is.
- **No third-party citations on any user-visible surface.** The banner, `--help`,
  the About box, the Settings footer, the HTML report and the landing page no
  longer name or cite other projects. The Settings line claiming the defaults
  matched another tool is gone; it now simply states the reference defaults.
  The copyright notices the licence requires remain in the source headers and in
  NOTICE, which travels inside the download.

## [1.2.0] - 2026-09-12

### Added

- **A Tree page.** Once two or more isolates share a scheme, WMLST draws a minimum
  spanning tree over their allele profiles. Two views: *Clonality*, with dots
  coloured by clonal group and sized by how many isolates share the sequence
  type, and *Labelled tree*, with every isolate named and its ST shown, coloured
  by ST. Edges carry the allelic distance. Zoom, pan, click an isolate to read
  its full profile and nearest neighbours, and save the picture as a PNG.
  Distances are computed only between isolates of the same scheme, ignoring loci
  missing in either, and the view states how many loci were compared.
- Results filter (Ctrl+F), row counts, and a View menu with Ctrl+1-4.

### Changed

- **The product is now IOWA-BioTech.** An existing installation keeps using the
  BLAST+ it already downloaded: the previous vendor folder is still searched
  before the new one, so upgrading never re-downloads 137 MB.
- **A failing scheme no longer aborts a database update.** The other schemes
  still update, the failed ones are left untouched on disk, and the run ends with
  a summary naming each failure and why.
- The animation is smoother and carries the IOWA-BioTech wordmark; panels are
  outlined so regions of the window read as distinct.
- README is now a user guide: download, three steps, what the output looks like
  on screen and as files. The copyright notices live in NOTICE, where the licence
  expects them.
- Releases publish the portable zip only for now; the installer step is paused
  behind a switch rather than removed, and v1.1.1's installer remains available.

## [1.1.1] - 2026-09-12

### Fixed

- **The Database tab showed no sequence-type or allele counts.** `SchemeInfo`
  had no fields for them, so the deep catalogue scan the GUI already ran filled
  in nothing and both columns stayed empty forever. They are now populated for
  160 of the 162 schemes (two ship no profile table), and `wmlst --info` remains
  byte-identical to the reference tool across all 163 rows.
- The tie notice on the status bar ran off the right edge and lost its ending.
  The status line now names both schemes and the file; the full explanation
  stays on the result card and in the report, where there is room for it.
- Two tests used POSIX-only APIs (`os.getuid`, and `os.chmod` on a directory)
  and could not run on Windows. They now assert the same product behaviour on
  every platform using a path whose parent is a regular file, and keep the
  mode-bit assertions where mode bits exist.
- The shutdown tests reported a killed child as surviving on Windows.
  `os.kill(pid, 0)` is not a liveness probe there — CPython maps it to
  `TerminateProcess` for every signal but CTRL_C/CTRL_BREAK, so it kills what it
  touches, and an exited process keeps an openable PID while anyone holds a
  handle to it. The probe now asks whether the process object is signalled.
  The shutdown behaviour itself was correct; only the measurement was wrong.

## [1.1.0] - 2026-09-12

### Changed - this release can alter a typing result

- **Scheme ties are no longer broken alphabetically.** When two schemes score
  identically, WMLST now prefers the one whose called alleles sit lowest in each
  locus' own allele registry, measured as a percentile so that a scheme's age and
  size do not decide the call. Scheme name remains the final key, so selection
  stays fully deterministic.

  This fixes real misidentifications. Measured on 210 labelled NCBI RefSeq
  genomes, 30 each across seven ESKAPE organisms: schemes matching the known
  organism went from 193/210 to **210/210**, and top-score ties resolved to the
  correct species went from **0/17 to 17/17**. Five *Klebsiella pneumoniae*
  genomes - two of them ST 258, the dominant carbapenem-resistant clone - were
  being reported as *Escherichia coli* ST 14464, because the Achtman scheme uses
  housekeeping genes conserved across Enterobacteriaceae. No non-tied call
  changed, and all 17 golden outputs remain byte-identical.

### Added

- **The organism is now the headline.** Results name the genus and species in
  italic binomial form, the scheme and its description, the locus count, and the
  primary reference with a PubMed link and a link to the authoritative PubMLST or
  Pasteur record. A new bundled table covers all 162 schemes: 158 carry a genus,
  and 75 carry a verified citation. A cryptic "abaumannii_2" now reads
  *Acinetobacter baumannii*.
- **Ties are shown, not silently resolved.** Two schemes fitting equally well is a
  real ambiguity; both are named, with both STs, on the card, in the table and in
  the HTML report.
- **Automatic performance tuning.** At start-up WMLST sizes itself to the machine:
  one file at a time per four cores, four BLAST threads each - 8 cores analyse 2
  files at once, 16 cores 4, 32 cores 8. Overridable in Settings. CLI defaults are
  unchanged for upstream parity.
- **Portable database.** A checkbox keeps the database beside WMLST.exe so a
  portable copy updates itself in place, with the resolved path always shown.
- **A redesigned interface.** New light, dark and high-contrast palettes, a real
  type scale, an 8px spacing grid, flat buttons with proper states, hairline cards
  with drawn rounded corners, and custom drawn checkboxes. The drop zone carries a
  circular chromosome whose seven locus arcs illuminate in turn, with drifting
  cocci and rods; it pauses when idle or running and honours reduced-motion.
- **Database tab controls**: Select all, Deselect all, a live selection count, and
  larger checkboxes that respond to click, Space and keyboard focus.
- **A one-click download page** that always resolves to the newest portable zip.

### Fixed

- **WMLST can always be closed.** Closing the window during a run left the process
  alive and only Task Manager could end it: ThreadPoolExecutor workers are
  non-daemon, concurrent.futures joins them at exit, and the BLAST children were
  never killed because terminate_all() was defined but never called. Closing now
  stops the engine, kills every child and exits - measured at 0.05s to kill the
  child and 1.04s to exit, with nothing left behind.
- Status is drawn as a shape rather than a font glyph, which rendered as an empty
  box on systems without the dingbat; the status word and colour are unchanged, so
  status is still never conveyed by colour alone.
- Several labels rendered as mojibake under some locales.
- db/scheme_refs.tsv is packaged in the wheel and sdist; without it the organism
  names silently disappeared.
- Shutting down no longer tells a user who just closed the window to install BLAST+.
- The CLI executable no longer carries the application icon, and the installer no
  longer creates a command-line Start Menu entry that a novice might click.

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
  needed into `%LOCALAPPDATA%\IOWA-BioTech\WMLST\blast\`. No administrator rights.
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
