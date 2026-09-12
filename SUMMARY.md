# WMLSTudio 0.2 — native research workbench

The revision implements native paired-read assembly, conservative automatic
MLST discovery, additional cgMLST/wgMLST profiles and local ad-hoc schemes,
actual HYDRA assembly analysis, managed storage, saved-library reuse, explicit
comparison/report cohorts and a dark native Qt interface.

No web server, WSL, Docker or end-user Python installation is required by the
portable Windows package. Non-commercial research is the intended use.
cgMLST.org database contents and local sequence inputs are not bundled.

## Observed scientific and engineering checks

- Final hosted [Windows and Linux build](https://github.com/iowa69/WMLSTudio/actions/runs/34713168070):
  **466/466 tests passed on each platform**, zero failures or skips; 86.7 s on
  Windows Server 2025 and 46.9 s on Linux. Both frozen builds passed their gates.
  The Windows build also passed native SKESA, app-local DLL closure, actual
  Segoe UI glyph checks and the isolated bundled-CDS ground-truth probe.
- Final Linux suites: **466 tests passed**, none skipped: 386 core/build checks,
  including real BLAST integration, and 80 native UI/workflow checks. Evidence:
  `artifacts/workbench-final-core.xml` and `artifacts/workbench-final-ui.xml`.
- S. epidermidis: ST184; three P. aeruginosa assemblies: ST155, ST2952 and
  ST1858, with all 21 allele calls matching the preserved WMLST comparator.
- Real E. faecium 1,423-locus control: 1,410 exact, one CDS-validated novel,
  seven missing, four ambiguous and one mixed locus; mixed result excluded
  from comparison. Novel sequence independently checked against coordinates.
- HYDRA real-assembly control: all 14 complete hit records agreed between
  Linux and native Windows Python/BLAST under Wine.
- Full 322,172-pair S. aureus dataset through Windows Python/SKESA 2.4.0 under
  Wine: ST20, exact 7/7 alleles, 20 contigs, 2,758,099 bp and N50 345,017 bp.
  Canonical contig sequences matched the Linux SKESA baseline; originals
  retained their hashes. Full pair validation and execution provenance retained.
- Native Windows SKESA compilation, DLL closure and adapter smoke passed
  [hosted Windows CI](https://github.com/iowa69/WMLSTudio/actions/runs/34710155127).
- Linux and native Windows/Wine surface audits exercised 42 page transitions
  each at 1380×940 and 1080×720. The audit also led to a reference-label cache
  fix: routine progress refresh no longer repeatedly scans every scheme folder.
- Targeted regressions cover stale AMR evidence, failed/stale profile reuse,
  lossless-but-unverified SHA allele tables and SHA-guarded input relinking.
- Packaging gates check app-local PE dependencies, readable Windows system
  fonts, dynamically loaded CDS implementations and complete attribution files.
  Windows uses native Schannel for Qt TLS; Python retains its own SSL libraries.

These are bounded software and scientific controls, not estimates of clinical
sensitivity/specificity, a population-level validation or SeqSphere+ equivalence.
The old 0.1 package and its 198-test record are historical; they do not describe
the current runtime. Final package gates and SHA-256 appear in the
[0.2 release notes](https://github.com/iowa69/WMLSTudio/releases/tag/v0.2.0-workbench.1).

## Portable artifact identity

The release uses the unmodified Windows artifact from successful run
`34713168070`, source/build revision `e5d5f726360de732a766550f06486a18d1cc2edc`.
`WMLSTudio-Windows-x64.zip`: **147,572,214 bytes** compressed; 462,042,768 bytes
of uncompressed files. SHA-256:

```text
695818d74f1021abf6c3691f4c0150748235d8b84529be2639acc864203b0923
```

ZIP integrity and embedded source-archive CRC checks passed. All packaged
sequence files belong to the declared classical-scheme and NCBI starter
reference directories; no local sample inputs or cgMLST.org snapshots are
included. The canonical reference fingerprint agrees across both platforms:
`65293f60de5dc9d55e3a8e20177a05f7220337888e517ab1a0fca6c919e6150d`.

After downloading that exact ZIP, native Windows binaries under Wine passed
QScreen captures at 1380×940 and 1080×720, real ST184 typing (1.855 s), and
HYDRA execution reproducing all 14 complete baseline hits (42.223 s). The input
assembly and all 202 starter-reference files retained their hashes. Independent
inspection resolved 400 dependency links and 44,364 imported symbols across
217 x64 binaries, without looking up build-host redistributables. The ZIP was
not repackaged or modified by these tests.

## Acceptance boundaries

Clean Windows 11/high-DPI/device acceptance is still required. Native hosted
Windows CI and Windows binaries under Wine are separate evidence, not a
substitute for that desktop test. The complete Kleborate/Kaptive,
AMRFinderPlus, agr/SCCmec/spa, MOB-recon and abricate execution stack, independent
species-complex confirmation, contamination assessment and broad clinical
deployment remain unfinished.

See [the capability audit](docs/STUDIO_CAPABILITIES.md),
[microbiology workflow contract](docs/MICROBIOLOGY_WORKFLOWS.md),
[build instructions](docs/STUDIO_WINDOWS.md) and [methods](tasks/METHODS.md).
