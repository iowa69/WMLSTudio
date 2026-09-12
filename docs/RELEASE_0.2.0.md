# WMLSTudio 0.2 — native research workbench

This revision replaces the earlier preview's disconnected workflows with a
sample-centred native desktop. It remains research software, not a validated
diagnostic device or an established SeqSphere+ equivalent.

## Retest

1. Extract the entire Windows ZIP into a writable folder. Keep the executables
   and `_internal` directory together, then start `WMLSTudio.exe`.
2. Use **Help → Practice project** to test the interface with synthetic samples.
3. Open a **copy** of an older project before importing your own samples: the
   project schema is upgraded and cannot be reopened with the old application.

The previous release remains available for rollback. The repository and release
assets remain private; sign in with an account that has access.

## What changed

- Dark native Qt interface, corrected page repainting, screen-aware startup,
  compact-layout checks, menus, contextual dialogs and pre-run review.
- Per-file or batch organism/workflow assignment: automatic provisional MLST,
  manual organism/scheme, or unknown. Full reviewed read pairing and native
  SKESA assembly retain both original read records and derived-assembly evidence.
- Optional managed copies, organism/ST folders and ST filename suffixes without
  changing originals; same-byte hash-guarded relinking for moved inputs.
- Additional cgMLST/wgMLST snapshots and local ad-hoc schemes without replacing
  classical MLST; exact-first matching and CDS-guarded full-SHA novel identifiers.
- Explicit comparison cohorts, saved-library profile reuse, missing-locus guards,
  metadata/manual colors, editable graph layout, group highlights and exports.
- Actual pinned HYDRA assembly execution with nucleotide/protein/mutation
  references, linked feature tables and primary-aware AMR matrices.
- Selected report cohorts and highlighted groups in PDF/HTML/CSV/TSV/JSON;
  complete sequence-free profile bundles for richer evidence exchange.
- Online reference catalogs and explicit immutable updates. Failed/stale profiles
  cannot silently become current comparison evidence; old AMR results remain
  archived and hashless imported mappings remain explicitly unverified.

## Evidence and boundaries

The binary is the unmodified artifact from
[successful build 34713168070](https://github.com/iowa69/WMLSTudio/actions/runs/34713168070),
revision `e5d5f726360de732a766550f06486a18d1cc2edc`. All **466 tests passed on
both hosted Windows and Linux**, with no failures or skips. Frozen checks cover
typing, HYDRA nucleotide/protein evidence, desktop startup, and native Windows
SKESA; Windows additionally passed app-local DLL, real Segoe UI glyph and bundled
CDS execution checks. ZIP: **147,572,214 bytes**; SHA-256:

```text
695818d74f1021abf6c3691f4c0150748235d8b84529be2639acc864203b0923
```

The exact downloaded ZIP was then tested under Wine: native desktop captures at
1380×940 and 1080×720, real ST184 typing, and all 14 complete HYDRA baseline hits
passed. Its input/reference hashes were unchanged. Independent inspection found
no missing dependency or exported symbol across the 217 packaged x64 binaries.

See [the build/test summary](../SUMMARY.md) and
[the capability audit](STUDIO_CAPABILITIES.md) for exact observed gates and
scientific controls. The release page records its artifact checksum and build
revision. Native hosted Windows tests, Windows binaries under Wine and a clean
Windows 11 desktop acceptance test are different evidence; the last remains
necessary, including high-DPI and device-specific behavior.

The complete Kleborate/Kaptive, AMRFinderPlus, agr/SCCmec/spa, MOB-recon and
abricate execution stack is not included. Provisional MLST lineage is not
independent species confirmation; genotypic AMR evidence is not measured
susceptibility; a close MST edge does not establish transmission.

Non-commercial research is the intended use. **cgMLST.org reference contents
are not bundled.** Provider terms and redistribution rights still apply.
Component notices and corresponding source/build materials accompany the ZIP.
The unsigned application does not require disabling Defender or SmartScreen;
review its source, provenance and checksum before deciding to run it.
