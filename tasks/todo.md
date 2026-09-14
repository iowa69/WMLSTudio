# Active work

## Investigation-first revision: the 100-isolate outbreak scenario

- [ ] Translate the 15 user questions into a source-backed problem-to-result workflow contract.
- [ ] Connect native workspace pages through a persistent investigation map, clear scope and actionable evidence states.
- [ ] Add independent species evidence and linked virulence/plasmid/drug-association interpretation without unsupported phenotype or transmission claims.
- [ ] Save named investigations, all-pair threshold clusters, review groups and immutable incremental snapshots.
- [ ] Make graph selection, rich labels, cluster halos and proximity/cohort reports a continuous workflow.
- [ ] Allocate CPU and RAM automatically across cancellable sample jobs and record actual resource plans.
- [ ] Review filename-based read association and explicitly attach validated FASTQ pairs to existing assemblies without replacing their evidence.
- [ ] Evaluate and implement defensible optional sequence-level refinement and missing-locus second-pass workflows with actual runtime checks.
- [ ] Test representative local Klebsiella, Acinetobacter and Enterococcus inputs; preserve originals and record reproducible commands and limitations.
- [ ] Publish a native in-app workflow guide and problem-to-solution documentation with verified screenshots and capability boundaries.
- [ ] Run regression, real-data, native Windows/package and usability acceptance gates before release handoff.

This revision starts from `2135c3e` on `investigation-v0.3`. Raw local sequence
data remain private and unchanged; only source, tests and sequence-free validation
summaries may be versioned. The existing 0.2 release remains available unchanged.

## Acceptance revision: the twelve requests, with honest status

Status words mean exactly this: **interface** — a user can do it in the running
application; **engine** — implemented and tested as a module, with no interface
surface; **open** — not started. An engine-only capability is not a delivered one.

| # | Request | Status | Where |
| --- | --- | --- | --- |
| 1 | Dedicated, interconnected tabs | interface | Seven keyed tabs, orientation strips, shared focus with visible provenance, per-tab cohorts |
| 2 | Auto-organised folders for inputs and databases | interface | Identify before copy, reviewed proposals, `Genus/species` folders, a six-bucket `_Unresolved` tree, re-filing |
| 3 | Right-click menu in every tab | interface | 29 handlers across all 15 registered views |
| 4 | Dual MST — original beside current | interface | Baseline pointer, snapshot replay beside the current tree, change summary with not-assessed states |
| 5 | Highlight / delete / add samples inline | interface | Archive as the default, restorable removals, Recently removed, focus accent distinct from saved highlights |
| 6 | High-resolution screen support the user can adjust | interface | Settings scales and `--display-scale`; graph text size redraws both trees without moving a node |
| 7 | Simple report: MST picture, resistance, proximity | interface | Reports tab action plus the `one_page` preset in the export path |
| 8 | Thresholds pre-set for common pathogens with references | interface | Catalogue and dialog; coverage audited honestly in `docs/THRESHOLDS.md` |
| 9 | Practice datasets, 10 single-species and 20 mixed-genus | interface | Data menu download with caveats, plus `studio_scripts/fetch_practice_cohort.py` |
| 10 | Manual genus/species override | interface | Assign from the context menu, re-file on change, CSV assignment import |
| 11 | Organism-specific tools (SCCmec, Kleborate-style) | interface | Plan-dialog selection, results column, drill-down, report section, CLI flags |
| 12 | Friendly to non-bioinformaticians | interface | Orientation strips, plain-language guide, evidence wording, the simple summary |

### Carried into the next batch

- [ ] Re-stage the bundled characterization starter to manifest format 2 so the
      organism modules stop reporting `not_run` on the shipped snapshot. CI stages
      it fresh, so this is a developer-tree gap; check `frozen-check.json` to see
      which state a given build actually shipped.
- [x] Regenerate `uv.lock` for the new `pyyaml` development dependency, with
      isolated uv 0.10.11 (archive SHA-256 verified against its published sum).
      `uv lock --check` resolves; only pyyaml was added.
- [x] `WorkbenchMixin.set_graph_text_scale` and the `widgets.py` scale it calls.
      Lettering only: no node moves and no edge length changes, so a rescaled
      tree is the same tree.
- [x] `organism_typing` as a flat export field, gated on `current_characterization`
      so typing from an earlier assembly leaves the column empty.
- [x] Delete the dead `navigate` override and the temporary shim in `app.py`.
- [x] Settle the organism-module registry order. It was decided by whichever
      assay module was imported first, so a table and its export could reorder
      their columns between runs.
- [ ] Widen the threshold catalogue, or record why each remaining organism cannot
      be curated. Twelve listed organisms still have no cutoff at all, including
      *S. pneumoniae*, *S. enterica* and *S. capitis*.
- [ ] Guard the per-window `setStyleSheet` re-polish that makes repeated window
      construction grow superlinearly in one process.
- [ ] Run the frozen SCCmec self-comparison for real: it needs a format 2 bundled
      panel, which only a fresh CI build currently produces.
- [ ] Clean Windows 11 acceptance, now including high-DPI and per-monitor scaling.

## Workbench revision requested after user acceptance feedback

- [x] Replace page-opacity rendering and implement a consistent dark native theme.
- [x] Add organism-aware import, batch assignment and pre-run review.
- [x] Add managed input storage and non-destructive ST organisation.
- [x] Automate conservative scheme discovery and retain organism evidence.
- [x] Add online organism/scheme catalog and versioned reference downloads.
- [x] Add explicit comparison cohort and secondary cgMLST/wgMLST analysis selection.
- [x] Add graph styling, persisted layout and sample/cluster highlights.
- [x] Join HYDRA AMR metadata to stable sample identities and feature tables.
- [x] Export selected report cohorts and highlighted sample groups.
- [x] Execute native Windows SKESA and HYDRA, with real-data Linux/Wine agreement.
- [x] Guard stale AMR, failed profile reuse and unverified external allele tokens.
- [x] Validate native workflows, rebuild Windows package and push the source.
- [x] Publish and download-verify the private Windows retest release.

Published `v0.2.0-workbench.1` with the unchanged CI Windows ZIP and checksum.
A fresh authenticated release download passed SHA-256 and full ZIP CRC checks;
the repository and release remain private. Artifact SHA-256:
`695818d74f1021abf6c3691f4c0150748235d8b84529be2639acc864203b0923`.

Design and acceptance contract: [WORKBENCH_DESIGN](../docs/WORKBENCH_DESIGN.md).

## Earlier implementation (not a production acceptance claim)

- [x] Native Qt desktop, accessible navigation and background jobs.
- [x] Streaming sequence validation and QC; exact MLST/cgMLST calling.
- [x] Local scheme import, portable resource resolution and bundled demo.
- [x] Persistent projects, provenance, trustworthy distance/MST and exports.
- [x] HYDRA JSON adapter, native evidence views and persistence.
- [x] Scientific edge-case tests, GUI tests and real-data validation.
- [x] Windows portable packaging recipe and CI; truthful capability matrix.
- [x] Execute Windows binary under Wine, including native Windows Qt platform and real typing.
- [ ] Windows 11 clean-machine acceptance, including USB/non-ASCII paths and cancellation.
- [x] Approximate/novel allele engine and per-target QC; bounded real cgMLST control.
- [x] Native Windows raw-read assembly and HYDRA engine execution.
- [ ] Independently labelled, multi-species cgMLST sensitivity/specificity benchmarks.
- [ ] Reproduce broad comparator benchmarks before competitiveness/production claims.
