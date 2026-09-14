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
