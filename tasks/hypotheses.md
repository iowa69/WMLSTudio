# Validation hypotheses

## Acceptance revision: alternatives introduced by this batch

| Explanation | Status | Discriminating check |
| --- | --- | --- |
| A frozen build silently omits an assay module the registry imports by computed name, so a portable install fails on every characterization run | **Addressed** | Disassembly shows `_load` emits no `IMPORT_NAME`; an AST scan finds no static reference; the spec derives hidden imports from `runner.__module__`; `check_frozen` asserts the registry's evidence blocks are present in a frozen report |
| A missing or out-of-date reference panel section is read as a negative assay result | Addressed by construction | A format 1 snapshot yields `not_run` naming the missing section; `check_frozen` reports `staged`/`absent` with the reason rather than a silent pass |
| An automatically created organism folder is taken for a laboratory identification | Unresolved | Folder text, quarantine README and cohort caveats all state the distinction; needs observation with a real non-specialist user |
| Filing quarantines the wrong isolates because "nothing matched" is read as "not this organism" | Unresolved | Quarantine reason is derived from the engine's own numbers, not from its prose; `Not_in_reference_panel` is an explicit statement about the installed panel |
| A published threshold is applied to a scheme it was never measured on | Addressed | `record_decision` refuses without an exact scheme key, the full target count, the reference fingerprint, the caller, the missing-data policy and a written justification; seventeen catalogue entries are refusable by that rule today |
| A practice cohort is mistaken for a validation set, or its agreement with a published study is read as software validation | Addressed in text, unresolved in use | Five caveats ship inside each cohort manifest and are shown at download; no expected ST, cluster or threshold is distributed |
| A re-rendered baseline snapshot reads as 100% locus completeness because only called alleles were stored | Addressed | `snapshot_graph` documents that `len(alleles)` is the callable count; completeness must come from the stored `callable_loci`/`total_loci` scalars |
| An incomparable pair of snapshots is summarised as "no change" | Addressed | Every not-assessed section and count is `None`, never `[]` or `0`; consumers must branch on `None` |
| A whole-interface scale leaves the window larger than the screen and unrecoverable | Addressed | Only scales the screen can display the whole window at are offered; `--display-scale` is a documented recovery path |
| Sequence data reach the repository or the portable ZIP through a practice cohort or the species panel | Addressed | Gitignore covers the download directories; the build refuses to package either, and `check_frozen` re-checks the built bundle |
| An SCCmec type call is read as an MRSA determination or a methicillin susceptibility result | Unresolved | `official_type` is unconditionally `None`, `mecC` is always `not_assayed`, and the limitation is repeated in every drill-down and report section; needs user observation |

### Recorded deferral: whole K and O locus typing

Not implemented, and recorded rather than promised. Kaptive v2.0.9 (commit
`b3856eac6e76b3017aa993319da2a8ea967a1ba0`, GPL-3.0-or-later) ships
`Klebsiella_k_locus_primary_reference.gbk` (8,325,855 B),
`Klebsiella_o_locus_primary_reference.gbk` (325,387 B),
`Klebsiella_k_locus_variant_reference.gbk` (1,303,472 B) and
`Klebsiella_o_locus_primary_reference.logic` (591 B). Kaptive 3 master no longer
carries the databases in-repo, so v2.0.9 is the pinnable artefact.

Two blockers. **Format**: GenBank, with no runtime parser and no justification
for adding one; it needs a build-time converter to whole-locus FASTA, per-CDS
gene FASTA and a locus-to-gene JSON, with Biopython in the development group only
and the manifest recording upstream `.gbk` SHA-256, derived SHA-256 and converter
version. **Algorithm**: Kaptive's output is a locus assignment with a confidence
grade derived from expected-gene coverage, missing and extra genes and locus
contiguity, plus `.logic` O-locus special rules. Reimplementing that and calling
the result a K or O locus type would be precisely the overstatement this project
refuses.

If pursued, the honest shape is `capsule_locus_candidate` by interval-union
coverage with per-expected-gene presence and a contig count, status capped at
`provisional_reference_match`, `official_locus` unconditionally `None`, and a
limitation naming what is not implemented. Bundle cost about **+9.9 MB**.

## Investigation revision: active scientific and usability alternatives

| Explanation | Status | Discriminating check |
| --- | --- | --- |
| A presumptive K. pneumoniae collection includes other species-complex members | Unresolved | Independent genome-reference evidence, runner-up margin and panel coverage; MLST is not independent taxonomy |
| Apparent proximity is a missing-locus, reference-version or mixed-sample artifact | Unresolved | Shared-locus denominators, quality gates, identical scheme snapshots and ambiguous-call exclusions |
| A threshold cluster is a single-link chain rather than uniformly close isolates | Unresolved | All comparable pair edges, cluster diameter, nearest neighbours and explicit chaining warnings |
| Shared resistance/plasmid markers reflect common mobile elements rather than recent isolate transmission | Unresolved | Keep core relatedness, accessory content and contig co-location separate; require epidemiology/long-read confirmation for stronger claims |
| Genotype-to-drug association is overinterpreted as a measured susceptible/resistant phenotype | Unresolved | Versioned determinant annotations; AST metadata kept separate; unknown/not-tested are never susceptible |
| New sampling or database updates alter apparent cluster membership without biological change | Unresolved | Immutable investigation snapshots, stable sample/profile keys and merge/split/addition audit |
| File stems pair the wrong reads or duplicate lanes | Unresolved | Conservative candidate matching, withheld ambiguity, explicit confirmation, full pair validation and input hashes |
| Parallel jobs oversubscribe RAM/threads or stall the interface | Unresolved | Resource-admission controls, actual child thread budgets, cancellation, GUI heartbeat and representative workload benchmarks |

Profile: computational validation and exploratory comparative genomics, not a
clinical intervention trial. Results must not establish clinical validity from
regression agreement alone. Revisit these alternatives after each real-data gate.

## Workbench revision hypotheses

| Explanation | Evidence needed | Status |
| --- | --- | --- |
| Full-page opacity effects retain stale native child frames | Effects removed; Linux and Windows/Wine native-surface navigation at two sizes | Addressed; clean Windows 11 acceptance remains |
| Scheme similarity is mistaken for species confirmation | Unique, weak, complex-level and competing-scheme controls | Tested; labels remain provisional |
| Sample selection is lost on refresh or sorting | Stable-ID cohort tests across page changes and worker completion | Tested |
| Managed copies overwrite or detach originals | Collision/cancellation/reopen and same-byte relink tests | Tested; originals unchanged in real read assembly |
| Missing/partial references create misleading distances | Same-snapshot selection and explicit shared-locus denominator tests | Tested |
| Name-based AMR linking merges unrelated samples | Duplicate-name, explicit mapping and input-hash evidence-state tests | Tested; hashless mappings remain unverified |
| Public API updates silently omit restricted reference records | Live catalog response audit and preserved restriction notices | Restricted access disclosed; no completeness claim |
| Changed inputs leave old AMR evidence looking current | Re-analysis with a different input hash, historical evidence retention and report regression | Tested |
| Failed results or external hash tokens become validated alleles on reuse | Failed/stale library rejection and lossless-but-unverified table round trips | Tested |

## Previous revision

These are software validation alternatives, not biological study hypotheses.

| Hypothesis | Test evidence required | Status |
| --- | --- | --- |
| Biological variation changes an allele | Practice ST1→ST2 differs at one locus; exact-caller and GUI tests | Tested |
| Parsing or reverse-strand handling loses loci | Wrapped, compressed, reverse-strand and chunk-boundary tests agree | Tested |
| Null: identical profiles have zero distance | Two identical practice profiles yield 0/7 differences | Tested |
| Sampling overstates read QC | Read-prefix limit tests and real FASTQ report disclose sampled=true | Tested |
| Database mismatch creates false type calls | Changed fingerprint/mixed-scheme comparison tests reject pairing | Tested |
| Missing evidence masquerades as similarity | Practice partial profile is disconnected at 0.95 shared coverage | Tested |

Reflection: real reference data contained incomplete profiles and ambiguous allele
sequences. The loader reports exclusions explicitly; it never matches ambiguity
codes as nucleotide wildcards. All 162 staged schemes load with those disclosures.
Agreement with the older caller on four assemblies is useful regression evidence,
not a sensitivity/specificity study or cgMLST equivalence claim.
