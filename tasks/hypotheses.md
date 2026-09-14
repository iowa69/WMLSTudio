# Validation hypotheses

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
