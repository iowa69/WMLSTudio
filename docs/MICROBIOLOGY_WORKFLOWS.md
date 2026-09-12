# What WMLSTudio should deliver to a microbiology laboratory

This is the scientific product contract behind the native workbench, informed by
an audit of [MLSTudio](https://github.com/iowa69/mlstudio), its local 2.0 checkout,
and actual calling/storage code. It is not a claim that every item below is
clinically validated. See the [implementation audit](STUDIO_CAPABILITIES.md).

## The result a user should take home

A defined cohort with defensible isolate identities, quality limitations,
classical STs, compatible genome-wide profiles, linked resistance/virulence
evidence, editable epidemiological metadata, an interpretable comparison, and
a self-contained report that another analyst can audit and reproduce.

The deliverable is not merely a colored tree. Every node needs an answer to:
which isolate, which input bytes, which reference snapshot, which method and
parameters, which loci failed, and which evidence supports the interpretation?

## Everyday workflows

1. **Routine intake.** Import files and review read pairs; assign a known organism,
   request provisional MLST discovery, or retain unknown. Store protected copies
   if requested. Keep sample identity separate from its filename.
2. **Type and characterize.** Review QC, call classical MLST, then apply the
   appropriate additional scheme. Run requested AMR methods against an explicit
   database snapshot. Preserve conflicting/inconclusive evidence for review.
3. **Investigate a cluster.** Search the saved library by organism, ST, gene,
   metadata or collection. Select the intended cohort and scheme. Inspect missing
   loci and comparison denominators before applying a locally justified group
   threshold. Keep epidemiological annotations distinct from genetic evidence.
4. **Communicate.** Export selected isolates and highlight those under
   investigation. Include QC exceptions, method/database versions, allele
   denominators and genotypic AMR evidence. Do not silently infer phenotype or
   transmission.
5. **Revisit.** Reuse frozen profiles from earlier projects or a collaborator's
   bundle. Retain old analyses when a new scheme/database snapshot is used.
   An update should never silently rewrite the original report.

## Required scientific controls

| Question | Required behavior |
| --- | --- |
| Is the identity supported? | Distinguish user assignment, provisional MLST lineage and independent species confirmation; do not force a species within an unresolved complex |
| Are the data suitable? | Report pair validation, sequence metrics, ambiguous bases and target completeness; successful execution is not a purity/completeness pass |
| Is the locus trustworthy? | Separate exact, validated local novel, missing, partial, duplicate/paralog, mixed and ambiguous evidence |
| Are profiles comparable? | Require reference identity and report shared loci/missingness; never convert missing tokens to matches or novel alleles |
| What does AMR mean? | Preserve identity, coverage, method, database, primary/secondary status and mutation evidence; absent results remain unknown |
| Does proximity prove transmission? | No. A spanning tree visualizes distances; sampling, dates, epidemiology and organism-specific validation remain necessary |
| Can this be reproduced? | Preserve hashes, software/tool versions, reference manifests, parameters, history and original inputs |
| Can this be shared? | Confirm sample-data permissions and database/software rights; make the export cohort explicit |

## Decisions implemented in this revision

- Primary MLST and additional cgMLST snapshots are separate.
- Novel sequences use full hashes, not collision-prone short suffixes.
- Ambiguous organism/scheme evidence stays unresolved.
- HYDRA is an actual optional engine, with linked provenance and primary-aware counts.
- Read assembly retains both read records and its derived assembly.
- Graph presentation is stored independently of scientific distances.
- Library/profile-bundle reuse works without FASTA availability.
- ST organization affects managed copies; reference updates create new snapshots.

## Work that cannot be replaced by interface polish

The complete organism-module stack, validated species-complex discrimination,
contamination assessment, platform-specific limits, formal AMR interpretation
rules and multi-laboratory concordance testing remain necessary before calling
this a production competitor. They require bounded reference datasets,
independent truth, reproducible methods and documented discrepancy review—not
invented success states or a checklist of buttons.
