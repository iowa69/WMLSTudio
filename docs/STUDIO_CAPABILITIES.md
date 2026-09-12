# WMLSTudio 0.2 capability and acceptance audit

Audit date: 2026-09-12. The current workbench implements the following workflows.
This is not a production, clinical, or SeqSphere+ equivalence claim. See the
[scientific product contract](MICROBIOLOGY_WORKFLOWS.md).

| Workflow | Current implementation | Boundary |
| --- | --- | --- |
| Native desktop | Dark Qt widgets, menus, contextual dialogs, animated helix, direct viewport painting | No WebUI/browser server or full-page opacity cache |
| Import / storage | Per-file/batch automatic, manual or unknown labels; hashed copies; optional organism/ST folders and ST suffixes | Original files never renamed or deleted |
| Automatic MLST | Installed classical-panel search with full exact verification and tied/weak evidence handling | Provisional lineage, not independent species/purity confirmation |
| Classical typing | Both-strand exact alleles and complete-profile lookup | No registered ST for missing/mixed loci or an unregistered combination |
| cgMLST/wgMLST | Exact-first calling, native BLAST+, complete-CDS checks and full-SHA novel IDs | Local novel IDs are not centrally registered allele numbers |
| Local ad-hoc scheme | Reference-anchored unique complete CDSs; explicit prevalence/identity/coverage | Cohort-defined research nomenclature |
| Paired short reads | Reviewed pairing, complete pair/QC validation, cancellable native SKESA, atomic assembly/provenance | Native tool must pass its platform gate; no fastp/SPAdes/long-read pipeline |
| Unassembled reads | Labelled bounded-prefix Phred+33 QC and full file-byte hash | Not complete-file quality or contamination validation |
| Library | SQLite projects, collections/history, organism/ST/AMR and metadata research search; cross-project frozen-profile reuse | Index represents saved evidence, not automatic reanalysis |
| Interoperability | Profile-only tables with missing-token validation; complete profile export and sequence-free bundles | Unverified external tables remain distinct snapshots |
| Compare | Explicit cohort/scheme, shared-locus denominators, missing-data guards, minimum spanning forest | Same fingerprint required; 95% default overlap is not a universally validated cutoff |
| Graph | Metadata/manual colors, editable labels/layout, identical-genotype pies, halos and persisted styles | Styling does not change distances; MST Newick is not a phylogenetic tree |
| HYDRA | Pinned upstream assembly engine in a dedicated native process, BLAST+, NCBI starter | Not the AMRFinderPlus executable; direct-read/pileup and full lineage stack disabled |
| Linked features | Stable-ID mapping, primary-aware AMR matrix, QC/typing/annotations | No report means unknown; gene evidence is not measured susceptibility |
| Reports | Explicit cohort/highlights, PDF/HTML/CSV/TSV/JSON, additional-profile summaries and complete JSON evidence | User-defined groups do not establish transmission |
| References | PubMLST/Pasteur catalogs, rights-gated cgMLST.org access, local import and immutable updates | Authentication, submission-date restrictions and provider rights limit availability |

## Current observed validation

- S. epidermidis positive-control assembly: ST184 and expected provisional lineage.
- Three P. aeruginosa assemblies: ST155, ST2952 and ST1858; all 21 alleles agreed
  with WMLST 1.2.1 / BLAST 2.12 using the same reference snapshot.
- E. faecium GCA_043869095.1 with a 1,423-locus research snapshot: 1,410 exact,
  one validated novel, seven missing, four ambiguous and one mixed locus. The
  mixed result is conservatively excluded from comparison. The 1,614-nt novel
  EFAU004_00165 CDS was independently checked against raw assembly coordinates.
- HYDRA: all 14 complete hit records from one real S. epidermidis assembly were
  identical on Linux and native Windows Python/BLAST under Wine.
- Ad-hoc truth test: two assemblies differing at one internal codon yielded one
  retained complete locus with distinct full-SHA alleles and exact recall.
- Native SKESA 2.4.0 passed its [hosted Windows build and validation
  gate](https://github.com/iowa69/WMLSTudio/actions/runs/34710155127): a 12,000-bp
  synthetic reference yielded a matching 11,928-bp contig. The real Python
  adapter also passed with Unicode/comma paths, and the executable launched
  without MSYS2 on PATH.
- Real S. aureus SRR12343864: Windows SKESA 2.4.0 and Windows Python under
  Wine 11.17 completed all 322,172 read pairs, producing 20 contigs,
  2,758,099 bp and N50 345,017 bp. Windows-native MLST calling returned ST20
  with all seven loci matched. The strand-normalized contig sequences and
  seven-locus profile were identical to the Linux SKESA 2.5.1 baseline;
  original compressed-input SHA-256 hashes remained unchanged. The Linux
  launcher was externally interrupted; its completed child output was
  independently recovered and validated. The Windows adapter completed normally.

These are bounded controls, not independent population-level sensitivity or
specificity estimates. Local sequence-bearing artifacts remain ignored. The
release notes record final observed build/test gates; a CI recipe alone is not
evidence that a run passed. Clean Windows 11/high-DPI/device acceptance and
independently labelled multi-species cohorts remain necessary.

The hosted SKESA gate validates that component, not the final complete desktop
release. Its success and the separate real-data Wine run do not constitute a
clean Windows 11 acceptance test.

## Remaining full-suite work and reference rights

The complete Kleborate/Kaptive, AMRFinderPlus, agr/SCCmec/spa, MOB-recon and
abricate execution stack is not included. Species-complex resolution,
contamination quantification, long-read assembly, fastp preprocessing, validated
phenotype flags and multi-user clinical deployment remain separate work.

cgMLST.org data are **not bundled**. Its
[server policy](https://www.cgmlst.org/serverpolicy.html) restricts use and requires
permission for database-driven products/services. PubMLST has separate
[submission-date-dependent terms](https://pubmlst.org/terms-conditions).
Software licenses do not grant reference database rights.

---

## Historical 0.1 audit — superseded implementation scope

The remainder records the earlier release for traceability. Statements below
about absent execution/assembly/novel calling describe **0.1 only**, not 0.2.

Audit date: 2026-09-12. This native desktop implementation is an initial,
testable replacement foundation. It does not yet reproduce every MLSTudio
workflow or establish scientific equivalence to Ridom SeqSphere+.

## Implemented scope

| Workflow | Current behavior | Meaning of the result |
| --- | --- | --- |
| Native desktop | Qt widgets, local project storage, guided input and analysis | No local web server or embedded browser engine |
| Assembly input | Streaming FASTA, including gzip/bzip2; content hashes and assembly QC | Original sequence files remain unchanged |
| MLST | Exact nucleotide allele matching on both strands against local schemes | ST is assigned only from complete unambiguous profiles present in that snapshot |
| Imported cgMLST | Same exact engine applied to user-supplied locus FASTA directories | Exact known alleles only; no novel-allele inference or nomenclature registration |
| Bundled references | 162 cached PubMLST MLST schemes, staged explicitly during development/build | Static snapshot, not every available cgMLST scheme |
| FASTQ QC | Bounded prefix of records with Phred+33 assumption, full input-byte hash | No read assembly, direct read MLST, contamination inference, or complete-file QC claim |
| Allele comparison | Explicit missing-data handling and minimum spanning tree | Comparisons require compatible scheme content; no universal outbreak threshold |
| HYDRA reports | Validated JSON import of per-sample AMR, virulence, mutation and typing evidence | External results retain original methods and provenance; no HYDRA execution bundled |
| Portable build | One-folder desktop and CLI, Python/Qt/compiled engine included | Native Windows output requires a Windows build and observed Windows validation |

## Gaps against the requested full suite

[MLSTudio's documented workflows](https://github.com/iowa69/mlstudio) also include
read assembly, AMRFinderPlus, Kleborate, organism-specific modules and broad
database connectors. Native execution of those workflows is not included in
this version. HYDRA report import covers exchange and inspection of already
computed evidence; its execution dependencies need separate integration work.

Ridom documents automated assembly pipelines and per-target quality checks,
including coverage and frame-shift criteria. Its target states distinguish
quality failure from failure to locate a target. Exact matching here does not
recover approximate targets or diagnose the biological reason for a missing
match. Competitive parity requires validated approximate/novel allele calling,
target QC, contamination assessment, organism-specific interpretation and stable
nomenclature integration. See the [Ridom pipeline guide](https://www.ridom.de/seqsphere/ug/v90/Tutorial_for_SeqSphere%2B_Pipeline.html).

An MST is a compact representation of allele differences, not a reconstructed
transmission chain. Low shared-locus coverage can hide differences. The
application must preserve missing/ambiguous calls and report comparison coverage.
Distance agreement across tools is meaningful only when scheme versions, locus
sets, allele definitions and missing-data policies agree.

## Acceptance evidence

`studio_tests` exercises sequence parsing and scientific edge cases, project
storage, comparisons, exports, native widgets, database staging and HYDRA import.
`studio_scripts/validate_real_data.py` provides opt-in local validation with input
and code SHA-256, scheme digest, measured timing, sampled-QC disclosure,
checkpointed JSON, and optional expected-ST assertions. Resume requires matching
input bytes, code, scheme, limits and expectation; failed records are retried.

The local positive-control assembly is
`../wmlst/tests/data/example.fna`, expected *Staphylococcus epidermidis* ST184
using `../wmlst/db/pubmlst/sepidermidis`. Example invocation:

```bash
uv run python studio_scripts/validate_real_data.py \
  --assembly ../wmlst/tests/data/example.fna \
  --scheme src/wmlstudio/resources/schemes/sepidermidis \
  --expected-st 184 --output artifacts/sepidermidis-validation.json
```

For a local cohort, supply `--assembly-root PATH --scheme PATH --max-assemblies N`.
For sampled read QC, supply `--fastq PATH --max-reads 10000`. Inputs are read only;
only the specified result report is written. A cohort without independent labels
is a robustness/performance exercise, not an accuracy benchmark.

Windows CI is configured; its existence is not evidence that a Windows run has
passed. Public release requires observed clean-machine tests and reference-data
redistribution review. Production claims require independently labelled,
multi-species cohorts, agreement analysis against validated comparators, and
explicit investigation of every discrepant and failed call.

## Observed local checks, 2026-09-12

| Real input | Observed result | Scope of evidence |
| --- | --- | --- |
| S. epidermidis `example.fna` | ST184; about 0.50 s | Matched the declared positive-control ST; frozen CLI also returned ST184 |
| P. aeruginosa `GCF059737285v1` SPAdes assembly | ST155; about 1.75 s | Matched WMLST 1.2.1 with BLAST 2.12.0 |
| P. aeruginosa `GCF059544455v1` SPAdes assembly | ST2952; about 1.33 s | Matched the same comparator |
| P. aeruginosa `GCF049679795v1` SPAdes assembly | ST1858; about 1.38 s | Matched the same comparator |
| S. aureus `SRR12343864_1.fastq.gz` | 10,000 reads; Q30 98.924%; about 1.46 s | Prefix QC only; Phred+33 assumption, complete compressed-file SHA-256 |

The three P. aeruginosa STs and all 21 individual allele assignments agreed with
the comparator using the same local reference snapshot. These are three real
assemblies, not an independent population-level accuracy benchmark. Their JSON
reports, original input hashes and comparator outputs are under `artifacts/`.
Timing includes validation overhead and is specific to this Linux host.
