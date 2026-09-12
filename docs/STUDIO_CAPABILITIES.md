# WMLSTudio capability and acceptance audit

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
