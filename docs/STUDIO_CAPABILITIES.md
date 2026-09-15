# WMLSTudio 0.4 pipeline revision — capability and reach

Audit date: 2026-09-15. This revision adds read trimming and assembly metrics at
the front of the pipeline, a SNP-distance tree beside the two allele trees, a
contig-level plasmid evidence screen, and an AMR reference catalogue that names
every database the engine can use instead of silently assuming two.

The **Reach** column says where each capability can actually be used, because a
capability that exists only in `studio_tests` is not one a microbiologist has.

| Capability | What is implemented | Boundary | Reach |
| --- | --- | --- | --- |
| Read QC — read trimming | Bundled fastp 1.3.7, pinned to its upstream SHA-256, run on a reviewed pair; fastp's own JSON report is displayed, never recomputed; the trimmed pair is recorded beside the originals, which are never replaced | Adapter trimming and quality filtering. Not an isolate validation, a purity check, a species assignment or a clinical result. Upstream publishes no Windows binary, so a Windows package carries fastp only when a reviewed native artifact was staged, and otherwise says trimming is unavailable and leaves the reads usable untrimmed | Where staged |
| Assembly metrics and read source | Contigs, total length, N50, largest/smallest contig, GC with its denominator, N and ambiguous bases, and the assembler's own reported depth; every assembly records whether it was built from trimmed or original reads | Assembler-reported depth is the assembler's number, not an independent coverage measurement; metrics are not a quality verdict | Interface |
| SNP tree — split k-mer distances | SKA2 split k-mer distances over the same cohort, with shared split k-mers as the denominator on every pair, an adjustable comparability floor and its basis, refused pairs listed with their reason, and a separate tree that names its own quantity | A separate quantity from MLST and cgMLST allele distances: no shared scale, axis, column or threshold. Insufficient overlap is a refused pair, never a zero distance. No curated SNP cutoff is applied, and a tree is not a phylogeny | Engine and payload; the tab is the last wiring step |
| Plasmid evidence | Replicon markers placed on the contigs they sit on, determinants reported as co-located or not, assembler-declared closure and coverage departure from the chromosomal backbone, and cohort replicon co-occurrence | **Not MOB-suite and not equivalent to it**: no relaxase or MPF type, no oriT, no plasmid reconstruction, no mobility prediction, no plasmid count. Contig co-location is an assembly observation, not proof two genes travel together | Interface |
| AMR reference catalogue | Every reference set the pinned engine can use is listed by name with its provider, licence, citation, upstream address and installed state; the NCBI core is bundled, anything else is a per-set download | A set with a non-open or unrecorded licence is listed and downloadable but never swept into "install and update everything"; provider terms are shown before any download | Interface |
| Bundled offline core | AMRFinderPlus nucleotide and protein reference data **plus** point mutations — 13 DNA catalogues, 30 organisms with curated protein mutations, 31 in total — bundled so the first run works with no network | The engine accepts 32 organism names, so one accepted organism has no catalogue at all. An organism outside the list is screened for acquired genes only: an absent catalogue, stated as such, never a negative mutation result | Interface and command line |
| Size budget | The portable archive is measured per component and the build fails above 1,000,000,000 bytes, warning from 80% | A size gate is a distribution promise, not a scientific control | Build gate |

MLST, cgMLST and SNP distances are three different quantities. They are reported
on three separate scales, in separate columns, with separate denominators, and no
threshold from one is offered for another.

## What the portable package costs

Measured from the staged trees this revision builds, compressed as the archive
stores them:

| Component | Bundled or downloaded | Compressed in the archive |
| --- | --- | --- |
| FastQC and its private Java runtime | bundled | 173.2 MB |
| BLAST+ 2.17.0 | bundled | 41.4 MB |
| SKA2 0.5.1 | bundled | 40.2 MB |
| Species, virulence and organism-module panels | bundled | 25.2 MB |
| MLST and cgMLST scheme snapshot | bundled | 7.0 MB |
| AMRFinderPlus core **and point mutations** | bundled | 5.8 MB |
| **fastp 1.3.7** | bundled where staged | **4,981,685 B** |
| SKESA | bundled | 3.3 MB |
| Read trimming, assembly metrics and read source (Python) | bundled | 25,011 B |
| SNP tree, including the SKA2 adapter's growth (Python) | bundled | 19,787 B — the tool itself was already staged |
| Database catalogue and install-everything (Python) | bundled | 22,561 B |
| Plasmid evidence screen (Python) | bundled | 10,128 B — no database, no binary, no new licence |
| Every other AMR database the catalogue lists | downloaded on request | 0 |
| MOB-suite database | not available | 0 — its 473 MB alone would have taken the archive past 845 MB |
| Practice cohorts and the broad species panel | downloaded on request | 0 |

The Python figures are each module compiled at optimize level 2 and compressed as
PyInstaller stores it, measured against the previous revision rather than
estimated: **77,487 bytes for everything this revision adds in code**. fastp is
the only new binary, at 4,981,685 bytes compressed — 0.5% of the budget — and it
adds nothing at all to a Windows package today, because upstream publishes no
Windows binary and the package ships without it.

The archive gate refuses any build over 1,000,000,000 bytes and warns from
800,000,000, with the table above printed so the growth has a name; the same
measurement runs on the frozen folder before the archive exists, as a
deliberately conservative prediction.

## 0.3 investigation revision

| Capability | What is implemented | Boundary | Reach |
| --- | --- | --- | --- |
| Interconnected tabs | Seven keyed tabs, each with a purpose sentence, a next-step button and a per-tab help entry | Tab navigation is not a workflow guarantee; a tab can be opened out of order | Interface |
| Shared focus and per-tab cohorts | Selection anywhere sets a visible focus with its origin and time; a tab adopts it only on an explicit button press, and states the cohort's provenance | Focus never writes a cohort; two tabs may legitimately hold different cohorts | Interface |
| Display scaling | Text size 80–150% applied live; whole-interface scale applied before the application starts, offering only sizes the screen can show the whole window at; `--display-scale` recovery | High-DPI behaviour is not verified on a clean Windows 11 host. The graph text size is saved but not yet applied to the drawn graph | Interface, except the graph text scale |
| Right-click actions | Twenty-nine handlers across fifteen views, with counts shown, single-selection actions disabled with a reason, and project-writing actions disabled while a job runs | Menu entries appear only where their handler exists | Interface |
| Archive instead of delete | Archiving retains the row, results, analyses, history and the user's file; removal writes the full record and every analysis into history first and can be restored | An archived isolate is hidden from working views, not deleted, and is still carried by explicit project exports | Interface |
| Automatic organism filing | Layered identification (broad ANI panel, focused complex panel, MLST panel compatibility) run before any copy is written, a reviewed proposal per file, `Genus/species` folders, a six-bucket `_Unresolved` tree, re-filing, manual override and CSV assignment import | A folder is a filing decision, not a laboratory identification; nothing below genus is ever auto-confirmed; originals are never moved | Interface |
| Organism-specific modules | Registry with a three-state match rule; SCCmec typing, *Klebsiella* locus STs and *wzi*/*wzc* capsule markers; multi-source pinned staging with manifest format 2; plan-dialog selection, results column, drill-down, report section and CLI flags | Not equivalent to Kleborate, Kaptive, staphopia-sccmec or SCCmecFinder; mecC not assayed; no K locus inferred; a format 1 snapshot reports `not_run`, never a negative | Interface and command line |
| Dual-tree change comparison | Baseline snapshot pointer, exact replay of a stored snapshot beside the current one, and a diff carrying distance, denominator, edge and cluster changes | Incomparable snapshots report "not assessed", never "no change"; an MST is not a phylogeny | Interface |
| Simple summary report | Five-section plain-language layout with an embedded tree image, resistance genes, closest matches and a fixed limitations block | Two to four pages, not one; the susceptibility caveat cannot be switched off | Interface |
| Practice cohorts and species panel | Two pinned public cohorts (10 single-species, 20 mixed-genus) and a 17-taxon species panel, downloaded and checksum-verified on request | No sequence enters the repository or the ZIP; no expected answer ships with either cohort | Interface and command line |
| Threshold catalogue | 29 organisms listed, 28 entries, 7 cited sources; 8 organisms with a scheme-bound cutoff | No cutoff is ever auto-applied; 12 listed organisms have no curated cutoff; no bundled scheme can bind a cgMLST cutoff | Interface and [audit](THRESHOLDS.md) |
| Flat exports | `organism_typing` as a single per-isolate field beside the AMR fields | Not implemented yet: the field is specified and not present in the CSV/TSV/JSON column set | Not yet reachable |

Nothing in this table is a clinical validity claim, a benchmark result or a clean
Windows 11 acceptance claim. The frozen-build checks described below are
self-consistency controls against bundled references, not independent biological
validation.

## Frozen-build verification

`studio_packaging/check_frozen.py` exercises the built bundle rather than assuming
staging succeeded. Added in 0.4:

- The bundled fastp is re-verified against its manifest and then actually run
  inside the frozen bundle on four synthetic pairs with Unicode paths. Its own
  report is read, the original read files are hashed before and after, and the
  executable is proved to resolve inside the package. A package with no fastp
  reports `not_bundled` **with the reason the application itself shows**, which is
  a legitimate build state and not a silent pass.
- The package is measured per component and gated against the size budget, so a
  bundle that outgrew the download promise fails the build instead of shipping.
- Packaging refuses an AMR core store without its point-mutation catalogues, so
  the offline first run cannot quietly become genes-only.

Added in 0.3:

- The bundled characterization snapshot is re-hashed file by file in the frozen
  bundle and its manifest fingerprint recomputed, so a panel truncated by
  packaging fails the build instead of later reading as a negative assay result.
- The report records the panel's format version, its per-source pinned revisions
  and licences, and whether the organism-module sections are `staged` or `absent`
  **with the reason**.
- The organism-module registry is proved to have loaded inside the frozen process.
  Its assay modules are reached through a computed import that the packaging
  module scan cannot follow, so they are declared explicitly; without that, every
  characterization run in the portable build would fail to import.
- When the panel and the command line both support it, a SCCmec self-comparison
  runs against the panel's own type IVa cassette reference. It is labelled
  *bundled reference self-comparison; not an independent biological validation*,
  exactly as the species control is.
- Practice-cohort genomes and the broad species panel are rejected at packaging
  time and their absence is re-checked in the built bundle. Both are the user's
  own downloads and never travel in the ZIP.

See [organism modules](ORGANISM_MODULES.md), [thresholds](THRESHOLDS.md) and
[practice cohorts](TEST_DATASETS.md).

---

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
| Paired short reads | Reviewed pairing, complete pair/QC validation, optional fastp trimming, cancellable native SKESA, atomic assembly/provenance | Native tool must pass its platform gate; trimming is optional and does not validate an isolate; no SPAdes or long-read pipeline |
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

The complete Kleborate/Kaptive, AMRFinderPlus, agr/SCCmec/spa, MOB-suite and
abricate execution stack is not included. Species-complex resolution,
contamination quantification, long-read assembly, validated phenotype flags and
multi-user clinical deployment remain separate work.

Read preprocessing is no longer absent: fastp 1.3.7 is bundled where a verified
native artifact exists for the platform, and its report is shown as fastp's own
numbers. It remains **optional** — an assembly records whether it used the
trimmed pair or the originals — and trimming an isolate's reads validates
nothing about that isolate.

MOB-suite specifically cannot be bundled or run in this portable application, and
that was verified rather than assumed: its reference database alone is a 473 MB
download, and its pipeline needs a second Python runtime with SciPy, pandas,
PyTables, ete3, pycurl, Biopython and a Windows Mash binary that does not exist.
The plasmid evidence screen reports contig-level replicon and determinant
placement instead, and says in the application that it is not MOB-suite: no
relaxase or MPF type, no oriT, no plasmid reconstruction, no mobility prediction.

The 0.3 organism modules narrow that gap without closing it, and the distinction
matters. WMLSTudio runs its own exact-allele and BLAST+ screens against *pinned
public reference data* from those projects; it does not execute those tools and
does not reproduce their outputs. Kleborate's aggregate virulence and resistance
scores are not computed. Kaptive's whole K and O locus references, its
match-confidence grading and its O-locus special logic are not implemented, and
the deferral is recorded with its reason in
[organism modules](ORGANISM_MODULES.md). No SCCmec result is an MRSA designation,
and mecC is outside the pinned panel entirely. spa typing, agr typing and
MOB-recon remain absent.

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
