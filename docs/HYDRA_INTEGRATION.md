# HYDRA integration assessment

HYDRA is a useful scientific engine candidate for WMLSTudio. Its scope covers
acquired AMR and virulence genes, mutations, read-derived allele fractions,
MLST and lineage typing. The native workspace imports its JSON reports and
now executes the pinned upstream **assembly** pipeline through a dedicated
native worker and bundled BLAST+ tools. Direct-read alignment/pileup and the
complete upstream typing stack are not enabled by this integration.

## Source and implemented contract

This assessment inspected HYDRA 1.4.0 at commit
[`6d36c109491c16544e8919fe6962b4b62e97d3d7`](https://github.com/iowa69/hydra/tree/6d36c109491c16544e8919fe6962b4b62e97d3d7),
specifically its [output writer](https://github.com/iowa69/hydra/blob/6d36c109491c16544e8919fe6962b4b62e97d3d7/src/hydra_amr/report/writer.py)
and [record definitions](https://github.com/iowa69/hydra/blob/6d36c109491c16544e8919fe6962b4b62e97d3d7/src/hydra_amr/records.py).
The JSON shape is:

```text
report
  hydra_version, command, databases[], parameters{}
  samples[]
    sample, input_type, inputs[], runtime_seconds
    species{}, mlst{}, typing[], scores{}, qc{}, warnings[]
    hits[]: gene, database, element_type, class, method, resolution,
            primary, identity_pct, coverage_pct, depth, allele_fraction, ...
```

`wmlstudio.hydra.load_hydra_report(path)` validates that structure, preserves the
upstream fields, computes a source-file SHA-256 and adds import provenance.
Each sample receives primary-aware summary counts; secondary cross-database hits
remain visible without inflating gene counts. Catalogued `POINTR` evidence is
kept distinct from uncatalogued `VARIANTR` observations. Upstream warnings, QC,
methods, resolution and allele fractions are retained. `runtime_seconds` means
the upstream whole-batch runtime, not the processing time of one sample.

Unknown versions matching the supported structure receive an import warning.
Malformed JSON, duplicate keys or samples, conflicting hit owners and impossible
percentages/fractions are rejected with readable errors. Paths embedded in
reports remain text and are never opened. Import does not execute reported
commands, contact a database, recalculate a call or overwrite native MLST evidence.

## Execution dependencies and Windows feasibility

| Component | Observed upstream dependency | Integration consequence |
| --- | --- | --- |
| Python orchestration | pandas, NumPy; optional SciPy/openpyxl | Package/version and startup-size testing needed |
| Assembly nucleotide/protein search | BLAST+ executables, including translated searches | Official Windows/Linux BLAST+ 2.17.0 staged with pinned hashes; real assembly parity checked under Wine |
| Read mode | minimap2, samtools, pysam/HTSlib | A Python wheel alone does not bundle these tools or establish a Windows read pipeline |
| Species sketch evidence | Optional Mash plus sketches | Must package/test a supported executable or visibly report that evidence unavailable |
| Reference databases | HYDRA database manager and explicit downloads | NCBI starter staged; explicit updates publish immutable snapshots with content hashes |

The dependency evidence is in
[pyproject.toml](https://github.com/iowa69/hydra/blob/6d36c109491c16544e8919fe6962b4b62e97d3d7/pyproject.toml),
[the reads engine](https://github.com/iowa69/hydra/blob/6d36c109491c16544e8919fe6962b4b62e97d3d7/src/hydra_amr/engines/reads.py),
[the BLAST wrapper](https://github.com/iowa69/hydra/blob/6d36c109491c16544e8919fe6962b4b62e97d3d7/src/hydra_amr/engines/blast.py),
and [species typing](https://github.com/iowa69/hydra/blob/6d36c109491c16544e8919fe6962b4b62e97d3d7/src/hydra_amr/typing/species.py).
The read engine pipes minimap2 output into samtools sort, indexes the BAM and
reads evidence through pysam. Each part requires native packaging and failure,
cancellation, path and parity tests. This is an engineering assessment, not a
claim that those components are impossible to run on Windows.

[NCBI documents native Windows BLAST installation](https://www.ncbi.nlm.nih.gov/books/NBK52637/).
[pysam's installation documentation](https://pysam.readthedocs.io/en/latest/installation.html)
explains its compiled HTSlib dependency. A source checkout of HYDRA by itself does
not deliver an all-in-one native executable. WSL would change the user experience
specified for this project and is not used by this import integration.

`wmlstudio.hydra_runtime.run_assemblies` launches HYDRA with argument arrays in
an isolated temporary workspace, captures diagnostics, validates the resulting
JSON and preserves input, reference and tool provenance. Cancellation terminates
the owned worker and its BLAST subprocesses; partial reports are not imported.
Analysis never downloads references. The frozen worker is `WMLSTudio-HYDRA.exe`.
`update_databases` downloads only after explicit user action, stages a new
selected-provider snapshot, validates it and atomically publishes it under a
new versioned directory. Existing snapshots are retained on success or failure.

## The reference catalogue, and what is bundled

The engine's own registry knows fifteen reference sets. Until this revision the
application silently used two of them and named none of the rest, which made a
missing database indistinguishable from a database that does not exist.
`hydra_runtime.database_catalogue` now returns every set the pinned engine can
use, installed or not, with the provider, title, purpose, licence, citation and
upstream address the registry records, plus whether the application knows how to
fetch it automatically or the user must obtain it by hand.

Only the NCBI sets are bundled, and they are bundled **whole**: nucleotide,
protein and point mutations. Packaging refuses a starter snapshot that lacks
`AMRProt-mutation.tsv` or the per-organism DNA catalogues under `mutation/dna`,
because a package that could screen only for acquired genes would report nothing
for point mutations, and a reader cannot tell that from a negative result. The
bundled release carries 13 DNA catalogues; 30 organisms have curated protein
mutations; their union is 31. The engine accepts 32 organism names, so one
accepted organism — *Burkholderia mallei* — has no mutation catalogue at all and
is reported as exactly that, not as an organism with no mutations.

`element_counts` reads the bundled protein table directly: 10,078 records, of
which 8,794 are AMR, 1,025 virulence, 259 stress and 288 point-mutation entries.
Those counts are what decides whether offering a virulence search would offer
anything; a store with no protein reference returns zeros and the caller says so
rather than letting an unperformed search resemble a clean result.

"Install and update everything" plans the whole store before touching it and
excludes every set whose licence is not an open one — CARD's academic licence, and
any provider that records no licence. Those remain individually downloadable after
their terms are shown. An automatic action must not accept licences on a user's
behalf.

Protein is imported before nucleotide references so HYDRA's own importer can
transfer curated family/class annotations. The wrapper reads NCBI's
`version.txt` before and after download and refuses to publish if it changes.
The inspected upstream downloader omits that file, so the wrapper preserves its
original `unknown` version as `upstream_reported_version` and adds the
independently observed release and source URL. The staged NCBI starter is
release `2026-08-07.1`: 9,786 nucleotide and 10,078 protein sequences, plus
organism-specific DNA mutation companions. Every companion file and index
participates in the reference hash inventory.

Native WMLSTudio MLST results are not replaced by HYDRA's MLST/lineage modules.
No organism means no automatic mutation-catalog selection; the user can provide
an organism explicitly. The organism also decides whether virulence and stress
elements are searched: the default searches them where the isolate's organism is
established, and the two overrides — always, or never — are recorded in the run's
provenance together with the reason, so a report always says what was not looked
for as well as what was found. Upstream nucleotide defaults are 80% identity and 60%
reference coverage; translated defaults are 90% identity and 90% complete-hit
coverage. The latter is **not** a hard hit-exclusion floor: upstream can retain
partial protein evidence down to its 50% partial-coverage threshold. Explicit
thresholds override provider-specific defaults and must be reported.

## Validation actually performed

On 2026-09-12, the local S. epidermidis `example.fna` assembly (SHA-256
`a18f883d81d4a09350526706698a9e2a0637e3240f2acd638a60fa854c8ded39`)
was screened with the same staged references and HYDRA commit on Linux and
native Windows Python/BLAST under Wine 11.17. All **14 complete hit records
matched exactly**, including primary/secondary flags, percentages, methods,
classes and element types. This is a bounded cross-platform engine check, not
clinical sensitivity/specificity. Runtime was 29.8 s on Linux and 45.1 s under
Wine on this host. Actual Windows 11 hardware validation remains pending.
The JSON evidence is under `artifacts/` as
`hydra-real-sepidermidis-linux-verified.json` and
`hydra-real-sepidermidis-windows-wine.json`.

Direct read mode should follow only after the complete alignment and pileup
stack passes Windows validation on paired real data. A native read-assembly
workflow, where available, is separate from direct HYDRA read evidence.

## Scientific interpretation and licensing

HYDRA's README reports upstream benchmark results. Those are upstream claims,
not WMLSTudio validation, and have not been reproduced by this integration.
Import tests establish parsing, provenance and count semantics; they do not
establish sensitivity, specificity, phenotype accuracy or read-call concordance.
A mutation/gene detection is molecular evidence and must not be silently
converted into a susceptible/resistant phenotype in the desktop UI.

HYDRA's software package declares MIT licensing. Its
[database registry](https://github.com/iowa69/hydra/blob/6d36c109491c16544e8919fe6962b4b62e97d3d7/src/hydra_amr/db/registry.py)
lists separate providers and license descriptions, including CARD academic
licensing. Those descriptions should be checked against each provider's current
terms and the exact snapshot before redistribution. The MIT software license
does not transfer ownership of CARD, VFDB, PubMLST or other reference data.
The importer does not execute reported commands. The assembly worker bundles
the pinned MIT engine and NCBI-only starter; it does not implicitly package
CARD, VFDB or other provider resources. Other components have different terms,
notably GPL-3.0-or-later Pyrodigal; see the portable package notices and included
corresponding-source archive.
