# HYDRA integration assessment

HYDRA is a useful scientific engine candidate for WMLSTudio. Its scope covers
acquired AMR and virulence genes, mutations, read-derived allele fractions,
MLST and lineage typing. The integration implemented here imports its JSON
reports into the native workspace. Running the complete engine inside a
portable Windows distribution needs additional work.

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
| Assembly nucleotide/protein search | BLAST+ executables, including translated searches | Native Windows BLAST exists; pinned binaries and engine-level validation are feasible next work |
| Read mode | minimap2, samtools, pysam/HTSlib | A Python wheel alone does not bundle these tools or establish a Windows read pipeline |
| Species sketch evidence | Optional Mash plus sketches | Must package/test a supported executable or visibly report that evidence unavailable |
| Reference databases | HYDRA database manager and explicit downloads | Need versioned manifests, licensing review and reproducible indexes |

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

The next execution milestone should be a pinned native assembly mode worker,
launched with argument arrays and isolated working directories, with explicit
capability detection, cancellable processes, captured diagnostics, and hashed
database manifests. Direct read mode should follow only after the complete
alignment and pileup stack passes Windows validation on paired real data.

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
WMLSTudio's JSON importer neither copies HYDRA source nor bundles those databases.
