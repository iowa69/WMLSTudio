# WMLSTudio dependency and reference-data notices

This desktop distribution includes Python, PySide6 Essentials (Qt for Python),
Shiboken, and pyahocorasick. Installed distribution metadata and applicable
packaged license files are retained under `_internal`. Version-matched Python,
PySide and QtBase license texts are staged explicitly under
`_internal/notices/licenses`, with source URLs and SHA-256 hashes. Python uses the PSF
license; pyahocorasick uses BSD-3-Clause. PySide6 and Shiboken provide LGPLv3,
GPLv3, and commercial licensing options. Qt components have their own notices.
See the authoritative [Qt licensing page](https://www.qt.io/licensing/) and the
license files accompanying the exact binary versions in this package. The
Windows target interpreter's complete `LICENSE.txt` is preserved separately
as `Python-Windows-runtime-LICENSE.txt`, including its bundled OpenSSL, bzip2,
libffi and native runtime notices; the shorter CPython source LICENSE alone
does not replace those notices. Windows Qt uses its native Schannel TLS
backend. The optional Qt OpenSSL backend and ambient runner-discovered
`libcrypto-3-x64.dll` / `libssl-3-x64.dll` are excluded. Python retains its own
official-interpreter `libcrypto-3.dll` / `libssl-3.dll` for explicit downloads.
The build preserves license/NOTICE files from all installed isolated build
distributions, including incidental transitive helper modules collected by
PyInstaller. The one-folder layout keeps Qt libraries separate. Distribution maintainers must
retain notices and satisfy the applicable license obligations, including
corresponding-source and replacement/relinking requirements where applicable.

The application implementation does not incorporate WMLST source code. The
build script copies only reference FASTA/profile/metadata files from the selected
WMLST database snapshot, crediting [WMLST](https://github.com/iowa69/WMLST),
Giovanni Lorenzin / IOWA-Tech, as the immediate source of that snapshot.
Software licensing does not establish the redistribution rights of databases.

The bundled allele sequences and profiles originate from
[PubMLST](https://pubmlst.org/), University of Oxford, and their organism-specific
scheme curators. Cite the original scheme publication and Jolley KA, Bray JE,
Maiden MCJ (2018), *Open-access bacterial population genomics: BIGSdb software,
the PubMLST.org website and their applications*, Wellcome Open Research 3:124,
[doi:10.12688/wellcomeopenres.14826.1](https://doi.org/10.12688/wellcomeopenres.14826.1).
Publications should include PubMLST's required acknowledgement from its terms.

The [PubMLST terms](https://pubmlst.org/terms-conditions), checked 2026-09-12,
distinguish submission date: pre-2025 records may be used and redistributed with
scientific citation and acknowledgement; post-2024 downloaded data have use and
redistribution restrictions. rMLST has separate terms. The cache includes
download/authentication metadata, but lacks independently verified submission
dates for individual records. An anonymous-download flag or a download date
does not independently prove a record's redistribution rights. Before public
distribution, establish the rights for the actual snapshot, retain evidence,
and use a rights-cleared snapshot where necessary. This build does not relabel
database content as MIT, BSD, or public domain.

`_internal/wmlstudio/resources/schemes/manifest.json` records the immediate
source, repository revision when available, per-file SHA-256 checksums, preserved
scheme metadata, snapshot SHA-256, and archive SHA-256 for downloaded builds.
The snapshot is static: the application makes no automatic database downloads.

## Added native analysis components

SKA2 0.5.1 is Apache-2.0, with original LICENSE/NOTICE retained in `Tools/ska2`.
The Windows build uses pinned source `fcf9413d2768dc6538d31f664a4bf310651449e7`,
Rust 1.90.0 and a static Microsoft C runtime; its manifest, Cargo.lock,
vendored dependency sources/notices and exact source archive travel together.
The Linux tool is the separate official upstream Linux release, hash-pinned
with its original notices and source snapshot; the Windows dependency lock
is not represented as the Linux build's lock. SKA2 split-kmer SNP distances
are a separate assay, not cgMLST allele distances or a transmission verdict.

FastQC 0.12.1 is the original Babraham distribution, pinned to the binary
archive checksum in `Tools/fastqc/manifest.json`. It is GPL-3.0-or-later;
its GPL and component license texts remain in that directory. The exact
upstream source revision `e7ef390bf10382f60786bdd0cf28abd4f8683ffd` is included
under `Tools/fastqc/sources`. Bundled JAR archives retain their component
license/NOTICE resources; FastQC is not represented as MIT software.
FastQC runs in its own process with app-local Eclipse Temurin OpenJDK
17.0.20.1+1, without requiring a system Java installation. The original
JRE legal tree, NOTICE, release information and complete matching OpenJDK
source archive accompany it. OpenJDK is GPL-2.0 with the Classpath Exception
and applicable additional component terms; consult `Tools/fastqc/jre/legal`.
No custom QC calculation is labeled FastQC, and FastQC reports do not imply
automatic trimming or clinical suitability. The full source archives are
intentionally retained despite the additional portable archive size.

HYDRA 1.4.0 is bundled from the pinned source revision
[`6d36c109491c16544e8919fe6962b4b62e97d3d7`](https://github.com/iowa69/hydra/tree/6d36c109491c16544e8919fe6962b4b62e97d3d7)
under its MIT license. Its original assembly calling pipeline runs in the
dedicated `WMLSTudio-HYDRA` process. pandas, NumPy, python-dateutil, six and
archspec license texts, including applicable bundled-component notices from
the target-platform wheels, are retained in `notices/licenses/packages`.

Pyrodigal 3.7.1 is **GPL-3.0-or-later**, not MIT. Its exact provider-verified
source distribution, including the underlying Prodigal sources and build
configuration, is included at `notices/licenses/sources/pyrodigal-3.7.1.tar.gz`;
the complete GPL text is retained with its package license. WMLSTudio's own
source files retain their stated license, but the combined native distribution
must not be represented as an MIT-only product. Preserve and provide the
corresponding application/source/build materials and satisfy all applicable
GPL requirements before further distribution. Public distribution still
requires a review of the complete bundle and reference-data rights.
The exact application/build sources are included in
`notices/licenses/sources/wmlstudio-application-source.zip`; the pinned HYDRA
source archive is beside it. These archives include source and build recipes,
not raw sample files.

Windows SKESA 2.4.0 is a native UCRT64 build of pinned NCBI source commit
`c1413581e4f37211892d3c4310d01f3d9a9b3490`, with the reviewed portability patch.
Its NCBI public-domain portions coexist with AGPL-3.0-or-later GATB portions;
the upstream license and its stated exceptions remain authoritative. The
complete modified corresponding source, build recipes, GNU AGPL text and
native dependency licenses accompany `skesa.exe` in
`wmlstudio/resources/tools/skesa`. Preserve this complete directory, not only
the executable. Its manifest records compiler, build flags, source revision,
CPU requirements, DLL dependencies and file checksums.

The native BLAST+ 2.17.0 files are from the official
[NCBI distribution](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.17.0/).
Their archive and extracted-file SHA-256 values are recorded in
`Tools/blast/manifest.json`; `Tools/blast/LICENSE` and `BLAST_PRIVACY` travel
with the binaries. NCBI's public-domain notice does not erase the separate
licenses of bundled third-party components.
BLAST's required Microsoft VC++ runtime DLLs are copied app-locally from the
official target-version PySide6 Windows wheel, not from a developer machine's
System32 or Wine libraries. Their exact provider and hashes are recorded in
`notices/licenses/manifest.json`. The dependency-closure gate rejects missing
redistributable libraries even when the build runner has them installed.

The NCBI-only HYDRA starter contains the AMRFinderPlus nucleotide, protein and
mutation reference data, with attribution to NCBI and Feldgarden et al.
(2021), *Scientific Reports* 11:12728. It is not the AMRFinderPlus executable
or a claim of identical AMRFinderPlus results. HYDRA's provider registry
describes these NCBI resources as public domain; provider terms remain
authoritative. No CARD, VFDB or other independently licensed provider data are
silently included by the HYDRA starter-staging command. Its reference files,
mutation companions, manifest and observed provider release are hashed in
`wmlstudio/resources/hydra/starter/snapshot_provenance.json`.

The bundled starter is pinned to NCBI AMRFinderPlus reference release
2026-08-07.1, staged 2026-09-12. The application reports that release, the day it
was staged and its age in days wherever the store is shown — the AMR reference
databases window, `studio_scripts/check_setup.py`, and every HYDRA report's
execution provenance — so a user can see how stale the bundled evidence is
before reporting from it. Updating publishes a new snapshot beside the old one:
nothing is replaced, and analyses already recorded keep the reference snapshot
they were run against.

## Bundled organism-module reference panels

The independent species/virulence starter and the organism-specific typing
panels are staged by `studio_scripts/stage_characterization.py` from three
public repositories, each pinned to one commit and each fetched together with
its own `LICENSE`. The staged snapshot records every upstream URL, byte count
and SHA-256, and its `sources` block repeats the repository, revision and
licence identifier for each source. Deriving a panel (splitting a multi-record
FASTA, rewriting headers, parsing a rules table) happens at staging time only;
the upstream files are retained verbatim beside the derived artefacts so the
derivation is auditable against the bytes that were actually downloaded.

[Kleborate](https://github.com/klebgenomics/Kleborate) is **GPL-3.0-or-later**.
Its species-reference accession set, virulence-locus allele FASTAs and
`profiles.tsv` locus-ST tables are pinned to commit
`550ce22a2c01c76064f4dabf403704ee2293356e`. The fetched `LICENSE` (35,141 bytes,
SHA-256 `589ed823e9a84c56feb95ac58e7cf384626b9cbf4fda2a907bc36e103de1bad2`)
travels in the snapshot as `source-LICENSE` (`source-LICENSE-kleborate` in a
format 2 snapshot). WMLSTudio does not incorporate Kleborate source code, does
not run Kleborate, and does not compute Kleborate's aggregate virulence or
resistance scores. A locus ST reported here is an exact-allele lookup against
that pinned profile table, and a lineage string is that table's own value,
reported verbatim.

[rpetit3/sccmec](https://github.com/rpetit3/sccmec) v1.2.0 is **MIT**,
"Copyright (c) 2024 Robert A. Petit III", pinned to commit
`b901cc618be8eb17284ccb0cf6ef9ee428d909c3`. Its `LICENSE` is 1,076 bytes,
SHA-256 `5545cae984ae5abd56b68d66ddf1811be5219ec3804c237da14a9eee067d0e82`,
staged as `source-LICENSE-sccmec`. The staged panel derives from
`data/sccmec-targets.fasta` (107,660 bytes, SHA-256
`4b18b4f97651345b389f26dd8e91c95b792d32841dd2aee64bc7842884aa6c4d`),
`data/sccmec-targets.yaml` (3,638 bytes, SHA-256
`feec4363b42076c4c6187b7abeae33860ea91f9e34393024db67e899e2bba03d`),
`data/sccmec-regions.fasta` (1,121,499 bytes, SHA-256
`8bce1de540374f487d3f2144292b5f32872cdd839acb86bad2dc84eb0cee789f`) and the two
companion TSV tables, all retained verbatim. The IWG ccr/mec type definitions
live in the staged manifest, not in WMLSTudio source, so the rules in force are
provable from the recorded `reference_digest`. The cassette references are
public GenBank records whose accessions remain in their FASTA headers. This is a
BLAST+ marker screen against that pinned panel: it is not staphopia-sccmec or
SCCmecFinder output, is not equivalent to them, and assigns no MRSA/MSSA
designation or methicillin susceptibility. `mecC` is absent from the upstream
20-target set and is recorded as `not_assayed`, never as absent.

[Kaptive](https://github.com/klebgenomics/Kaptive) v2.0.9 is
**GPL-3.0-or-later** (GNU GPL v3), pinned to commit
`b3856eac6e76b3017aa993319da2a8ea967a1ba0`. Only
`reference_database/wzi_wzc_db.fasta` (246,938 bytes, SHA-256
`5349423a9cbeedbce35ea499b441a23f1a965d64d265bdc29c96713e775e820d`) is staged;
the K and O locus GenBank databases are not. The fetched `LICENSE` is staged as
`source-LICENSE-kaptive` with its byte count and SHA-256 recorded in the
snapshot manifest's `files` list. WMLSTudio reports a `wzi` or `wzc` allele
number only. No wzi-allele-to-K-type mapping is shipped or applied, no K or O
locus is assigned, and Kaptive's match-confidence grading and O-locus logic are
not implemented. This is not Kaptive output and is not equivalent to it.

A snapshot staged before these panels existed is format version 1. It keeps
validating unchanged, and the organism-specific assays then report `not_run`
naming the missing manifest section. That is a diagnosable reference state, not
a negative result.

## Public genome references downloaded on request

Neither the practice cohorts nor the broad species panel is included in the
portable ZIP. Both are fetched from NCBI to the user's own adjacent `Data`
directory after an explicit request, verified against pinned checksums, and the
build refuses to package a copy that was staged into the source tree. Only
accession and checksum tables are held in WMLSTudio; no sequence bytes are
redistributed by this project.

> Assemblies are public NCBI RefSeq records retrieved from
> https://ftp.ncbi.nlm.nih.gov/genomes/all/. NCBI places no restrictions on the
> use or distribution of the data it hosts, but it does not hold their copyright
> and cannot grant rights on behalf of the depositing submitters; individual
> submitters may assert terms. Cite the assembly accession and the originating
> submitters, not WMLSTudio, when reusing these sequences.

The broad species panel is pinned as revision
`ncbi-refseq-species-panel-2026-09-14.2`: 18 RefSeq assemblies, about 16.5 MiB,
each row pinning the accession, assembly directory, compressed byte count,
compressed SHA-256 and decompressed-FASTA SHA-256. It is a triage panel, not a
representation of within-species diversity, and it does not distinguish
*Escherichia coli* from *Shigella*. Most taxa carry one assembly;
*Listeria monocytogenes* carries two, because its lineages straddle the 95% ANI
species line and a single reference left common genomes unresolved.

The practice cohorts are pinned as content digests over their accession tables:
`kpneumoniae-10` (10 assemblies, 16,746,561 bytes, digest
`238aae60e7f48449aec656cf4e4a500c6ad0136d6e9dc71ecf4fee3c545e6dbc`) and
`mixed-genus-20` (20 assemblies, 21,703,421 bytes, digest
`459071d13922ec70677618e35c5ddab0b1629cf46a8fec6b53424e3b8a5b5b97`). Each
download re-checks the pinned size, MD5 and SHA-256 against the checksums NCBI
publishes today and aborts on any disagreement. Organism and strain labels are
the ones NCBI records for those assemblies; they were not independently
verified here, and no expected ST, cluster or threshold is shipped with either
cohort. Full provenance and per-accession tables are in
[docs/TEST_DATASETS.md](../docs/TEST_DATASETS.md).

Database downloads/updates occur only after explicit user action and publish
a new versioned snapshot. Existing reference snapshots are retained. Database
access and software licensing must not be confused with clinical validation:
gene/mutation evidence does not establish a susceptible/resistant phenotype.
