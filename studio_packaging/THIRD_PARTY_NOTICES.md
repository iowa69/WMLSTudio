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

Database downloads/updates occur only after explicit user action and publish
a new versioned snapshot. Existing reference snapshots are retained. Database
access and software licensing must not be confused with clinical validation:
gene/mutation evidence does not establish a susceptible/resistant phenotype.
