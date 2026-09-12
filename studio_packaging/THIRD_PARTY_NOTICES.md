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
one-folder layout keeps Qt libraries separate. Distribution maintainers must
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
