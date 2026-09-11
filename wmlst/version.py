# SPDX-License-Identifier: GPL-2.0-only
# Copyright (C) 2025-2026 IOWA-Tech - Giovanni Lorenzin
# Copyright (C) Torsten Seemann (upstream `mlst`, from which WMLST is ported)
"""Single source of version truth for WMLST.

WMLST is a Windows-native port of `mlst` by Torsten Seemann.
Keeping UPSTREAM_MLST_VERSION here documents exactly which upstream release
the byte-compatibility contract in docs/ARCHITECTURE.md was verified against.
"""

__version__ = "1.0.2"

#: The tseemann/mlst release whose output WMLST reproduces byte-for-byte.
UPSTREAM_MLST_VERSION = "2.35.0"

#: NCBI BLAST+ release that the Windows bootstrapper installs.
BLAST_VERSION = "2.17.0"

#: Bundled PubMLST database snapshot date (mirrors db/VERSION.txt).
BUNDLED_DB_VERSION = "2025-12-29"

__all__ = [
    "BLAST_VERSION",
    "BUNDLED_DB_VERSION",
    "UPSTREAM_MLST_VERSION",
    "__version__",
]
