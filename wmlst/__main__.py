# SPDX-License-Identifier: GPL-2.0-only
# Copyright (C) 2025 IOWA-BioTech - Giovanni Lorenzin
# Copyright (C) Torsten Seemann (upstream `mlst`, GPL-2.0-only)
"""``python -m wmlst`` entry point (docs/ARCHITECTURE.md section 2.1)."""

import sys

from .cli import main

if __name__ == "__main__":
    sys.exit(main())
