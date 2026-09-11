# SPDX-License-Identifier: GPL-2.0-only
# Copyright (C) 2025-2026 IOWA-Tech - Giovanni Lorenzin
# Copyright (C) Torsten Seemann (upstream `mlst`, from which WMLST is ported)
"""Anchor package for the bundled PubMLST database (docs/ARCHITECTURE.md 13.1).

Mapped to the importable name ``wmlst_db`` by pyproject.toml, so the 162
schemes can be located with importlib.resources from a wheel, an editable
install or a frozen bundle alike.
"""

__all__ = []
