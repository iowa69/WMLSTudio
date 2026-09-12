# SPDX-License-Identifier: GPL-2.0-only
# Copyright (C) 2025 IOWA-BioTech - Giovanni Lorenzin
# Copyright (C) Torsten Seemann (upstream `mlst`, GPL-2.0-only)
"""WMLST - native Windows MLST typing from assembled contigs.

A Python port of `mlst` 2.35.0 by Torsten Seemann. See docs/ARCHITECTURE.md.

This package exports only ``__version__``; everything else lives in the
individual modules so that importing :mod:`wmlst` costs nothing.
"""

from .version import __version__

__all__ = ["__version__"]
