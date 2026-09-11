# SPDX-License-Identifier: GPL-2.0-only
# Copyright (C) 2025-2026 IOWA-Tech - Giovanni Lorenzin
# Copyright (C) Torsten Seemann (upstream `mlst`, from which WMLST is ported)
"""PyInstaller runtime hook: survive having no standard streams.

Implements docs/ARCHITECTURE.md section 13.1.

Under ``--windowed`` (``console=False``) a frozen Windows application is started
with no console, so ``sys.stdout``, ``sys.stderr`` and ``sys.stdin`` are all
``None``. Any library that calls ``print()``, ``sys.stderr.write()`` or
``logging.StreamHandler()`` then dies with
``AttributeError: 'NoneType' object has no attribute 'write'`` -- typically deep
inside a worker thread, where the traceback goes nowhere because there is no
stderr to write it to. Substituting real file-like objects is the only reliable
cure; the GUI reads them back when the user opens the log pane.

It also forces UTF-8 on the child-process environment. BLAST+ writes UTF-8 and
the database carries non-ASCII scheme descriptions (``MLST#2 (Bekő)``), while a
frozen process on an Italian or Polish Windows install otherwise inherits
cp1252 and mangles them.

This module runs before any application code, so it must import nothing beyond
the standard library and must never raise.
"""

import io
import os
import sys


class _NullStream(io.TextIOWrapper):
    """A real text stream over an in-memory buffer, standing in for a missing one.

    Subclassing TextIOWrapper (rather than handing back a bare StringIO) means
    ``.buffer``, ``.encoding``, ``.reconfigure()`` and the absence of a usable
    ``.fileno()`` all behave the way ``cli.reconfigure_streams`` expects
    (section 4.8). ``name`` is a read-only descriptor on the base class, so it
    is re-exposed as a property rather than assigned.
    """

    def __init__(self, name, write_through=True):
        super().__init__(
            io.BytesIO(), encoding="utf-8", errors="replace", newline="\n",
            write_through=write_through,
        )
        self._stream_name = name

    @property
    def name(self):
        return self._stream_name

    def isatty(self):
        return False


def _install():
    os.environ.setdefault("PYTHONIOENCODING", "utf-8")
    os.environ.setdefault("PYTHONUTF8", "1")
    # NCBI's opt-out for the BLAST+ usage-reporting beacon. Set here as well as
    # in blastbin so a frozen GUI never makes an unexpected outbound request.
    os.environ.setdefault("BLAST_USAGE_REPORT", "false")

    for name in ("stdout", "stderr", "stdin"):
        if getattr(sys, name, None) is None:
            stream = _NullStream("<%s>" % name)
            setattr(sys, name, stream)
            setattr(sys, "__%s__" % name, stream)


try:
    _install()
except Exception:  # pragma: no cover - a hook must never break start-up
    pass
