#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-2.0-only
# Copyright (C) 2025-2026 IOWA-BioTech - Giovanni Lorenzin
# Copyright (C) Torsten Seemann (upstream `mlst`, from which WMLST is ported)
"""Fail the build if any module under ``wmlst/`` imports a non-stdlib package.

Implements docs/ARCHITECTURE.md section 13.1 and the hard rule of section 13.3:
WMLST is GPL-2.0-only, and a GPLv2 work may not link an Apache-2.0 library
(the patent-termination clause is a further restriction GPLv2 section 6
forbids). Zero mandatory runtime dependencies is therefore a *licensing*
invariant, not a packaging preference, and it is enforced here rather than
trusted to reviewers.

Exactly one structural exception is allowed: ``tkinterdnd2`` (MIT), and only
when its import sits inside a ``try: ... except ImportError: ...`` guard so the
GUI degrades to a Browse button when the package is absent.

Usage::

    python3 scripts/check_stdlib_only.py [path ...]     # default: wmlst

Exit code 0 when clean, 1 when a violation is found.
"""

from __future__ import annotations

import argparse
import ast
import os
import sys

#: The single sanctioned optional third-party import (section 2.1).
OPTIONAL_ALLOWED = frozenset({"tkinterdnd2"})

#: Stdlib packages that exist only on a NEWER Python than the 3.9 floor.
#:
#: ``sys.stdlib_module_names`` describes the interpreter running this script, so
#: a module added after that version looks third-party here even though it ships
#: with CPython. Checking on 3.14 would pass and on 3.12 would fail, which is how
#: this rule first broke CI. These are stdlib, but because they are absent on the
#: floor they must still sit behind ``try/except ImportError``.
VERSIONED_STDLIB = {
    "compression": (3, 14),   # compression.zstd
    "tomllib": (3, 11),
    "zoneinfo": (3, 9),
    "graphlib": (3, 9),
}

#: Package roots that are first-party and therefore always fine.
FIRST_PARTY = frozenset({"wmlst", "wmlst_db", "db"})

# sys.stdlib_module_names landed in 3.10. On the 3.9 floor we fall back to a
# list derived from the same source, plus the handful of private helpers the
# stdlib itself re-exports.
_EXTRA_STDLIB = frozenset(
    {
        "__future__", "_abc", "_ast", "_bootlocale", "_collections_abc",
        "_compat_pickle", "_compression", "_csv", "_ctypes", "_datetime",
        "_decimal", "_frozen_importlib", "_functools", "_hashlib", "_heapq",
        "_imp", "_io", "_json", "_locale", "_lzma", "_markupbase", "_md5",
        "_multiprocessing", "_opcode", "_operator", "_osx_support", "_pickle",
        "_posixsubprocess", "_py_abc", "_pydecimal", "_pyio", "_queue",
        "_random", "_sha1", "_sha256", "_sha512", "_signal", "_sitebuiltins",
        "_socket", "_sqlite3", "_sre", "_ssl", "_stat", "_string", "_strptime",
        "_struct", "_thread", "_threading_local", "_tkinter", "_tracemalloc",
        "_warnings", "_weakref", "_weakrefset", "_winapi", "_bz2",
    }
)


def _stdlib_names() -> frozenset:
    """Return the set of standard-library top-level module names."""
    names = getattr(sys, "stdlib_module_names", None)
    if names is None:  # pragma: no cover - only on Python 3.9
        import distutils.sysconfig as _sc  # noqa: F401  (presence probe only)

        names = frozenset(sys.builtin_module_names) | _PY39_STDLIB
    return frozenset(names) | _EXTRA_STDLIB


# Python 3.9's stdlib top-level names, for the one interpreter that lacks
# sys.stdlib_module_names. Generated from the 3.9 module index.
_PY39_STDLIB = frozenset(
    """abc aifc antigravity argparse array ast asynchat asyncio asyncore atexit
    audioop base64 bdb binascii binhex bisect builtins bz2 cProfile calendar cgi
    cgitb chunk cmath cmd code codecs codeop collections colorsys compileall
    concurrent configparser contextlib contextvars copy copyreg crypt csv ctypes
    curses dataclasses datetime dbm decimal difflib dis distutils doctest email
    encodings ensurepip enum errno faulthandler fcntl filecmp fileinput fnmatch
    formatter fractions ftplib functools gc genericpath getopt getpass gettext
    glob graphlib grp gzip hashlib heapq hmac html http idlelib imaplib imghdr
    imp importlib inspect io ipaddress itertools json keyword lib2to3 linecache
    locale logging lzma mailbox mailcap marshal math mimetypes mmap modulefinder
    msilib msvcrt multiprocessing netrc nis nntplib nt ntpath nturl2path numbers
    opcode operator optparse os ossaudiodev parser pathlib pdb pickle pickletools
    pipes pkgutil platform plistlib poplib posix posixpath pprint profile pstats
    pty pwd py_compile pyclbr pydoc pydoc_data pyexpat queue quopri random re
    readline reprlib resource rlcompleter runpy sched secrets select selectors
    shelve shlex shutil signal site smtpd smtplib sndhdr socket socketserver
    spwd sqlite3 sre_compile sre_constants sre_parse ssl stat statistics string
    stringprep struct subprocess sunau symbol symtable sys sysconfig syslog
    tabnanny tarfile telnetlib tempfile termios textwrap this threading time
    timeit tkinter token tokenize trace traceback tracemalloc tty turtle
    turtledemo types typing unicodedata unittest urllib uu uuid venv warnings
    wave weakref webbrowser winreg winsound wsgiref xdrlib xml xmlrpc zipapp
    zipfile zipimport zlib""".split()
)

STDLIB = _stdlib_names()


class Violation:
    """One offending import: the file, the line, the module and why."""

    __slots__ = ("lineno", "module", "path", "reason")

    def __init__(self, path, lineno, module, reason):
        self.path = path
        self.lineno = lineno
        self.module = module
        self.reason = reason

    def __str__(self):
        return "{0}:{1}: {2} -- {3}".format(
            self.path, self.lineno, self.module, self.reason
        )


def _root_of(module_name):
    """Return the top-level package name of a dotted module path."""
    return (module_name or "").split(".", 1)[0]


def _guarded_nodes(tree):
    """Return the set of ``ast`` import nodes that sit inside a try/except ImportError.

    Only such imports may name an optional third-party package (section 2.1).
    """
    guarded = set()
    for node in ast.walk(tree):
        if not isinstance(node, ast.Try):
            continue
        catches_import_error = False
        for handler in node.handlers:
            exc = handler.type
            candidates = []
            if isinstance(exc, ast.Tuple):
                candidates = list(exc.elts)
            elif exc is not None:
                candidates = [exc]
            elif exc is None:
                catches_import_error = True  # bare except catches everything
            for cand in candidates:
                name = getattr(cand, "id", None) or getattr(cand, "attr", None)
                if name in ("ImportError", "ModuleNotFoundError", "Exception"):
                    catches_import_error = True
        if not catches_import_error:
            continue
        for stmt in node.body:
            for sub in ast.walk(stmt):
                if isinstance(sub, (ast.Import, ast.ImportFrom)):
                    guarded.add(sub)
    return guarded


def check_file(path):
    """Return a list of Violation for one ``.py`` file. Implements section 13.1."""
    with open(path, "rb") as fh:
        source = fh.read()
    try:
        tree = ast.parse(source, filename=path)
    except SyntaxError as exc:
        return [Violation(path, exc.lineno or 0, "<syntax>", "cannot parse: %s" % exc)]

    guarded = _guarded_nodes(tree)
    out = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            modules = [alias.name for alias in node.names]
            relative = False
        elif isinstance(node, ast.ImportFrom):
            relative = (node.level or 0) > 0
            modules = [node.module or ""]
        else:
            continue

        if relative:
            continue  # `from .engine import ...` is always first-party

        for module in modules:
            root = _root_of(module)
            if not root or root in STDLIB or root in FIRST_PARTY:
                continue
            if root in VERSIONED_STDLIB:
                if node in guarded:
                    continue
                major, minor = VERSIONED_STDLIB[root]
                out.append(
                    Violation(
                        path,
                        node.lineno,
                        module,
                        "stdlib only from Python %d.%d; must sit inside "
                        "try/except ImportError to keep the 3.9 floor working"
                        % (major, minor),
                    )
                )
                continue
            if root in OPTIONAL_ALLOWED:
                if node in guarded:
                    continue
                out.append(
                    Violation(
                        path,
                        node.lineno,
                        module,
                        "optional dependency must sit inside try/except ImportError",
                    )
                )
                continue
            out.append(
                Violation(
                    path,
                    node.lineno,
                    module,
                    "third-party import; WMLST has zero mandatory runtime deps "
                    "(GPL-2.0-only, see NOTICE section 4)",
                )
            )
    return out


def iter_python_files(targets):
    """Yield every ``.py`` path under the given files or directories, sorted."""
    for target in targets:
        if os.path.isfile(target):
            if target.endswith(".py"):
                yield target
            continue
        for dirpath, dirnames, filenames in os.walk(target):
            dirnames[:] = sorted(
                d for d in dirnames if d not in ("__pycache__", ".git", "build", "dist")
            )
            for name in sorted(filenames):
                if name.endswith(".py"):
                    yield os.path.join(dirpath, name)


def main(argv=None):
    """Command-line entry point. Returns a process exit code."""
    parser = argparse.ArgumentParser(
        prog="check_stdlib_only.py",
        description="Assert that wmlst/ imports nothing outside the standard library.",
    )
    parser.add_argument(
        "targets",
        nargs="*",
        default=None,
        help="files or directories to scan (default: the wmlst package)",
    )
    parser.add_argument("-q", "--quiet", action="store_true", help="only report failures")
    ns = parser.parse_args(argv)

    targets = ns.targets
    if not targets:
        here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        targets = [os.path.join(here, "wmlst")]

    files = list(iter_python_files(targets))
    if not files:
        sys.stderr.write("check_stdlib_only: no Python files found in %r\n" % (targets,))
        return 1

    violations = []
    for path in files:
        violations.extend(check_file(path))

    if violations:
        sys.stderr.write("check_stdlib_only: FAIL\n")
        for violation in violations:
            sys.stderr.write("  %s\n" % violation)
        sys.stderr.write(
            "\n%d violation(s). WMLST is GPL-2.0-only with zero mandatory runtime\n"
            "dependencies (docs/ARCHITECTURE.md 13.1, 13.3). The only sanctioned\n"
            "optional import is tkinterdnd2 behind try/except ImportError.\n"
            % len(violations)
        )
        return 1

    if not ns.quiet:
        sys.stdout.write(
            "check_stdlib_only: OK - %d file(s) scanned, stdlib only\n" % len(files)
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
