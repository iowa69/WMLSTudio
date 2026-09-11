#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-2.0-only
# Copyright (C) 2025-2026 IOWA-Tech - Giovanni Lorenzin
# Copyright (C) Torsten Seemann (upstream `mlst`, from which WMLST is ported)
"""Enforce the per-file licence header mandated by docs/ARCHITECTURE.md 13.3.

WMLST is a derivative work of ``mlst`` 2.35.0 by Torsten Seemann, which is
GPL-2.0-only with no "or any later version" clause. GPLv2 section 2(b) licenses
the derivative as a whole, so every source file must carry:

1. ``SPDX-License-Identifier: GPL-2.0-only``
2. a copyright line naming IOWA-Tech / Giovanni Lorenzin
3. a copyright line naming Torsten Seemann

and, for the modules that are direct line-by-line translations of upstream Perl,
4. a provenance reference naming the upstream file and line range it came from,
   e.g. ``Ported from bin/mlst:370-395``.

Usage::

    python3 scripts/check_headers.py [path ...]
    python3 scripts/check_headers.py --fix        # insert missing headers

Exit code 0 when every file complies, 1 otherwise.
"""

from __future__ import annotations

import argparse
import os
import re
import sys

SPDX = "GPL-2.0-only"

#: Number of leading lines that count as "the header region".
HEADER_LINES = 60

RE_SPDX = re.compile(r"SPDX-License-Identifier:\s*(?P<id>[^\s*#]+)")
RE_VENDOR = re.compile(r"IOWA[- ]?Tech", re.IGNORECASE)
RE_AUTHOR = re.compile(r"Giovanni\s+Lorenzin", re.IGNORECASE)
RE_UPSTREAM = re.compile(r"Torsten\s+Seemann", re.IGNORECASE)

#: `path:LINE-LINE` or `path:LINE`, e.g. bin/mlst:370-395, MLST/Scheme.pm:12
RE_PROVENANCE = re.compile(
    r"(?:bin/mlst|perl5/MLST/[A-Za-z]+\.pm|MLST/[A-Za-z]+\.pm|"
    r"scripts/mlst-[A-Za-z_]+|any2fasta|test/test\.sh)"
    r"(?::\d+(?:\s*-\s*\d+)?)?"
)

#: Modules that are direct translations and therefore need rule 4.
#: Keys are paths relative to the repository root, in POSIX form.
TRANSLATED_MODULES = frozenset(
    {
        "wmlst/engine.py",
        "wmlst/schemes.py",
        "wmlst/any2fasta.py",
        "wmlst/updatedb.py",
        "wmlst/cli.py",
        "wmlst/report.py",
    }
)

#: Trivial files that need the SPDX tag but no provenance prose.
SKIP_DIRS = frozenset({"__pycache__", ".git", "build", "dist", ".venv", "venv",
                       ".tox", ".mypy_cache", ".pytest_cache", "db"})

HEADER_TEMPLATE = (
    "# SPDX-License-Identifier: GPL-2.0-only\n"
    "# Copyright (C) 2025-2026 IOWA-Tech - Giovanni Lorenzin\n"
    "# Copyright (C) Torsten Seemann (upstream `mlst`, from which WMLST is ported)\n"
)


def header_region(path):
    """Return the first HEADER_LINES lines of a file as one string."""
    with open(path, encoding="utf-8", errors="replace") as fh:
        lines = []
        for i, line in enumerate(fh):
            if i >= HEADER_LINES:
                break
            lines.append(line)
    return "".join(lines)


def check_file(path, relpath):
    """Return a list of human-readable problems for one file (section 13.3)."""
    text = header_region(path)
    problems = []

    match = RE_SPDX.search(text)
    if match is None:
        problems.append("missing 'SPDX-License-Identifier: %s'" % SPDX)
    elif match.group("id") != SPDX:
        problems.append(
            "SPDX tag is %r, must be %r - WMLST cannot be relicensed, it is a "
            "GPLv2 derivative (13.3)" % (match.group("id"), SPDX)
        )

    if not (RE_VENDOR.search(text) and RE_AUTHOR.search(text)):
        problems.append("missing the 'IOWA-Tech - Giovanni Lorenzin' copyright line")

    if not RE_UPSTREAM.search(text):
        problems.append(
            "missing the upstream copyright line naming Torsten Seemann "
            "(attribution is mandatory and comes first, 13.3)"
        )

    if relpath in TRANSLATED_MODULES and not RE_PROVENANCE.search(text):
        problems.append(
            "direct translation: the header must name the upstream file and line "
            "range it was ported from, e.g. 'Ported from bin/mlst:370-395'"
        )

    return problems


def insert_header(path):
    """Insert the canonical three-line header at the top of a file.

    Preserves a shebang and a PEP 263 coding cookie, which must stay on lines
    1 and 2. Returns True when the file was modified.
    """
    with open(path, encoding="utf-8") as fh:
        text = fh.read()
    if RE_SPDX.search(text[:4096]):
        # An SPDX tag is already present; only the prose lines are missing and
        # rewriting them automatically would be guesswork.
        return False
    lines = text.splitlines(True)
    at = 0
    if lines and lines[0].startswith("#!"):
        at = 1
    if len(lines) > at and re.match(r"^#.*coding[:=]", lines[at]):
        at += 1
    lines[at:at] = [HEADER_TEMPLATE]
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("".join(lines))
    return True


def iter_python_files(targets, root):
    """Yield (abspath, repo-relative posix path) for every .py under targets."""
    for target in targets:
        if os.path.isfile(target):
            if target.endswith(".py"):
                yield os.path.abspath(target), _rel(target, root)
            continue
        for dirpath, dirnames, filenames in os.walk(target):
            dirnames[:] = sorted(d for d in dirnames if d not in SKIP_DIRS)
            for name in sorted(filenames):
                if name.endswith(".py"):
                    full = os.path.join(dirpath, name)
                    yield os.path.abspath(full), _rel(full, root)


def _rel(path, root):
    """Repository-relative POSIX path, for stable comparison on Windows."""
    return os.path.relpath(os.path.abspath(path), root).replace(os.sep, "/")


def main(argv=None):
    """Command-line entry point. Returns a process exit code."""
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    parser = argparse.ArgumentParser(
        prog="check_headers.py",
        description="Assert the GPL-2.0-only SPDX header on every source file.",
    )
    parser.add_argument("targets", nargs="*", help="files or directories (default: wmlst scripts packaging tests)")
    parser.add_argument("--fix", action="store_true", help="insert the header where it is absent")
    parser.add_argument("-q", "--quiet", action="store_true")
    ns = parser.parse_args(argv)

    targets = ns.targets or [
        os.path.join(root, d) for d in ("wmlst", "scripts", "packaging", "tests")
    ]
    targets = [t for t in targets if os.path.exists(t)]

    failures = []
    fixed = []
    checked = 0
    for path, relpath in iter_python_files(targets, root):
        checked += 1
        problems = check_file(path, relpath)
        if problems and ns.fix:
            if insert_header(path):
                fixed.append(relpath)
                problems = check_file(path, relpath)
        if problems:
            failures.append((relpath, problems))

    for relpath in fixed:
        sys.stdout.write("check_headers: inserted header into %s\n" % relpath)

    if failures:
        sys.stderr.write("check_headers: FAIL\n")
        for relpath, problems in failures:
            sys.stderr.write("  %s\n" % relpath)
            for problem in problems:
                sys.stderr.write("      - %s\n" % problem)
        sys.stderr.write(
            "\n%d of %d file(s) non-compliant. Run:\n"
            "    python3 scripts/check_headers.py --fix\n"
            "then add the upstream provenance line by hand where it is asked for.\n"
            "See docs/ARCHITECTURE.md 13.3 and the NOTICE file.\n"
            % (len(failures), checked)
        )
        return 1

    if not ns.quiet:
        sys.stdout.write("check_headers: OK - %d file(s) carry the GPL-2.0-only header\n" % checked)
    return 0


if __name__ == "__main__":
    sys.exit(main())
