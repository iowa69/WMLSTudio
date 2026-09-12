# SPDX-License-Identifier: GPL-2.0-only
# Copyright (C) 2025 IOWA-BioTech - Giovanni Lorenzin
# Copyright (C) Torsten Seemann (upstream `mlst`, from which WMLST is ported)
"""`--help` snapshot and structure tests (docs/ARCHITECTURE.md section 13.4).

The help screen is a user-facing contract: it reproduces the six section
headings of `bin/mlst:485-513` in order, with the upstream misspellings intact,
and is snapshot-compared against ``tests/golden/help.txt``.

Runs under ``python3 -m pytest tests/test_cli_help.py`` and as a plain script.
"""

import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
GOLDEN = os.path.join(HERE, "golden", "help.txt")

if REPO not in sys.path:
    sys.path.insert(0, REPO)

from wmlst import branding, cli

#: The six upstream headings, in the order bin/mlst emits them.
UPSTREAM_SECTIONS = ("GENERAL", "SCHEME", "OUTPUT FORMAT", "OUTPUT FILES",
                     "SCORING", "DATABASE")

#: Every option upstream's @Options table declares (bin/mlst:440-473).
UPSTREAM_OPTIONS = (
    "help", "version", "check", "skipcheck", "quiet", "threads", "debug",
    "fofn", "scheme", "info", "list", "longlist", "exclude", "full", "legacy",
    "csv", "label", "nopath", "outfile", "novel", "json", "minid", "mincov",
    "minscore", "blastdb", "datadir",
)

#: The WMLST-only switches this port adds (divergence D14).
WMLST_OPTIONS = ("html", "evidence-tsv", "jobs", "gui")


def _run(*args):
    env = dict(os.environ)
    env["PYTHONPATH"] = REPO + os.pathsep + env.get("PYTHONPATH", "")
    return subprocess.run([sys.executable, "-m", "wmlst", *list(args)],
                          cwd=REPO, env=env, capture_output=True,
                          text=True)


def test_help_matches_golden_snapshot():
    with open(GOLDEN, encoding="utf-8", newline="") as fh:
        expected = fh.read()
    assert cli.usage_text() == expected, (
        "--help drifted from tests/golden/help.txt; regenerate deliberately")


def test_golden_snapshot_has_no_crlf():
    with open(GOLDEN, "rb") as fh:
        assert b"\r" not in fh.read()


def test_six_upstream_sections_appear_in_order():
    text = cli.usage_text()
    positions = []
    for name in UPSTREAM_SECTIONS:
        line = "\n%s\n" % name
        assert line in text, "missing help section %s" % name
        positions.append(text.index(line))
    assert positions == sorted(positions), "help sections are out of order"


def test_wmlst_section_comes_after_the_upstream_ones():
    text = cli.usage_text()
    assert text.index("\nWMLST EXTRAS\n") > text.index("\nDATABASE\n")


def test_every_upstream_option_is_offered():
    text = cli.usage_text()
    for name in UPSTREAM_OPTIONS:
        assert "  --%s" % name in text, "missing option --%s" % name


def test_wmlst_extensions_are_offered():
    text = cli.usage_text()
    for name in WMLST_OPTIONS:
        assert "  --%s" % name in text, "missing option --%s" % name


def test_upstream_misspellings_are_preserved():
    text = cli.usage_text()
    for typo in ("allelles", "Minumum", "%identity"):
        assert typo in text, "lost upstream spelling %r" % typo


def test_synopsis_and_usage_lines_match_upstream():
    text = cli.usage_text()
    assert text.startswith(
        "SYNOPSIS\n  Automatic MLST calling from assembled contigs\nUSAGE\n")
    assert "# list known schemes" in text
    assert "# auto-detect scheme" in text
    assert "# force a scheme" in text


def test_help_carries_the_vendor_branding():
    text = cli.usage_text()
    assert branding.VENDOR in text
    assert branding.AUTHOR in text
    assert branding.HOMEPAGE in text


def test_boolean_defaults_render_like_perl():
    text = cli.usage_text()
    # Perl renders a defined-but-false DEFAULT as " (default '0')", which is
    # truthy, so every negatable flag prints "(default OFF)" (bin/mlst:499).
    assert "--quiet           Quiet - no stderr output (default OFF)" in text
    assert "(default ON)" not in text
    assert "--threads INT" in text and "(default '1')" in text


def test_help_goes_to_stdout_with_exit_zero():
    proc = _run("--help")
    assert proc.returncode == 0
    assert proc.stdout.startswith("SYNOPSIS")
    assert proc.stderr == ""


def test_help_mentions_threads():
    # bats case 3: [[ "$output" =~ "threads" ]]
    assert "threads" in _run("--help").stdout


def test_help_is_reachable_through_dash_h():
    assert _run("-h").stdout == cli.usage_text()


def test_version_prints_wmlst_and_version_only():
    from wmlst.version import __version__

    proc = _run("--version")
    assert proc.returncode == 0
    assert proc.stdout == "wmlst %s\n" % __version__
    assert "2.35.0" not in proc.stdout  # C15: compat version never on this line


def test_dash_dash_h_is_help_not_html():
    """Finding 15: `--h` is unambiguous upstream, so --html must not shadow it.

    Getopt::Long's auto_abbrev resolves `--h` to `--help` in `mlst` 2.35.0; the
    WMLST-only `--html`/`--html-evidence` may not take that prefix away.
    """
    proc = _run("--h")
    assert proc.returncode == 0
    assert proc.stdout == cli.usage_text()


def test_upstream_long_names_holds_exactly_the_upstream_table():
    """Finding 15: the prefix-preference set must not creep to the extras."""
    assert cli.UPSTREAM_LONG_NAMES >= set(UPSTREAM_OPTIONS)
    for name in (*WMLST_OPTIONS, "blast-timeout", "html-evidence",
                 "repair-locus-ids", "tsv-evidence"):
        assert name not in cli.UPSTREAM_LONG_NAMES, name


def _main():
    failures = 0
    for name, fn in sorted(globals().items()):
        if not name.startswith("test_") or not callable(fn):
            continue
        try:
            fn()
        except BaseException as exc:
            failures += 1
            print("FAIL %s: %s" % (name, exc))
        else:
            print("ok   %s" % name)
    print("%d failure(s)" % failures)
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(_main())


# ---------------------------------------------------------------------------
# The maintenance flags README and .github/workflows both invoke
# ---------------------------------------------------------------------------

MAINTENANCE_FLAGS = ("--update-db", "--check-only", "--make-blast-db",
                     "--bootstrap-blast")


def test_maintenance_flags_are_accepted():
    """These four are promised by README and used by CI, so they must exist.

    They were documented and wired into .github/workflows before cli.py grew
    them, which turned every CI job red; this pins the three together.
    """
    help_text = _run("--help").stdout
    for flag in MAINTENANCE_FLAGS:
        assert flag in help_text, "%s is missing from --help" % flag


def test_every_wmlst_flag_used_by_ci_and_readme_exists():
    """Scan the workflows and README for `wmlst --flag` and check each one."""
    import re

    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    sources = [os.path.join(root, "README.md")]
    wf = os.path.join(root, ".github", "workflows")
    if os.path.isdir(wf):
        sources += [os.path.join(wf, n) for n in sorted(os.listdir(wf))
                    if n.endswith((".yml", ".yaml"))]

    help_text = _run("--help").stdout
    pattern = re.compile(r"wmlst(?:-cli\.exe)?\s+((?:--[a-z][a-z0-9-]*\s*)+)")
    missing = set()
    for path in sources:
        with open(path, encoding="utf-8") as fh:
            text = fh.read()
        for run in pattern.findall(text):
            for flag in run.split():
                if flag.startswith("--") and flag not in help_text:
                    missing.add((os.path.basename(path), flag))
    assert not missing, "documented/invoked but not implemented: %s" % sorted(missing)
