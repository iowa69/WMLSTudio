# SPDX-License-Identifier: GPL-2.0-only
# Copyright (C) 2025-2026 IOWA-Tech - Giovanni Lorenzin
# Copyright (C) Torsten Seemann (upstream `mlst`, from which WMLST is ported)
"""Packaging, licensing and repository-hygiene tests (docs/ARCHITECTURE.md 13.4).

Everything here guards a promise that is invisible at runtime and expensive to
discover later: that the version is the same number in five files, that the
icon really is a multi-size ICO, that the 116 MB allele database is still under
version control, that the derived BLAST index is still not, that the licence is
GPL-2.0-only everywhere it is spelled out, and that every console entry point
points at something real.

Runnable two ways::

    python -m pytest tests/test_packaging.py
    python tests/test_packaging.py
"""

from __future__ import annotations

import ast
import glob
import os
import re
import subprocess
import sys
import zipfile

import pytest

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

PACKAGING = os.path.join(REPO_ROOT, "packaging")
ICON = os.path.join(PACKAGING, "wmlst.ico")

#: The sizes a Windows icon must carry to look right everywhere from the
#: taskbar (16) to the 250%-DPI Alt-Tab switcher (256).
REQUIRED_ICON_SIZES = (16, 24, 32, 48, 64, 128, 256)

#: Section 13.4: the eleven upstream booleans that must accept --no-<flag>.
NEGATABLE_BOOLEANS = (
    "quiet", "full", "legacy", "csv", "nopath", "debug",
    "check", "skipcheck", "info", "list", "longlist",
)

#: Section 13.1: four console scripts plus one GUI script.
CONSOLE_SCRIPTS = {
    "wmlst": "wmlst.cli:main",
    "wmlst-update-db": "wmlst.cli:main_update_db",
    "wmlst-make-blast-db": "wmlst.cli:main_make_blast_db",
    "wmlst-show-seqs": "wmlst.cli:main_show_seqs",
}
GUI_SCRIPTS = {"wmlst-gui": "wmlst.gui:main"}


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
def _read(relpath, binary=False):
    path = os.path.join(REPO_ROOT, relpath)
    mode = "rb" if binary else "r"
    kwargs = {} if binary else {"encoding": "utf-8"}
    with open(path, mode, **kwargs) as handle:
        return handle.read()


def _iss_defines():
    """-> {NAME: value} for every ``#define NAME "value"`` in wmlst.iss."""
    out = {}
    for line in _read(os.path.join("packaging", "wmlst.iss")).splitlines():
        match = re.match(r'^#define\s+(\w+)\s+"(.*)"\s*$', line.strip())
        if match:
            out[match.group(1)] = match.group(2)
    return out


def _iss_expand(line):
    """Substitute ``{#Name}`` preprocessor references, as ISCC would."""
    defines = _iss_defines()
    for _ in range(4):  # defines may nest; four passes is far more than enough
        expanded = re.sub(
            r"\{#(\w+)\}", lambda m: defines.get(m.group(1), m.group(0)), line
        )
        if expanded == line:
            break
        line = expanded
    return line


def _iss_section(name):
    """-> the directive lines of one wmlst.iss section, comments stripped.

    Inno treats a leading ``;`` as a comment, and the [Icons] section carries a
    long explanation of why there is no command-line shortcut. Grepping the raw
    text for "cmd.exe" would find that explanation and fail; only the directive
    lines are assertions about the installer.
    """
    text = _read(os.path.join("packaging", "wmlst.iss"))
    lines = []
    current = None
    pending = ""
    for raw in text.splitlines():
        line = raw.strip()
        match = re.fullmatch(r"\[([A-Za-z]+)\]", line)
        if match:
            current = match.group(1)
            pending = ""
            continue
        if current != name:
            continue
        if not pending and (not line or line.startswith(";")):
            continue
        # Inno continues a directive onto the next line with a trailing "\".
        if line.endswith("\\"):
            pending += line[:-1].rstrip() + " "
            continue
        lines.append(_iss_expand(pending + line))
        pending = ""
    if pending:
        lines.append(_iss_expand(pending))
    assert lines or ("[%s]" % name) in text, "wmlst.iss has no [%s] section" % name
    return lines


def _load_toml():
    try:
        import tomllib
    except ImportError:  # Python 3.9/3.10
        try:
            import tomli as tomllib  # type: ignore[no-redef]
        except ImportError:
            pytest.skip("no TOML parser available (pip install 'wmlst[test]')")
    with open(os.path.join(REPO_ROOT, "pyproject.toml"), "rb") as handle:
        return tomllib.load(handle)


def _git(*args):
    """Run git in the repository, returning stdout, or None when unusable."""
    try:
        out = subprocess.run(
            ["git", *args],
            cwd=REPO_ROOT,
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
            check=False,
        )
    except (OSError, subprocess.SubprocessError):
        return None
    return out.stdout.decode("utf-8", "replace")


def _repo_has_commits():
    """True when HEAD resolves to a real object.

    On a repository with no commits `git rev-parse HEAD` prints the literal
    string "HEAD" on stdout and the error on stderr, so a truthiness test on
    the output is not enough.
    """
    out = (_git("rev-parse", "HEAD") or "").strip()
    return bool(re.fullmatch(r"[0-9a-f]{40}", out))


def _declared_version():
    """Read __version__ out of wmlst/version.py without importing it."""
    tree = ast.parse(_read(os.path.join("wmlst", "version.py")))
    for node in tree.body:
        if isinstance(node, ast.Assign):
            for target in node.targets:
                if getattr(target, "id", None) == "__version__":
                    return node.value.value
    raise AssertionError("wmlst/version.py does not define __version__")


VERSION = _declared_version()


# ---------------------------------------------------------------------------
# 13.4 -- version sync across version.py / CITATION.cff / wmlst.iss
# ---------------------------------------------------------------------------
def test_version_is_a_release_number():
    """__version__ must be a plain dotted release number (section 13.1)."""
    assert re.match(r"^\d+\.\d+\.\d+([.-]?(a|b|rc)\d+)?$", VERSION), VERSION


def test_citation_cff_version_matches():
    """CITATION.cff must carry the same version as wmlst/version.py (13.3)."""
    text = _read("CITATION.cff")
    versions = re.findall(r'^\s*version:\s*"?([^"\n]+)"?\s*$', text, re.M)
    assert versions, "CITATION.cff declares no version"
    for found in versions:
        if found == "2.35.0":
            continue  # the upstream mlst reference entry, deliberately pinned
        assert found == VERSION, "CITATION.cff says %r, version.py says %r" % (found, VERSION)


def test_inno_setup_version_matches():
    """packaging/wmlst.iss #define MyAppVersion must match version.py (13.4)."""
    text = _read(os.path.join("packaging", "wmlst.iss"))
    found = re.search(r'#define\s+MyAppVersion\s+"([^"]+)"', text)
    assert found, "wmlst.iss defines no MyAppVersion"
    assert found.group(1) == VERSION


def test_windows_version_resource_matches():
    """The PE VERSIONINFO resource must match version.py in all four places."""
    text = _read(os.path.join("packaging", "version_info.txt"))
    parts = (*(int(x) for x in VERSION.split(".")[:3]), 0)
    expected_tuple = "(%d, %d, %d, %d)" % parts
    assert "filevers=%s" % expected_tuple in text
    assert "prodvers=%s" % expected_tuple in text
    dotted = "%d.%d.%d.%d" % parts
    assert "'FileVersion', '%s'" % dotted in text
    assert "'ProductVersion', '%s'" % dotted in text


def test_bundled_db_version_matches_the_tree():
    """version.BUNDLED_DB_VERSION must equal db/VERSION.txt (section 6.1)."""
    from wmlst.version import BUNDLED_DB_VERSION

    on_disk = _read(os.path.join("db", "VERSION.txt")).strip()
    assert on_disk == BUNDLED_DB_VERSION, (on_disk, BUNDLED_DB_VERSION)


# ---------------------------------------------------------------------------
# 13.4 -- the icon
# ---------------------------------------------------------------------------
def _make_icon_module():
    import importlib.util

    spec = importlib.util.spec_from_file_location(
        "wmlst_make_icon", os.path.join(PACKAGING, "make_icon.py")
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_icon_exists_and_is_committed():
    """packaging/wmlst.ico is a build input, so it lives in the repository."""
    assert os.path.isfile(ICON), "run: python packaging/make_icon.py"
    assert os.path.getsize(ICON) > 1024


def test_icon_is_a_wellformed_ico_with_every_required_size():
    """The .ico must parse and contain {16,24,32,48,64,128,256} (section 13.4)."""
    module = _make_icon_module()
    sizes = module.verify_ico(ICON)
    assert sizes == tuple(sorted(REQUIRED_ICON_SIZES)), sizes


def test_icon_header_is_an_icon_not_a_cursor():
    """ICONDIR.type must be 1. A 2 there is a cursor and Windows silently
    refuses to use it as an application icon."""
    blob = _read(os.path.join("packaging", "wmlst.ico"), binary=True)
    assert blob[:2] == b"\x00\x00", "ICONDIR.reserved must be zero"
    assert blob[2:4] == b"\x01\x00", "ICONDIR.type must be 1 (icon)"


def test_icon_is_reproducible():
    """Regenerating the icon must yield the committed bytes exactly.

    The generator is pure arithmetic, so a difference means someone edited the
    .ico by hand or changed the artwork without committing the result.
    """
    module = _make_icon_module()
    assert module.build_ico() == _read(os.path.join("packaging", "wmlst.ico"), binary=True)


# ---------------------------------------------------------------------------
# 13.4 -- all eleven boolean flags accept --no-X
# ---------------------------------------------------------------------------
def test_all_eleven_boolean_flags_accept_a_no_prefix():
    """bats 21/63/96: every upstream boolean must have a --no-<flag> form."""
    import argparse

    from wmlst import cli

    parser = cli.build_parser()
    options = set()
    for action in parser._actions:
        options.update(action.option_strings)

    missing = []
    for name in NEGATABLE_BOOLEANS:
        if "--%s" % name not in options:
            missing.append("--%s" % name)
        if "--no-%s" % name not in options:
            missing.append("--no-%s" % name)
    assert not missing, "missing options: %s" % ", ".join(missing)

    negatable = {
        action.dest
        for action in parser._actions
        if isinstance(action, argparse.BooleanOptionalAction)
    }
    assert set(NEGATABLE_BOOLEANS).issubset(negatable), sorted(
        set(NEGATABLE_BOOLEANS) - negatable
    )


# ---------------------------------------------------------------------------
# 13.2 -- repository hygiene
# ---------------------------------------------------------------------------
def test_the_allele_database_is_present_and_the_right_shape():
    """162 scheme directories and 1,108 locus files (section 6.1)."""
    pubmlst = os.path.join(REPO_ROOT, "db", "pubmlst")
    assert os.path.isdir(pubmlst), "db/pubmlst is missing from this checkout"
    # __pycache__ appears here the moment anything imports wmlst_db.pubmlst,
    # and a dot-directory is a backup or a staging area. Neither is a scheme.
    schemes = [
        name
        for name in os.listdir(pubmlst)
        if os.path.isdir(os.path.join(pubmlst, name))
        and not name.startswith((".", "__"))
    ]
    assert len(schemes) == 162, len(schemes)
    tfa = glob.glob(os.path.join(pubmlst, "*", "*.tfa"))
    assert len(tfa) == 1108, len(tfa)


def test_db_is_tracked():
    """db/pubmlst must never be gitignored (section 13.2).

    116 MB raw becomes ~5.8 MB of git objects, and it is the corpus the golden
    outputs were generated against. Dropping it from version control silently
    changes every result the next contributor gets.
    """
    ignored = _git("check-ignore", "-v", "db/pubmlst/abaumannii/abaumannii.txt")
    if ignored is None:
        pytest.skip("git is not available")
    assert ignored.strip() == "", "db/pubmlst is gitignored: %s" % ignored.strip()

    tracked = _git("ls-files", "db/pubmlst")
    if tracked is None:
        pytest.skip("git is not available")
    count = len([line for line in tracked.splitlines() if line.strip()])
    if count == 0 and not _repo_has_commits():
        pytest.skip("repository has no commits yet; nothing is tracked at all")
    assert count >= 1108, "only %d database files are tracked" % count


def test_derived_blast_index_is_gitignored():
    """db/blast is 180 MB of rebuildable index and must stay out of git."""
    ignored = _git("check-ignore", "-v", "db/blast/mlst.fa")
    if ignored is None:
        pytest.skip("git is not available")
    assert "db/blast" in ignored, "db/blast is NOT ignored: %r" % ignored


def test_gitignore_covers_the_usual_build_noise():
    text = _read(".gitignore")
    for needle in ("db/blast/", "__pycache__/", "dist/", "build/", "*.egg-info/", ".pytest_cache/"):
        assert needle in text, "missing from .gitignore: %s" % needle
    rules = [
        line.strip()
        for line in text.splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ]
    # `db/pubmlst/*/database_version.txt` is legitimate: that file is derived and
    # only written with --write-version-files. A rule that swallows the alleles
    # themselves is not.
    for rule in rules:
        assert rule.rstrip("/") not in ("db/pubmlst", "db", "/db/pubmlst", "/db"), (
            "the allele database must never be gitignored: %r" % rule
        )
        assert not rule.startswith("db/pubmlst/*.") , rule
        assert "*.tfa" not in rule, "an allele file pattern is gitignored: %r" % rule


def test_gitattributes_freezes_the_byte_identical_paths():
    """core.autocrlf=true must not be able to rewrite an allele file (13.2).

    Changing one line ending changes that allele's MD5, which changes every
    --novel identifier and every BLAST hit it produces.
    """
    text = _read(".gitattributes")
    for path in ("db/pubmlst/**", "tests/data/**", "tests/golden/**"):
        pattern = re.compile(r"^%s\s+.*-text" % re.escape(path), re.M)
        assert pattern.search(text), "%s is not marked -text" % path
    for ext in ("*.tfa", "*.ico", "*.gz", "*.bz2", "*.zip"):
        assert re.search(r"^%s\s+binary" % re.escape(ext), text, re.M), ext


def test_git_reports_the_frozen_paths_as_not_text():
    out = _git("check-attr", "text", "eol", "--", "db/pubmlst/abaumannii/abaumannii.txt")
    if out is None:
        pytest.skip("git is not available")
    assert "text: unset" in out, out


# ---------------------------------------------------------------------------
# 13.3 -- licensing
# ---------------------------------------------------------------------------
def test_license_is_gpl_v2_only():
    """LICENSE must be the GPLv2 text and must not be GPLv3 (section 13.3)."""
    text = _read("LICENSE")
    assert "GNU GENERAL PUBLIC LICENSE" in text
    assert "Version 2, June 1991" in text
    assert "Version 3" not in text
    assert "TERMS AND CONDITIONS FOR COPYING, DISTRIBUTION AND MODIFICATION" in text


def test_spdx_identifier_is_consistent_everywhere():
    """GPL-2.0-only in pyproject, CITATION.cff and the installer (13.3)."""
    data = _load_toml()
    assert data["project"]["license"] == "GPL-2.0-only", data["project"]["license"]
    assert "license: GPL-2.0-only" in _read("CITATION.cff")
    assert "GPL-2.0-only" in _read(os.path.join("packaging", "wmlst.iss"))


def test_license_files_ship_with_the_distribution():
    """GPLv2 section 1 requires the licence to travel with every copy."""
    data = _load_toml()
    assert "LICENSE" in data["project"]["license-files"]
    manifest = _read("MANIFEST.in")
    assert "include LICENSE" in manifest
    assert "include NOTICE" in manifest
    iss = _read(os.path.join("packaging", "wmlst.iss"))
    assert "LicenseFile=" in iss
    assert 'DestName: "LICENSE.txt"' in iss, "the installer must put LICENSE in {app}"
    assert 'DestName: "NOTICE.txt"' in iss


def test_notice_credits_upstream_pubmlst_and_ncbi():
    """Attribution is mandatory and upstream comes first (section 13.3)."""
    text = _read("NOTICE")
    for needle in (
        "Torsten Seemann",
        "https://github.com/tseemann/mlst",
        "Jolley",
        "30345391",
        "Wellcome Open Research",
        "PubMLST",
        "BLAST+",
        "public domain",
        "IOWA-Tech",
        "Giovanni Lorenzin",
        "GPL-2.0-only",
    ):
        assert needle in text, "NOTICE does not mention %r" % needle


def test_readme_carries_the_attribution_and_the_branding():
    """Acceptance check 12 and added bats case 51."""
    text = _read("README.md")
    for needle in (
        "IOWA-Tech",
        "Giovanni Lorenzin",
        "Torsten Seemann",
        "tseemann/mlst",
        "GPL-2.0-only",
        "30345391",
        "SmartScreen",
    ):
        assert needle in text, "README.md does not mention %r" % needle


def test_readme_does_not_claim_the_build_is_signed():
    """Section 13.1: be honest about code signing, in every document."""
    text = _read("README.md").lower()
    assert "not code-signed" in text or "not code signed" in text
    for lie in ("digitally signed by", "signed with an authenticode", "verified publisher"):
        assert lie not in text, "README claims a signature it does not have: %r" % lie


def test_packaging_sources_carry_the_spdx_header():
    """Every file this module owns carries the GPL-2.0-only header (13.3).

    Scoped to the files this module owns, because the whole-repository run is
    the CI gate (.github/workflows/ci.yml), where a missing header in someone
    else's module is that owner's failure to fix rather than a reason for the
    packaging suite to go red.
    """
    result = subprocess.run(
        [sys.executable, os.path.join("scripts", "check_headers.py"),
         "scripts", "packaging",
         os.path.join("tests", "conftest.py"),
         os.path.join("tests", "test_packaging.py")],
        cwd=REPO_ROOT,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )
    assert result.returncode == 0, result.stdout.decode("utf-8", "replace")


def test_check_headers_is_wired_into_ci_for_the_whole_repository():
    """The repository-wide enforcement lives in CI (section 13.3)."""
    assert "check_headers.py" in _read(os.path.join(".github", "workflows", "ci.yml"))
    assert "check_headers.py" in _read(os.path.join(".github", "workflows", "release.yml"))


# ---------------------------------------------------------------------------
# 13.1 -- dependencies and entry points
# ---------------------------------------------------------------------------
def test_there_are_no_mandatory_runtime_dependencies():
    """A GPLv2 work cannot link Apache-2.0 code; zero deps is the guard (13.3)."""
    data = _load_toml()
    assert data["project"]["dependencies"] == [], data["project"]["dependencies"]


def test_stdlib_only_check_passes():
    """scripts/check_stdlib_only.py must be green (acceptance check 6)."""
    result = subprocess.run(
        [sys.executable, os.path.join("scripts", "check_stdlib_only.py"), "wmlst"],
        cwd=REPO_ROOT,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )
    assert result.returncode == 0, result.stdout.decode("utf-8", "replace")


def test_python_floor_is_three_nine():
    data = _load_toml()
    assert data["project"]["requires-python"] == ">=3.9"


def test_entry_points_are_declared():
    """Four console scripts and one GUI script (section 13.1)."""
    data = _load_toml()
    assert data["project"]["scripts"] == CONSOLE_SCRIPTS
    assert data["project"]["gui-scripts"] == GUI_SCRIPTS


def test_the_gui_is_a_gui_script_not_a_console_script():
    """setuptools binds [project.scripts] to python.exe and [project.gui-scripts]
    to pythonw.exe. WMLST.exe must never flash a console window, so wmlst-gui
    has to be in the second table (section 13.1)."""
    data = _load_toml()
    assert "wmlst-gui" not in data["project"]["scripts"]
    assert "wmlst-gui" in data["project"]["gui-scripts"]


@pytest.mark.parametrize(
    "name,target", sorted(dict(CONSOLE_SCRIPTS, **GUI_SCRIPTS).items())
)
def test_entry_point_target_exists(name, target):
    """Every console script must resolve to a real callable.

    A declared entry point that does not resolve produces an exe that fails with
    an ImportError the first time a user runs it, which is the worst possible
    place to discover it. Set WMLST_STRICT_PACKAGING=1 to make an unimplemented
    target a hard failure (the release workflow does).
    """
    import importlib

    module_name, _, attr = target.partition(":")
    module = importlib.import_module(module_name)
    func = getattr(module, attr, None)
    if func is None:
        message = (
            "%s -> %s does not exist yet. wmlst/cli.py must define %s() "
            "(docs/ARCHITECTURE.md 13.1 names five entry points)." % (name, target, attr)
        )
        if os.environ.get("WMLST_STRICT_PACKAGING") == "1":
            pytest.fail(message)
        pytest.skip(message)
    assert callable(func)


# ---------------------------------------------------------------------------
# 13.1 -- the wheel carries the database
# ---------------------------------------------------------------------------
def test_package_data_maps_the_database_into_the_wheel():
    data = _load_toml()
    setuptools = data["tool"]["setuptools"]
    assert setuptools["package-dir"]["wmlst_db"] == "db"
    assert "wmlst_db" in setuptools["packages"]
    patterns = setuptools["package-data"]["wmlst_db"]
    for pattern in ("VERSION.txt", "pubmlst/*/*.tfa", "pubmlst/*/*.txt", "pubmlst/*/*.json"):
        assert pattern in patterns, pattern
    # db/pubmlst must NOT be a declared package: importing it writes a
    # __pycache__ directory that every scheme enumerator would read as a 163rd
    # scheme and emit in --list, --longlist and --info.
    assert "wmlst_db.pubmlst" not in setuptools["packages"]


def test_manifest_ships_the_database_and_prunes_the_index():
    text = _read("MANIFEST.in")
    assert "recursive-include db/pubmlst *.tfa" in text
    assert "prune db/blast" in text


@pytest.mark.slow
def test_built_wheel_contains_the_alleles_and_not_the_index(tmp_path):
    """Build a real wheel and look inside it (section 13.1).

    Marked slow: it copies 116 MB. Skipped when `build` is not installed.
    """
    # Guard on build.__main__, not on `build`: a local `build/` directory (pip
    # leaves one behind after an editable install, and it is gitignored) is
    # importable as a namespace package, so importorskip("build") succeeds and
    # the run then dies on "'build' is a package and cannot be directly executed".
    pytest.importorskip("build.__main__", reason="pip install 'wmlst[build]'")
    out = tmp_path / "dist"
    result = subprocess.run(
        [sys.executable, "-m", "build", "--wheel", "--outdir", str(out), REPO_ROOT],
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )
    assert result.returncode == 0, result.stdout.decode("utf-8", "replace")[-4000:]
    wheels = list(out.glob("*.whl"))
    assert len(wheels) == 1, wheels
    names = zipfile.ZipFile(str(wheels[0])).namelist()
    assert len([n for n in names if n.endswith(".tfa")]) == 1108
    assert not [n for n in names if n.startswith("wmlst_db/blast/")]
    assert "wmlst_db/VERSION.txt" in names
    assert "wmlst_db/pubmlst/__init__.py" in names


# ---------------------------------------------------------------------------
# 13.1 -- the PyInstaller spec and the installer
# ---------------------------------------------------------------------------
def test_pyinstaller_spec_is_valid_python():
    ast.parse(_read(os.path.join("packaging", "wmlst.spec")))


def test_pyinstaller_spec_builds_both_executables_and_merges_them():
    text = _read(os.path.join("packaging", "wmlst.spec"))
    assert 'name="WMLST"' in text and "console=False" in text
    assert 'name="wmlst-cli"' in text and "console=True" in text
    assert "MERGE(" in text, "without MERGE the 116 MB database is stored twice"
    assert "COLLECT(" in text, "one-folder only; one-file unpacks the DB on every launch"
    assert "upx=False" in text, "UPX corrupts the MSVC runtime DLLs"
    assert "runtime_hook_noconsole.py" in text
    assert "wmlst.ico" in text


def test_pyinstaller_spec_excludes_the_heavy_libraries():
    text = _read(os.path.join("packaging", "wmlst.spec"))
    for name in ("numpy", "scipy", "pandas", "matplotlib", "PIL", "pytest", "setuptools"):
        assert '"%s"' % name in text, "%s is not in the spec's excludes" % name


def test_pyinstaller_spec_does_not_bundle_the_derived_index():
    text = _read(os.path.join("packaging", "wmlst.spec"))
    assert 'excludes=["blast"' in text.replace("'", '"')


def test_inno_setup_is_configured_the_way_13_1_requires():
    text = _read(os.path.join("packaging", "wmlst.iss"))
    assert "Compression=lzma2/max" in text
    assert "SolidCompression=yes" in text
    assert "ArchitecturesAllowed=x64compatible" in text
    assert "MinVersion=10.0.17763" in text
    assert "PrivilegesRequired=lowest" in text
    assert "IOWA-Tech" in text and "Giovanni Lorenzin" in text
    assert "[UninstallDelete]" in text


def test_inno_setup_creates_a_start_menu_and_a_default_desktop_shortcut():
    """Start Menu always; desktop shortcut offered and ticked by default.

    Double-clicking a desktop icon is how most people will start WMLST, so the
    task is pre-selected -- but it stays a [Tasks] entry the user can untick,
    never an unconditional [Icons] line they cannot refuse.
    """
    text = _read(os.path.join("packaging", "wmlst.iss"))
    assert "{group}\\{#MyAppName}" in text, "no Start Menu entry"
    assert "{autodesktop}\\{#MyAppName}" in text, "no desktop shortcut"
    task = re.search(r'^Name:\s*"desktopicon".*$', text, re.M)
    assert task, "the desktop shortcut must remain a deselectable [Tasks] entry"
    assert "unchecked" not in task.group(0), (
        "the desktop shortcut should be ticked by default"
    )
    assert re.search(r'Tasks:\s*desktopicon', text), (
        "the desktop [Icons] line must be gated on the desktopicon task"
    )


# ---------------------------------------------------------------------------
# C21 -- drag-and-drop is optional at runtime, mandatory in a frozen build
# ---------------------------------------------------------------------------
# wmlst/gui.py imports tkinterdnd2 inside a try/except and degrades silently to
# click-to-browse (ARCHITECTURE.md C21). Silent is right for `pip install wmlst`
# and wrong for a release: dropping a FASTA on the window is the primary
# interaction, and a build that lost it looks exactly like one that did not.
# The spec therefore refuses to build without it, and these tests execute the
# spec -- with PyInstaller stubbed out -- to prove the refusal actually happens.
# ---------------------------------------------------------------------------

#: What collect_all("tkinterdnd2") returns on a healthy Windows build machine:
#: the two Python modules, the tkdnd Tcl scripts, and the native extension that
#: TkinterDnD._require() loads out of <package>/tkdnd/<platform>/.
FAKE_DND_OK = (
    [
        ("/site/tkinterdnd2/tkdnd/win-x64/tkdnd.tcl", "tkinterdnd2/tkdnd/win-x64"),
        ("/site/tkinterdnd2/tkdnd/win-x64/pkgIndex.tcl", "tkinterdnd2/tkdnd/win-x64"),
    ],
    [("/site/tkinterdnd2/tkdnd/win-x64/libtkdnd2.10.2.dll", "tkinterdnd2/tkdnd/win-x64")],
    ["tkinterdnd2", "tkinterdnd2.TkinterDnD"],
)

#: The package imports but its Tcl payload was not collected. The frozen GUI
#: would import tkinterdnd2 happily and then die on `package require tkdnd`.
FAKE_DND_NO_TCL = (
    [("/site/tkinterdnd2/__init__.py", "tkinterdnd2")],
    [],
    ["tkinterdnd2", "tkinterdnd2.TkinterDnD"],
)

#: tkinterdnd2 is not installed at all.
FAKE_DND_ABSENT = ([], [], [])


def _stub_pyinstaller(monkeypatch, collect_all_result, record):
    """Put a minimal fake PyInstaller in sys.modules.

    The spec is executed for real, so every name it imports has to exist; none
    of them has to build anything. `record` collects the EXE()/COLLECT()/...
    calls so the caller can inspect what the spec asked for.
    """
    import types

    def maker(kind):
        def factory(*args, **kwargs):
            record.append((kind, args, kwargs))
            stub = types.SimpleNamespace(**kwargs)
            for attr in ("scripts", "pure", "zipped_data", "binaries", "datas"):
                setattr(stub, attr, [])
            return stub

        return factory

    root = types.ModuleType("PyInstaller")
    building = types.ModuleType("PyInstaller.building")
    api = types.ModuleType("PyInstaller.building.api")
    build_main = types.ModuleType("PyInstaller.building.build_main")
    datastruct = types.ModuleType("PyInstaller.building.datastruct")
    utils = types.ModuleType("PyInstaller.utils")
    hooks = types.ModuleType("PyInstaller.utils.hooks")

    for name in ("COLLECT", "EXE", "MERGE", "PYZ"):
        setattr(api, name, maker(name))
    build_main.Analysis = maker("Analysis")
    datastruct.Tree = maker("Tree")
    if collect_all_result is not None:
        hooks.collect_all = lambda package: collect_all_result

    modules = {
        "PyInstaller": root,
        "PyInstaller.building": building,
        "PyInstaller.building.api": api,
        "PyInstaller.building.build_main": build_main,
        "PyInstaller.building.datastruct": datastruct,
        "PyInstaller.utils": utils,
        "PyInstaller.utils.hooks": hooks,
    }
    for name, module in modules.items():
        monkeypatch.setitem(sys.modules, name, module)


def _exec_spec(tmp_path, monkeypatch, collect_all_result, *, allow_no_dnd=None,
               platform="win32"):
    """Run packaging/wmlst.spec end to end against the stubs.

    The spec is copied into `tmp_path` first: it writes generated entry-point
    launchers next to itself, and a test must not leave files in the checkout.
    -> (spec globals, recorded calls). Raises SystemExit when the spec refuses.
    """
    import shutil

    staged = tmp_path / "packaging"
    staged.mkdir(parents=True)
    for name in (
        "wmlst.spec", "wmlst.ico", "version_info.txt", "runtime_hook_noconsole.py",
    ):
        shutil.copy2(os.path.join(PACKAGING, name), staged / name)

    record = []
    _stub_pyinstaller(monkeypatch, collect_all_result, record)
    monkeypatch.setattr(sys, "platform", platform, raising=False)
    if allow_no_dnd is None:
        monkeypatch.delenv("WMLST_ALLOW_NO_DND", raising=False)
    else:
        monkeypatch.setenv("WMLST_ALLOW_NO_DND", allow_no_dnd)

    source = (staged / "wmlst.spec").read_text(encoding="utf-8")
    namespace = {"__name__": "__main__", "SPECPATH": str(staged)}
    exec(compile(source, "wmlst.spec", "exec"), namespace)
    return namespace, record


def test_spec_bundles_tkinterdnd2_when_it_is_available(tmp_path, monkeypatch):
    """The happy path: the Tcl payload and the hidden imports reach the bundle."""
    namespace, _record = _exec_spec(tmp_path, monkeypatch, FAKE_DND_OK)

    assert "tkinterdnd2" in namespace["HIDDEN"]
    assert "tkinterdnd2.TkinterDnD" in namespace["HIDDEN"]

    dests = {dest for _src, dest in namespace["DATAS"]}
    assert "tkinterdnd2/tkdnd/win-x64" in dests, (
        "the tkdnd Tcl scripts must be collected, not just the Python shim"
    )
    assert namespace["BINARIES"], "the native tkdnd extension was not collected"
    assert any(
        "libtkdnd" in os.path.basename(src) for src, _dest in namespace["BINARIES"]
    )


@pytest.mark.parametrize(
    ("payload", "because"),
    [
        (FAKE_DND_ABSENT, "tkinterdnd2 is not installed"),
        (FAKE_DND_NO_TCL, "the tkdnd Tcl extension was not collected"),
        (None, "PyInstaller.utils.hooks.collect_all is unavailable"),
    ],
)
def test_spec_refuses_to_build_without_drag_and_drop(
    tmp_path, monkeypatch, payload, because
):
    """Fail loudly, not quietly.

    A frozen WMLST without drag-and-drop is a regression the build log must
    scream about, because nothing downstream of the build can detect it.
    """
    with pytest.raises(SystemExit) as excinfo:
        _exec_spec(tmp_path, monkeypatch, payload)

    message = str(excinfo.value)
    assert "tkinterdnd2" in message, because
    # The message has to tell the person reading a red CI log what to type.
    assert "pip install" in message
    assert "WMLST_ALLOW_NO_DND" in message, (
        "the refusal must document its own escape hatch"
    )


def test_spec_refuses_when_tkdnd_has_no_windows_build(tmp_path, monkeypatch):
    """A Linux-only tkdnd payload is useless to a Windows release."""
    linux_only = (
        [("/site/tkinterdnd2/tkdnd/linux-x64/tkdnd.tcl", "tkinterdnd2/tkdnd/linux-x64")],
        [("/site/tkinterdnd2/tkdnd/linux-x64/libtkdnd.so", "tkinterdnd2/tkdnd/linux-x64")],
        ["tkinterdnd2"],
    )
    with pytest.raises(SystemExit) as excinfo:
        _exec_spec(tmp_path, monkeypatch, linux_only, platform="win32")
    assert "win-x64" in str(excinfo.value)

    # ... and the same payload is fine when the build host is not Windows.
    namespace, _ = _exec_spec(tmp_path / "linux", monkeypatch, linux_only,
                              platform="linux")
    assert "tkinterdnd2" in namespace["HIDDEN"]


def test_spec_escape_hatch_builds_without_drag_and_drop(tmp_path, monkeypatch):
    """WMLST_ALLOW_NO_DND=1 is a deliberate, explicit downgrade -- and it works.

    It exists so a developer without the extra can still produce a build. It is
    never set in CI; test_workflows_install_the_gui_extra pins that.
    """
    namespace, record = _exec_spec(
        tmp_path, monkeypatch, FAKE_DND_ABSENT, allow_no_dnd="1"
    )
    assert "tkinterdnd2" not in namespace["HIDDEN"]
    assert namespace["BINARIES"] == []
    assert [kind for kind, _a, _k in record].count("EXE") == 2


def test_workflows_install_the_gui_extra_wherever_they_freeze_the_app():
    """Every workflow job that runs PyInstaller must install tkinterdnd2.

    The spec now fails without it, so a workflow that forgot the extra would be
    a broken release pipeline rather than a silently worse product -- but the
    cheapest place to notice is here.
    """
    for relpath in (
        os.path.join(".github", "workflows", "ci.yml"),
        os.path.join(".github", "workflows", "release.yml"),
    ):
        text = _read(relpath)
        if "pyinstaller --noconfirm" not in text:
            continue
        assert re.search(r'pip install -e "\.\[[^"]*\bgui\b[^"]*\]"', text), (
            "%s freezes the app but never installs the gui extra" % relpath
        )
        assert "WMLST_ALLOW_NO_DND" not in text, (
            "%s must not disable the drag-and-drop guard" % relpath
        )


def test_gui_extra_provides_tkinterdnd2():
    """The `gui` extra is what the spec's error message tells people to install."""
    extras = _load_toml()["project"]["optional-dependencies"]
    assert any(
        req.split(">")[0].split("=")[0].split("[")[0].strip() == "tkinterdnd2"
        for req in extras["gui"]
    ), extras["gui"]


# ---------------------------------------------------------------------------
# 13.1 -- only ONE thing in this product is meant to be clicked
# ---------------------------------------------------------------------------
# WMLST.exe (windowed) and wmlst-cli.exe (console) live in the same folder. For
# a novice, two files wearing the same artwork is a coin flip, and losing it
# means a black console window that appears and vanishes. The rule that follows
# from that: the branded icon marks the GUI and nothing else, and the Start Menu
# offers the GUI and nothing else. wmlst-cli.exe stays installed and stays
# usable from a terminal -- it is simply never advertised as clickable.
# ---------------------------------------------------------------------------
def _spec_exe_calls():
    """-> {exe name: {keyword: unparsed source}} for every EXE() in wmlst.spec.

    Parsed rather than grepped: ``icon=None`` and ``icon = None`` and a keyword
    inherited from ``**_exe_common`` are three different facts, and only an AST
    tells them apart.
    """
    tree = ast.parse(_read(os.path.join("packaging", "wmlst.spec")))
    found = {}
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        if not (isinstance(node.func, ast.Name) and node.func.id == "EXE"):
            continue
        kwargs = {}
        starred = []
        for keyword in node.keywords:
            if keyword.arg is None:
                starred.append(ast.unparse(keyword.value))
            else:
                kwargs[keyword.arg] = ast.unparse(keyword.value)
        kwargs["**"] = starred
        found[ast.literal_eval(kwargs["name"])] = kwargs
    return found


def _spec_exe_common_keys():
    """-> the keyword names of the shared ``_exe_common = dict(...)`` literal."""
    tree = ast.parse(_read(os.path.join("packaging", "wmlst.spec")))
    for node in ast.walk(tree):
        if not isinstance(node, ast.Assign):
            continue
        targets = [t.id for t in node.targets if isinstance(t, ast.Name)]
        if "_exe_common" not in targets:
            continue
        if isinstance(node.value, ast.Call) and isinstance(node.value.func, ast.Name):
            return {k.arg for k in node.value.keywords if k.arg}
    raise AssertionError("wmlst.spec no longer defines _exe_common = dict(...)")


def test_only_the_windowed_exe_carries_the_branded_icon():
    """WMLST.exe gets wmlst.ico; wmlst-cli.exe gets no custom icon at all.

    The console executable keeps the default console icon on purpose, so that
    Explorer, the taskbar and Alt-Tab all show it as what it is. If this test
    fails because someone put ``icon=ICON`` back on the CLI, the fix is to take
    it off again, not to update the test.
    """
    exes = _spec_exe_calls()
    assert set(exes) == {"WMLST", "wmlst-cli"}, sorted(exes)

    assert "icon" not in _spec_exe_common_keys(), (
        "_exe_common must not set `icon`: sharing it gives both executables the "
        "same artwork, which is exactly what this change removed"
    )

    assert exes["WMLST"]["icon"] == "ICON"
    assert exes["WMLST"]["console"] == "False"

    assert "icon" in exes["wmlst-cli"], (
        "wmlst-cli must pass icon explicitly; omitting the keyword would silently "
        "inherit one again the next time _exe_common changes"
    )
    assert exes["wmlst-cli"]["icon"] == "None"
    assert exes["wmlst-cli"]["console"] == "True"


def test_inno_setup_has_no_start_menu_shortcut_for_the_command_line():
    """No 'WMLST command line' entry, and nothing in [Icons] launches a shell."""
    icons = _iss_section("Icons")
    for line in icons:
        lowered = line.lower()
        assert "cmd.exe" not in lowered, (
            "an [Icons] entry opens a console window: %s" % line
        )
        assert "powershell" not in lowered, line
        assert "wmlst-cli" not in lowered, (
            "the command-line tool must not have a Start Menu shortcut: %s" % line
        )
        assert "command line" not in lowered, line

    filenames = [
        re.search(r'Filename:\s*"([^"]+)"', line).group(1)
        for line in icons
        if re.search(r'Filename:\s*"([^"]+)"', line)
    ]
    assert filenames, "[Icons] has no entries at all"
    exe_targets = [f for f in filenames if f.lower().endswith(".exe")]
    assert exe_targets, "[Icons] points at no executable"
    for target in exe_targets:
        assert os.path.basename(target.replace("\\", "/")) == "WMLST.exe", (
            "the only clickable executable is WMLST.exe: %s" % target
        )


def test_inno_setup_still_installs_the_cli_and_still_offers_the_path_task():
    """Removing the shortcut must not remove the tool or the way to reach it."""
    text = _read(os.path.join("packaging", "wmlst.iss"))
    assert '#define MyAppCliName     "wmlst-cli.exe"' in text, (
        "wmlst-cli.exe is still part of the product"
    )

    # The whole dist\WMLST tree is copied verbatim, and wmlst-cli.exe is in it.
    files = _iss_section("Files")
    assert any(
        "dist\\WMLST\\*" in line and "recursesubdirs" in line for line in files
    ), "the installer no longer copies the PyInstaller one-folder tree"

    tasks = _iss_section("Tasks")
    addtopath = [line for line in tasks if 'Name: "addtopath"' in line]
    assert addtopath, "the optional add-to-PATH task was removed"
    assert "unchecked" in addtopath[0], (
        "add-to-PATH must stay opt-in; it edits the user's environment"
    )
    assert "wmlst-cli.exe" in addtopath[0], (
        "the task must still name the tool it adds to PATH: %s" % addtopath[0]
    )

    registry = _iss_section("Registry")
    path_rows = [line for line in registry if "Tasks: addtopath" in line]
    assert path_rows, "no [Registry] row gated on the addtopath task"
    assert any("Environment" in line and "Path" in line for line in path_rows)


def test_uninstall_deletes_only_derived_data():
    """Never delete a scientist's results on uninstall (section 13.1)."""
    block = _read(os.path.join("packaging", "wmlst.iss")).split("[UninstallDelete]")[1]
    block = block.split("[Code]")[0]
    # Strip Inno comments (`;` to end of line) so the explanation of what must
    # never be deleted is not mistaken for a deletion rule.
    block = "\n".join(
        line for line in block.splitlines() if not line.lstrip().startswith(";")
    )
    assert "db\\blast" in block
    for forbidden in ("{userdocs}", "{userdesktop}", "{userprofile}", "results"):
        assert forbidden not in block, "uninstall would remove user data: %s" % forbidden


def test_smartscreen_is_documented_honestly():
    """Section 13.1: say plainly that it is unsigned, and do not suggest a fix
    that is really a security downgrade."""
    for relpath in (
        os.path.join("packaging", "wmlst.iss"),
        os.path.join("packaging", "SMARTSCREEN.txt"),
        os.path.join("docs", "WINDOWS.md"),
        "README.md",
    ):
        text = _read(relpath)
        assert "SmartScreen" in text, relpath
        assert "Run anyway" in text, relpath
        lowered = text.lower()
        assert "not code-signed" in lowered or "not code signed" in lowered, relpath
        # Telling a user to switch SmartScreen off is a security downgrade
        # dressed up as support. If the phrase appears at all it must appear in
        # a sentence promising we will never ask for it.
        for phrase in ("disable smartscreen", "turn off smartscreen", "disable defender"):
            if phrase in lowered:
                assert "never" in lowered, "%s tells the user to %s" % (relpath, phrase)


# ---------------------------------------------------------------------------
# CI
# ---------------------------------------------------------------------------
def _workflow_files():
    return sorted(glob.glob(os.path.join(REPO_ROOT, ".github", "workflows", "*.yml")))


def test_workflows_exist():
    names = [os.path.basename(p) for p in _workflow_files()]
    assert "ci.yml" in names
    assert "release.yml" in names


def test_workflow_yaml_parses():
    yaml = pytest.importorskip("yaml", reason="pip install 'wmlst[test]'")
    for path in _workflow_files():
        with open(path, encoding="utf-8") as handle:
            document = yaml.safe_load(handle)
        assert "jobs" in document, path
        for name, job in document["jobs"].items():
            assert "runs-on" in job, "%s: job %s has no runs-on" % (path, name)


def test_workflows_cover_ubuntu_and_windows():
    text = "".join(_read(os.path.relpath(p, REPO_ROOT)) for p in _workflow_files())
    assert "ubuntu-latest" in text
    assert "windows-latest" in text
    assert "macos" not in text, "macOS runners are billed 10x; section 13.4 says no"


def test_every_action_is_pinned_to_a_major_tag():
    """An unpinned `uses:` is a supply-chain hole and a reproducibility hole."""
    unpinned = []
    for path in _workflow_files():
        for lineno, line in enumerate(_read(os.path.relpath(path, REPO_ROOT)).splitlines(), 1):
            match = re.search(r"^\s*-?\s*uses:\s*(\S+)", line)
            if match and "@" not in match.group(1):
                unpinned.append("%s:%d %s" % (os.path.basename(path), lineno, match.group(1)))
    assert not unpinned, unpinned


def test_ci_runs_the_licence_and_dependency_gates():
    text = _read(os.path.join(".github", "workflows", "ci.yml"))
    assert "check_stdlib_only.py" in text
    assert "check_headers.py" in text
    assert "make_icon.py --verify" in text


def test_release_workflow_gates_on_the_test_suite_and_writes_hashes():
    text = _read(os.path.join(".github", "workflows", "release.yml"))
    assert "pytest" in text
    assert "SHA256SUMS.txt" in text
    assert "pyinstaller" in text
    assert "wmlst.iss" in text
    assert "twine upload" not in text, "PyPI upload stays a manual step (section 13.4)"


# ---------------------------------------------------------------------------
# stand-alone runner
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    raise SystemExit(pytest.main([os.path.abspath(__file__), "-v"]))
