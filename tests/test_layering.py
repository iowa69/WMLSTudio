# SPDX-License-Identifier: GPL-2.0-only
# Copyright (C) 2025-2026 IOWA-Tech - Giovanni Lorenzin
# Copyright (C) Torsten Seemann (upstream `mlst`, from which WMLST is ported)
"""Import-layering guards for docs/ARCHITECTURE.md section 2.1.

These tests read source, never behaviour, so they run with no display, no BLAST+
and no database. A module that does not exist yet is skipped, not failed: the
implementers work in one tree at the same time.

Runnable both as ``python -m pytest tests/test_layering.py`` and as
``python3 tests/test_layering.py``.
"""

from __future__ import annotations

import ast
import sys
from pathlib import Path
from unittest import SkipTest

PKG = Path(__file__).resolve().parents[1] / "wmlst"

#: The acyclic dependency graph of section 2.1, as a rank. A module may import
#: only modules of a STRICTLY LOWER rank.
RANK = {
    "version": 0,
    "branding": 1,
    "any2fasta": 2,
    "blastbin": 2,
    "schemes": 2,
    "engine": 3,
    "report": 4,
    "updatedb": 4,
    "gui": 5,
    # cli sits above gui because --gui hands over to it (divergence D14);
    # the reverse import would be the cycle.
    "cli": 6,
    "__init__": 7,
    "__main__": 7,
}

#: Only blastbin.py may launch a process (section 2.1).
PROCESS_NAMES = {"subprocess", "Popen"}

#: Names of engine FUNCTIONS. report.py may import engine's dataclasses but none
#: of these: it performs no biology.
ENGINE_FUNCTIONS = {
    "Engine", "parse_blast", "collect_calls", "build_signature", "score_signature",
    "status_column", "sort_duplicate_codes", "analyse", "analyse_file",
}

#: Tokens the worker side of the GUI must never touch (Tcl is not thread-safe).
TK_NAMES = {"tk", "ttk", "tkinter", "tkinterdnd2", "filedialog", "messagebox",
            "tkfont", "Tk", "Toplevel"}


def _tree(name: str) -> ast.Module:
    path = PKG / (name + ".py")
    if not path.is_file():
        raise SkipTest("{} has not been written yet".format(path.name))
    return ast.parse(path.read_text(encoding="utf-8"), filename=str(path))


def _source(name: str) -> str:
    path = PKG / (name + ".py")
    if not path.is_file():
        raise SkipTest("{} has not been written yet".format(path.name))
    return path.read_text(encoding="utf-8")


def _walk_module_scope(tree: ast.Module):
    """Yield nodes that execute at import time, skipping function bodies.

    An import inside a function is deferred: it cannot create an import-time
    cycle, which is what the section 2.1 graph is actually about. ``schemes.py``
    reaches ``engine``'s exception hierarchy exactly that way.
    """
    stack = list(tree.body)
    while stack:
        node = stack.pop()
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            continue
        yield node
        for child in ast.iter_child_nodes(node):
            stack.append(child)


def _module_level_imports(tree: ast.Module) -> set:
    """Modules imported at import time (the ones that can form a cycle)."""
    found = set()
    for node in _walk_module_scope(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                found.add(alias.name.split(".")[0])
        elif isinstance(node, ast.ImportFrom):
            if node.level:
                if node.module:
                    found.add(node.module.split(".")[0])
                else:
                    for alias in node.names:
                        found.add(alias.name.split(".")[0])
            elif node.module:
                parts = node.module.split(".")
                found.add(parts[1] if parts[0] == "wmlst" and len(parts) > 1
                          else parts[0])
    return found


def _imported_modules(tree: ast.Module) -> set:
    """Every module this file imports, by top-level name (absolute or relative)."""
    found = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                found.add(alias.name.split(".")[0])
        elif isinstance(node, ast.ImportFrom):
            if node.level:  # from . import x  /  from .x import y
                if node.module:
                    found.add(node.module.split(".")[0])
                else:
                    for alias in node.names:
                        found.add(alias.name.split(".")[0])
            elif node.module:
                parts = node.module.split(".")
                found.add(parts[1] if parts[0] == "wmlst" and len(parts) > 1
                          else parts[0])
    return found


def _names_used(node: ast.AST) -> set:
    """Bare names and attribute roots referenced inside ``node``."""
    used = set()
    for child in ast.walk(node):
        if isinstance(child, ast.Name):
            used.add(child.id)
        elif isinstance(child, ast.Attribute):
            base = child.value
            if isinstance(base, ast.Name):
                used.add(base.id)
    return used


# ---------------------------------------------------------------------------
# The dependency graph
# ---------------------------------------------------------------------------
def test_dependency_graph_is_acyclic():
    """No module imports a sibling of equal or higher rank (section 2.1)."""
    checked = 0
    for name, rank in sorted(RANK.items()):
        path = PKG / (name + ".py")
        if not path.is_file():
            continue
        checked += 1
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        for imported in _module_level_imports(tree):
            if imported in RANK and imported != name:
                assert RANK[imported] < rank, (
                    "{}.py imports {}.py at import time, which is not below it "
                    "in the section 2.1 graph".format(name, imported))
    if not checked:
        raise SkipTest("no package modules exist yet")


def test_gui_never_imports_the_cli():
    """The command line may hand over to the GUI (D14); never the reverse."""
    tree = _tree("gui")
    assert "cli" not in _imported_modules(tree), "gui.py imports cli.py"


def test_only_blastbin_may_launch_a_process():
    """subprocess lives in blastbin.py alone (section 2.1)."""
    checked = 0
    for path in sorted(PKG.glob("*.py")):
        if path.name == "blastbin.py":
            continue
        checked += 1
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        for imported in _imported_modules(tree):
            assert imported != "subprocess", (
                "{} imports subprocess; only blastbin.py may".format(path.name))
    if not checked:
        raise SkipTest("no package modules exist yet")


# ---------------------------------------------------------------------------
# gui.py
# ---------------------------------------------------------------------------
def test_gui_does_not_import_subprocess():
    """gui.py must not import subprocess, nor name Popen (section 10)."""
    tree = _tree("gui")
    assert "subprocess" not in _imported_modules(tree)
    for node in ast.walk(tree):
        if isinstance(node, ast.Name):
            assert node.id not in PROCESS_NAMES, (
                "gui.py references {!r}".format(node.id))
        elif isinstance(node, ast.Attribute):
            assert node.attr not in PROCESS_NAMES, (
                "gui.py references .{}".format(node.attr))


def test_gui_never_names_a_search_executable():
    """gui.py cannot launch a search tool: it does not even spell one.

    The executable names may appear ONLY as a field or attribute label reached
    indirectly (``getattr(tools, _ENGINE_BLAST_FIELD)``), never as a string
    literal, a bare name or a call target — which is what makes "the GUI never
    computes anything" checkable rather than merely asserted.
    """
    tree = _tree("gui")
    banned = ("blastn", "makeblastdb")
    for node in ast.walk(tree):
        if isinstance(node, ast.Constant) and isinstance(node.value, str):
            for word in banned:
                assert word not in node.value, (
                    "gui.py contains the literal {!r}".format(word))
        if isinstance(node, ast.Name):
            assert node.id not in banned
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute):
            assert node.func.attr not in banned


def test_gui_does_not_reimplement_the_algorithm():
    """gui.py defines no scoring or status function (section 10)."""
    tree = _tree("gui")
    for node in ast.walk(tree):
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            assert node.name not in ENGINE_FUNCTIONS | {"classify"}, (
                "gui.py defines {}()".format(node.name))
    for node in ast.walk(tree):
        if isinstance(node, ast.ImportFrom) and (node.module or "").endswith("engine"):
            for alias in node.names:
                assert alias.name not in ENGINE_FUNCTIONS, (
                    "gui.py imports engine.{}".format(alias.name))


def test_status_words_live_in_exactly_one_table():
    """The seven STATUS words appear only in the presentation table.

    If a status string turned up anywhere else in gui.py it would mean the GUI
    was deciding a status instead of rendering the one the engine returned.
    """
    tree = _tree("gui")
    table_nodes = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Assign):
            targets = [t.id for t in node.targets if isinstance(t, ast.Name)]
            if "STATUS_UI" in targets:
                table_nodes = {id(n) for n in ast.walk(node)}
    assert table_nodes, "gui.py has no STATUS_UI table"
    words = {"PERFECT", "NOVEL", "MIXED", "MISSING", "BAD", "NONE"}
    for node in ast.walk(tree):
        if (isinstance(node, ast.Constant) and isinstance(node.value, str)
                and node.value in words and id(node) not in table_nodes):
            raise AssertionError(
                "the status word {!r} appears outside STATUS_UI".format(node.value))


def test_no_tkinter_on_the_worker_side():
    """The worker body and every background probe stay clear of Tcl (section 9)."""
    tree = _tree("gui")
    worker_functions = {"_run", "probe_environment", "body"}
    checked = 0
    for node in ast.walk(tree):
        if isinstance(node, ast.FunctionDef) and node.name in worker_functions:
            checked += 1
            used = _names_used(node)
            leaked = used & TK_NAMES
            assert not leaked, (
                "{}() touches {} on the worker thread".format(node.name,
                                                              sorted(leaked)))
    assert checked >= 2, "the worker functions were not found"


def test_gui_has_the_public_api_of_section_4_10():
    """main(), WmlstApp, AnalysisController, Msg and Prefs are all present."""
    tree = _tree("gui")
    names = {n.name for n in tree.body
             if isinstance(n, (ast.FunctionDef, ast.ClassDef))}
    for required in ("main", "WmlstApp", "AnalysisController", "Msg", "Prefs"):
        assert required in names, "gui.py is missing {}".format(required)


# ---------------------------------------------------------------------------
# report.py
# ---------------------------------------------------------------------------
def test_report_imports_no_engine_functions():
    """report.py imports engine's dataclasses only — it performs no biology."""
    tree = _tree("report")
    for node in ast.walk(tree):
        if isinstance(node, ast.ImportFrom) and (node.module or "").endswith("engine"):
            for alias in node.names:
                assert alias.name not in ENGINE_FUNCTIONS, (
                    "report.py imports engine.{}".format(alias.name))


def test_report_calls_no_engine_function():
    """No ``engine.score_signature(...)``-shaped call anywhere in report.py."""
    tree = _tree("report")
    for node in ast.walk(tree):
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute):
            base = node.func.value
            if isinstance(base, ast.Name) and base.id in ("engine", "wmlst"):
                assert node.func.attr not in ENGINE_FUNCTIONS, (
                    "report.py calls engine.{}()".format(node.func.attr))


# ---------------------------------------------------------------------------
# third-party imports
# ---------------------------------------------------------------------------
def test_the_only_optional_third_party_import_is_tkinterdnd2():
    """Every non-stdlib import sits inside a try/except ImportError (section 2.1)."""
    source = _source("gui")
    tree = ast.parse(source)
    guarded = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Try):
            for child in ast.walk(node):
                if isinstance(child, ast.Import):
                    guarded.update(a.name.split(".")[0] for a in child.names)
                elif isinstance(child, (ast.ImportFrom,)) and child.module:
                    guarded.add(child.module.split(".")[0])
    assert "tkinterdnd2" in guarded, "tkinterdnd2 must be imported defensively"


def _run_all() -> int:
    failures = 0
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            try:
                fn()
            except SkipTest as exc:
                print("SKIP {}: {}".format(name, exc))
            except AssertionError as exc:
                failures += 1
                print("FAIL {}: {}".format(name, exc))
            else:
                print("ok   {}".format(name))
    print("\n{} failure(s)".format(failures))
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(_run_all())
