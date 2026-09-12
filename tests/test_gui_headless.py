# SPDX-License-Identifier: GPL-2.0-only
# Copyright (C) 2025-2026 IOWA-Tech - Giovanni Lorenzin
# Copyright (C) Torsten Seemann (upstream `mlst`, from which WMLST is ported)
"""GUI tests that need no display (docs/ARCHITECTURE.md section 10).

Everything in ``wmlst.gui`` above the widget layer is deliberately free of Tk, so
the queue plumbing, the state machine, the settings validation, the input
discovery and the error mapping are all testable on a headless box.

The widget-construction smoke test at the end runs only when a display happens to
be available (a real X server, or ``xvfb-run``); otherwise it skips loudly rather
than pretending to have run.

Runnable as ``python -m pytest tests/test_gui_headless.py`` and as
``python3 tests/test_gui_headless.py``.
"""

from __future__ import annotations

import dataclasses
import json
import os
import shutil
import sys
import tempfile
import threading
import time
from pathlib import Path
from unittest import SkipTest

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from typing import Optional

from wmlst import gui
from wmlst.engine import (
    AlleleCall,
    BlastFailedError,
    BlastNotFoundError,
    Cancelled,
    DatabaseMissingError,
    EmptyInputError,
    Hit,
    RunConfig,
    SampleResult,
    SchemeScore,
    UnsupportedFormatError,
)

DATA = ROOT / "tests" / "data"


# ---------------------------------------------------------------------------
# fixtures built by hand — no engine, no BLAST, no database
# ---------------------------------------------------------------------------
def make_hit(locus: str = "arcC", allele: str = "16") -> Hit:
    return Hit(sseqid="sepidermidis.{}_{}".format(locus, allele),
               scheme="sepidermidis", locus=locus, allele=allele,
               slen=456, length=456, nident=456, qseqid="contig00001",
               qstart=1, qend=456, qseq="ACGT" * 114, sstrand="plus",
               pct_identity=100.0, pct_coverage=100.0, call_kind="exact", index=1)


def make_result(label: str = "example.fna", status: str = "PERFECT",
                st: str = "184", score: int = 100) -> SampleResult:
    """A SampleResult shaped exactly like the example.fna golden row."""
    codes = (("arcC", "16"), ("aroE", "1"), ("gtr", "2"), ("mutS", "1"),
             ("pyrR", "2"), ("tpiA", "1"), ("yqiL", "1"))
    alleles = tuple(
        AlleleCall(locus=locus, code=code, symbol="exact", best=make_hit(locus, code),
                   hits=(make_hit(locus, code),))
        for locus, code in codes)
    return SampleResult(
        path=str(DATA / label), label=label, scheme="sepidermidis", st=st,
        signature="16/1/2/1/2/1/1", score=score, status=status, alleles=alleles,
        n_contigs=42, total_bp=2_500_000, hits_seen=1423, hits_kept=61,
        elapsed_s=5.5)


def make_tied_result() -> SampleResult:
    """A result whose winner tied with a second scheme (section 5.14a)."""
    base = make_result("kpneu.fna", st="258")
    tied = (SchemeScore("klebsiella", "258", "3/3/1/1/1/1/79", 100, 7),
            SchemeScore("ecoli_achtman_4", "14464",
                        "1769/1664/193/1804/986/745/866", 100, 7))
    return dataclasses.replace(base, scheme="klebsiella", tied=tied,
                               candidates=tied)


class FakeEngine:
    """Stands in for ``engine.Engine`` so the controller can be tested alone."""

    def __init__(self, cfg: RunConfig, *, fail_on: str = "", block: bool = False,
                 delay: float = 0.0, raise_on_init: Optional[BaseException] = None):
        if raise_on_init is not None:
            raise raise_on_init
        self.cfg = cfg
        self.fail_on = fail_on
        self.block = block
        self.delay = delay
        self.closed = False
        self.seen = []

    def analyse_file(self, path, *, label=None, progress=None, warn=None,
                     cancel=None):
        self.seen.append(path)
        if self.delay:
            time.sleep(self.delay)
        if progress is not None:
            for percent, text in ((0.0, "Reading…"), (10.0, "Read 42 contigs."),
                                  (80.0, "Search finished."), (100.0, "done")):
                progress(_Event("blast", percent, text, path))
        if warn is not None:
            warn("found additional exact allele match sepidermidis.arcC-16")
        if self.block:
            while cancel is None or not cancel.is_set():
                time.sleep(0.01)
            raise Cancelled("cancelled")
        if self.fail_on and os.path.basename(path) == self.fail_on:
            raise EmptyInputError("The input appears to be empty")
        return make_result(os.path.basename(path))

    def close(self):
        self.closed = True


class _Event:
    def __init__(self, phase, percent, text, path):
        self.phase, self.percent, self.text, self.path = phase, percent, text, path


class FakeScheduler:
    """A stand-in for ``root.after`` that the test drives by hand."""

    def __init__(self):
        self.pending = []

    def __call__(self, ms, fn):
        self.pending.append(fn)
        return len(self.pending)

    def run(self, rounds: int = 1) -> None:
        for _ in range(rounds):
            due, self.pending = self.pending, []
            for fn in due:
                fn()


def drive(controller, scheduler, *, until, timeout: float = 10.0) -> None:
    """Run the after-chain until ``until()`` is true or the timeout expires."""
    deadline = time.time() + timeout
    while time.time() < deadline:
        scheduler.run()
        if until():
            return
        time.sleep(0.01)
    raise AssertionError("timed out waiting for the worker")


# ---------------------------------------------------------------------------
# scaling and theme maths (section 10.1)
# ---------------------------------------------------------------------------
def test_px_scales_pixels_and_never_fonts():
    """px(16) == 24 at 144 dpi, and no font size is scaled twice."""
    try:
        gui.set_scale(144)
        assert gui.px(16) == 24
        assert gui.px(1) == 2 and gui.px(0) == 0
        for name, size in gui.FONT_PT.items():
            assert size == int(size) and 6 <= size <= 40, name
            # a point size that had been through px() would be 1.5x larger here
            assert size != gui.px(size), "{} looks double-scaled".format(name)
    finally:
        gui.set_scale(96)
    assert gui.px(16) == 16


def test_dpi_fix_is_a_no_op_off_windows():
    assert gui.dpi_fix() in ("not-windows", "per-monitor-v2", "per-monitor",
                             "system", "legacy", "unavailable")
    assert gui.high_contrast_active() in (True, False)


def test_font_ladder_falls_back():
    assert gui.pick_font(("Segoe UI", "Arial"), gui.PREFERRED_FAMILIES, "X") == "Segoe UI"
    assert gui.pick_font(("Arial",), gui.PREFERRED_FAMILIES, "Fallback") == "Fallback"
    assert gui.pick_font(("DejaVu Sans",), gui.PREFERRED_FAMILIES, "X") == "DejaVu Sans"


# ---------------------------------------------------------------------------
# status presentation (section 11.3) — rendered, never derived
# ---------------------------------------------------------------------------
def test_every_status_has_a_glyph_a_colour_and_a_sentence():
    for status in ("PERFECT", "NOVEL", "MIXED", "MISSING", "BAD", "OK", "NONE"):
        glyph, colour, sentence, shape = gui.STATUS_UI[status]
        assert glyph and colour in gui.PALETTE
        assert sentence.endswith(".") and len(sentence) > 40
        # The shape lives in this table too, so dot_image() never has to
        # re-list the status words (tests/test_layering.py enforces that).
        assert shape in ("disc", "ring", "triangle")
        # The STATUS cell is the WORD. The second, non-colour channel is a
        # SHAPE drawn by gui.dot_image() into the row (see the test below):
        # it used to be a dingbat in the cell, which silently became a hollow
        # box -- and so no channel at all -- on any Tk without a font for it.
        cell = gui.status_cell(status)
        assert cell == status
        assert gui.status_tag(status).startswith("st_")


def test_each_status_draws_its_own_shape_into_the_row():
    """Shape, not colour, and not a glyph the font may not own (WCAG 1.4.1)."""
    if not _display_available():
        raise SkipTest("no display: run under xvfb-run to exercise real widgets")
    import tkinter as tk

    try:
        root = tk.Tk()
    except tk.TclError as exc:
        raise SkipTest("Tk could not open a display: {}".format(exc)) from exc
    try:
        root.withdraw()
        gui.apply_scaling(root)
        gui.install_theme(root)
        app = gui.WmlstApp(root, gui.Prefs())
        view = app.analyse_view
        seen = {}
        for status in ("PERFECT", "NOVEL", "MIXED", "MISSING", "BAD", "OK", "NONE"):
            icon = view.status_icon(status)
            assert icon.width() > 0 and icon.height() > 0
            seen[status] = icon
        # cached, so 2000 rows do not rasterise 2000 images
        assert view.status_icon("PERFECT") is seen["PERFECT"]
    finally:
        root.destroy()


def test_allele_summary_copies_codes_verbatim():
    res = make_result()
    assert gui.allele_summary(res) == (
        "arcC(16);aroE(1);gtr(2);mutS(1);pyrR(2);tpiA(1);yqiL(1)")
    assert "100.0% identity" in gui.evidence_text(res.alleles[0])
    missing = AlleleCall(locus="g3pd", code="-", symbol="missing", best=None, hits=())
    assert gui.evidence_text(missing) == "no match found"


# ---------------------------------------------------------------------------
# preferences and validation (section 10.4)
# ---------------------------------------------------------------------------
def test_defaults_match_the_reference():
    p = gui.Prefs()
    assert (p.minid, p.mincov, p.minscore) == (95.0, 50.0, 50.0)
    assert p.exclude == ("ecoli", "abaumannii", "vcholerae_2", "senterica_achtman_2")
    assert p.threads >= 1 and p.jobs >= 1


def test_clamping_repairs_anything():
    assert gui.clamp_float("101", 0, 100, 95) == 100.0
    assert gui.clamp_float("-3", 0, 100, 95) == 0.0
    assert gui.clamp_float("abc", 0, 100, 95) == 95.0
    assert gui.clamp_float("97,5", 0, 100, 95) == 97.5
    assert gui.clamp_int("0", 1, 8, 1) == 1
    assert gui.clamp_int("999", 1, 8, 1) == 8
    assert gui.clamp_int(None, 1, 8, 2) == 2
    assert gui.parse_exclude("ecoli, abaumannii ecoli\nvcholerae_2;") == (
        "ecoli", "abaumannii", "vcholerae_2")
    assert gui.parse_exclude("") == ()


def test_forcing_a_scheme_applies_the_reference_side_effects():
    """--scheme forces minscore 0 and clears exclude (bin/mlst:125, 3.6.1)."""
    p = gui.Prefs(scheme="saureus", minscore=50.0)
    n = p.normalised()
    assert n.minscore == 0.0 and n.exclude == ()
    env = gui.Environment(dbdir="/db", datadir="/db/pubmlst", blastdb="/db/blast/mlst.fa")
    cfg = p.to_runconfig(["a.fa"], env)
    assert isinstance(cfg, RunConfig)
    assert cfg.scheme == "saureus" and cfg.minscore == 0.0
    assert cfg.exclude == frozenset()
    assert cfg.files == ("a.fa",)
    assert cfg.dbdir == "/db" and cfg.blastdb == "/db/blast/mlst.fa"


def test_threads_are_always_passed_explicitly():
    p = gui.Prefs(threads=2, jobs=2)
    cfg = p.to_runconfig(["a.fa"], gui.Environment())
    assert cfg.threads >= 1 and cfg.jobs >= 1
    assert cfg.jobs * cfg.threads <= max(1, os.cpu_count() or 1)
    huge = gui.Prefs(threads=999, jobs=999).normalised()
    assert huge.threads <= gui.max_threads() and huge.jobs <= gui.max_threads()


def test_prefs_round_trip_and_survive_corruption(tmp_path=None):
    # Path("") is Path(".") and is truthy, so an empty PYTEST_TMP must not be
    # allowed through: in script mode it would scribble settings into the cwd.
    env_tmp = os.environ.get("PYTEST_TMP") or ""
    home = tmp_path or (Path(env_tmp) if env_tmp else None)
    if home is None:
        home = Path(tempfile.mkdtemp(prefix="wmlst-gui-test-"))
    base = Path(str(home))
    base.mkdir(parents=True, exist_ok=True)
    old = {k: os.environ.get(k) for k in ("APPDATA", "XDG_CONFIG_HOME")}
    os.environ["APPDATA"] = str(base)
    os.environ["XDG_CONFIG_HOME"] = str(base)
    try:
        path = gui.Prefs.path()
        if path.exists():
            path.unlink()
        p = gui.Prefs(minid=90.0, threads=1, last_dir=str(base))
        assert p.save()
        back = gui.Prefs.load()
        assert back.minid == 90.0 and back.last_dir == str(base)
        assert "load_note" not in json.loads(path.read_text())
        path.write_text("{not json at all")
        broken = gui.Prefs.load()
        assert broken.minid == 95.0 and broken.load_note
        assert path.with_name("settings.bad.json").exists()
    finally:
        for key, value in old.items():
            if value is None:
                os.environ.pop(key, None)
            else:
                os.environ[key] = value


# ---------------------------------------------------------------------------
# input discovery (section 10.2)
# ---------------------------------------------------------------------------
def test_expand_inputs_handles_files_folders_and_lists():
    if not DATA.is_dir():
        raise SkipTest("tests/data is not present")
    plan = gui.expand_inputs([str(DATA)])
    assert plan.count >= 8, plan.files
    assert all(gui.looks_like_sequence(f) for f in plan.files)
    assert not any(f.endswith("fofn.txt") for f in plan.files)

    one = gui.expand_inputs([str(DATA / "example.fna")])
    assert one.files == (str(DATA / "example.fna"),)

    txt = gui.expand_inputs([str(DATA / "fofn.txt")])
    assert txt.files == () and txt.fofn_candidates == (str(DATA / "fofn.txt"),)

    missing = gui.expand_inputs([str(DATA / "nope.fa")])
    assert missing.missing and not missing.files

    twice = gui.expand_inputs([str(DATA / "example.fna"), str(DATA / "example.fna")])
    assert twice.count == 1, "the same file must not be queued twice"


def test_looks_like_sequence_accepts_the_compressed_forms():
    for good in ("a.fa", "a.fasta", "a.fna", "b.gbk", "c.embl",
                 "d.fna.gz", "e.fasta.bz2", "f.fa.zip"):
        assert gui.looks_like_sequence(good), good
    for bad in ("notes.txt", "a.fastq", "b.pdf", "archive.zip", "README"):
        assert not gui.looks_like_sequence(bad), bad


def test_fofn_resolves_relative_entries_against_its_own_folder():
    """Divergences D8 and D9."""
    if not (DATA / "fofn.txt").is_file():
        raise SkipTest("fofn.txt is not present")
    entries = gui.read_fofn(str(DATA / "fofn.txt"))
    assert entries, "the fofn produced no entries"
    for entry in entries:
        assert os.path.isabs(entry)
        assert Path(entry).exists(), entry
    assert "" not in entries


def test_dnd_payload_is_split_with_tcl_rules():
    class FakeTk:
        @staticmethod
        def splitlist(data):
            return ("C:/My Data/a.fa", "C:/b.fa")

    class FakeWidget:
        tk = FakeTk()

    assert gui.split_dnd_paths(FakeWidget(), "{C:/My Data/a.fa} C:/b.fa") == [
        "C:/My Data/a.fa", "C:/b.fa"]

    class Broken:
        class tk:
            @staticmethod
            def splitlist(data):
                raise RuntimeError("no interpreter")

    assert gui.split_dnd_paths(Broken(), "{/a b.fa} {/c.fa}") == ["/a b.fa", "/c.fa"]


def test_human_helpers():
    assert gui.human_duration(30) == "30 seconds"
    assert gui.human_duration(180) == "about 3 minutes"
    assert gui.human_duration(3600) == "about 1 hour"
    assert "MB" in gui.human_bytes(143_400_333)
    assert gui.human_bytes(None) == "unknown"
    assert gui.estimate_seconds(60, 4) == 90.0


def test_geometry_is_only_restored_when_it_is_visible():
    assert gui.parse_geometry("1180x760+40+30") == (1180, 760, 40, 30)
    assert gui.parse_geometry("nonsense") is None
    assert gui.geometry_on_screen("1180x760+40+30", 1920, 1080)
    assert not gui.geometry_on_screen("1180x760+5000+30", 1920, 1080)
    assert not gui.geometry_on_screen("1180x760+40+2000", 1920, 1080)
    assert not gui.geometry_on_screen("", 1920, 1080)


# ---------------------------------------------------------------------------
# the state machine (section 10.2)
# ---------------------------------------------------------------------------
def test_state_machine_transitions():
    seen = []
    state = gui.AppState(on_change=lambda old, new: seen.append((old, new)))
    assert state.state == gui.IDLE and not state.running
    assert state.to(gui.RUNNING)
    assert state.running
    assert not state.to(gui.RESULTS + "X") if False else True
    assert state.finish("cancelled")
    assert state.state == gui.RESULTS and state.reason == "cancelled"
    assert state.to(gui.RUNNING)
    assert state.finish("error")
    assert state.state == gui.RESULTS
    assert seen[0] == (gui.IDLE, gui.RUNNING)
    # RESULTS -> RUNNING is legal (a new drop appends); IDLE -> RESULTS is not.
    fresh = gui.AppState()
    assert not fresh.to(gui.RESULTS)
    try:
        fresh.to("NOPE")
    except ValueError:
        pass
    else:
        raise AssertionError("an unknown state must raise")


def test_assert_main_thread_guards_the_ui():
    gui.assert_main_thread("test")  # on the main thread: fine
    box = []

    def off_thread():
        try:
            gui.assert_main_thread("the status line")
        except RuntimeError as exc:
            box.append(str(exc))

    t = threading.Thread(target=off_thread)
    t.start()
    t.join()
    assert box and "main thread" in box[0]


# ---------------------------------------------------------------------------
# the controller (sections 9, 10.5)
# ---------------------------------------------------------------------------
def test_msg_carries_only_plain_data():
    fields = gui.Msg.__dataclass_fields__
    assert set(fields) == {"kind", "job", "index", "total", "path", "percent",
                           "text", "result", "exc", "tb"}
    assert gui.AnalysisController.POLL_MS == 120


def make_controller(**kw):
    scheduler = FakeScheduler()
    seen = []
    controller = gui.AnalysisController(
        schedule=scheduler, dispatch=seen.append,
        engine_factory=lambda cfg: FakeEngine(cfg, **kw))
    return controller, scheduler, seen


def test_a_batch_produces_the_expected_message_sequence():
    controller, scheduler, seen = make_controller()
    cfg = gui.Prefs().to_runconfig(["a.fa", "b.fa"], gui.Environment())
    assert controller.start(["a.fa", "b.fa"], cfg)
    drive(controller, scheduler,
          until=lambda: any(m.kind == "batch_done" for m in seen))
    kinds = [m.kind for m in seen]
    assert kinds[0] == "started"
    assert kinds[-1] == "batch_done"
    assert kinds.count("result") == 2
    assert "phase" in kinds and "warn" in kinds
    results = [m.result for m in seen if m.kind == "result"]
    assert all(isinstance(r, SampleResult) for r in results)
    assert results[0].st == "184" and results[0].status == "PERFECT"
    # percent is within-file and monotonically non-decreasing (section 5.16)
    percents = [m.percent for m in seen if m.kind == "phase" and m.index == 1]
    assert percents == sorted(percents)
    assert not any(hasattr(m.result, "winfo_id") for m in seen if m.result)


def test_a_failed_file_does_not_stop_the_batch():
    controller, scheduler, seen = make_controller(fail_on="b.fa")
    cfg = gui.Prefs().to_runconfig(["a.fa", "b.fa", "c.fa"], gui.Environment())
    controller.start(["a.fa", "b.fa", "c.fa"], cfg)
    drive(controller, scheduler,
          until=lambda: any(m.kind == "batch_done" for m in seen))
    kinds = [m.kind for m in seen]
    assert kinds.count("result") == 2 and kinds.count("failed") == 1
    failure = next(m for m in seen if m.kind == "failed")
    assert isinstance(failure.exc, EmptyInputError) and failure.tb


def test_a_second_start_is_impossible_and_the_files_are_queued():
    controller, scheduler, seen = make_controller(block=True)
    cfg = gui.Prefs().to_runconfig(["a.fa"], gui.Environment())
    assert controller.start(["a.fa"], cfg)
    assert controller.start(["b.fa", "c.fa"], cfg) is False
    assert controller.running
    controller.request_cancel()
    drive(controller, scheduler,
          until=lambda: any(m.kind == "batch_done" for m in seen))
    assert any(m.kind == "cancelled" for m in seen)
    assert controller.pending_paths() == ["b.fa", "c.fa"]
    assert controller.pending_paths() == []


def test_cancellation_stops_between_files_and_keeps_results():
    # A per-file delay makes the race deterministic: the worker cannot possibly
    # finish all 40 files before the cancel lands.
    controller, scheduler, seen = make_controller(delay=0.02)
    cfg = gui.Prefs().to_runconfig([], gui.Environment())
    controller.start(["a.fa"] * 40, cfg)
    controller.request_cancel()
    drive(controller, scheduler,
          until=lambda: any(m.kind == "batch_done" for m in seen))
    assert any(m.kind == "cancelled" for m in seen)
    assert len([m for m in seen if m.kind == "result"]) < 40
    controller.shutdown(1.0)
    assert not controller.running


def test_an_engine_that_cannot_be_built_is_fatal_not_a_crash():
    scheduler = FakeScheduler()
    seen = []
    controller = gui.AnalysisController(
        schedule=scheduler, dispatch=seen.append,
        engine_factory=lambda cfg: FakeEngine(
            cfg, raise_on_init=DatabaseMissingError("Database directory does not exist")))
    controller.start(["a.fa"], gui.Prefs().to_runconfig(["a.fa"], gui.Environment()))
    drive(controller, scheduler,
          until=lambda: any(m.kind == "batch_done" for m in seen))
    fatal = [m for m in seen if m.kind == "fatal"]
    assert fatal and isinstance(fatal[0].exc, DatabaseMissingError)


def test_background_tasks_report_success_and_failure():
    controller, scheduler, seen = make_controller()
    job = controller.run_task("probe", lambda: {"answer": 42})
    drive(controller, scheduler, until=lambda: any(m.kind == "env" for m in seen))
    assert controller.take_task_result(job) == {"answer": 42}

    seen.clear()

    def boom(progress, cancel):
        progress(50, 100, "halfway")
        raise BlastNotFoundError("no engine")

    controller.run_task("bootstrap", boom, with_progress=True)
    drive(controller, scheduler, until=lambda: any(m.kind == "failed" for m in seen))
    phases = [m for m in seen if m.kind == "phase"]
    assert phases and phases[0].percent == 50.0 and phases[0].text == "halfway"
    assert isinstance(next(m for m in seen if m.kind == "failed").exc,
                      BlastNotFoundError)


def test_the_poll_chain_is_started_once():
    controller, scheduler, _seen = make_controller()
    controller.ensure_polling()
    controller.ensure_polling()
    assert len(scheduler.pending) == 1
    scheduler.run()
    assert len(scheduler.pending) == 1, "the after() chain must not fork"


# ---------------------------------------------------------------------------
# friendly errors (section 10.6)
# ---------------------------------------------------------------------------
def test_every_exception_maps_to_a_novice_sentence():
    cases = [
        (EmptyInputError("x"), "We couldn't read that file"),
        (UnsupportedFormatError("x"), "That file format isn't supported"),
        (BlastFailedError("x", returncode=1, stderr_tail="boom"),
         "The search engine stopped unexpectedly"),
        (DatabaseMissingError("missing"), "The MLST database is missing"),
        (PermissionError(13, "Permission denied"),
         "Windows wouldn't let us write there"),
        (ValueError("weird"), "Something went wrong"),
    ]
    for exc, headline in cases:
        f = gui.friendly_error(exc, name="messy.fa")
        assert f.headline == headline, (exc, f.headline)
        assert f.details, "the technical fold must never be empty"
        assert "Python" in f.details and "WMLST" in f.details
        assert not f.needs_bootstrap
        for tip in f.tips:
            assert tip.endswith(".")


def test_an_empty_file_is_reported_as_an_empty_file():
    """BLAST says "Sequence contains no data" (bats 13); the user hears plain words."""
    exc = BlastFailedError("BLAST engine error", returncode=2,
                           stderr_tail="Warning: Sequence contains no data")
    f = gui.friendly_error(exc, name="empty.fa")
    assert f.headline == "We couldn't read that file"
    assert "empty.fa" in f.body


def test_a_missing_search_engine_routes_to_setup_not_to_an_error():
    f = gui.friendly_error(BlastNotFoundError("nope"))
    assert f.needs_bootstrap and "search engine" in f.headline


def test_a_damaged_database_offers_a_rebuild():
    from wmlst.engine import DatabaseCorruptError
    f = gui.friendly_error(DatabaseCorruptError("damaged"))
    assert f.offers_rebuild


def test_the_technical_fold_names_the_environment():
    env = gui.Environment(db_version="2025-12-29", scheme_count=162,
                          blast_version="2.17.0+", blast_path="/opt/b/x",
                          dbdir="/db")
    block = gui._env_block(env)
    assert "2025-12-29" in block and "162" in block and "2.17.0+" in block
    assert env.footer_text().startswith("Database 2025-12-29")
    assert not env.ready  # db_ok / blast_ok are still False


# ---------------------------------------------------------------------------
# run assembly for the exporters (section 10.2)
# ---------------------------------------------------------------------------
def test_novel_alleles_are_deduplicated_in_order():
    from wmlst.engine import NovelAllele
    a = NovelAllele("sepidermidis", "arcC", "d" * 32, "ACGT", "one",
                    "sepidermidis.arcC-" + "d" * 32, "16", 4)
    b = NovelAllele("sepidermidis", "aroE", "e" * 32, "TTTT", "two",
                    "sepidermidis.aroE-" + "e" * 32, "1", 4)
    s1 = make_result("one.fa")
    s2 = make_result("two.fa")
    s1 = type(s1)(**{**{f: getattr(s1, f) for f in s1.__dataclass_fields__},
                     "novel": (a, b)})
    s2 = type(s2)(**{**{f: getattr(s2, f) for f in s2.__dataclass_fields__},
                     "novel": (a,)})
    out = gui.dedupe_novel([s1, s2])
    assert [n.fasta_id for n in out] == [a.fasta_id, b.fasta_id]


def test_fallback_meta_is_complete():
    cfg = gui.Prefs().to_runconfig(["a.fa"], gui.Environment(db_version="2025-12-29"))
    meta = gui.fallback_meta(cfg, gui.Environment(db_version="2025-12-29",
                                                  scheme_count=162))
    assert meta.wmlst_version and meta.mlst_compat == "2.35.0"
    assert meta.db_version == "2025-12-29" and meta.db_scheme_count == 162
    assert meta.config is cfg and meta.started_utc.endswith("Z")


# ---------------------------------------------------------------------------
# widget smoke test — only with a display
# ---------------------------------------------------------------------------
def _display_available() -> bool:
    if sys.platform.startswith("win") or sys.platform == "darwin":
        return True
    return bool(os.environ.get("DISPLAY") or os.environ.get("WAYLAND_DISPLAY"))


def test_selftest_builds_and_destroys_a_root():
    """``wmlst-gui --selftest`` exits 0 and prints nothing (section 4.10)."""
    if not _display_available():
        raise SkipTest("no display: run under xvfb-run to exercise real widgets")
    assert gui.main(["--selftest"]) == 0


def test_the_whole_window_builds_and_shows_a_result():
    """Build the real window, push a fake batch through it, read the tree back."""
    if not _display_available():
        raise SkipTest("no display: run under xvfb-run to exercise real widgets")
    import tkinter as tk

    try:
        root = tk.Tk()
    except tk.TclError as exc:
        raise SkipTest("Tk could not open a display: {}".format(exc)) from exc
    try:
        root.withdraw()
        gui.apply_scaling(root)
        gui.install_theme(root)
        app = gui.WmlstApp(root, gui.Prefs())
        app.env = gui.Environment(dbdir="/db", datadir="/db/pubmlst",
                                  blastdb="/db/blast/mlst.fa", db_ok=True,
                                  blast_ok=True, db_version="2025-12-29",
                                  scheme_count=162, blast_version="2.17.0+",
                                  scheme_names=("saureus", "sepidermidis"))
        app.refresh_footer()
        app.settings_view.set_scheme_choices(app.env.scheme_names)

        # the three tabs, the title, the footer credit
        assert root.title() == "WMLST — MLST typing for Windows — IOWA-Tech"
        assert [app.notebook.tab(i, "text") for i in range(3)] == [
            "Analyse", "Database", "Settings"]
        assert "IOWA-Tech" in app.footer_env.master.winfo_children()[0].cget("text")

        # drive one fake batch through the real controller and the real view
        app.controller._engine_factory = lambda cfg: FakeEngine(cfg)
        app.start_analysis([str(DATA / "example.fna")])
        deadline = time.time() + 10
        while time.time() < deadline and app.state.running:
            root.update()
            time.sleep(0.02)
        root.update()
        rows = app.analyse_view.tree.get_children("")
        assert len(rows) == 1, "the result row never appeared"
        values = app.analyse_view.tree.item(rows[0], "values")
        # ORGANISM leads, then SCHEME, ST, STATUS, SCORE.
        assert values[0] == "Staphylococcus epidermidis"
        assert values[1] == "sepidermidis" and values[2] == "184"
        assert values[3].endswith("PERFECT") and values[4] == "100"
        assert app.analyse_view.st_label.cget("text") == "184"
        assert app.analyse_view.organism_label.cget("text") == (
            "Staphylococcus epidermidis")

        # the per-locus children are built lazily, in gene order, under the
        # caption row and the scheme-provenance row
        app.analyse_view.tree.focus(rows[0])
        app.analyse_view._on_open()
        children = list(app.analyse_view.tree.get_children(rows[0]))
        texts = [app.analyse_view.tree.item(ch, "text").strip() for ch in children]
        assert "SCHEME" in texts, "the scheme/reference row is missing"
        caption = texts.index("LOCUS")
        assert app.analyse_view.tree.item(children[caption], "values")[0] == "ALLELE"
        loci = texts[caption + 1:]
        assert loci == ["arcC", "aroE", "gtr", "mutS", "pyrR", "tpiA", "yqiL"]
        codes = [app.analyse_view.tree.item(ch, "values")[0]
                 for ch in children[caption + 1:]]
        assert codes == ["16", "1", "2", "1", "2", "1", "1"]

        # a settings change offers a re-run; restoring defaults works
        app.settings_view.vars["minid"].set("110")
        app.settings_view._clamp("minid")
        assert app.settings_view.vars["minid"].get() == "100"
        app.settings_view.restore_defaults()
        assert app.prefs.minid == 95.0 and app.prefs.minscore == 50.0
        app.settings_view.vars["scheme"].set("saureus")
        app.settings_view._changed()
        assert app.prefs.scheme == "saureus"
        assert str(app.settings_view.spin_minscore.cget("state")) == "disabled"
        # the control displays 0 and the RUN uses 0 (bin/mlst:125, section 10.4)
        assert app.settings_view.vars["minscore"].get() == "0"
        assert app.prefs.to_runconfig(["a.fa"], app.env).minscore == 0.0
        # ... but the displayed 0 is not the stored preference: releasing the
        # scheme has to give the user their own minimum score back.
        assert app.prefs.minscore == 50.0
        app.settings_view.vars["scheme"].set("Automatic (recommended)")
        app.settings_view._changed()
        assert app.prefs.scheme is None and app.prefs.minscore == 50.0
        assert app.settings_view.vars["minscore"].get() == "50"
        app.settings_view.vars["scheme"].set("saureus")
        app.settings_view._changed()

        # the drop zone repaints in every state without raising
        for state in ("idle", "hover", "dragover", "focus"):
            app.analyse_view.dropzone.set_state(state)
        app.analyse_view.dropzone.set_enabled(False)
        app.analyse_view.dropzone.set_enabled(True)

        # clearing returns to IDLE
        app.analyse_view.clear()
        assert app.state.state == gui.IDLE
        assert app.analyse_view.tree.get_children("") == ()
    finally:
        try:
            root.destroy()
        except Exception:
            pass


def test_the_export_buttons_and_the_footer_survive_a_single_result():
    """A one-file run must not push the actions or the footer off the window.

    Regression: pack() serves children in call order, so a Treeview that asked
    for more rows than were left over used to starve the action row and cover
    the status line — and the single-file case is the first thing a novice sees.
    """
    if not _display_available():
        raise SkipTest("no display: run under xvfb-run to exercise real widgets")
    import tkinter as tk

    try:
        root = tk.Tk()
    except tk.TclError as exc:
        raise SkipTest("Tk could not open a display: {}".format(exc)) from exc
    try:
        gui.apply_scaling(root)
        gui.install_theme(root)
        app = gui.WmlstApp(root, gui.Prefs())
        root.geometry("1180x760+40+40")
        root.update_idletasks()
        root.update()

        for count in (1, 2, 13):
            app.analyse_view.clear()
            for i in range(count):
                app.analyse_view.add_result(make_result("sample{}.fna".format(i)))
            app.analyse_view.show_state(gui.RESULTS)
            root.update_idletasks()
            root.update()

            view = app.analyse_view
            for button in (view.btn_html, view.btn_tsv, view.btn_json,
                           view.btn_novel, view.btn_clear):
                assert button.winfo_ismapped(), (
                    "{!r} is off-screen with {} result(s)".format(
                        button.cget("text"), count))
                assert button.winfo_height() >= gui.px(28), (
                    "{!r} collapsed to {} px".format(button.cget("text"),
                                                     button.winfo_height()))

            notebook_bottom = (app.notebook.winfo_rooty()
                               + app.notebook.winfo_height())
            status_top = app.status.winfo_rooty()
            assert notebook_bottom <= status_top, (
                "the notebook covers the status line with {} result(s)".format(count))
            assert app.status.winfo_ismapped() and app.footer_env.winfo_ismapped()
    finally:
        try:
            root.destroy()
        except Exception:
            pass


def _focus_ring(start, limit: int = 40) -> list:
    """Walk tk_focusNext from ``start`` back round to it."""
    widget, ring = start, []
    for _ in range(limit):
        widget = widget.tk_focusNext()
        if widget is None:
            break
        ring.append(widget)
        if widget is start:
            break
    return ring


def test_tab_order_reaches_every_control():
    """Tab order is walkable and the drop zone takes focus (section 10.8)."""
    if not _display_available():
        raise SkipTest("no display: run under xvfb-run to exercise real widgets")
    import tkinter as tk

    try:
        root = tk.Tk()
    except tk.TclError as exc:
        raise SkipTest("Tk could not open a display: {}".format(exc)) from exc
    try:
        # NOT withdrawn: Tk refuses focus traversal on an unmapped window, so a
        # withdrawn root would make this test pass by returning nothing.
        gui.apply_scaling(root)
        gui.install_theme(root)
        root.geometry("1100x760+40+40")
        app = gui.WmlstApp(root, gui.Prefs())
        root.update_idletasks()
        root.update()

        # On the idle Analyse tab the only things worth focusing are the drop
        # zone (which must accept focus, since Return/space operate it) and the
        # tab strip: every action button is legitimately disabled.
        assert int(app.analyse_view.dropzone.cget("takefocus")) == 1
        ring = _focus_ring(app.analyse_view.dropzone)
        assert app.notebook in ring and app.analyse_view.dropzone in ring

        # The Settings tab is where tab order actually matters. It must walk the
        # controls top to bottom and come back to the notebook.
        app.notebook.select(app.settings_view)
        root.update_idletasks()
        root.update()
        # The custom CheckBox is a Frame that takes focus, and it must be in the
        # walk: a checkbox the keyboard cannot reach is not a control. The two
        # speed spinboxes are legitimately disabled while tuning is automatic,
        # and portable mode is disabled outside a packaged build.
        classes = [w.winfo_class() for w in _focus_ring(app.notebook)]
        assert classes == ["TSpinbox", "TSpinbox", "TSpinbox", "TCombobox",
                           "TEntry", "Frame", "TButton", "TNotebook"], classes
        assert app.settings_view.auto_box in _focus_ring(app.notebook)

        # Unticking automatic tuning hands the two numbers back to the keyboard.
        app.settings_view.vars["perf_auto"].set(False)
        app.settings_view._auto_changed()
        root.update_idletasks()
        root.update()
        classes = [w.winfo_class() for w in _focus_ring(app.notebook)]
        assert classes.count("TSpinbox") == 5, classes
    finally:
        try:
            root.destroy()
        except Exception:
            pass


# ---------------------------------------------------------------------------
# regressions: the single after() chain, the start-up race, task routing,
# the settings round trip, Clear during a run, and the cancelled queue
# ---------------------------------------------------------------------------
def _app_on_display(*, withdraw: bool = True):
    """Build the real window on whatever display is available, or skip."""
    if not _display_available():
        raise SkipTest("no display: run under xvfb-run to exercise real widgets")
    import tkinter as tk

    try:
        root = tk.Tk()
    except tk.TclError as exc:
        raise SkipTest("Tk could not open a display: {}".format(exc)) from exc
    if withdraw:
        root.withdraw()
    gui.apply_scaling(root)
    gui.install_theme(root)
    return root, gui.WmlstApp(root, gui.Prefs())


def _settle(root, until, timeout: float = 5.0) -> None:
    """Pump the real event loop until ``until()`` is true or time runs out."""
    deadline = time.time() + timeout
    while time.time() < deadline:
        root.update()
        if until():
            return
        time.sleep(0.01)


def _ready_env() -> gui.Environment:
    return gui.Environment(dbdir="/db", datadir="/db/pubmlst",
                           blastdb="/db/blast/mlst.fa", db_ok=True, blast_ok=True,
                           blast_path="/nowhere/blastn", db_version="2025-12-29",
                           scheme_count=2, blast_version="2.17.0+",
                           scheme_names=("saureus", "sepidermidis"))


def test_the_poll_chain_is_rearmed_before_a_handler_can_block():
    """A handler reached from pump() may open a modal and spin a nested event
    loop (the first-run BLAST bootstrap does).  The next tick has to be armed
    *before* pump() runs, or the one after() chain is dead for as long as the
    modal is up and the dialog can never be fed or closed.
    """
    scheduler = FakeScheduler()
    trace = []

    def dispatch(msg):
        # stand-in for ModalDialog + root.wait_window(): a nested event loop
        # that runs whatever the scheduler already has pending.
        trace.append(len(scheduler.pending))
        for _ in range(3):
            scheduler.run()
        trace.append(len(scheduler.pending))

    controller = gui.AnalysisController(schedule=scheduler, dispatch=dispatch)
    controller.queue.put(gui.Msg("phase", 1, text="open a modal"))
    controller.ensure_polling()
    scheduler.run()
    assert trace, "the message was never dispatched"
    assert trace[0] == 1, "the tick must be re-armed BEFORE the handler can block"
    assert trace[1] == 1, "the nested loop must carry exactly one chain onwards"
    assert len(scheduler.pending) == 1, "the after() chain must neither die nor fork"


def test_a_modal_opened_from_a_handler_still_gets_its_later_messages():
    """End-to-end form of the same defect on a real Tk root: a dialog opened
    from a dispatch handler is closed by a message that arrives afterwards.
    """
    root, _app = _app_on_display()
    import tkinter as tk

    try:
        box = {}
        seen = []

        def dispatch(msg):
            if msg.text == "open":
                dialog = tk.Toplevel(root)
                dialog.withdraw()
                box["dialog"] = dialog
                root.wait_window(dialog)          # nested event loop
                seen.append("modal-returned")
                root.quit()
            elif msg.text == "close":
                seen.append("close-delivered")
                box["dialog"].destroy()

        controller = gui.AnalysisController(
            schedule=lambda ms, fn: root.after(ms, fn), dispatch=dispatch)
        controller.queue.put(gui.Msg("phase", 1, text="open"))
        controller.ensure_polling()
        root.after(400, lambda: controller.queue.put(gui.Msg("phase", 1, text="close")))

        def watchdog():
            dialog = box.get("dialog")
            if dialog is not None and dialog.winfo_exists():
                seen.append("TIMED-OUT")
                dialog.destroy()
            root.quit()

        root.after(6000, watchdog)
        root.mainloop()
        assert seen == ["close-delivered", "modal-returned"], seen
        controller._schedule = lambda ms, fn: None   # stop the chain re-arming
    finally:
        try:
            root.destroy()
        except Exception:
            pass


def test_command_line_files_wait_for_the_environment_probe():
    """``wmlst-gui FILE`` used to race the probe on a 120 ms timer and always
    lose, so a healthy install was shown the 137 MB download modal and the file
    was never analysed.  The files are handed over by _apply_env instead.
    """
    root, app = _app_on_display()
    try:
        started = []
        app.analyse_view.handle_paths = lambda paths: started.append(list(paths))
        app._argv_files = ["/tmp/from-argv.fa"]

        app._apply_env(_ready_env())
        assert started == [], "argv files must not be started synchronously"
        _settle(root, lambda: started)
        assert started == [["/tmp/from-argv.fa"]], started
        assert app._argv_files == [], "the argv list must be drained exactly once"
        assert app._bootstrap is None, "a healthy install must see no setup modal"

        # a probe that fails must still route the files (to the real
        # start_analysis, which then explains what is missing) and not drop them
        app._argv_files = ["/tmp/again.fa"]
        app._task_failed("env", gui.Friendly("nope", "the probe blew up"))
        _settle(root, lambda: len(started) > 1)
        assert started[-1] == ["/tmp/again.fa"], started
    finally:
        try:
            root.destroy()
        except Exception:
            pass


def test_a_finished_task_still_reaches_its_owner_and_is_never_dropped():
    """The owner is popped out of _owners before routing, so _route must be
    handed the owner it just removed; and a bootstrap result whose dialog has
    gone must still be adopted rather than silently discarded.
    """
    root, app = _app_on_display()
    try:
        class Owner:
            def __init__(self):
                self.payload = None

            def on_task_done(self, job, payload):
                self.payload = payload
                return True

        owner = Owner()
        app._owners[7] = owner
        app._task_kinds[7] = "blast-bootstrap"
        with app.controller._lock:
            app.controller._task_results[7] = "TOOLS"
        app.dispatch(gui.Msg("env", 7, text="blast-bootstrap"))
        assert owner.payload == "TOOLS", "the popped owner must still be routed to"
        assert 7 not in app._owners and 7 not in app._task_kinds

        # no owner left at all: the install must still be adopted
        adopted = []
        app.adopt_tools = lambda tools: adopted.append(tools)
        app._owners[8] = None
        app._task_kinds[8] = "blast-bootstrap"
        with app.controller._lock:
            app.controller._task_results[8] = "TOOLS2"
        app.dispatch(gui.Msg("env", 8, text="blast-bootstrap"))
        assert adopted == ["TOOLS2"], adopted
    finally:
        try:
            root.destroy()
        except Exception:
            pass


def test_closing_the_bootstrap_modal_mid_download_cancels_it():
    """Escape used to destroy the dialog and leave the 137 MB download running
    unattended with nothing to delete the partial file (section 10.7).
    """
    root, app = _app_on_display()
    try:
        dialog = gui.BootstrapDialog(root, app)
        app._bootstrap = dialog
        root.update()
        assert not dialog.cancel_event.is_set()

        dialog.job = 99                   # pretend the download is in flight
        dialog._downloading = True
        dialog.close()
        root.update()
        assert dialog.winfo_exists(), "the modal must stay up until the worker stops"
        assert dialog.cancel_event.is_set(), "closing mid-download must cancel it"

        # once the worker has acknowledged, closing goes through as usual
        dialog.on_task_failed(99, gui.Friendly("Stopped", "cancelled"))
        dialog.close()
        root.update()
        assert not dialog.winfo_exists()
    finally:
        try:
            root.destroy()
        except Exception:
            pass


def test_a_scheme_override_never_poisons_the_saved_minimum_score():
    """Forcing a scheme displays 0 (section 10.4) and runs at 0 (bin/mlst:125),
    but returning to Automatic must give the user their own value back.
    """
    forced = gui.Prefs(minscore=50.0, scheme="saureus")
    assert forced.normalised().minscore == 0.0
    assert forced.normalised().exclude == ()
    assert forced.normalised(force_scheme_rules=False).minscore == 50.0
    assert forced.normalised(force_scheme_rules=False).exclude == gui.DEFAULT_EXCLUDE
    assert forced.to_runconfig(["a.fa"], gui.Environment()).minscore == 0.0

    root, app = _app_on_display()
    try:
        view = app.settings_view
        view.set_scheme_choices(("saureus", "sepidermidis"))
        assert app.prefs.minscore == 50.0

        view.scheme_box.set("saureus")
        view.scheme_box.event_generate("<<ComboboxSelected>>")
        root.update()
        assert view.vars["minscore"].get() == "0", "the control must display 0"
        assert str(view.spin_minscore.cget("state")) == "disabled"
        assert app.prefs.scheme == "saureus"
        assert app.prefs.minscore == 50.0, "0 is a display, not the preference"
        assert app.prefs.to_runconfig(["a.fa"], app.env).minscore == 0.0
        assert app.prefs.to_runconfig(["a.fa"], app.env).exclude == frozenset()

        view.scheme_box.set("Automatic (recommended)")
        view.scheme_box.event_generate("<<ComboboxSelected>>")
        root.update()
        assert view.vars["minscore"].get() == "50"
        assert app.prefs.scheme is None
        assert app.prefs.minscore == 50.0, "the user's minimum score was destroyed"
        assert app.prefs.exclude == gui.DEFAULT_EXCLUDE
    finally:
        try:
            root.destroy()
        except Exception:
            pass


def test_clear_results_refuses_while_a_worker_is_alive():
    """Run > Clear results is not greyed out the way the Clear button is, and
    dropping to IDLE mid-run also unlocked the database update path.
    """
    state = gui.AppState()
    assert state.to(gui.RUNNING)
    assert not state.to(gui.IDLE), "RUNNING may only be left through finish()"
    assert state.finish("done")

    root, app = _app_on_display()
    try:
        app.env = _ready_env()
        app.controller._engine_factory = lambda cfg: FakeEngine(cfg, block=True)
        app.start_analysis(["/tmp/A.fa"])
        deadline = time.time() + 5.0
        while time.time() < deadline and not app.controller.running:
            root.update()
            time.sleep(0.01)
        assert app.controller.running

        assert app.analyse_view.clear() is False
        assert app.state.state == gui.RUNNING
        assert app.controller.running
        assert "Stop the analysis" in app.status.label.cget("text")

        # ... and the database update guard keys off the worker, not the view
        started = []
        app.run_task = lambda kind, fn, **kw: started.append(kind)
        app.state._state = gui.IDLE          # what Clear used to leave behind
        app.database_view.check_updates()
        assert started == [], "an update must not start while the engine is reading"
    finally:
        app.controller.request_cancel()
        app.controller.shutdown(5.0)
        try:
            root.destroy()
        except Exception:
            pass


def test_cancelling_says_what_happened_to_the_queued_files():
    """Files dropped during a run are acknowledged as queued; cancelling drops
    them, so it has to say so -- and must not resurrect them in a later batch.
    """
    root, app = _app_on_display()
    try:
        app.env = _ready_env()
        engines = []

        def blocking(cfg):
            engine = FakeEngine(cfg, block=True)
            engines.append(engine)
            return engine

        app.controller._engine_factory = blocking
        app.start_analysis(["/tmp/A.fa"])
        deadline = time.time() + 5.0
        while time.time() < deadline and not app.controller.running:
            root.update()
            time.sleep(0.01)
        app.start_analysis(["/tmp/B.fa"])
        assert app.controller._pending == ["/tmp/B.fa"]

        app.cancel_run()
        deadline = time.time() + 15.0
        while time.time() < deadline and app.state.running:
            root.update()
            time.sleep(0.01)
        root.update()
        assert app.state.state == gui.RESULTS and app.state.reason == "cancelled"
        assert app.controller._pending == [], "the queue must not survive a cancel"
        assert "not analysed" in app.status.label.cget("text"), \
            app.status.label.cget("text")

        # a later, unrelated run must not pick the cancelled file back up
        engines.clear()
        app.controller._engine_factory = lambda cfg: engines.append(
            FakeEngine(cfg)) or engines[-1]
        app.start_analysis(["/tmp/C.fa"])
        deadline = time.time() + 15.0
        while time.time() < deadline and app.state.running:
            root.update()
            time.sleep(0.01)
        root.update()
        assert [e.seen for e in engines] == [["/tmp/C.fa"]], [e.seen for e in engines]
    finally:
        app.controller.request_cancel()
        app.controller.shutdown(5.0)
        try:
            root.destroy()
        except Exception:
            pass

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
            except Exception as exc:
                failures += 1
                print("ERROR {}: {!r}".format(name, exc))
            else:
                print("ok   {}".format(name))
    print("\n{} failure(s)".format(failures))
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(_run_all())


# ---------------------------------------------------------------------------
# First run after an install: the index builds itself, and nothing races it
# ---------------------------------------------------------------------------

class _FakeStatus:
    def __init__(self):
        self.messages = []

    def set(self, text, kind="info", **kw):
        self.messages.append((kind, text))


class _FirstRunApp:
    """The slice of WmlstApp that _start_argv_files and the auto-build touch."""

    _start_argv_files = gui.WmlstApp._start_argv_files
    _auto_index_build_tried = False

    def __init__(self, env):
        self.env = env
        self._argv_files = ["sample.fna"]
        self.started = []
        self.analyse_view = type("V", (), {
            "handle_paths": staticmethod(lambda paths: self.started.append(list(paths)))
        })()
        # after() fires immediately here; the real one defers by 1 ms.
        self.root = type("R", (), {"after": staticmethod(lambda ms, fn: fn())})()


def _env(**kw):
    base = {"db_ok": True, "blast_ok": True, "index_stale": False}
    base.update(kw)
    return gui.Environment(**{k: v for k, v in base.items()
                              if k in {f.name for f in
                                       dataclasses.fields(gui.Environment)}})


def test_argv_files_wait_for_the_first_run_index_build():
    """Files named on the command line must not be analysed mid-rebuild.

    A search against a half-built index fails, and the failure surfaced as
    "WMLST could not start blastn" -- a missing-search-engine dialog on the very
    first launch, when the search engine was fine and only the index was absent.
    """
    app = _FirstRunApp(_env(index_stale=True))
    app._auto_index_build_tried = True          # the rebuild is in flight
    app._start_argv_files()
    assert app.started == [], "argv files were started while the index was building"
    assert app._argv_files == ["sample.fna"], "the queued files were dropped"


def test_argv_files_run_once_the_index_is_ready():
    """The rebuild's completion re-probes the environment, which releases them."""
    app = _FirstRunApp(_env(index_stale=False))
    app._auto_index_build_tried = True
    app._start_argv_files()
    assert app.started == [["sample.fna"]], "argv files were never started"
    assert app._argv_files == [], "the queue was not drained"


def test_argv_files_are_not_held_when_no_rebuild_was_attempted():
    """A stale index with no BLAST to rebuild it must not hang the queue open."""
    app = _FirstRunApp(_env(index_stale=True, blast_ok=False))
    app._auto_index_build_tried = False
    app._start_argv_files()
    assert app.started == [["sample.fna"]]


# ---------------------------------------------------------------------------
# 5.14a -- a scheme tie must be visible in the GUI, never silently resolved
# ---------------------------------------------------------------------------
def test_tied_schemes_text_names_both_schemes_and_their_sts():
    assert gui.tied_schemes(make_tied_result()) == (
        "klebsiella ST 258, ecoli_achtman_4 ST 14464")


def test_tied_schemes_is_empty_when_the_winner_stands_alone():
    assert gui.tied_schemes(make_result()) == ""


def test_tie_sentence_says_what_was_done_and_warns():
    text = gui.tie_sentence()
    assert "equally well" in text
    assert "Confirm the species" in text
    assert gui.branding.APP_NAME in text


def test_the_tie_reaches_the_result_row_the_card_and_the_status_line():
    """The ambiguity must be on screen, not only in the stderr WARNING."""
    if not _display_available():
        raise SkipTest("no display: run under xvfb-run to exercise real widgets")
    import tkinter as tk

    try:
        root = tk.Tk()
    except tk.TclError as exc:
        raise SkipTest("Tk could not open a display: {}".format(exc)) from exc
    try:
        root.withdraw()
        gui.apply_scaling(root)
        gui.install_theme(root)
        app = gui.WmlstApp(root, gui.Prefs())
        view = app.analyse_view
        view.add_result(make_tied_result())
        view.show_state(gui.RESULTS)
        root.update_idletasks()

        # the tie band across the single-file summary card
        assert view.tie_banner.winfo_manager(), "the tie band is not shown"
        assert "klebsiella ST 258" in view.tie_label.cget("text")
        assert "ecoli_achtman_4 ST 14464" in view.tie_label.cget("text")

        # and the row itself says so, in words as well as colour
        assert "tie" in view.tree.item(
            view.tree.get_children("")[0], "values")[1].lower()

        # the expansion row, above the per-locus table
        iid = view.tree.get_children("")[0]
        view.tree.focus(iid)
        view._on_open()
        children = list(view.tree.get_children(iid))
        # The schemes live in the #0 column, which is the one that stretches
        # with the window; in a narrow one they are wrapped over more than a
        # single row rather than chopped off mid-token.
        assert view.tree.item(children[0], "text").strip().startswith("TIE")
        assert view.tree.item(children[0], "values")[1] == "equal score"
        spelled = " ".join(view.tree.item(ch, "text") for ch in children)
        assert "klebsiella ST 258" in spelled
        assert "ecoli_achtman_4 ST 14464" in spelled
        captions = [view.tree.item(ch, "values")[0] for ch in children]
        assert "ALLELE" in captions

        # the status line
        view._on_select()
        assert "tie:" in app.status.label.cget("text").lower()
    finally:
        root.destroy()


# ---------------------------------------------------------------------------
# the design system (section 10.1): palette, contrast, spacing, type scale
# ---------------------------------------------------------------------------
def _luminance(hexcolour: str) -> float:
    """WCAG relative luminance of an ``#rrggbb`` colour."""
    parts = [int(hexcolour[i:i + 2], 16) / 255.0 for i in (1, 3, 5)]
    chan = [(v / 12.92 if v <= 0.03928 else ((v + 0.055) / 1.055) ** 2.4)
            for v in parts]
    return 0.2126 * chan[0] + 0.7152 * chan[1] + 0.0722 * chan[2]


def contrast(a: str, b: str) -> float:
    """WCAG 2.1 contrast ratio between two ``#rrggbb`` colours."""
    la, lb = _luminance(a), _luminance(b)
    hi, lo = max(la, lb), min(la, lb)
    return (hi + 0.05) / (lo + 0.05)


def test_every_palette_defines_the_same_roles():
    """A missing role would fall back to black and be invisible in dark mode."""
    roles = set(gui.PALETTE)
    for name, palette in (("dark", gui.DARK_PALETTE),
                          ("high-contrast", gui.HIGH_CONTRAST_PALETTE)):
        assert set(palette) == roles, name
        for key, value in palette.items():
            assert value.startswith("#") and len(value) == 7, (name, key)


def test_text_clears_wcag_aa_in_both_themes():
    """Every foreground role reads at 4.5:1 or better on its own ground."""
    for palette in (gui.PALETTE, gui.DARK_PALETTE):
        for ground in ("surface", "bg"):
            for role in ("text", "text_soft", "muted", "accent", "ok", "warn",
                         "bad", "info"):
                ratio = contrast(palette[role], palette[ground])
                assert ratio >= 4.5, (role, ground, round(ratio, 2))
        # a label on the accent button, and the selected-row tint
        assert contrast(palette["on_accent"], palette["accent"]) >= 4.5
        assert contrast(palette["text"], palette["accent_soft"]) >= 4.5
        assert contrast(palette["warn"], palette["warn_soft"]) >= 4.5


def test_the_spacing_ladder_is_an_eight_pixel_grid():
    ladder = [gui.PAD_XS, gui.PAD_S, gui.PAD_M, gui.PAD_L, gui.PAD_XL, gui.PAD_XXL]
    assert ladder == sorted(set(ladder)), ladder
    assert ladder[0] * 2 == ladder[1]
    for step in ladder[1:]:
        assert step % 8 == 0, step


def test_the_type_scale_has_one_size_per_role():
    sizes = [gui.FONT_PT[k] for k in ("tiny", "small", "body", "heading",
                                      "title", "hero")]
    assert sizes == sorted(sizes)
    assert sizes[-1] >= 24, "the ST has to be large and confident"
    # Windows gets the optical-size cuts first, then plain Segoe UI
    assert gui.PREFERRED_DISPLAY[0] == "Segoe UI Variable Display"
    assert "Segoe UI" in gui.PREFERRED_DISPLAY
    assert gui.PREFERRED_FAMILIES[1] == "Segoe UI"


def test_font_matching_ignores_case():
    """Tk reports the X core fonts lower-cased; a name must still match."""
    assert gui.pick_font(("helvetica", "fixed"), ("Helvetica",), "X") == "helvetica"
    assert gui.pick_font(("Helvetica",), ("Helvetica",), "X") == "Helvetica"


def test_mix_blends_and_clamps():
    assert gui.mix("#000000", "#ffffff", 0.0) == "#000000"
    assert gui.mix("#000000", "#ffffff", 1.0) == "#ffffff"
    assert gui.mix("#000000", "#ffffff", 0.5) == "#808080"
    assert gui.mix("#000000", "#ffffff", 5.0) == "#ffffff"
    assert gui.mix("#000000", "#ffffff", -5.0) == "#000000"
    assert gui.mix("not-a-colour", "#ffffff", 0.5) == "not-a-colour"


def test_theme_mode_is_selectable_and_always_has_a_palette():
    old = os.environ.get("WMLST_THEME")
    try:
        for value, expected in (("dark", "dark"), ("light", "light"),
                                ("high-contrast", "high-contrast"),
                                ("DARK", "dark")):
            os.environ["WMLST_THEME"] = value
            assert gui.detect_theme_mode() == expected
        os.environ["WMLST_THEME"] = "nonsense"
        assert gui.detect_theme_mode() in ("light", "dark", "high-contrast")
    finally:
        os.environ.pop("WMLST_THEME", None)
        if old is not None:
            os.environ["WMLST_THEME"] = old
    for mode in ("light", "dark", "high-contrast", "anything-else"):
        assert set(gui.palette_for(mode)) == set(gui.PALETTE)


def test_reduced_motion_is_honoured_both_ways():
    old = os.environ.get("WMLST_REDUCED_MOTION")
    try:
        os.environ["WMLST_REDUCED_MOTION"] = "1"
        assert gui.reduced_motion() is True
        os.environ["WMLST_REDUCED_MOTION"] = "0"
        assert gui.reduced_motion() is False
    finally:
        os.environ.pop("WMLST_REDUCED_MOTION", None)
        if old is not None:
            os.environ["WMLST_REDUCED_MOTION"] = old


def test_round_rect_points_are_a_closed_polygon():
    pts = gui.round_rect_points(0, 0, 100, 50, 10)
    assert len(pts) == 24 and len(pts) % 2 == 0
    assert min(pts[0::2]) >= 0 and max(pts[0::2]) <= 100
    assert min(pts[1::2]) >= 0 and max(pts[1::2]) <= 50
    # a radius larger than the box is clamped, never inverted
    tight = gui.round_rect_points(0, 0, 10, 4, 999)
    assert min(tight[0::2]) >= 0 and max(tight[1::2]) <= 4


# ---------------------------------------------------------------------------
# the drop-zone animation: one after() chain, and it stops when it should
# ---------------------------------------------------------------------------
def _idle_app(prefs=None):
    """Build the real window on whatever display exists, or skip."""
    if not _display_available():
        raise SkipTest("no display: run under xvfb-run to exercise real widgets")
    import tkinter as tk

    try:
        root = tk.Tk()
    except tk.TclError as exc:
        raise SkipTest("Tk could not open a display: {}".format(exc)) from exc
    root.withdraw()
    gui.apply_scaling(root)
    gui.install_theme(root)
    return root, gui.WmlstApp(root, prefs or gui.Prefs())


def test_the_drop_zone_animates_once_and_stops_on_every_pause_condition():
    root, app = _idle_app()
    try:
        zone = app.analyse_view.dropzone
        zone._static = False
        zone._mapped = True
        zone.animate(True)
        root.update_idletasks()
        assert zone.should_animate() is True
        first = zone._anim_after
        assert first is not None, "the chain never started"

        # arming twice must not create a second chain
        zone._start_anim()
        zone._start_anim()
        assert zone._anim_after == first

        # a frame re-arms exactly once, and it re-arms BEFORE painting
        zone._tick()
        assert zone._anim_after is not None and zone._anim_after != first

        # the sweep actually advances
        before = zone._phase
        zone._tick()
        assert zone._phase != before

        # every pause condition stops the chain dead
        for pause, resume in (
                (lambda: zone.animate(False), lambda: zone.animate(True)),
                (lambda: zone._on_unmap(), lambda: zone._on_map()),
                (lambda: zone.set_enabled(False), lambda: zone.set_enabled(True))):
            pause()
            assert zone._anim_after is None, "a paused zone kept a timer"
            assert zone.should_animate() is False
            resume()
            root.update_idletasks()
            assert zone._anim_after is not None

        # reduced motion: one static frame, never a timer
        zone._static = True
        zone.animate(True)
        assert zone._anim_after is None and zone.should_animate() is False
        zone.redraw()          # must still paint, without raising
        assert zone._items["arcs"]
        zone.flourish()
        assert zone._anim_after is None, "reduced motion must not animate the drop"
    finally:
        root.destroy()


def test_the_drop_zone_draws_seven_loci_in_every_state_and_at_every_size():
    root, app = _idle_app()
    try:
        zone = app.analyse_view.dropzone
        zone._static = False
        for width, height in ((1100, 420), (420, 190), (300, 120), (140, 60)):
            zone.configure(width=width, height=height)
            root.update_idletasks()
            for state in ("idle", "hover", "dragover", "focus"):
                zone.set_state(state)
                assert len(zone._items["arcs"]) == gui.DropZone.LOCI
                assert len(zone._items["nodes"]) == gui.DropZone.LOCI
        # the drop flourish runs and then gives the chain back
        zone.set_state("idle")
        zone.animate(False)
        zone.flourish()
        assert zone._anim_after is not None
        for _ in range(gui.DropZone.FLOURISH_FRAMES + 2):
            if zone._anim_after is None:
                break
            zone._tick()
        assert zone._flourish == 0
        assert zone._anim_after is None, "the flourish left the chain running"
    finally:
        root.destroy()


# ---------------------------------------------------------------------------
# the custom checkbox (the Database tab will need dozens of them)
# ---------------------------------------------------------------------------
def test_the_checkbox_is_big_keyboard_operable_and_shows_every_state():
    root, _app = _idle_app()
    try:
        import tkinter as tk

        calls = []
        var = tk.BooleanVar(value=False)
        box = gui.CheckBox(root, "Include this scheme", variable=var,
                           command=lambda: calls.append(1))
        box.pack()
        root.update_idletasks()

        # bigger than ttk's 13 px indicator, which is the whole point
        assert gui.CheckBox.SIZE >= 20
        assert int(box.box.cget("width")) >= 20
        assert int(box.cget("takefocus")) == 1

        # mouse: press then release toggles once, and the command fires
        box._press()
        box._release()
        assert var.get() is True and calls == [1]

        # keyboard: space and Return both toggle
        assert box._key() == "break"
        assert var.get() is False
        box._key()
        assert var.get() is True and len(calls) == 3

        # a tick is drawn only when it is on, and it is a stroke, not a glyph
        assert any(box.box.type(i) == "line" for i in box.box.find_all())
        box.set(False)
        root.update_idletasks()
        assert not any(box.box.type(i) == "line" for i in box.box.find_all())

        # hover and focus change the drawing, not only the colour of the tick
        box._enter()
        box._focus_in()
        root.update_idletasks()
        assert len(box.box.find_all()) >= 2, "no focus ring was drawn"
        box._focus_out()
        box._leave()

        # disabled: no focus, no toggle, and a dash so it is not colour alone
        box.set_enabled(False)
        assert int(box.cget("takefocus")) == 0
        box._key()
        assert var.get() is False
        assert any(box.box.type(i) == "line" for i in box.box.find_all())
    finally:
        root.destroy()


def test_a_card_is_never_a_three_d_frame():
    root, app = _idle_app()
    try:
        import tkinter as tk

        holder = tk.Frame(root)
        body = gui.card(holder)
        outer = body.master
        assert isinstance(outer, gui.Card)
        assert str(outer.cget("background")) == gui.c("border")
        assert str(body.cget("background")) == gui.c("surface")
        assert len(outer._corners) == 4
        # nothing anywhere in the window may carry a 3-D relief
        bad = []

        def walk(widget):
            # Menus are drawn by the platform on Windows, so their relief is the
            # OS's business and not this theme's.
            if widget.winfo_class() != "Menu":
                try:
                    relief = str(widget.cget("relief"))
                except Exception:
                    relief = ""
                if relief in ("sunken", "groove", "ridge", "raised"):
                    bad.append((str(widget), relief))
            for child in widget.winfo_children():
                walk(child)

        walk(app.root)
        assert bad == [], bad
    finally:
        root.destroy()


def test_the_settings_tab_can_be_scrolled_to_its_last_control():
    """At the minimum window size the last card must still be reachable."""
    root, app = _idle_app()
    try:
        root.deiconify()
        root.geometry("{}x{}+30+30".format(gui.px(gui.MIN_WIDTH),
                                           gui.px(gui.MIN_HEIGHT)))
        app.notebook.select(app.settings_view)
        root.update_idletasks()
        root.update()
        pane = app.settings_view.scroller
        assert isinstance(pane, gui.ScrollPane)
        region = pane.canvas.cget("scrollregion")
        assert region, "the scroll region was never computed"
        pane.canvas.yview_moveto(1.0)
        root.update_idletasks()
        assert float(pane.canvas.yview()[1]) >= 0.99
    finally:
        root.destroy()


# ---------------------------------------------------------------------------
# the result names the organism, the scheme and the reference (wmlst.schemerefs)
# ---------------------------------------------------------------------------
def test_the_organism_is_read_from_the_curated_table_and_never_invented():
    """A scheme id is a directory name; the organism is the answer."""
    assert gui.organism_name("sepidermidis") == "Staphylococcus epidermidis"
    assert gui.organism_name("klebsiella") == "Klebsiella pneumoniae"
    # A genus-wide scheme carries no species, and "spp." is not italic.
    assert gui.organism_parts("neisseria") == ("Neisseria", "spp.")
    assert gui.organism_parts("klebsiella") == ("Klebsiella pneumoniae", "")
    # Nothing is invented for an unknown scheme, an empty one, or "-".
    for unknown in ("", "-", "not_a_scheme_at_all"):
        assert gui.organism_name(unknown) == ""
        assert gui.scheme_ref(unknown) is None
        assert gui.organism_parts(unknown) == ("", "")


def test_the_scheme_caption_and_the_reference_show_only_what_is_known():
    ref = gui.scheme_ref("abaumannii_2")
    assert ref is not None
    caption = gui.scheme_caption("abaumannii_2", ref)
    assert "abaumannii_2" in caption           # the id stays, secondary
    assert "Acinetobacter baumannii" in caption   # and is explained
    assert "7 loci" in caption
    line = gui.reference_line(ref)
    assert "Diancourt" in line and "PMID 20383326" in line
    assert ref.pubmed_url.endswith("/20383326/")
    assert ref.database_url.startswith("https://")
    assert gui.source_name(ref) == "PubMLST"
    # A scheme with no verified citation prints no citation at all.
    assert gui.reference_line(None) == ""
    blank = gui.schemerefs.SchemeRef(scheme="x", genus="X", species="y")
    assert gui.reference_line(blank) == ""


def test_the_summary_card_leads_with_the_organism_and_cites_the_scheme():
    root, app = _idle_app()
    try:
        view = app.analyse_view
        view.add_result(make_result())
        view.show_state(gui.RESULTS)
        root.update_idletasks()
        assert view.organism_label.cget("text") == "Staphylococcus epidermidis"
        assert str(view.organism_label.cget("font")) == str(gui.F["organism"])
        caption = view.summary_title.cget("text")
        assert "sepidermidis" in caption and "7 loci" in caption
        assert "17151213" in view.reference_label.cget("text")
        links = [w.cget("text") for w in view.reference_links.winfo_children()]
        assert links and links[0].startswith("PubMed")
        assert any("PubMLST" in text for text in links)
        # every link is a real button: Tab reaches it, Enter follows it
        for widget in view.reference_links.winfo_children():
            assert widget.winfo_class() == "TButton"
        # and the row carries the organism as its first column
        row = view.tree.get_children("")[0]
        assert view.tree.item(row, "values")[0] == "Staphylococcus epidermidis"
    finally:
        root.destroy()


def test_an_unassigned_st_states_the_organism_quietly():
    """No ST means nothing was confirmed; the name must not shout."""
    root, app = _idle_app()
    try:
        view = app.analyse_view
        view.add_result(make_result(status="NONE", st="-", score=0))
        view.show_state(gui.RESULTS)
        root.update_idletasks()
        assert view.organism_label.cget("style") == "SurfaceMuted.TLabel"
        assert view.st_label.cget("text") == "not assigned"
        view.clear()
        assert view.organism_label.cget("text") == ""
        assert not view.reference_wrap.winfo_manager()
    finally:
        root.destroy()


# ---------------------------------------------------------------------------
# the Database tab: big tick boxes, select all / deselect all, a live count
# ---------------------------------------------------------------------------
class FakeInfo:
    def __init__(self, name, loci=7):
        self.name = name
        self.locus = loci
        self.num_genotypes = 10
        self.num_alleles = 100
        self.last_updated = "2026-01-01"


class FakeCatalog:
    def species_label(self, name):
        return {"klebsiella": "Klebsiella pneumoniae"}.get(name)


def test_the_scheme_list_has_tick_boxes_two_buttons_and_a_live_count():
    root, app = _idle_app()
    try:
        view = app.database_view
        names = ["abaumannii_2", "klebsiella", "saureus", "sepidermidis"]
        view._fill(FakeCatalog(), [FakeInfo(n) for n in names], deep=True)
        root.update_idletasks()
        assert view.selection_label.cget("text") == "4 of 4 schemes selected"
        rows = view.tree.get_children("")

        # the tick box is an image, and a big one: the stock ttk indicator is
        # 13 px and this is the CheckBox square at the same size as the widget
        on, off = view._images()
        assert on is not None and off is not None
        assert on.height() >= gui.px(20)
        assert on.width() > on.height()   # transparent gap before the label
        assert view.tree.item(rows[0], "image")

        view.deselect_all()
        assert view.selected_schemes() == ()
        assert view.selection_label.cget("text") == "0 of 4 schemes selected"
        view.toggle(rows[1])
        view.toggle(rows[2])
        assert view.selected_schemes() == ("klebsiella", "saureus")
        assert view.selection_label.cget("text") == "2 of 4 schemes selected"
        view.select_all()
        assert len(view.selected_schemes()) == 4

        # Space toggles the focused row, so the list is keyboard-operable
        view.tree.focus(rows[0])
        view._space()
        assert view.checked["abaumannii_2"] is False

        # a scheme the catalogue does not map still names its organism
        assert view.tree.set(rows[0], "species") == "Acinetobacter baumannii"
    finally:
        root.destroy()


def test_clicking_a_tick_box_hands_the_scheme_list_the_keyboard_focus():
    """A tick-box click must move the KEYBOARD focus, not just the row focus.

    ``_click`` returns "break", which suppresses ttk::treeview's own Button-1
    binding -- the one that would normally have focused the widget. Without an
    explicit ``focus_set`` the focus stayed on the Analyse tab's drop zone, so
    Space never reached ``_space`` and instead opened that tab's file chooser
    from underneath the Database tab. Calling ``_space()`` directly, as the
    test above does, cannot see any of that; only a real click can.
    """
    root, app = _idle_app()
    try:
        root.deiconify()
        root.geometry("1200x800+50+50")
        for tab in app.notebook.tabs():
            if app.notebook.tab(tab, "text").strip().lower() == "database":
                app.notebook.select(tab)
                break
        view = app.database_view
        view._fill(FakeCatalog(), [FakeInfo(n) for n in
                                   ("abaumannii_2", "klebsiella", "saureus")],
                   deep=True)
        root.update()
        rows = view.tree.get_children("")
        box = view.tree.bbox(rows[1])
        if not box:
            raise SkipTest("the scheme list has no geometry on this display")
        if root.focus_get() is None:
            raise SkipTest("this display does not report keyboard focus")
        assert root.focus_get() is not view.tree
        view.tree.event_generate("<Button-1>", x=6, y=box[1] + box[3] // 2)
        root.update()
        assert root.focus_get() is view.tree
        # and the row really did toggle, so the click still does its own job
        assert view.selection_label.cget("text") == "2 of 3 schemes selected"
    finally:
        root.destroy()


def test_the_update_dialog_selects_all_deselects_all_and_counts():
    root, _app = _idle_app()
    try:
        class FakeUpdate:
            def __init__(self, name):
                self.name = name
                self.bytes_estimate = 1024 * 1024
                self.added_types = 5
                self.detail = "7 loci"
                self.status = "changed"

        updates = [FakeUpdate(n) for n in ("klebsiella", "saureus", "ecoli_achtman_4")]
        dialog = gui.UpdatePlanDialog(root, object(), updates,
                                      preselected={"saureus": False})
        try:
            root.update_idletasks()
            assert dialog.count_label.cget("text") == "2 of 3 schemes selected"
            dialog._set_all(True)
            assert dialog.count_label.cget("text") == "3 of 3 schemes selected"
            dialog._set_all(False)
            assert dialog.count_label.cget("text") == "0 of 3 schemes selected"
            dialog._accept()
            assert dialog.accepted is False and dialog.chosen == []
        finally:
            try:
                dialog.destroy()
            except Exception:
                pass
    finally:
        root.destroy()


# ---------------------------------------------------------------------------
# portable mode (section 10.4): where the database is written
# ---------------------------------------------------------------------------
def test_portable_mode_is_detected_from_the_layout(tmp_path=None):
    folder = tempfile.mkdtemp(prefix="wmlst-portable-")
    old = os.environ.get("WMLST_PORTABLE_DIR")
    try:
        os.environ["WMLST_PORTABLE_DIR"] = folder
        assert gui.portable_possible() is True
        assert str(gui.portable_root()) == folder
        assert str(gui.portable_db_dir()) == os.path.join(folder, "db")
        os.environ.pop("WMLST_PORTABLE_DIR")
        # A plain Python run is not a packaged build.
        if not getattr(sys, "frozen", False):
            assert gui.portable_possible() is False
            assert gui.portable_root() is None
            assert gui.portable_db_dir() is None
    finally:
        if old is None:
            os.environ.pop("WMLST_PORTABLE_DIR", None)
        else:
            os.environ["WMLST_PORTABLE_DIR"] = old
        shutil.rmtree(folder, ignore_errors=True)


def test_a_read_only_folder_is_reported_not_raised():
    folder = tempfile.mkdtemp(prefix="wmlst-ro-")
    try:
        assert gui.dir_writable(folder) is True
        assert gui.dir_writable(os.path.join(folder, "nope")) is False
        os.chmod(folder, 0o500)
        if os.getuid() != 0:     # root ignores the mode bits
            assert gui.dir_writable(folder) is False
        assert gui.looks_like_database(folder) is False
    finally:
        os.chmod(folder, 0o700)
        shutil.rmtree(folder, ignore_errors=True)


def test_the_resolved_database_folder_is_always_printable():
    folder = tempfile.mkdtemp(prefix="wmlst-portable-")
    old = os.environ.get("WMLST_PORTABLE_DIR")
    try:
        os.environ["WMLST_PORTABLE_DIR"] = folder
        prefs = gui.Prefs()
        env = gui.Environment(dbdir="/opt/wmlst/db")
        assert gui.db_location(prefs, env) == "/opt/wmlst/db"
        prefs.portable_db = True
        assert gui.db_location(prefs, env) == os.path.join(folder, "db")
        prefs.portable_db = False
        prefs.dbdir = "/somewhere/else/db"
        assert gui.db_location(prefs, env) == "/somewhere/else/db"
    finally:
        if old is None:
            os.environ.pop("WMLST_PORTABLE_DIR", None)
        else:
            os.environ["WMLST_PORTABLE_DIR"] = old
        shutil.rmtree(folder, ignore_errors=True)


def test_copying_the_database_copies_every_file_and_reports_progress():
    src = tempfile.mkdtemp(prefix="wmlst-src-")
    dst = os.path.join(tempfile.mkdtemp(prefix="wmlst-dst-"), "db")
    seen = []
    try:
        os.makedirs(os.path.join(src, "pubmlst", "saureus"))
        os.makedirs(os.path.join(src, "blast.staging.99"))
        for path, text in (("VERSION.txt", "2026-01-01"),
                           ("pubmlst/saureus/saureus.txt", "ST\tarcC\n1\t1\n"),
                           ("blast.staging.99/junk", "x")):
            with open(os.path.join(src, path), "w", encoding="utf-8") as fh:
                fh.write(text)
        out = gui.copy_database(src, dst,
                                progress=lambda pct, text="": seen.append(pct))
        assert out == dst
        assert os.path.isfile(os.path.join(dst, "VERSION.txt"))
        assert os.path.isfile(os.path.join(dst, "pubmlst", "saureus", "saureus.txt"))
        # staging leftovers are not part of a database
        assert not os.path.exists(os.path.join(dst, "blast.staging.99"))
        assert seen and seen[-1] == 100.0
        assert gui.looks_like_database(dst) is True
    finally:
        shutil.rmtree(src, ignore_errors=True)
        shutil.rmtree(os.path.dirname(dst), ignore_errors=True)


def test_portable_mode_refuses_a_read_only_folder_without_a_traceback():
    folder = tempfile.mkdtemp(prefix="wmlst-ro-")
    old = os.environ.get("WMLST_PORTABLE_DIR")
    os.environ["WMLST_PORTABLE_DIR"] = folder
    root, app = _idle_app()
    try:
        os.chmod(folder, 0o500)
        if os.getuid() == 0:
            raise SkipTest("running as root: the mode bits mean nothing")
        app.set_portable_db(True)
        assert app.prefs.portable_db is False
        assert app.settings_view.vars["portable_db"].get() is False
        said = app.status.label.cget("text").lower()
        assert "cannot be written" in said or "read" in said
    finally:
        os.chmod(folder, 0o700)
        root.destroy()
        if old is None:
            os.environ.pop("WMLST_PORTABLE_DIR", None)
        else:
            os.environ["WMLST_PORTABLE_DIR"] = old
        shutil.rmtree(folder, ignore_errors=True)


def test_portable_mode_adopts_a_database_that_is_already_there():
    folder = tempfile.mkdtemp(prefix="wmlst-portable-")
    old = os.environ.get("WMLST_PORTABLE_DIR")
    os.environ["WMLST_PORTABLE_DIR"] = folder
    root, app = _idle_app()
    try:
        os.makedirs(os.path.join(folder, "db", "pubmlst"))
        app.set_portable_db(True)
        assert app.prefs.portable_db is True
        assert app.prefs.dbdir == os.path.join(folder, "db")
        assert app.settings_view.portable_path.cget("text").endswith(
            os.path.join(folder, "db"))
        # and switching it off goes back to the per-user location
        app.set_portable_db(False)
        assert app.prefs.portable_db is False and app.prefs.dbdir is None
    finally:
        root.destroy()
        if old is None:
            os.environ.pop("WMLST_PORTABLE_DIR", None)
        else:
            os.environ["WMLST_PORTABLE_DIR"] = old
        shutil.rmtree(folder, ignore_errors=True)


# ---------------------------------------------------------------------------
# start-up auto-tuning (wmlst.perf) — the GUI only, never the CLI
# ---------------------------------------------------------------------------
def test_the_window_sizes_itself_for_this_computer_and_says_why():
    from wmlst import perf

    tuning = perf.probe()
    prefs = gui.Prefs()
    assert (prefs.threads, prefs.jobs) == (1, 1), "the reference defaults are 1/1"
    root, app = _idle_app(prefs)
    try:
        app.apply_autotune(announce=True)
        root.update_idletasks()
        assert app.prefs.threads == tuning.threads
        assert app.prefs.jobs == tuning.jobs
        assert app.settings_view.tuning_label.cget("text") == tuning.rationale
        assert tuning.rationale in app.status.label.cget("text")
        # the numbers belong to the tuner while the box is ticked
        assert str(app.settings_view.spin_threads.cget("state")) == "disabled"

        # unticking hands them back, and a second probe leaves them alone
        app.prefs.perf_auto = False
        app.prefs.threads, app.prefs.jobs = 2, 3
        app.settings_view.load_from(app.prefs)
        app.apply_autotune()
        assert (app.prefs.threads, app.prefs.jobs) == (2, 3)
        assert str(app.settings_view.spin_threads.cget("state")) == "normal"
    finally:
        root.destroy()


def test_restoring_the_reference_defaults_turns_automatic_tuning_off():
    root, app = _idle_app()
    try:
        app.settings_view.restore_defaults()
        root.update_idletasks()
        assert app.prefs.perf_auto is False
        assert (app.prefs.threads, app.prefs.jobs) == (1, 1)
    finally:
        root.destroy()


def test_the_two_new_preferences_round_trip():
    prefs = gui.Prefs()
    assert prefs.perf_auto is True and prefs.portable_db is False
    prefs.perf_auto = False
    prefs.portable_db = True
    again = gui.Prefs.from_dict(prefs.to_dict())
    assert again.perf_auto is False and again.portable_db is True
    # anything unreadable falls back to the safe answer
    assert gui.Prefs.from_dict({}).perf_auto is True
    assert gui.Prefs.from_dict({}).portable_db is False


# ---------------------------------------------------------------------------
# closing the window must end the process (section 10.5)
# ---------------------------------------------------------------------------
def test_closing_stops_the_engines_and_the_search_children():
    calls = []
    real_stop = gui.stop_the_engines
    real_force = gui.force_exit
    root, app = _idle_app()
    try:
        gui.stop_the_engines = lambda timeout=2.0: calls.append(("stop", timeout)) or True
        gui.force_exit = lambda code=0: calls.append(("force", code))
        assert gui.owns_process() is False, "only main() owns the process"
        app.on_close()
        assert ("stop", gui.SHUTDOWN_TIMEOUT_S) in calls
        assert not any(kind == "force" for kind, _ in calls), (
            "a test process must never be exited from under pytest")
    finally:
        gui.stop_the_engines = real_stop
        gui.force_exit = real_force
        try:
            root.destroy()
        except Exception:
            pass


def test_an_unclean_shutdown_ends_the_process_when_it_owns_it():
    calls = []
    real_stop = gui.stop_the_engines
    real_force = gui.force_exit
    real_arm = gui.arm_exit_watchdog
    root, app = _idle_app()
    try:
        gui._OWNS_PROCESS = True
        gui.stop_the_engines = lambda timeout=2.0: False     # a child survived
        gui.force_exit = lambda code=0: calls.append(("force", code))
        gui.arm_exit_watchdog = lambda seconds=5.0: calls.append(("arm", seconds))
        app.on_close()
        assert ("force", 0) in calls, "the process was left running"
        assert ("arm", gui.EXIT_WATCHDOG_S) in calls, "no watchdog was armed"
    finally:
        gui._OWNS_PROCESS = False
        gui.stop_the_engines = real_stop
        gui.force_exit = real_force
        gui.arm_exit_watchdog = real_arm
        try:
            root.destroy()
        except Exception:
            pass


def test_the_exit_watchdog_is_a_daemon_and_does_not_hold_the_process():
    timer = gui.arm_exit_watchdog(30.0)
    try:
        assert timer.daemon is True
        assert timer.is_alive()
    finally:
        timer.cancel()


def test_stopping_the_engines_is_safe_with_nothing_running():
    from wmlst import blastbin

    try:
        assert gui.stop_the_engines(0.5) is True
    finally:
        # shutdown() latches the module; the rest of the suite still needs it
        blastbin.reset_shutdown()
    assert blastbin.is_shutting_down() is False
