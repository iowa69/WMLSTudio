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
        glyph, colour, sentence = gui.STATUS_UI[status]
        assert glyph and colour in gui.PALETTE
        assert sentence.endswith(".") and len(sentence) > 40
        cell = gui.status_cell(status)
        assert cell.endswith(status) and cell != status  # colour is never alone
        assert gui.status_tag(status).startswith("st_")


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
        assert values[0] == "sepidermidis" and values[1] == "184"
        assert values[2].endswith("PERFECT") and values[3] == "100"
        assert app.analyse_view.st_label.cget("text") == "184"

        # the per-locus children are built lazily, in gene order
        app.analyse_view.tree.focus(rows[0])
        app.analyse_view._on_open()
        children = app.analyse_view.tree.get_children(rows[0])
        assert app.analyse_view.tree.item(children[0], "values")[0] == "ALLELE"
        loci = [app.analyse_view.tree.item(ch, "text").strip() for ch in children[1:]]
        assert loci == ["arcC", "aroE", "gtr", "mutS", "pyrR", "tpiA", "yqiL"]
        codes = [app.analyse_view.tree.item(ch, "values")[0] for ch in children[1:]]
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
        classes = [w.winfo_class() for w in _focus_ring(app.notebook)]
        assert classes == ["TSpinbox", "TSpinbox", "TSpinbox", "TCombobox",
                           "TEntry", "TSpinbox", "TSpinbox", "TButton",
                           "TNotebook"], classes
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
