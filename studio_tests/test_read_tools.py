"""fastp adapter tests use an explicit process stub; a genuine run is an opt-in.

A stub proves the adapter's contract. It never proves that fastp trims anything,
and a Linux run never validates the Windows executable.
"""

import bz2
import gzip
import hashlib
import importlib.util
import io
import json
import os
import sys
import tarfile
import threading
from pathlib import Path

import pytest

from wmlstudio import read_tools
from wmlstudio.project import Project
from wmlstudio.sequence import AnalysisCancelled, SequenceError

_spec = importlib.util.spec_from_file_location(
    "stage_read_tools", Path(__file__).resolve().parents[1] / "studio_packaging/stage_read_tools.py")
staging = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(staging)


def reads(tmp_path, count=4, stem="reads, λ", length=12):
    Path(tmp_path).mkdir(parents=True, exist_ok=True)
    paths = [tmp_path / f"{stem}_R1.fastq", tmp_path / f"{stem}_R2.fastq"]
    sequence = ("ACGTTGCAACGT" * (length // 12 + 1))[:length]
    for mate, path in enumerate(paths, 1):
        path.write_text("".join(f"@read{index}/{mate}\n{sequence}\n+\n{'I' * length}\n"
                                for index in range(count)), encoding="ascii")
    return paths


def gzipped(path, records):
    with gzip.open(path, "wb") as handle:
        handle.write("".join(f"@trimmed{index}/1\nACGTTGCA\n+\nIIIIIIII\n"
                             for index in range(records)).encode("ascii"))


def fastp_json(directory, *, before_reads=1000, after_reads=900, version=read_tools.VERSION):
    report = {
        "summary": {
            "fastp_version": version, "sequencing": "paired end (100 cycles + 100 cycles)",
            "before_filtering": {"total_reads": before_reads, "total_bases": before_reads * 100,
                                 "q20_bases": before_reads * 90, "q30_bases": before_reads * 80,
                                 "q20_rate": 0.9, "q30_rate": 0.8, "gc_content": 0.51,
                                 "read1_mean_length": 100, "read2_mean_length": 100},
            "after_filtering": {"total_reads": after_reads, "total_bases": after_reads * 98,
                                "q20_bases": after_reads * 95, "q30_bases": after_reads * 88,
                                "q20_rate": 0.97, "q30_rate": 0.9, "gc_content": 0.51,
                                "read1_mean_length": 98, "read2_mean_length": 98}},
        "filtering_result": {"passed_filter_reads": after_reads,
                             "low_quality_reads": before_reads - after_reads, "too_many_N_reads": 0,
                             "adapter_dimer_reads": 0, "too_short_reads": 0, "too_long_reads": 0},
        "duplication": {"rate": 0.031},
        "insert_size": {"peak": 210, "unknown": 12},
        "adapter_cutting": {"adapter_trimmed_reads": 42, "adapter_trimmed_bases": 900,
                            "read1_adapter_sequence": "AGATCGGAAGAGC",
                            "read2_adapter_sequence": "AGATCGGAAGAGC"},
    }
    (Path(directory) / "fastp.json").write_text(json.dumps(report), encoding="utf-8")
    return report


@pytest.fixture
def process_stub(tmp_path, monkeypatch):
    """Only adapter contract tests; this fixture does not validate a read trimmer."""
    root = tmp_path / "tools"
    root.mkdir()
    binary = root / ("fastp.exe" if os.name == "nt" else "fastp")
    binary.write_bytes(b"adapter fixture, not a trimmer")
    monkeypatch.setattr(read_tools, "runtime_capabilities", lambda root=None: {
        "available": True, "binary": str(binary), "expected_sha256": "stub-only",
        "version": read_tools.VERSION, "platform": "test-only",
        "source_commit": read_tools.SOURCE_COMMIT, "reason": ""})
    monkeypatch.setattr(read_tools, "_tool_info", lambda capability, directory, cancelled: {
        "name": "fastp", "version": read_tools.VERSION, "executable": str(binary),
        "executable_sha256": "stub-only", "source_commit": read_tools.SOURCE_COMMIT,
        "platform": "test-only", "reported_version": f"fastp {read_tools.VERSION}"})
    state = {"after_reads": 900, "unpaired": 0, "exit": 0, "calls": []}

    def run(command, directory, log_path, cancelled=None, progress=None, timeout=None):
        state["calls"].append(command)
        plan = json.loads((directory / "run-plan.json").read_text())
        assert plan["state"] == "prepared_not_completed"
        assert plan["inputs"][0]["sha256"] and plan["inputs"][1]["sha256"]
        Path(log_path).write_text("TEST ONLY subprocess fixture\n")
        fastp_json(directory, after_reads=state["after_reads"])
        (directory / "fastp.html").write_text("<h1>TEST ONLY fastp report fixture</h1>")
        if state["exit"]:
            return state["exit"]
        for name in ("trimmed_R1.fastq.gz", "trimmed_R2.fastq.gz"):
            gzipped(directory / name, max(1, state["after_reads"] // 2))
        for name in ("unpaired_R1.fastq.gz", "unpaired_R2.fastq.gz"):
            gzipped(directory / name, state["unpaired"])
        return 0

    monkeypatch.setattr(read_tools, "_run", run)
    return state


def trim(tmp_path, state, **options):
    paths = reads(tmp_path)
    return paths, read_tools.run_fastp(*paths, tmp_path / "trimmed", **options)


def test_trimming_publishes_the_tools_own_report_and_never_touches_the_originals(tmp_path, process_stub):
    paths = reads(tmp_path)
    original = [path.read_bytes() for path in paths]
    messages = []
    result = read_tools.run_fastp(*paths, tmp_path / "trimmed",
                                  progress=lambda done, total, message: messages.append(message))
    assert [path.read_bytes() for path in paths] == original
    assert result["report"]["before"]["total_reads"] == 1000
    assert result["report"]["after"]["total_reads"] == 900
    assert result["report"]["reads_removed"] == 100
    assert result["report"]["reads_removed_percent"] == pytest.approx(10.0)
    assert result["report"]["duplication_rate"] == 0.031
    assert result["report"]["adapter"]["trimmed_reads"] == 42
    link = result["provenance"]["link"]
    assert [item["original_sha256"] for item in link] == [hashlib.sha256(data).hexdigest() for data in original]
    assert all(Path(item["trimmed_path"]).is_file() for item in link)
    assert [item["trimmed_path"] for item in link] == [result["read1_path"], result["read2_path"]]
    assert all(Path(item["original_path"]) in paths for item in link)
    assert json.loads((tmp_path / "trimmed/read-trimming.json").read_text()) == result
    assert (tmp_path / "trimmed/fastp.json").is_file() and (tmp_path / "trimmed/fastp.html").is_file()
    assert not list(tmp_path.glob(".wmlstudio-trim-*"))
    assert any("Cancel" in message or "Validating" in message for message in messages)


def test_trimming_never_claims_it_validated_the_isolate(tmp_path, process_stub):
    _paths, result = trim(tmp_path, process_stub)
    assert read_tools.TRIMMING_DISCLAIMER in result["notes"]
    assert any("fastp's own output" in note for note in result["notes"])
    assert any("the originals still hold them" in note for note in result["notes"])
    assert "Phred+33" in result["provenance"]["quality_encoding"]


def test_an_empty_unpaired_file_is_not_reported_as_a_rescued_read(tmp_path, process_stub):
    _paths, result = trim(tmp_path, process_stub)
    assert result["provenance"]["unpaired_outputs"] == []
    assert not (tmp_path / "trimmed/unpaired_R1.fastq.gz").exists()
    assert not any("unpaired" in note for note in result["notes"])
    process_stub["unpaired"] = 3
    rescued = read_tools.run_fastp(*reads(tmp_path / "again"), tmp_path / "second")
    assert [item["mate"] for item in rescued["provenance"]["unpaired_outputs"]] == [1, 2]
    assert (tmp_path / "second/unpaired_R2.fastq.gz").is_file()
    assert any("unpaired" in note for note in rescued["notes"])


def test_a_run_that_keeps_no_reads_publishes_nothing_but_keeps_its_evidence(tmp_path, process_stub):
    process_stub["after_reads"] = 0
    with pytest.raises(read_tools.ReadTrimError, match="kept none"):
        trim(tmp_path, process_stub)
    assert not (tmp_path / "trimmed").exists()
    plan = json.loads(next(tmp_path.glob("trimmed.failed-*.json")).read_text())
    assert plan["state"] == "failed"
    assert plan["report"]["before"]["total_reads"] == 1000
    assert not list(tmp_path.glob(".wmlstudio-trim-*"))


def test_a_failing_process_preserves_its_log_without_publishing_partial_output(tmp_path, process_stub):
    process_stub["exit"] = 9
    with pytest.raises(read_tools.ReadTrimError, match="code 9"):
        trim(tmp_path, process_stub)
    assert not (tmp_path / "trimmed").exists()
    assert "TEST ONLY subprocess fixture" in next(tmp_path.glob("trimmed.failed-*.log")).read_text()


def test_cancellation_before_work_publishes_nothing(tmp_path, process_stub):
    with pytest.raises(AnalysisCancelled):
        trim(tmp_path, process_stub, cancelled=lambda: True)
    assert not (tmp_path / "trimmed").exists()
    assert not process_stub["calls"]


def test_a_real_subprocess_is_cancellable_and_keeps_its_log(tmp_path):
    event = threading.Event()
    timer = threading.Timer(0.5, event.set)
    timer.start()
    try:
        with pytest.raises(AnalysisCancelled):
            read_tools._run([sys.executable, "-u", "-c",
                             "import time; print('started', flush=True); time.sleep(60)"],
                            tmp_path, tmp_path / "child.log", event.is_set)
    finally:
        timer.cancel()
    assert "started" in (tmp_path / "child.log").read_text()


@pytest.mark.parametrize("mutation,error", [
    ("same_file", "different files"),
    ("fasta", "two FASTQ"),
    ("bzip2", "bzip2"),
])
def test_unsuitable_inputs_are_refused_before_any_process_runs(tmp_path, process_stub, mutation, error):
    left, right = reads(tmp_path)
    if mutation == "same_file":
        os.link(left, tmp_path / "same.fastq")
        right = tmp_path / "same.fastq"
    elif mutation == "fasta":
        right.write_text(">contig\nACGT\n")
    else:
        right = tmp_path / "compressed.fastq.bz2"
        right.write_bytes(bz2.compress(b"@r/2\nACGT\n+\nIIII\n"))
    with pytest.raises(SequenceError, match=error):
        read_tools.run_fastp(left, right, tmp_path / "trimmed")
    assert not process_stub["calls"]
    assert not (tmp_path / "trimmed").exists()


def test_an_existing_output_directory_is_never_replaced(tmp_path, process_stub):
    left, right = reads(tmp_path)
    (tmp_path / "trimmed").mkdir()
    with pytest.raises(FileExistsError):
        read_tools.run_fastp(left, right, tmp_path / "trimmed")
    assert not process_stub["calls"]


@pytest.mark.parametrize("options", [
    {"threads": 0}, {"threads": True}, {"min_length": 0}, {"quality_threshold": 61},
    {"unqualified_percent_limit": 101}, {"compression": 0}, {"detect_adapter": "yes"},
])
def test_out_of_range_settings_are_refused_before_the_inputs_are_even_read(tmp_path, process_stub, options):
    with pytest.raises(ValueError):
        read_tools.run_fastp(tmp_path / "missing_R1.fastq", tmp_path / "missing_R2.fastq",
                             tmp_path / "trimmed", **options)


def test_a_thread_request_above_the_tools_own_limit_is_recorded_not_hidden(tmp_path, process_stub):
    _paths, result = trim(tmp_path, process_stub, threads=64)
    parameters = result["provenance"]["parameters"]
    assert parameters["threads_requested"] == 64
    assert parameters["threads_applied"] == read_tools.MAX_THREADS
    command = result["provenance"]["command"]
    assert command[command.index("--thread") + 1] == str(read_tools.MAX_THREADS)


def test_a_non_ascii_input_is_linked_rather_than_copied_or_renamed(tmp_path, process_stub):
    paths = reads(tmp_path)
    before = [(path.stat().st_ino, path.read_bytes()) for path in paths]
    result = read_tools.run_fastp(*paths, tmp_path / "trimmed")
    assert [(path.stat().st_ino, path.read_bytes()) for path in paths] == before
    assert all(path.is_file() for path in paths)
    access = result["provenance"]["tool_input_access"]
    assert all(description != "original path" for description in access)
    command = result["provenance"]["command"]
    assert "λ" not in command[command.index("--in1") + 1]


def test_an_absent_runtime_degrades_honestly_and_names_no_substitute(tmp_path, monkeypatch):
    capability = read_tools.runtime_capabilities(tmp_path / "nothing-staged-here")
    assert capability["available"] is False
    assert "No download happens during analysis" in capability["reason"]
    assert "unavailable" in capability["reason"]
    monkeypatch.setattr(read_tools, "runtime_capabilities", lambda root=None: capability)
    left, right = reads(tmp_path)
    with pytest.raises(read_tools.ReadTrimError, match="not staged"):
        read_tools.run_fastp(left, right, tmp_path / "trimmed")


def test_a_staged_package_with_the_wrong_version_is_not_used(tmp_path):
    root = tmp_path / "tools"
    root.mkdir()
    binary = root / ("fastp.exe" if os.name == "nt" else "fastp")
    binary.write_bytes(b"not fastp")
    (root / "manifest.json").write_text(json.dumps({
        "tool": "fastp", "version": "0.0.1", "platform": "linux-x64",
        "source_commit": read_tools.SOURCE_COMMIT,
        "files": [{"path": binary.name, "sha256": "x"}]}))
    capability = read_tools.runtime_capabilities(root)
    assert capability["available"] is False
    assert "does not match this build" in capability["reason"]


@pytest.mark.parametrize("mutation,error", [
    ("version", "staged engine"),
    ("impossible", "more reads or bases after filtering"),
    ("truncated", "complete fastp JSON report"),
])
def test_a_report_that_cannot_be_trusted_is_refused(tmp_path, mutation, error):
    if mutation == "version":
        fastp_json(tmp_path, version="0.23.4")
    elif mutation == "impossible":
        fastp_json(tmp_path, before_reads=100, after_reads=200)
    else:
        (tmp_path / "fastp.json").write_text(json.dumps({"summary": {}}))
    with pytest.raises(read_tools.ReadTrimError, match=error):
        read_tools.parse_report(tmp_path / "fastp.json")


def test_every_reported_rate_carries_the_denominator_it_was_measured_over(tmp_path):
    fastp_json(tmp_path)
    report = read_tools.parse_report(tmp_path / "fastp.json")
    assert "not of reads" in report["denominators"]["q20_rate/q30_rate"]
    assert "before filtering" in report["denominators"]["reads_removed_percent"]
    assert "estimate" in report["denominators"]["duplication_rate"]
    assert report["insert_size"]["undetermined_pairs"] == 12


def test_recording_a_trimming_run_keeps_the_originals_as_the_records_inputs(tmp_path, process_stub):
    paths = reads(tmp_path)
    result = read_tools.run_fastp(*paths, tmp_path / "trimmed")
    with Project(tmp_path / "study.wmlstudio") as project:
        first, second = [project.add_sample(path) for path in paths]
        project.set_result(first, {"kind": "fastq", "qc": {"records": 4}})
        assert read_tools.record_trimming(project, first, second, result) == [first, second]
        samples = [project.get_sample(identifier) for identifier in (first, second)]
        assert [sample["input_path"] for sample in samples] == [str(path) for path in paths]
        assert project.get_sample(first)["result"]["kind"] == "fastq"
        record = samples[0]["metadata"]["read_trimming"]
        assert record["status"] == "completed" and record["mate"] == 1
        assert record["paired_with"] == second
        assert samples[1]["metadata"]["read_trimming"]["paired_with"] == first
        assert record["report"]["after"]["total_reads"] == 900
        assert any(entry["action"] == "read_trimming_recorded" for entry in project.history(first))
    assert [path.is_file() for path in paths] == [True, True]


def test_recording_refuses_a_run_that_belongs_to_other_read_records(tmp_path, process_stub):
    paths = reads(tmp_path)
    result = read_tools.run_fastp(*paths, tmp_path / "trimmed")
    other = reads(tmp_path / "other")
    with Project(tmp_path / "study.wmlstudio") as project:
        first, second = [project.add_sample(path) for path in other]
        before = project.samples()
        with pytest.raises(ValueError, match="does not match"):
            read_tools.record_trimming(project, first, second, result)
        assert project.samples() == before


def test_the_assembler_is_offered_the_trimmed_pair_and_told_which_one_it_is(tmp_path, process_stub):
    paths = reads(tmp_path)
    result = read_tools.run_fastp(*paths, tmp_path / "trimmed")
    with Project(tmp_path / "study.wmlstudio") as project:
        first, second = [project.add_sample(path) for path in paths]
        samples = [project.get_sample(identifier) for identifier in (first, second)]
        assert read_tools.reads_for_assembly(*samples) == (paths[0], paths[1], None)
        read_tools.record_trimming(project, first, second, result)
        samples = [project.get_sample(identifier) for identifier in (first, second)]
        read1, read2, source = read_tools.reads_for_assembly(*samples)
        assert [str(read1), str(read2)] == [result["read1_path"], result["read2_path"]]
        assert source["kind"] == "fastp_trimmed"
        assert source["read1_sha256"] == result["provenance"]["outputs"][0]["sha256"]
        assert [item["path"] for item in source["originals"]] == [str(path) for path in paths]
        assert read_tools.reads_for_assembly(*samples, prefer_trimmed=False) == (paths[0], paths[1], None)
        Path(result["read1_path"]).unlink()
        assert read_tools.reads_for_assembly(*samples) == (paths[0], paths[1], None)


def test_the_read_source_a_run_hands_the_assembler_matches_its_own_outputs(tmp_path, process_stub):
    _paths, result = trim(tmp_path, process_stub)
    source = read_tools.read_source_for(result)
    assert source["kind"] == "fastp_trimmed"
    assert source["read2_sha256"] == result["provenance"]["outputs"][1]["sha256"]
    assert source["reads_removed"] == 100
    assert "unchanged" in source["description"]


def test_the_read_tab_says_plainly_when_a_record_was_not_trimmed(tmp_path, process_stub):
    view = read_tools.trimming_view({"id": "s1", "name": "isolate", "input_path": "R1.fastq"})
    assert view["status"] == "not_trimmed" and view["report"] is None
    assert any("were not trimmed" in note for note in view["notes"])
    assert read_tools.TRIMMING_DISCLAIMER in view["notes"]
    assembled = read_tools.trimming_view({"id": "s2", "name": "isolate", "input_path": "contigs.fasta",
                                          "metadata": {"workflow": {"source_kind": "assembly"}}})
    assert assembled["status"] == "not_reads"


def test_the_read_tab_shows_the_recorded_run_without_inventing_numbers(tmp_path, process_stub):
    paths = reads(tmp_path)
    result = read_tools.run_fastp(*paths, tmp_path / "trimmed")
    with Project(tmp_path / "study.wmlstudio") as project:
        first, second = [project.add_sample(path) for path in paths]
        read_tools.record_trimming(project, first, second, result)
        view = read_tools.trimming_view(project.get_sample(first))
    assert view["status"] == "trimmed" and view["mate"] == 1
    assert view["report"] == result["report"]
    assert view["original_path"] == str(paths[0])
    assert view["trimmed_path"] == result["read1_path"]
    assert any("original FASTQ is unchanged" in note for note in view["notes"])


def fake_upstream(tmp_path, monkeypatch, *, binary=b"TEST ONLY fastp fixture"):
    cache = tmp_path / "cache"
    cache.mkdir()
    (cache / f"fastp.{staging.VERSION}").write_bytes(binary)
    archive = cache / staging.SOURCE_ARCHIVE
    with tarfile.open(archive, "w:gz") as bundle:
        for name, content in (("LICENSE", b"MIT License test fixture"),
                              ("README.md", b"# fastp test fixture"),
                              ("src/fastp.cpp", b"int main() { return 0; }")):
            info = tarfile.TarInfo(f"fastp-{staging.SOURCE_COMMIT}/{name}")
            info.size = len(content)
            bundle.addfile(info, io.BytesIO(content))
    monkeypatch.setattr(staging, "LINUX_SHA256", hashlib.sha256(binary).hexdigest())
    monkeypatch.setattr(staging, "SOURCE_SHA256", staging.digest(archive))
    return cache


def test_staging_verifies_every_file_and_reuses_a_valid_snapshot(tmp_path, monkeypatch):
    cache = fake_upstream(tmp_path, monkeypatch)
    destination = tmp_path / "fastp"
    manifest = staging.stage(destination, "linux-x64", cache=cache)
    assert manifest["license"].startswith("MIT")
    assert (destination / "LICENSE").read_bytes() == b"MIT License test fixture"
    assert (destination / staging.SOURCE_ARCHIVE).is_file()
    assert not (destination / "src").exists()
    names = {entry["path"] for entry in manifest["files"]}
    assert {"fastp", "LICENSE", staging.SOURCE_ARCHIVE, "fastp-README.md"} <= names
    assert staging.stage(destination, "linux-x64", cache=cache)["files"] == manifest["files"]
    (destination / "LICENSE").write_text("tampered")
    with pytest.raises(ValueError, match="checksum mismatch"):
        staging.verify(destination, "linux-x64")


def test_staging_refuses_a_download_that_does_not_match_its_pin(tmp_path, monkeypatch):
    cache = fake_upstream(tmp_path, monkeypatch)
    monkeypatch.setattr(staging, "LINUX_SHA256", "0" * 64)
    destination = tmp_path / "fastp"
    with pytest.raises(ValueError, match="checksum mismatch"):
        staging.stage(destination, "linux-x64", cache=cache)
    assert not destination.exists()
    assert not list(tmp_path.glob(".fastp-stage-*"))


def test_windows_staging_refuses_to_invent_an_upstream_binary(tmp_path, monkeypatch):
    cache = fake_upstream(tmp_path, monkeypatch)
    with pytest.raises(ValueError, match="no official Windows binary"):
        staging.stage(tmp_path / "windows", "windows-x64", cache=cache)
    assert not (tmp_path / "windows").exists()
    with pytest.raises(ValueError, match="Unsupported fastp platform"):
        staging.stage(tmp_path / "other", "darwin-arm64", cache=cache)


def test_a_windows_artifact_without_smoke_evidence_is_not_staged(tmp_path, monkeypatch):
    source = tmp_path / "artifact"
    source.mkdir()
    files = {"fastp.exe": b"MZ test fixture", "LICENSE": b"MIT License test fixture",
             staging.SOURCE_ARCHIVE: b"source archive fixture",
             "smoke-test.json": json.dumps({"passed": False}).encode()}
    for name, content in files.items():
        (source / name).write_bytes(content)
    (source / "manifest.json").write_text(json.dumps({
        "tool": "fastp", "version": staging.VERSION, "platform": "windows-x64",
        "source_commit": staging.SOURCE_COMMIT,
        "files": [{"path": name, "sha256": hashlib.sha256(content).hexdigest()}
                  for name, content in files.items()]}))
    with pytest.raises(ValueError, match="smoke evidence"):
        staging.stage(tmp_path / "windows", "windows-x64", source=source)
    assert not (tmp_path / "windows").exists()


@pytest.mark.skipif(not os.environ.get("WMLSTUDIO_TEST_FASTP_ROOT"),
                    reason="Explicit genuine fastp integration opt-in")
def test_genuine_fastp_trims_a_unicode_named_pair_without_changing_it(tmp_path):
    paths = reads(tmp_path, count=200, length=100)
    before = [path.read_bytes() for path in paths]
    result = read_tools.run_fastp(*paths, tmp_path / "trimmed é",
                                  root=os.environ["WMLSTUDIO_TEST_FASTP_ROOT"], threads=2)
    assert result["report"]["version"] == read_tools.VERSION
    assert result["report"]["before"]["total_reads"] == 400
    assert [path.read_bytes() for path in paths] == before
    assert all(Path(item["trimmed_path"]).is_file() for item in result["provenance"]["link"])
    source = read_tools.read_source_for(result)
    assert source["read1_sha256"] != source["originals"][0]["sha256"]
