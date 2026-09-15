"""Adapter tests use explicit process stubs; real assembly is gated in Windows CI."""

import hashlib
import json
import os
import sys
import threading
from pathlib import Path

import pytest

from wmlstudio import assembly
from wmlstudio.project import Project
from wmlstudio.sequence import AnalysisCancelled, SequenceError
from wmlstudio.storage import organize_sample


def reads(tmp_path, count=3, marked=True):
    paths = [tmp_path / "reads, λ_R1.fastq", tmp_path / "reads, λ_R2.fastq"]
    for mate, path in enumerate(paths, 1):
        path.write_text("".join(
            f"@read{i}{'/' + str(mate) if marked else ''}\nACGTTGCAACGT\n+\nIIIIIIIIIIII\n"
            for i in range(count)), encoding="ascii")
    return paths


@pytest.fixture
def process_stub(tmp_path, monkeypatch):
    """Only adapter contract tests; this fixture does not validate an assembler."""
    executable = tmp_path / "test-only-tool"
    executable.write_bytes(b"adapter fixture, not an assembler")
    monkeypatch.setattr(assembly, "_tool_info", lambda exe, directory, cancel: {
        "name": "TEST-ONLY adapter stub", "executable_sha256": "stub"})
    calls = []

    def run(command, directory, log_path, cancelled=None, progress=None, timeout=None):
        calls.append(command)
        plan = json.loads((directory / "run-plan.json").read_text())
        assert plan["state"] == "prepared_not_completed"
        assert plan["inputs"][0]["sha256"]
        assert command[command.index("--reads") + 1] == "mate1.fastq,mate2.fastq"
        assert (directory / "mate1.fastq").read_text().startswith("@pair1/1\n")
        assert (directory / "mate2.fastq").read_text().startswith("@pair1/2\n")
        (directory / "contigs.fasta").write_text(">test_fixture_contig\nACGTTGCAACGT\n")
        Path(log_path).write_text("TEST ONLY subprocess fixture\n")
        return 0

    monkeypatch.setattr(assembly, "_run_process", run)
    return executable, calls


def test_adapter_publishes_validated_evidence_and_never_modifies_originals(tmp_path, process_stub):
    paths = reads(tmp_path)
    original = [p.read_bytes() for p in paths]
    output = tmp_path / "assembly"
    result = assembly.run_skesa(*paths, output, executable=process_stub[0], min_contig=1)
    assert [p.read_bytes() for p in paths] == original
    assert Path(result["assembly_path"]).is_file()
    assert result["pairing"]["records_checked"] == 3
    assert result["pairing"]["complete_file"] and not result["pairing"]["sampled"]
    assert result["read_qc"][0]["records"] == 3
    assert result["qc"]["n50"] == 12
    assert result["provenance"]["inputs"][0]["sha256"] == hashlib.sha256(original[0]).hexdigest()
    assert result["provenance"]["biological_qc"]["contamination"] == "not assessed"
    assert json.loads((output / "assembly.json").read_text()) == result
    assert not list(output.glob("*.fastq"))
    assert not list(tmp_path.glob(".wmlstudio-assembly-*"))


def test_unmarked_explicitly_chosen_pairs_are_not_claimed_biologically_verified(tmp_path, process_stub):
    result = assembly.run_skesa(*reads(tmp_path, marked=False), tmp_path / "out", executable=process_stub[0])
    assert result["pairing"]["status"] == "unmarked"
    assert any("user-assigned" in note for note in result["notes"])


@pytest.mark.parametrize("mutation,error", [
    ("mismatch_last", "identifiers differ"),
    ("truncated_last", "truncated"),
    ("extra_record", "different numbers"),
    ("reversed_mate", "reversed"),
    ("fasta", "two FASTQ"),
])
def test_full_pair_validation_precedes_process_launch(tmp_path, process_stub, mutation, error):
    left, right = reads(tmp_path)
    text = right.read_text()
    if mutation == "mismatch_last":
        text = text.replace("@read2/2", "@different/2")
    elif mutation == "truncated_last":
        text = text.rsplit("\n", 3)[0] + "\n+\nI\n"
    elif mutation == "extra_record":
        text += "@read3/2\nACGT\n+\nIIII\n"
    elif mutation == "reversed_mate":
        text = text.replace("/2", "/1")
    else:
        text = ">contig\nACGT\n"
    right.write_text(text)
    with pytest.raises(SequenceError, match=error):
        assembly.run_skesa(left, right, tmp_path / "out", executable=process_stub[0])
    assert not process_stub[1]
    assert not (tmp_path / "out").exists()


def test_existing_output_and_hardlinked_mates_cannot_be_overwritten(tmp_path, process_stub):
    left, right = reads(tmp_path)
    output = tmp_path / "existing"
    output.mkdir()
    with pytest.raises(FileExistsError):
        assembly.run_skesa(left, right, output, executable=process_stub[0])
    alias = tmp_path / "same-data.fastq"
    os.link(left, alias)
    with pytest.raises(SequenceError, match="different files"):
        assembly.run_skesa(left, alias, tmp_path / "new", executable=process_stub[0])


def test_failure_preserves_log_without_publishing_partial_assembly(tmp_path, process_stub, monkeypatch):
    def failure(command, directory, log_path, *args, **kwargs):
        Path(log_path).write_text("real diagnostic example: allocation failed")
        (directory / "contigs.fasta").write_text(">partial\nACGT\n")
        return 7
    monkeypatch.setattr(assembly, "_run_process", failure)
    with pytest.raises(assembly.AssemblyError, match="code 7"):
        assembly.run_skesa(*reads(tmp_path), tmp_path / "out", executable=process_stub[0])
    assert not (tmp_path / "out").exists()
    assert "allocation failed" in next(tmp_path.glob("out.failed-*.log")).read_text()
    plan = json.loads(next(tmp_path.glob("out.failed-*.json")).read_text())
    assert plan["state"] == "failed" and plan["inputs"][0]["sha256"]
    assert not list(tmp_path.glob(".wmlstudio-assembly-*"))


def test_success_exit_without_contigs_is_not_success(tmp_path, process_stub, monkeypatch):
    def empty(command, directory, log_path, *args, **kwargs):
        Path(log_path).write_text("no adequate coverage")
        return 0
    monkeypatch.setattr(assembly, "_run_process", empty)
    with pytest.raises(assembly.AssemblyError, match="no contigs"):
        assembly.run_skesa(*reads(tmp_path), tmp_path / "out", executable=process_stub[0])
    assert not (tmp_path / "out").exists()


def test_real_subprocess_cancellation_terminates_child_and_retains_log(tmp_path):
    event = threading.Event()
    timer = threading.Timer(0.5, event.set)
    timer.start()
    try:
        with pytest.raises(AnalysisCancelled):
            assembly._run_process(
                [sys.executable, "-u", "-c", "import time; print('started', flush=True); time.sleep(60)"],
                tmp_path, tmp_path / "child.log", event.is_set)
    finally:
        timer.cancel()
    assert "started" in (tmp_path / "child.log").read_text()


def test_tool_manifest_detects_tampered_executable(tmp_path, monkeypatch):
    executable = tmp_path / "skesa"
    executable.write_bytes(b"tampered")
    (tmp_path / "manifest.json").write_text(json.dumps({
        "tool": "SKESA", "files": [{"path": "skesa", "sha256": "wrong"}]}))
    def version(command, directory, log_path, *args, **kwargs):
        Path(log_path).write_text("SKESA 2.4.0\n")
        return 0
    monkeypatch.setattr(assembly, "_run_process", version)
    with pytest.raises(assembly.AssemblyError, match="integrity"):
        assembly._tool_info(executable, tmp_path, None)


def test_association_preserves_both_read_records_history_and_fasta_extension(tmp_path, process_stub):
    paths = reads(tmp_path)
    originals = [p.read_bytes() for p in paths]
    result = assembly.run_skesa(*paths, tmp_path / "out", executable=process_stub[0])
    with Project(tmp_path / "study.wmlstudio") as project:
        first, second = [project.add_sample(path) for path in paths]
        project.set_result(first, {"kind": "fastq", "qc": {"records": 3}})
        project.set_metadata(first, {"organism": {"genus": "Escherichia", "species": "coli"},
                                    "workflow": {"source_path": str(paths[0])}})
        assert assembly.associate_assembly(project, first, second, result) == first
        assert len(project.samples()) == 2
        primary, mate = project.get_sample(first), project.get_sample(second)
        assert primary["input_path"] == result["assembly_path"] and primary["status"] == "queued"
        assert primary["result"] is None
        assert mate["input_path"] == str(paths[1])
        assert mate["metadata"]["workflow"] == {
            "source_kind": "read_mate", "paired_with": first, "assembly_path": result["assembly_path"]}
        assert any(h["action"] == "result_invalidated" for h in project.history(first))
        project.set_result(first, {"kind": "fasta", "status": "complete", "st": "131",
                                  "input_sha256": result["provenance"]["assembly_sha256"]})
        path = organize_sample(project, first, tmp_path / "managed", append_st=True)
        assert path.suffix == ".fasta" and "_ST_131" in path.name
        assert project.get_sample(first)["metadata"]["assembly"]["assembly_path"] == str(path)
        assert [p.read_bytes() for p in paths] == originals


def test_association_mapping_failure_is_atomic(tmp_path, process_stub):
    paths = reads(tmp_path)
    result = assembly.run_skesa(*paths, tmp_path / "out", executable=process_stub[0])
    with Project(tmp_path / "study.wmlstudio") as project:
        first, second = [project.add_sample(path) for path in paths]
        before = project.samples()
        with pytest.raises(ValueError, match="does not match"):
            assembly.associate_assembly(project, second, first, result)
        assert project.samples() == before


def test_cancel_before_work_cannot_publish(tmp_path, process_stub):
    with pytest.raises(AnalysisCancelled):
        assembly.run_skesa(*reads(tmp_path), tmp_path / "out", executable=process_stub[0],
                           cancelled=lambda: True)
    assert not (tmp_path / "out").exists()


def test_upstream_memory_reserve_is_explained_before_launch(tmp_path, process_stub):
    with pytest.raises(ValueError, match="reserves 2 GB"):
        assembly.run_skesa(*reads(tmp_path), tmp_path / "out", executable=process_stub[0], memory_gb=2)


def trimmed_source(paths, **overrides):
    """A trimming record that honestly describes these exact files."""
    digests = [hashlib.sha256(path.read_bytes()).hexdigest() for path in paths]
    return {"kind": "fastp_trimmed", "tool": "fastp", "version": "1.3.7",
            "description": "fastp 1.3.7 trimmed and quality-filtered reads.",
            "read1_path": str(paths[0]), "read2_path": str(paths[1]),
            "read1_sha256": digests[0], "read2_sha256": digests[1],
            "originals": [{"mate": 1, "path": "raw_R1.fastq", "sha256": "a" * 64},
                          {"mate": 2, "path": "raw_R2.fastq", "sha256": "b" * 64}],
            "reads_removed": 1200, "bases_removed": 60000, **overrides}


def test_an_untrimmed_run_says_so_and_claims_no_preprocessing(tmp_path, process_stub):
    result = assembly.run_skesa(*reads(tmp_path), tmp_path / "out", executable=process_stub[0], min_contig=1)
    assert result["read_source"] == assembly.RAW_READ_SOURCE
    assert result["provenance"]["read_preprocessing"].startswith("Format normalization only")
    assert any("not a species" in note for note in result["notes"])


def test_an_assembly_records_the_trimmed_reads_it_actually_used(tmp_path, process_stub):
    paths = reads(tmp_path)
    result = assembly.run_skesa(*paths, tmp_path / "out", executable=process_stub[0], min_contig=1,
                                read_source=trimmed_source(paths))
    assert result["read_source"]["kind"] == "fastp_trimmed"
    assert result["read_source"]["reads_removed"] == 1200
    assert "fastp 1.3.7" in result["provenance"]["read_preprocessing"]
    assert "1,200 reads were removed" in result["provenance"]["read_preprocessing"]
    assert result["provenance"]["read_source"] == result["read_source"]
    assert json.loads((tmp_path / "out/assembly.json").read_text())["read_source"] == result["read_source"]


def test_an_assembly_cannot_borrow_a_trimming_record_from_other_files(tmp_path, process_stub):
    paths = reads(tmp_path)
    borrowed = trimmed_source(paths, read1_sha256="c" * 64)
    with pytest.raises(assembly.AssemblyError, match="different files"):
        assembly.run_skesa(*paths, tmp_path / "out", executable=process_stub[0], read_source=borrowed)
    assert not (tmp_path / "out").exists()
    with pytest.raises(ValueError, match="either fastp-trimmed reads or nothing"):
        assembly.run_skesa(*paths, tmp_path / "other", executable=process_stub[0],
                           read_source={"kind": "bbduk_trimmed"})


def test_metrics_and_depth_come_from_the_assembly_the_run_produced(tmp_path, process_stub):
    result = assembly.run_skesa(*reads(tmp_path), tmp_path / "out", executable=process_stub[0], min_contig=1)
    metrics = result["metrics"]
    assert metrics["contigs"] == 1 and metrics["total_length"] == 12
    assert metrics["n50"] == 12 and metrics["largest_contig"] == 12
    assert metrics["assembler_depth"]["reported"] is False
    assert "does not carry SKESA contig names" in metrics["assembler_depth"]["basis"]
    coverage = result["coverage"]
    assert coverage["input_read_bases"] == 72
    assert coverage["read_bases_per_assembled_base"] == pytest.approx(6.0)
    assert "Not genome coverage" in coverage["basis"]


def test_skesa_contig_names_are_read_as_the_assemblers_own_depth_estimate(tmp_path):
    path = tmp_path / "contigs.fasta"
    path.write_text(">Contig_1_35.5\n" + "ACGT" * 25 + "\n"
                    ">Contig_2_12.25_Circ\n" + "ACGT" * 5 + "\n"
                    ">Contig_3_8\n" + "ACGTN" * 4 + "\n")
    metrics = assembly.assembly_metrics(path)
    assert metrics["contigs"] == 3 and metrics["total_length"] == 140
    assert metrics["largest_contig"] == 100 and metrics["smallest_contig"] == 20
    assert metrics["n50"] == 100
    depth = metrics["assembler_depth"]
    assert depth["reported"] is True and depth["contigs_with_depth"] == 3
    assert depth["minimum"] == 8 and depth["maximum"] == 35.5
    assert depth["circular_contigs"] == 1
    assert depth["length_weighted_mean"] == pytest.approx((35.5 * 100 + 12.25 * 20 + 8 * 20) / 140)
    assert metrics["n_bases"] == 4 and metrics["ambiguous_bases"] == 4
    assert "not a species" in metrics["interpretation"]


def test_metrics_refuse_reads_and_report_no_depth_for_a_plain_assembly(tmp_path):
    left, _right = reads(tmp_path)
    with pytest.raises(assembly.AssemblyError, match="not raw reads"):
        assembly.assembly_metrics(left)
    plain = tmp_path / "plain.fasta"
    plain.write_text(">contig1\nACGTACGTAC\n")
    assert assembly.assembly_metrics(plain)["assembler_depth"]["reported"] is False


def test_the_assembly_tab_shows_one_samples_recorded_evidence(tmp_path, process_stub):
    paths = reads(tmp_path)
    result = assembly.run_skesa(*paths, tmp_path / "out", executable=process_stub[0], min_contig=1,
                                read_source=trimmed_source(paths))
    with Project(tmp_path / "study.wmlstudio") as project:
        first, second = [project.add_sample(path) for path in paths]
        assembly.associate_assembly(project, first, second, result)
        primary = assembly.assembly_view(project.get_sample(first))
        mate = assembly.assembly_view(project.get_sample(second))
    assert primary["status"] == "assembled"
    assert primary["from_trimmed_reads"] is True
    assert primary["read_source"]["kind"] == "fastp_trimmed"
    assert primary["metrics"]["n50"] == 12
    assert primary["coverage"]["read_bases_per_assembled_base"] == pytest.approx(6.0)
    assert primary["assembly_path"] == result["assembly_path"]
    assert any("not a species" in note for note in primary["notes"])
    assert mate["status"] == "read_mate"
    assert any("second mate" in note for note in mate["notes"])


def test_the_assembly_tab_never_pretends_an_unassembled_record_is_assembled(tmp_path, process_stub):
    paths = reads(tmp_path)
    contigs = tmp_path / "supplied.fasta"
    contigs.write_text(">contig1\nACGTACGTAC\n")
    with Project(tmp_path / "study.wmlstudio") as project:
        read_id = project.add_sample(paths[0])
        supplied_id = project.add_sample(contigs)
        unread_id = project.add_sample(paths[1])
        paths[1].unlink()
        views = {name: assembly.assembly_view(project.get_sample(identifier))
                 for name, identifier in (("reads", read_id), ("supplied", supplied_id),
                                          ("missing", unread_id))}
    assert views["reads"]["status"] == "reads_not_assembled"
    assert views["reads"]["metrics"] is None and views["reads"]["from_trimmed_reads"] is None
    assert views["supplied"]["status"] == "assembly_supplied"
    assert any("supplied, not produced here" in note for note in views["supplied"]["notes"])
    assert views["missing"]["status"] == "unknown"
    assert any("could not be read" in note for note in views["missing"]["notes"])


def test_an_assembly_recorded_before_metrics_existed_admits_what_is_missing():
    view = assembly.assembly_view({"id": "s1", "name": "isolate", "metadata": {"assembly": {
        "assembly_path": "/tmp/contigs.fasta",
        "provenance": {"engine": {"name": "SKESA", "version": "2.4.0"},
                       "read_preprocessing": "Format normalization only."}}}})
    assert view["status"] == "assembled" and view["metrics"] is None
    assert view["from_trimmed_reads"] is None
    assert view["read_source"]["kind"] == "unrecorded"
    assert any("before per-contig metrics were kept" in note for note in view["notes"])
    assert any("No depth was recorded" in note for note in view["notes"])
