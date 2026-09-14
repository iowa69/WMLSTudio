import json
import os
import zipfile
from pathlib import Path

import pytest

from wmlstudio import fastqc
from wmlstudio.fastqc_dialog import reads_for_sample, run_project_fastqc
from wmlstudio.scheduler import GIB, HardwareSnapshot, plan_resources
from wmlstudio.sequence import AnalysisCancelled, file_sha256


def report(path, version="0.12.1", status="WARN"):
    with zipfile.ZipFile(path, "w") as archive:
        archive.writestr("reads_fastqc/summary.txt", f"PASS\tBasic Statistics\treads é.fastq\n{status}\tAdapter Content\treads é.fastq\n")
        archive.writestr("reads_fastqc/fastqc_data.txt", f"##FastQC\t{version}\n>>Basic Statistics\tpass\n#Measure\tValue\nTotal Sequences\t2\n>>END_MODULE\n")


def reads(path):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("@r1\nACGT\n+\nIIII\n@r2\nTGCA\n+\n!!!!\n", encoding="ascii")
    return path


@pytest.fixture
def fake_runtime(tmp_path, monkeypatch):
    root = tmp_path / "tools"
    root.mkdir()
    (root / "manifest.json").write_text("{}")
    monkeypatch.setattr(fastqc, "runtime_capabilities", lambda value=None: {
        "available": True, "java": str(root / "java.exe"), "fastqc_root": str(root),
        "version": "0.12.1", "java_version": "17"})
    commands = []
    def run(command, directory, log_path, cancelled, progress):
        commands.append(command)
        output = Path(next(value.split("=", 1)[1] for value in command if value.startswith("-Dfastqc.output_dir=")))
        report(output / "reads_fastqc.zip")
        (output / "reads_fastqc.html").write_text("<h1>FastQC</h1>")
        log_path.write_text("Complete")
        return 0
    monkeypatch.setattr(fastqc, "_run", run)
    return commands


def test_parse_original_module_flags_and_version(tmp_path):
    path = tmp_path / "report.zip"
    report(path)
    parsed = fastqc.parse_report(path)
    assert parsed["qc_status"] == "WARN"
    assert parsed["metrics"]["Total Sequences"] == "2"
    assert parsed["modules"][0]["filename"] == "reads é.fastq"
    report(path, version="0.11.9")
    with pytest.raises(fastqc.FastQCError, match="version"):
        fastqc.parse_report(path)
    report(path, status="MAYBE")
    with pytest.raises(fastqc.FastQCError, match="status"):
        fastqc.parse_report(path)


def test_run_preserves_inputs_and_separates_same_named_files(tmp_path, fake_runtime):
    inputs = [reads(tmp_path / name / "reads é.fastq") for name in ("one", "two")]
    before = [file_sha256(path) for path in inputs]
    result = fastqc.run_fastqc(inputs, tmp_path / "reports", threads=2, memory_gb=1)
    assert result["status"] == "completed"
    assert [item["qc_status"] for item in result["reports"]] == ["WARN", "WARN"]
    assert len({item["report_zip"] for item in result["reports"]}) == 2
    assert all(Path(item["report_html"]).is_file() for item in result["reports"])
    assert [file_sha256(path) for path in inputs] == before
    assert result["provenance"]["files_parallel"] == 1
    for command in fake_runtime:
        assert "-Dfile.encoding=UTF-8" in command
        assert "-XX:+UseSerialGC" in command
        assert "-XX:ActiveProcessorCount=2" in command
        assert "-Dfastqc.threads=1" in command
    assert json.loads((tmp_path / "reports/fastqc-result.json").read_text())["inputs"] == result["inputs"]
    with pytest.raises(ValueError, match="already exists"):
        fastqc.run_fastqc(inputs, tmp_path / "reports")


def test_failed_or_cancelled_fastqc_publishes_no_partial_report(tmp_path, fake_runtime, monkeypatch):
    source = reads(tmp_path / "input.fastq")
    def fail(*args):
        args[2].write_text("Real engine error")
        return 1
    monkeypatch.setattr(fastqc, "_run", fail)
    with pytest.raises(fastqc.FastQCError, match="Real engine error"):
        fastqc.run_fastqc([source], tmp_path / "failure")
    assert not (tmp_path / "failure").exists()
    with pytest.raises(AnalysisCancelled):
        fastqc.run_fastqc([source], tmp_path / "cancelled", cancelled=lambda: True)
    assert not (tmp_path / "cancelled").exists()
    assert not list(tmp_path.glob(".fastqc-run-*"))


def test_runtime_unavailable_and_invalid_input_rejected(tmp_path, fake_runtime, monkeypatch):
    source = reads(tmp_path / "input.fastq")
    with pytest.raises(ValueError, match="more than once"):
        fastqc.run_fastqc([source, source], tmp_path / "duplicate")
    for options in ({"threads": True}, {"threads": 0}, {"memory_gb": 0}):
        with pytest.raises(ValueError):
            fastqc.run_fastqc([source], tmp_path / "invalid", **options)
    monkeypatch.setattr(fastqc, "runtime_capabilities", lambda value=None: {"available": False, "message": "Unavailable"})
    with pytest.raises(fastqc.FastQCError, match="Unavailable"):
        fastqc.run_fastqc([source], tmp_path / "missing")


def test_project_links_reports_without_replacing_assembly_or_typing(tmp_path, fake_runtime):
    from wmlstudio.project import Project
    first, second = [reads(tmp_path / f"R{index}.fastq") for index in (1, 2)]
    assembly = tmp_path / "isolate.fasta"
    assembly.write_text(">c\nACGT\n")
    with Project(tmp_path / "project.wmlstudio") as project:
        identifier = project.add_sample(assembly)
        project.set_result(identifier, {"kind": "fasta", "status": "complete", "st": "20"})
        project.update_metadata(identifier, {"workflow": {"read1_path": str(first), "read2_path": str(second)}})
        before = project.get_sample(identifier)
        allocation = plan_resources(threads_per_sample=1, memory_gb=1, max_parallel=1,
                                    hardware=HardwareSnapshot(2, 8 * GIB))
        result = run_project_fastqc(project, [before], tmp_path / "reports", allocation)
        after = project.get_sample(identifier)
        assert result[0]["sample_id"] == identifier
        assert after["input_path"] == before["input_path"]
        assert after["result"] == before["result"]
        assert after["metadata"]["fastqc"]["engine"] == "FastQC"
        assert len(after["metadata"]["fastqc"]["reports"]) == 2


def test_attached_original_reads_override_primary_assembly():
    sample = {"input_path": "assembly.fasta", "metadata": {"reads": {
        "reads": [{"mate": 2, "path": "R2.fastq"}, {"mate": 1, "path": "R1.fastq"}]}}}
    assert reads_for_sample(sample) == ["R1.fastq", "R2.fastq"]
    assert reads_for_sample({"input_path": "assembly.fasta"}) == []


@pytest.mark.skipif(not os.environ.get("WMLSTUDIO_TEST_FASTQC_ROOT"), reason="Explicit genuine FastQC integration opt-in")
def test_genuine_fastqc_on_unicode_input(tmp_path):
    source = reads(tmp_path / "original reads é.fastq")
    before = file_sha256(source)
    result = fastqc.run_fastqc([source], tmp_path / "native reports é", root=os.environ["WMLSTUDIO_TEST_FASTQC_ROOT"])
    assert result["version"] == "0.12.1"
    assert result["reports"][0]["metrics"]["Total Sequences"] == "2"
    assert result["reports"][0]["metrics"]["Filename"] == source.name
    assert file_sha256(source) == before
