"""Runtime isolation, safe native invocation, immutable database publication."""

import json
import sys
import time
from pathlib import Path

import pytest

from wmlstudio import hydra_runtime as runtime
from wmlstudio.sequence import AnalysisCancelled


def database(root):
    target = root / "nucl/ncbi"
    target.mkdir(parents=True)
    (target / "sequences.fna").write_text(">ref\nACGTACGT\n")
    (target / "meta.tsv").write_text("seqid\tgene\nref\ttest\n")
    (root / "manifest.json").write_text(json.dumps({"hydra_db_version": 1,
        "databases": {"ncbi": {"path": "nucl/ncbi", "version": "test-snapshot", "kind": "nucl"}}}))
    return root


@pytest.fixture
def ready(monkeypatch, tmp_path):
    root = database(tmp_path / "reference")
    monkeypatch.setattr(runtime, "_ncbi_version", lambda: "2026-08-07.1")
    monkeypatch.setattr(runtime, "runtime_capabilities", lambda db_root=None: {
        "assembly_available": True, "tools": {name: f"/tools/{name}" for name in runtime.TOOLS},
        "hydra_version": runtime.HYDRA_VERSION,
        "databases": runtime.installed_databases(db_root) if db_root else {}, "limitations": ["Assembly evidence only."]})
    return root


def test_capabilities_and_missing_store_are_read_only(tmp_path):
    root = tmp_path / "absent"
    info = runtime.runtime_capabilities(root)
    assert "available" in info and "message" in info and "blastn" in info
    assert info["reads_available"] is False
    assert info["databases"] == {}
    assert not root.exists()


def test_manifest_rejects_external_paths_and_invalid_shapes(tmp_path):
    for entry in [{"path": "../outside"}, "not a dictionary", {"path": 1}]:
        (tmp_path / "manifest.json").write_text(json.dumps({"databases": {"ncbi": entry}}))
        with pytest.raises(runtime.HydraRuntimeError):
            runtime.installed_databases(tmp_path)


def test_assembly_runs_real_upstream_argument_contract_and_preserves_inputs(ready, tmp_path, monkeypatch):
    assembly = tmp_path / "isolate space é.fasta"
    assembly.write_text(">contig\nACGTACGT\n")
    original = assembly.read_bytes()
    captured = []

    def child(arguments, directory, cancelled=None, progress=None):
        captured.extend(arguments)
        output = Path(arguments[arguments.index("--outdir") + 1])
        output.mkdir()
        (output / "hydra.json").write_text(json.dumps({
            "hydra_version": "1.4.0", "command": "hydra run", "databases": ["ncbi"],
            "parameters": {}, "samples": [{"sample": "sample_a", "input_type": "assembly",
                                             "inputs": [str(assembly)], "hits": []}]}))
        return ["worker", *arguments], "Completed original HYDRA pipeline"

    monkeypatch.setattr(runtime, "_run_child", child)
    result = runtime.run_assemblies([assembly], ready, ["ncbi"], sample_names=["sample_a"],
                                    work_root=tmp_path / "jobs")
    assert captured[0] == "run"
    assert "--no-mlst" in captured and "--no-typing" in captured
    assert "--no-auto-organism" in captured and "--no-heteroresistance" in captured
    assert captured[captured.index("--assembly") + 1] == str(assembly)
    assert captured[captured.index("--name") + 1] == "sample_a"
    provenance = result["execution_provenance"]
    assert provenance["source_commit"] == runtime.HYDRA_COMMIT
    assert provenance["inputs"][0]["sha256"] == runtime.file_sha256(assembly)
    assert "ncbi" in provenance["reference_snapshot"]["databases"]
    assert assembly.read_bytes() == original
    assert list((tmp_path / "jobs").iterdir()) == []


def test_explicit_organism_and_options_are_forwarded_without_guessing(ready, tmp_path, monkeypatch):
    assembly = tmp_path / "isolate.fasta"
    assembly.write_text(">contig\nACGT\n")
    captured = []

    def child(arguments, directory, *args):
        captured.extend(arguments)
        raise runtime.HydraRuntimeError("stopped after capturing argv")

    monkeypatch.setattr(runtime, "_run_child", child)
    with pytest.raises(runtime.HydraRuntimeError, match="capturing"):
        runtime.run_assemblies([assembly], ready, ["ncbi"], organism="Staphylococcus_aureus",
                               protein=False, point_mutations=False)
    assert captured[captured.index("--organism") + 1] == "Staphylococcus_aureus"
    assert "--no-auto-organism" not in captured
    assert "--no-protein" in captured and "--no-point-mutations" in captured


def test_reads_and_duplicate_inputs_are_not_silently_analyzed(ready, tmp_path):
    reads = tmp_path / "reads.fasta"
    reads.write_text("@read\nACGT\n+\nIIII\n")
    with pytest.raises(runtime.HydraRuntimeError, match="FASTA assemblies only"):
        runtime.run_assemblies([reads], ready)
    with pytest.raises(runtime.HydraRuntimeError, match="distinct"):
        runtime.run_assemblies([reads, reads], ready)
    with pytest.raises(runtime.HydraRuntimeError, match="Threads"):
        runtime.run_assemblies([reads], ready, threads=0)


@pytest.mark.parametrize("value", [-1, 101, float("nan"), float("inf"), True, "90"])
def test_invalid_scientific_thresholds_are_rejected(ready, tmp_path, value):
    assembly = tmp_path / "isolate.fa"
    assembly.write_text(">contig\nACGT\n")
    with pytest.raises(runtime.HydraRuntimeError, match="percentage"):
        runtime.run_assemblies([assembly], ready, min_identity=value)


def test_missing_database_does_not_trigger_network_or_child(ready, tmp_path, monkeypatch):
    assembly = tmp_path / "isolate.fa"
    assembly.write_text(">contig\nACGT\n")
    monkeypatch.setattr(runtime, "_run_child", lambda *args: pytest.fail("must not launch"))
    with pytest.raises(runtime.HydraRuntimeError, match="not installed"):
        runtime.run_assemblies([assembly], ready, ["card"])


def test_database_publication_is_atomic_and_preserves_old_snapshot(ready, tmp_path, monkeypatch):
    old_manifest = (ready / "manifest.json").read_bytes()

    def child(arguments, directory, *args):
        assert arguments[:2] == ["db", "download"]
        assert "--force" not in arguments
        destination = Path(arguments[arguments.index("--db-dir") + 1])
        assert destination != ready
        database(destination)
        return arguments, "Download complete"

    monkeypatch.setattr(runtime, "_run_child", child)
    result = runtime.update_databases(ready, ["ncbi"])
    published = Path(result["database_root"])
    assert published.is_dir() and published != ready
    assert published.parent == ready.parent
    assert (ready / "manifest.json").read_bytes() == old_manifest
    assert result["reference_snapshot"]["root"] == str(published)
    assert not list(tmp_path.glob(".wmlstudio-hydra-stage-*"))


def test_cancelled_database_build_publishes_nothing_and_preserves_old(ready, tmp_path, monkeypatch):
    original = sorted(tmp_path.iterdir())

    def child(arguments, directory, *args):
        database(Path(arguments[arguments.index("--db-dir") + 1]))
        raise AnalysisCancelled("User cancelled")

    monkeypatch.setattr(runtime, "_run_child", child)
    with pytest.raises(AnalysisCancelled):
        runtime.update_databases(ready, ["ncbi"])
    assert sorted(tmp_path.iterdir()) == original
    assert runtime.installed_databases(ready)["ncbi"]["version"] == "test-snapshot"


def test_child_failure_keeps_diagnostics_and_uses_argument_array(tmp_path, monkeypatch):
    monkeypatch.setattr(runtime, "_worker_command", lambda args: [
        sys.executable, "-c", "import sys; print('native search diagnostic'); sys.exit(7)"])
    with pytest.raises(runtime.HydraRuntimeError, match="native search diagnostic"):
        runtime._run_child([], tmp_path)


def test_child_cancellation_terminates_owned_process(tmp_path, monkeypatch):
    monkeypatch.setattr(runtime, "_worker_command", lambda args: [
        sys.executable, "-c", "import time; time.sleep(60)"])
    started = time.monotonic()
    with pytest.raises(AnalysisCancelled):
        runtime._run_child([], tmp_path, cancelled=lambda: time.monotonic() - started > 0.2)
    assert time.monotonic() - started < 15
