"""Native worker contracts exercised with real Qt threads and sequence engines."""

import gzip
import threading
from pathlib import Path

import pytest
from PySide6.QtTest import QSignalSpy

from wmlstudio import __version__, jobs
from wmlstudio.typing import load_scheme

ARC1 = "AACCGTACGTTAG"
ARC2 = "AACCGTTCGTTAG"
GYR1 = "TTGGCATACCTGA"


@pytest.fixture
def scheme_path(tmp_path):
    path = tmp_path / "known_scheme"
    path.mkdir()
    (path / "arcA.tfa").write_text(f">arcA_1\n{ARC1}\n>arcA_2\n{ARC2}\n")
    (path / "gyrB.tfa").write_text(f">gyrB_1\n{GYR1}\n")
    (path / "profiles.tsv").write_text("ST\tarcA\tgyrB\n7\t1\t1\n8\t2\t1\n")
    return path


def sample(tmp_path, name, data, identifier="sample-1"):
    path = tmp_path / name
    path.write_text(data)
    return {"id": identifier, "name": f"Display {identifier}", "input_path": str(path)}


def run_worker(qtbot, worker, signal_names, action=None):
    spies = {name: QSignalSpy(getattr(worker, name)) for name in signal_names}
    try:
        with qtbot.waitSignal(worker.finished, timeout=10000):
            worker.start()
            if action is not None:
                action(worker)
    finally:
        if worker.isRunning():
            worker.cancel()
            assert worker.wait(10000), "Worker did not stop after cancellation."
    return {name: [spy.at(index) for index in range(spy.count())]
            for name, spy in spies.items()}


ANALYSIS_SIGNALS = (
    "sample_started", "sample_finished", "sample_failed", "sample_cancelled", "progress",
)
IMPORT_SIGNALS = ("imported", "failed", "progress")


def test_analysis_dispatches_by_contents_and_preserves_display_identity(
    qtbot, tmp_path, scheme_path,
):
    assembly = sample(tmp_path, "assembly.fastq", f">a\n{ARC1}\n>b\n{GYR1}\n", "assembly")
    reads = sample(tmp_path, "reads.fasta", "@one/1\nACGT\n+\nIIII\n", "reads")
    worker = jobs.AnalysisWorker([assembly, reads], scheme_path=scheme_path)
    signals = run_worker(qtbot, worker, ANALYSIS_SIGNALS)
    assert signals["sample_failed"] == []
    assert signals["sample_cancelled"] == []
    results = {identifier: result for identifier, result in signals["sample_finished"]}
    assert results["assembly"]["kind"] == "fasta"
    assert results["assembly"]["status"] == "complete"
    assert results["assembly"]["st"] == "7"
    assert results["reads"]["kind"] == "fastq"
    assert results["reads"]["status"] == "qc_only"
    assert results["reads"]["st"] is None
    assert results["reads"]["alleles"] == {}
    assert any("Raw reads have not been assembled or typed" in note
               for note in results["reads"]["notes"])
    assert not any("Assembly quality only" in note for note in results["reads"]["notes"])
    for identifier, result in results.items():
        assert result["sample_id"] == identifier
        assert result["sample_name"] == f"Display {identifier}"
        assert result["software_version"] == __version__


def test_invalid_sample_does_not_stop_next_sample(qtbot, tmp_path):
    bad = sample(tmp_path, "bad.fastq", "@r\nACGT\n+\nII\n", "bad")
    good = sample(tmp_path, "good.fasta", ">a\nACGTNN\n", "good")
    signals = run_worker(qtbot, jobs.AnalysisWorker([bad, good]), ANALYSIS_SIGNALS)
    assert signals["sample_started"] == [["bad"], ["good"]]
    assert len(signals["sample_failed"]) == 1
    assert signals["sample_failed"][0][0] == "bad"
    assert "truncated FASTQ" in signals["sample_failed"][0][1]
    assert len(signals["sample_finished"]) == 1
    identifier, result = signals["sample_finished"][0]
    assert identifier == "good"
    assert result["qc"]["total_bases"] == 6
    assert result["status"] == "qc_only"
    assert signals["progress"][-1] == [100, "Analysis finished"]


def test_per_sample_unknown_and_manual_no_scheme_do_not_force_global_typing(qtbot, tmp_path, scheme_path):
    unknown = sample(tmp_path, 'unknown.fa', f'>a\n{ARC1}\n>b\n{GYR1}\n', 'unknown')
    unknown['metadata'] = {'workflow': {'typing_mode': 'unknown'}}
    manual = sample(tmp_path, 'manual.fa', f'>a\n{ARC1}\n>b\n{GYR1}\n', 'manual')
    manual['metadata'] = {'organism': {'genus': 'Staphylococcus', 'species': 'epidermidis'},
                          'workflow': {'typing_mode': 'manual', 'scheme_path': None}}
    signals = run_worker(qtbot, jobs.AnalysisWorker([unknown, manual], scheme_path), ANALYSIS_SIGNALS)
    assert signals['sample_failed'] == []
    results = dict(signals['sample_finished'])
    assert results['unknown']['status'] == 'qc_only'
    assert results['manual']['status'] == 'qc_only'
    assert results['manual']['organism_assignment'] == 'user supplied'
    assert results['manual']['organism']['species'] == 'epidermidis'


def test_per_sample_full_cds_and_exact_route_explicitly(qtbot, tmp_path, scheme_path, monkeypatch):
    import copy
    calls = []
    real = jobs.call_assembly
    def cg(path, scheme, cancelled, progress, genetic_code):
        calls.append('full_cds')
        assert callable(cancelled)
        assert genetic_code == 11
        return real(path, scheme, cancelled=cancelled, progress=progress)
    monkeypatch.setattr(jobs, 'call_cgassembly', cg)
    first = sample(tmp_path, 'sample.fa', f'>a\n{ARC1}\n>b\n{GYR1}\n', 'first')
    first['metadata'] = {'workflow': {'typing_mode': 'manual', 'scheme_path': str(scheme_path),
                                     'calling_mode': 'full_cds'}}
    second = copy.deepcopy(first)
    second['id'] = 'second'
    second['metadata']['workflow']['calling_mode'] = 'exact'
    signals = run_worker(qtbot, jobs.AnalysisWorker([first, second]), ANALYSIS_SIGNALS)
    assert signals['sample_failed'] == []
    assert calls == ['full_cds']
    assert [result['st'] for _, result in signals['sample_finished']] == ['7', '7']


def test_scheme_failure_does_not_prevent_read_quality(qtbot, tmp_path):
    assembly = sample(tmp_path, "assembly.fasta", ">a\nACGT\n", "assembly")
    reads = sample(tmp_path, "reads.fastq", "@r\nACGT\n+\nIIII\n", "reads")
    worker = jobs.AnalysisWorker([assembly, reads], scheme_path=tmp_path / "missing_scheme")
    signals = run_worker(qtbot, worker, ANALYSIS_SIGNALS)
    assert signals["sample_failed"][0][0] == "assembly"
    assert "Scheme could not be loaded" in signals["sample_failed"][0][1]
    assert signals["sample_finished"][0][0] == "reads"
    assert signals["sample_finished"][0][1]["qc"]["q30_percent"] == 100


def test_worker_reports_bounded_fastq_prefix(qtbot, tmp_path):
    reads = sample(
        tmp_path, "reads.fastq", "@r1\nACGT\n+\nIIII\n@r2\nTGCA\n+\n!!!!\n",
    )
    signals = run_worker(qtbot, jobs.AnalysisWorker([reads], max_reads=1), ANALYSIS_SIGNALS)
    result = signals["sample_finished"][0][1]
    assert result["qc"]["records"] == 1
    assert result["qc"]["sampled"] is True
    assert result["qc"]["q30_percent"] == 100
    assert any("not validated" in note for note in result["notes"])


def test_scheme_load_observes_worker_cancellation(qtbot, tmp_path, scheme_path, monkeypatch):
    entered, released = threading.Event(), threading.Event()
    real_load = jobs.load_scheme

    def gated_load(path, cancelled=None):
        entered.set()
        assert released.wait(5), "Test did not release the scheme loader."
        return real_load(path, cancelled=cancelled)

    monkeypatch.setattr(jobs, "load_scheme", gated_load)
    row = sample(tmp_path, "assembly.fasta", f">a\n{ARC1}\n", "unstarted")
    worker = jobs.AnalysisWorker([row], scheme_path=scheme_path)

    def cancel_during_load(active):
        qtbot.waitUntil(entered.is_set, timeout=3000)
        active.cancel()
        released.set()

    signals = run_worker(qtbot, worker, ANALYSIS_SIGNALS, cancel_during_load)
    assert signals["sample_started"] == []
    assert signals["sample_finished"] == []
    assert signals["sample_failed"] == []
    assert signals["progress"][-1] == [100, "Analysis cancelled"]


def test_cancelled_active_sample_never_completes_or_starts_next(
    qtbot, tmp_path, monkeypatch,
):
    entered, released = threading.Event(), threading.Event()
    real_inspect = jobs.inspect_sequence

    def gated_inspect(path, max_reads, cancelled):
        entered.set()
        assert released.wait(5), "Test did not release the sequence inspector."
        return real_inspect(path, max_reads=max_reads, cancelled=cancelled)

    monkeypatch.setattr(jobs, "inspect_sequence", gated_inspect)
    rows = [sample(tmp_path, f"{name}.fasta", ">a\nACGT\n", name) for name in ("active", "next")]
    worker = jobs.AnalysisWorker(rows)

    def cancel_during_inspection(active):
        qtbot.waitUntil(entered.is_set, timeout=3000)
        active.cancel()
        released.set()

    signals = run_worker(qtbot, worker, ANALYSIS_SIGNALS, cancel_during_inspection)
    assert signals["sample_started"] == [["active"]]
    assert signals["sample_cancelled"] == [["active"]]
    assert signals["sample_finished"] == []
    assert signals["sample_failed"] == []
    assert signals["progress"][-1] == [100, "Analysis cancelled"]


def test_import_validates_copied_scheme_and_preserves_source(qtbot, tmp_path, scheme_path):
    root = tmp_path / "workspace"
    root.mkdir()
    original = {path.name: path.read_bytes() for path in scheme_path.iterdir()}
    signals = run_worker(qtbot, jobs.SchemeImportWorker(scheme_path, root), IMPORT_SIGNALS)
    assert signals["failed"] == []
    assert len(signals["imported"]) == 1
    destination, locus_count = signals["imported"][0]
    assert locus_count == 2
    assert Path(destination).parent == root / "schemes"
    assert load_scheme(destination).digest == load_scheme(scheme_path).digest
    assert {path.name: path.read_bytes() for path in scheme_path.iterdir()} == original
    assert not list(root.glob("scheme-*"))
    repeated = run_worker(qtbot, jobs.SchemeImportWorker(scheme_path, root), IMPORT_SIGNALS)
    assert repeated["failed"] == []
    assert repeated["imported"] == signals["imported"]


def test_import_supports_compressed_allele_files(qtbot, tmp_path, scheme_path):
    allele = scheme_path / "arcA.tfa"
    compressed = scheme_path / "arcA.tfa.gz"
    compressed.write_bytes(gzip.compress(allele.read_bytes()))
    allele.unlink()
    root = tmp_path / "workspace"
    root.mkdir()
    signals = run_worker(qtbot, jobs.SchemeImportWorker(scheme_path, root), IMPORT_SIGNALS)
    assert signals["failed"] == []
    assert len(signals["imported"]) == 1
    destination = Path(signals["imported"][0][0])
    assert (destination / "arcA.tfa.gz").read_bytes() == compressed.read_bytes()
    assert load_scheme(destination).digest == load_scheme(scheme_path).digest


def test_import_rejects_symlinked_allele_file(qtbot, tmp_path, scheme_path):
    original = scheme_path / "arcA.tfa"
    external = tmp_path / "external.tfa"
    original.rename(external)
    try:
        original.symlink_to(external)
    except (OSError, NotImplementedError):
        pytest.skip("This environment cannot create symbolic links.")
    root = tmp_path / "workspace"
    root.mkdir()
    signals = run_worker(qtbot, jobs.SchemeImportWorker(scheme_path, root), IMPORT_SIGNALS)
    assert signals["imported"] == []
    assert len(signals["failed"]) == 1
    assert "symbolic" in signals["failed"][0][0].lower()
    assert not list(root.glob("scheme-*"))
    assert external.is_file() and original.is_symlink()


def test_import_rejects_source_changed_while_copying(
    qtbot, tmp_path, scheme_path, monkeypatch,
):
    real_copy = jobs.shutil.copy2

    def change_before_copy(source, destination):
        if Path(source).name == "arcA.tfa":
            source.write_text(source.read_text().replace(ARC1, "T" + ARC1[1:]))
        return real_copy(source, destination)

    monkeypatch.setattr(jobs.shutil, "copy2", change_before_copy)
    root = tmp_path / "workspace"
    root.mkdir()
    signals = run_worker(qtbot, jobs.SchemeImportWorker(scheme_path, root), IMPORT_SIGNALS)
    assert signals["imported"] == []
    assert len(signals["failed"]) == 1
    assert "changed" in signals["failed"][0][0].lower()
    assert not list(root.glob("scheme-*"))
    assert not (root / "schemes").exists()


def test_import_revalidates_preexisting_destination(qtbot, tmp_path, scheme_path):
    root = tmp_path / "workspace"
    root.mkdir()
    first = run_worker(qtbot, jobs.SchemeImportWorker(scheme_path, root), IMPORT_SIGNALS)
    destination = Path(first["imported"][0][0])
    reference = destination / "arcA.tfa"
    tampered = reference.read_text().replace(ARC1, "T" + ARC1[1:])
    reference.write_text(tampered)
    again = run_worker(qtbot, jobs.SchemeImportWorker(scheme_path, root), IMPORT_SIGNALS)
    assert again["imported"] == []
    assert len(again["failed"]) == 1
    assert reference.read_text() == tampered


def test_import_cancellation_during_copy_removes_staging_only(
    qtbot, tmp_path, scheme_path, monkeypatch,
):
    root = tmp_path / "workspace"
    root.mkdir()
    unrelated = root / "keep.txt"
    unrelated.write_text("existing project content")
    worker = jobs.SchemeImportWorker(scheme_path, root)
    real_copy = jobs.shutil.copy2

    def cancel_after_copy(source, destination):
        result = real_copy(source, destination)
        worker.cancel()
        return result

    monkeypatch.setattr(jobs.shutil, "copy2", cancel_after_copy)
    signals = run_worker(qtbot, worker, IMPORT_SIGNALS)
    assert signals["imported"] == []
    assert signals["failed"] == []
    assert signals["progress"][-1] == [100, "Scheme import cancelled"]
    assert list(root.iterdir()) == [unrelated]
    assert unrelated.read_text() == "existing project content"
    assert load_scheme(scheme_path).locus_count == 2
