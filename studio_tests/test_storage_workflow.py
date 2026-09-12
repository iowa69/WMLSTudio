import hashlib
from pathlib import Path

import pytest

from wmlstudio.project import Project
from wmlstudio.sequence import AnalysisCancelled
from wmlstudio.storage import assign_organism, import_samples, organize_sample, safe_component


def sequence(path, text="ACGTACGT"):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(f">contig\n{text}\n")
    return path


@pytest.fixture
def project(tmp_path):
    with Project(tmp_path / "study.wmlstudio") as project:
        yield project


def test_assigning_organism_does_not_discard_inputless_external_profile(project):
    profile = {"scheme": "external MLST", "scheme_digest": "snapshot",
               "alleles": {"arcA": "1"}, "st": "131", "status": "profile_imported"}
    sid = project.add_profile("external isolate", profile)
    before = project.get_sample(sid)["result"]
    assign_organism(project, [sid], "Escherichia", "coli")
    sample = project.get_sample(sid)
    assert sample["result"] == before and sample["status"] == "completed"
    assert sample["metadata"]["organism"]["genus"] == "Escherichia"
    assert sample["profile_only"] and not sample["input_path"]


def test_same_stems_import_without_collisions_or_original_changes(project, tmp_path):
    a = sequence(tmp_path / "run1" / "same.fasta", "ACGT")
    b = sequence(tmp_path / "run2" / "same.fasta", "AAAA")
    before = {a: a.read_bytes(), b: b.read_bytes()}
    identifiers = import_samples(project, [{"path": path, "typing_mode": "manual", "genus": "Klebsiella",
                                            "species": "pneumoniae"} for path in (a, b)], tmp_path / "managed")
    assert len(set(identifiers)) == 2
    samples = project.samples()
    assert samples[0]["input_path"] != samples[1]["input_path"]
    for sample, original in zip(samples, (a, b), strict=True):
        copied = Path(sample["input_path"])
        assert copied.read_bytes() == before[original] == original.read_bytes()
        assert "Klebsiella/pneumoniae/ST_unassigned" in copied.as_posix()
        assert sample["id"] in copied.parts
        assert sample["metadata"]["workflow"]["source_sha256"] == hashlib.sha256(before[original]).hexdigest()


def test_typing_organises_managed_copy_and_appends_st_only_to_copy(project, tmp_path):
    original = sequence(tmp_path / "original" / "isolate.fasta")
    sid = import_samples(project, [{"path": original, "typing_mode": "manual", "genus": "Escherichia",
                                    "species": "coli"}], tmp_path / "managed", append_st=True)[0]
    previous = Path(project.get_sample(sid)["input_path"])
    project.set_result(sid, {"status": "complete", "st": "131", "scheme": "MLST", "scheme_digest": "abc",
                             "alleles": {"adk": "1"}, "input_sha256": hashlib.sha256(original.read_bytes()).hexdigest()})
    destination = organize_sample(project, sid)
    assert "Escherichia/coli/ST_131" in destination.as_posix()
    assert destination.name == "isolate_ST_131.fasta"
    assert destination.read_bytes() == original.read_bytes()
    assert original.name == "isolate.fasta" and original.exists()
    assert not previous.exists()
    assert project.get_sample(sid)["result"]["st"] == "131"
    assert organize_sample(project, sid) == destination
    assert any(event["action"] == "input_relocated" for event in project.history(sid))


def test_auto_identification_controls_folders_without_rewriting_manual_assignment(project, tmp_path):
    original = sequence(tmp_path / "isolate.fasta")
    sid = import_samples(project, [{"path": original, "typing_mode": "auto"}], tmp_path / "managed")[0]
    project.set_result(sid, {"status": "complete", "st": "7", "alleles": {},
                             "identification": {"genus": "Listeria", "species": "monocytogenes"}})
    path = organize_sample(project, sid)
    assert "Listeria/monocytogenes/ST_7" in path.as_posix()
    assert project.get_sample(sid)["metadata"]["organism"] == {"genus": "", "species": ""}


def test_cancelled_copy_leaves_no_sample_or_partial_file(project, tmp_path):
    original = sequence(tmp_path / "large.fasta", "A" * (3 * 1024 * 1024))
    before = original.read_bytes()
    checks = 0

    def cancel():
        nonlocal checks
        checks += 1
        return checks > 3

    with pytest.raises(AnalysisCancelled):
        import_samples(project, [{"path": original}], tmp_path / "managed", cancelled=cancel)
    assert project.samples() == []
    assert original.read_bytes() == before
    assert not [path for path in (tmp_path / "managed").rglob("*") if path.is_file()]


def test_database_failure_rolls_back_batch_and_removes_only_new_copies(project, tmp_path, monkeypatch):
    originals = [sequence(tmp_path / f"{i}.fasta") for i in range(2)]
    setter = project.set_metadata
    calls = 0

    def fail_second(sample_id, metadata):
        nonlocal calls
        calls += 1
        if calls == 2:
            raise OSError("Database unavailable")
        setter(sample_id, metadata)

    monkeypatch.setattr(project, "set_metadata", fail_second)
    with pytest.raises(OSError):
        import_samples(project, [{"path": path} for path in originals], tmp_path / "managed")
    assert project.samples() == []
    assert all(path.exists() for path in originals)
    assert not [path for path in (tmp_path / "managed").rglob("*") if path.is_file()]


def test_notification_failure_after_commit_cannot_delete_saved_input(project, tmp_path):
    source = sequence(tmp_path / "original.fasta")

    def progress(done, total, message):
        if done == total:
            raise RuntimeError("UI notification failed")

    with pytest.raises(RuntimeError):
        import_samples(project, [{"path": source}], tmp_path / "managed", progress=progress)
    assert len(project.samples()) == 1
    assert Path(project.samples()[0]["input_path"]).read_bytes() == source.read_bytes()


def test_source_changed_during_copy_is_rejected(project, tmp_path):
    source = sequence(tmp_path / "source.fasta")
    checks = 0

    def change_source():
        nonlocal checks
        checks += 1
        if checks == 3:
            source.write_text(">changed\nAAAAA\n")
        return False

    with pytest.raises(ValueError, match="changed"):
        import_samples(project, [{"path": source}], tmp_path / "managed", cancelled=change_source)
    assert project.samples() == []
    assert source.read_text() == ">changed\nAAAAA\n"


@pytest.mark.parametrize("name", ["CON", "aux.txt", "NUL.", "LPT1", "../../escape", "a\\b", "foo:bar", "\x00evil", "..."])
def test_windows_safe_components(name):
    component = safe_component(name)
    assert component and component not in {".", ".."}
    assert not any(char in component for char in '<>:"/\\|?*\x00')
    assert component.split(".", 1)[0].upper() not in {"CON", "AUX", "NUL", "LPT1"}
    assert not component.endswith((" ", "."))


def test_symlinked_storage_subdirectory_cannot_escape_root(project, tmp_path):
    source = sequence(tmp_path / "input.fasta")
    root = tmp_path / "managed"
    root.mkdir()
    outside = tmp_path / "outside"
    outside.mkdir()
    try:
        (root / "Genus").symlink_to(outside, target_is_directory=True)
    except OSError:
        pytest.skip("Directory symlinks unavailable")
    with pytest.raises(ValueError, match="escapes"):
        import_samples(project, [{"path": source, "genus": "Genus", "typing_mode": "manual"}], root)
    assert list(outside.iterdir()) == []
    assert project.samples() == []


def test_assignment_retains_custom_metadata_and_archives_previous_result(project, tmp_path):
    sample_id = project.add_sample(sequence(tmp_path / "input.fasta"))
    project.set_metadata(sample_id, {"ward": "ICU", "workflow": {"source_path": "original"}})
    project.set_result(sample_id, {"st": "1", "status": "complete"})
    assign_organism(project, [sample_id], "Klebsiella", "pneumoniae", "scheme", "manual")
    sample = project.get_sample(sample_id)
    assert sample["metadata"]["ward"] == "ICU"
    assert sample["metadata"]["workflow"]["source_path"] == "original"
    assert sample["metadata"]["workflow"]["scheme_path"] == "scheme"
    assert sample["result"] is None and sample["status"] == "queued"
    event = next(event for event in project.history(sample_id) if event["action"] == "result_invalidated")
    assert event["details"]["result"]["st"] == "1"
