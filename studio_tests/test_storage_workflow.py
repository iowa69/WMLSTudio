import hashlib
from pathlib import Path

import pytest

from wmlstudio.project import Project
from wmlstudio.sequence import AnalysisCancelled
from wmlstudio.storage import (
    assign_organism,
    import_samples,
    organize_sample,
    plan_filing,
    safe_component,
)


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


def test_organising_by_st_prunes_the_directory_it_vacated(project, tmp_path):
    original = sequence(tmp_path / "original" / "isolate.fasta")
    root = tmp_path / "managed"
    sid = import_samples(project, [{"path": original, "typing_mode": "manual", "genus": "Escherichia",
                                    "species": "coli"}], root)[0]
    vacated = Path(project.get_sample(sid)["input_path"]).parent
    project.set_result(sid, {"status": "complete", "st": "131", "alleles": {},
                             "input_sha256": hashlib.sha256(original.read_bytes()).hexdigest()})
    destination = organize_sample(project, sid)
    assert "Escherichia/coli/ST_131" in destination.as_posix()
    assert not vacated.exists() and not vacated.parent.exists()
    assert (root / "Escherichia" / "coli").is_dir()


def test_recorded_organism_evidence_travels_with_the_import_without_changing_filing(project, tmp_path):
    original = sequence(tmp_path / "isolate.fasta")
    evidence = {"status": "confirmed", "basis": "genomic_ani", "confidence": "genomic_reference_supported",
                "proposed": {"genus": "Klebsiella", "species": "pneumoniae"},
                "accepted": {"genus": "Klebsiella", "species": "pneumoniae"}}
    sid = import_samples(project, [{"path": original, "typing_mode": "manual", "genus": "Klebsiella",
                                    "species": "pneumoniae", "organism_evidence": evidence}],
                         tmp_path / "managed")[0]
    stored = project.get_sample(sid)["metadata"]["organism_evidence"]
    assert stored["basis"] == "genomic_ani" and stored["format_version"] == 1
    assert "Klebsiella/pneumoniae/ST_unassigned" in Path(project.get_sample(sid)["input_path"]).as_posix()
    assert plan_filing(project, sid)["changed"] is False
    assert any(event["action"] == "organism_identified" for event in project.history(sid))


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


def typed(project, sample_id, digest):
    """One classical ST and one core-genome profile, stored the way the app stores them."""
    project.set_result(sample_id, {"scheme": "Klebsiella MLST", "scheme_digest": "mlst-v1",
                                   "st": "258", "status": "complete", "input_sha256": digest,
                                   "alleles": {f"gene{index}": "1" for index in range(7)}})
    project.set_analysis(sample_id, {"scheme": "Klebsiella cgMLST", "scheme_digest": "cg-v1",
                                     "st": None, "status": "complete", "input_sha256": digest,
                                     "alleles": {f"locus{index:04d}": "3" for index in range(400)},
                                     "parameters": {"method": "full-cds-cgmlst-v2"}})


def hydra_report(digest):
    return {"samples": [{"sample": "isolate", "summary": {"amr_genes": 2},
                         "hits": [{"gene": "blaKPC-2", "element_type": "AMR", "class": "CARBAPENEM",
                                   "primary": True}],
                         "mlst": {}, "species": {}, "input_sha256": digest}],
            "import_provenance": {"sha256": "a" * 64}, "execution_provenance": {}, "databases": []}


def test_a_saved_project_reopens_with_every_sample_result_and_evidence_still_usable(tmp_path):
    from wmlstudio.sample_workflow import hydra_evidence_status, link_hydra
    from wmlstudio.storage import confirm_organism
    path = tmp_path / "study.wmlstudio"
    original = sequence(tmp_path / "inbox" / "isolate.fasta", "ACGTACGTAC")
    root = tmp_path / "study.files"
    with Project(path) as project:
        sid = import_samples(project, [{"path": original, "typing_mode": "manual",
                                        "genus": "Klebsiella", "species": "pneumoniae"}], root)[0]
        copied = Path(project.get_sample(sid)["input_path"])
        digest = hashlib.sha256(copied.read_bytes()).hexdigest()
        typed(project, sid, digest)
        confirm_organism(project, [sid], "Klebsiella", "pneumoniae")
        link_hydra(project, hydra_report(digest), {"isolate": sid})
        project.create_collection("Ward A")
        before = project.get_sample(sid)
        planned = plan_filing(project, sid)

    with Project(path) as project:
        sample = project.get_sample(sid)
        assert sample["name"] == before["name"] and sample["status"] == "completed"
        assert sample["input_path"] == str(copied) and Path(sample["input_path"]).is_file()
        assert sample["missing_input"] is False
        assert sample["metadata"]["organism"] == {"genus": "Klebsiella", "species": "pneumoniae"}
        assert sample["metadata"]["organism_evidence"]["status"] == "confirmed"
        assert sample["metadata"]["workflow"]["managed"] is True
        assert project.latest_analysis(sid, "mlst")["st"] == "258"
        assert len(project.latest_analysis(sid, "cgmlst")["alleles"]) == 400
        assert hydra_evidence_status(sample)["status"] == "current"
        assert sample["metadata"]["hydra"]["summary"]["amr_genes"] == 2
        assert project.collections()[0]["name"] == "Ward A"
        # Immediately usable means the reopened project plans exactly what the open
        # one planned: the same managed root, the same destination, the same ST.
        assert plan_filing(project, sid) == planned
        assert planned["root"] == root and planned["st"] == "258"


def test_a_project_that_moved_with_its_folder_finds_its_managed_copies_again(tmp_path):
    from wmlstudio.storage import (
        missing_managed_copies,
        relocate_managed_storage,
    )
    original = sequence(tmp_path / "inbox" / "isolate.fasta", "ACGTACGTAC")
    first = tmp_path / "drive_e"
    first.mkdir()
    with Project(first / "study.wmlstudio") as project:
        sid = import_samples(project, [{"path": original, "typing_mode": "manual",
                                        "genus": "Klebsiella"}], first / "study.files")[0]
        digest = hashlib.sha256(Path(project.get_sample(sid)["input_path"]).read_bytes()).hexdigest()
        typed(project, sid, digest)
        relative = Path(project.get_sample(sid)["input_path"]).relative_to(first / "study.files")

    moved = tmp_path / "drive_f"
    first.rename(moved)
    with Project(moved / "study.wmlstudio") as project:
        outstanding = missing_managed_copies(project)
        assert [entry["sample_id"] for entry in outstanding] == [sid]
        assert project.get_sample(sid)["missing_input"] is True

        report = relocate_managed_storage(project)
        assert report["relocated"] == [(sid, str(moved / "study.files" / relative))]
        assert report["unresolved"] == []
        sample = project.get_sample(sid)
        assert sample["missing_input"] is False
        assert Path(sample["input_path"]).read_bytes() == original.read_bytes()
        assert sample["metadata"]["workflow"]["storage_root"] == str(moved / "study.files")
        assert project.latest_analysis(sid, "mlst")["st"] == "258"
        assert missing_managed_copies(project) == []
        event = next(e for e in project.history(sid) if e["action"] == "managed_storage_relocated")
        assert event["details"]["sha256"] == digest


def test_a_file_with_the_right_name_but_other_contents_is_never_adopted_as_the_copy(tmp_path):
    from wmlstudio.storage import relocate_managed_storage
    original = sequence(tmp_path / "inbox" / "isolate.fasta", "ACGTACGTAC")
    first = tmp_path / "before"
    first.mkdir()
    with Project(first / "study.wmlstudio") as project:
        sid = import_samples(project, [{"path": original, "typing_mode": "manual",
                                        "genus": "Klebsiella"}], first / "study.files")[0]
        relative = Path(project.get_sample(sid)["input_path"]).relative_to(first / "study.files")

    moved = tmp_path / "after"
    first.rename(moved)
    impostor = moved / "study.files" / relative
    impostor.write_text(">contig\nTTTTTTTTTT\n")
    with Project(moved / "study.wmlstudio") as project:
        stale = project.get_sample(sid)["input_path"]
        report = relocate_managed_storage(project)
        assert report["relocated"] == []
        assert "contents differ" in report["unresolved"][0][1]
        assert project.get_sample(sid)["input_path"] == stale
        assert impostor.read_text() == ">contig\nTTTTTTTTTT\n"
