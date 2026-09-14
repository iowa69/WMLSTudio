import hashlib
import json
import random
from pathlib import Path

import pytest

from wmlstudio.characterization_refs import reference_digest
from wmlstudio.organism_id import assignment_for, identify_batch
from wmlstudio.project import Project
from wmlstudio.sequence import file_sha256
from wmlstudio.storage import (
    QUARANTINE_ROOT,
    assign_organism,
    confirm_organism,
    delete_managed_copy,
    import_samples,
    managed_copy_for,
    organize_sample,
    plan_filing,
    preview_target,
    quarantine_samples,
    reassign_organism,
    refile_samples,
    remove_samples,
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


def quarantined(path, reason, **evidence):
    return {"path": path, "typing_mode": "auto", "quarantine": reason,
            "organism_evidence": {"status": "quarantined", "quarantine_reason": reason,
                                  "basis": "none", "confidence": "unresolved", **evidence}}


def test_unidentified_input_is_filed_for_review_and_never_into_a_guessed_genus(project, tmp_path):
    original = sequence(tmp_path / "inbox" / "mystery.fasta")
    root = tmp_path / "managed"
    sid = import_samples(project, [quarantined(original, "not_in_reference_panel")], root)[0]
    copied = Path(project.get_sample(sid)["input_path"])
    relative = copied.relative_to(root)
    assert relative.parts[:2] == (QUARANTINE_ROOT, "Not_in_reference_panel")
    assert sid in relative.parts
    assert copied.read_bytes() == original.read_bytes() and original.exists()
    assert not (root / "Unknown_genus").exists()
    notice = (root / QUARANTINE_ROOT / "README.txt").read_text(encoding="utf-8")
    assert "declined to decide" in notice and "not because the organism is new" in notice
    assert "Not_in_reference_panel" in notice
    event = next(e for e in project.history(sid) if e["action"] == "organism_quarantined")
    assert event["details"]["reason"] == "not_in_reference_panel"


def test_needs_review_folder_cannot_collide_with_an_organism_folder(project, tmp_path):
    original = sequence(tmp_path / "input.fasta")
    root = tmp_path / "managed"
    sid = import_samples(project, [{"path": original, "typing_mode": "manual",
                                    "genus": QUARANTINE_ROOT, "species": "pneumoniae"}], root)[0]
    filed = Path(project.get_sample(sid)["input_path"]).relative_to(root)
    assert filed.parts[0] != QUARANTINE_ROOT
    assert filed.parts[0].startswith(QUARANTINE_ROOT)


def test_confirming_an_organism_refiles_the_copy_and_prunes_the_empty_bucket(project, tmp_path):
    original = sequence(tmp_path / "isolate.fasta")
    root = tmp_path / "managed"
    sid = import_samples(project, [quarantined(original, "awaiting_identification")], root)[0]
    before = Path(project.get_sample(sid)["input_path"])
    bucket = before.parent.parent
    confirm_organism(project, [sid], "Klebsiella", "pneumoniae")
    plan = plan_filing(project, sid)
    assert plan["changed"] and plan["quarantine"] is None
    assert plan["relative"] == "Klebsiella/pneumoniae/ST_unassigned/%s/isolate.fasta" % sid
    report = refile_samples(project, [sid])
    assert [identifier for identifier, _ in report["moved"]] == [sid]
    after = Path(project.get_sample(sid)["input_path"])
    assert after == plan["destination"] and after.read_bytes() == original.read_bytes()
    assert not before.exists() and not before.parent.exists() and not bucket.exists()
    assert original.exists()
    assert project.get_sample(sid)["metadata"]["organism_evidence"]["status"] == "confirmed"


def test_explicit_unknown_mode_still_uses_unknown_genus_rather_than_the_review_tree(project, tmp_path):
    original = sequence(tmp_path / "unknown.fasta")
    root = tmp_path / "managed"
    sid = import_samples(project, [{"path": original, "typing_mode": "unknown"}], root)[0]
    filed = Path(project.get_sample(sid)["input_path"]).relative_to(root)
    assert filed.parts[:3] == ("Unknown_genus", "Unknown_species", "ST_unassigned")
    assert QUARANTINE_ROOT not in filed.parts


def test_sending_a_filed_sample_back_to_review_moves_it_and_records_why(project, tmp_path):
    original = sequence(tmp_path / "isolate.fasta")
    root = tmp_path / "managed"
    sid = import_samples(project, [{"path": original, "typing_mode": "manual",
                                    "genus": "Klebsiella", "species": "pneumoniae"}], root)[0]
    filed = Path(project.get_sample(sid)["input_path"])
    quarantine_samples(project, [sid], "user_deferred")
    refile_samples(project, [sid])
    moved = Path(project.get_sample(sid)["input_path"]).relative_to(root)
    assert moved.parts[:2] == (QUARANTINE_ROOT, "User_deferred")
    assert not filed.exists() and not (root / "Klebsiella").exists()


def test_duplicate_import_keeps_one_sample_and_leaves_no_orphan_copy(project, tmp_path):
    first = sequence(tmp_path / "run1" / "isolate.fasta", "ACGTACGTAA")
    second = sequence(tmp_path / "run2" / "copy_of_isolate.fasta", "ACGTACGTAA")
    root = tmp_path / "managed"
    notes = []
    identifiers = import_samples(
        project, [{"path": first, "typing_mode": "manual", "genus": "Klebsiella"},
                  {"path": second, "typing_mode": "manual", "genus": "Klebsiella"}],
        root, duplicates="skip", notes=notes)
    assert len(identifiers) == 1 and len(project.samples()) == 1
    assert notes[0]["duplicate_of"] == identifiers[0]
    assert notes[0]["sha256"] == hashlib.sha256(first.read_bytes()).hexdigest()
    copies = [path for path in root.rglob("*") if path.is_file()]
    assert len(copies) == 1
    assert second.exists() and first.exists()
    skipped = [e for e in project.history() if e["action"] == "duplicate_import_skipped"]
    assert len(skipped) == 1 and skipped[0]["sample_id"] is None


def test_allow_is_still_the_default_and_imports_both_identical_files(project, tmp_path):
    first = sequence(tmp_path / "run1" / "isolate.fasta", "ACGTACGTAA")
    second = sequence(tmp_path / "run2" / "isolate.fasta", "ACGTACGTAA")
    identifiers = import_samples(
        project, [{"path": path, "typing_mode": "manual", "genus": "Klebsiella"}
                  for path in (first, second)], tmp_path / "managed")
    assert len(identifiers) == 2 and len(project.samples()) == 2


def test_removing_a_sample_deletes_only_a_verified_managed_copy(project, tmp_path):
    original = sequence(tmp_path / "original" / "isolate.fasta")
    root = tmp_path / "managed"
    sid = import_samples(project, [{"path": original, "typing_mode": "manual",
                                    "genus": "Klebsiella"}], root)[0]
    copy = Path(project.get_sample(sid)["input_path"])
    report = remove_samples(project, [sid], delete_managed_copy=True)
    assert report["deleted"] == [str(copy)] and report["retained"] == []
    assert not copy.exists() and not copy.parent.exists()
    assert original.exists() and original.read_text().startswith(">contig")
    assert project.samples() == []


def test_a_changed_managed_copy_is_never_deleted_by_removal(project, tmp_path):
    original = sequence(tmp_path / "original" / "isolate.fasta")
    sid = import_samples(project, [{"path": original, "typing_mode": "manual",
                                    "genus": "Klebsiella"}], tmp_path / "managed")[0]
    copy = Path(project.get_sample(sid)["input_path"])
    copy.write_text(">edited\nAAAA\n")
    assert managed_copy_for(project, sid) is None
    report = remove_samples(project, [sid], delete_managed_copy=True)
    assert copy.exists() and report["deleted"] == []
    assert "no longer matches" in report["retained"][0][1]
    assert original.exists()


def test_an_unmanaged_original_is_never_deleted_by_removal(project, tmp_path):
    original = sequence(tmp_path / "original" / "isolate.fasta")
    sid = import_samples(project, [{"path": original, "typing_mode": "manual", "genus": "Klebsiella"}],
                         tmp_path / "managed", managed=False)[0]
    assert managed_copy_for(project, sid) is None
    report = remove_samples(project, [sid], delete_managed_copy=True)
    assert original.exists() and report["deleted"] == []
    assert "linked, not copied" in report["retained"][0][1]


def test_a_file_outside_the_storage_root_is_never_deleted_by_removal(project, tmp_path):
    original = sequence(tmp_path / "original" / "isolate.fasta")
    sid = import_samples(project, [{"path": original, "typing_mode": "manual", "genus": "Klebsiella"}],
                         tmp_path / "managed")[0]
    copy = Path(project.get_sample(sid)["input_path"])
    project.update_metadata(sid, {"workflow": {"storage_root": str(tmp_path / "elsewhere")}})
    assert managed_copy_for(project, sid) is None
    report = remove_samples(project, [sid], delete_managed_copy=True)
    assert copy.exists() and original.exists()
    assert report["deleted"] == [] and "outside the managed storage location" in report["retained"][0][1]


def test_deleting_a_managed_copy_relinks_to_the_original_and_refuses_without_one(project, tmp_path):
    original = sequence(tmp_path / "original" / "isolate.fasta")
    root = tmp_path / "managed"
    sid = import_samples(project, [{"path": original, "typing_mode": "manual",
                                    "genus": "Klebsiella"}], root)[0]
    copy = Path(project.get_sample(sid)["input_path"])
    assert managed_copy_for(project, sid) == copy
    report = delete_managed_copy(project, sid)
    assert report["deleted"] and not copy.exists()
    assert Path(project.get_sample(sid)["input_path"]) == original.resolve()
    assert project.get_sample(sid)["metadata"]["workflow"]["managed"] is False

    other = sequence(tmp_path / "second" / "isolate.fasta", "TTTT")
    second = import_samples(project, [{"path": other, "typing_mode": "manual",
                                       "genus": "Klebsiella"}], root)[0]
    other.unlink()
    kept = Path(project.get_sample(second)["input_path"])
    outcome = delete_managed_copy(project, second)
    assert outcome["deleted"] is False and kept.exists()
    assert "only remaining evidence" in outcome["reason"]


def test_correcting_a_label_keeps_a_result_earned_from_a_pinned_scheme(project, tmp_path):
    original = sequence(tmp_path / "pinned.fasta")
    sid = import_samples(project, [{"path": original, "typing_mode": "manual", "genus": "Escherichia",
                                    "species": "coli", "scheme_path": "scheme-A"}], tmp_path / "managed")[0]
    project.set_result(sid, {"status": "complete", "st": "131", "scheme_digest": "abc", "alleles": {}})
    assign_organism(project, [sid], "Klebsiella", "pneumoniae", "scheme-A", "manual")
    sample = project.get_sample(sid)
    assert sample["result"]["st"] == "131" and sample["status"] == "completed"
    assert sample["metadata"]["organism"]["genus"] == "Klebsiella"


def test_correcting_a_label_reruns_typing_when_the_scheme_follows_the_organism(project, tmp_path):
    original = sequence(tmp_path / "auto.fasta")
    sid = import_samples(project, [{"path": original, "typing_mode": "auto"}], tmp_path / "managed")[0]
    project.set_result(sid, {"status": "complete", "st": "7", "alleles": {}})
    assign_organism(project, [sid], "Listeria", "monocytogenes", None, "auto")
    assert project.get_sample(sid)["result"] is None
    assert project.get_sample(sid)["status"] == "queued"


def test_changing_the_pinned_scheme_still_invalidates_the_result(project, tmp_path):
    original = sequence(tmp_path / "pinned.fasta")
    sid = import_samples(project, [{"path": original, "typing_mode": "manual", "genus": "Escherichia",
                                    "scheme_path": "scheme-A"}], tmp_path / "managed")[0]
    project.set_result(sid, {"status": "complete", "st": "131", "alleles": {}})
    assign_organism(project, [sid], "Escherichia", "", "scheme-B", "manual")
    assert project.get_sample(sid)["result"] is None


def test_reassigning_an_organism_confirms_it_and_moves_the_file_in_one_call(project, tmp_path):
    original = sequence(tmp_path / "isolate.fasta")
    root = tmp_path / "managed"
    sid = import_samples(project, [quarantined(original, "low_confidence")], root)[0]
    report = reassign_organism(project, [sid], "Enterobacter", "", storage_root=root)
    assert len(report["moved"]) == 1
    filed = Path(project.get_sample(sid)["input_path"]).relative_to(root)
    assert filed.parts[:3] == ("Enterobacter", "Unknown_species", "ST_unassigned")
    evidence = project.get_sample(sid)["metadata"]["organism_evidence"]
    assert evidence["status"] == "confirmed" and evidence["basis"] == "user_assigned"
    assert evidence["confirmed_by"] == "user" and evidence["accepted"]["genus"] == "Enterobacter"


def test_the_plan_shown_before_a_move_is_the_path_the_move_creates(project, tmp_path):
    original = sequence(tmp_path / "isolate.fasta")
    root = tmp_path / "managed"
    sid = import_samples(project, [{"path": original, "typing_mode": "manual",
                                    "genus": "Klebsiella", "species": "pneumoniae"}], root)[0]
    project.set_result(sid, {"status": "complete", "st": "258", "alleles": {},
                             "input_sha256": hashlib.sha256(original.read_bytes()).hexdigest()})
    plan = plan_filing(project, sid)
    assert plan["eligible"] and plan["changed"] and plan["st"] == "258"
    assert organize_sample(project, sid) == plan["destination"]
    assert plan_filing(project, sid)["changed"] is False


def test_preview_target_matches_what_an_import_creates_for_every_tier(project, tmp_path):
    root = tmp_path / "managed"
    cases = [({"genus": "Klebsiella", "species": "pneumoniae"}, None),
             ({"genus": "Enterobacter", "species": ""}, None),
             ({"genus": "", "species": ""}, "conflicting_evidence"),
             ({"genus": "", "species": ""}, "reads_not_assembled")]
    for index, (organism, bucket) in enumerate(cases):
        original = sequence(tmp_path / f"case{index}.fasta", "ACGT" * (index + 1))
        assignment = {"path": original, "typing_mode": "auto", **organism}
        if bucket:
            assignment.update(quarantined(original, bucket))
        sid = import_samples(project, [assignment], root)[0]
        preview = preview_target(root, sid, original, organism["genus"], organism["species"],
                                 quarantine=bucket)
        assert Path(project.get_sample(sid)["input_path"]) == preview
        assert sid in preview.relative_to(root).parts


def test_a_quarantine_reason_the_storage_layer_does_not_know_is_refused(project, tmp_path):
    original = sequence(tmp_path / "isolate.fasta")
    with pytest.raises(ValueError, match="Unknown quarantine bucket"):
        import_samples(project, [{"path": original, "typing_mode": "auto",
                                  "quarantine": "probably_novel_species"}], tmp_path / "managed")
    assert project.samples() == []
    assert not list((tmp_path / "managed").rglob("*.fasta"))


def test_profile_only_samples_are_skipped_rather_than_failed_by_refiling(project):
    sid = project.add_profile("external isolate", {"scheme": "external", "scheme_digest": "x",
                                                   "alleles": {"adk": "1"}, "st": "7",
                                                   "status": "profile_imported"})
    report = refile_samples(project, [sid])
    assert report["moved"] == [] and report["skipped"][0][0] == sid
    assert "profile-only" in report["skipped"][0][1]


def test_every_managed_file_lives_under_one_of_the_three_known_shapes(project, tmp_path):
    root = tmp_path / "managed"
    assignments = [
        {"path": sequence(tmp_path / "a.fasta", "AAAA"), "typing_mode": "manual", "genus": "Klebsiella"},
        {"path": sequence(tmp_path / "b.fasta", "CCCC"), "typing_mode": "unknown"},
        quarantined(sequence(tmp_path / "c.fasta", "GGGG"), "low_confidence"),
    ]
    identifiers = import_samples(project, assignments, root)
    for sample in project.samples():
        parts = Path(sample["input_path"]).relative_to(root).parts
        assert sample["id"] in parts
        if parts[0] == QUARANTINE_ROOT:
            assert len(parts) == 4 and parts[1] in {"Low_confidence"}
        else:
            assert len(parts) == 5 and parts[2].startswith("ST_")
            assert parts[0] == safe_component(sample["metadata"]["organism"]["genus"], "Unknown_genus")
    assert len(identifiers) == 3


def ani_panel(root, genomes):
    """A panel in the shape identify_species already accepts, built from real bytes."""
    root.mkdir(parents=True, exist_ok=True)
    manifest = {"format_version": 1, "source_revision": "synthetic-truth",
                "source_repository": "synthetic-test-fixture", "species": [], "virulence": {},
                "files": [], "limitations": ["One reference per taxon is a triage panel."]}
    for identifier, genus, species, text in genomes:
        (root / f"{identifier}.fasta").write_text(f">{identifier}\n{text}\n")
        manifest["species"].append({"id": identifier, "path": f"{identifier}.fasta", "genus": genus,
                                    "species": species, "subspecies": "", "outgroup": True})
    for path in sorted(root.iterdir()):
        manifest["files"].append({"path": path.name, "bytes": path.stat().st_size,
                                  "sha256": file_sha256(path)})
    manifest["reference_digest"] = reference_digest(manifest)
    (root / "manifest.json").write_text(json.dumps(manifest))
    return root


def test_a_mixed_drop_is_identified_on_the_originals_then_filed_once(project, tmp_path):
    genomes = {seed: "".join(random.Random(seed).choices("ACGT", k=200_000)) for seed in (1, 2, 3)}
    panel = ani_panel(tmp_path / "panel", [("kp", "Klebsiella", "pneumoniae", genomes[1]),
                                           ("pa", "Pseudomonas", "aeruginosa", genomes[2])])
    inbox = tmp_path / "inbox"
    inbox.mkdir()
    inputs = [sequence(inbox / "first.fasta", genomes[1]),
              sequence(inbox / "second.fasta", genomes[2]),
              sequence(inbox / "stranger.fasta", genomes[3])]
    reads = inbox / "sample_R1.fastq"
    reads.write_text("@read1\nACGTACGT\n+\nIIIIIIII\n")
    inputs.append(reads)
    before = {path: path.read_bytes() for path in inputs}

    verdicts = identify_batch(inputs, species_panel_root=panel)
    assert [Path(verdict["input_path"]) for verdict in verdicts] == inputs
    # Identification ran on the originals; nothing has been copied yet.
    assert {path: path.read_bytes() for path in inputs} == before
    assert not (tmp_path / "managed").exists()
    assert project.samples() == []

    # The default policy reviews everything, so no proposal files itself.
    assert all(verdict["status"] != "confirmed" for verdict in verdicts)
    accepted = []
    for verdict in verdicts:
        if verdict["confidence"] == "genomic_reference_supported":
            verdict = dict(verdict, status="confirmed", confirmed_by="user",
                           accepted=dict(verdict["proposed"]))
        accepted.append(assignment_for(verdict))

    root = tmp_path / "managed"
    notes = []
    identifiers = import_samples(project, accepted, root, duplicates="skip", notes=notes)
    assert len(identifiers) == 4 and notes == []
    filed = {Path(sample["input_path"]).relative_to(root).parts[:2]
             for sample in project.samples()}
    assert ("Klebsiella", "pneumoniae") in filed and ("Pseudomonas", "aeruginosa") in filed
    assert (QUARANTINE_ROOT, "Not_in_reference_panel") in filed
    assert (QUARANTINE_ROOT, "Reads_not_assembled") in filed
    assert {path: path.read_bytes() for path in inputs} == before
    for sample in project.samples():
        evidence = sample["metadata"]["organism_evidence"]
        assert evidence["confirmed_by"] != "auto_policy"
        assert sample["id"] in Path(sample["input_path"]).parts
        if evidence["status"] == "quarantined":
            assert sample["metadata"]["organism"] == {"genus": "", "species": ""}

    # A quarantined isolate is a first-class sample that re-files on a later decision.
    stranger = next(sample for sample in project.samples()
                    if sample["metadata"]["organism_evidence"].get("quarantine_reason")
                    == "not_in_reference_panel")
    reassign_organism(project, [stranger["id"]], "Serratia", "marcescens", storage_root=root)
    moved = Path(project.get_sample(stranger["id"])["input_path"]).relative_to(root)
    assert moved.parts[:3] == ("Serratia", "marcescens", "ST_unassigned")
    assert not (root / QUARANTINE_ROOT / "Not_in_reference_panel").exists()
