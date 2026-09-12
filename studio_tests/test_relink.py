import hashlib
from pathlib import Path

import pytest

from wmlstudio.project import Project
from wmlstudio.sequence import AnalysisCancelled
from wmlstudio.storage import relink_input


def sample(project, tmp_path):
    original = tmp_path / "old" / "genome.fasta"
    original.parent.mkdir()
    original.write_text(">genome\nACGTTGCA\n")
    sid = project.add_sample(original)
    digest = hashlib.sha256(original.read_bytes()).hexdigest()
    project.set_result(sid, {"scheme": "MLST", "scheme_digest": "scheme",
                            "st": "20", "status": "complete", "alleles": {"arcC": "4"},
                            "input_sha256": digest})
    return sid, original, digest


def test_relink_moved_copy_keeps_profile_status_and_other_evidence(tmp_path):
    with Project(tmp_path / "study.wmlstudio") as project:
        sid, original, digest = sample(project, tmp_path)
        project.set_analysis(sid, {"scheme": "cgMLST", "scheme_digest": "secondary",
                                  "alleles": {"locus": "1"}, "input_sha256": digest})
        project.set_metadata(sid, {"hydra": {"report_sha256": "retained"},
                                  "workflow": {"managed": True, "storage_root": str(original.parent),
                                               "managed_sha256": digest, "source_path": str(original)}})
        before = project.get_sample(sid)
        analyses = project.analysis_results(sid)
        new = tmp_path / "moved.fasta"
        new.write_bytes(original.read_bytes())
        original.unlink()
        assert relink_input(project, sid, new) == new
        after = project.get_sample(sid)
        assert not after["missing_input"] and after["input_path"] == str(new)
        assert after["result"] == before["result"] and after["status"] == before["status"]
        assert project.analysis_results(sid) == analyses
        assert after["metadata"]["hydra"] == before["metadata"]["hydra"]
        assert after["metadata"]["workflow"]["managed"] is False
        assert after["metadata"]["workflow"]["storage_root"] is None
        assert any(row["action"] == "input_relinked" for row in project.history(sid))


def test_different_bytes_or_missing_hash_never_reassign_evidence(tmp_path):
    with Project(tmp_path / "study.wmlstudio") as project:
        sid, original, digest = sample(project, tmp_path)
        wrong = tmp_path / "wrong.fasta"
        wrong.write_text(">genome\nAAAAAAAA\n")
        before = project.get_sample(sid)
        with pytest.raises(ValueError, match="different SHA-256"):
            relink_input(project, sid, wrong)
        assert project.get_sample(sid) == before
        untyped = project.add_sample(wrong)
        with pytest.raises(ValueError, match="No complete input"):
            relink_input(project, untyped, original)


def test_assembly_hash_takes_priority_over_old_raw_source_hash(tmp_path):
    with Project(tmp_path / "study.wmlstudio") as project:
        sid, original, digest = sample(project, tmp_path)
        project.invalidate_result(sid, "Assembly awaits typing")
        project.set_metadata(sid, {"assembly": {"provenance": {"assembly_sha256": digest}},
                                  "workflow": {"source_kind": "assembly", "source_sha256": "1" * 64}})
        new = tmp_path / "assembly.fasta"
        new.write_bytes(original.read_bytes())
        relink_input(project, sid, new)
        result = project.get_sample(sid)
        assert result["status"] == "queued" and result["result"] is None
        assert result["metadata"]["assembly"]["assembly_path"] == str(new)
        assert result["metadata"]["assembly"]["provenance"]["assembly_sha256"] == digest


def test_untyped_managed_sample_uses_stored_full_input_hash(tmp_path):
    with Project(tmp_path / "study.wmlstudio") as project:
        path = tmp_path / "untyped.fasta"
        path.write_text(">a\nACGT\n")
        sid = project.add_sample(path)
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        root = tmp_path / "Sequences"
        destination = root / "Genus" / "species" / "ST_unassigned" / sid / path.name
        destination.parent.mkdir(parents=True)
        destination.write_bytes(path.read_bytes())
        project.set_metadata(sid, {"workflow": {"managed": True, "storage_root": str(root),
                                               "source_path": str(path), "managed_sha256": digest}})
        relink_input(project, sid, destination)
        assert project.get_sample(sid)["metadata"]["workflow"]["managed"] is True


def test_cancel_and_running_sample_do_not_change_path(tmp_path):
    with Project(tmp_path / "study.wmlstudio") as project:
        sid, original, _ = sample(project, tmp_path)
        with pytest.raises(AnalysisCancelled):
            relink_input(project, sid, original, cancelled=lambda: True)
        project.set_status(sid, "running")
        with pytest.raises(ValueError, match="Wait"):
            relink_input(project, sid, original)


def test_changed_replacement_file_is_detected_before_database_write(tmp_path, monkeypatch):
    import wmlstudio.storage as storage

    with Project(tmp_path / "study.wmlstudio") as project:
        sid, original, digest = sample(project, tmp_path)
        replacement = tmp_path / "new.fasta"
        replacement.write_bytes(original.read_bytes())
        def changed_hash(path, cancelled):
            Path(path).write_text(">changed\nTTTT\n")
            return digest
        monkeypatch.setattr(storage, "file_sha256", changed_hash)
        with pytest.raises(ValueError, match="changed during hashing"):
            relink_input(project, sid, replacement)
        assert project.get_sample(sid)["input_path"] == str(original)
