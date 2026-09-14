import copy
import gzip
from pathlib import Path

import pytest

from wmlstudio.project import Project
from wmlstudio.read_attachment_dialog import AttachReadsDialog
from wmlstudio.read_attachments import (
    attach_read_pair,
    suggest_read_attachments,
    validate_read_attachment,
)
from wmlstudio.sequence import AnalysisCancelled, SequenceError, file_sha256


def fixture(tmp_path, marked=True):
    assembly = tmp_path / "isolate_contigs.fasta"
    assembly.write_text(">genome\nACGTTGCAACGT\n")
    paths = [tmp_path / f"isolate_R{mate}_001.fastq.gz" for mate in (1, 2)]
    for mate, path in enumerate(paths, 1):
        text = "".join(f"@pair{index}{'/' + str(mate) if marked else ''}\nACGT\n+\nIIII\n" for index in range(3))
        path.write_bytes(gzip.compress(text.encode()))
    sample = {"id": "stable-assembly", "name": "isolate", "input_path": str(assembly), "metadata": {}}
    reads = [{"id": f"read-{i}", "name": path.name, "input_path": str(path)} for i, path in enumerate(paths, 1)]
    return sample, reads


def test_exact_filename_suggestions_withhold_ambiguous_and_shared_read_matches(tmp_path):
    sample, reads = fixture(tmp_path)
    result = suggest_read_attachments([sample], reads)
    assert result["suggestions"][0]["read1_id"] == "read-1"
    extra = dict(sample, id="alternate-assembly")
    result = suggest_read_attachments([sample, extra], reads)
    assert not result["suggestions"]
    assert set(result["ambiguous"]) == {sample["id"], extra["id"]}
    duplicate_read = dict(reads[0], id="duplicate-r1")
    assert not suggest_read_attachments([sample], [*reads, duplicate_read])["suggestions"]
    unrelated = dict(sample, id="unrelated", input_path=str(tmp_path / "isolate2.fasta"))
    assert suggest_read_attachments([unrelated], reads)["missing"] == ["unrelated"]


def test_full_pair_qc_hashes_are_validated_without_altering_any_original(tmp_path):
    sample, reads = fixture(tmp_path)
    paths = [Path(sample["input_path"]), *[Path(read["input_path"]) for read in reads]]
    original = [path.read_bytes() for path in paths]
    evidence = validate_read_attachment(sample, *paths[1:])
    assert [path.read_bytes() for path in paths] == original
    assert evidence["pairing"] == {"records_checked": 3, "complete_file": True, "sampled": False,
                                    "explicit_mates": True, "status": "verified"}
    assert [read["qc"]["q30_percent"] for read in evidence["reads"]] == [100, 100]
    assert [read["sha256"] for read in evidence["reads"]] == [file_sha256(path) for path in paths[1:]]
    assert evidence["attachment_id"] == validate_read_attachment(sample, *paths[1:])["attachment_id"]
    assert "not FastQC or fastp" in evidence["preprocessing"]
    assert "do not prove" in evidence["identity_basis"]


def test_unmarked_pairs_are_user_associated_not_biologically_verified(tmp_path):
    sample, reads = fixture(tmp_path, marked=False)
    evidence = validate_read_attachment(sample, *[read["input_path"] for read in reads])
    assert evidence["pairing"]["status"] == "unmarked"
    assert not evidence["pairing"]["explicit_mates"]


def test_last_pair_mismatch_and_cancel_are_not_silently_accepted(tmp_path):
    sample, reads = fixture(tmp_path)
    paths = [Path(read["input_path"]) for read in reads]
    second = gzip.decompress(paths[1].read_bytes()).replace(b"@pair2/2", b"@different/2")
    paths[1].write_bytes(gzip.compress(second))
    with pytest.raises(SequenceError, match="pair 3"):
        validate_read_attachment(sample, *paths)
    with pytest.raises(AnalysisCancelled):
        validate_read_attachment(sample, *paths, cancelled=lambda: True)


def test_attachment_preserves_stable_assembly_and_prior_evidence_and_links_source_records(tmp_path):
    sample, reads = fixture(tmp_path)
    with Project(tmp_path / "library.wmlstudio") as project:
        identifier = project.add_sample(sample["input_path"], "isolate", sample_id=sample["id"])
        for read in reads:
            project.add_sample(read["input_path"], read["name"], sample_id=read["id"])
        project.set_result(identifier, {"status": "complete", "st": "10", "alleles": {"abc": "1"}})
        project.update_metadata(identifier, {"hydra": {"report_sha256": "retained"}})
        before = project.get_sample(identifier)
        evidence = validate_read_attachment(before, *[read["input_path"] for read in reads],
                                            read_sample_ids=[read["id"] for read in reads])
        key = attach_read_pair(project, identifier, evidence)
        assert key == attach_read_pair(project, identifier, evidence)
        after = project.get_sample(identifier)
        assert after["id"] == before["id"] and after["input_path"] == before["input_path"]
        assert after["result"] == before["result"] and after["status"] == before["status"]
        assert after["metadata"]["hydra"] == before["metadata"]["hydra"]
        assert len(project.samples()) == 3
        for read in reads:
            retained = project.get_sample(read["id"])
            assert retained["input_path"] == read["input_path"]
            assert retained["metadata"]["workflow"]["source_kind"] == "read_mate"
            assert retained["metadata"]["workflow"]["paired_with"] == identifier
        assert sum(item["action"] == "reads_attached" for item in project.history(identifier)) == 1


def test_attachment_rejects_changed_input_and_wrong_stable_id_without_mutation(tmp_path):
    sample, reads = fixture(tmp_path)
    with Project(tmp_path / "library.wmlstudio") as project:
        sid = project.add_sample(sample["input_path"], "isolate", sample_id=sample["id"])
        evidence = validate_read_attachment(sample, *[read["input_path"] for read in reads])
        forged = copy.deepcopy(evidence)
        forged["sample_id"] = "another"
        with pytest.raises(ValueError, match="different stable"):
            attach_read_pair(project, sid, forged)
        Path(reads[0]["input_path"]).write_bytes(b"changed")
        with pytest.raises(ValueError, match="changed"):
            attach_read_pair(project, sid, evidence)
        assert "reads" not in project.get_sample(sid)["metadata"]


def test_native_dialog_suggests_and_requires_two_distinct_mates(qtbot, tmp_path):
    sample, reads = fixture(tmp_path)
    dialog = AttachReadsDialog([sample], reads)
    qtbot.addWidget(dialog)
    assert dialog.rows[0][0].currentData()["sample_id"] == "read-1"
    assert dialog.rows[0][1].currentData()["sample_id"] == "read-2"
    dialog.rows[0][1].setCurrentIndex(1)
    dialog.accept()
    assert not dialog.assignments
    assert "only one isolate" in dialog.feedback.text()
    dialog.rows[0][1].setCurrentIndex(2)
    dialog.accept()
    assert dialog.assignments[0]["sample_id"] == sample["id"]
    assert dialog.assignments[0]["read_sample_ids"] == ["read-1", "read-2"]
