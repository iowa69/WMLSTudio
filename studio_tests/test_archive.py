import pytest

from wmlstudio.archive import (
    active_samples,
    archive_record,
    archive_samples,
    archived_samples,
    is_archived,
    restore_samples,
)
from wmlstudio.project import Project


@pytest.fixture
def sequence(tmp_path):
    path = tmp_path / "isolate.fasta"
    path.write_text(">contig\nACGTACGT\n")
    return path


@pytest.fixture
def project(tmp_path):
    with Project(tmp_path / "study.wmlstudio") as project:
        yield project


def typed_sample(project, sequence, name="Isolate A"):
    sample_id = project.add_sample(sequence, name)
    project.set_metadata(sample_id, {"organism": {"genus": "Klebsiella", "species": "pneumoniae"},
                                     "annotations": {"ward": "ICU"}})
    project.set_result(sample_id, {"scheme": "MLST", "scheme_digest": "mlst-v1", "st": "258",
                                   "alleles": {"gapA": "3"}, "input_sha256": "a" * 64})
    project.set_analysis(sample_id, {"scheme": "cgMLST", "scheme_digest": "cgmlst-v1",
                                     "alleles": {"locus0001": "7"}})
    return sample_id


def test_archiving_hides_an_isolate_from_the_working_set_but_destroys_no_evidence(project, sequence):
    sample_id = typed_sample(project, sequence)
    before = project.get_sample(sample_id)
    analyses = project.analysis_results(sample_id)
    payload = sequence.read_bytes()

    assert archive_samples(project, [sample_id], reason="Wrong organism assigned at import") == [sample_id]

    sample = project.get_sample(sample_id)
    assert is_archived(sample) is True
    assert active_samples(project.samples()) == []
    assert [s["id"] for s in archived_samples(project.samples())] == [sample_id]
    assert len(project.samples()) == 1
    assert sample["result"] == before["result"]
    assert sample["status"] == before["status"]
    assert sample["metadata"]["organism"] == before["metadata"]["organism"]
    assert sample["metadata"]["annotations"] == before["metadata"]["annotations"]
    assert project.analysis_results(sample_id) == analyses
    assert sequence.is_file() and sequence.read_bytes() == payload


def test_archive_history_records_the_reason_and_the_time_it_was_hidden(project, sequence):
    sample_id = typed_sample(project, sequence)
    archive_samples(project, [sample_id], reason="  Contaminated culture  ")

    events = [event for event in project.history(sample_id) if event["action"] == "sample_archived"]
    assert len(events) == 1
    assert events[0]["details"]["reason"] == "Contaminated culture"
    assert events[0]["details"]["at"]
    assert archive_record(project.get_sample(sample_id))["reason"] == "Contaminated culture"


def test_restoring_removes_the_archive_note_instead_of_marking_it_false(project, sequence):
    sample_id = typed_sample(project, sequence)
    before = project.get_sample(sample_id)
    analyses = project.analysis_results(sample_id)
    archive_samples(project, [sample_id], reason="Filed by mistake")

    assert restore_samples(project, [sample_id]) == [sample_id]

    sample = project.get_sample(sample_id)
    assert "archived" not in sample["metadata"]
    assert is_archived(sample) is False
    assert sample["metadata"] == before["metadata"]
    assert sample["result"] == before["result"]
    assert sample["name"] == before["name"] and sample["input_path"] == before["input_path"]
    assert project.analysis_results(sample_id) == analyses
    assert [s["id"] for s in active_samples(project.samples())] == [sample_id]
    assert any(event["action"] == "sample_restored" for event in project.history(sample_id))


def test_archiving_and_restoring_are_idempotent_and_keep_the_first_archive_note(project, sequence):
    sample_id = typed_sample(project, sequence)
    archive_samples(project, [sample_id], reason="First reason")
    first = archive_record(project.get_sample(sample_id))

    assert archive_samples(project, [sample_id], reason="Second reason") == []
    assert archive_record(project.get_sample(sample_id)) == first

    restore_samples(project, [sample_id])
    assert restore_samples(project, [sample_id]) == []
    assert len([e for e in project.history(sample_id) if e["action"] == "sample_restored"]) == 1


def test_archiving_a_running_analysis_is_refused_and_writes_nothing(project, sequence):
    sample_id = typed_sample(project, sequence)
    project.set_status(sample_id, "running")
    before = project.get_sample(sample_id)["metadata"]

    with pytest.raises(ValueError, match="Wait for this sample's analysis"):
        archive_samples(project, [sample_id])

    assert project.get_sample(sample_id)["metadata"] == before
    assert not [e for e in project.history(sample_id) if e["action"] == "sample_archived"]


def test_an_unknown_identifier_rolls_back_the_whole_archive_batch(project, sequence):
    first = typed_sample(project, sequence, "Isolate A")
    second = typed_sample(project, sequence, "Isolate B")

    with pytest.raises(KeyError):
        archive_samples(project, [first, "does-not-exist", second])

    assert active_samples(project.samples()) == project.samples()
    assert not [e for e in project.history() if e["action"] == "sample_archived"]


def test_a_stray_archived_annotation_of_another_shape_never_hides_an_isolate(project, sequence):
    sample_id = project.add_sample(sequence)
    project.set_metadata(sample_id, {"archived": "no"})

    assert is_archived(project.get_sample(sample_id)) is False
    assert archive_record(project.get_sample(sample_id)) is None
    assert [s["id"] for s in active_samples(project.samples())] == [sample_id]


def test_filters_accept_plain_records_without_touching_a_project():
    records = [{"id": "a", "metadata": {}}, {"id": "b", "metadata": {"archived": {"at": "now"}}},
               {"id": "c"}]
    assert [record["id"] for record in active_samples(records)] == ["a", "c"]
    assert [record["id"] for record in archived_samples(records)] == ["b"]
    assert active_samples(records)[0] is records[0]
