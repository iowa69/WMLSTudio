import sqlite3

import pytest

from wmlstudio.project import APPLICATION_ID, Project, ProjectError


@pytest.fixture
def sequence(tmp_path):
    path = tmp_path / "sample.fasta"
    path.write_text(">contig\nACGTACGT\n")
    return path


def test_project_roundtrip_and_missing_source(tmp_path, sequence):
    path = tmp_path / "study.wmlstudio"
    with Project(path) as project:
        sample_id = project.add_sample(sequence, "Isolate α")
        project.set_metadata(sample_id, {"site": "Dublin", "year": 2026})
        project.set_result(sample_id, {"sample_name": "Isolate α", "alleles": {"abc": "7"}})
        project.set_setting("scheme", {"path": "/example/scheme"})
        before = project.get_sample(sample_id)
    sequence.unlink()
    with Project(path) as project:
        sample = project.samples()[0]
        assert sample["id"] == sample_id
        assert sample["status"] == "completed"
        assert sample["result"]["sample_id"] == sample_id
        assert sample["result"]["alleles"] == {"abc": "7"}
        assert sample["metadata"] == {"site": "Dublin", "year": 2026}
        assert sample["created_at"] == before["created_at"]
        assert sample["missing_input"] is True
        assert project.get_setting("scheme") == {"path": "/example/scheme"}
        assert project.get_setting("unconfigured", 9) == 9


def test_reopening_recovers_running_only(tmp_path, sequence):
    path = tmp_path / "study.wmlstudio"
    with Project(path) as project:
        queued = project.add_sample(sequence)
        running = project.add_sample(sequence)
        failed = project.add_sample(sequence)
        project.set_status(running, "running")
        project.set_status(failed, "failed", "Invalid FASTQ")
    with Project(path) as project:
        assert project.get_sample(queued)["status"] == "queued"
        assert project.get_sample(running)["status"] == "interrupted"
        assert "Rerun" in project.get_sample(running)["error"]
        assert project.get_sample(failed)["error"] == "Invalid FASTQ"


def test_transaction_rolls_back_multiple_writes(tmp_path, sequence):
    with Project(tmp_path / "study.wmlstudio") as project:
        sample_id = project.add_sample(sequence)
        with pytest.raises(RuntimeError, match="simulate"):
            with project.transaction():
                project.set_metadata(sample_id, {"site": "new"})
                project.add_sample(sequence)
                project.set_setting("threshold", 0.9)
                raise RuntimeError("simulate interrupted import")
        assert len(project.samples()) == 1
        assert project.get_sample(sample_id)["metadata"] == {}
        assert project.get_setting("threshold") is None


def test_invalid_data_does_not_destroy_result(tmp_path, sequence):
    with Project(tmp_path / "study.wmlstudio") as project:
        sample_id = project.add_sample(sequence)
        project.set_result(sample_id, {"st": "1"})
        with pytest.raises(ValueError):
            project.set_result(sample_id, {"st": float("nan")})
        assert project.get_sample(sample_id)["result"]["st"] == "1"
        with pytest.raises(ValueError):
            project.set_status(sample_id, "success-ish")
        with pytest.raises(KeyError):
            project.set_status("does-not-exist", "failed")
        with pytest.raises(FileNotFoundError):
            project.add_sample(tmp_path / "missing.fasta")
        with pytest.raises(ValueError):
            project.add_sample(sequence, " ")


def test_deleting_record_keeps_input(tmp_path, sequence):
    with Project(tmp_path / "study.wmlstudio") as project:
        sample_id = project.add_sample(sequence)
        project.remove_sample(sample_id)
        assert project.samples() == []
        assert sequence.is_file()
        with pytest.raises(KeyError):
            project.remove_sample(sample_id)


def test_rejects_foreign_database_without_modifying_it(tmp_path):
    path = tmp_path / "foreign.sqlite"
    with sqlite3.connect(path) as connection:
        connection.execute("CREATE TABLE valuable (value TEXT)")
        connection.execute("INSERT INTO valuable VALUES ('keep')")
    before = path.read_bytes()
    with pytest.raises(ProjectError, match="not a WMLSTudio"):
        Project(path)
    assert path.read_bytes() == before


def test_rejects_future_schema_without_recovering_jobs(tmp_path, sequence):
    path = tmp_path / "future.wmlstudio"
    with Project(path) as project:
        project.set_status(project.add_sample(sequence), "running")
    with sqlite3.connect(path) as connection:
        connection.execute("PRAGMA user_version = 999")
    before = path.read_bytes()
    with pytest.raises(ProjectError, match="version 999"):
        Project(path)
    assert path.read_bytes() == before


def test_rejects_corrupt_file_and_invalid_structure(tmp_path):
    corrupt = tmp_path / "corrupt.wmlstudio"
    corrupt.write_bytes(b"not sqlite")
    with pytest.raises(ProjectError, match="not a readable SQLite"):
        Project(corrupt)
    malformed = tmp_path / "malformed.wmlstudio"
    with sqlite3.connect(malformed) as connection:
        connection.execute(f"PRAGMA application_id = {APPLICATION_ID}")
        connection.execute("PRAGMA user_version = 1")
        connection.execute("CREATE TABLE samples (id TEXT)")
    with pytest.raises(ProjectError, match="missing required columns"):
        Project(malformed)
