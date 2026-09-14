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


def test_renaming_a_sample_keeps_its_input_result_and_identifier(tmp_path, sequence):
    with Project(tmp_path / "study.wmlstudio") as project:
        sample_id = project.add_sample(sequence, "Isolate α")
        project.set_result(sample_id, {"scheme": "MLST", "scheme_digest": "mlst-v1", "st": "258"})
        before = project.get_sample(sample_id)
        project.rename_sample(sample_id, "  Ward B isolate 7  ")
        sample = project.get_sample(sample_id)
        assert sample["name"] == "Ward B isolate 7"
        assert sample["id"] == sample_id
        assert sample["input_path"] == before["input_path"]
        assert sample["result"] == before["result"]
        assert sequence.is_file()
        renames = [event for event in project.history(sample_id) if event["action"] == "sample_renamed"]
        assert renames[0]["details"] == {"from": "Isolate α", "to": "Ward B isolate 7"}
        project.rename_sample(sample_id, "Ward B isolate 7")
        assert len([e for e in project.history(sample_id) if e["action"] == "sample_renamed"]) == 1


def test_renaming_refuses_a_blank_name_or_an_unknown_sample(tmp_path, sequence):
    with Project(tmp_path / "study.wmlstudio") as project:
        sample_id = project.add_sample(sequence, "Isolate α")
        with pytest.raises(ValueError, match="cannot be empty"):
            project.rename_sample(sample_id, "   ")
        assert project.get_sample(sample_id)["name"] == "Isolate α"
        with pytest.raises(KeyError):
            project.rename_sample("does-not-exist", "Anything")


def test_removal_records_the_analyses_it_deletes_so_they_can_be_restored(tmp_path, sequence):
    with Project(tmp_path / "study.wmlstudio") as project:
        sample_id = project.add_sample(sequence, "Isolate α")
        project.set_metadata(sample_id, {"organism": {"genus": "Klebsiella", "species": "pneumoniae"}})
        project.set_result(sample_id, {"scheme": "MLST", "scheme_digest": "mlst-v1", "st": "258",
                                       "alleles": {"gapA": "3"}})
        project.set_analysis(sample_id, {"scheme": "cgMLST", "scheme_digest": "cgmlst-v1",
                                         "alleles": {"locus0001": "7"}})
        before = project.get_sample(sample_id)
        analyses = project.analysis_results(sample_id)
        project.remove_sample(sample_id)

        removal = [event for event in project.history(sample_id) if event["action"] == "sample_removed"][0]
        assert removal["details"]["format_version"] == 2
        assert removal["details"]["sample"] == before
        assert removal["details"]["analyses"] == analyses
        assert sequence.is_file()

        assert project.restore_removed_sample(removal["id"]) == sample_id
        restored = project.get_sample(sample_id)
        assert restored["name"] == before["name"]
        assert restored["input_path"] == before["input_path"]
        assert restored["status"] == before["status"] and restored["error"] == before["error"]
        assert restored["metadata"] == before["metadata"]
        assert restored["result"] == before["result"]
        assert restored["created_at"] == before["created_at"]
        assert project.analysis_results(sample_id) == analyses
        assert any(event["action"] == "sample_restored_from_history"
                   for event in project.history(sample_id))


def test_restoring_refuses_a_live_identifier_or_an_unrelated_history_entry(tmp_path, sequence):
    with Project(tmp_path / "study.wmlstudio") as project:
        sample_id = project.add_sample(sequence, "Isolate α")
        project.remove_sample(sample_id)
        removal = [event for event in project.history(sample_id) if event["action"] == "sample_removed"][0]
        project.restore_removed_sample(removal["id"])

        with pytest.raises(ValueError, match="already in this project"):
            project.restore_removed_sample(removal["id"])
        imported = [event for event in project.history(sample_id) if event["action"] == "sample_imported"][0]
        with pytest.raises(ValueError, match="not hold a restorable"):
            project.restore_removed_sample(imported["id"])
        with pytest.raises(KeyError):
            project.restore_removed_sample(99999)
        assert len(project.samples()) == 1


def test_restoring_a_sample_removed_mid_analysis_never_claims_it_is_still_running(tmp_path, sequence):
    with Project(tmp_path / "study.wmlstudio") as project:
        sample_id = project.add_sample(sequence)
        project.set_status(sample_id, "running")
        project.remove_sample(sample_id)
        removal = [event for event in project.history(sample_id) if event["action"] == "sample_removed"][0]
        project.restore_removed_sample(removal["id"])
        restored = project.get_sample(sample_id)
        assert restored["status"] == "interrupted"
        assert "Rerun" in restored["error"]


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


def test_version_one_project_migrates_without_losing_original_mlst(tmp_path, sequence):
    path = tmp_path / "legacy.wmlstudio"
    with Project(path) as project:
        sid = project.add_sample(sequence)
        project.set_metadata(sid, {"ward": "A"})
        project.set_result(sid, {"scheme": "MLST", "scheme_digest": "mlst-v1", "st": "131", "alleles": {"adk": "1"}})
        project.set_setting("favorite", "retained")
    with sqlite3.connect(path) as connection:
        connection.execute("DROP TABLE history")
        connection.execute("DROP TABLE analyses")
        connection.execute("PRAGMA user_version = 1")
    with Project(path) as project:
        assert project.get_sample(sid)["metadata"] == {"ward": "A"}
        assert project.get_setting("favorite") == "retained"
        assert project.analysis_results(sid)[0]["st"] == "131"
        project.set_analysis(sid, {"scheme": "cgMLST", "scheme_digest": "cgmlst-v1", "st": None,
                                   "alleles": {"locus0001": "7"}})
        assert project.get_sample(sid)["result"]["st"] == "131"
        assert len(project.analysis_results(sid)) == 2
        assert any(event["action"] == "schema_migrated" for event in project.history())
    with Project(path) as project:
        assert len(project.analysis_results(sid)) == 2
        assert project.get_sample(sid)["result"]["st"] == "131"


def test_profile_only_samples_need_no_fasta_and_keep_stable_ids(tmp_path):
    path = tmp_path / "profiles.wmlstudio"
    with Project(path) as project:
        sid = project.add_profile("External A", {"scheme": "cgMLST", "scheme_digest": "external",
                                                "kind": "profile", "status": "profile_imported", "alleles": {"l1": "7"}})
        assert project.get_sample(sid)["profile_only"] is True
        assert project.get_sample(sid)["missing_input"] is False
        assert project.get_sample(sid)["input_path"] == ""
    with Project(path) as project:
        assert project.analysis_results(sid)[0]["alleles"] == {"l1": "7"}


def test_collections_persist_and_do_not_delete_samples(tmp_path, sequence):
    path = tmp_path / "study.wmlstudio"
    with Project(path) as project:
        ids = [project.add_sample(sequence) for _ in range(2)]
        collection = project.create_collection("Ward A")
        assert project.create_collection("ward a") == collection
        project.set_collection_members(collection, ids)
        project.set_collection_members(collection, [ids[0]], add=False)
        project.rename_collection(collection, "Investigation 2026")
        assert project.collections()[0]["sample_ids"] == [ids[1]]
        with pytest.raises(KeyError):
            project.set_collection_members(collection, ["missing"])
    with Project(path) as project:
        assert project.collections()[0]["name"] == "Investigation 2026"
        project.delete_collection(collection)
        assert len(project.samples()) == 2
        assert project.collections() == []


MLST = {"scheme": "Klebsiella MLST", "scheme_digest": "mlst-v1", "st": "258", "status": "complete",
        "alleles": {f"gene{index}": "1" for index in range(7)}}
CGMLST = {"scheme": "Klebsiella cgMLST", "scheme_digest": "cg-v1", "st": None, "status": "complete",
          "alleles": {f"locus{index:04d}": "3" for index in range(400)},
          "parameters": {"method": "full-cds-cgmlst-v2"}}


def test_a_stored_profile_says_which_typing_it_is_and_on_what_basis():
    from wmlstudio.project import classify_typing
    assert classify_typing(MLST) == {"kind": "mlst", "basis": "7 loci, at classical scheme size"}
    assert classify_typing(CGMLST)["kind"] == "cgmlst"
    assert "full-cds-cgmlst-v2" in classify_typing(CGMLST)["basis"]
    # A caller that already knows is believed, and a declared scheme type outranks size.
    assert classify_typing({**MLST, "analysis_kind": "cgmlst"})["kind"] == "cgmlst"
    small = {"scheme_digest": "x", "alleles": {"a": "1"}, "scheme_metadata": {"type": "cgMLST"}}
    assert classify_typing(small)["kind"] == "cgmlst"
    assert "declares type" in classify_typing(small)["basis"]
    # Quality-only work is not a profile of either kind and is never sorted as one.
    quality = {"status": "qc_only", "st": None, "alleles": {}, "scheme_digest": None}
    assert classify_typing(quality) == {"kind": "unclassified",
                                        "basis": "no allelic profile is stored"}


def test_classical_st_and_core_genome_profiles_are_retrievable_independently(tmp_path, sequence):
    path = tmp_path / "study.wmlstudio"
    with Project(path) as project:
        sample_id = project.add_sample(sequence, "Isolate 1")
        project.set_result(sample_id, MLST)
        project.set_analysis(sample_id, CGMLST)

        assert project.latest_analysis(sample_id, "mlst")["st"] == "258"
        assert project.latest_analysis(sample_id, "cgmlst")["scheme_digest"] == "cg-v1"
        assert len(project.latest_analysis(sample_id, "cgmlst")["alleles"]) == 400
        grouped = project.analyses_by_kind(sample_id)
        assert [result["scheme_digest"] for result in grouped["mlst"]] == ["mlst-v1"]
        assert [result["scheme_digest"] for result in grouped["cgmlst"]] == ["cg-v1"]
        assert grouped["unclassified"] == []
        assert project.typing_kind_index()[sample_id] == {"mlst": 1, "cgmlst": 1,
                                                          "unclassified": 0}
        kinds = {row["scheme_digest"]: row["typing_kind"]
                 for row in project.analysis_summaries(sample_id)}
        assert kinds == {"mlst-v1": "mlst", "cg-v1": "cgmlst"}
        # The stored scientific payload is exactly what was handed in; the kind is
        # recorded beside it, never injected into the result.
        assert "typing_kind" not in project.latest_analysis(sample_id, "mlst")

    with Project(path) as project:
        assert project.latest_analysis(sample_id, "mlst")["st"] == "258"
        assert len(project.latest_analysis(sample_id, "cgmlst")["alleles"]) == 400


def test_a_core_genome_run_does_not_take_over_the_classical_st(tmp_path, sequence):
    with Project(tmp_path / "study.wmlstudio") as project:
        sample_id = project.add_sample(sequence)
        project.set_result(sample_id, MLST)
        project.set_result(sample_id, CGMLST)
        # The headline result is whatever ran last, but the ST is still the ST.
        assert project.get_sample(sample_id)["result"]["scheme_digest"] == "cg-v1"
        assert project.latest_analysis(sample_id, "mlst")["st"] == "258"
        assert project.latest_analysis(sample_id, "cgmlst")["st"] is None
        event = next(e for e in project.history(sample_id)
                     if e["action"] == "analysis_completed" and e["details"]["scheme_digest"] == "cg-v1")
        assert event["details"]["typing_kind"] == "cgmlst"


def test_each_typing_kind_can_be_collected_across_the_project_without_the_other(tmp_path, sequence):
    with Project(tmp_path / "study.wmlstudio") as project:
        first = project.add_sample(sequence, "One")
        second = project.add_sample(sequence, "Two")
        project.set_result(first, MLST)
        project.set_analysis(first, CGMLST)
        project.set_result(second, MLST)

        classical = project.results_of_kind("mlst")
        assert sorted(row["sample_name"] for row in classical) == ["One", "Two"]
        assert {row["st"] for row in classical} == {"258"}
        core = project.results_of_kind("cgmlst")
        assert [row["sample_id"] for row in core] == [first]
        assert core[0]["job_status"] == "completed"
        # A sample with no core-genome profile is absent from that view, never
        # represented there by its ST.
        assert second not in {row["sample_id"] for row in core}
        assert project.latest_analysis(second, "cgmlst") is None
        with pytest.raises(ValueError, match="Unknown typing kind"):
            project.results_of_kind("snp")


def test_a_schema_three_project_gains_typing_kinds_without_rewriting_its_results(tmp_path, sequence):
    path = tmp_path / "legacy.wmlstudio"
    with Project(path) as project:
        sample_id = project.add_sample(sequence)
        project.set_result(sample_id, MLST)
        project.set_analysis(sample_id, CGMLST)
        stored = project.analysis_results(sample_id)
    with sqlite3.connect(path) as connection:
        connection.execute("CREATE TABLE legacy AS SELECT sample_id, scheme_digest, result, "
                           "updated_at FROM analyses")
        connection.execute("DROP TABLE analyses")
        connection.execute("ALTER TABLE legacy RENAME TO analyses")
        connection.execute("PRAGMA user_version = 3")
    with Project(path) as project:
        assert project.analysis_results(sample_id) == stored
        assert project.latest_analysis(sample_id, "mlst")["st"] == "258"
        assert project.latest_analysis(sample_id, "cgmlst")["scheme_digest"] == "cg-v1"
        assert any(event["action"] == "schema_migrated" and event["details"]["to"] == 4
                   for event in project.history())
