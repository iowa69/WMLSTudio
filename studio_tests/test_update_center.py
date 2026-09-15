"""The Update page's own words: measured sizes, and the sentences a probe may say.

These are the pieces that decide what a row claims, so they are tested without a
window: a size is measured or it is not published, "not checked" never becomes "up
to date", and a licence nobody read is reported rather than quietly dropped.
"""

from wmlstudio.update_center import (
    check_record,
    decisions_needed,
    describe_set_update,
    describe_size,
    describe_update,
    describe_version,
    human_size,
    installed_path,
    measure_tree,
    probe,
    published_size,
    size_note,
)

FAILED = {"ok": False, "when": "2026-09-15 10:00", "error": "no route to host",
          "last_success": "2026-09-01 08:30"}
NEVER = {"ok": False, "when": "2026-09-15 10:00", "error": "no route to host",
         "last_success": ""}
REACHED = {"ok": True, "when": "2026-09-15 10:00", "error": "",
           "message": "The installed reference release 2026-08-07.1 is the one NCBI publishes."}


def tree(root, files=3, size=1024):
    root.mkdir(parents=True, exist_ok=True)
    for index in range(files):
        (root / f"file{index}.bin").write_bytes(b"x" * size)
    return root


def test_a_size_is_measured_on_this_computer_or_it_is_not_claimed(tmp_path):
    measured, complete = measure_tree(tree(tmp_path / "store", files=4, size=2048))
    assert (measured, complete) == (4 * 2048, True)
    # A folder nobody has is not a folder of size zero.
    assert measure_tree(tmp_path / "absent") == (None, True)
    assert measure_tree("") == (None, True)
    single = tmp_path / "one.bin"
    single.write_bytes(b"y" * 10)
    assert measure_tree(single) == (10, True)


def test_a_measurement_that_stopped_early_says_so_instead_of_undercounting(tmp_path):
    """A cgMLST library is hundreds of thousands of files; "at least" is the honest word."""
    measured, complete = measure_tree(tree(tmp_path / "library", files=9, size=100), limit=3)
    assert complete is False
    assert 0 < measured < 9 * 100
    item = {"size_bytes": measured, "size_complete": False}
    assert describe_size(item).startswith("at least ")
    assert "stopped early" in size_note(item) and "larger than the number shown" in size_note(item)


def test_sizes_read_as_room_taken_here_or_as_a_download_never_as_a_guess():
    assert human_size(0) == ""
    assert human_size(2048) == "2.0 KB"
    assert human_size(15 * 1024 ** 2) == "15 MB"
    assert describe_size({"size_bytes": 3 * 1024 ** 2, "size_complete": True}) == "3.0 MB here"
    # Not installed: the only honest number is the one the item's own action states.
    published = {"action": {"detail": "Downloads 18 pinned genomes (about 17 MB) into the folder."}}
    assert published_size(published) == "17 MB"
    assert describe_size(published) == "about 17 MB to download"
    assert describe_size({"action": {"detail": "no size here"}}) == "Not published"


def test_installed_path_reads_the_probe_rather_than_inventing_a_folder():
    assert installed_path({"detail": {"path": "/a/b"}}) == "/a/b"
    assert installed_path({"detail": {"roots": ["/c", "/d"]}}) == "/c"
    assert installed_path({"detail": {}}) == ""
    assert installed_path({}) == ""
    # Path("") is the working directory, so a thing nobody has installed used to
    # be given a "last changed" date of today, read from whatever folder the
    # application happened to be started in.
    assert describe_version({"key": "assembly_runtime", "detail": {}}) == "Not recorded"
    assert "last changed" not in describe_version({"key": "blast_tools", "detail": {}})


def test_a_version_column_never_shows_what_a_download_would_have_brought():
    """"Not installed · 18 reference genomes" read as though eighteen were here."""
    missing = {"key": "species_panel", "state": "missing", "ready": False,
               "detail": {"expected_revision": "panel-2026-09-14.2", "reference_count": 18}}
    said = describe_version(missing)
    assert said == "Nothing installed · this build expects revision panel-2026-09-14.2"
    assert "18" not in said
    installed = {"key": "species_panel", "state": "ready", "ready": True,
                 "detail": {"source_revision": "panel-2026-09-14.2", "reference_count": 18}}
    assert "18 reference genomes" in describe_version(installed)


def test_a_missing_thing_is_reported_from_this_computer_not_from_a_provider():
    item = {"key": "species_panel", "ready": False, "detail": {},
            "action": {"label": "Install the species reference panel"}}
    assert describe_update(item) == "Install the species reference panel"
    assert describe_update(item, REACHED) == "Install the species reference panel"


def test_a_revision_this_build_expects_is_never_presented_as_an_online_check():
    item = {"key": "species_panel", "ready": True,
            "detail": {"source_revision": "2026-01-01", "expected_revision": "2026-09-01"}}
    said = describe_update(item)
    assert "Update available" in said and "2026-09-01" in said
    assert "not with a provider" in said
    matching = {"key": "species_panel", "ready": True,
                "detail": {"source_revision": "2026-09-01", "expected_revision": "2026-09-01"}}
    assert "No provider was asked" in describe_update(matching)
    assert "up to date" not in describe_update(matching).casefold()
    # A part that ships inside the application, and a scheme library nobody
    # version-compares, each say what "an update" means for them rather than
    # borrowing the sentence that belongs to a set with a published release.
    assert "not from here" in describe_update({"key": "blast_tools", "ready": True, "detail": {}})
    assert "one scheme at a time" in describe_update(
        {"key": "cgmlst_schemes", "ready": True, "detail": {}})
    assert "cannot be confirmed here" in describe_update(
        {"key": "characterization", "ready": True, "detail": {}})


def test_a_failed_check_never_reads_as_up_to_date_and_names_the_last_one_that_worked():
    item = {"key": "hydra_database", "ready": True, "detail": {}}
    assert "Not yet checked" in describe_update(item)
    failed = describe_update(item, FAILED)
    assert "The check failed" in failed and "no route to host" in failed
    assert "2026-09-01 08:30" in failed
    assert "up to date" not in failed.casefold()
    assert "No check has ever reached a provider" in describe_update(item, NEVER)
    assert describe_update(item, REACHED) == REACHED["message"]


def test_a_reference_set_whose_terms_nobody_read_is_reported_not_skipped():
    restricted = {"name": "vfdb_full", "title": "VFDB full dataset", "installed": False,
                  "open_licence": False, "licence": "academic use only",
                  "provider": "VFDB", "download": "automatic", "purpose": "virulence factors"}
    said = describe_set_update(restricted)
    assert "Needs your decision" in said
    assert "academic use only" in said and "VFDB" in said
    assert "will not accept those terms for you" in said
    assert "never part of" in said


def test_a_set_this_engine_cannot_fetch_says_where_it_comes_from():
    manual = {"name": "card", "title": "CARD", "installed": False, "open_licence": True,
              "licence": "CC-BY", "provider": "McMaster", "download": "by hand",
              "url": "https://card.mcmaster.ca", "purpose": "resistance genes"}
    said = describe_set_update(manual)
    assert "cannot fetch it" in said and "https://card.mcmaster.ca" in said


def test_an_installed_set_says_which_question_was_answered_about_it():
    installed = {"name": "ncbi", "title": "NCBI", "installed": True, "open_licence": True,
                 "licence": "public domain", "provider": "NCBI", "download": "automatic",
                 "version": "2026-08-07.1", "purpose": "resistance genes"}
    assert "Not yet checked" in describe_set_update(installed)
    failed = describe_set_update(installed, None, FAILED)
    assert "The check failed" in failed and "unknown" in failed
    step = {"action": "update", "reason": "NCBI publishes 2026-09-01.1 today."}
    assert describe_set_update(installed, step, REACHED) == step["reason"]
    # A provider with no published version cannot confirm anything, and says so
    # instead of borrowing the phrasing of a set that was really checked.
    other = dict(installed, name="vfdb_core", provider="VFDB")
    said = describe_set_update(other, {"action": "skip", "reason": "already current"}, REACHED)
    assert "publishes no version" in said and "cannot be confirmed here" in said


def test_install_everything_names_the_sets_it_refuses_to_accept_terms_for():
    plan = {"databases": ["ncbi"], "steps": [
        {"kind": "requirement", "key": "species_panel", "action": "install"},
        {"kind": "database", "name": "ncbi", "action": "update"},
        {"kind": "database", "name": "vfdb_full", "action": "skip"},
        {"kind": "database", "name": "resfinder", "action": "skip"},
    ], "catalogue": {"entries": [
        {"name": "ncbi", "title": "NCBI", "installed": True, "open_licence": True,
         "licence": "public domain", "provider": "NCBI", "url": "", "size": "14 MB",
         "purpose": "resistance genes"},
        {"name": "vfdb_full", "title": "VFDB full dataset", "installed": False,
         "open_licence": False, "licence": "academic use only", "provider": "VFDB",
         "url": "http://vfdb", "size": "about 70 MB", "purpose": "virulence factors"},
        {"name": "resfinder", "title": "ResFinder", "installed": False, "open_licence": True,
         "licence": "Apache 2.0", "provider": "DTU", "url": "", "size": "",
         "purpose": "acquired resistance"},
    ]}}
    waiting = decisions_needed(plan)
    assert [row["name"] for row in waiting] == ["vfdb_full"]
    assert waiting[0]["licence"] == "academic use only"
    # A set somebody asked for by name is being installed, not left out, and saying
    # it was left out would be its own untruth.
    assert decisions_needed(dict(plan, databases=["ncbi", "vfdb_full"])) == []
    assert decisions_needed({}) == []


def test_the_last_successful_check_outlives_the_window(tmp_path):
    assert check_record(tmp_path) is None
    stored = check_record(tmp_path, {"when": "2026-09-01 08:30", "latest": "2026-08-07.1"})
    assert check_record(tmp_path) == stored
    assert check_record(tmp_path)["when"] == "2026-09-01 08:30"


def test_the_probe_reads_this_computer_and_asks_nobody_anything(tmp_path, monkeypatch):
    """Every sentence on the page starts here, and this step must never open a socket."""
    import urllib.request

    def refuse(*args, **options):
        raise AssertionError("the probe must not contact a server")

    monkeypatch.setattr(urllib.request, "urlopen", refuse)
    data = probe(str(tmp_path / "data"), [], "")
    assert data["items"] and data["verified"] is False
    assert "catalogue" in data and data["catalogue"]["entries"]
    for item in data["items"]:
        assert "size_bytes" in item and "size_complete" in item
        assert describe_size(item)
        # Nothing here has been compared with a provider, and no row pretends it has.
        assert "up to date" not in describe_update(item).casefold()
