"""The installation report must diagnose a silent gap, never assume a working one."""

import hashlib
import importlib.util
import json
import sys
from pathlib import Path

import pytest

from wmlstudio import characterization_refs, hydra_runtime, provisioning
from wmlstudio.organism_panel import PANEL_DIRECTORY, PANEL_PREFIX, panel_digest
from wmlstudio.project import Project

CLI = Path(__file__).resolve().parents[1] / "studio_scripts" / "check_setup.py"
_SPEC = importlib.util.spec_from_file_location("check_setup_test", CLI)
check_setup = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(check_setup)


def _publish(snapshot, manifest):
    """Hash every staged file into the manifest exactly as provisioning does."""
    for entry in manifest["files"]:
        target = snapshot / entry["path"]
        entry["bytes"] = target.stat().st_size
        entry["sha256"] = hashlib.sha256(target.read_bytes()).hexdigest()
    manifest["files"].sort(key=lambda item: item["path"])
    manifest["reference_digest"] = panel_digest(manifest)
    (snapshot / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return snapshot


def species_panel(data_root, taxa=(("Klebsiella", "pneumoniae"), ("Escherichia", "coli"))):
    """A real, verifiable species panel snapshot in the folder the application reads."""
    store = Path(data_root) / PANEL_DIRECTORY
    manifest = {"format_version": 1, "source_revision": "test-species-panel", "species": [],
                "virulence": {}, "files": [], "panel_kind": "species_panel"}
    for index, (genus, species) in enumerate(taxa):
        accession = f"GCF_0000000{index:02d}.1"
        relative = f"species/{accession}.fna.gz"
        manifest["species"].append({"id": accession, "path": relative, "genus": genus,
                                    "species": species, "subspecies": "", "outgroup": True,
                                    "taxonomy_basis": "synthetic fixture"})
        manifest["files"].append({"path": relative})
    snapshot = store / (PANEL_PREFIX + "0" * 20)
    for entry in manifest["files"]:
        target = snapshot / entry["path"]
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(b"synthetic reference bytes " + entry["path"].encode())
    published = _publish(snapshot, manifest)
    # The folder name carries the digest, the way an installed panel's does.
    final = store / (PANEL_PREFIX + manifest["reference_digest"][:20])
    published.rename(final)
    return final


def characterization_snapshot(root, *, format_version=1):
    """A starter snapshot; format 2 additionally carries the organism-module sections."""
    root = Path(root)
    manifest = {"format_version": format_version, "source_revision": "test-characterization",
                "species": [{"id": "GCF_000016305.1", "path": "species/kp.fna.gz",
                             "genus": "Klebsiella", "species": "pneumoniae", "subspecies": "",
                             "outgroup": False, "taxonomy_basis": "synthetic fixture"}],
                "virulence": {"ybt": [{"gene": "ybtS", "path": "virulence/ybt/ybtS.fasta"}]},
                "files": [{"path": "species/kp.fna.gz"}, {"path": "virulence/ybt/ybtS.fasta"}]}
    if format_version >= 2:
        manifest.update({
            "sources": {"kleborate": {"repository": "https://example.invalid", "revision": "x"}},
            "locus_profiles": {"ybt": {"module": "klebsiella__ybst", "st_field": "ST",
                                       "path": "virulence/ybt/profiles.tsv",
                                       "lineage_field": "lineage"}},
            "sccmec": {"targets": [{"gene": "mecA", "path": "sccmec/targets/mecA.fasta"}],
                       "regions": [], "rules": {"targets": ["mecA"], "aliases": [], "types": []}},
            "capsule": {"loci": [{"gene": "wzi", "path": "capsule/wzi.fasta"}]}})
        manifest["files"] += [{"path": "virulence/ybt/profiles.tsv"},
                              {"path": "sccmec/targets/mecA.fasta"}, {"path": "capsule/wzi.fasta"}]
    for entry in manifest["files"]:
        target = root / entry["path"]
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(b"synthetic " + entry["path"].encode())
    return _publish(root, manifest)


def hydra_store(root, names=("ncbi", "protein"),
                organisms=("Escherichia", "Klebsiella_pneumoniae"), release="2026-08-07.1",
                staged="2026-09-12 18:41:37"):
    """A HYDRA database store in the shape hydra_runtime.installed_databases reads."""
    root = Path(root)
    entries = {}
    for name in names:
        target = root / "stores" / name
        target.mkdir(parents=True, exist_ok=True)
        (target / "sequences.fna").write_text(">ref\nACGT\n")
        entries[name] = {"path": f"stores/{name}", "version": release, "kind": "nucl",
                         "installed": staged}
        if name == "protein":
            entries[name]["organisms"] = list(organisms)
    (root / "manifest.json").write_text(json.dumps({"hydra_db_version": 1, "databases": entries}))
    return root


def scheme_folder(root, name, loci, *, declared=""):
    folder = Path(root) / name
    folder.mkdir(parents=True)
    for index in range(loci):
        (folder / f"{name}_locus{index}.tfa").write_text(f">{name}_locus{index}_1\nACGTACGT\n")
    metadata = {"name": name, "organism": "Testella testis"}
    if declared:
        metadata["type"] = declared
    (folder / f"{name}_info.json").write_text(json.dumps(metadata))
    return folder


@pytest.fixture
def bare(monkeypatch):
    """An installation that staged nothing: no bundled AMR store, no characterization."""
    monkeypatch.setattr(hydra_runtime, "bundled_database_root", lambda: None)
    monkeypatch.setattr(characterization_refs, "bundled_reference_root", lambda: None)


def probe(tmp_path, **options):
    options.setdefault("scheme_paths", [])
    options.setdefault("hydra_database_root", str(tmp_path / "no-database"))
    return provisioning.report(data_root=tmp_path / "data", **options)


def item_of(report, key):
    return next(item for item in report["items"] if item["key"] == key)


def test_an_unprovisioned_installation_names_the_reason_samples_are_never_filed(tmp_path, bare):
    report = probe(tmp_path)
    panel = item_of(report, "species_panel")
    assert panel["state"] == "missing" and panel["required"] is True
    assert "unresolved" in panel["consequence"] and "organism folder" in panel["consequence"]
    assert "no genomic comparison" in panel["reason"]
    assert report["ready"] is False and report["blocking"][0] == "species_panel"
    assert report["capabilities"]["organism_filing"]["ready"] is False
    # The single next thing to do is an action a caller can actually carry out.
    assert report["next_action"]["entry_point"] == "wmlstudio.provisioning:install_species_panel"
    assert report["next_action"]["automatic"] is True
    assert report["summary"].startswith("Not yet")
    assert "install the species reference panel" in report["summary"]


def test_the_species_panel_is_ready_only_once_it_is_installed_and_verified(tmp_path, bare):
    assert item_of(probe(tmp_path), "species_panel")["state"] == "missing"
    installed = species_panel(tmp_path / "data")
    panel = item_of(probe(tmp_path), "species_panel")
    assert panel["state"] == "ready" and panel["detail"]["path"] == str(installed)
    assert panel["detail"]["taxa"] == ["Escherichia coli", "Klebsiella pneumoniae"]
    assert "re-hashed" in panel["probe"]


def test_a_tampered_panel_is_refused_and_a_quick_check_says_it_did_not_hash(tmp_path, bare):
    installed = species_panel(tmp_path / "data")
    reference = next(installed.glob("species/*.fna.gz"))
    # Same length, different bytes: only a content check can see this.
    reference.write_bytes(b"X" * reference.stat().st_size)
    verified = item_of(probe(tmp_path), "species_panel")
    assert verified["state"] == "unusable" and "did not pass validation" in verified["reason"]
    quick = item_of(probe(tmp_path, verify=False), "species_panel")
    assert quick["state"] == "ready"
    assert "recorded size" in quick["probe"] and "re-hashed" not in quick["probe"]


def test_a_quick_check_still_refuses_a_panel_whose_files_are_gone(tmp_path, bare):
    installed = species_panel(tmp_path / "data")
    next(installed.glob("species/*.fna.gz")).unlink()
    quick = item_of(probe(tmp_path, verify=False), "species_panel")
    assert quick["state"] == "unusable" and "missing or a different size" in quick["reason"]


def test_classical_mlst_and_cgmlst_are_reported_as_two_separate_capabilities(tmp_path, bare):
    schemes = tmp_path / "schemes"
    scheme_folder(schemes, "sevenlocus", 7)
    report = probe(tmp_path, scheme_paths=sorted(schemes.iterdir()))
    classical, core = item_of(report, "mlst_schemes"), item_of(report, "cgmlst_schemes")
    assert classical["state"] == "ready" and classical["detail"]["count"] == 1
    assert core["state"] == "missing" and core["detail"]["count"] == 0
    assert core["required"] is False
    assert "classical 7-locus MLST is unaffected" in core["consequence"]
    assert report["capabilities"]["mlst"]["ready"] is True
    assert report["capabilities"]["cgmlst"]["ready"] is False
    assert "cgmlst_schemes" in report["capabilities"]["cgmlst"]["blocked_by"]

    scheme_folder(schemes, "coregenome", 40)
    scheme_folder(schemes, "declared", 5, declared="cgMLST")
    report = probe(tmp_path, scheme_paths=sorted(schemes.iterdir()))
    assert item_of(report, "mlst_schemes")["detail"]["count"] == 1
    assert item_of(report, "cgmlst_schemes")["detail"]["count"] == 2
    assert item_of(report, "cgmlst_schemes")["state"] == "ready"


def test_missing_blast_names_every_executable_and_every_place_it_looked():
    tools = dict.fromkeys(hydra_runtime.TOOLS, None)
    absent = provisioning.blast_requirement(tools, verify=False)
    assert absent.state == "missing" and absent.required is True
    assert all(name in absent.reason for name in hydra_runtime.TOOLS)
    assert absent.detail["searched_directories"] and "on PATH" in absent.reason
    assert "HYDRA" in absent.consequence and "cgMLST" in absent.consequence

    partial = provisioning.blast_requirement({**tools, "blastn": "/somewhere/blastn"}, verify=False)
    assert partial.state == "partial" and partial.detail["found"] == {"blastn": "/somewhere/blastn"}


def test_present_blast_binaries_are_not_called_ready_until_one_is_actually_run():
    # sys.executable exists and refuses "-version": a resolvable tool that fails.
    tools = dict.fromkeys(hydra_runtime.TOOLS, sys.executable)
    assumed = provisioning.blast_requirement(tools, verify=False)
    assert assumed.state == "ready" and "were not run" in assumed.reason
    executed = provisioning.blast_requirement(tools, verify=True)
    assert executed.state == "unusable" and "-version" in executed.probe


def test_a_missing_assembler_blocks_reads_without_touching_imported_assemblies(monkeypatch):
    from wmlstudio import assembly
    def refuse(*args):
        raise assembly.AssemblyError("no SKESA here")

    monkeypatch.setattr(assembly, "resolve_skesa", refuse)
    absent = provisioning.assembly_runtime_requirement()
    assert absent.state == "missing" and absent.reason == "no SKESA here"
    assert "never typed" in absent.consequence and "unaffected" in absent.consequence
    monkeypatch.setattr(assembly, "resolve_skesa", lambda *args: Path("/tools/skesa"))
    assert provisioning.assembly_runtime_requirement().state == "ready"


def test_a_portable_build_missing_its_hydra_worker_is_not_called_ready(monkeypatch, tmp_path):
    capabilities = {"hydra_version": hydra_runtime.HYDRA_VERSION}
    monkeypatch.setattr(provisioning, "hydra_worker_path", lambda: tmp_path / "WMLSTudio-HYDRA")
    broken = provisioning.hydra_engine_requirement(capabilities)
    assert broken.state == "unusable" and "worker program is missing" in broken.reason
    assert broken.action.automatic is False
    (tmp_path / "WMLSTudio-HYDRA").write_text("#!/bin/sh\n")
    assert provisioning.hydra_engine_requirement(capabilities).state == "ready"


def test_the_hydra_engine_is_reported_apart_from_its_reference_data():
    missing = provisioning.hydra_engine_requirement({"hydra_version": None, "tools": {}})
    assert missing.state == "missing"
    wrong = provisioning.hydra_engine_requirement({"hydra_version": "0.0.1", "tools": {}})
    assert wrong.state == "unusable" and hydra_runtime.HYDRA_VERSION in wrong.reason
    ready = provisioning.hydra_engine_requirement({"hydra_version": hydra_runtime.HYDRA_VERSION})
    assert ready.state == "ready"


def test_an_installed_amr_store_lists_the_organisms_its_catalogs_cover(tmp_path, bare):
    store = hydra_store(tmp_path / "amr")
    report = probe(tmp_path, hydra_database_root=str(store))
    database = item_of(report, "hydra_database")
    assert database["state"] == "ready"
    assert sorted(database["detail"]["databases"]) == ["ncbi", "protein"]
    # This list is what an organism chooser may offer for HYDRA; nothing else.
    assert database["detail"]["organisms"] == ["Escherichia", "Klebsiella_pneumoniae"]


def test_a_store_without_the_protein_set_is_partial_and_says_what_it_cannot_report(tmp_path, bare):
    store = hydra_store(tmp_path / "amr", names=("ncbi",))
    database = item_of(probe(tmp_path, hydra_database_root=str(store)), "hydra_database")
    assert database["state"] == "partial"
    assert "protein" in database["reason"] and "point mutation" in database["reason"]


def test_an_empty_selection_is_pointed_at_the_store_already_installed(tmp_path, monkeypatch):
    staged = hydra_store(tmp_path / "staged")
    monkeypatch.setattr(hydra_runtime, "bundled_database_root", lambda: staged)
    monkeypatch.setattr(characterization_refs, "bundled_reference_root", lambda: None)
    empty = tmp_path / "chosen"
    empty.mkdir()
    database = item_of(probe(tmp_path, hydra_database_root=str(empty)), "hydra_database")
    assert database["state"] == "missing"
    assert str(staged) in database["reason"] and "not the one selected" in database["reason"]
    assert database["action"]["entry_point"] == "wmlstudio.provisioning:select_hydra_database"
    assert database["action"]["argument"] == str(staged)


def test_a_broken_database_manifest_is_unusable_rather_than_silently_empty(tmp_path, bare):
    store = tmp_path / "amr"
    store.mkdir()
    (store / "manifest.json").write_text(json.dumps({"databases": {"ncbi": {"path": "../escape"}}}))
    database = item_of(probe(tmp_path, hydra_database_root=str(store)), "hydra_database")
    assert database["state"] == "unusable" and "leaves the selected store" in database["reason"]


def test_selecting_a_store_records_only_one_that_actually_holds_a_database(tmp_path):
    project = Project(tmp_path / "investigation.wmlstudio")
    empty = tmp_path / "empty"
    empty.mkdir()
    with pytest.raises(ValueError, match="no readable HYDRA database manifest"):
        provisioning.select_hydra_database(project, empty)
    assert project.get_setting("hydra_database_root") is None
    store = hydra_store(tmp_path / "amr")
    recorded = provisioning.select_hydra_database(project, store)
    assert project.get_setting("hydra_database_root") == recorded == str(store.resolve())
    project.close()


def test_a_format_one_snapshot_is_partial_and_names_the_assays_it_cannot_run(tmp_path, bare):
    root = characterization_snapshot(tmp_path / "characterization", format_version=1)
    item = item_of(probe(tmp_path, characterization_root=root), "characterization")
    assert item["state"] == "partial" and item["required"] is False
    assert item["detail"]["manifest_format_version"] == 1
    assert "sccmec.targets" in item["reason"]
    unavailable = [module for module in item["detail"]["modules"] if not module["available"]]
    assert {module["key"] for module in unavailable} >= {"sccmec", "klebsiella_capsule"}


def test_a_format_two_snapshot_reports_every_organism_specific_assay_as_available(tmp_path, bare):
    root = characterization_snapshot(tmp_path / "characterization", format_version=2)
    item = item_of(probe(tmp_path, characterization_root=root), "characterization")
    assert item["state"] == "ready" and item["detail"]["manifest_format_version"] == 2
    assert all(module["available"] for module in item["detail"]["modules"])
    assert all(module["taxa"] for module in item["detail"]["modules"])


def test_an_absent_characterization_snapshot_never_blocks_typing(tmp_path, bare):
    report = probe(tmp_path)
    item = item_of(report, "characterization")
    assert item["state"] == "missing" and item["required"] is False
    assert "characterization" not in report["blocking"]


def test_the_report_is_plain_data_a_user_interface_can_show_verbatim(tmp_path, bare):
    report = probe(tmp_path)
    assert json.loads(json.dumps(report)) == report
    assert [item["key"] for item in report["items"]] == list(provisioning.ORDER)
    for item in report["items"]:
        assert item["state"] in provisioning.STATES
        assert item["reason"] and item["probe"] and item["consequence"].startswith("Without this,")
        assert item["ready"] == (item["state"] == "ready")
        if not item["ready"]:
            assert item["action"]["label"]


def requirement(key, state, *, required=True, label="Install it"):
    return provisioning.Requirement(
        key=key, title=key.replace("_", " ").title(), capability="Some capability",
        blocks="nothing at all happens", required=required, state=state,
        reason="because", probe="looked", action=provisioning.Action(label=label))


def test_the_one_sentence_summary_answers_can_i_work_now_and_what_to_do_first():
    ready = [requirement("species_panel", "ready"), requirement("cgmlst_schemes", "ready",
                                                                required=False)]
    assert provisioning.summarize(ready).startswith("Yes, everything this application needs")

    optional = [requirement("species_panel", "ready"),
                requirement("cgmlst_schemes", "missing", required=False,
                            label="Download a cgMLST scheme")]
    sentence = provisioning.summarize(optional)
    assert sentence.startswith("Yes, you can import samples")
    assert "download a cgMLST scheme" in sentence

    blocked = [requirement("species_panel", "partial"),
               requirement("blast_tools", "missing")]
    assert provisioning.summarize(blocked).startswith("Not yet: nothing at all happens;")
    # One sentence, so a status bar or a banner can show it unedited.
    assert provisioning.summarize(blocked).count(".") == 1


def test_a_requirement_cannot_be_given_a_state_the_interface_does_not_understand():
    with pytest.raises(ValueError, match="Unknown provisioning state"):
        requirement("species_panel", "probably fine")


def test_installing_the_panel_targets_the_folder_the_report_named(tmp_path, monkeypatch):
    from wmlstudio import organism_panel
    seen = {}
    monkeypatch.setattr(organism_panel, "provision_species_panel",
                        lambda root, **options: seen.setdefault("root", Path(root)))
    provisioning.install_species_panel(tmp_path / "data")
    assert seen["root"] == tmp_path / "data" / PANEL_DIRECTORY


def test_the_command_line_check_reports_the_same_gaps_and_exits_non_zero(tmp_path, bare, capsys):
    code = check_setup.main(["--json", "--data-root", str(tmp_path / "data"),
                             "--hydra-database", str(tmp_path / "absent"), "--quick"])
    printed = json.loads(capsys.readouterr().out)
    assert code == 1
    assert printed["ready"] is False and printed["blocking"][0] == "species_panel"
    assert printed["verified"] is False


def test_the_command_line_check_prints_the_next_step_for_a_non_expert(tmp_path, bare, capsys):
    assert check_setup.main(["--data-root", str(tmp_path / "data"),
                             "--hydra-database", str(tmp_path / "absent"), "--quick"]) == 1
    printed = capsys.readouterr().out
    assert "[ MISSING ] Species reference panel (required)" in printed
    assert "Next step: Install the species reference panel" in printed
    assert "Samples > Reference data > Species panel" in printed
    assert "Can you run an investigation right now?" in printed


def test_the_command_line_check_exits_zero_when_nothing_required_is_missing(tmp_path, bare,
                                                                            monkeypatch, capsys):
    ready = {**probe(tmp_path), "ready": True, "blocking": [], "next_action": None}
    monkeypatch.setattr(check_setup, "report", lambda **options: ready)
    assert check_setup.main([]) == 0
    assert "installation check" in capsys.readouterr().out


def test_the_installed_amr_release_is_reported_with_its_date_and_its_age(tmp_path, bare):
    current = item_of(probe(tmp_path, hydra_database_root=str(hydra_store(tmp_path / "now"))),
                      "hydra_database")
    assert current["state"] == "ready"
    assert current["detail"]["release"] == "2026-08-07.1"
    assert current["detail"]["staged"].startswith("2026-09-12")
    assert current["detail"]["stale"] is False and current["detail"]["age_days"] >= 0
    assert "reference release 2026-08-07.1" in current["reason"]
    assert "on 2026-09-12" in current["reason"] and "days old" in current["reason"]
    # Each store says what it is searched for, so an absent one is an explainable gap.
    assert "blastn" in current["detail"]["purposes"]["ncbi"]
    assert "point-mutation" in current["detail"]["purposes"]["protein"]

    old = item_of(probe(tmp_path, hydra_database_root=str(
        hydra_store(tmp_path / "old", release="2019-01-01.1"))), "hydra_database")
    assert old["state"] == "ready" and old["detail"]["stale"] is True
    assert "determinants named after that release are not in it" in old["reason"].lower()
    # The age is whatever it is today; the point is that the number is in the sentence.
    assert f"{old['detail']['age_days']} days old" in old["reason"]
    assert old["detail"]["age_days"] > hydra_runtime.REFERENCE_AGE_DAYS
    # A stale-but-working store is offered an update, not a first-time install.
    assert old["action"]["entry_point"] == "wmlstudio.provisioning:update_hydra_databases"


def test_a_hydra_run_is_refused_with_the_missing_piece_and_the_action_that_fixes_it(tmp_path,
                                                                                    monkeypatch):
    monkeypatch.setattr(hydra_runtime, "bundled_database_root", lambda: None)
    monkeypatch.setattr(hydra_runtime, "runtime_capabilities", lambda db_root=None: {
        "hydra_version": hydra_runtime.HYDRA_VERSION, "assembly_available": True, "databases": {},
        "tools": {name: f"/tools/{name}" for name in hydra_runtime.TOOLS}, "limitations": []})
    empty = tmp_path / "store"
    empty.mkdir()
    check = provisioning.hydra_prerequisites(tmp_path / "data", str(empty))
    assert check["ready"] is False
    assert "hydra_database" in check["blocking"]
    assert "every gene and mutation this screen can name" in check["message"]
    assert check["action"]["automatic"] is True
    assert check["action"]["entry_point"] == "wmlstudio.provisioning:install_hydra_databases"
    assert json.loads(json.dumps(check)) == check


def test_a_ready_hydra_gate_still_says_what_the_run_will_not_be_able_to_report(tmp_path,
                                                                               monkeypatch):
    monkeypatch.setattr(hydra_runtime, "bundled_database_root", lambda: None)
    monkeypatch.setattr(hydra_runtime, "runtime_capabilities", lambda db_root=None: {
        "hydra_version": hydra_runtime.HYDRA_VERSION, "assembly_available": True, "databases": {},
        "tools": {name: f"/tools/{name}" for name in hydra_runtime.TOOLS},
        "limitations": ["Gene and mutation evidence is not a susceptibility phenotype."]})
    store = hydra_store(tmp_path / "amr")
    check = provisioning.hydra_prerequisites(tmp_path / "data", str(store),
                                             organism="Klebsiella pneumoniae")
    assert check["ready"] is True and check["blocking"] == []
    assert check["databases"] == ["ncbi", "protein"]
    assert check["organism"]["resolved"] == "Klebsiella_pneumoniae"
    assert check["database"]["release"] == "2026-08-07.1"
    assert check["limitations"] and "susceptibility phenotype" in check["limitations"][0]

    outside = provisioning.hydra_prerequisites(tmp_path / "data", str(store),
                                               organism="Listeria monocytogenes")
    assert outside["ready"] is True and outside["organism"]["resolved"] == ""
    assert any("not an absence of mutations" in text for text in outside["warnings"])


def test_updating_downloads_and_starts_using_the_new_store_in_one_action(tmp_path, monkeypatch):
    """install publishes a snapshot; nobody discovers they must also select it."""
    project = Project(tmp_path / "investigation.wmlstudio")
    published = hydra_store(tmp_path / "published")
    monkeypatch.setattr(provisioning, "install_hydra_databases",
                        lambda root, names, **options: {"database_root": str(published),
                                                        "previous_database_root": str(root)})
    result = provisioning.update_hydra_databases(tmp_path / "downloads", project=project)
    assert result["selected"] is True
    assert project.get_setting("hydra_database_root") == str(published.resolve())
    assert "reference release 2026-08-07.1" in result["summary"]
    assert "keep the snapshot they were run against" in result["summary"]

    # Without a project there is nothing to record it on, and it says so.
    unrecorded = provisioning.update_hydra_databases(tmp_path / "downloads")
    assert unrecorded["selected"] is False
    assert "Select it for a project" in unrecorded["summary"]
    project.close()


def test_the_whole_reference_catalogue_is_listed_not_only_what_is_installed(tmp_path, bare):
    """HYDRA needs inputs it never names; a user cannot choose from a list they cannot see."""
    store = hydra_store(tmp_path / "amr")
    catalogue = provisioning.hydra_database_catalogue(tmp_path / "data", str(store))
    rows = {entry["name"]: entry for entry in catalogue["entries"]}
    assert catalogue["installed"] == ["ncbi", "protein"]
    assert len(rows) > 10 and set(catalogue["available"]) == set(rows) - {"ncbi", "protein"}
    assert all(entry["purpose"] and entry["provider"] and entry["licence"]
               for entry in catalogue["entries"])
    assert "2 of" in catalogue["summary"] and "downloaded only when you ask" in catalogue["summary"]
    # The same count reaches the installation report, so the Update page can show it.
    detail = item_of(probe(tmp_path, hydra_database_root=str(store)), "hydra_database")["detail"]
    assert detail["catalogue"]["installed"] == ["ncbi", "protein"]
    assert detail["catalogue"]["total"] == len(rows)
    assert set(detail["catalogue"]["downloadable"]) <= set(rows)


def test_the_plan_says_what_is_missing_what_is_stale_and_what_is_left_alone(tmp_path, bare,
                                                                            monkeypatch):
    store = hydra_store(tmp_path / "amr")
    monkeypatch.setattr(provisioning, "hydra_update_available", lambda *args, **options: {
        "installed": "2026-08-07.1", "latest": "2026-09-01.1", "newer_available": True,
        "message": "NCBI publishes 2026-09-01.1 today.", "error": ""})
    plan = provisioning.installation_plan(tmp_path / "data", str(store))
    steps = {step["key"]: step for step in plan["steps"]}
    # Installed but superseded: an update, named with both releases.
    assert steps["database:ncbi"]["action"] == "update"
    assert "2026-09-01.1" in steps["database:ncbi"]["reason"]
    # Not installed and nobody asked for it: listed, and explicitly not fetched.
    assert steps["database:card"]["action"] == "skip"
    assert steps["database:card"]["state"] == "not selected"
    assert "no isolate is screened against it" in steps["database:card"]["reason"]
    # Missing, and something in this application can install it.
    assert steps["species_panel"]["action"] == "install"
    assert plan["databases"] == ["ncbi", "protein"]
    assert plan["work"] is True
    assert "published beside" in plan["summary"]
    assert json.loads(json.dumps(plan)) == plan


def test_a_reference_set_nobody_asked_for_is_listed_but_never_fetched(tmp_path, bare, monkeypatch):
    store = hydra_store(tmp_path / "amr")
    monkeypatch.setattr(provisioning, "hydra_update_available", lambda *args, **options: {
        "installed": "2026-08-07.1", "latest": "2026-08-07.1", "newer_available": False,
        "message": "", "error": ""})
    quiet = provisioning.installation_plan(tmp_path / "data", str(store))
    assert quiet["databases"] == []
    asked = provisioning.installation_plan(tmp_path / "data", str(store), databases=["card"])
    assert asked["databases"] == ["card"]
    # The snapshot that would be published still carries everything already installed:
    # fetching only the missing set would publish a store narrower than the one in use.
    assert asked["database_snapshot_contents"] == ["card", "ncbi", "protein"]


def test_a_set_whose_provider_publishes_no_version_is_judged_by_age_not_assumed_current(
        tmp_path, bare, monkeypatch):
    """Only NCBI publishes a release string; for the rest, "unknown" is not "current"."""
    store = hydra_store(tmp_path / "amr", names=("ncbi", "protein", "vfdb"),
                        staged="2019-01-04 09:00:00")
    monkeypatch.setattr(provisioning, "hydra_update_available", lambda *args, **options: {
        "installed": "2026-08-07.1", "latest": "2026-08-07.1", "newer_available": False,
        "message": "", "error": ""})
    steps = {step["key"]: step for step in
             provisioning.installation_plan(tmp_path / "data", str(store))["steps"]}
    aged = steps["database:vfdb"]
    assert aged["state"] == "stale" and aged["action"] == "update"
    assert "publishes no version this application can compare against" in aged["reason"]
    assert "days ago" in aged["reason"]
    # The NCBI sets are judged by the release NCBI publishes, not by their age.
    assert steps["database:ncbi"]["action"] == "skip"
    assert "is the one NCBI publishes today" in steps["database:ncbi"]["reason"]

    # When NCBI could not be reached, an old copy is not quietly called current.
    monkeypatch.setattr(provisioning, "hydra_update_available", lambda *args, **options: {
        "installed": "2026-08-07.1", "latest": "", "newer_available": False,
        "message": "", "error": "Could not reach NCBI."})
    offline = {step["key"]: step for step in
               provisioning.installation_plan(tmp_path / "data", str(store))["steps"]}
    assert offline["database:ncbi"]["action"] == "update"
    assert "current release could not be checked" in offline["database:ncbi"]["reason"]


def test_a_hand_installed_set_is_not_handed_to_a_downloader_that_cannot_fetch_it(tmp_path, bare,
                                                                                  monkeypatch):
    """Asking the engine for a set it has no fetcher for would fail the whole batch."""
    store = hydra_store(tmp_path / "amr", names=("ncbi", "protein", "vfdb_full"))
    monkeypatch.setattr(provisioning, "hydra_update_available", lambda *args, **options: {
        "installed": "2026-08-07.1", "latest": "2026-09-01.1", "newer_available": True,
        "message": "", "error": ""})
    plan = provisioning.installation_plan(tmp_path / "data", str(store))
    assert plan["installed_by_hand"] == ["vfdb_full"]
    assert "vfdb_full" not in plan["databases"]
    assert "vfdb_full" not in plan["database_snapshot_contents"]
    # And its absence from the new snapshot is stated rather than discovered later.
    assert "installed by hand" in plan["summary"] and "vfdb_full" in plan["summary"]


def test_a_licence_that_has_to_be_read_is_never_swept_into_install_everything(tmp_path, bare,
                                                                              monkeypatch):
    store = hydra_store(tmp_path / "amr")
    monkeypatch.setattr(provisioning, "hydra_update_available", lambda *args, **options: {
        "installed": "2026-08-07.1", "latest": "2026-08-07.1", "newer_available": False,
        "message": "", "error": ""})
    plan = provisioning.installation_plan(tmp_path / "data", str(store))
    restricted = [step for step in plan["steps"] if step["kind"] == "database"
                  and step.get("licence") and "academic" in step["licence"]]
    assert restricted and all(step["action"] == "skip" for step in restricted)
    assert all(step["name"] not in plan["databases"] for step in restricted)


def test_the_update_page_shows_the_plan_before_anything_is_downloaded(tmp_path, bare, monkeypatch):
    """"Update everything" is a leap of faith unless the list comes first."""
    from wmlstudio.update_center import UpdateCenter

    store = hydra_store(tmp_path / "amr")
    monkeypatch.setattr(provisioning, "hydra_update_available", lambda *args, **options: {
        "installed": "2026-08-07.1", "latest": "2026-09-01.1", "newer_available": True,
        "message": "NCBI publishes 2026-09-01.1 today.", "error": ""})
    text = UpdateCenter.plan_text(provisioning.installation_plan(tmp_path / "data", str(store)))
    assert "Will be installed:" in text and "Will be updated:" in text
    assert "Species reference panel" in text and "NCBI AMRFinderPlus" in text
    assert "Cannot be installed from here:" in text
    assert "NCBI publishes 2026-09-01.1 today." in text
    assert "sequences are never uploaded" in text and "published beside it" in text


def _no_network(monkeypatch, *, newer=False):
    monkeypatch.setattr(provisioning, "hydra_update_available", lambda *args, **options: {
        "installed": "2026-08-07.1", "latest": "2026-09-01.1" if newer else "2026-08-07.1",
        "newer_available": newer, "message": "", "error": ""})


def test_install_and_update_everything_is_safe_to_press_twice(tmp_path, bare, monkeypatch):
    """The second press must find nothing to do rather than download the same data again."""
    store = hydra_store(tmp_path / "amr")
    published = hydra_store(tmp_path / "published")
    calls = []
    _no_network(monkeypatch, newer=True)
    monkeypatch.setattr(provisioning, "install_species_panel",
                        lambda root=None, **options: calls.append("panel"))
    monkeypatch.setattr(provisioning, "update_hydra_databases",
                        lambda root, names, **options: calls.append(list(names)) or {
                            "database_root": str(published), "selected": False})
    first = provisioning.install_everything(tmp_path / "data", selected=str(store))
    assert "panel" in calls and ["ncbi", "protein"] in calls
    assert first["failed"] == [] and first["updated"]
    assert "Previous snapshots are kept" in first["summary"]

    calls.clear()
    species_panel(tmp_path / "data")
    _no_network(monkeypatch, newer=False)
    second = provisioning.install_everything(tmp_path / "data", selected=str(store))
    assert calls == [], "a second press must not re-download what is already current"
    assert second["installed"] == [] and second["updated"] == []
    assert "Nothing needed installing or updating" in second["summary"]


def test_install_everything_publishes_beside_the_store_an_analysis_already_used(tmp_path, bare,
                                                                                monkeypatch):
    project = Project(tmp_path / "investigation.wmlstudio")
    store = hydra_store(tmp_path / "amr")
    manifest = (store / "manifest.json").read_bytes()
    published = hydra_store(tmp_path / "published")
    _no_network(monkeypatch, newer=True)
    monkeypatch.setattr(provisioning, "install_species_panel", lambda root=None, **options: None)
    # Only the download itself is replaced: publishing beside, recording the new
    # snapshot and saying so are the behaviour under test.
    monkeypatch.setattr(provisioning, "install_hydra_databases",
                        lambda root, names, **options: {"database_root": str(published),
                                                        "previous_database_root": str(root)})
    result = provisioning.install_everything(tmp_path / "data", project=project,
                                             selected=str(store))
    assert result["database_root"] == str(published)
    assert project.get_setting("hydra_database_root") == str(published.resolve())
    # The store the previous analyses ran against is byte-for-byte where it was.
    assert (store / "manifest.json").read_bytes() == manifest
    assert "keep the snapshot they were run against" in result["summary"]
    assert json.loads(json.dumps(result)) == result
    project.close()


def test_a_failed_download_is_named_and_leaves_everything_else_that_worked(tmp_path, bare,
                                                                           monkeypatch):
    store = hydra_store(tmp_path / "amr")
    manifest = (store / "manifest.json").read_bytes()
    _no_network(monkeypatch, newer=True)
    monkeypatch.setattr(provisioning, "install_species_panel", lambda root=None, **options: None)

    def refuse(root, names, **options):
        raise ValueError("Downloading protein, ncbi failed; the store you are using is unchanged.")

    monkeypatch.setattr(provisioning, "update_hydra_databases", refuse)
    result = provisioning.install_everything(tmp_path / "data", selected=str(store))
    assert result["installed"] == ["species_panel"], "the step that worked is still recorded"
    assert [failure["key"] for failure in result["failed"]] == ["hydra_database"]
    assert "the store you are using is unchanged" in result["failed"][0]["error"]
    assert "unchanged" in result["failed"][0]["state"]
    assert "failed" in result["summary"]
    assert result["database_root"] == str(store.resolve())
    assert (store / "manifest.json").read_bytes() == manifest


def test_checking_for_a_newer_release_never_reports_up_to_date_on_a_failed_request(tmp_path,
                                                                                   monkeypatch):
    store = hydra_store(tmp_path / "amr")
    monkeypatch.setattr(hydra_runtime, "bundled_database_root", lambda: None)
    monkeypatch.setattr(hydra_runtime, "latest_release",
                        lambda **options: (_ for _ in ()).throw(OSError("no network")))
    offline = provisioning.hydra_update_available(tmp_path / "data", str(store))
    assert offline["newer_available"] is False and offline["latest"] == ""
    assert "Could not reach NCBI" in offline["message"]
    assert offline["installed"] == "2026-08-07.1"

    monkeypatch.setattr(hydra_runtime, "latest_release", lambda **options: "2026-09-01.1")
    newer = provisioning.hydra_update_available(tmp_path / "data", str(store))
    assert newer["newer_available"] is True
    assert "2026-09-01.1" in newer["message"] and "2026-08-07.1" in newer["message"]
    assert "not in what you have" in newer["message"]

    monkeypatch.setattr(hydra_runtime, "latest_release", lambda **options: "2026-08-07.1")
    current = provisioning.hydra_update_available(tmp_path / "data", str(store))
    assert current["newer_available"] is False
    assert "is the one NCBI publishes today" in current["message"]
