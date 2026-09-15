"""Runtime isolation, safe native invocation, immutable database publication."""

import json
import sys
import time
from pathlib import Path

import pytest
from PySide6.QtCore import Qt

from wmlstudio import hydra_runtime as runtime
from wmlstudio.sequence import AnalysisCancelled


def database(root):
    target = root / "nucl/ncbi"
    target.mkdir(parents=True)
    (target / "sequences.fna").write_text(">ref\nACGTACGT\n")
    (target / "meta.tsv").write_text("seqid\tgene\nref\ttest\n")
    (root / "manifest.json").write_text(json.dumps({"hydra_db_version": 1,
        "databases": {"ncbi": {"path": "nucl/ncbi", "version": "test-snapshot", "kind": "nucl"}}}))
    return root


def complete_database(root, *, release="2026-08-07.1", separator="/",
                      organisms=("Escherichia", "Klebsiella_pneumoniae", "Staphylococcus_aureus"),
                      protein_mutations=None, virulence=3, stress=2):
    """A store in the shape a real staged snapshot has: both sets, plus the catalogues.

    ``separator`` reproduces a snapshot staged on Windows, whose manifest records
    ``nucl\\ncbi`` and which must still read correctly wherever it is opened.
    ``organisms`` have a DNA mutation catalogue; ``protein_mutations`` defaults to
    the same names and is given separately where the two lists must differ, as
    they do in the real release.
    """
    protein_mutations = organisms if protein_mutations is None else protein_mutations
    for relative in ("nucl/ncbi", "prot/protein"):
        (root / relative).mkdir(parents=True)
        (root / relative / "sequences.fna").write_text(">ref\nACGTACGT\n")
    (root / "prot/protein/taxgroup.tsv").write_text(
        "#taxgroup\tgpipe_taxgroup\tnumber_of_nucl_ref_genes\n"
        + "".join(f"{name}\t{name}\t1\n" for name in (*organisms, "Pseudomonas_aeruginosa")))
    (root / "prot/protein/AMRProt-mutation.tsv").write_text(
        "#taxgroup\taccession_version\tmutation_position\n"
        + "".join(f"{name}\tWP_00{index}.1\t{index}\n"
                  for index, name in enumerate(protein_mutations)))
    (root / "prot/protein/AMRProt-suppress.tsv").write_text(
        "#taxgroup\tprotein_accession\n"
        + "".join(f"{name}\tWP_9{index}.1\n" for index, name in enumerate(organisms)))
    rows = [("AMR", "AMR")] * 5 + [("VIRULENCE", "VIRULENCE")] * virulence
    rows += [("STRESS", "METAL")] * stress + [("AMR", "POINT")]
    (root / "prot/protein/meta.tsv").write_text(
        "seqid\tgene\telement_type\telement_subtype\n"
        + "".join(f"afp_{index:05}\tgene{index}\t{kind}\t{subtype}\n"
                  for index, (kind, subtype) in enumerate(rows)))
    (root / "mutation/dna").mkdir(parents=True)
    for name in organisms:
        (root / "mutation/dna" / f"{name}.fna").write_text(">locus\nACGT\n")
    (root / "manifest.json").write_text(json.dumps({"hydra_db_version": 1, "databases": {
        "ncbi": {"path": separator.join(("nucl", "ncbi")), "version": release, "kind": "nucl",
                 "installed": "2026-09-12 18:41:37"},
        "protein": {"path": separator.join(("prot", "protein")), "version": release, "kind": "prot",
                    "installed": "2026-09-12 18:41:33", "organisms": list(organisms)}}}))
    return root


def capabilities(tools=None, version=None):
    return {"assembly_available": True, "hydra_version": version or runtime.HYDRA_VERSION,
            "tools": tools if tools is not None else {name: f"/tools/{name}"
                                                      for name in runtime.TOOLS},
            "databases": {}, "limitations": ["Assembly evidence only."]}


@pytest.fixture
def ready(monkeypatch, tmp_path):
    root = database(tmp_path / "reference")
    monkeypatch.setattr(runtime, "_ncbi_version", lambda: "2026-08-07.1")
    monkeypatch.setattr(runtime, "runtime_capabilities", lambda db_root=None: {
        "assembly_available": True, "tools": {name: f"/tools/{name}" for name in runtime.TOOLS},
        "hydra_version": runtime.HYDRA_VERSION,
        "databases": runtime.installed_databases(db_root) if db_root else {}, "limitations": ["Assembly evidence only."]})
    return root


def test_capabilities_and_missing_store_are_read_only(tmp_path):
    root = tmp_path / "absent"
    info = runtime.runtime_capabilities(root)
    assert "available" in info and "message" in info and "blastn" in info
    assert info["reads_available"] is False
    assert info["databases"] == {}
    assert not root.exists()


def test_manifest_rejects_external_paths_and_invalid_shapes(tmp_path):
    for entry in [{"path": "../outside"}, "not a dictionary", {"path": 1}]:
        (tmp_path / "manifest.json").write_text(json.dumps({"databases": {"ncbi": entry}}))
        with pytest.raises(runtime.HydraRuntimeError):
            runtime.installed_databases(tmp_path)


def test_assembly_runs_real_upstream_argument_contract_and_preserves_inputs(ready, tmp_path, monkeypatch):
    assembly = tmp_path / "isolate space é.fasta"
    assembly.write_text(">contig\nACGTACGT\n")
    original = assembly.read_bytes()
    captured = []

    def child(arguments, directory, cancelled=None, progress=None):
        captured.extend(arguments)
        output = Path(arguments[arguments.index("--outdir") + 1])
        output.mkdir()
        (output / "hydra.json").write_text(json.dumps({
            "hydra_version": "1.4.0", "command": "hydra run", "databases": ["ncbi"],
            "parameters": {}, "samples": [{"sample": "sample_a", "input_type": "assembly",
                                             "inputs": [str(assembly)], "hits": []}]}))
        return ["worker", *arguments], "Completed original HYDRA pipeline"

    monkeypatch.setattr(runtime, "_run_child", child)
    result = runtime.run_assemblies([assembly], ready, ["ncbi"], sample_names=["sample_a"],
                                    work_root=tmp_path / "jobs")
    assert captured[0] == "run"
    assert "--no-mlst" in captured and "--no-typing" in captured
    assert "--no-auto-organism" in captured and "--no-heteroresistance" in captured
    assert captured[captured.index("--assembly") + 1] == str(assembly)
    assert captured[captured.index("--name") + 1] == "sample_a"
    provenance = result["execution_provenance"]
    assert provenance["source_commit"] == runtime.HYDRA_COMMIT
    assert provenance["inputs"][0]["sha256"] == runtime.file_sha256(assembly)
    assert "ncbi" in provenance["reference_snapshot"]["databases"]
    assert assembly.read_bytes() == original
    assert list((tmp_path / "jobs").iterdir()) == []


def test_explicit_organism_and_options_are_forwarded_without_guessing(ready, tmp_path, monkeypatch):
    assembly = tmp_path / "isolate.fasta"
    assembly.write_text(">contig\nACGT\n")
    captured = []

    def child(arguments, directory, *args):
        captured.extend(arguments)
        raise runtime.HydraRuntimeError("stopped after capturing argv")

    monkeypatch.setattr(runtime, "_run_child", child)
    with pytest.raises(runtime.HydraRuntimeError, match="capturing"):
        runtime.run_assemblies([assembly], ready, ["ncbi"], organism="Staphylococcus_aureus",
                               protein=False, point_mutations=False)
    assert captured[captured.index("--organism") + 1] == "Staphylococcus_aureus"
    assert "--no-auto-organism" not in captured
    assert "--no-protein" in captured and "--no-point-mutations" in captured


def test_reads_and_duplicate_inputs_are_not_silently_analyzed(ready, tmp_path):
    reads = tmp_path / "reads.fasta"
    reads.write_text("@read\nACGT\n+\nIIII\n")
    with pytest.raises(runtime.HydraRuntimeError, match="FASTA assemblies only"):
        runtime.run_assemblies([reads], ready)
    with pytest.raises(runtime.HydraRuntimeError, match="distinct"):
        runtime.run_assemblies([reads, reads], ready)
    with pytest.raises(runtime.HydraRuntimeError, match="Threads"):
        runtime.run_assemblies([reads], ready, threads=0)


@pytest.mark.parametrize("value", [-1, 101, float("nan"), float("inf"), True, "90"])
def test_invalid_scientific_thresholds_are_rejected(ready, tmp_path, value):
    assembly = tmp_path / "isolate.fa"
    assembly.write_text(">contig\nACGT\n")
    with pytest.raises(runtime.HydraRuntimeError, match="percentage"):
        runtime.run_assemblies([assembly], ready, min_identity=value)


def test_missing_database_does_not_trigger_network_or_child(ready, tmp_path, monkeypatch):
    assembly = tmp_path / "isolate.fa"
    assembly.write_text(">contig\nACGT\n")
    monkeypatch.setattr(runtime, "_run_child", lambda *args: pytest.fail("must not launch"))
    with pytest.raises(runtime.HydraRuntimeError, match="not installed"):
        runtime.run_assemblies([assembly], ready, ["card"])


def test_database_publication_is_atomic_and_preserves_old_snapshot(ready, tmp_path, monkeypatch):
    old_manifest = (ready / "manifest.json").read_bytes()

    def child(arguments, directory, *args):
        assert arguments[:2] == ["db", "download"]
        assert "--force" not in arguments
        destination = Path(arguments[arguments.index("--db-dir") + 1])
        assert destination != ready
        database(destination)
        return arguments, "Download complete"

    monkeypatch.setattr(runtime, "_run_child", child)
    result = runtime.update_databases(ready, ["ncbi"])
    published = Path(result["database_root"])
    assert published.is_dir() and published != ready
    assert published.parent == ready.parent
    assert (ready / "manifest.json").read_bytes() == old_manifest
    assert result["reference_snapshot"]["root"] == str(published)
    assert not list(tmp_path.glob(".wmlstudio-hydra-stage-*"))


def test_cancelled_database_build_publishes_nothing_and_preserves_old(ready, tmp_path, monkeypatch):
    original = sorted(tmp_path.iterdir())

    def child(arguments, directory, *args):
        database(Path(arguments[arguments.index("--db-dir") + 1]))
        raise AnalysisCancelled("User cancelled")

    monkeypatch.setattr(runtime, "_run_child", child)
    with pytest.raises(AnalysisCancelled):
        runtime.update_databases(ready, ["ncbi"])
    assert sorted(tmp_path.iterdir()) == original
    assert runtime.installed_databases(ready)["ncbi"]["version"] == "test-snapshot"


def test_child_failure_keeps_diagnostics_and_uses_argument_array(tmp_path, monkeypatch):
    monkeypatch.setattr(runtime, "_worker_command", lambda args: [
        sys.executable, "-c", "import sys; print('native search diagnostic'); sys.exit(7)"])
    with pytest.raises(runtime.HydraRuntimeError, match="native search diagnostic"):
        runtime._run_child([], tmp_path)


def test_child_cancellation_terminates_owned_process(tmp_path, monkeypatch):
    monkeypatch.setattr(runtime, "_worker_command", lambda args: [
        sys.executable, "-c", "import time; time.sleep(60)"])
    started = time.monotonic()
    with pytest.raises(AnalysisCancelled):
        runtime._run_child([], tmp_path, cancelled=lambda: time.monotonic() - started > 0.2)
    assert time.monotonic() - started < 15


def test_a_run_with_no_reference_data_is_refused_before_anything_is_opened(tmp_path, monkeypatch):
    """The silence the user reported as "HYDRA does not work" must become a sentence."""
    assembly = tmp_path / "isolate.fasta"
    assembly.write_text(">contig\nACGTACGT\n")
    empty = tmp_path / "store"
    empty.mkdir()
    monkeypatch.setattr(runtime, "runtime_capabilities", lambda db_root=None: capabilities())
    monkeypatch.setattr(runtime, "_run_child", lambda *args, **options: pytest.fail("must not launch"))
    with pytest.raises(runtime.HydraRuntimeError) as failure:
        runtime.run_assemblies([assembly], empty)
    message = str(failure.value)
    assert "HYDRA cannot start" in message and "nothing was run" in message
    assert "every gene and mutation this screen can name" in message
    assert str(empty) in message and "no reference database" in message
    assert "Nothing is ever downloaded during a run." in message


def test_a_missing_blast_tool_names_the_tool_and_the_search_it_performs(tmp_path, monkeypatch):
    assembly = tmp_path / "isolate.fasta"
    assembly.write_text(">contig\nACGTACGT\n")
    store = complete_database(tmp_path / "store")
    tools = {name: (None if name == "blastx" else f"/tools/{name}") for name in runtime.TOOLS}
    monkeypatch.setattr(runtime, "runtime_capabilities",
                        lambda db_root=None: capabilities(tools=tools))
    monkeypatch.setattr(runtime, "_run_child", lambda *args, **options: pytest.fail("must not launch"))
    with pytest.raises(runtime.HydraRuntimeError) as failure:
        runtime.run_assemblies([assembly], store)
    message = str(failure.value)
    assert "blastx" in message and "point mutations are read" in message
    assert "beside the application or on PATH" in message
    # The tools that are present are not listed as problems.
    assert "blastn —" not in message


def test_an_unpinned_engine_is_refused_by_name_rather_than_run_anyway(tmp_path, monkeypatch):
    store = complete_database(tmp_path / "store")
    monkeypatch.setattr(runtime, "runtime_capabilities",
                        lambda db_root=None: capabilities(version="0.9.0"))
    checked = runtime.preflight(store)
    assert checked["ready"] is False
    assert f"HYDRA {runtime.HYDRA_VERSION}" in checked["message"]
    assert "0.9.0 is installed" in checked["message"]


def test_an_absent_reference_set_is_a_recorded_warning_not_a_silent_gap(tmp_path, monkeypatch):
    """A store with no protein set still runs; what it could not look for is stated."""
    monkeypatch.setattr(runtime, "runtime_capabilities", lambda db_root=None: capabilities())
    checked = runtime.preflight(database(tmp_path / "partial"))
    assert checked["ready"] is True and checked["databases"] == ["ncbi"]
    warning = next(text for text in checked["warnings"] if "'protein'" in text)
    assert "point-mutation catalogue" in warning
    assert "that is not a negative result" in warning


def test_a_store_staged_on_one_platform_is_readable_on_the_other(tmp_path):
    """A Windows-staged manifest records nucl\\ncbi; the same snapshot must still read."""
    store = complete_database(tmp_path / "windows", separator="\\")
    assert sorted(runtime.installed_databases(store)) == ["ncbi", "protein"]
    assert runtime.database_status(store)["release"] == "2026-08-07.1"
    # Normalising separators must not weaken the containment check.
    (store / "manifest.json").write_text(json.dumps(
        {"databases": {"ncbi": {"path": "..\\outside"}}}))
    with pytest.raises(runtime.HydraRuntimeError, match="leaves the selected store"):
        runtime.installed_databases(store)


def test_the_installed_release_is_labelled_with_its_version_date_and_age(tmp_path):
    fresh = runtime.database_status(complete_database(tmp_path / "fresh"))
    assert fresh["release"] == "2026-08-07.1" and fresh["staged"].startswith("2026-09-12")
    assert fresh["age_days"] is not None and fresh["stale"] is False
    assert "NCBI AMRFinderPlus reference release 2026-08-07.1" in fresh["label"]
    assert "ncbi, protein" in fresh["label"] and "days old" in fresh["label"]

    old = runtime.database_status(complete_database(tmp_path / "old", release="2019-01-01.1"))
    assert old["stale"] is True and old["age_days"] > runtime.REFERENCE_AGE_DAYS

    absent = runtime.database_status(tmp_path / "nothing-here")
    assert absent["installed"] == {} and absent["age_days"] is None
    assert absent["label"] == "No AMR reference database is installed."


def test_the_organism_catalogue_separates_accepted_names_from_curated_ones(tmp_path):
    catalogue = runtime.organism_catalogue(complete_database(tmp_path / "store"))
    # Pseudomonas is an accepted taxgroup with no DNA catalogue staged for it.
    assert "Pseudomonas_aeruginosa" in catalogue["accepted"]
    assert "Pseudomonas_aeruginosa" not in catalogue["point_mutations"]
    assert "Staphylococcus_aureus" in catalogue["point_mutations"]


@pytest.mark.parametrize("assigned,expected", [
    ("Staphylococcus aureus", "Staphylococcus_aureus"),
    ("Klebsiella pneumoniae", "Klebsiella_pneumoniae"),
    # Upstream groups every Escherichia at genus level; the species resolves to it.
    ("Escherichia coli", "Escherichia"),
    ("Listeria monocytogenes", None),
    ("", None),
])
def test_an_assigned_organism_is_resolved_against_the_installed_catalogue(tmp_path, assigned,
                                                                          expected):
    store = complete_database(tmp_path / "store")
    assert runtime.match_organism(assigned, store) == expected


def test_a_spaced_organism_reaches_the_engine_in_the_form_it_accepts(tmp_path, monkeypatch):
    """The assignment reads "Staphylococcus aureus"; the engine only takes the taxgroup."""
    store = complete_database(tmp_path / "store")
    assembly = tmp_path / "isolate.fasta"
    assembly.write_text(">contig\nACGTACGT\n")
    monkeypatch.setattr(runtime, "runtime_capabilities", lambda db_root=None: capabilities())
    captured = []

    def child(arguments, directory, cancelled=None, progress=None):
        captured.extend(arguments)
        output = Path(arguments[arguments.index("--outdir") + 1])
        output.mkdir()
        (output / "hydra.json").write_text(json.dumps({
            "hydra_version": "1.4.0", "command": "hydra run", "databases": ["ncbi", "protein"],
            "parameters": {}, "samples": [{"sample": "iso", "hits": []}]}))
        return ["worker", *arguments], "done"

    monkeypatch.setattr(runtime, "_run_child", child)
    report = runtime.run_assemblies([assembly], store, sample_names=["iso"],
                                    organism="Staphylococcus aureus")
    assert captured[captured.index("--organism") + 1] == "Staphylococcus_aureus"
    organism = report["execution_provenance"]["organism"]
    assert organism["requested"] == "Staphylococcus aureus"
    assert organism["resolved"] == "Staphylococcus_aureus"
    assert report["execution_provenance"]["reference_release"]["release"] == "2026-08-07.1"


def test_an_organism_with_no_catalogue_is_reported_not_replaced(tmp_path, monkeypatch):
    store = complete_database(tmp_path / "store")
    assembly = tmp_path / "isolate.fasta"
    assembly.write_text(">contig\nACGTACGT\n")
    monkeypatch.setattr(runtime, "runtime_capabilities", lambda db_root=None: capabilities())
    captured = []

    def child(arguments, directory, cancelled=None, progress=None):
        captured.extend(arguments)
        output = Path(arguments[arguments.index("--outdir") + 1])
        output.mkdir()
        (output / "hydra.json").write_text(json.dumps({
            "hydra_version": "1.4.0", "command": "hydra run", "databases": ["ncbi"],
            "parameters": {}, "samples": [{"sample": "iso", "hits": []}]}))
        return ["worker", *arguments], "done"

    monkeypatch.setattr(runtime, "_run_child", child)
    report = runtime.run_assemblies([assembly], store, sample_names=["iso"],
                                    organism="Listeria monocytogenes")
    # Not passed at all, and no neighbouring organism's catalogue substituted.
    assert "--organism" not in captured and "--no-auto-organism" in captured
    warning = next(text for text in report["import_warnings"] if "Listeria" in text)
    assert "not an absence of mutations" in warning
    assert report["execution_provenance"]["organism"]["resolved"] == ""


def test_a_store_that_publishes_no_catalogue_does_not_second_guess_the_caller(tmp_path):
    """With nothing to check a name against, the caller's organism is passed through."""
    store = database(tmp_path / "reference")
    assert runtime.organism_catalogue(store)["accepted"] == []
    checked = runtime.preflight(store, organism="Staphylococcus_aureus",
                                capabilities=capabilities())
    assert checked["organism"]["resolved"] == "Staphylococcus_aureus"


def test_a_reference_release_older_than_the_threshold_says_so_before_the_run(tmp_path,
                                                                             monkeypatch):
    monkeypatch.setattr(runtime, "runtime_capabilities", lambda db_root=None: capabilities())
    store = complete_database(tmp_path / "old", release="2019-01-01.1")
    checked = runtime.preflight(store, organism="Staphylococcus aureus")
    assert checked["ready"] is True
    stale = next(text for text in checked["warnings"] if "days old" in text)
    assert "Determinants named after it are not in it" in stale


def test_the_database_page_lets_a_user_see_and_choose_from_the_whole_catalogue(tmp_path, qtbot):
    """The user's complaint verbatim: HYDRA needs inputs it never names."""
    from wmlstudio.amr_databases import COLUMNS, AMRDatabaseDialog

    store = complete_database(tmp_path / "store")
    dialog = AMRDatabaseDialog(store)
    qtbot.addWidget(dialog)
    names = [dialog.table.item(row, 0).text() for row in range(dialog.table.rowCount())]
    assert set(names) == {entry["name"] for entry in dialog.catalogue["entries"]}
    assert len(names) > 2, "the page must list more than what happens to be installed"
    rows = {dialog.table.item(row, 0).text(): row for row in range(dialog.table.rowCount())}
    state, purpose = COLUMNS.index("Installed"), COLUMNS.index("What it is searched for")
    assert "Installed" in dialog.table.item(rows["ncbi"], state).text()
    assert "Not installed" in dialog.table.item(rows["card"], state).text()
    assert "nothing is reported for it" in dialog.table.item(rows["card"], state).text()
    assert dialog.table.item(rows["card"], purpose).text()

    # Installed sets are ticked, so the button means "update these"; a set the
    # engine cannot fetch cannot be ticked at all.
    assert dialog.selected_names() == ["ncbi", "protein"]
    manual = next(entry["name"] for entry in dialog.catalogue["entries"]
                  if entry["download"] == "by hand")
    assert not (dialog.table.item(rows[manual], 0).flags()
                & Qt.ItemFlag.ItemIsUserCheckable)

    # A tick survives re-reading the store from disk.
    dialog.table.item(rows["card"], 0).setCheckState(Qt.CheckState.Checked)
    assert "card" in dialog.selected_names()
    dialog.refresh()
    assert "card" in dialog.selected_names()
    # And the confirmation names the licence before anything is downloaded.
    text = dialog.confirmation(dialog.selected_names())
    assert "CARD academic licence" in text and "read them before installing" in text
    assert "samples are not uploaded" in text
    assert dialog.organism_sentence(runtime.database_status(store)).count("catalogue") >= 2
    assert "no other organism's catalogue is ever substituted" in dialog.organism_sentence(
        runtime.database_status(store))


def test_the_database_surface_says_what_an_empty_store_would_report(tmp_path):
    """A table with no rows reads as "nothing found"; the gap must be said in words."""
    from wmlstudio.amr_databases import AMRDatabaseDialog

    nothing = AMRDatabaseDialog.gap_sentence(runtime.database_status(tmp_path / "absent"))
    assert "would have nothing to search" in nothing
    assert "That is not a negative result." in nothing

    partial = AMRDatabaseDialog.gap_sentence(runtime.database_status(database(tmp_path / "part")))
    assert "Missing: protein" in partial and "point-mutation catalogue" in partial

    complete = AMRDatabaseDialog.gap_sentence(
        runtime.database_status(complete_database(tmp_path / "full")))
    assert "Missing:" not in complete and "screened for genes only" in complete

    old = AMRDatabaseDialog.gap_sentence(
        runtime.database_status(complete_database(tmp_path / "old", release="2019-01-01.1")))
    assert "days old" in old and "Determinants named after it are not in it" in old


def test_every_reference_set_the_engine_can_search_is_named_whether_or_not_it_is_here(tmp_path):
    """HYDRA's silence about its inputs is the complaint; the whole list answers it."""
    catalogue = runtime.database_catalogue(complete_database(tmp_path / "store"))
    rows = {entry["name"]: entry for entry in catalogue["entries"]}
    # Every set the engine knows, not only the two that are installed.
    assert set(catalogue["installed"]) == {"ncbi", "protein"}
    assert len(rows) > len(catalogue["installed"])
    assert {"card", "vfdb", "plasmidfinder"} <= set(catalogue["available"])
    for entry in catalogue["entries"]:
        assert entry["purpose"] and entry["title"] and entry["licence"]
        assert entry["state"] in {"installed", "not installed"}
        assert entry["download"] in {"automatic", "by hand"}
        # A size is either measured, estimated with its basis named, or admitted absent.
        assert entry["size_basis"]
    assert rows["ncbi"]["installed"] and rows["ncbi"]["version"] == "2026-08-07.1"
    assert "measured" in rows["ncbi"]["size_basis"] and rows["ncbi"]["size"].endswith("B")
    assert not rows["card"]["installed"] and rows["card"]["size"] == ""
    # A licence nobody can act on without reading it is never presented as open.
    assert rows["ncbi"]["open_licence"] is True and rows["ncbi"]["licence_note"] == ""
    assert rows["card"]["open_licence"] is False and "terms apply" in rows["card"]["licence_note"]
    assert json.loads(json.dumps(catalogue)) == catalogue


def test_a_set_the_engine_cannot_fetch_is_listed_rather_than_promised(tmp_path):
    catalogue = runtime.database_catalogue(complete_database(tmp_path / "store"))
    rows = {entry["name"]: entry for entry in catalogue["entries"]}
    manual = [name for name, entry in rows.items() if entry["download"] == "by hand"]
    assert manual and not set(manual) & set(catalogue["automatic"])
    # It is still named, still says what it is for, and still says who publishes it.
    assert all(rows[name]["purpose"] and rows[name]["url"] for name in manual)


def test_a_reference_set_name_the_engine_does_not_know_is_refused_before_any_download(
        ready, monkeypatch):
    monkeypatch.setattr(runtime, "_run_child", lambda *args: pytest.fail("must not launch"))
    with pytest.raises(runtime.HydraRuntimeError) as failure:
        runtime.update_databases(ready, ["ncbi", "carrd"])
    assert "no reference set called carrd" in str(failure.value)
    assert "card" in str(failure.value)


def test_a_failed_download_says_what_failed_and_leaves_the_old_snapshot_alone(ready, tmp_path,
                                                                              monkeypatch):
    before = sorted(tmp_path.iterdir())
    manifest = (ready / "manifest.json").read_bytes()

    def child(arguments, directory, *args):
        raise runtime.HydraRuntimeError("HYDRA exited with code 1.\nprovider unreachable")

    monkeypatch.setattr(runtime, "_run_child", child)
    with pytest.raises(runtime.HydraRuntimeError) as failure:
        runtime.update_databases(ready, ["ncbi", "protein"])
    message = str(failure.value)
    assert "Downloading protein, ncbi failed" in message
    assert "reference store you are using is unchanged" in message
    assert "provider unreachable" in message
    assert sorted(tmp_path.iterdir()) == before
    assert (ready / "manifest.json").read_bytes() == manifest


def test_dna_and_protein_point_mutation_catalogues_are_counted_separately(tmp_path):
    """Two different catalogues answer two different questions; one number hides that."""
    store = complete_database(tmp_path / "store",
                              organisms=("Escherichia", "Staphylococcus_aureus"),
                              protein_mutations=("Escherichia", "Staphylococcus_aureus",
                                                 "Pseudomonas_aeruginosa"))
    catalogue = runtime.organism_catalogue(store)
    assert catalogue["dna_point_mutations"] == ["Escherichia", "Staphylococcus_aureus"]
    assert "Pseudomonas_aeruginosa" in catalogue["protein_point_mutations"]
    assert "Pseudomonas_aeruginosa" not in catalogue["dna_point_mutations"]
    assert "Pseudomonas_aeruginosa" in catalogue["point_mutations"]
    status = runtime.database_status(store)
    assert status["dna_point_mutation_organisms"] == catalogue["dna_point_mutations"]
    assert status["protein_point_mutation_organisms"] == catalogue["protein_point_mutations"]


def test_a_genus_level_catalogue_is_recognised_as_covering_its_species(tmp_path, monkeypatch):
    """The engine applies a parent taxgroup to a child; saying otherwise understates cover."""
    monkeypatch.setattr(runtime, "runtime_capabilities", lambda db_root=None: capabilities())
    store = complete_database(tmp_path / "store",
                              organisms=("Klebsiella_pneumoniae",),
                              protein_mutations=("Klebsiella",))
    assert runtime.catalogue_covers("Klebsiella_pneumoniae", ["Klebsiella"]) is True
    assert runtime.catalogue_covers("Klebsiella_pneumoniae", ["Klebsiella_oxytoca"]) is False
    checked = runtime.preflight(store, organism="Klebsiella pneumoniae")
    assert checked["organism"]["resolved"] == "Klebsiella_pneumoniae"
    assert checked["organism"]["point_mutation_level"] == "dna_and_protein"
    assert checked["organism"]["reason"] == ""


@pytest.mark.parametrize("organism,level,phrase", [
    ("Staphylococcus aureus", "dna_and_protein", ""),
    ("Pseudomonas aeruginosa", "protein_only", "no DNA catalogue"),
    ("Listeria monocytogenes", "none", "not an absence of mutations"),
])
def test_what_a_point_mutation_catalogue_covers_is_stated_per_isolate(tmp_path, monkeypatch,
                                                                      organism, level, phrase):
    monkeypatch.setattr(runtime, "runtime_capabilities", lambda db_root=None: capabilities())
    store = complete_database(tmp_path / "store",
                              organisms=("Escherichia", "Staphylococcus_aureus"),
                              protein_mutations=("Escherichia", "Staphylococcus_aureus",
                                                 "Pseudomonas_aeruginosa"))
    checked = runtime.preflight(store, organism=organism)
    assert checked["organism"]["point_mutation_level"] == level
    if phrase:
        assert phrase in checked["organism"]["reason"]
        assert checked["organism"]["reason"] in checked["warnings"]
    else:
        assert checked["organism"]["reason"] == ""


def test_virulence_is_searched_for_an_established_organism_and_not_for_an_unknown_one(
        tmp_path, monkeypatch):
    """The organism decides whose curation applies, so an unknown one gets neither."""
    monkeypatch.setattr(runtime, "runtime_capabilities", lambda db_root=None: capabilities())
    store = complete_database(tmp_path / "store")
    known = runtime.preflight(store, organism="Staphylococcus aureus")
    assert known["virulence"]["enabled"] is True and known["virulence"]["requested"] == "auto"
    assert known["virulence"]["organism_curated"] is True
    assert known["virulence"]["virulence_elements"] == 3
    assert known["virulence"]["stress_elements"] == 2
    assert "not a validated virulence prediction" in known["virulence"]["reason"]

    unknown = runtime.preflight(store, organism="Listeria monocytogenes")
    assert unknown["virulence"]["enabled"] is False
    assert unknown["virulence"]["available"] is True
    assert "no organism curation applied" in unknown["virulence"]["reason"]
    assert unknown["virulence"]["reason"] in unknown["warnings"]


def test_virulence_asked_for_without_an_organism_runs_uncurated_and_says_so(tmp_path, monkeypatch):
    monkeypatch.setattr(runtime, "runtime_capabilities", lambda db_root=None: capabilities())
    store = complete_database(tmp_path / "store")
    forced = runtime.preflight(store, virulence=True)
    assert forced["virulence"]["enabled"] is True
    assert forced["virulence"]["organism_curated"] is False
    assert forced["virulence"]["reason"] in forced["warnings"]

    off = runtime.preflight(store, organism="Staphylococcus aureus", virulence=False)
    assert off["virulence"]["enabled"] is False and off["virulence"]["requested"] is False
    # Turning it off narrows the protein search only: the nucleotide catalogues
    # carry virulence-typed genes and report them either way, so the sentence must
    # not claim that nothing virulent was looked for.
    assert "governs the protein search only" in off["virulence"]["reason"]
    assert "limited to acquired resistance" in off["virulence"]["reason"]
    assert off["virulence"]["reason"] not in off["warnings"]


def test_a_store_with_no_protein_set_cannot_offer_virulence_at_all(tmp_path, monkeypatch):
    monkeypatch.setattr(runtime, "runtime_capabilities", lambda db_root=None: capabilities())
    checked = runtime.preflight(database(tmp_path / "genes-only"), virulence=True)
    assert checked["virulence"]["available"] is False
    assert checked["virulence"]["enabled"] is False
    assert "not among the databases being searched" in checked["virulence"]["reason"]
    assert runtime.element_counts(tmp_path / "genes-only")["read"] is False


def test_the_recorded_command_says_which_of_the_two_searches_was_run(tmp_path, monkeypatch):
    """--plus and --no-plus are different searches; a report must not be ambiguous."""
    store = complete_database(tmp_path / "store")
    assembly = tmp_path / "isolate.fasta"
    assembly.write_text(">contig\nACGTACGT\n")
    monkeypatch.setattr(runtime, "runtime_capabilities", lambda db_root=None: capabilities())
    captured = []

    def child(arguments, directory, cancelled=None, progress=None):
        captured.append(list(arguments))
        output = Path(arguments[arguments.index("--outdir") + 1])
        output.mkdir()
        (output / "hydra.json").write_text(json.dumps({
            "hydra_version": "1.4.0", "command": "hydra run", "databases": ["ncbi", "protein"],
            "parameters": {}, "samples": [{"sample": "iso", "hits": []}]}))
        return ["worker", *arguments], "done"

    monkeypatch.setattr(runtime, "_run_child", child)
    report = runtime.run_assemblies([assembly], store, sample_names=["iso"],
                                    organism="Staphylococcus aureus")
    assert "--plus" in captured[-1] and "--no-plus" not in captured[-1]
    assert report["execution_provenance"]["virulence"]["enabled"] is True

    unknown = runtime.run_assemblies([assembly], store, sample_names=["iso"],
                                     organism="Listeria monocytogenes")
    assert "--no-plus" in captured[-1] and "--plus" not in captured[-1]
    assert unknown["execution_provenance"]["virulence"]["enabled"] is False
    assert any("no organism curation applied" in text
               for text in unknown["import_warnings"])


def test_a_named_database_that_is_absent_stops_the_run_rather_than_quietly_narrowing_it(
        tmp_path, monkeypatch):
    """Asking for ncbi and protein and getting only ncbi answers a different question."""
    monkeypatch.setattr(runtime, "runtime_capabilities", lambda db_root=None: capabilities())
    partial = database(tmp_path / "partial")
    asked = runtime.preflight(partial, databases=["ncbi", "protein"])
    assert asked["ready"] is False and "protein" in asked["message"]
    # The same store, with nothing named, runs on what it has and says what it lost.
    defaulted = runtime.preflight(partial)
    assert defaulted["ready"] is True and defaulted["databases"] == ["ncbi"]
