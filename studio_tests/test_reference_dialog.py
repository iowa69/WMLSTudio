import json
import threading
import time

from PySide6.QtCore import QThread
from PySide6.QtTest import QSignalSpy
from PySide6.QtWidgets import QMessageBox

from wmlstudio import cgmlst_schemes, reference_dialog
from wmlstudio.reference_dialog import ReferenceManagerDialog
from wmlstudio.sequence import check_cancelled

KLEBSIELLA = "cgmlst.org:kpneumoniae-2358"
# The two organisms for which the catalogue pins more than one target set: a core
# set with an accessory set beside it, and a family of four including a pan-genome
# set. They are what makes "cgMLST or cgMLST plus accessory genes" a real choice.
ANTHRACIS_CORE = "pubmlst:banthracis-cgmlst-3803"
ANTHRACIS_ACCESSORY = "pubmlst:banthracis-accessory-1263"
GONOCOCCUS_PANGENOME = "pubmlst:ngonorrhoeae-pgmlst-1907"


def install_scheme(folder, *, targets=40, organism="Klebsiella pneumoniae", name="Kp cgMLST",
                   api="https://www.cgmlst.org/ncs/schema/Kpneumoniae_complex/", source="cgMLST.org"):
    """A scheme folder shaped exactly like a finished download: alleles and scheme.json."""
    folder.mkdir(parents=True, exist_ok=True)
    for index in range(targets):
        (folder / f"target{index:04d}.fasta").write_text(f">target{index:04d}_1\nACGT\n")
    (folder / "scheme.json").write_text(json.dumps({
        "name": name, "organism": organism, "description": name, "type": "cgMLST",
        "source": source, "API": api, "locus_count": targets, "last_updated": "2026-09-14"}))
    return folder


class Catalog:
    """A public catalog stand-in: no sockets, every call recorded with its thread."""

    TERMS_URL = "https://example.invalid/terms"
    TERMS_NOTICE = "Stand-in terms for a provider that does not grant redistribution."

    def __init__(self, schemes=None, installs=None):
        self.calls = []
        self.schemes = schemes if schemes is not None else [
            {"id": "PubMLST:example:1", "name": "Example MLST", "organism": "Examplegenus species",
             "locus_count": 7, "type": "MLST", "last_updated": "2026-01-01",
             "url": "https://rest.pubmlst.org/db/example/schemes/1",
             "access_notice": "Anonymous access excludes recent alleles."}]
        self.installs = installs or {}

    def search_schemes(self, query, *, cancelled, progress, **filters):
        self.calls.append(("search", QThread.currentThread()))
        check_cancelled(cancelled)
        return {"schemes": list(self.schemes), "organisms_searched": 2,
                "errors": [{"organism": "Unavailable genus", "error": "offline"}]}

    def download_scheme(self, entry, root, *, cancelled, progress, **options):
        self.calls.append(("download", QThread.currentThread()))
        self.last_download = (dict(entry), str(root), dict(options))
        check_cancelled(cancelled)
        folder = self.installs.get(entry.get("id"))
        if folder is not None:
            install_scheme(folder)
            return {"path": str(folder), "created": True, "locus_count": 40, "notes": []}
        return {"path": str(root / "verified"), "created": True, "locus_count": 7, "notes": []}


def open_dialog(qtbot, root, catalog=None):
    dialog = ReferenceManagerDialog(root, catalog=catalog)
    qtbot.addWidget(dialog)
    return dialog


def rows_for(page):
    return {(row["organism"], row["count"], row["status"]) for row in page.rows}


def test_catalog_dialog_does_not_network_until_search_and_workers_stay_off_gui(qtbot, tmp_path, qapp):
    client = Catalog()
    dialog = open_dialog(qtbot, tmp_path, client)
    assert client.calls == []
    dialog.query.setText("Examplegenus")
    dialog.search()
    qtbot.waitUntil(lambda: bool(dialog.pages["mlst"].online) and dialog.download_button.isEnabled(),
                    timeout=3000)
    # A seven-locus scheme belongs to the classical library, so that is the tab it
    # lands in and the tab the dialog brings forward.
    assert dialog.page is dialog.pages["mlst"]
    assert dialog.table.rowCount() == 1
    assert "Anonymous access" in dialog.notice.toPlainText()
    assert "Unavailable genus: offline" in dialog.notice.toPlainText()
    assert dialog.download_button.isEnabled()
    installed = QSignalSpy(dialog.schemeInstalled)
    dialog.download()
    qtbot.waitUntil(lambda: installed.count() == 1 and not dialog.worker.isRunning(), timeout=3000)
    assert installed.at(0)[0] == str(tmp_path / "schemes" / "verified")
    assert [name for name, _ in client.calls] == ["search", "download"]
    assert all(thread is not qapp.thread() for _, thread in client.calls)
    dialog.close()


def test_dialog_close_cancels_active_network_work_without_publishing(qtbot, tmp_path, monkeypatch):
    entered = threading.Event()

    class BlockingCatalog(Catalog):
        def search_schemes(self, query, *, cancelled, progress, **filters):
            entered.set()
            while not cancelled():
                threading.Event().wait(.005)
            check_cancelled(cancelled)

    dialog = open_dialog(qtbot, tmp_path, BlockingCatalog())
    dialog.show()
    dialog.query.setText("Examplegenus")
    dialog.search()
    qtbot.waitUntil(entered.is_set, timeout=3000)
    installed = QSignalSpy(dialog.schemeInstalled)
    # Closing now asks before it stops a download; this test is about what happens
    # once that is answered, so answer it rather than waiting on a modal forever.
    monkeypatch.setattr(dialog, "confirm_stop", lambda: True)
    dialog.reject()
    qtbot.waitUntil(lambda: not dialog.worker.isRunning() and not dialog.isVisible(), timeout=3000)
    assert installed.count() == 0
    assert "cancelled" in dialog.status.text().lower()


def test_the_two_libraries_are_separate_surfaces_counted_in_their_own_units(qtbot, tmp_path):
    """A seven-locus distance and a 2,358-target distance never share one list."""
    install_scheme(tmp_path / "cgmlst" / "Klebsiella_pneumoniae__cgmlst_org_2358")
    dialog = open_dialog(qtbot, tmp_path)
    cgmlst, mlst = dialog.pages["cgmlst"], dialog.pages["mlst"]

    assert [cgmlst.table.horizontalHeaderItem(index).text() for index in range(3)] == \
        ["Organism", "Scheme", "Targets"]
    assert [mlst.table.horizontalHeaderItem(index).text() for index in range(3)] == \
        ["Organism", "Scheme", "Loci"]
    installed_cgmlst = [row for row in cgmlst.rows if row["state"] == "installed"]
    assert len(installed_cgmlst) == 1 and installed_cgmlst[0]["count"] == 40
    assert not [row for row in mlst.rows if row["count"] > 30], \
        "a gene-by-gene scheme is never offered where a classical one is expected"
    assert all(row["count"] <= 30 for row in mlst.rows)
    assert {row["state"] for row in mlst.rows} == {"installed"}
    assert dialog.tabs.tabText(0) == "cgMLST schemes", "the cgMLST library comes first"


def test_a_downloaded_cgmlst_scheme_is_listed_installed_under_a_readable_name(
        qtbot, tmp_path, monkeypatch):
    """The reported bug, seen from the interface: found after the download, not lost."""
    monkeypatch.setattr(reference_dialog.QMessageBox, "question",
                        lambda *args, **kwargs: QMessageBox.StandardButton.Yes)
    destination = tmp_path / "cgmlst" / "Klebsiella_pneumoniae__cgmlst_org_2358"
    entry = {"id": KLEBSIELLA, "name": "Kp cgMLST", "organism": "Klebsiella pneumoniae",
             "locus_count": 2358, "type": "cgMLST", "last_updated": None,
             "url": "https://www.cgmlst.org/ncs/schema/Kpneumoniae_complex/", "slug": "Kpneumoniae_complex",
             "access_notice": "", "terms_notice": "Confirm your use is permitted.",
             "terms_url": "https://www.cgmlst.org/serverpolicy.html"}
    client = Catalog(schemes=[entry], installs={KLEBSIELLA: destination})
    dialog = open_dialog(qtbot, tmp_path, client)
    dialog.show_page("cgmlst")
    # Searched for first, so a leftover filter is in force when the download lands.
    dialog.page.query.setText("Klebsiella pneumoniae/variicola")
    assert dialog.page.select_row_id("key:" + KLEBSIELLA)

    installed = QSignalSpy(dialog.schemeInstalled)
    dialog.download()
    qtbot.waitUntil(lambda: installed.count() == 1 and not dialog.worker.isRunning(), timeout=5000)

    assert dialog.page is dialog.pages["cgmlst"]
    row = dialog.page.selected_row()
    assert row is not None and row["path"] == str(destination)
    assert row["state"] == "installed" and row["count"] == 40
    assert "Klebsiella pneumoniae" in row["title"] and "40 targets" in row["title"]
    assert "cgMLST.org" in row["title"]
    assert str(destination) in dialog.status.text()
    assert "cgMLST library" in dialog.status.text()
    # And it is not offered in the classical library under any name.
    assert not [item for item in dialog.pages["mlst"].rows if item["path"] == str(destination)]


def test_a_cgmlst_scheme_downloaded_into_the_classical_folder_is_moved_and_named(qtbot, tmp_path):
    """An earlier release filed it between two seven-locus schemes under a digest name."""
    misfiled = install_scheme(tmp_path / "schemes" / "cgmlst_org_kpneumoniae_abcdef0123456789")
    dialog = open_dialog(qtbot, tmp_path)

    assert not misfiled.exists(), "the bytes were moved, not copied"
    # Into the labelled slot the pinned catalogue already describes, not beside it.
    moved = tmp_path / "cgmlst" / "Klebsiella_pneumoniae__cgmlst_org_2358"
    assert (moved / "scheme.json").is_file()
    assert "Moved 1 cgMLST scheme" in dialog.status.text()
    installed = [row for row in dialog.pages["cgmlst"].rows if row["state"] == "installed"]
    assert [row["path"] for row in installed] == [str(moved)]
    assert installed[0]["organism"] == "Klebsiella pneumoniae"
    assert not [row for row in dialog.pages["mlst"].rows if "kpneumoniae" in row["location"]]


def test_a_core_and_an_accessory_set_are_an_explicit_choice_and_one_set_is_stated_plainly(
        qtbot, tmp_path):
    """One set pinned must read as an answer, two as a choice a person can reach.

    An organism with a single catalogued set said so in words. An organism with a
    core set and an accessory set beside it must offer both, because running
    cgMLST and running cgMLST plus accessory genes are two different runs.
    """
    dialog = open_dialog(qtbot, tmp_path)
    page = dialog.pages["cgmlst"]
    page.query.setText("Klebsiella pneumoniae")
    page.table.selectRow(0)

    # Only a core set is pinned for this organism: that is said in words, and the
    # choice is not offered as an empty menu.
    assert page.variant.count() == 1 and not page.variant.isEnabled()
    assert "Only a core target set is catalogued" in page.variant_note.text()

    page.query.setText("Bacillus anthracis")
    assert page.select_row_id("key:" + ANTHRACIS_CORE)

    assert page.variant.isEnabled() and page.variant.count() == 2
    assert page.variant.itemText(0) == ("Core target set — 3803 targets — B. anthracis cgMLST · "
                                        "PubMLST (University of Oxford)")
    assert page.variant.itemText(1).startswith("Accessory target set — 1263 targets — ")
    assert page.variant.currentData() == ANTHRACIS_CORE
    assert [page.variant.itemData(index) for index in range(2)] == \
        [ANTHRACIS_CORE, ANTHRACIS_ACCESSORY]

    note = page.variant_note.text()
    assert "1 core and 1 accessory target set(s) are catalogued" in note
    assert "DIFFERENT QUANTITIES" in note
    assert "never shares a scale, an axis or a threshold" in note
    assert "not offered for a core-plus-accessory run" in note
    assert "not comparable at all" in note
    assert "Type every isolate of one comparison against the set you choose here." in note
    # And what each set is for, so the choice is made on what they measure.
    assert "expected in every isolate" in note
    assert "biology and not a failed call" in note


def test_choosing_the_accessory_set_in_the_menu_selects_that_scheme_and_not_the_core_one(
        qtbot, tmp_path):
    """The menu was inert: it listed the sets but choosing one changed nothing."""
    dialog = open_dialog(qtbot, tmp_path)
    page = dialog.pages["cgmlst"]
    page.query.setText("Bacillus anthracis")
    assert page.select_row_id("key:" + ANTHRACIS_CORE)

    page.choose_variant(page.variant.findData(ANTHRACIS_ACCESSORY))

    row = page.selected_row()
    assert row["key"] == ANTHRACIS_ACCESSORY and row["count"] == 1263
    assert row["target_set_detail"] == "accessory"


def test_every_catalogued_target_set_for_one_organism_is_offered_with_its_own_target_count(
        qtbot, tmp_path):
    """Four gonococcal sets: two cores, an accessory set and a pan-genome set.

    The reported bug was that the choice could not be reached at all. A menu that
    lists every set, each with the number of targets it holds, is the difference
    between choosing a run and discovering afterwards which one was made.
    """
    dialog = open_dialog(qtbot, tmp_path)
    page = dialog.pages["cgmlst"]
    page.query.setText("Neisseria gonorrhoeae")
    assert page.select_row_id("key:" + GONOCOCCUS_PANGENOME)

    assert page.variant.isEnabled()
    assert [page.variant.itemText(index) for index in range(page.variant.count())] == [
        "Core target set — 1430 targets — N. gonorrhoeae cgMLST v2 · PubMLST (University of Oxford)",
        "Core target set — 1649 targets — N. gonorrhoeae cgMLST v1.0 · "
        "PubMLST (University of Oxford)",
        "Accessory target set — 251 targets — N. gonorrhoeae agMLST v1.0 · "
        "PubMLST (University of Oxford)",
        "Whole-genome set (core + accessory) — 1907 targets — N. gonorrhoeae pgMLST v1.0 · "
        "PubMLST (University of Oxford)"]
    assert page.variant.currentData() == GONOCOCCUS_PANGENOME
    note = page.variant_note.text()
    assert "2 core, 1 accessory and 1 whole-genome target set(s) are catalogued" in note
    assert "Run it instead of a core scheme, never beside one" in note


def test_a_pan_genome_scheme_is_labelled_whole_genome_and_never_core(qtbot, tmp_path):
    """It read 'Not recorded' beside two core sets, which invited reading it as one.

    A 1,907-target pan-genome distance is not a cgMLST distance. The column that
    names the target set must say which set each row is, in the same words the
    catalogue uses, or the one row that must not be mistaken for a core scheme is
    the one row that carries no label.
    """
    dialog = open_dialog(qtbot, tmp_path)
    page = dialog.pages["cgmlst"]
    page.query.setText("Neisseria gonorrhoeae")

    column = page.COLUMNS["cgmlst"].index("Target set")
    labelled = {int(page.table.item(index, 2).text()): page.table.item(index, column).text()
                for index in range(page.table.rowCount())}
    assert labelled == {1430: "Core", 1649: "Core", 251: "Accessory", 1907: "Whole genome"}


def test_a_core_set_and_a_core_plus_accessory_set_stay_separate_all_the_way_down(qtbot, tmp_path):
    """Two isolates typed against different target sets are not comparable.

    So the two sets must stay distinguishable everywhere a result can be traced
    back to a scheme: a different catalogue key, a different install folder, a
    different pinned target list, and no cutoff carried across.
    """
    dialog = open_dialog(qtbot, tmp_path)
    page = dialog.pages["cgmlst"]
    page.query.setText("Bacillus anthracis")
    rows = {row["key"]: row for row in page.visible}

    core, accessory = rows[ANTHRACIS_CORE], rows[ANTHRACIS_ACCESSORY]
    assert core["row_id"] != accessory["row_id"]
    assert core["location"] != accessory["location"], "one folder could not hold both"
    assert (core["count"], accessory["count"]) == (3803, 1263)

    pins = {cgmlst_schemes.entry_for(key)["target_list_sha256"]
            for key in (ANTHRACIS_CORE, ANTHRACIS_ACCESSORY)}
    assert len(pins) == 2, "the two target lists are pinned apart, not by count alone"

    # The published five-allele cutoff was derived on the core set. It is not
    # offered for the accessory set, and the accessory row says why in so many words.
    assert cgmlst_schemes.threshold_for(ANTHRACIS_ACCESSORY)["threshold"] is None
    page.select_row_id(accessory["row_id"])
    detail = page.notice.toPlainText()
    assert "Target set: Accessory · 1263 targets." in detail
    assert "the core cutoff is not offered for it" in detail
    assert "never added together into one number" in detail
    assert "ACCESSORY target set, not a cgMLST scheme" in detail


def test_a_scheme_found_online_reports_its_target_set_as_unrecorded_rather_than_core(
        qtbot, tmp_path):
    """A search result says nothing about its target set, so neither may the column.

    Rendering an unlabelled set as 'Core' would turn "we do not know" into "this is
    the core genome", which is the one reading that lets a pan-genome scheme be
    compared against a core cutoff.
    """
    found = {"id": "PubMLST:ng:81", "name": "pgMLST v1.0", "organism": "Neisseria gonorrhoeae",
             "locus_count": 1907, "type": "cgMLST", "last_updated": "2026-09-15",
             "url": "https://rest.pubmlst.org/db/pubmlst_neisseria_seqdef/schemes/81",
             "access_notice": ""}
    dialog = open_dialog(qtbot, tmp_path, Catalog(schemes=[found]))
    dialog.query.setText("Neisseria gonorrhoeae")
    dialog.search()
    qtbot.waitUntil(lambda: bool(dialog.pages["cgmlst"].online), timeout=3000)

    page = dialog.pages["cgmlst"]
    online = [row for row in page.rows if row["state"] == "online"]
    assert [row["target_set"] for row in online] == [""]
    assert reference_dialog._target_set_words(online[0]) == "Not recorded"


def test_another_providers_scheme_for_the_same_organism_is_named_but_never_offered_as_a_variant(
        qtbot, tmp_path):
    dialog = open_dialog(qtbot, tmp_path)
    page = dialog.pages["cgmlst"]
    page.query.setText("Acinetobacter baumannii")
    counts = sorted(row["count"] for row in page.visible)
    page.select_row_id("key:cgmlst.org:abaumannii-2390")

    assert counts == [2133, 2390], "both catalogued schemes are listed"
    assert [page.variant.itemData(index) for index in range(page.variant.count())] == \
        ["cgmlst.org:abaumannii-2390"]
    note = page.variant_note.text()
    assert "1 other cgMLST scheme(s) are catalogued for this organism" in note
    assert "2133 targets" in note and "PubMLST" in note
    assert "never share a distance, an axis or a cutoff" in note


def test_the_library_summary_names_the_organisms_that_have_more_than_one_target_set(
        qtbot, tmp_path):
    """The choice was invisible until the right row happened to be selected.

    A person who does not already know that B. anthracis has an accessory set has
    no way of discovering that the greyed-out menu ever lights up.
    """
    dialog = open_dialog(qtbot, tmp_path)
    summary = dialog.pages["cgmlst"].summary.text()

    assert "More than one target set is catalogued for " in summary
    assert "Bacillus anthracis" in summary and "Neisseria gonorrhoeae" in summary
    # Two core schemes from one provider are two schemes, not two target sets.
    assert "Salmonella enterica" not in summary
    assert "the accessory or whole-genome set beside it" in summary
    assert "More than one target set" not in dialog.pages["mlst"].summary.text()


def test_an_installed_scheme_offers_its_own_providers_target_sets_and_no_one_elses(qtbot, tmp_path):
    """An installed scheme recorded no scheme family, so the menu fell back to the
    organism and offered a different provider's core scheme as if it were the
    accessory half of this one. Two providers' core schemes for one organism are
    two schemes, not two target sets of one."""
    install_scheme(tmp_path / "cgmlst" / "Acinetobacter_baumannii__pubmlst_2133",
                   organism="Acinetobacter baumannii", name="cgMLST v1", source="PubMLST",
                   api="https://rest.pubmlst.org/db/pubmlst_abaumannii_seqdef/schemes/3")
    dialog = open_dialog(qtbot, tmp_path)
    page = dialog.pages["cgmlst"]
    page.query.setText("Acinetobacter baumannii")
    installed = next(row for row in page.visible if row["state"] == "installed")
    page.select_row_id(installed["row_id"])

    assert installed["key"] == "pubmlst:abaumannii-cgmlst-2133"
    assert [page.variant.itemData(index) for index in range(page.variant.count())] == \
        ["pubmlst:abaumannii-cgmlst-2133"]
    assert not page.variant.isEnabled(), "there is one set, so there is nothing to choose"
    note = page.variant_note.text()
    assert "1 other cgMLST scheme(s) are catalogued for this organism" in note
    assert "2390 targets · Core · cgMLST.org" in note
    assert "never share a distance, an axis or a cutoff" in note


def test_an_installed_scheme_that_records_no_target_set_says_so_instead_of_claiming_core(
        qtbot, tmp_path):
    """A folder that never said which set it holds must not be read as a core one.

    "We did not look" rendering as "this is the core genome" is what would let a
    hand-copied accessory or pan-genome folder be compared against a core cutoff.
    """
    folder = tmp_path / "cgmlst" / "Homebrew__unknown_40"
    folder.mkdir(parents=True)
    for index in range(40):
        (folder / f"target{index:04d}.fasta").write_text(f">target{index:04d}_1\nACGT\n")
    (folder / "scheme.json").write_text(json.dumps(
        {"name": "homebrew", "organism": "Mystery organism", "locus_count": 40}))
    dialog = open_dialog(qtbot, tmp_path)
    page = dialog.pages["cgmlst"]
    page.select_path(folder)
    row = page.selected_row()

    assert row["target_set"] == "" and row["target_set_detail"] == ""
    column = page.COLUMNS["cgmlst"].index("Target set")
    assert page.table.item(page.table.currentRow(), column).text() == "Not recorded"
    assert "Target set: not recorded." in page.notice.toPlainText()
    assert "cannot be read as a cgMLST core distance" in page.notice.toPlainText()


def test_nothing_is_downloaded_until_the_providers_terms_are_confirmed(qtbot, tmp_path, monkeypatch):
    client = Catalog()
    dialog = open_dialog(qtbot, tmp_path, client)
    dialog.show_page("cgmlst")
    assert dialog.page.select_row_id("key:" + KLEBSIELLA)
    assert dialog.download_button.text() == "Download under the provider's terms…"

    monkeypatch.setattr(reference_dialog.QMessageBox, "question",
                        lambda *args, **kwargs: QMessageBox.StandardButton.No)
    dialog.download()
    assert client.calls == []
    assert "Nothing was downloaded" in dialog.status.text()

    monkeypatch.setattr(reference_dialog.QMessageBox, "question",
                        lambda *args, **kwargs: QMessageBox.StandardButton.Yes)
    dialog.download()
    qtbot.waitUntil(lambda: not dialog.worker.isRunning(), timeout=3000)
    entry, root, options = client.last_download
    assert [name for name, _ in client.calls] == ["download"]
    assert entry["slug"] == "Kpneumoniae_complex" and entry["locus_count"] == 2358
    assert entry["url"] == "https://www.cgmlst.org/ncs/schema/Kpneumoniae_complex/"
    assert root == str(tmp_path / "schemes"), "the installer redirects a cgMLST scheme itself"


def test_the_pinned_cgmlst_catalogue_says_where_each_scheme_installs_before_it_is_downloaded(
        qtbot, tmp_path):
    dialog = open_dialog(qtbot, tmp_path)
    page = dialog.pages["cgmlst"]
    page.select_row_id("key:" + KLEBSIELLA)
    row, detail = page.selected_row(), page.notice.toPlainText()

    assert row["status"] == "Not installed · download under the provider's terms"
    assert row["location"] == str(tmp_path / "cgmlst" / "Klebsiella_pneumoniae__cgmlst_org_2358")
    assert "2358 targets" in row["title"] and row["target_set"] == "core"
    assert row["location"] in detail
    assert "Ridom GmbH" in detail and "https://www.cgmlst.org/serverpolicy.html" in detail
    assert "never shares a scale with a seven-locus MLST distance" in detail
    # Counted from the catalogue rather than written down here, so pinning another
    # scheme is a catalogue change and not a test failure. What matters is that
    # every pinned scheme is offered, not that there are exactly N of them.
    catalogued = len(cgmlst_schemes.catalog_entries())
    assert f"{catalogued} catalogued and ready to download" in page.summary.text()
    assert str(tmp_path / "cgmlst") in page.summary.text()


def test_an_online_search_files_every_result_into_the_library_it_would_install_into(qtbot, tmp_path):
    classical = {"id": "PubMLST:kp:1", "name": "MLST", "organism": "Klebsiella pneumoniae",
                 "locus_count": 7, "type": "MLST", "last_updated": "2026-01-01",
                 "url": "https://rest.pubmlst.org/db/pubmlst_klebsiella_seqdef/schemes/1",
                 "access_notice": ""}
    gene_by_gene = {**classical, "id": "PubMLST:kp:2", "name": "cgMLST", "locus_count": 2358,
                    "type": "cgMLST",
                    "url": "https://rest.pubmlst.org/db/pubmlst_klebsiella_seqdef/schemes/2"}
    dialog = open_dialog(qtbot, tmp_path, Catalog(schemes=[classical, gene_by_gene]))
    dialog.query.setText("Klebsiella")
    dialog.search()
    qtbot.waitUntil(lambda: all(page.online for page in dialog.pages.values()), timeout=3000)

    online = {kind: [row for row in page.rows if row["state"] == "online"]
              for kind, page in dialog.pages.items()}
    assert [row["count"] for row in online["mlst"]] == [7]
    assert [row["count"] for row in online["cgmlst"]] == [2358]
    assert "1 listed under cgMLST schemes and 1 under classical MLST schemes" in dialog.status.text()
    assert online["cgmlst"][0]["title"].endswith("2358 targets")
    assert online["mlst"][0]["title"].endswith("7 loci")


def test_a_long_step_says_how_far_through_it_is_and_what_is_left():
    """A step naming only the locus it finished looks the same at 10 and at 2,000.

    That is why a working cgMLST download read as a frozen one: the theme draws
    the bar as a seven-pixel sliver with transparent text, so the status line was
    the only thing a person could read, and it carried no position at all.
    """
    describe = reference_dialog.describe_progress
    assert describe(1204, 2358, "Checking allele files", 240) == (
        "Checking allele files · 1,204 of 2,358 (51%) · about 3 min 50 s left")
    # Finished work states its position and stops estimating.
    assert describe(2358, 2358, "Checking allele files", 470) == (
        "Checking allele files · 2,358 of 2,358 (100%)")
    # Too early to estimate honestly: report position only.
    assert "left" not in describe(3, 2358, "Checking allele files", 4)
    # A step with no countable total still shows it is moving.
    assert describe(0, 0, "Downloading: 12 MB received", 30) == (
        "Downloading: 12 MB received · 30 seconds so far")
    assert describe(0, 0, "Starting", 1) == "Starting"
    # A total that is overrun is reported as complete, never as more than all.
    assert describe(2400, 2358, "Checking", 10).endswith("2,358 of 2,358 (100%)")


def test_closing_during_a_download_asks_first_and_only_stops_when_told_to(qtbot, tmp_path,
                                                                            monkeypatch):
    """Closing used to cancel silently, losing a download that was nearly done."""
    dialog = ReferenceManagerDialog(tmp_path)
    qtbot.addWidget(dialog)

    class Running:
        def __init__(self):
            self.cancelled = False

        def isRunning(self):
            return True

        def cancel(self):
            self.cancelled = True

    dialog.worker = Running()

    # Declining leaves the download alone and the window open.
    monkeypatch.setattr(dialog, "confirm_stop", lambda: False)
    dialog.reject()
    assert dialog.worker.cancelled is False, "a declined close must not stop the download"
    assert dialog._close_pending is False

    # Accepting stops it, and only then.
    monkeypatch.setattr(dialog, "confirm_stop", lambda: True)
    dialog.reject()
    assert dialog.worker.cancelled is True and dialog._close_pending is True
    dialog.worker = None


def test_the_stop_question_names_what_happens_to_what_was_already_downloaded(qtbot, tmp_path,
                                                                            monkeypatch):
    dialog = ReferenceManagerDialog(tmp_path)
    qtbot.addWidget(dialog)
    shown = {}

    def question(parent, title, text, buttons, default):
        shown["title"], shown["text"] = title, text
        return QMessageBox.StandardButton.Cancel

    monkeypatch.setattr(QMessageBox, "question", question)
    assert dialog.confirm_stop() is False, "the safe answer is the default"
    assert "kept" in shown["text"] and "continues from where it stopped" in shown["text"]


def test_each_stage_of_a_download_is_timed_from_its_own_start(qtbot, tmp_path):
    """A stage restarting its count must not inherit the previous stage's clock.

    Extraction reports 2,358 loci over twenty minutes; the check that follows
    counts the same 2,358 again. Dividing twenty minutes by that stage's first
    locus would promise hours left on work that takes a couple of minutes — a
    wrong number, which is worse than the silence it replaced.
    """
    dialog = ReferenceManagerDialog(tmp_path)
    qtbot.addWidget(dialog)
    dialog._started = time.monotonic() - 1200
    dialog._progress_at, dialog._progress_total = 2358, 2358
    # The message that announces the next stage restarts the count at zero.
    dialog._progress(0, 2358, "Every allele file is now read back")
    assert time.monotonic() - dialog._started < 5, "the new stage starts its own clock"
    dialog._progress(600, 2358, "Checking what was downloaded: Reading locus KP1_RS00100")
    assert "600 of 2,358" in dialog.status.text()
    assert "2 h" not in dialog.status.text()
    # A stage with a different total is a different stage too.
    dialog._started = time.monotonic() - 1200
    dialog._progress(1, 2360, "Checking what was downloaded: Fingerprinting arc.fasta")
    assert time.monotonic() - dialog._started < 5
