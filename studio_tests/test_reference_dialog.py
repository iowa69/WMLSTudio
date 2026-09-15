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
        qtbot, tmp_path, monkeypatch):
    dialog = open_dialog(qtbot, tmp_path)
    page = dialog.pages["cgmlst"]
    page.query.setText("Klebsiella pneumoniae")
    page.table.selectRow(0)

    # Only a core set is pinned for this organism: that is said in words, and the
    # choice is not offered as an empty menu.
    assert page.variant.count() == 1 and not page.variant.isEnabled()
    assert "Only a core target set is catalogued" in page.variant_note.text()

    core = cgmlst_schemes.entry_for(KLEBSIELLA)
    accessory = {**core, "key": "cgmlst.org:kpneumoniae-accessory-900",
                 "target_set": "accessory", "locus_count": 900,
                 "scheme_name": "Klebsiella pneumoniae accessory set",
                 "title": "Klebsiella pneumoniae · accessory set · 900 targets · cgMLST.org"}
    monkeypatch.setattr(cgmlst_schemes, "scheme_variants", lambda *args, **kwargs: {
        "organism": "Klebsiella pneumoniae", "group": core["scheme_group"], "core": [core],
        "accessory": [accessory], "has_core": True, "has_accessory": True,
        "message": "They are different quantities: a distance from one never shares a scale."})
    page.show_entry()

    assert page.variant.isEnabled() and page.variant.count() == 2
    assert page.variant.itemText(0).startswith("Core target set — ")
    assert page.variant.itemText(1).startswith("Accessory / whole-genome set — ")
    assert page.variant.currentData() == KLEBSIELLA
    assert "different quantities" in page.variant_note.text()


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
