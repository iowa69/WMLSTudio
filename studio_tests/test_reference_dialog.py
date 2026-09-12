import threading

from PySide6.QtCore import QThread
from PySide6.QtTest import QSignalSpy

from wmlstudio.reference_dialog import ReferenceManagerDialog
from wmlstudio.sequence import check_cancelled


class Catalog:
    def __init__(self):
        self.calls = []

    def search_schemes(self, query, *, cancelled, progress, **filters):
        self.calls.append(('search', QThread.currentThread()))
        check_cancelled(cancelled)
        return {'schemes': [{'name': 'Example MLST', 'organism': 'Examplegenus species',
                             'locus_count': 7, 'type': 'MLST', 'last_updated': '2026-01-01',
                             'url': 'https://rest.pubmlst.org/db/example/schemes/1',
                             'access_notice': 'Anonymous access excludes recent alleles.'}],
                'organisms_searched': 2,
                'errors': [{'organism': 'Unavailable genus', 'error': 'offline'}]}

    def download_scheme(self, entry, root, *, cancelled, progress):
        self.calls.append(('download', QThread.currentThread()))
        check_cancelled(cancelled)
        return {'path': str(root / 'verified'), 'created': True, 'locus_count': 7, 'notes': []}


def test_catalog_dialog_does_not_network_until_search_and_workers_stay_off_gui(qtbot, tmp_path, qapp):
    client = Catalog()
    dialog = ReferenceManagerDialog(tmp_path, catalog=client)
    qtbot.addWidget(dialog)
    assert client.calls == []
    dialog.query.setText('Examplegenus')
    dialog.search()
    qtbot.waitUntil(lambda: dialog.table.rowCount() == 1 and not dialog.worker.isRunning(), timeout=3000)
    assert 'Anonymous access' in dialog.notice.toPlainText()
    assert 'Unavailable genus: offline' in dialog.notice.toPlainText()
    assert dialog.download_button.isEnabled()
    installed = QSignalSpy(dialog.schemeInstalled)
    dialog.download()
    qtbot.waitUntil(lambda: installed.count() == 1 and not dialog.worker.isRunning(), timeout=3000)
    assert installed.at(0)[0] == str(tmp_path / 'schemes' / 'verified')
    assert [name for name, _ in client.calls] == ['search', 'download']
    assert all(thread is not qapp.thread() for _, thread in client.calls)
    dialog.close()


def test_dialog_close_cancels_active_network_work_without_publishing(qtbot, tmp_path):
    entered = threading.Event()
    class BlockingCatalog(Catalog):
        def search_schemes(self, query, *, cancelled, progress, **filters):
            entered.set()
            while not cancelled():
                threading.Event().wait(.005)
            check_cancelled(cancelled)
    dialog = ReferenceManagerDialog(tmp_path, catalog=BlockingCatalog())
    qtbot.addWidget(dialog)
    dialog.show()
    dialog.query.setText('Examplegenus')
    dialog.search()
    qtbot.waitUntil(entered.is_set, timeout=3000)
    installed = QSignalSpy(dialog.schemeInstalled)
    dialog.reject()
    qtbot.waitUntil(lambda: not dialog.worker.isRunning() and not dialog.isVisible(), timeout=3000)
    assert installed.count() == 0
    assert 'cancelled' in dialog.status.text().lower()
