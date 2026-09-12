"""Cancellable desktop tasks with result delivery confined to the Qt GUI thread."""

import threading

from PySide6.QtCore import QThread, Signal

from wmlstudio.sequence import AnalysisCancelled


class FunctionWorker(QThread):
    completed = Signal(object)
    failed = Signal(str)
    progress = Signal(int, str)

    def __init__(self, operation, parent=None):
        super().__init__(parent)
        self.operation = operation
        self.cancel_event = threading.Event()

    def cancel(self):
        self.cancel_event.set()

    def report(self, done, total, message):
        self.progress.emit(int(done / max(1, total) * 100), str(message))

    def run(self):
        try:
            result = self.operation(self.cancel_event.is_set, self.report)
            self.completed.emit(result)
        except AnalysisCancelled:
            self.progress.emit(100, "Operation cancelled; original inputs are unchanged.")
        except Exception as exc:
            self.failed.emit(str(exc))
