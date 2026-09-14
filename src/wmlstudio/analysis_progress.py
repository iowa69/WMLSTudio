"""Focused modal progress without blocking Qt's worker or cancellation signals."""

from PySide6.QtCore import QElapsedTimer, Qt, QTimer
from PySide6.QtWidgets import QDialog, QHBoxLayout, QPlainTextEdit, QProgressBar, QVBoxLayout

from wmlstudio.widgets import button, label


class AnalysisProgressDialog(QDialog):
    def __init__(self, window):
        super().__init__(window)
        self.window = window
        self.finished_safely = False
        self._last_message = None
        self.setWindowTitle("Analysis in progress")
        self.setWindowModality(Qt.WindowModality.WindowModal)
        self.resize(640, 270)
        layout = QVBoxLayout(self)
        layout.addWidget(label("Your analysis is running", "title", True))
        self.message = label("Preparing the reviewed workflow…", "muted", True)
        layout.addWidget(self.message)
        self.progress = QProgressBar()
        self.progress.setRange(0, 100)
        layout.addWidget(self.progress)
        self.elapsed = label("", "small")
        layout.addWidget(self.elapsed)
        layout.addWidget(label("Workspace changes are paused until this stage finishes. Completed evidence is retained if you cancel; original inputs stay unchanged.", "small", True))
        self.logs = QPlainTextEdit()
        self.logs.setReadOnly(True)
        self.logs.setMaximumBlockCount(1000)
        self.logs.hide()
        layout.addWidget(self.logs, 1)
        controls = QHBoxLayout()
        self.details = button("Show progress log", self.toggle_details)
        controls.addWidget(self.details)
        controls.addStretch()
        self.cancel = button("Cancel safely", self.request_cancel)
        controls.addWidget(self.cancel)
        layout.addLayout(controls)
        self.clock = QElapsedTimer()
        self.clock.start()
        self.timer = QTimer(self)
        self.timer.timeout.connect(self.update_elapsed)
        self.timer.start(1000)
        self.update_elapsed()

    def update_elapsed(self):
        seconds = self.clock.elapsed() // 1000
        self.elapsed.setText(f"Elapsed {seconds // 60:02d}:{seconds % 60:02d} · progress is stage-specific")

    def update_progress(self, percent, message):
        self.progress.setValue(max(0, min(100, int(percent))))
        self.message.setText(str(message))
        if self._last_message != str(message):
            self.logs.appendPlainText(str(message))
            self._last_message = str(message)

    def toggle_details(self):
        show = not self.logs.isVisible()
        self.logs.setVisible(show)
        self.details.setText("Hide progress log" if show else "Show progress log")
        self.resize(self.width(), 520 if show else 270)

    def request_cancel(self):
        self.cancel.setEnabled(False)
        self.message.setText("Cancelling safely; waiting for active processes to finish…")
        self.window.cancel_analysis()

    def reject(self):
        if self.finished_safely:
            super().reject()
        else:
            self.request_cancel()

    def finish(self):
        self.finished_safely = True
        self.timer.stop()
        self.accept()
