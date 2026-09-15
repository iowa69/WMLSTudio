"""Focused modal progress without blocking Qt's worker or cancellation signals.

The dialog says what is actually running, because it was reported twice that it
did not. A person pressed Install on the Update page, was told in bold that
"your analysis is running", and concluded that the update button runs an
analysis instead of fetching anything. Nothing was wrong underneath — the
download was working — but every word on the dialog belonged to a different kind
of work, so the caption comes from whoever starts the task rather than from here.

It also answers "is this getting anywhere". A bar with no number beside it and a
line reading "progress is stage-specific" cannot distinguish a step that is
half-done from one that is stuck, which is the same complaint in its other form.
"""

from dataclasses import dataclass

from PySide6.QtCore import QElapsedTimer, Qt, QTimer
from PySide6.QtWidgets import QDialog, QHBoxLayout, QPlainTextEdit, QProgressBar, QVBoxLayout

from wmlstudio.progress_words import remaining_words
from wmlstudio.widgets import button, label


@dataclass(frozen=True)
class TaskCaption:
    """What this particular task is, in the words the person who started it used."""

    title: str
    heading: str
    note: str
    starting: str
    cancel: str = "Cancel safely"
    cancelling: str = "Cancelling safely; waiting for active processes to finish…"


#: Computing something from the user's own sequences. The original, and still the
#: default for every caller that does not say otherwise.
ANALYSIS = TaskCaption(
    title="Analysis in progress",
    heading="Your analysis is running",
    note="Workspace changes are paused until this stage finishes. Completed evidence is "
         "retained if you cancel; original inputs stay unchanged.",
    starting="Preparing the reviewed workflow…")


def installing(what: str) -> TaskCaption:
    """The caption for fetching reference data. This computes nothing from anybody's samples.

    Deliberately makes no promise about resuming: a scheme download continues from
    where it stopped, the species panel starts again, and a dialog shared by both
    may only say what is true of both.
    """
    return TaskCaption(
        title="Installing reference data",
        heading=f"Installing {what}",
        note="This is a download, not an analysis: nothing is being computed from your samples, "
             "and no sequence of yours leaves this computer. What is already installed stays "
             "exactly as it is until the new copy has been fetched and checked.",
        starting="Contacting the provider…",
        cancel="Stop the download",
        cancelling="Stopping. Nothing that is already installed is changed or removed.")


class AnalysisProgressDialog(QDialog):
    def __init__(self, window, caption=None):
        super().__init__(window)
        self.window = window
        self.finished_safely = False
        self._last_message = None
        self._percent = 0
        self.caption = caption or ANALYSIS
        self.setWindowModality(Qt.WindowModality.WindowModal)
        self.resize(640, 270)
        layout = QVBoxLayout(self)
        self.heading = label(self.caption.heading, "title", True)
        layout.addWidget(self.heading)
        self.message = label(self.caption.starting, "muted", True)
        layout.addWidget(self.message)
        self.progress = QProgressBar()
        self.progress.setRange(0, 100)
        layout.addWidget(self.progress)
        self.elapsed = label("", "small")
        layout.addWidget(self.elapsed)
        self.note = label(self.caption.note, "small", True)
        layout.addWidget(self.note)
        self.logs = QPlainTextEdit()
        self.logs.setReadOnly(True)
        self.logs.setMaximumBlockCount(1000)
        self.logs.hide()
        layout.addWidget(self.logs, 1)
        controls = QHBoxLayout()
        self.details = button("Show progress log", self.toggle_details)
        controls.addWidget(self.details)
        controls.addStretch()
        self.cancel = button(self.caption.cancel, self.request_cancel)
        controls.addWidget(self.cancel)
        layout.addLayout(controls)
        self.clock = QElapsedTimer()
        self.clock.start()
        self.timer = QTimer(self)
        self.timer.timeout.connect(self.update_elapsed)
        self.timer.start(1000)
        self.set_caption(self.caption)

    def set_caption(self, caption):
        """Re-label a dialog that is being shown again for a different kind of work."""
        self.caption = caption or ANALYSIS
        self.setWindowTitle(self.caption.title)
        self.heading.setText(self.caption.heading)
        self.note.setText(self.caption.note)
        self.cancel.setText(self.caption.cancel)
        self.cancel.setEnabled(True)
        self._last_message = None
        self._percent = 0
        self.progress.setValue(0)
        self.message.setText(self.caption.starting)
        self.clock.restart()
        self.update_elapsed()
        return self.caption

    def describe_elapsed(self) -> str:
        """How long this has been running, and — once it can be told — how long is left."""
        seconds = self.clock.elapsed() // 1000
        clock = f"Elapsed {seconds // 60:02d}:{seconds % 60:02d}"
        left = remaining_words(self._percent, seconds)
        if self._percent <= 0:
            # No step has reported a fraction yet, so there is nothing to divide:
            # saying so is the honest version of "progress is stage-specific".
            return f"{clock} · this step has not reported how far along it is yet"
        return f"{clock} · {self._percent}%" + (f" · {left}" if left else "")

    def update_elapsed(self):
        self.elapsed.setText(self.describe_elapsed())

    def update_progress(self, percent, message):
        self._percent = max(0, min(100, int(percent)))
        self.progress.setValue(self._percent)
        self.message.setText(str(message))
        self.update_elapsed()
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
        self.message.setText(self.caption.cancelling)
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
