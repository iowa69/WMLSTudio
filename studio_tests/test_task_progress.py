"""What a long task tells you while it runs, and what it is allowed to call itself.

Two reports, the same complaint in two places. Pressing Install on the Update
page raised a dialog headed "Your analysis is running", so the update button was
read as running analyses instead of fetching anything. And a cgMLST download that
reached 2,358 of 2,358 sat there with a full bar for several minutes while it
re-read every allele file, reported as hung.

Neither was broken underneath. Both were the words.
"""

import pytest
from PySide6.QtWidgets import QWidget

from wmlstudio.analysis_progress import ANALYSIS, AnalysisProgressDialog, installing
from wmlstudio.progress_words import describe_progress, duration_words, remaining_words


class Host(QWidget):
    """Just enough window for the dialog: it only ever calls cancel_analysis."""

    def __init__(self):
        super().__init__()
        self.cancelled = 0

    def cancel_analysis(self):
        self.cancelled += 1


@pytest.fixture
def host(qtbot):
    widget = Host()
    qtbot.addWidget(widget)
    return widget


def test_seconds_are_said_the_way_a_person_says_them():
    assert duration_words(0) == "0 seconds"
    assert duration_words(45) == "45 seconds"
    assert duration_words(60) == "1 min"
    assert duration_words(150) == "2 min 30 s"
    assert duration_words(7400) == "2 h 3 min"
    assert duration_words(-5) == "0 seconds"


def test_an_estimate_is_offered_only_once_it_would_mean_something():
    # Two percent in, "about 4 h left" is noise dressed as information.
    assert remaining_words(2, 30) == ""
    assert remaining_words(50, 2) == ""
    assert remaining_words(100, 600) == ""
    assert remaining_words(50, 60) == "about 1 min left"
    assert remaining_words(25, 60) == "about 3 min left"
    assert remaining_words("nonsense", 60) == ""


def test_a_counted_step_says_where_it_has_got_to():
    """The locus name alone reads identically at locus 10 and at locus 2,300."""
    said = describe_progress(1204, 2358, "Checking allele files", 240)
    assert "1,204 of 2,358 (51%)" in said
    assert "about" in said and "left" in said
    early = describe_progress(2, 2358, "Checking allele files", 1)
    assert "2 of 2,358" in early and "left" not in early
    assert describe_progress(0, 0, "Starting", 1) == "Starting"
    assert "so far" in describe_progress(0, 0, "Starting", 30)


def test_a_download_is_never_announced_as_an_analysis(host, qtbot):
    """The reported bug, in one assertion: these words belong to different work."""
    caption = installing("the species reference panel")
    dialog = AnalysisProgressDialog(host, caption)
    qtbot.addWidget(dialog)
    assert dialog.windowTitle() == "Installing reference data"
    assert dialog.heading.text() == "Installing the species reference panel"
    for text in (dialog.windowTitle(), dialog.heading.text(), dialog.note.text()):
        assert "analysis" not in text.casefold() or "not an analysis" in text.casefold()
    assert "no sequence of yours leaves this computer" in dialog.note.text()
    # Stopping a download is not "cancelling safely" in the sense an analysis
    # means it, and it promises nothing about resuming: one of these downloads
    # continues from where it stopped and the other starts again.
    assert dialog.cancel.text() == "Stop the download"
    assert "resum" not in dialog.caption.cancelling.casefold()
    assert "continue" not in dialog.caption.cancelling.casefold()


def test_the_analysis_wording_is_still_what_an_analysis_gets(host, qtbot):
    dialog = AnalysisProgressDialog(host)
    qtbot.addWidget(dialog)
    assert dialog.windowTitle() == "Analysis in progress"
    assert dialog.heading.text() == "Your analysis is running"
    assert dialog.caption is ANALYSIS
    assert dialog.cancel.text() == "Cancel safely"


def test_a_dialog_shown_again_is_relabelled_rather_than_keeping_the_last_task(host, qtbot):
    """It is reused between tasks, so a stale caption would be a lie about this one."""
    dialog = AnalysisProgressDialog(host)
    qtbot.addWidget(dialog)
    dialog.update_progress(40, "Typing isolate 4")
    dialog.set_caption(installing("the NCBI reference set"))
    assert dialog.windowTitle() == "Installing reference data"
    assert dialog.heading.text() == "Installing the NCBI reference set"
    # And it starts from nothing: the previous task's bar and message would read
    # as progress this one has not made.
    assert dialog.progress.value() == 0
    assert dialog.message.text() == dialog.caption.starting
    assert "40" not in dialog.elapsed.text()
    dialog.set_caption(None)
    assert dialog.windowTitle() == "Analysis in progress"


def test_the_elapsed_line_answers_is_this_getting_anywhere(host, qtbot):
    dialog = AnalysisProgressDialog(host)
    qtbot.addWidget(dialog)
    # Nothing has reported a fraction yet, and the old line ("progress is
    # stage-specific") could not tell a half-finished step from a stuck one.
    assert "has not reported how far along it is yet" in dialog.describe_elapsed()
    dialog.update_progress(50, "Halfway")
    dialog._percent, dialog.clock = 50, dialog.clock
    assert "50%" in dialog.describe_elapsed()
    assert dialog.progress.value() == 50


def test_cancelling_says_what_stopping_means_for_this_kind_of_work(host, qtbot):
    dialog = AnalysisProgressDialog(host, installing("the species reference panel"))
    qtbot.addWidget(dialog)
    dialog.request_cancel()
    assert host.cancelled == 1
    assert "already installed is changed or removed" in dialog.message.text()
    assert not dialog.cancel.isEnabled()
