"""Searchable offline problem-to-solution handbook, rendered by native Qt."""

from __future__ import annotations

import sys
from pathlib import Path

from PySide6.QtCore import Qt, Signal
from PySide6.QtGui import QDesktopServices
from PySide6.QtWidgets import (
    QDialog, QDialogButtonBox, QHBoxLayout, QLineEdit, QListWidget,
    QListWidgetItem, QSplitter, QTextBrowser, QVBoxLayout,
)

GUIDE_ACTIONS = frozenset({"import", "samples", "review", "analyse", "characterize",
                           "features", "compare", "reports"})


def guide_path() -> Path:
    root = Path(getattr(sys, "_MEIPASS", Path(__file__).resolve().parents[2]))
    path = root / "docs" / "WORKFLOW_GUIDE.md"
    if not path.is_file():
        raise FileNotFoundError("The offline workflow guide is missing from this application bundle.")
    return path


def guide_chapters(text: str) -> list[tuple[str, str]]:
    chapters = []
    title, lines = "Start here", []
    for line in text.splitlines():
        if line.startswith("## "):
            if lines:
                chapters.append((title, "\n".join(lines).strip()))
            title, lines = line[3:].strip(), [line]
        else:
            lines.append(line)
    if lines:
        chapters.append((title, "\n".join(lines).strip()))
    return chapters


class WorkflowGuide(QDialog):
    actionRequested = Signal(str)

    def __init__(self, parent=None, topic=None):
        super().__init__(parent)
        self.setWindowTitle("WMLSTudio · Problem → solution guide")
        self.resize(1040, 740)
        self.setMinimumSize(720, 480)
        self.chapters = guide_chapters(guide_path().read_text(encoding="utf-8"))
        layout = QVBoxLayout(self)
        self.search = QLineEdit()
        self.search.setPlaceholderText("Search a question: plasmids, next week, missing loci, resistance…")
        self.search.setClearButtonEnabled(True)
        self.search.textChanged.connect(self.filter_chapters)
        layout.addWidget(self.search)
        splitter = QSplitter(Qt.Orientation.Horizontal)
        self.contents = QListWidget()
        self.contents.setAccessibleName("Workflow guide topics")
        self.contents.setWordWrap(True)
        self.contents.currentItemChanged.connect(self.show_chapter)
        splitter.addWidget(self.contents)
        self.browser = QTextBrowser()
        self.browser.setOpenLinks(False)
        self.browser.setOpenExternalLinks(False)
        self.browser.anchorClicked.connect(self.follow_link)
        self.browser.setAccessibleName("Workflow explanation and evidence boundaries")
        splitter.addWidget(self.browser)
        splitter.setSizes([270, 720])
        layout.addWidget(splitter, 1)
        bottom = QHBoxLayout()
        buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Close)
        buttons.rejected.connect(self.reject)
        bottom.addStretch()
        bottom.addWidget(buttons)
        layout.addLayout(bottom)
        self.filter_chapters(topic or "")

    def filter_chapters(self, query):
        self.contents.clear()
        query = query.casefold().strip()
        for index, (title, body) in enumerate(self.chapters):
            if query and query not in (title + "\n" + body).casefold():
                continue
            item = QListWidgetItem(title)
            item.setData(Qt.ItemDataRole.UserRole, index)
            self.contents.addItem(item)
        if self.contents.count():
            self.contents.setCurrentRow(0)
        else:
            self.browser.setPlainText("No matching topics. Try a shorter term or clear the search.")

    def show_chapter(self, current, previous=None):
        if current is not None:
            self.browser.setMarkdown(self.chapters[current.data(Qt.ItemDataRole.UserRole)][1])

    def follow_link(self, url):
        if url.scheme() == "wmlstudio":
            action = url.path().strip("/")
            if action in GUIDE_ACTIONS:
                self.accept()
                self.actionRequested.emit(action)
        elif url.scheme() == "https":
            QDesktopServices.openUrl(url)
        elif not url.scheme() and url.fragment():
            self.browser.scrollToAnchor(url.fragment())
