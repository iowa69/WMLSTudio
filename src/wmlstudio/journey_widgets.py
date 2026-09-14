"""Responsive native investigation map; actions are routed by the main window."""

from PySide6.QtCore import Signal
from PySide6.QtWidgets import QGridLayout, QHBoxLayout, QVBoxLayout, QWidget

from wmlstudio.widgets import button, card, label


class InvestigationMap(QWidget):
    actionRequested = Signal(str)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setAccessibleName("Investigation workflow map")
        self._cards = {}
        self._columns = 0
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 14, 0, 0)
        self.scope = label("No isolates imported yet", "cardTitle", True)
        layout.addWidget(self.scope)
        self.context = label("One project holds your evidence. Selections define what you investigate and report.", "muted", True)
        layout.addWidget(self.context)
        self.grid = QGridLayout()
        self.grid.setSpacing(12)
        layout.addLayout(self.grid)
        controls = QHBoxLayout()
        controls.addWidget(button("Import / add next batch…", lambda: self.actionRequested.emit("import"), True))
        controls.addWidget(button("Review flagged isolates", lambda: self.actionRequested.emit("review")))
        controls.addWidget(button("Problem → solution guide", lambda: self.actionRequested.emit("guide")))
        controls.addStretch()
        # The shared wrapping layout keeps all actions reachable at compact sizes.
        from wmlstudio.ui_common import FlowLayout
        wrapped = FlowLayout()
        while controls.count():
            item = controls.takeAt(0)
            if item.widget():
                wrapped.addWidget(item.widget())
        layout.addLayout(wrapped)
        layout.addStretch()

    def set_summary(self, summary, investigation=None):
        counts = summary["counts"]
        self.scope.setText(f"{counts['total']} isolates in project · Choose isolates separately for each analysis or report")
        if investigation:
            self.context.setText(str(investigation))
        else:
            self.context.setText("Start with your research question. Each step links to the same saved isolates and evidence.")
        for step in summary["steps"]:
            if step.key not in self._cards:
                frame, content = card()
                frame.setMinimumWidth(260)
                frame.setAccessibleName(step.title)
                title = label(step.title, "cardTitle", True)
                question = label(step.question, "muted", True)
                description = label(step.summary, "small", True)
                action = button(step.action_label, lambda checked=False, key=step.action: self.actionRequested.emit(key))
                content.addWidget(title)
                content.addWidget(question)
                content.addWidget(description)
                content.addStretch()
                content.addWidget(action)
                self._cards[step.key] = (frame, description)
                self._columns = 0
            frame, description = self._cards[step.key]
            description.setText(step.summary)
            frame.setToolTip("Needs evidence review" if step.state == "review" else step.question)
        self._arrange()

    def _arrange(self):
        columns = 2 if self.width() >= 740 else 1
        if columns == self._columns:
            return
        self._columns = columns
        for index, (frame, _) in enumerate(self._cards.values()):
            self.grid.addWidget(frame, index // columns, index % columns)
        self.grid.setColumnStretch(0, 1)
        self.grid.setColumnStretch(1, 1 if columns == 2 else 0)

    def resizeEvent(self, event):
        super().resizeEvent(event)
        self._arrange()
