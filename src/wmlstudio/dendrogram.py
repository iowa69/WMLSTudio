"""A dendrogram on screen: what :mod:`wmlstudio.hierarchical` computed, and nothing else.

This module draws. It decides no distance, no height and no cluster: every
number it paints comes from :func:`wmlstudio.hierarchical.build_tree` and
:func:`wmlstudio.hierarchical.cut_tree`, and every caveat those two record
travels to the notes beside the picture rather than being dropped because it did
not fit. The drawing is plain Qt, so a dendrogram costs no new dependency.

Three things about this picture are easy to misread, and the widget is built
around saying them out loud:

* **The height axis has a quantity.** A seven-locus MLST height, a cgMLST target
  height and a SNP height are different numbers that never share an axis, so the
  axis is captioned with the ``scale_caption`` the tree was built with. A tree
  built without one says the quantity was not recorded rather than leaving the
  axis bare.
* **A cut is a line drawn on a picture.** It groups isolates the tree has already
  placed; it measures nothing new, and it is never evidence on its own.
* **A cluster can be held together by a chain.** Under single linkage two
  isolates far apart can sit in one group because a run of near neighbours joins
  them. Where the core reports that, the join is marked in the drawing and the
  note is printed under it.
"""

from __future__ import annotations

import html
import math

from PySide6.QtCore import QEvent, QPointF, QRectF, Qt, Signal
from PySide6.QtGui import QColor, QFontMetricsF, QPainter, QPen, QPolygonF
from PySide6.QtWidgets import (
    QComboBox,
    QDoubleSpinBox,
    QHBoxLayout,
    QSizePolicy,
    QSplitter,
    QTextBrowser,
    QToolTip,
    QVBoxLayout,
    QWidget,
)

from wmlstudio.hierarchical import (
    CHAINING_RULE,
    LINKAGES,
    MISSING_DATA_RULE,
    TOLERANCE,
    cut_tree,
    dendrogram_layout,
)
from wmlstudio.theme import BORDER, INK, MUTED, PALETTE
from wmlstudio.widgets import button, label

# The three linkages the core offers, in the order a reader meets them: the one
# that chains, the one that cannot, and the mean of the two.
LINKAGE_ORDER = ("single", "complete", "average")

CUT_RULE = (
    "A cut is a line drawn on a picture. It groups the isolates this tree already placed and it "
    "measures nothing new, so it is never evidence of an outbreak on its own."
)
HEIGHT_UNRECORDED = "Height — the quantity this tree was built from was not recorded"
# How many blocked pairs and chained joins are printed in full before the notes
# start counting instead. Everything is still counted; nothing is dropped.
NOTE_LIMIT = 8


def _tick_step(span: float, target: int = 5) -> float:
    """A round step near ``span / target``, so the axis reads 0, 5, 10 and not 0, 4.7."""
    if span <= 0:
        return 1.0
    raw = span / max(1, target)
    power = 10.0 ** math.floor(math.log10(raw))
    for multiple in (1, 2, 2.5, 5, 10):
        if raw <= multiple * power:
            return multiple * power
    return 10 * power


def _format(value: float) -> str:
    return f"{value:g}"


class DendrogramView(QWidget):
    """A rectangular cladogram with a height axis and one draggable cut line.

    Leaves sit at the left edge and heights grow to the right, so an isolate's
    name reads horizontally beside its own tip. The clusters below the cut are
    coloured; everything above it is drawn in one neutral colour, because a
    branch above the cut belongs to no group the cut made.
    """

    cutChanged = Signal(float)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setMinimumHeight(220)
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        self.setMouseTracking(True)
        self.setAutoFillBackground(False)
        self._tree = {}
        self._layout = {"leaves": [], "nodes": [], "max_height": 0.0, "height_caption": ""}
        self._names = {}
        self._members = {}
        self._cut = 0.0
        self._groups = {}
        self._message = "No comparison has been built yet, so there is no tree to draw."
        self._plan = {}
        self._dragging = False

    # --- what is drawn ---------------------------------------------------
    def set_tree(self, tree, names=None, cut=None):
        """Draw one tree. ``names`` maps a leaf label to the name to print for it."""
        self._tree = dict(tree or {})
        self._layout = dendrogram_layout(self._tree) if self._tree else {
            "leaves": [], "nodes": [], "max_height": 0.0, "height_caption": ""}
        self._names = dict(names or {})
        self._members = {merge["id"]: list(merge["members"])
                         for merge in self._tree.get("merges", [])}
        self._message = "" if self._layout["leaves"] else (
            "This cohort produced no tree to draw. Nothing was clustered, and that is not a "
            "finding that the isolates are alike.")
        self.set_cut(self._cut if cut is None else cut, notify=False)

    def clear(self, message=""):
        self._tree, self._members, self._groups = {}, {}, {}
        self._layout = {"leaves": [], "nodes": [], "max_height": 0.0, "height_caption": ""}
        self._message = message or "No comparison has been built yet, so there is no tree to draw."
        self._plan = {}
        self.update()

    def set_cut(self, height, notify=True):
        try:
            height = max(0.0, float(height))
        except (TypeError, ValueError):
            return
        changed = abs(height - self._cut) > TOLERANCE
        self._cut = height
        self._groups = cut_tree(self._tree, height) if self._tree else {}
        self._plan = {}
        self.update()
        if changed and notify:
            self.cutChanged.emit(height)

    def cut(self) -> float:
        return self._cut

    def groups(self) -> dict:
        """The cut as the core reported it, for whoever has to describe it in words."""
        return self._groups

    def axis_caption(self) -> str:
        """Never empty: an unrecorded quantity says so rather than leaving the axis bare."""
        return self._layout.get("height_caption") or HEIGHT_UNRECORDED

    # --- geometry, separate from painting so it can be read in a test -----
    def plan(self, rect=None) -> dict:
        """Where every leaf, branch, tick and the cut line go inside ``rect``."""
        rect = QRectF(rect) if rect is not None else QRectF(self.rect())
        leaves = self._layout.get("leaves", [])
        plan = {"leaves": [], "branches": [], "ticks": [], "cut": None, "labels_drawn": False,
                "caption": self.axis_caption(), "message": self._message, "rect": rect}
        if not leaves or rect.width() < 60 or rect.height() < 90:
            self._plan = plan
            return plan
        metrics = QFontMetricsF(self.font())
        rows = rect.height() - 64
        row_height = rows / len(leaves)
        labels_drawn = row_height >= metrics.height() + 1
        names = [self._names.get(leaf["label"], leaf["label"]) for leaf in leaves]
        gutter = 10.0
        if labels_drawn:
            widest = max((metrics.horizontalAdvance(name) for name in names), default=0.0)
            gutter = min(rect.width() * 0.4, widest + 14.0)
        left = rect.left() + gutter
        right = rect.right() - 14
        top = rect.top() + 26
        bottom = top + rows
        span = max(self._layout.get("max_height", 0.0), self._cut, TOLERANCE)
        span *= 1.08
        width = max(right - left, 1.0)

        def position(height):
            return left + (float(height) / span) * width

        def line(index):
            return top + (index + 0.5) * row_height

        group_of, colour_of = self._group_colours()
        for index, leaf in enumerate(leaves):
            group = group_of.get(leaf["label"], {})
            plan["leaves"].append({
                "label": leaf["label"], "name": names[index],
                "x": position(0.0), "y": line(index),
                "colour": colour_of.get(group.get("id"), MUTED),
                "status": group.get("status", "singleton"),
                "group": group.get("name", ""),
            })
        positions = {"leaf:" + leaf["label"]: (position(0.0), line(index))
                     for index, leaf in enumerate(leaves)}
        for node in self._layout.get("nodes", []):
            point = (position(node["y"]), line(node["x"]))
            positions[node["id"]] = point
            members = self._members.get(node["id"], [])
            ids = {group_of.get(member, {}).get("id") for member in members}
            inside = (len(ids) == 1 and members
                      and group_of.get(members[0], {}).get("status") == "cluster"
                      and node["y"] <= self._cut + TOLERANCE)
            colour = colour_of.get(group_of[members[0]]["id"], MUTED) if inside else BORDER
            left_point = positions[node["left"]]
            right_point = positions[node["right"]]
            plan["branches"].append({
                "id": node["id"], "x": point[0], "y": point[1],
                "colour": colour, "inside": bool(inside),
                "provisional": bool(node["provisional"]), "chained": bool(node["chained"]),
                "height": node["y"],
                "segments": [(left_point[0], left_point[1], point[0], left_point[1]),
                             (right_point[0], right_point[1], point[0], right_point[1]),
                             (point[0], left_point[1], point[0], right_point[1])],
            })
        step = _tick_step(span)
        value = 0.0
        while value <= span + TOLERANCE:
            plan["ticks"].append({"x": position(value), "text": _format(value)})
            value += step
        plan["cut"] = {"x": position(self._cut), "height": self._cut,
                       "top": top - 18, "bottom": bottom + 6}
        plan.update(labels_drawn=labels_drawn, left=left, right=right, top=top, bottom=bottom,
                    row_height=row_height, span=span)
        self._plan = plan
        return plan

    def _group_colours(self):
        """One colour per cut cluster; anything the cut did not cluster stays muted."""
        group_of, colour_of = {}, {}
        for group in self._groups.get("groups", []):
            for member in group["members"]:
                group_of[member] = group
            if group["status"] == "cluster" and group.get("number"):
                colour_of[group["id"]] = PALETTE[(group["number"] - 1) % len(PALETTE)]
        return group_of, colour_of

    # --- painting ---------------------------------------------------------
    def paintEvent(self, event):
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing, True)
        rect = QRectF(self.rect())
        plan = self.plan(rect)
        if not plan["leaves"]:
            painter.setPen(QPen(QColor(MUTED)))
            painter.drawText(rect.adjusted(16, 16, -16, -16),
                             int(Qt.AlignmentFlag.AlignTop | Qt.AlignmentFlag.AlignLeft)
                             | int(Qt.TextFlag.TextWordWrap), plan["message"])
            painter.end()
            return
        for branch in plan["branches"]:
            pen = QPen(QColor(branch["colour"]), 2.0 if branch["inside"] else 1.2)
            if branch["provisional"]:
                # A provisional join was measured over fewer pairs than it spans,
                # so it is drawn as a broken line rather than a solid claim.
                pen.setStyle(Qt.PenStyle.DashLine)
            painter.setPen(pen)
            for x1, y1, x2, y2 in branch["segments"]:
                painter.drawLine(QPointF(x1, y1), QPointF(x2, y2))
            if branch["chained"]:
                self._draw_chain_mark(painter, branch)
        metrics = QFontMetricsF(self.font())
        for leaf in plan["leaves"]:
            colour = QColor(leaf["colour"])
            painter.setPen(QPen(colour, 1.2))
            painter.setBrush(colour if leaf["status"] == "cluster" else Qt.BrushStyle.NoBrush)
            painter.drawEllipse(QPointF(leaf["x"], leaf["y"]), 3.0, 3.0)
            if plan["labels_drawn"]:
                painter.setBrush(Qt.BrushStyle.NoBrush)
                painter.setPen(QPen(QColor(INK if leaf["status"] == "cluster" else MUTED)))
                box = QRectF(plan["rect"].left(), leaf["y"] - metrics.height() / 2,
                             leaf["x"] - plan["rect"].left() - 6, metrics.height())
                painter.drawText(box, int(Qt.AlignmentFlag.AlignRight
                                          | Qt.AlignmentFlag.AlignVCenter),
                                 metrics.elidedText(leaf["name"], Qt.TextElideMode.ElideRight,
                                                    int(max(box.width(), 1))))
        if not plan["labels_drawn"]:
            # A name that was dropped for want of room is said to be dropped. A
            # reader must never take an unlabelled tip for an unnamed isolate.
            painter.setPen(QPen(QColor(MUTED)))
            note = (f"{len(plan['leaves'])} isolate names do not fit at this size — "
                    "hover a leaf or a join to name it.")
            painter.drawText(QRectF(plan["rect"].left() + 4, plan["rect"].top() + 2,
                                    plan["rect"].width() * 0.6, 18),
                             int(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter),
                             metrics.elidedText(note, Qt.TextElideMode.ElideRight,
                                                int(plan["rect"].width() * 0.6)))
        self._paint_axis(painter, plan)
        self._paint_cut(painter, plan)
        painter.end()

    def _draw_chain_mark(self, painter, branch):
        """Mark a join whose furthest compared pair is wider than the join itself."""
        painter.setBrush(QColor("#FFB976"))
        painter.setPen(QPen(QColor("#FFB976"), 1.0))
        x, y = branch["x"], branch["y"]
        painter.drawPolygon(QPolygonF([QPointF(x, y - 5), QPointF(x + 5, y + 4),
                                       QPointF(x - 5, y + 4)]))
        painter.setBrush(Qt.BrushStyle.NoBrush)

    def _paint_axis(self, painter, plan):
        metrics = QFontMetricsF(self.font())
        baseline = plan["bottom"] + 8
        painter.setPen(QPen(QColor(BORDER), 1.0))
        painter.drawLine(QPointF(plan["left"], baseline), QPointF(plan["right"], baseline))
        painter.setPen(QPen(QColor(MUTED)))
        for tick in plan["ticks"]:
            painter.drawLine(QPointF(tick["x"], baseline), QPointF(tick["x"], baseline + 4))
            painter.drawText(QRectF(tick["x"] - 30, baseline + 5, 60, metrics.height()),
                             int(Qt.AlignmentFlag.AlignHCenter), tick["text"])
        # The axis never goes out without its quantity: a height of 5 on seven loci
        # and a height of 5 on 2 358 targets are different numbers.
        painter.drawText(QRectF(plan["left"], baseline + 5 + metrics.height(),
                                max(plan["right"] - plan["left"], 1), metrics.height() + 2),
                         int(Qt.AlignmentFlag.AlignHCenter),
                         metrics.elidedText(plan["caption"], Qt.TextElideMode.ElideRight,
                                            int(max(plan["right"] - plan["left"], 1))))

    def _paint_cut(self, painter, plan):
        cut = plan["cut"]
        pen = QPen(QColor("#FF86A4"), 1.6)
        pen.setStyle(Qt.PenStyle.DashLine)
        painter.setPen(pen)
        painter.drawLine(QPointF(cut["x"], cut["top"]), QPointF(cut["x"], cut["bottom"]))
        painter.setPen(QPen(QColor("#FF86A4")))
        painter.drawText(QRectF(cut["x"] + 4, plan["rect"].top() + 2, 220, 18),
                         int(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter),
                         f"cut at {_format(cut['height'])}")

    # --- dragging the cut -------------------------------------------------
    def _height_at(self, x):
        plan = self._plan or self.plan()
        if not plan.get("leaves"):
            return self._cut
        width = max(plan["right"] - plan["left"], 1.0)
        return max(0.0, (x - plan["left"]) / width * plan["span"])

    def _near_cut(self, point):
        plan = self._plan or self.plan()
        cut = plan.get("cut")
        return bool(cut) and abs(point.x() - cut["x"]) <= 6

    def _inside_plot(self, point):
        plan = self._plan or self.plan()
        if not plan.get("leaves"):
            return False
        return (plan["left"] - 8 <= point.x() <= plan["right"] + 8
                and plan["cut"]["top"] <= point.y() <= plan["cut"]["bottom"])

    def mousePressEvent(self, event):
        # A click inside the plot moves the cut, which is how the line is found
        # at all; outside it nothing moves, so a stray click changes no grouping.
        if event.button() == Qt.MouseButton.LeftButton and self._inside_plot(event.position()):
            self._dragging = True
            self.set_cut(self._height_at(event.position().x()))
            event.accept()
            return
        super().mousePressEvent(event)

    def mouseMoveEvent(self, event):
        if self._dragging:
            self.set_cut(self._height_at(event.position().x()))
            event.accept()
            return
        self.setCursor(Qt.CursorShape.SplitHCursor if self._near_cut(event.position())
                       else Qt.CursorShape.ArrowCursor)
        super().mouseMoveEvent(event)

    def mouseReleaseEvent(self, event):
        self._dragging = False
        super().mouseReleaseEvent(event)

    def event(self, incoming):
        if incoming.type() == QEvent.Type.ToolTip:
            text = self.describe(incoming.pos())
            if text:
                QToolTip.showText(incoming.globalPos(), text, self)
            else:
                QToolTip.hideText()
            return True
        return super().event(incoming)

    def describe(self, point) -> str:
        """What is under the cursor, in the words the core used for it."""
        plan = self._plan or self.plan()
        if not plan.get("leaves"):
            return ""
        merges = {merge["id"]: merge for merge in self._tree.get("merges", [])}
        for branch in plan["branches"]:
            if abs(branch["x"] - point.x()) <= 5 and \
                    min(branch["segments"][2][1], branch["segments"][2][3]) - 4 <= point.y() \
                    <= max(branch["segments"][2][1], branch["segments"][2][3]) + 4:
                merge = merges.get(branch["id"], {})
                lines = [f"Joined at {_format(merge.get('height', branch['height']))} "
                         f"({self._tree.get('linkage', '')} linkage) · {merge.get('size', 0)} isolates",
                         self.axis_caption()]
                if merge.get("chain_note"):
                    lines.append(merge["chain_note"])
                if merge.get("provisional"):
                    lines.append(f"Provisional: {merge.get('uncompared_within_pairs', 0)} pair(s) "
                                 "inside this cluster were never compared.")
                return "\n".join(filter(None, lines))
        for leaf in plan["leaves"]:
            if abs(leaf["y"] - point.y()) <= max(plan["row_height"] / 2, 3):
                return "\n".join(filter(None, [
                    leaf["name"], leaf["group"] or "No group at this cut", self.axis_caption()]))
        return ""


class DendrogramPanel(QWidget):
    """The dendrogram, its controls, and every caveat the core recorded about it.

    The panel never recomputes a tree; it is handed one. What it does own is the
    sentence beside the picture, and that sentence has to hold two things at
    once: which linkage and which cut produced the groups drawn here, and how
    those relate to the single-linkage "Group ≤" list the page already shows
    elsewhere. Two cluster lists with no explanation is worse than one.
    """

    linkageChanged = Signal(str)
    cutChanged = Signal(float)
    cutAdopted = Signal(int)

    def __init__(self, parent=None):
        super().__init__(parent)
        self._tree = {}
        self._comparison = None
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(6)
        controls = QHBoxLayout()
        controls.setContentsMargins(0, 0, 0, 0)
        controls.addWidget(label("Linkage", "small"))
        self.linkage = QComboBox()
        self.linkage.setProperty("compactCharacters", 10)
        for name in LINKAGE_ORDER:
            self.linkage.addItem(name.capitalize() + " linkage", name)
        self.linkage.setToolTip("\n\n".join(LINKAGES[name] for name in LINKAGE_ORDER))
        self.linkage.currentIndexChanged.connect(
            lambda: self.linkageChanged.emit(self.linkage.currentData()))
        controls.addWidget(self.linkage)
        controls.addWidget(label("Cut ≤", "small"))
        self.cut = QDoubleSpinBox()
        self.cut.setDecimals(2)
        self.cut.setRange(0, 1000000)
        self.cut.setSingleStep(1)
        self.cut.setToolTip(CUT_RULE)
        self.cut.valueChanged.connect(self._cut_typed)
        controls.addWidget(self.cut)
        self.adopt = button("Use this cut as Group ≤", self._adopt_cut)
        self.adopt.setEnabled(False)
        controls.addWidget(self.adopt)
        controls.addStretch(1)
        layout.addLayout(controls)
        self.caption = label("", "small")
        # Ignored width for the same reason as the graph captions: a long caption
        # must never be what pushes this page into a horizontal scroll.
        self.caption.setSizePolicy(QSizePolicy.Policy.Ignored, QSizePolicy.Policy.Preferred)
        layout.addWidget(self.caption)
        split = QSplitter(Qt.Orientation.Vertical)
        split.setChildrenCollapsible(False)
        self.view = DendrogramView()
        self.view.cutChanged.connect(self._cut_dragged)
        split.addWidget(self.view)
        self.notes = QTextBrowser()
        self.notes.setMinimumHeight(90)
        split.addWidget(self.notes)
        split.setStretchFactor(0, 3)
        split.setStretchFactor(1, 1)
        layout.addWidget(split, 1)
        self.show_message("No comparison has been built yet, so there is no tree to draw.")

    # --- what the page hands in ------------------------------------------
    def selected_linkage(self) -> str:
        return self.linkage.currentData() or LINKAGE_ORDER[0]

    def set_linkage(self, name):
        index = self.linkage.findData(str(name))
        if index >= 0 and index != self.linkage.currentIndex():
            self.linkage.blockSignals(True)
            self.linkage.setCurrentIndex(index)
            self.linkage.blockSignals(False)

    def show_message(self, text):
        """Say why there is no tree. An analysis that was not run must say so."""
        self._tree, self._comparison = {}, None
        self.view.clear(text)
        self.caption.setText(text)
        self.caption.setToolTip(text)
        self.adopt.setEnabled(False)
        self.notes.setHtml("<p>" + html.escape(text) + "</p>")

    def show_tree(self, tree, *, names=None, comparison=None, cut=None):
        """Draw ``tree``; ``comparison`` describes the Groups tab list to reconcile with."""
        if not tree:
            return self.show_message(
                "Nothing was clustered for this cohort. That is an analysis that did not run, "
                "not a finding that these isolates are unrelated.")
        self._tree = dict(tree)
        self._comparison = dict(comparison) if comparison else None
        height = self.cut.value() if cut is None else float(cut)
        self.cut.blockSignals(True)
        self.cut.setValue(height)
        self.cut.blockSignals(False)
        self.view.set_tree(self._tree, names=names, cut=height)
        self._refresh_caption()
        self._refresh_notes()

    def tree(self) -> dict:
        """The tree as the core returned it, for whoever has to report or check it."""
        return self._tree

    def cut_height(self) -> float:
        return self.view.cut()

    def notes_text(self) -> str:
        return self.notes.toPlainText()

    # --- the cut ----------------------------------------------------------
    def _cut_typed(self, value):
        self.view.set_cut(float(value), notify=False)
        self._refresh_caption()
        self._refresh_notes()
        self.cutChanged.emit(float(value))

    def _cut_dragged(self, height):
        self.cut.blockSignals(True)
        self.cut.setValue(height)
        self.cut.blockSignals(False)
        self._refresh_caption()
        self._refresh_notes()
        self.cutChanged.emit(float(height))

    def _adopt_cut(self):
        height = self.view.cut()
        if self._can_adopt():
            self.cutAdopted.emit(int(round(height)))

    def _can_adopt(self):
        """Only a single-linkage cut at a whole number is the same rule as Group ≤.

        Group ≤ keeps every pair of isolates that differ by at most that many
        alleles in one component, which is exactly a single-linkage cut at that
        height and is not what a complete or average cut produces. Offering the
        button anywhere else would quietly hand one method's number to another.
        """
        height = self.view.cut()
        return bool(self._tree) and self._tree.get("linkage") == "single" \
            and abs(height - round(height)) <= TOLERANCE

    # --- the words beside the picture -------------------------------------
    def _refresh_caption(self):
        if not self._tree:
            return
        cut = self.view.groups()
        clusters = cut.get("clusters", 0)
        caption = " · ".join(filter(None, [
            self._tree.get("scale_caption", "") or "quantity not recorded",
            f"{self._tree.get('linkage', '')} linkage",
            f"{len(self._tree.get('labels', []))} isolates",
            f"cut at {_format(self.view.cut())}",
            f"{clusters} group(s) of 2 or more"]))
        self.caption.setText(caption)
        self.caption.setToolTip("\n".join(filter(None, [
            caption, self._tree.get("linkage_caption", ""), CUT_RULE, MISSING_DATA_RULE])))
        self.adopt.setEnabled(self._can_adopt())
        self.adopt.setToolTip(
            "Set the page's Group ≤ threshold to this height. A single-linkage cut and the "
            "Group ≤ grouping are the same rule, so the two lists stay one list."
            if self._can_adopt() else
            "Only a single-linkage cut at a whole number is the same rule as Group ≤. A "
            "complete or average cut is a different method, and its height is not a link "
            "threshold.")

    def _refresh_notes(self):
        if not self._tree:
            return
        tree = self._tree
        parts = [f"<h3>{html.escape(str(tree.get('scale_caption') or 'Quantity not recorded'))}</h3>",
                 "<p>" + html.escape(str(tree.get("linkage_caption", ""))) + "</p>",
                 "<p>" + html.escape(CUT_RULE) + "</p>",
                 "<p><b>Uncompared pairs:</b> "
                 + html.escape(str(tree.get("missing_policy_note", ""))) + " "
                 + html.escape(MISSING_DATA_RULE) + "</p>"]
        parts.append(self._blocked_html())
        parts.append(self._provisional_html())
        parts.append(self._chaining_html())
        parts.append(self._agreement_html())
        lonely = tree.get("isolates_without_a_comparable_pair", [])
        if lonely:
            parts.append("<p><b>" + str(len(lonely)) + " isolate(s) had no accepted comparison "
                         "at all:</b> " + html.escape(", ".join(map(str, lonely[:NOTE_LIMIT])))
                         + ("…" if len(lonely) > NOTE_LIMIT else "")
                         + ". They are unassessed, not unrelated.</p>")
        if tree.get("inversions"):
            parts.append("<p><b>" + str(len(tree["inversions"])) + " join(s) sit lower than a join "
                         "beneath them,</b> which happens when a provisional height was measured "
                         "over fewer pairs than it spans. Read those heights as lower bounds.</p>")
        self.notes.setHtml("".join(parts))

    def _blocked_html(self):
        blocked = self._tree.get("blocked_merges", [])
        if not blocked:
            return ""
        rows = []
        for entry in blocked[:NOTE_LIMIT]:
            rows.append("<li>" + html.escape(", ".join(map(str, entry["left_members"][:4])))
                        + (" …" if len(entry["left_members"]) > 4 else "")
                        + " and " + html.escape(", ".join(map(str, entry["right_members"][:4])))
                        + (" …" if len(entry["right_members"]) > 4 else "")
                        + f" — {entry['uncompared_pairs']} pair(s) never compared: "
                        + html.escape("; ".join(entry["reasons"][:2])) + "</li>")
        more = ("<li>… and " + str(len(blocked) - NOTE_LIMIT) + " more.</li>"
                if len(blocked) > NOTE_LIMIT else "")
        return ("<p><b>" + str(len(blocked)) + " pair(s) of clusters could not be joined at any "
                "height,</b> because pairs across them were never compared. They are apart for "
                "want of evidence, not because they were measured as far apart.</p><ul>"
                + "".join(rows) + more + "</ul>")

    def _provisional_html(self):
        merges = [merge for merge in self._tree.get("merges", []) if merge["provisional"]]
        groups = [group for group in self.view.groups().get("groups", []) if group["provisional"]]
        if not merges and not groups:
            return ""
        names = ", ".join(html.escape(str(group["name"])) for group in groups[:NOTE_LIMIT])
        return ("<p><b>Provisional joins:</b> " + str(len(merges)) + " join(s) span pairs that were "
                "never compared, so their heights were measured over fewer pairs than they cover "
                "and are lower bounds. "
                + (f"At this cut that affects: {names}." if groups else
                   "At this cut no group rests on one.") + "</p>")

    def _chaining_html(self):
        chained = [group for group in self.view.groups().get("groups", []) if group["chained"]]
        tree_chained = self._tree.get("chained_merges", [])
        if not chained and not tree_chained:
            return ""
        rows = "".join("<li>" + html.escape(str(group["name"])) + " — "
                       + html.escape(str(group["chain_note"])) + "</li>"
                       for group in chained[:NOTE_LIMIT])
        return ("<p><b>Chaining:</b> " + html.escape(CHAINING_RULE) + "</p>"
                + ("<ul>" + rows + "</ul>" if rows else
                   "<p>" + str(len(tree_chained)) + " join(s) in this tree are chained, but no group "
                   "at this cut is held together by one.</p>"))

    def _agreement_html(self):
        """Name the method behind each cluster list, and say plainly where they differ."""
        cut = self.view.groups()
        mine = [group for group in cut.get("groups", []) if group["status"] == "cluster"]
        here = (f"<p><b>These groups:</b> {self._tree.get('linkage', '')} linkage, cut at "
                f"{_format(self.view.cut())} — {len(mine)} group(s) of 2 or more, over "
                f"{html.escape(str(self._tree.get('scale_caption') or 'an unrecorded quantity'))}.</p>")
        if not self._comparison:
            return here + ("<p>No Group ≤ list was supplied for comparison, so this tab states its "
                           "own method and claims nothing about any other list.</p>")
        threshold = self._comparison.get("threshold")
        theirs = [sorted(map(str, members)) for members in self._comparison.get("groups", [])
                  if len(members) > 1]
        text = here + (f"<p><b>The Groups tab:</b> single linkage, link ≤ {html.escape(str(threshold))}"
                       f" — {len(theirs)} group(s) of 2 or more.</p>")
        placement = {}
        for group in cut.get("groups", []):
            for member in group["members"]:
                placement[str(member)] = frozenset(map(str, group["members"]))
        other = {}
        for members in self._comparison.get("groups", []):
            for member in members:
                other[str(member)] = frozenset(map(str, members))
        moved = [name for name in set(placement) | set(other)
                 if placement.get(name, frozenset([name])) != other.get(name, frozenset([name]))]
        if not moved:
            return text + "<p>Both lists place every isolate the same way.</p>"
        reason = ("The two lists are the same rule at the same height, so they differ only where "
                  "this tree refused to join clusters across a pair that was never compared. That "
                  "refusal is the point: the Groups tab has no edge for such a pair and can still "
                  "reach around it."
                  if self._tree.get("linkage") == "single"
                  and abs(float(self.view.cut()) - float(threshold or 0)) <= TOLERANCE else
                  "They are two different methods, so a difference is expected rather than a fault: "
                  "neither list is the corrected version of the other.")
        return text + ("<p><b>" + str(len(moved)) + " isolate(s) are placed differently by the two "
                       "lists.</b> " + reason + "</p>")
