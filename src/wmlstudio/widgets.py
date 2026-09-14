"""Small native visual components for the desktop workspace."""

import json
import math
from collections import Counter, defaultdict
from pathlib import Path
from xml.etree import ElementTree

from PySide6.QtCore import QPointF, QRectF, QSize, Qt, QTimer, Signal
from PySide6.QtGui import QColor, QFont, QImage, QPainter, QPen
from PySide6.QtWidgets import (
    QColorDialog,
    QFrame,
    QGraphicsEllipseItem,
    QGraphicsItem,
    QGraphicsScene,
    QGraphicsSimpleTextItem,
    QGraphicsView,
    QHBoxLayout,
    QInputDialog,
    QLabel,
    QMenu,
    QPushButton,
    QVBoxLayout,
    QWidget,
)

from wmlstudio.theme import BACKGROUND, INK, MUTED, PALETTE

# Printed under every exported forest, including each half of a side-by-side pair.
GRAPH_SUBTITLE = ("Not a phylogeny or transmission tree · positions are editable "
                  "· edge labels are allele differences")
SIDE_BY_SIDE_NOTE = ("Both panels are layouts of allele differences between profiles. Neither is a "
                     "phylogeny, a time line or a transmission chain, and a line between two isolates "
                     "is similarity, not a proven link.")


def label(text: str, style: str = "", wrap: bool = False) -> QLabel:
    item = QLabel(text)
    item.setTextFormat(Qt.TextFormat.PlainText)
    item.setObjectName(style)
    item.setWordWrap(wrap)
    return item


def button(text: str, callback=None, primary=False) -> QPushButton:
    item = QPushButton(text)
    item.setCursor(Qt.CursorShape.PointingHandCursor)
    if primary:
        item.setObjectName("primary")
    if callback:
        item.clicked.connect(callback)
    return item


def card() -> tuple[QFrame, QVBoxLayout]:
    frame = QFrame()
    frame.setObjectName("card")
    layout = QVBoxLayout(frame)
    layout.setContentsMargins(21, 18, 21, 18)
    layout.setSpacing(12)
    return frame, layout


class Metric(QFrame):
    def __init__(self, title, hint, parent=None):
        super().__init__(parent)
        self.setObjectName("card")
        layout = QVBoxLayout(self)
        layout.setContentsMargins(20, 15, 20, 15)
        layout.addWidget(label(title, "muted", wrap=True))
        self.value = label("0", "metric")
        layout.addWidget(self.value)
        self.hint = label(hint, "small", wrap=True)
        layout.addWidget(self.hint)


class Helix(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setMinimumSize(230, 170)
        self.phase = 0.0
        self.timer = QTimer(self)
        self.timer.timeout.connect(self.advance)
        self.timer.start(60)

    def advance(self):
        if self.isVisible():
            self.phase += 0.024
            self.update()

    def paintEvent(self, event):
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        painter.translate(self.width() / 2, self.height() / 2)
        painter.rotate(-28)
        for i in range(17):
            y = (i - 8) * 11
            phase = i * 0.40 + self.phase
            x = math.sin(phase) * 57
            depth = (math.cos(phase) + 1) / 2
            painter.setPen(QPen(QColor("#365F69"), 2))
            painter.drawLine(QPointF(-x, y), QPointF(x, y))
            for sign, color in [(1, "#48DCC0"), (-1, "#FFB976")]:
                painter.setPen(Qt.PenStyle.NoPen)
                c = QColor(color)
                c.setAlphaF(0.55 + depth * 0.45)
                painter.setBrush(c)
                painter.drawEllipse(QPointF(sign * x, y), 4.8, 4.8)
        painter.end()


class DropZone(QFrame):
    filesDropped = Signal(list)
    browseRequested = Signal()

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setObjectName("drop")
        self.setAcceptDrops(True)
        self.setCursor(Qt.CursorShape.PointingHandCursor)
        layout = QHBoxLayout(self)
        layout.setContentsMargins(24, 20, 24, 20)
        icon = label("＋")
        self._icon = icon
        icon.setStyleSheet("font-size: 29px; color: #83E9D3; background: #204449; border-radius: 22px;")
        icon.setFixedSize(46, 46)
        icon.setAlignment(Qt.AlignmentFlag.AlignCenter)
        layout.addWidget(icon)
        text = QVBoxLayout()
        text.addWidget(label("Bring your sequences into focus", "cardTitle", wrap=True))
        text.addWidget(label("Drop FASTA, FASTQ, or a folder here · .gz and .bz2 supported", "small", wrap=True))
        layout.addLayout(text, 1)
        self._browse = button("Browse files", self.browseRequested.emit)
        layout.addWidget(self._browse)

    def resizeEvent(self, event):
        self._icon.setVisible(self.width() >= 570)
        self._browse.setText("Browse files" if self.width() >= 620 else "Browse")
        super().resizeEvent(event)

    def mousePressEvent(self, event):
        if event.button() == Qt.MouseButton.LeftButton:
            self.browseRequested.emit()

    def dragEnterEvent(self, event):
        if event.mimeData().hasUrls() and all(u.isLocalFile() for u in event.mimeData().urls()):
            event.acceptProposedAction()

    def dropEvent(self, event):
        self.filesDropped.emit([url.toLocalFile() for url in event.mimeData().urls() if url.isLocalFile()])
        event.acceptProposedAction()


def _components(keys, edges):
    """Stable connected components, independent of input ordering."""
    adjacent = {key: set() for key in keys}
    for edge in edges:
        a, b = str(edge["source"]), str(edge["target"])
        if a in adjacent and b in adjacent and a != b:
            adjacent[a].add(b)
            adjacent[b].add(a)
    remaining = set(keys)
    groups = []
    while remaining:
        pending = [min(remaining)]
        group = set()
        while pending:
            key = pending.pop()
            if key not in group:
                group.add(key)
                pending.extend(sorted(adjacent[key] - group, reverse=True))
        remaining.difference_update(group)
        groups.append(sorted(group))
    return groups


def forest_layout(keys, edges):
    """Deterministic force-directed layout; geometric length has no biological meaning."""
    components = _components(sorted(keys), edges)
    positions = {}
    offset_x = offset_y = row_height = 0.0
    # A small disconnected isolate belongs beside the main component when it
    # fits a landscape plot. Starting another tall row needlessly shrank every
    # label in compact native windows.
    row_width = max(1100, math.sqrt(len(keys)) * 210)
    for group in sorted(components, key=lambda items: (-len(items), items)):
        count = len(group)
        index = {key: i for i, key in enumerate(group)}
        links = sorted({tuple(sorted((index[str(e["source"])], index[str(e["target"])])))
                        for e in edges if str(e["source"]) in index and str(e["target"]) in index
                        and e["source"] != e["target"]})
        points = [[math.cos(i * 2.399963) * 80 * math.sqrt(i + 1),
                   math.sin(i * 2.399963) * 80 * math.sqrt(i + 1)] for i in range(count)]
        iterations = 100 if count < 60 else 65
        for step in range(iterations if count > 1 else 0):
            motion = [[0.0, 0.0] for _ in points]
            for i in range(count):
                for j in range(i):
                    dx, dy = points[i][0] - points[j][0], points[i][1] - points[j][1]
                    distance2 = max(dx * dx + dy * dy, 1.0)
                    force = 9500 / distance2
                    motion[i][0] += dx * force
                    motion[i][1] += dy * force
                    motion[j][0] -= dx * force
                    motion[j][1] -= dy * force
            for i, j in links:
                dx, dy = points[i][0] - points[j][0], points[i][1] - points[j][1]
                distance = max(math.hypot(dx, dy), 1.0)
                force = distance / 110
                motion[i][0] -= dx * force
                motion[i][1] -= dy * force
                motion[j][0] += dx * force
                motion[j][1] += dy * force
            temperature = 28 * (1 - step / iterations) + 0.4
            for i, delta in enumerate(motion):
                norm = max(math.hypot(*delta), 0.01)
                points[i][0] += delta[0] / norm * min(norm, temperature)
                points[i][1] += delta[1] / norm * min(norm, temperature)
        # Align each component's principal axis horizontally. Small chains then
        # use the wide native plot panel rather than shrinking into a tall strip.
        mean_x = sum(p[0] for p in points) / count
        mean_y = sum(p[1] for p in points) / count
        centered = [(p[0] - mean_x, p[1] - mean_y) for p in points]
        xx = sum(x * x for x, y in centered)
        yy = sum(y * y for x, y in centered)
        xy = sum(x * y for x, y in centered)
        angle = math.atan2(2 * xy, xx - yy) / 2
        cosine, sine = math.cos(angle), math.sin(angle)
        points = [[x * cosine + y * sine, -x * sine + y * cosine] for x, y in centered]
        left = min(p[0] for p in points)
        top = min(p[1] for p in points)
        width = max(p[0] for p in points) - left + 180
        height = max(p[1] for p in points) - top + 170
        if offset_x and offset_x + width > row_width:
            offset_x, offset_y, row_height = 0, offset_y + row_height + 50, 0
        for key, (x, y) in zip(group, points, strict=True):
            positions[key] = (x - left + offset_x + 80, y - top + offset_y + 65)
        offset_x += width + 55
        row_height = max(row_height, height)
    return positions


_GRAPH_TEXT_SCALE = 100


def set_graph_text_scale(percent):
    """Resize graph labels for a high-resolution screen. Positions never move.

    Only the lettering changes; node coordinates, edge lengths and the distances
    they carry are untouched, so a rescaled tree is the same tree. Exported
    images use the same size as the screen, so a picture matches what was read.
    """
    global _GRAPH_TEXT_SCALE
    _GRAPH_TEXT_SCALE = max(80, min(150, int(percent)))
    return _GRAPH_TEXT_SCALE


def graph_text_scale():
    return _GRAPH_TEXT_SCALE


def _graph_font(size, weight=QFont.Weight.DemiBold):
    font = QFont()
    font.setFamilies(["Segoe UI", "Inter", "DejaVu Sans"])
    font.setPointSize(max(6, round(size * _GRAPH_TEXT_SCALE / 100)))
    font.setWeight(weight)
    return font


class GraphLabel(QGraphicsSimpleTextItem):
    def __init__(self, text, edge=False, callback=None):
        super().__init__(text)
        self.edge = edge
        self.callback = callback
        self.node_key = None
        self.setFont(_graph_font(11 if edge else 14))
        self.setBrush(QColor("#C7D8E8" if edge else INK))
        self.setAcceptedMouseButtons(Qt.MouseButton.LeftButton if callback else Qt.MouseButton.NoButton)
        self.setZValue(1 if edge else 3)

    def paint(self, painter, option, widget=None):
        if self.edge:
            painter.setPen(Qt.PenStyle.NoPen)
            painter.setBrush(QColor(BACKGROUND))
            painter.drawRoundedRect(self.boundingRect(), 5, 5)
        super().paint(painter, option, widget)

    def boundingRect(self):
        bounds = super().boundingRect()
        return bounds.adjusted(-5, -2, 5, 2) if self.edge else bounds

    def mousePressEvent(self, event):
        if self.callback:
            self.callback(event.modifiers())
            event.accept()
        else:
            super().mousePressEvent(event)


class TreeNode(QGraphicsEllipseItem):
    def __init__(self, name, subtitle, color, callback, *, key="", count=1):
        radius = 24 + min(14, 5 * math.log2(count))
        super().__init__(-radius, -radius, radius * 2, radius * 2)
        self.key = key
        self.count = count
        self.highlighted = False
        self.slices = [(str(color), count)]
        self.setBrush(QColor(color))
        self.setPen(QPen(QColor("#0B1220"), 3))
        self.setFlags(QGraphicsItem.GraphicsItemFlag.ItemIsMovable |
                      QGraphicsItem.GraphicsItemFlag.ItemIsSelectable |
                      QGraphicsItem.GraphicsItemFlag.ItemIsFocusable |
                      QGraphicsItem.GraphicsItemFlag.ItemSendsGeometryChanges)
        self.setCursor(Qt.CursorShape.OpenHandCursor)
        self.setAcceptHoverEvents(True)
        self.setToolTip(f"{name}\n{subtitle}\nDrag to arrange · Ctrl-click to select · right-click to edit")
        self.callback = callback
        self.setZValue(2)

    def boundingRect(self):
        return super().boundingRect().adjusted(-8, -8, 8, 8)

    def paint(self, painter, option, widget=None):
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        if self.isSelected() or self.highlighted:
            painter.setBrush(Qt.BrushStyle.NoBrush)
            painter.setPen(QPen(QColor("#FFD58A" if self.highlighted else "#E9F7FF"), 2.5))
            painter.drawEllipse(self.rect().adjusted(-5, -5, 5, 5))
        painter.setPen(QPen(QColor(BACKGROUND), 2))
        start = 90 * 16
        for color, count in self.slices:
            span = round(count / self.count * 360 * 16)
            painter.setBrush(QColor(color))
            if len(self.slices) == 1:
                painter.drawEllipse(self.rect())
            else:
                painter.drawPie(self.rect(), start, span)
            start += span
        if self.count > 1:
            painter.setPen(QColor(BACKGROUND))
            painter.setFont(_graph_font(10, QFont.Weight.Bold))
            painter.drawText(self.rect(), Qt.AlignmentFlag.AlignCenter, str(self.count))

    def itemChange(self, change, value):
        if change == QGraphicsItem.GraphicsItemChange.ItemPositionHasChanged:
            self.callback()
        return super().itemChange(change, value)


class TreeView(QGraphicsView):
    """Editable *presentation* of an allele-distance MST/forest, never a phylogeny.

    Source results and allele calls are not modified by arranging, recoloring,
    relabeling, highlighting, or merging nodes. Stable sample IDs carry edits.
    """

    selectionChanged = Signal(list)
    nodeActivated = Signal(str)
    colorsChanged = Signal(dict)
    layoutChanged = Signal(dict)
    legendChanged = Signal(dict)
    labelsChanged = Signal(dict)
    reportRequested = Signal(list)
    proximityRequested = Signal(str)

    def __init__(self, parent=None):
        self.canvas = QGraphicsScene()
        super().__init__(self.canvas, parent)
        self.canvas.setParent(self)
        self.setRenderHints(QPainter.RenderHint.Antialiasing | QPainter.RenderHint.TextAntialiasing)
        self.setBackgroundBrush(QColor(BACKGROUND))
        self.canvas.setBackgroundBrush(QColor(BACKGROUND))
        self.setFrameShape(QFrame.Shape.NoFrame)
        # Native Windows redraws must invalidate the opaque viewport, including after navigation.
        self.setViewportUpdateMode(QGraphicsView.ViewportUpdateMode.FullViewportUpdate)
        self.setTransformationAnchor(QGraphicsView.ViewportAnchor.AnchorUnderMouse)
        self.setResizeAnchor(QGraphicsView.ViewportAnchor.AnchorViewCenter)
        self.nodes, self.edges, self.labels = {}, [], {}
        self._results, self._source_edges, self._display_edges = {}, [], []
        self._members, self._manual_colors, self._aliases = {}, {}, {}
        self._cluster_groups, self._halos = [], []
        self._group_definitions, self._halo_labels = [], {}
        self.label_fields = ["sample_name", "primary_st"]
        self._label_guides = {}
        self._legend, self._positions = {}, {}
        # A colour agreed with another view, so two forests shown side by side do
        # not give one category two colours. Empty means "decide my own colours".
        self._pinned_legend = {}
        # Optional window hook: context_extension(view, menu, node_key) may append
        # shared workspace actions and returns {QAction: callable} for dispatch.
        self.context_extension = None
        self._palette = [QColor(color).name() for color in PALETTE]
        self.color_by = "cluster"
        self.show_labels = self.show_st = self.show_edge_labels = True
        self.merge_identical = False
        self.show_halos = True
        self.cluster_threshold = 1
        self._highlight = ""
        self._drawing = False
        self._panning = None
        self.set_interaction_mode("select")
        self.canvas.selectionChanged.connect(self._emit_selection)
        self.setToolTip("Allele-distance minimum spanning forest; not a phylogeny.\n"
                        "Drag nodes · Ctrl-click / rectangle-select · middle-drag to pan · scroll to zoom")

    def _emit_selection(self):
        if not self._drawing:
            self.selectionChanged.emit(self.selected_ids())

    def selected_ids(self):
        return sorted(member for key, node in self.nodes.items() if node.isSelected()
                      for member in self._members[key])

    def select_ids(self, sample_ids):
        selected = set(map(str, sample_ids))
        self.canvas.blockSignals(True)
        for key, node in self.nodes.items():
            node.setSelected(bool(selected.intersection(self._members[key])))
        self.canvas.blockSignals(False)
        self._emit_selection()

    def set_interaction_mode(self, mode):
        if mode not in {"select", "pan"}:
            raise ValueError("Interaction mode must be 'select' or 'pan'.")
        self.interaction_mode = mode
        self.setDragMode(QGraphicsView.DragMode.RubberBandDrag if mode == "select"
                         else QGraphicsView.DragMode.ScrollHandDrag)

    def wheelEvent(self, event):
        factor = 1.12 if event.angleDelta().y() > 0 else 1 / 1.12
        if 0.05 <= self.transform().m11() * factor <= 8:
            self.scale(factor, factor)
        event.accept()

    def mousePressEvent(self, event):
        if event.button() == Qt.MouseButton.MiddleButton:
            self._panning = event.position()
            self.setCursor(Qt.CursorShape.ClosedHandCursor)
            event.accept()
        else:
            super().mousePressEvent(event)

    def mouseMoveEvent(self, event):
        if self._panning is not None:
            delta = event.position() - self._panning
            self._panning = event.position()
            self.horizontalScrollBar().setValue(self.horizontalScrollBar().value() - round(delta.x()))
            self.verticalScrollBar().setValue(self.verticalScrollBar().value() - round(delta.y()))
            event.accept()
        else:
            super().mouseMoveEvent(event)

    def mouseReleaseEvent(self, event):
        if event.button() == Qt.MouseButton.MiddleButton and self._panning is not None:
            self._panning = None
            self.unsetCursor()
            event.accept()
        else:
            super().mouseReleaseEvent(event)

    def mouseDoubleClickEvent(self, event):
        node = self.itemAt(event.position().toPoint())
        if isinstance(node, GraphLabel) and node.node_key in self.nodes:
            node = self.nodes[node.node_key]
        if isinstance(node, TreeNode):
            self.nodeActivated.emit(self._members[node.key][0])
        else:
            self.fit_tree()
        event.accept()

    def keyPressEvent(self, event):
        if event.key() == Qt.Key.Key_A and event.modifiers() & Qt.KeyboardModifier.ControlModifier:
            self.select_ids(self._results)
            event.accept()
        elif event.key() == Qt.Key.Key_Escape:
            self.canvas.clearSelection()
            self.highlight("")
            event.accept()
        else:
            super().keyPressEvent(event)

    def contextMenuEvent(self, event):
        node = self.itemAt(event.pos())
        if isinstance(node, GraphLabel) and node.node_key in self.nodes:
            node = self.nodes[node.node_key]
        if isinstance(node, TreeNode) and not node.isSelected():
            self.select_ids(self._members[node.key])
        menu = QMenu(self)
        cluster = menu.addAction("Select this threshold group")
        cluster.setEnabled(isinstance(node, TreeNode))
        report = menu.addAction("Report selected isolates…")
        report.setEnabled(bool(self.selected_ids()))
        proximity = menu.addAction("Inspect nearest neighbours…")
        proximity.setEnabled(isinstance(node, TreeNode))
        menu.addSeparator()
        color = menu.addAction("Set selected node color…")
        color.setEnabled(bool(self.selected_ids()))
        rename = menu.addAction("Edit display label…")
        rename.setEnabled(isinstance(node, TreeNode))
        reset = menu.addAction("Reset selected colors")
        reset.setEnabled(bool(self.selected_ids()))
        menu.addSeparator()
        arrange = menu.addAction("Reset automatic layout")
        fit = menu.addAction("Fit graph")
        extension = {}
        if callable(self.context_extension):
            key = node.key if isinstance(node, TreeNode) else None
            extension = dict(self.context_extension(self, menu, key) or {})
        action = menu.exec(event.globalPos())
        menu.deleteLater()
        if action in extension:
            extension[action]()
        elif action == cluster:
            self.select_cluster(self._members[node.key][0])
        elif action == report:
            self.reportRequested.emit(self.selected_ids())
        elif action == proximity:
            self.proximityRequested.emit(self._members[node.key][0])
        elif action == color:
            chosen = QColorDialog.getColor(QColor(PALETTE[0]), self, "Node color",
                                            QColorDialog.ColorDialogOption.DontUseNativeDialog)
            if chosen.isValid():
                self.set_selected_color(chosen.name())
        elif action == rename:
            key = self._members[node.key][0]
            value, okay = QInputDialog.getText(self, "Display label", "Presentation only; original sample name is retained:",
                                               text=self._aliases.get(key, self._results[key]["sample_name"]))
            if okay:
                self.set_node_label(key, value)
        elif action == reset:
            self.reset_colors(self.selected_ids())
        elif action == arrange:
            self.reset_layout()
        elif action == fit:
            self.fit_tree()

    @staticmethod
    def _merge_signature(result):
        # Do not equate zero differences over shared loci with identical complete genotypes.
        alleles, calls = result.get("alleles", {}), result.get("calls", [])
        digest = result.get("scheme_digest")
        if not digest or not alleles or not calls or result.get("status") == "mixed":
            return None
        if any(call.get("status") != "exact" for call in calls):
            return None
        if len(calls) != len(alleles) or {c.get("locus") for c in calls} != set(alleles):
            return None
        if any(str(call.get("allele")) != str(alleles[call["locus"]]) for call in calls):
            return None
        if any(value is None or str(value) in {"", "0", "-", "?"} for value in alleles.values()):
            return None
        return str(digest), tuple(sorted((str(k), str(v)) for k, v in alleles.items()))

    def draw_results(self, results, edges, cluster_threshold=1, *, groups=None):
        records = {}
        for result in results:
            key = str(result.get("sample_id", result.get("id", result["sample_name"])))
            if key in records:
                raise ValueError(f"Duplicate sample identifier in graph: {key}")
            records[key] = dict(result)
        selected = self.selected_ids()
        self._positions.update({key: (node.pos().x(), node.pos().y()) for key, node in self.nodes.items()})
        self._results = dict(sorted(records.items()))
        self._source_edges = [dict(e) for e in edges]
        self.cluster_threshold = cluster_threshold
        self._group_definitions = [dict(group) for group in groups] if groups is not None else []
        self._drawing = True
        self.canvas.blockSignals(True)
        self.nodes, self.edges, self.labels, self._halos, self._label_guides = {}, [], {}, [], {}
        self._halo_labels = {}
        self.canvas.clear()
        self._members = {}
        signatures = {}
        mapping = {}
        for key, result in self._results.items():
            signature = self._merge_signature(result) if self.merge_identical else None
            representative = signatures.setdefault(signature, key) if signature is not None else key
            mapping[key] = representative
            self._members.setdefault(representative, []).append(key)
        self._display_edges = []
        seen = set()
        for edge in sorted(self._source_edges, key=lambda e: (e["distance"], str(e["source"]), str(e["target"]))):
            a, b = mapping.get(str(edge["source"])), mapping.get(str(edge["target"]))
            pair = tuple(sorted((a, b))) if a is not None and b is not None else None
            if a is None or b is None or a == b or pair in seen:
                continue
            seen.add(pair)
            self._display_edges.append(dict(edge, source=a, target=b))
        self._cluster_groups = _components(self._members, [e for e in self._display_edges
                                                           if e["distance"] <= cluster_threshold])
        if self._group_definitions:
            self._cluster_groups = [sorted({mapping[sid] for sid in group["members"] if sid in mapping})
                                    for group in self._group_definitions]
            self._cluster_groups = [group for group in self._cluster_groups if group]
        else:
            self._group_definitions = [{"name": f"Group {index + 1}", "number": index + 1,
                                       "members": [member for key in group for member in self._members[key]],
                                       "status": "cluster" if len(group) > 1 else "singleton"}
                                      for index, group in enumerate(self._cluster_groups)]
        # Laying out a 500-node forest costs seconds. When every node already has
        # a position (a redraw after a threshold, merge or label change) reuse it;
        # "Reset automatic layout" still calls forest_layout explicitly.
        fresh_positions = (dict(self._positions)
                           if all(key in self._positions for key in self._members)
                           else forest_layout(self._members, self._display_edges))
        for key, members in self._members.items():
            result = self._results[key]
            names = "\n".join(self._results[member]["sample_name"] for member in members)
            node = TreeNode(names, f"ST {result.get('st') or 'unassigned'}", PALETTE[0],
                            self.update_edges, key=key, count=len(members))
            self.nodes[key] = node
            self.canvas.addItem(node)
            node.setPos(*self._positions.get(key, fresh_positions[key]))
            title = GraphLabel("", callback=lambda modifiers, k=key: self._select_label(k, modifiers))
            title.node_key = key
            self.canvas.addItem(title)
            self.labels[key] = title
            guide_pen = QPen(QColor("#36516C"), 0.8, Qt.PenStyle.DotLine)
            guide_pen.setCosmetic(True)
            guide = self.canvas.addLine(0, 0, 0, 0, guide_pen)
            guide.setZValue(0.1)
            guide.setAcceptedMouseButtons(Qt.MouseButton.NoButton)
            guide.setToolTip("Display-label guide, not an allele-distance edge.")
            guide.hide()
            self._label_guides[key] = guide
        for edge in self._display_edges:
            pen = QPen(QColor("#5F819B"), 2)
            pen.setCosmetic(True)
            if edge["distance"] > cluster_threshold:
                pen.setStyle(Qt.PenStyle.DashLine)
            line = self.canvas.addLine(0, 0, 0, 0, pen)
            text = GraphLabel(str(edge["distance"]), edge=True)
            tooltip = (f"{edge['distance']} differing / {edge['shared_loci']} shared loci\n"
                       "Line length is layout only, not evolutionary time or transmission.")
            text.setToolTip(tooltip)
            line.setToolTip(tooltip)
            self.canvas.addItem(text)
            text.setVisible(self.show_edge_labels)
            self.edges.append((edge["source"], edge["target"], line, text))
        for index, group in enumerate(self._cluster_groups):
            definition = self._group_definitions[index]
            fill = QColor(self._palette[index % len(self._palette)])
            border = QColor(fill)
            fill.setAlpha(15)
            border.setAlpha(65)
            halo = self.canvas.addEllipse(QRectF(), QPen(border, 1), fill)
            halo.setZValue(-2)
            halo.setAcceptedMouseButtons(Qt.MouseButton.NoButton)
            halo.setToolTip(f"Single-link group at ≤ {cluster_threshold} allele differences.\n"
                            "A visual grouping, not evidence of an outbreak or transmission.")
            halo.setVisible(self.show_halos)
            self._halos.append((group, halo))
            caption = GraphLabel(str(definition["name"]) + (" · chaining" if definition.get("chained") else ""), edge=True)
            caption.setZValue(-1)
            caption.setVisible(self.show_halos)
            caption.setToolTip(halo.toolTip())
            self.canvas.addItem(caption)
            self._halo_labels[index] = caption
        if not records:
            message = GraphLabel("Analyse at least two assemblies with the same scheme to compare profiles.")
            message.setBrush(QColor(MUTED))
            self.canvas.addItem(message)
        self._drawing = False
        self.canvas.blockSignals(False)
        if self.color_by not in self.available_color_fields():
            self.color_by = "cluster"
        self._refresh_labels()
        self._apply_colors()
        self.highlight(self._highlight)
        self.select_ids(selected)
        self.update_edges()
        self.fit_tree()

    def _redraw(self):
        self.draw_results(list(self._results.values()), self._source_edges, self.cluster_threshold,
                          groups=self._group_definitions)

    def _select_label(self, key, modifiers):
        ids = set(self.selected_ids()) if modifiers & Qt.KeyboardModifier.ControlModifier else set()
        members = set(self._members[key])
        ids = ids - members if members <= ids else ids | members
        self.select_ids(ids)

    def select_cluster(self, sample_id):
        group = next((g for g in self._group_definitions if sample_id in g["members"]), None)
        self.select_ids(group["members"] if group else [sample_id])

    def set_merge_identical(self, enabled):
        self.merge_identical = bool(enabled)
        self._redraw()

    def set_halos_visible(self, visible):
        self.show_halos = bool(visible)
        for _group, halo in self._halos:
            halo.setVisible(self.show_halos)
        for caption in self._halo_labels.values():
            caption.setVisible(self.show_halos)
        self.viewport().update()

    def available_color_fields(self):
        fields = {key for result in self._results.values() for key in self._metadata_values(result)}
        return ["cluster", "st"] + [f"metadata:{key}" for key in sorted(fields)]

    @staticmethod
    def _metadata_values(result):
        values = {}

        def collect(mapping, prefix="", depth=0):
            for key, value in mapping.items():
                if not prefix and key in {'hydra', 'characterization', 'assembly'}:
                    continue
                field = prefix + str(key)
                if isinstance(value, dict) and value and depth < 4:
                    collect(value, field + ".", depth + 1)
                elif isinstance(value, (str, int, float, bool)) or value is None:
                    values[field] = value

        collect(result.get("metadata") or {})
        if not (result.get('metadata') or {}).get('characterization'):
            return values
        from wmlstudio.characterization import current_characterization
        state = current_characterization(result)
        values['characterization.state'] = state['status']
        evidence = state.get('evidence') or {}
        species = evidence.get('species_evidence') or {}
        organism = ' '.join(str(species.get(k) or '') for k in ('genus', 'species', 'subspecies')).strip()
        values['characterization.species'] = (organism + ' · ' + species.get('status', 'not_run')) if organism else state['status']
        for key, source, field in [('virulence', 'hits', 'gene'), ('plasmid_hypotheses', 'replicons', 'gene'),
                                    ('drug_associations', 'associations', 'class')]:
            module = evidence.get(key) or {}
            names = sorted({str(row[field]) for row in module.get(source, []) if isinstance(row, dict) and row.get(field)})
            values['characterization.' + key] = (', '.join(names[:12]) if module.get('status') in {'completed', 'detected', 'not_detected'} and names
                                                  else module.get('status', state['status']))
        return values

    def set_color_by(self, field):
        if field not in self.available_color_fields():
            raise ValueError(f"Unavailable coloring field: {field}")
        self.color_by = field
        self._apply_colors()

    def set_palette(self, colors):
        if not colors or any(not QColor(color).isValid() for color in colors):
            raise ValueError("A palette needs at least one valid color.")
        self._palette = [QColor(color).name() for color in colors]
        self._apply_colors()

    def _category(self, key, clusters):
        result = self._results[key]
        if self.color_by == "cluster":
            return clusters[key]
        if self.color_by == "st":
            st = result.get('primary_st') or result.get('st')
            return f"ST {st}" if st is not None else "ST unassigned"
        value = self._metadata_values(result).get(self.color_by.removeprefix("metadata:"))
        return "Not recorded" if value is None or value == "" else (
            json.dumps(value, ensure_ascii=False, sort_keys=True) if isinstance(value, (dict, list)) else str(value))

    def color_categories(self):
        """The colour categories this view currently shows, for a shared legend."""
        clusters = {member: group["name"] for group in self._group_definitions for member in group["members"]}
        return sorted({self._category(key, clusters) for key in self._results})

    def set_pinned_legend(self, mapping):
        """Force given categories to given colours so two views agree.

        Only categories this view actually shows are affected; an unknown
        category is ignored rather than invented. Recolouring happens only when
        the mapping really changed, which is what keeps a shared legend from
        looping through ``legendChanged``.
        """
        pinned = {str(key): QColor(str(value)).name() for key, value in dict(mapping or {}).items()
                  if QColor(str(value)).isValid()}
        if pinned == self._pinned_legend:
            return False
        self._pinned_legend = pinned
        self._apply_colors()
        return True

    def _apply_colors(self):
        clusters = {member: group["name"] for group in self._group_definitions for member in group["members"]}
        categories = {key: self._category(key, clusters) for key in self._results}
        self._legend = {category: self._palette[index % len(self._palette)]
                        for index, category in enumerate(sorted(set(categories.values())))}
        if self.color_by == "cluster":
            self._legend = {g["name"]: (self._palette[((g.get("number") or 1) - 1) % len(self._palette)]
                                       if g.get("status") == "cluster" else "#73869A")
                            for g in self._group_definitions}
        if self._pinned_legend:
            self._legend = {category: self._pinned_legend.get(category, color)
                            for category, color in self._legend.items()}
        for key, members in self._members.items():
            colors = Counter(self._manual_colors.get(member, self._legend[categories[member]]) for member in members)
            self.nodes[key].slices = sorted(colors.items())
            self.nodes[key].setBrush(QColor(self.nodes[key].slices[0][0]))
            self.nodes[key].update()
        for index, (group, halo) in enumerate(self._halos):
            definition = self._group_definitions[index]
            color = QColor(self._palette[((definition.get("number") or 1) - 1) % len(self._palette)]
                           if definition.get("status") == "cluster" else "#73869A")
            border = QColor(color)
            color.setAlpha(15)
            border.setAlpha(65)
            halo.setBrush(color)
            halo.setPen(QPen(border, 1))
        self.legendChanged.emit(self.legend())
        self.viewport().update()

    def legend(self):
        return dict(self._legend)

    def set_selected_color(self, color):
        self.set_node_colors({key: color for key in self.selected_ids()})

    def set_node_colors(self, colors):
        if any(key not in self._results or not QColor(color).isValid() for key, color in colors.items()):
            raise ValueError("Node colors require known sample IDs and valid colors.")
        self._manual_colors.update({key: QColor(color).name() for key, color in colors.items()})
        self._apply_colors()
        self.colorsChanged.emit(dict(self._manual_colors))

    def reset_colors(self, sample_ids=None):
        if sample_ids is None:
            self._manual_colors.clear()
        else:
            for key in sample_ids:
                self._manual_colors.pop(str(key), None)
        self._apply_colors()
        self.colorsChanged.emit(dict(self._manual_colors))

    def set_node_label(self, sample_id, text):
        key = str(sample_id)
        if key not in self._results:
            raise ValueError(f"Unknown sample ID: {key}")
        text = str(text).strip()
        if len(text) > 120 or "\n" in text or "\r" in text:
            raise ValueError("Display labels must be one line of at most 120 characters.")
        if text:
            self._aliases[key] = text
        else:
            self._aliases.pop(key, None)
        self._refresh_labels()
        self.labelsChanged.emit(dict(self._aliases))

    def _refresh_labels(self):
        for key, item in self.labels.items():
            members = self._members[key]
            result = self._results[key]
            lines = []
            for field in self.label_fields:
                if field == "sample_name":
                    value = self._aliases.get(key, result["sample_name"])
                elif field in {"st", "primary_st"}:
                    if not self.show_st:
                        continue
                    value = "ST " + str(result.get("primary_st", result.get("st")) or "unassigned")
                elif field == "amr_genes":
                    state = result.get('amr_evidence_status')
                    value = "AMR" + (f" ({state})" if state and state != 'current' else '') + ": " + (", ".join(result.get("amr_genes") or []) or "not recorded")
                elif field.startswith("metadata:"):
                    name = field.removeprefix("metadata:")
                    value = name.rsplit(".", 1)[-1] + ": " + str(self._metadata_values(result).get(name) or "not recorded")
                else:
                    value = str(result.get(field) or "not recorded")
                value = str(value).replace("\n", " ")
                lines.append(value if len(value) <= 64 else value[:61] + "…")
            text = "\n".join(lines)
            if len(members) > 1:
                text += f" +{len(members) - 1}"
            item.setText(text)
            item.setVisible(self.show_labels)
            item.setToolTip("\n".join(self._results[member]["sample_name"] for member in members))
        self.update_edges()

    def available_label_fields(self):
        return ["sample_name", "primary_st", "scheme", "amr_genes"] + [
            field for field in self.available_color_fields() if field.startswith("metadata:")]

    def set_label_fields(self, fields):
        fields = list(dict.fromkeys(map(str, fields)))
        if not fields or len(fields) > 8:
            raise ValueError("Choose one to eight graph label fields.")
        if any(field not in {"sample_name", "primary_st", "st", "scheme", "amr_genes"}
               and not field.startswith("metadata:") for field in fields):
            raise ValueError("Unknown graph label field.")
        self.label_fields = fields
        self._refresh_labels()

    def set_labels_visible(self, visible):
        self.show_labels = bool(visible)
        self._refresh_labels()

    def set_show_st(self, visible):
        self.show_st = bool(visible)
        self._refresh_labels()

    def set_edge_labels_visible(self, visible):
        self.show_edge_labels = bool(visible)
        for _a, _b, _line, text in self.edges:
            text.setVisible(self.show_edge_labels)

    def highlight(self, text):
        self._highlight = str(text).strip().casefold()
        matches = []
        for key, members in self._members.items():
            matched = bool(self._highlight) and any(
                self._highlight in " ".join((self._results[member]["sample_name"], self._aliases.get(member, ""),
                                             f"ST {self._results[member].get('st') or 'unassigned'}",
                                             json.dumps(self._results[member].get("metadata") or {}, ensure_ascii=False))).casefold()
                for member in members)
            self.nodes[key].highlighted = matched
            self.nodes[key].update()
            if matched:
                matches.extend(members)
        return sorted(matches)

    def update_edges(self):
        if self._drawing:
            return
        occupied = []
        # The item's bookkeeping bounds include an additional paint-invalidation
        # margin. Treating that margin as the visible node sent ordinary labels
        # several rows away even in a two-isolate graph.
        node_bounds = {key: node.mapRectToScene(node.rect()).adjusted(-7, -7, 7, 7)
                       for key, node in self.nodes.items()}
        for key, title in self.labels.items():
            node = self.nodes[key]
            pos = node.pos()
            bounds = title.boundingRect()
            radius = node.rect().height() / 2
            x = pos.x() - bounds.width() / 2
            options = []
            for row in range(9):
                for shift in (0, -0.55 * bounds.width(), 0.55 * bounds.width()):
                    options.extend([QPointF(x + shift, pos.y() + radius + 15 + row * (bounds.height() + 8)),
                                    QPointF(x + shift, pos.y() - radius - bounds.height() - 15 - row * (bounds.height() + 8))])
            options.sort(key=lambda point: abs(point.x() - x) * 1.2 + abs(point.y() + bounds.height() / 2 - pos.y()))
            chosen = options[0]
            if title.isVisible():
                other_nodes = [rectangle for other_key, rectangle in node_bounds.items() if other_key != key]
                for candidate in options:
                    rectangle = bounds.translated(candidate).adjusted(-4, -3, 4, 3)
                    if not any(rectangle.intersects(other) for other in occupied + other_nodes):
                        chosen = candidate
                        break
                occupied.append(bounds.translated(chosen).adjusted(-4, -3, 4, 3))
            title.setPos(chosen)
            guide = self._label_guides[key]
            rectangle = bounds.translated(chosen)
            target = QPointF(max(rectangle.left(), min(pos.x(), rectangle.right())),
                             max(rectangle.top(), min(pos.y(), rectangle.bottom())))
            guide.setLine(pos.x(), pos.y(), target.x(), target.y())
            guide.setVisible(title.isVisible() and (abs(chosen.x() - x) > 1
                             or abs(target.y() - pos.y()) > radius + 16))
        for a, b, line, text in self.edges:
            start, end = self.nodes[a].pos(), self.nodes[b].pos()
            line.setLine(start.x(), start.y(), end.x(), end.y())
            middle = (start + end) / 2
            text.setPos(middle.x() - text.boundingRect().width() / 2,
                        middle.y() - text.boundingRect().height() / 2)
        for index, (group, halo) in enumerate(self._halos):
            bounds = QRectF()
            for key in group:
                rect = self.nodes[key].sceneBoundingRect().united(self.labels[key].sceneBoundingRect())
                bounds = bounds.united(rect)
            halo.setRect(bounds.adjusted(-28, -24, 28, 24))
            caption = self._halo_labels[index]
            caption.setPos(halo.rect().left() + 12, halo.rect().top() - caption.boundingRect().height() / 2)
        self._positions.update({key: (node.pos().x(), node.pos().y()) for key, node in self.nodes.items()})
        self.layoutChanged.emit({key: list(value) for key, value in self._positions.items() if key in self.nodes})
        self.viewport().update()

    def reset_layout(self):
        positions = forest_layout(self.nodes, self._display_edges)
        self._drawing = True
        for key, node in self.nodes.items():
            node.setPos(*positions[key])
        self._drawing = False
        self.update_edges()
        self.fit_tree()

    def fit_tree(self):
        # Hidden halos and labels still contribute to itemsBoundingRect(). Fit
        # only what the user actually sees, including every visible label.
        bounds = QRectF()
        for item in self.canvas.items():
            if item.isVisible():
                bounds = bounds.united(item.sceneBoundingRect())
        bounds = bounds.adjusted(-25, -22, 25, 22)
        self.canvas.setSceneRect(bounds)
        self.fitInView(bounds, Qt.AspectRatioMode.KeepAspectRatio)
        self.viewport().update()

    def export_state(self):
        return {"version": 1, "positions": {key: list(value) for key, value in self._positions.items()
                                             if key in self._results},
                "colors": {key: value for key, value in self._manual_colors.items() if key in self._results},
                "labels": {key: value for key, value in self._aliases.items() if key in self._results},
                "color_by": self.color_by, "palette": list(self._palette), "show_labels": self.show_labels,
                "show_st": self.show_st, "show_edge_labels": self.show_edge_labels,
                "merge_identical": self.merge_identical, "show_halos": self.show_halos,
                "label_fields": list(self.label_fields)}

    def restore_state(self, state):
        if state == {}:
            state = {"version": 1}
        if not isinstance(state, dict) or state.get("version") != 1:
            raise ValueError("Unsupported graph presentation state version.")
        positions = state.get("positions", {})
        for key, value in positions.items():
            if (not isinstance(value, (list, tuple)) or len(value) != 2
                    or any(not isinstance(v, (int, float)) or not math.isfinite(v) or abs(v) > 1e7 for v in value)):
                raise ValueError(f"Invalid node position: {key}")
        self._positions.update({key: tuple(value) for key, value in positions.items() if key in self._results})
        self.show_labels = bool(state.get("show_labels", True))
        self.show_st = bool(state.get("show_st", True))
        self.show_edge_labels = bool(state.get("show_edge_labels", True))
        self.show_halos = bool(state.get("show_halos", True))
        self.label_fields = list(state.get("label_fields", ["sample_name", "primary_st"]))
        self.merge_identical = bool(state.get("merge_identical", False))
        if state.get("palette"):
            self.set_palette(state["palette"])
        self._manual_colors = {}
        self._aliases = {}
        self.set_node_colors({key: color for key, color in state.get("colors", {}).items() if key in self._results})
        for key, value in state.get("labels", {}).items():
            if key in self._results:
                self.set_node_label(key, value)
        self.color_by = state.get("color_by", "cluster")
        # draw_results snapshots live positions; apply saved coordinates before redrawing.
        self._drawing = True
        for key, node in self.nodes.items():
            if key in positions:
                node.setPos(*positions[key])
        self._drawing = False
        self._redraw()

    def _render(self, painter, width, height, *, title=None, subtitle=None):
        # A saved or reported picture keeps the standard lettering, so the same
        # comparison reads the same on every computer whatever the reader set
        # their own screen to. Only export paths reach here, never a repaint.
        if graph_text_scale() != 100:
            reader_scale = graph_text_scale()
            set_graph_text_scale(100)
            self._redraw()
            try:
                return self._render(painter, width, height, title=title, subtitle=subtitle)
            finally:
                set_graph_text_scale(reader_scale)
                self._redraw()
        bounds = self.canvas.itemsBoundingRect().adjusted(-35, -35, 35, 35)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        painter.fillRect(QRectF(0, 0, width, height), QColor(BACKGROUND))
        painter.setPen(QColor(INK))
        painter.setFont(_graph_font(16))
        painter.drawText(QPointF(35, 38), title or "Allele-distance minimum spanning forest")
        painter.setPen(QColor(MUTED))
        painter.setFont(_graph_font(10, QFont.Weight.Normal))
        painter.drawText(QPointF(35, 63), subtitle or GRAPH_SUBTITLE)
        self.canvas.render(painter, QRectF(0, 90, width, height - 205), bounds)
        painter.setFont(_graph_font(10))
        painter.setPen(QColor(MUTED))
        painter.drawText(QPointF(35, height - 82),
                         f"Color by: {self.color_by} · single-link threshold: {self.cluster_threshold} differences"
                         + (" · manual color overrides present" if self._manual_colors else ""))
        x, y = 35, height - 54
        for category, color in self._legend.items():
            text = category if len(category) <= 40 else category[:37] + "…"
            space = painter.fontMetrics().horizontalAdvance(text) + 43
            if x + space > width - 35:
                x, y = 35, y + 25
            if y > height - 10:
                painter.setPen(QColor(MUTED))
                painter.drawText(QPointF(width - 270, height - 8), "More categories in GraphML export")
                break
            painter.setPen(Qt.PenStyle.NoPen)
            painter.setBrush(QColor(color))
            painter.drawEllipse(QPointF(x + 6, y - 5), 6, 6)
            painter.setPen(QColor(INK))
            painter.drawText(QPointF(x + 21, y), text)
            x += space

    def save_png(self, path):
        return self.save_image(path)

    def save_image(self, path, *, format=None, title=None, subtitle=None):
        image = QImage(1800, 1200, QImage.Format.Format_ARGB32)
        image.fill(QColor(BACKGROUND))
        painter = QPainter(image)
        self._render(painter, 1800, 1200, title=title, subtitle=subtitle)
        painter.end()
        if not image.save(str(path), format):
            raise OSError(f"Could not save image: {path}")

    def save_svg(self, path, *, title=None, subtitle=None):
        from PySide6.QtSvg import QSvgGenerator
        generator = QSvgGenerator()
        generator.setFileName(str(path))
        generator.setSize(QSize(1800, 1200))
        generator.setViewBox(QRectF(0, 0, 1800, 1200))
        generator.setTitle("WMLSTudio allele-distance minimum spanning forest")
        generator.setDescription("Presentation layout only; not a phylogeny or transmission tree.")
        painter = QPainter(generator)
        if not painter.isActive():
            raise OSError(f"Could not save SVG: {path}")
        self._render(painter, 1800, 1200, title=title, subtitle=subtitle)
        painter.end()

    def save_graphml(self, path):
        namespace = "http://graphml.graphdrawing.org/xmlns"
        ElementTree.register_namespace("", namespace)
        def tag(name):
            return f"{{{namespace}}}{name}"
        root = ElementTree.Element(tag("graphml"))
        for name, target, kind in [("label", "node", "string"), ("members", "node", "string"),
                                   ("metadata", "node", "string"), ("x", "node", "double"),
                                   ("y", "node", "double"), ("colors", "node", "string"),
                                   ("allele_differences", "edge", "int"), ("shared_loci", "edge", "int"),
                                   ("interpretation", "graph", "string")]:
            ElementTree.SubElement(root, tag("key"), {"id": name, "for": target, "attr.name": name, "attr.type": kind})
        graph = ElementTree.SubElement(root, tag("graph"), {"id": "allele-distance-forest", "edgedefault": "undirected"})
        ElementTree.SubElement(graph, tag("data"), {"key": "interpretation"}).text = (
            "Allele-distance minimum spanning forest; layout coordinates are presentation only. Not a phylogeny or transmission tree.")
        for key, node in self.nodes.items():
            item = ElementTree.SubElement(graph, tag("node"), {"id": key})
            values = {"label": self._aliases.get(key, self._results[key]["sample_name"]),
                      "members": json.dumps(self._members[key], ensure_ascii=False),
                      "metadata": json.dumps({member: {"sample_name": self._results[member]["sample_name"],
                                                        "st": self._results[member].get("st"),
                                                        "scheme_digest": self._results[member].get("scheme_digest"),
                                                        "metadata": self._results[member].get("metadata", {})}
                                               for member in self._members[key]}, ensure_ascii=False, sort_keys=True),
                      "x": node.pos().x(), "y": node.pos().y(), "colors": json.dumps(node.slices)}
            for name, value in values.items():
                ElementTree.SubElement(item, tag("data"), {"key": name}).text = str(value)
        for index, edge in enumerate(self._display_edges):
            item = ElementTree.SubElement(graph, tag("edge"), {"id": f"e{index}", "source": edge["source"], "target": edge["target"]})
            for name, value in (("allele_differences", edge["distance"]), ("shared_loci", edge["shared_loci"])):
                ElementTree.SubElement(item, tag("data"), {"key": name}).text = str(value)
        ElementTree.ElementTree(root).write(path, encoding="utf-8", xml_declaration=True)

    def save_newick(self, path):
        """Export rooted MST topology for interchange, explicitly not a phylogeny.

        Every isolate is a labeled zero-length leaf attached to its observed
        genotype vertex. MST edge lengths are raw allele differences. Distinct
        disconnected components are separate Newick trees, never invented links.
        """
        adjacent = defaultdict(list)
        for edge in self._display_edges:
            adjacent[edge["source"]].append((edge["target"], edge["distance"]))
            adjacent[edge["target"]].append((edge["source"], edge["distance"]))
        visited = set()

        def quote(value):
            return "'" + str(value).replace("'", "''").replace("\n", " ").replace("\r", " ") + "'"

        def subtree(key):
            visited.add(key)
            branches = [quote(self._results[member]["sample_name"] + " [" + member + "]") + ":0"
                        for member in self._members[key]]
            for other, distance in sorted(adjacent[key]):
                if other not in visited:
                    branches.append(subtree(other) + ":" + str(distance))
            return "(" + ",".join(branches) + ")"

        trees = ["[Allele-distance MST topology; not a phylogeny; root arbitrary; disconnected components separate]"]
        for key in sorted(self.nodes):
            if key not in visited:
                trees.append(subtree(key) + ";")
        Path(path).write_text("\n".join(trees) + "\n", encoding="utf-8")


def render_side_by_side(left, right, *, left_title, right_title, left_subtitle=None,
                        right_subtitle=None, headline, note=None, width=1800, height=1200,
                        header=76):
    """Draw two forests into one image, each keeping its own title and legend.

    The two panels are separate measurements, so neither borrows the other's
    caption: the headline states whether they share a reference, a link threshold
    and a minimum shared-locus fraction, and ``note`` carries the interpretation
    limit into the picture itself, where it cannot be cropped away from the trees
    the way a surrounding caption can.
    """
    image = QImage(int(width) * 2, int(height) + int(header), QImage.Format.Format_ARGB32)
    image.fill(QColor(BACKGROUND))
    painter = QPainter(image)
    try:
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        painter.setPen(QColor(INK))
        painter.setFont(_graph_font(18))
        painter.drawText(QPointF(35, 34), str(headline))
        painter.setPen(QColor(MUTED))
        painter.setFont(_graph_font(10, QFont.Weight.Normal))
        painter.drawText(QPointF(35, 58), str(note or SIDE_BY_SIDE_NOTE))
        for index, (view, title, subtitle) in enumerate(
                ((left, left_title, left_subtitle), (right, right_title, right_subtitle))):
            painter.save()
            painter.translate(index * width, header)
            view._render(painter, width, height, title=title, subtitle=subtitle)
            painter.restore()
        painter.setPen(QPen(QColor(MUTED), 2))
        painter.drawLine(int(width), int(header), int(width), int(header) + int(height))
    finally:
        painter.end()
    return image
