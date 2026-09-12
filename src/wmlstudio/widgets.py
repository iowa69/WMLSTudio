"""Small native visual components for the desktop workspace."""

import math

from PySide6.QtCore import QPointF, QRectF, Qt, QTimer, Signal
from PySide6.QtGui import QColor, QFont, QPainter, QPen
from PySide6.QtWidgets import (
    QFrame,
    QGraphicsEllipseItem,
    QGraphicsItem,
    QGraphicsScene,
    QGraphicsView,
    QHBoxLayout,
    QLabel,
    QPushButton,
    QVBoxLayout,
    QWidget,
)

from wmlstudio.theme import PALETTE


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
        layout.addWidget(label(title, "muted"))
        self.value = label("0", "metric")
        layout.addWidget(self.value)
        self.hint = label(hint, "small")
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
            painter.setPen(QPen(QColor("#BBD9C8"), 2))
            painter.drawLine(QPointF(-x, y), QPointF(x, y))
            for sign, color in [(1, "#238775"), (-1, "#D4A071")]:
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
        icon.setStyleSheet("font-size: 29px; color: #237E6D; background: #EAF4EE; border-radius: 22px;")
        icon.setFixedSize(46, 46)
        icon.setAlignment(Qt.AlignmentFlag.AlignCenter)
        layout.addWidget(icon)
        text = QVBoxLayout()
        text.addWidget(label("Bring your sequences into focus", "cardTitle"))
        text.addWidget(label("Drop FASTA, FASTQ, or a folder here · .gz and .bz2 supported", "small"))
        layout.addLayout(text, 1)
        layout.addWidget(button("Browse files", self.browseRequested.emit))

    def mousePressEvent(self, event):
        if event.button() == Qt.MouseButton.LeftButton:
            self.browseRequested.emit()

    def dragEnterEvent(self, event):
        if event.mimeData().hasUrls() and all(u.isLocalFile() for u in event.mimeData().urls()):
            event.acceptProposedAction()

    def dropEvent(self, event):
        self.filesDropped.emit([url.toLocalFile() for url in event.mimeData().urls() if url.isLocalFile()])
        event.acceptProposedAction()


class TreeNode(QGraphicsEllipseItem):
    def __init__(self, name, subtitle, color, callback):
        super().__init__(-19, -19, 38, 38)
        self.setBrush(QColor(color))
        self.setPen(QPen(QColor("white"), 3))
        self.setFlags(QGraphicsItem.GraphicsItemFlag.ItemIsMovable |
                      QGraphicsItem.GraphicsItemFlag.ItemSendsGeometryChanges)
        self.setCursor(Qt.CursorShape.OpenHandCursor)
        self.setToolTip(f"{name}\n{subtitle}\nDrag to arrange · scroll to zoom")
        self.callback = callback
        self.setZValue(2)

    def itemChange(self, change, value):
        if change == QGraphicsItem.GraphicsItemChange.ItemPositionHasChanged:
            self.callback()
        return super().itemChange(change, value)


class TreeView(QGraphicsView):
    """Interactive minimum spanning forest, including explicit disconnected nodes."""

    def __init__(self, parent=None):
        self.canvas = QGraphicsScene()
        super().__init__(self.canvas, parent)
        self.setRenderHints(QPainter.RenderHint.Antialiasing | QPainter.RenderHint.TextAntialiasing)
        self.setBackgroundBrush(QColor("#FBFDFB"))
        self.setFrameShape(QFrame.Shape.NoFrame)
        self.setDragMode(QGraphicsView.DragMode.ScrollHandDrag)
        self.setTransformationAnchor(QGraphicsView.ViewportAnchor.AnchorUnderMouse)
        self.nodes = {}
        self.edges = []
        self.labels = {}

    def wheelEvent(self, event):
        factor = 1.12 if event.angleDelta().y() > 0 else 1 / 1.12
        scale = self.transform().m11() * factor
        if 0.1 <= scale <= 5:
            self.scale(factor, factor)
        event.accept()

    def draw_results(self, results, edges, cluster_threshold=1):
        self.canvas.clear()
        self.nodes, self.edges, self.labels = {}, [], {}
        if not results:
            message = self.canvas.addText("Analyse at least two assemblies with the same scheme to compare profiles.")
            message.setDefaultTextColor(QColor("#71877A"))
            self.fit_tree()
            return
        labels = [str(r.get("sample_id", r.get("id", r["sample_name"]))) for r in results]
        parent = {key: key for key in labels}

        def find(key):
            while parent[key] != key:
                key = parent[key]
            return key

        for edge in edges:
            a, b = str(edge["source"]), str(edge["target"])
            if a in parent and b in parent and edge["distance"] <= cluster_threshold:
                parent[find(b)] = find(a)
        groups = list(dict.fromkeys(find(key) for key in labels))
        for i, (key, result) in enumerate(zip(labels, results, strict=True)):
            color = PALETTE[groups.index(find(key)) % len(PALETTE)]
            node = TreeNode(result["sample_name"], f"ST {result.get('st') or 'unassigned'}", color, self.update_edges)
            self.nodes[key] = node
            self.canvas.addItem(node)
            # Fixed deterministic radial layout, freely rearrangeable by the user.
            angle = 2 * math.pi * i / len(results) - math.pi / 2
            radius = max(125, min(420, 34 * len(results)))
            node.setPos(math.cos(angle) * radius * 1.4, math.sin(angle) * radius)
            title = self.canvas.addText(result["sample_name"], QFont("Segoe UI", 10))
            title.setDefaultTextColor(QColor("#2B5041"))
            title.setZValue(3)
            self.labels[key] = title
        for edge in edges:
            a, b = str(edge["source"]), str(edge["target"])
            if a not in self.nodes or b not in self.nodes:
                continue
            pen = QPen(QColor("#ADC2B7"), 2)
            if edge["distance"] > cluster_threshold:
                pen.setStyle(Qt.PenStyle.DashLine)
            line = self.canvas.addLine(0, 0, 0, 0, pen)
            text = self.canvas.addText(str(edge["distance"]), QFont("Segoe UI", 10, QFont.Weight.Bold))
            text.setDefaultTextColor(QColor("#4D7261"))
            text.setToolTip(f"{edge['distance']} differing / {edge['shared_loci']} shared loci")
            self.edges.append((a, b, line, text))
        self.update_edges()
        self.fit_tree()

    def update_edges(self):
        for key, title in self.labels.items():
            pos = self.nodes[key].pos()
            title.setPos(pos.x() - title.boundingRect().width() / 2, pos.y() + 24)
        for a, b, line, text in self.edges:
            start, end = self.nodes[a].pos(), self.nodes[b].pos()
            line.setLine(start.x(), start.y(), end.x(), end.y())
            text.setPos((start + end) / 2 - QPointF(5, 15))

    def fit_tree(self):
        bounds = self.canvas.itemsBoundingRect().adjusted(-45, -35, 45, 35)
        self.canvas.setSceneRect(bounds)
        self.fitInView(bounds, Qt.AspectRatioMode.KeepAspectRatio)

    def save_png(self, path):
        from PySide6.QtGui import QImage
        bounds = self.canvas.itemsBoundingRect().adjusted(-35, -35, 35, 35)
        image = QImage(1800, 1200, QImage.Format.Format_ARGB32)
        image.fill(QColor("#FBFDFB"))
        painter = QPainter(image)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        self.canvas.render(painter, QRectF(0, 0, 1800, 1200), bounds)
        painter.end()
        if not image.save(str(path)):
            raise OSError(f"Could not save image: {path}")
