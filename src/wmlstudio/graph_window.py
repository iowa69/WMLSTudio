"""One distance forest in its own window: movable, editable, self-describing.

A window draws a *copy* of a panel's forest — the same records, edges, threshold,
groups and scale in a view of its own — so dragging, recoloring or relabelling
here never moves a node on the page behind it while the window is open. What
happens to that arrangement afterwards is the host's decision: the window hands
it out on :attr:`GraphWindow.stateChanged` and a host that keeps it says so
through :meth:`GraphWindow.keep_arrangement_in`, so the window never implies a
save nobody is performing. Several windows can stand side by side, which is why each
one states its own typing kind, reference, target count, cohort, link threshold
and the time it was opened: a "3" on an edge of a classical MLST forest, a "3" on
an edge of a cgMLST forest and a "3" on an edge of a SKA2 split k-mer SNP forest
are three different quantities, and two windows must never be mistaken for one
another. A window takes those words from the scale it was opened on rather than
assuming allele typing, so a SNP window never says "alleles" anywhere.

Nothing here reads a sequence, recalculates a distance or edits a call. Arranging
and coloring are presentation, and the window says so where it cannot be cropped
away from the picture.
"""

from dataclasses import dataclass, replace
from datetime import datetime
from pathlib import Path

from PySide6.QtCore import QSize, Qt, QTimer, Signal
from PySide6.QtGui import QColor
from PySide6.QtWidgets import (
    QCheckBox,
    QColorDialog,
    QComboBox,
    QFileDialog,
    QHBoxLayout,
    QInputDialog,
    QMainWindow,
    QVBoxLayout,
    QWidget,
)

from wmlstudio import display
from wmlstudio.theme import PALETTE
from wmlstudio.ui_common import FlowLayout
from wmlstudio.widgets import TreeView, button, card, label

# A detached graph may be much smaller than the workspace window: it holds one
# picture, not seven tabs. The ceiling is the shared one, so a graph window and
# the workspace agree about what this screen can show.
MINIMUM_SIZE = (760, 520)
DEFAULT_SIZE = (1180, 820)

# The three statements that must travel with the picture, in the window and in
# every image saved from it.
_OVERLAP = "Too little overlap is an excluded pair, never a zero distance."
HONESTY = (
    "Groups are single linkage: a chain of small differences can join two isolates whose own "
    "distance is larger than the threshold.",
    "Every distance is counted over the loci the two isolates actually share. " + _OVERLAP,
    "This is a layout of allele differences. It is not a phylogeny, not a time line and not a "
    "transmission chain, and line length carries no meaning.",
)


def honesty_lines(identity=None):
    """The same three statements, in the words of whatever this window measures.

    An identity that records no vocabulary of its own is an allele-distance forest
    and gets :data:`HONESTY` unchanged. One that names its own difference word and
    its own denominator states those instead, because "a layout of allele
    differences" is a false description of a SNP forest and the sentence that
    keeps unshared sequence from reading as zero has to be true of the quantity it
    is printed beside.
    """
    note = str(getattr(identity, "denominator_note", "") or "")
    difference = str(getattr(identity, "difference_word", "") or "allele differences")
    return (HONESTY[0], f"{note} {_OVERLAP}" if note else HONESTY[1],
            f"This is a layout of {difference}. It is not a phylogeny, not a time line and not a "
            "transmission chain, and line length carries no meaning.")
EDIT_NOTE = ("Arranging, coloring and renaming change this picture only. Allele calls, distances "
             "and group membership are untouched.")
# Said only by a host that has connected stateChanged and really does keep what it is
# handed. The signal is emitted whether or not anything listens, so a window that
# promised "kept" on its own would be claiming a save that may never happen. The
# sentence names what travels and what does not, because a reader who has carefully
# zoomed in would otherwise expect the zoom back.
KEEP_NOTE = ("Arranging, coloring and renaming are kept: they are applied to {where} as you make "
             "them and once more when this window closes, and saved with the project. Node "
             "positions, colors, display labels and which details are shown travel back; this "
             "window's own zoom and pan do not. Allele calls, distances and group membership are "
             "untouched.")
EXPORT_FORMATS = (("PNG image", "png"), ("JPEG image", "jpg"), ("SVG vector", "svg"),
                  ("GraphML nodes and edges", "graphml"), ("Newick topology", "nwk"))

# Remembered window sizes for this session when no settings file was handed in,
# keyed by typing kind so a cgMLST window and an MLST window keep their own size.
_SESSION_SIZES = {}


def _size_key(kind):
    return f"graph_window/{kind or 'unclassified'}/size"


def _parse_size(token):
    try:
        width, height = (int(part) for part in str(token).lower().split("x"))
    except (AttributeError, TypeError, ValueError):
        return None
    return (width, height)


def _sides(size):
    """(width, height) from a QSize or a plain pair."""
    width, height = getattr(size, "width", None), getattr(size, "height", None)
    return (int(width() if callable(width) else size[0]),
            int(height() if callable(height) else size[1]))


def clamp_size(width, height, available_size=None):
    """A graph-window size inside the supported range, and inside this screen when known."""
    width = max(MINIMUM_SIZE[0], min(display.MAXIMUM_WINDOW_WIDTH, int(width)))
    height = max(MINIMUM_SIZE[1], min(display.MAXIMUM_WINDOW_HEIGHT, int(height)))
    if available_size is not None:
        screen_width, screen_height = _sides(available_size)
        if screen_width >= MINIMUM_SIZE[0]:
            width = min(width, screen_width)
        if screen_height >= MINIMUM_SIZE[1]:
            height = min(height, screen_height)
    return width, height


def _clock(value=None):
    """A plain minute stamp for a window title, from a datetime, an ISO string or text."""
    if value is None:
        return datetime.now().strftime("%Y-%m-%d %H:%M")
    if isinstance(value, datetime):
        return value.strftime("%Y-%m-%d %H:%M")
    text = str(value).strip()
    try:
        return datetime.fromisoformat(text).astimezone().strftime("%Y-%m-%d %H:%M")
    except ValueError:
        return text


@dataclass(frozen=True)
class GraphIdentity:
    """What one window is showing, in the words the window and its exports use.

    Everything here is copied from the comparison that produced the forest. None
    of it is inferred: an unrecorded scheme stays "reference not recorded" and an
    unrecorded target count stays absent rather than being guessed from how many
    loci happened to be called.
    """

    kind: str = "unclassified"
    title: str = "Unclassified typing"
    scheme: str = ""
    targets: int = 0
    target_word: str = "loci"
    cohort: str = ""
    threshold: int | None = None
    created: str = ""
    separation: str = ""
    note: str = ""
    # What this window's numbers are, in its own words. The defaults are the
    # allele-typing vocabulary every window used before a SNP forest could be
    # opened, so a scale that records none of these reads exactly as it did.
    difference_word: str = "allele differences"
    distance_phrase: str = "allele-distance"
    caption_text: str = ""
    denominator_note: str = ""

    @classmethod
    def from_scale(cls, scale=None, *, kind=None, cohort="", threshold=None, created=None, note=""):
        """Build an identity from :func:`wmlstudio.investigation.typing_scale` output."""
        if isinstance(scale, cls):
            return replace(scale, cohort=cohort or scale.cohort,
                           threshold=scale.threshold if threshold is None else int(threshold),
                           created=_clock(created) if created is not None else scale.created,
                           note=note or scale.note)
        scale = dict(scale or {})
        return cls(kind=str(kind or scale.get("kind") or "unclassified"),
                   title=str(scale.get("title") or "Unclassified typing"),
                   scheme=str(scale.get("scheme") or ""),
                   targets=int(scale.get("targets") or 0),
                   target_word=str(scale.get("target_word") or "loci"),
                   cohort=str(cohort or ""),
                   threshold=None if threshold is None else int(threshold),
                   created=_clock(created), separation=str(scale.get("separation") or ""),
                   note=str(note or ""),
                   difference_word=str(scale.get("difference_word") or "allele differences"),
                   distance_phrase=str(scale.get("distance_phrase") or "allele-distance"),
                   caption_text=str(scale.get("caption") or ""),
                   denominator_note=str(scale.get("denominator_note") or ""))

    def caption(self):
        """'cgMLST · kpneumoniae · 2358 targets', naming the quantity this window shows.

        A scale that wrote its own caption is trusted with it: a forest whose
        denominator is per pair has no cohort target count, and "target count not
        recorded" would describe that as a gap rather than as the truth.
        """
        return self.caption_text or " · ".join(
            [self.title, self.scheme or "reference not recorded",
             f"{self.targets} {self.target_word}" if self.targets else "target count not recorded"])

    def threshold_words(self):
        if self.threshold is None:
            return "link threshold not recorded"
        # A negative threshold is this application's way of saying no cutoff was
        # justified. It groups nothing, and must never be printed as a number.
        if self.threshold < 0:
            return "no link threshold set"
        return (f"link ≤ {self.threshold} of {self.targets} {self.target_word}" if self.targets
                else f"link ≤ {self.threshold} {self.difference_word}")

    def cohort_words(self):
        return self.cohort or "cohort not named"

    def window_title(self):
        """The distinguishing facts first, so two open windows differ in the task bar."""
        return " · ".join(filter(None, [f"{self.title} forest", self.threshold_words(),
                                        self.cohort_words(), self.scheme, self.created]))

    def summary(self):
        return " · ".join(filter(None, [self.cohort_words(), self.threshold_words(),
                                        f"opened {self.created}" if self.created else "",
                                        self.note]))

    def export_title(self):
        return f"{self.caption()} · {self.distance_phrase} minimum spanning forest"

    def export_subtitle(self):
        """The line printed under every picture saved from this window."""
        return " · ".join(filter(None, [
            self.cohort_words(),
            self.threshold_words() + (", single linkage" if (self.threshold or 0) >= 0 else ""),
            f"opened {self.created}" if self.created else "",
            "not a phylogeny or transmission chain"]))

    def as_dict(self):
        return {"kind": self.kind, "title": self.title, "scheme": self.scheme,
                "targets": self.targets, "target_word": self.target_word, "cohort": self.cohort,
                "threshold": self.threshold, "created": self.created, "note": self.note,
                "difference_word": self.difference_word, "distance_phrase": self.distance_phrase}


def _field_title(field):
    """A coloring field in the words a microbiologist reading the window would use."""
    if field == "cluster":
        return "Single-link groups"
    if field == "st":
        return "Sequence type"
    name = str(field).removeprefix("metadata:").rsplit(".", 1)[-1].replace("_", " ")
    return (name[:1].upper() + name[1:] if name else field) + " (recorded detail)"


class GraphWindow(QMainWindow):
    """A detached, resizable forest: the same evidence, its own arrangement.

    The window owns its view and its copy of the data, so opening two windows on
    one comparison — at two thresholds, or one baseline beside one current tree —
    gives two pictures that cannot overwrite each other's layout, and no evidence
    on the page it came from is changed by anything done here.
    """

    closed = Signal()
    stateChanged = Signal(dict)
    selectionChanged = Signal(list)
    nodeActivated = Signal(str)
    reportRequested = Signal(list)
    proximityRequested = Signal(str)

    def __init__(self, identity=None, *, view=None, parent=None, preferences=None):
        super().__init__(parent)
        self.identity = (identity if isinstance(identity, GraphIdentity)
                         else GraphIdentity.from_scale(identity))
        self.preferences = preferences
        self.view = view if view is not None else TreeView()
        # Assigned by the host so an export can never overwrite an input file.
        self.check_output = None
        # The forest is fitted to the window while the reader has not placed the
        # view themselves; their own zoom or pan ends that, and "Fit" resumes it.
        self.follow_window = True
        self._contents = {"results": [], "edges": []}
        # What the status bar falls back to between messages. A host that keeps this
        # window's arrangement replaces it through keep_arrangement_in().
        self._edit_note = EDIT_NOTE
        self.setWindowFlag(Qt.WindowType.Window, True)
        self.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose, False)
        self.setMinimumSize(QSize(*MINIMUM_SIZE))
        self.resize(*self.remembered_size())
        body = QWidget()
        body.setObjectName("page")
        self.setCentralWidget(body)
        layout = QVBoxLayout(body)
        layout.setContentsMargins(16, 14, 16, 12)
        layout.setSpacing(10)
        self._build_header(layout)
        self._build_controls(layout)
        # The picture keeps a usable height however small the window is made; the
        # honesty lines sit in one block so they cost one gap, not three.
        self.view.setMinimumHeight(220)
        layout.addWidget(self.view, 1)
        self.honesty = label("\n".join(honesty_lines(self.identity)), "small", wrap=True)
        layout.addWidget(self.honesty)
        self._fit_timer = QTimer(self)
        self._fit_timer.setSingleShot(True)
        self._fit_timer.timeout.connect(self._fit_to_window)
        self._state_timer = QTimer(self)
        self._state_timer.setSingleShot(True)
        self._state_timer.timeout.connect(self._emit_state)
        self.view.viewAdjusted.connect(self._stop_following)
        self.view.selectionChanged.connect(self._selection_changed)
        self.view.nodeActivated.connect(self.nodeActivated)
        self.view.reportRequested.connect(self.reportRequested)
        self.view.proximityRequested.connect(self.proximityRequested)
        for signal in (self.view.layoutChanged, self.view.colorsChanged,
                       self.view.labelsChanged):
            signal.connect(lambda _value: self._state_timer.start(120))
        self.apply_identity(self.identity)
        self.refresh_controls()
        self.say(self._edit_note)

    # --- construction --------------------------------------------------------
    def _build_header(self, layout):
        frame, inner = card()
        inner.setContentsMargins(16, 12, 16, 12)
        inner.setSpacing(6)
        top = QHBoxLayout()
        self.kind_badge = label(self.identity.title, "badge")
        top.addWidget(self.kind_badge)
        self.headline = label(self.identity.caption(), "cardTitle", wrap=True)
        top.addWidget(self.headline, 1)
        inner.addLayout(top)
        self.summary = label(self.identity.summary(), "muted", wrap=True)
        inner.addWidget(self.summary)
        self.separation = label(self.identity.separation, "small", wrap=True)
        self.separation.setVisible(bool(self.identity.separation))
        inner.addWidget(self.separation)
        # Where an arrangement made here ends up, in the window rather than only in
        # a status message a later click would replace. Hidden until a host says so.
        self.keep_note = label("", "small", wrap=True)
        self.keep_note.setVisible(False)
        inner.addWidget(self.keep_note)
        layout.addWidget(frame)

    def _build_controls(self, layout):
        strip = FlowLayout()
        strip.addWidget(button("Fit whole forest", self.fit))
        strip.addWidget(button("Reset arrangement", self.reset_layout))
        self.mode_button = button("Pan with the mouse", self.toggle_mode)
        self.mode_button.setToolTip("Switch between selecting isolates and panning the picture. "
                                    "Nodes can always be dragged; middle-drag always pans.")
        strip.addWidget(self.mode_button)
        self.color_field = QComboBox()
        self.color_field.setToolTip("Which recorded detail decides node color.")
        self.color_field.currentIndexChanged.connect(self._color_field_chosen)
        strip.addWidget(self.color_field)
        strip.addWidget(button("Color selected…", self.color_selected))
        strip.addWidget(button("Color whole group…", self.color_group))
        strip.addWidget(button("Edit display label…", self.rename_selected))
        strip.addWidget(button("Reset colors", self.reset_colors))
        strip.addWidget(button("Save picture…", self.export_dialog, primary=True))
        layout.addLayout(strip)
        toggles = FlowLayout()
        self.toggles = {}
        for name, text, tip in (
                ("labels", "Isolate labels", "Show the display label beside each node."),
                ("st", "Sequence types", "Include the sequence type in each label."),
                ("edges", "Allele differences on edges",
                 "Show how many alleles differ, over the shared loci in the tooltip."),
                ("halos", "Single-link groups",
                 "Outline the groups formed at this window's threshold. A grouping, not an outbreak."),
                ("merge", "Merge identical complete profiles",
                 "Draw one node for isolates whose complete profiles are identical. "
                 "Incomplete or ambiguous profiles are never merged.")):
            box = QCheckBox(text)
            box.setToolTip(tip)
            box.toggled.connect(lambda checked, key=name: self._toggled(key, checked))
            self.toggles[name] = box
            toggles.addWidget(box)
        layout.addLayout(toggles)

    # --- contents ------------------------------------------------------------
    @classmethod
    def from_view(cls, source, identity=None, *, parent=None, preferences=None, view=None,
                  cohort="", created=None, note=""):
        """Open a window on the same forest a panel is showing, arrangement included.

        The data and the presentation state are copied; the panel's own view is
        never reparented, so the page keeps its graph while the window is open.
        """
        contents = source.graph_contents()
        if identity is None:
            identity = GraphIdentity.from_scale(contents.get("scale"),
                                                threshold=contents.get("cluster_threshold"),
                                                cohort=cohort, created=created, note=note)
        window = cls(identity, view=view, parent=parent, preferences=preferences)
        window.set_contents(contents)
        window.view.restore_state(source.export_state())
        window.refresh_controls()
        return window

    def set_contents(self, contents):
        """Draw a copy of a forest, as :meth:`wmlstudio.widgets.TreeView.graph_contents` gives it."""
        self._contents = {"results": list(dict(contents or {}).get("results") or []),
                          "edges": list(dict(contents or {}).get("edges") or [])}
        self.follow_window = True
        self.view.show_contents(contents)
        self.refresh_controls()
        self.say(self._edit_note)
        return self

    def apply_identity(self, identity):
        """Restate what this window shows; nothing is redrawn and no number changes."""
        self.identity = (identity if isinstance(identity, GraphIdentity)
                         else GraphIdentity.from_scale(identity))
        difference = self.identity.difference_word
        self.setWindowTitle(self.identity.window_title())
        self.kind_badge.setText(self.identity.title)
        self.headline.setText(self.identity.caption())
        self.summary.setText(self.identity.summary())
        self.separation.setText(self.identity.separation)
        self.separation.setVisible(bool(self.identity.separation))
        # The controls name this window's own quantity too: a checkbox offering
        # "allele differences on edges" over a SNP forest is simply false.
        self.toggles["edges"].setText(f"{difference[:1].upper()}{difference[1:]} on edges")
        self.toggles["edges"].setToolTip(f"Show how many {difference} separate two isolates, over "
                                         "the sequence they share, which the tooltip states.")
        self.honesty.setText("\n".join(honesty_lines(self.identity)))
        return self.identity

    def refresh_controls(self):
        """Match the controls to the view, without letting that count as an edit."""
        fields = self.view.available_color_fields()
        self.color_field.blockSignals(True)
        self.color_field.clear()
        for field in fields:
            self.color_field.addItem(_field_title(field), field)
        index = self.color_field.findData(self.view.color_by)
        self.color_field.setCurrentIndex(max(0, index))
        self.color_field.blockSignals(False)
        for key, checked in (("labels", self.view.show_labels), ("st", self.view.show_st),
                             ("edges", self.view.show_edge_labels), ("halos", self.view.show_halos),
                             ("merge", self.view.merge_identical)):
            box = self.toggles[key]
            box.blockSignals(True)
            box.setChecked(bool(checked))
            box.blockSignals(False)
        self.mode_button.setText("Select isolates" if self.view.interaction_mode == "pan"
                                 else "Pan with the mouse")

    def refresh_theme(self):
        """Repaint the forest in the theme that is active now.

        The view notices a new theme by itself on its next repaint; this is for a
        host that has just changed the theme and wants every open window to follow
        at once, including one that is behind another and not repainting. Nothing
        is recalculated: the same nodes, edges and numbers on another ground.
        """
        self.view.refresh_theme()

    # --- arrangement ---------------------------------------------------------
    def fit(self):
        self.follow_window = True
        self.view.fit_tree()

    def reset_layout(self):
        self.view.reset_layout()
        self.say("Automatic arrangement restored. The distances it lays out did not change.")

    def toggle_mode(self):
        self.view.set_interaction_mode("select" if self.view.interaction_mode == "pan" else "pan")
        self.refresh_controls()

    def _stop_following(self):
        self.follow_window = False

    def _fit_to_window(self):
        if self.follow_window:
            self.view.fit_tree()

    def resizeEvent(self, event):
        super().resizeEvent(event)
        # The view is refitted after the resize settles, not during it: fitting on
        # every intermediate width makes a large forest crawl.
        if getattr(self, "_fit_timer", None) is not None:
            self._fit_timer.start(90)

    def _toggled(self, key, checked):
        {"labels": self.view.set_labels_visible, "st": self.view.set_show_st,
         "edges": self.view.set_edge_labels_visible, "halos": self.view.set_halos_visible,
         "merge": self.view.set_merge_identical}[key](checked)
        self._state_timer.start(120)

    def _color_field_chosen(self, index):
        field = self.color_field.itemData(index)
        if field and field in self.view.available_color_fields():
            self.view.set_color_by(field)
            self._state_timer.start(120)

    # --- editing -------------------------------------------------------------
    def selected_ids(self):
        return self.view.selected_ids()

    def _selection_changed(self, ids):
        self.selectionChanged.emit(list(ids))
        self.say((f"{len(ids)} isolate(s) selected. " + self._edit_note) if ids else self._edit_note)

    def _chosen_color(self):
        color = QColorDialog.getColor(QColor(PALETTE[0]), self, "Node color",
                                      QColorDialog.ColorDialogOption.DontUseNativeDialog)
        return color.name() if color.isValid() else ""

    def color_selected(self):
        ids = self.selected_ids()
        if not ids:
            return self.say("Select one or more isolates first.")
        color = self._chosen_color()
        if color:
            self.view.set_node_colors({key: color for key in ids})
            self.say(f"{len(ids)} isolate(s) recolored. Color is presentation, never evidence.")

    def color_group(self):
        """Color every isolate in the single-link group holding the selection."""
        ids = self.selected_ids()
        if not ids:
            return self.say("Select an isolate in the group you want to color.")
        members = sorted({member for key in ids for member in self.view.group_members(key)})
        color = self._chosen_color()
        if color:
            self.view.set_node_colors({key: color for key in members})
            self.say(f"{len(members)} isolate(s) in {len(ids)} selected group(s) recolored. "
                     "A group is a single-link grouping, not an outbreak assignment.")

    def rename_selected(self):
        # One node, not one isolate: a node that merges identical complete
        # profiles carries several isolates and still has one display label.
        keys = sorted(key for key, node in self.view.nodes.items() if node.isSelected())
        if len(keys) != 1:
            return self.say("Select exactly one isolate to relabel.")
        text, okay = QInputDialog.getText(
            self, "Display label", "Presentation only; the original sample name is retained:",
            text=self.view.node_label(keys[0]))
        if okay:
            self.view.set_node_label(keys[0], text)
            self.say("Display label changed. The stored sample name is unchanged.")

    def reset_colors(self):
        ids = self.selected_ids()
        self.view.reset_colors(ids or None)
        self.say(f"Colors reset for {len(ids)} isolate(s)." if ids else "All manual colors reset.")

    # --- saved presentation --------------------------------------------------
    def presentation_state(self):
        """This window's arrangement, colors and labels, with what it is showing."""
        return {"identity": self.identity.as_dict(), "view": self.view.export_state()}

    def restore_presentation(self, state):
        """Put a saved arrangement back; unknown isolates are ignored, never invented."""
        state = dict(state or {})
        self.view.restore_state(state.get("view", state))
        self.refresh_controls()
        self.say("Saved arrangement restored.")

    def _emit_state(self):
        self.stateChanged.emit(self.view.export_state())

    def keep_arrangement_in(self, where):
        """Say in the window that a host is keeping what is arranged here.

        Only the host that connected :attr:`stateChanged` knows whether the state it
        is handed is actually applied and stored, so the promise is the host's to
        make and its words for where the arrangement lands go in the sentence.
        """
        self._edit_note = KEEP_NOTE.format(where=str(where).strip() or "the page it came from")
        self.keep_note.setText(self._edit_note)
        self.keep_note.setVisible(True)
        self.say(self._edit_note)
        return self._edit_note

    # --- exports -------------------------------------------------------------
    def export_dialog(self):
        """Ask where to save, then write it with this window's own identity on it."""
        filters = ";;".join(f"{name} (*.{suffix})" for name, suffix in EXPORT_FORMATS)
        suggested = f"{self.identity.kind}-forest.png"
        path, chosen = QFileDialog.getSaveFileName(self, "Save this forest", suggested, filters)
        if not path:
            return None
        # A typed name without an extension gets the format that was chosen in the
        # dialog, rather than a file the computer cannot open again.
        picked = next((suffix for name, suffix in EXPORT_FORMATS
                       if chosen and chosen.startswith(name)), "png")
        target = Path(path)
        target = target if target.suffix else target.with_suffix("." + picked)
        try:
            return self.export(target)
        except Exception as error:  # a full disk or a read-only folder, said plainly
            self.say(f"Could not save: {error}")
            return None

    def export(self, path, *, suffix=None):
        """Save the picture or the graph, at the standard lettering size.

        Text size is a per-reader comfort setting, so exports use the standard
        size whatever this reader chose: the same comparison then reads the same
        on every computer. The title and the line under it carry this window's
        typing kind, reference, target count, cohort and threshold into the file.
        """
        path = Path(path)
        suffix = str(suffix or path.suffix.lstrip(".")).lower()
        if callable(self.check_output):
            self.check_output(path)
        if suffix in {"png", "jpg", "jpeg"}:
            self.view.save_image(path, title=self.identity.export_title(),
                                 subtitle=self.identity.export_subtitle())
        elif suffix == "svg":
            self.view.save_svg(path, title=self.identity.export_title(),
                               subtitle=self.identity.export_subtitle())
        elif suffix == "graphml":
            self.view.save_graphml(path)
        elif suffix == "nwk":
            self.view.save_newick(path)
        else:
            raise ValueError(f"Unsupported graph export format: {suffix or 'none given'}")
        self.say(f"Saved {path.name} · {self.identity.caption()} · {self.identity.threshold_words()}.")
        return path

    # --- window --------------------------------------------------------------
    def say(self, message):
        self.statusBar().showMessage(str(message))

    def available_screen_size(self):
        screen = self.screen()
        return screen.availableGeometry().size() if screen is not None else None

    def remembered_size(self):
        """The size this kind of graph window was last left at, inside this screen."""
        stored = None
        if self.preferences is not None:
            stored = _parse_size(self.preferences.value(_size_key(self.identity.kind), ""))
        if stored is None:
            stored = _SESSION_SIZES.get(self.identity.kind)
        return clamp_size(*(stored or DEFAULT_SIZE), self.available_screen_size())

    def store_size(self):
        """Remember the size of a normal window; a maximised one keeps its own size."""
        geometry = self.normalGeometry() if self.isMaximized() or self.isFullScreen() else self.geometry()
        width, height = (geometry.width(), geometry.height()) if geometry.isValid() else (self.width(), self.height())
        width, height = clamp_size(width, height)
        _SESSION_SIZES[self.identity.kind] = (width, height)
        if self.preferences is not None:
            self.preferences.setValue(_size_key(self.identity.kind), f"{width}x{height}")
            self.preferences.sync()
        return width, height

    def closeEvent(self, event):
        self.store_size()
        # An arrangement made in the last moment before closing is still an
        # arrangement the host may want to keep.
        if self._state_timer.isActive():
            self._state_timer.stop()
            self._emit_state()
        super().closeEvent(event)
        if event.isAccepted():
            self.closed.emit()
