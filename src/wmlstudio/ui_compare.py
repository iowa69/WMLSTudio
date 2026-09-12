"""Explicit cohort/scheme selection, reusable profile views and graph presentation."""

import html
import json
import threading
from collections import Counter
from pathlib import Path
from types import SimpleNamespace

from PySide6.QtCore import QAbstractTableModel, QModelIndex, Qt, QThread, QTimer, Signal
from PySide6.QtGui import QColor
from PySide6.QtWidgets import (
    QCheckBox,
    QColorDialog,
    QComboBox,
    QDoubleSpinBox,
    QFileDialog,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMenu,
    QMessageBox,
    QSizePolicy,
    QSpinBox,
    QSplitter,
    QTableView,
    QTabWidget,
    QTextBrowser,
    QToolButton,
    QVBoxLayout,
    QWidget,
)

from wmlstudio.comparison import forest_from_distances, pairwise_distances
from wmlstudio.identification import scheme_organism
from wmlstudio.jobs import AnalysisWorker
from wmlstudio.sequence import AnalysisCancelled, check_cancelled
from wmlstudio.ui_common import FlowLayout, cell, make_table, organism_for
from wmlstudio.widgets import TreeView, button, label


def _available_profiles(project, sample):
    results = project.analysis_results(sample['id']) if hasattr(project, 'analysis_results') else []
    primary = sample.get('result') or {}
    if primary.get('alleles') and sample['status'] == 'completed':
        results = [primary, *[r for r in results if r.get('scheme_digest') != primary.get('scheme_digest')]]
    if sample['status'] in {'failed', 'interrupted', 'running'}:
        results = [r for r in results if r.get('scheme_digest') != primary.get('scheme_digest')]
    current_hash = primary.get('input_sha256') or sample.get('metadata', {}).get('assembly', {}).get('provenance', {}).get('assembly_sha256')
    if current_hash:
        results = [r for r in results if not r.get('input_sha256') or r['input_sha256'] == current_hash]
    return [r for r in results if r.get('alleles') and r.get('scheme_digest')]


def _requested_results(project, request, cancelled=None):
    """Read profile evidence without touching Qt widgets; safe in worker threads."""
    samples = project.samples()
    chosen = request['chosen'] if request['chosen'] is not None else {s['id'] for s in samples}
    available = []
    for sample in samples:
        check_cancelled(cancelled)
        genus, species, _ = organism_for(sample)
        if (sample['id'] not in chosen or (request['genus'] and genus != request['genus'])
                or (request['species'] and species != request['species'])):
            continue
        for result in _available_profiles(project, sample):
            available.append(dict(result, sample_id=sample['id'], sample_name=sample['name'],
                                  metadata=sample.get('metadata', {})))
    choice = request['choice']
    digest = choice.removeprefix('digest:') if choice and choice.startswith('digest:') else None
    if choice and not choice.startswith('digest:'):
        matching = [r for r in available if Path(str(r.get('scheme_path', ''))).resolve() == Path(choice).resolve()]
        if not matching:
            return [], len(chosen)
        digest = matching[0]['scheme_digest']
    if digest is None and available:
        digest = Counter(r['scheme_digest'] for r in available).most_common(1)[0][0]
    seen, results = set(), []
    for result in available:
        if result['scheme_digest'] == digest and result['sample_id'] not in seen:
            seen.add(result['sample_id'])
            results.append(result)
    return results, len(chosen)


def _calculate_comparison(project, request, cancelled=None):
    revision = project.comparison_revision()
    if request.get('revision', revision) != revision:
        raise ValueError('Stored profiles changed before comparison began. Build the comparison again.')
    results, total = _requested_results(project, request, cancelled)
    if len(results) > 500:
        return {'results': [], 'rows': [], 'edges': [], 'total': total, 'too_large': True}
    rows = pairwise_distances(results, request['overlap'], cancelled=cancelled)
    edges = forest_from_distances(rows, cancelled=cancelled)
    check_cancelled(cancelled)
    if project.comparison_revision() != revision:
        raise ValueError('Stored profiles changed during comparison. No mixed-version comparison was displayed; build again.')
    return {'results': results, 'rows': rows, 'edges': edges, 'total': total, 'too_large': False}


class ComparisonWorker(QThread):
    succeeded = Signal(int, object)
    failed = Signal(int, str)

    def __init__(self, project, request, generation, parent=None):
        super().__init__(parent)
        self.project, self.request, self.generation = project, request, generation
        self.cancel_event = threading.Event()

    def cancel(self):
        self.cancel_event.set()

    def run(self):
        try:
            payload = _calculate_comparison(self.project, self.request, self.cancel_event.is_set)
            self.succeeded.emit(self.generation, payload)
        except AnalysisCancelled:
            pass
        except Exception as error:
            self.failed.emit(self.generation, str(error))


class ComparisonTreeView(TreeView):
    """Fit all labels after a real viewport resize, without an expanding page hint."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self._resize_fit_timer = QTimer(self)
        self._resize_fit_timer.setSingleShot(True)
        self._resize_fit_timer.timeout.connect(self.fit_tree)

    def resizeEvent(self, event):
        super().resizeEvent(event)
        if hasattr(self, '_resize_fit_timer'):
            self._resize_fit_timer.start(80)


class EvidenceMatrixModel(QAbstractTableModel):
    """Read-only large matrices: cells are rendered on demand, not preallocated."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.headers = []
        self.rows = []

    def replace(self, headers, rows):
        self.beginResetModel()
        self.headers, self.rows = list(headers), list(rows)
        self.endResetModel()

    def rowCount(self, parent=QModelIndex()):
        return 0 if parent.isValid() else len(self.rows)

    def columnCount(self, parent=QModelIndex()):
        return 0 if parent.isValid() else len(self.headers)

    def data(self, index, role=Qt.ItemDataRole.DisplayRole):
        if index.isValid() and role in {Qt.ItemDataRole.DisplayRole, Qt.ItemDataRole.ToolTipRole}:
            value = self.rows[index.row()].get(self.headers[index.column()])
            return "—" if value is None else str(value)
        return None

    def headerData(self, section, orientation, role=Qt.ItemDataRole.DisplayRole):
        if role == Qt.ItemDataRole.DisplayRole:
            return self.headers[section] if orientation == Qt.Orientation.Horizontal else section + 1
        return None

    def sort(self, column, order=Qt.SortOrder.AscendingOrder):
        if not 0 <= column < len(self.headers):
            return
        self.layoutAboutToBeChanged.emit()
        key = self.headers[column]
        self.rows.sort(key=lambda row: str(row.get(key) or "").casefold(),
                       reverse=order == Qt.SortOrder.DescendingOrder)
        self.layoutChanged.emit()


class ComparisonWorkspaceMixin:
    def build_compare(self):
        _, layout = self.page()
        layout.setSpacing(8)
        self.heading(layout, "Compare a cohort", "Choose isolates and one reference snapshot. Reuse saved profiles without reading the genomes again.")
        filters = QHBoxLayout()
        self.compare_genus = QComboBox()
        self.compare_genus.setProperty('compactCharacters', 9)
        self.compare_genus.addItem("All genera", "")
        self.compare_species = QComboBox()
        self.compare_species.setProperty('compactCharacters', 9)
        self.compare_species.addItem("All species", "")
        self.compare_scheme = QComboBox()
        self.compare_scheme.addItem("Existing profile snapshot", None)
        self.compare_scheme.setMinimumWidth(190)
        self.compare_scheme.setSizeAdjustPolicy(QComboBox.SizeAdjustPolicy.AdjustToMinimumContentsLengthWithIcon)
        self.compare_scheme.setMinimumContentsLength(18)
        self.compare_genus.currentIndexChanged.connect(self.refresh_cohort_table)
        self.compare_species.currentIndexChanged.connect(self.refresh_cohort_table)
        self.compare_scheme.currentIndexChanged.connect(self.comparison_scheme_changed)
        filters.addWidget(self.compare_genus, 1)
        filters.addWidget(self.compare_species, 1)
        filters.addWidget(self.compare_scheme, 2)
        filters.addWidget(button("Call scheme…", self.type_comparison_scheme))
        layout.addLayout(filters)
        controls = FlowLayout()
        controls.addWidget(button("Use selection", self.compare_selected))
        self.overlap = QDoubleSpinBox()
        self.overlap.setRange(0.01, 1.0)
        self.overlap.setDecimals(2)
        self.overlap.setSingleStep(0.05)
        self.overlap.setValue(0.95)
        self.overlap.setToolTip("Minimum shared called loci divided by all loci. Missing data never becomes a matching allele.")
        self.cluster_threshold = QSpinBox()
        self.cluster_threshold.setRange(0, 100000)
        self.cluster_threshold.setValue(1)
        self.cluster_threshold.valueChanged.connect(self.refresh_comparison)
        for text, control in [('Shared loci ≥', self.overlap), ('Group threshold', self.cluster_threshold)]:
            group = QWidget()
            row = QHBoxLayout(group)
            row.setContentsMargins(0, 0, 0, 0)
            row.addWidget(label(text, 'small'))
            row.addWidget(control)
            controls.addWidget(group)
        controls.addWidget(button("Build comparison", self.refresh_comparison, True))
        layout.addLayout(controls)
        splitter = QSplitter(Qt.Orientation.Horizontal)
        splitter.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Ignored)
        splitter.setMinimumHeight(220)
        splitter.setChildrenCollapsible(True)
        chooser = QWidget()
        left = QVBoxLayout(chooser)
        left.setContentsMargins(0, 0, 8, 0)
        self.cohort_count = label("Choose a cohort", "cardTitle")
        left.addWidget(self.cohort_count)
        self.cohort_search = QLineEdit()
        self.cohort_search.setPlaceholderText("Find isolates…")
        self.cohort_search.textChanged.connect(self.refresh_cohort_table)
        left.addWidget(self.cohort_search)
        self.show_all_comparison_schemes = QCheckBox('Show all schemes')
        self.show_all_comparison_schemes.setToolTip('Include schemes for other organisms. Schemes with unknown organism metadata remain available in either mode.')
        self.show_all_comparison_schemes.toggled.connect(self.refresh_cohort_table)
        left.addWidget(self.show_all_comparison_schemes)
        cohort_buttons = FlowLayout()
        cohort_buttons.addWidget(button('All visible', lambda: self.set_cohort_visible(True)))
        cohort_buttons.addWidget(button('Clear', lambda: self.set_cohort_visible(False)))
        left.addLayout(cohort_buttons)
        self.cohort_table = make_table(["Use", "Sample", "Organism", "Stored profile"])
        self.cohort_table.setColumnWidth(0, 48)
        self.cohort_table.itemChanged.connect(self.cohort_item_changed)
        left.addWidget(self.cohort_table, 1)
        left.addWidget(button("Report this cohort", self.report_comparison))
        splitter.addWidget(chooser)
        right = QWidget()
        content = QVBoxLayout(right)
        content.setContentsMargins(0, 0, 0, 0)
        content.setSpacing(6)
        toolbar = FlowLayout()
        self.color_by = QComboBox()
        self.color_by.setProperty('compactCharacters', 12)
        self.color_by.addItem("Colour: cluster", "cluster")
        self.color_by.addItem("Colour: ST", "st")
        self.color_by.currentIndexChanged.connect(self.set_graph_color_by)
        toolbar.addWidget(self.color_by)
        graph_options = QToolButton()
        graph_options.setText('Graph options')
        graph_options.setPopupMode(QToolButton.ToolButtonPopupMode.InstantPopup)
        menu = QMenu(graph_options)
        menu.addAction('Fit graph', lambda: self.tree.fit_tree())
        focus = menu.addAction('Focus graph — hide cohort list')
        focus.setCheckable(True)
        focus.toggled.connect(lambda checked: chooser.setVisible(not checked))
        menu.addAction('Colour selected nodes…', self.color_graph_selection)
        menu.addAction('Highlight selected nodes…', self.highlight_graph_selection)
        menu.addSeparator()
        graph_options.setMenu(menu)
        toolbar.addWidget(graph_options)
        export = QComboBox()
        export.setProperty('compactCharacters', 12)
        export.addItems(["Export graph…", "PNG image", "SVG vector", "GraphML", "Newick (MST topology)", "Distance JSON", "Distance matrix TSV"])
        export.activated.connect(lambda index: self.export_graph_action(index, export))
        toolbar.addWidget(export)
        content.addLayout(toolbar)
        for text, method, default in [("Labels", "set_labels_visible", True), ("ST", "set_show_st", True),
                                      ("Edge distances", "set_edge_labels_visible", True),
                                      ("Merge identical", "set_merge_identical", False), ("Cluster halos", "set_halos_visible", False)]:
            check = menu.addAction(text)
            check.setCheckable(True)
            check.setChecked(default)
            check.toggled.connect(lambda enabled, name=method: self.graph_option(name, enabled))
        tabs = QTabWidget()
        tabs.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        self.tree = ComparisonTreeView()
        tabs.addTab(self.tree, "Graph")
        self.profile_table = QTableView()
        self.profile_model = EvidenceMatrixModel(self.profile_table)
        self.profile_table.setModel(self.profile_model)
        self.profile_table.setSortingEnabled(True)
        self.profile_table.setAlternatingRowColors(True)
        tabs.addTab(self.profile_table, "Profiles")
        self.statistics_view = QTextBrowser()
        tabs.addTab(self.statistics_view, "Quality / statistics")
        content.addWidget(tabs, 1)
        self.graph_legend = QLabel()
        self.graph_legend.setWordWrap(True)
        self.graph_legend.setTextFormat(Qt.TextFormat.RichText)
        content.addWidget(self.graph_legend)
        splitter.addWidget(right)
        splitter.setSizes([260, 800])
        layout.addWidget(splitter, 1)
        self.tree_status = label("Choose a cohort and a compatible scheme snapshot. Groups indicate allele similarity, not proven transmission.", "small", True)
        layout.addWidget(self.tree_status)
        self.distance_rows = []
        self._cohort_filling = False
        self._last_comparison = []
        self.comparison_worker = None
        self._comparison_generation = 0
        self._comparison_pending = None
        self._comparison_cache = None
        self._comparison_closing = False
        self._scheme_label_cache = {}
        self._comparison_timer = QTimer(self)
        self._comparison_timer.setSingleShot(True)
        self._comparison_timer.timeout.connect(self._start_pending_comparison)
        for name, callback in [("colorsChanged", self.persist_graph_style), ("layoutChanged", self.persist_graph_style),
                               ("legendChanged", self.update_graph_legend), ("nodeActivated", self.inspect_graph_sample)]:
            if hasattr(self.tree, name):
                getattr(self.tree, name).connect(callback)

    def available_profiles(self, sample):
        return _available_profiles(self.project, sample)

    def invalidate_comparison_scheme_labels(self, path=None):
        """Explicit reference-revision hook; sample/progress refreshes must not clear it.

        Bundled references and installed snapshots are immutable during normal
        use. Import, update and the manual installed-reference refresh all call
        ``populate_schemes``; callers replacing metadata directly must invoke
        this hook (optionally for just the replaced path) before refreshing.
        """
        if not hasattr(self, '_scheme_label_cache') or path is None:
            self._scheme_label_cache = {}
        else:
            self._scheme_label_cache.pop(str(Path(path)), None)

    def populate_schemes(self):
        self.invalidate_comparison_scheme_labels()
        return super().populate_schemes()

    def comparison_scheme_label(self, path):
        """Read bounded metadata once per path/reference revision, never alleles."""
        path = Path(path)
        # Do not resolve, stat or glob before this lookup: even metadata-only
        # directory scans stall repeated progress refreshes on Windows.
        key = str(path)
        if key in self._scheme_label_cache:
            return self._scheme_label_cache[key]
        try:
            candidates = [path / 'scheme.json'] if (path / 'scheme.json').is_file() else sorted(path.glob('*_info.json'))
            candidate = candidates[0] if len(candidates) == 1 else None
            stat = candidate.stat() if candidate else None
            metadata = json.loads(candidate.read_text(encoding='utf-8')) if candidate and stat.st_size <= 2 * 1024 * 1024 else {}
            if not isinstance(metadata, dict):
                metadata = {}
        except (OSError, ValueError):
            metadata = {}
        name = str(metadata.get('name') or path.name.replace('_', ' '))
        _, organism = scheme_organism(SimpleNamespace(metadata=metadata, name=name))
        value = (name, organism.get('genus', ''), organism.get('species', ''))
        self._scheme_label_cache[key] = value
        return value

    def refresh_cohort_table(self):
        if not hasattr(self, "cohort_table") or self._cohort_filling:
            return
        self._cohort_filling = True
        try:
            samples = self.project.samples()
            taxa = [organism_for(s) for s in samples]
            self.fill_filter(self.compare_genus, [t[0] for t in taxa], "All genera")
            selected_genus = self.compare_genus.currentData()
            self.fill_filter(self.compare_species, [t[1] for t in taxa if not selected_genus or t[0] == selected_genus], "All species")
            selected_scheme = self.compare_scheme.currentData()
            self.compare_scheme.blockSignals(True)
            self.compare_scheme.clear()
            self.compare_scheme.addItem("Use a stored compatible snapshot", None)
            fingerprints = {}
            for sample in samples:
                genus, species, _ = organism_for(sample)
                if (selected_genus and genus and genus != selected_genus) or (
                        self.compare_species.currentData() and species and species != self.compare_species.currentData()):
                    continue
                for result in self.available_profiles(sample):
                    fingerprints[result["scheme_digest"]] = result.get("scheme") or "Unnamed scheme"
            for digest, name in sorted(fingerprints.items(), key=lambda p: (p[1], p[0])):
                self.compare_scheme.addItem(f"Saved: {name} · {digest[:8]}", "digest:" + digest)
            for path in self.scheme_paths:
                name, genus, species = self.comparison_scheme_label(path)
                if not self.show_all_comparison_schemes.isChecked() and (
                    (selected_genus and genus and genus != selected_genus) or
                    (self.compare_species.currentData() and species and species != self.compare_species.currentData())
                ):
                    continue
                self.compare_scheme.addItem('Call: ' + name + (' · organism unknown' if not genus else ''), str(path))
            self.compare_scheme.setCurrentIndex(max(0, self.compare_scheme.findData(selected_scheme)))
            self.compare_scheme.blockSignals(False)
            chosen = self.cohort_ids if self.cohort_ids is not None else {s["id"] for s in samples}
            query = self.cohort_search.text().strip().casefold()
            visible = []
            for sample in samples:
                genus, species, _ = organism_for(sample)
                if self.compare_genus.currentData() and genus != self.compare_genus.currentData():
                    continue
                if self.compare_species.currentData() and species != self.compare_species.currentData():
                    continue
                if query and query not in sample["name"].casefold():
                    continue
                visible.append(sample)
            self.cohort_table.blockSignals(True)
            self.cohort_table.setSortingEnabled(False)
            self.cohort_table.setRowCount(len(visible))
            for row, sample in enumerate(visible):
                genus, species, _ = organism_for(sample)
                values = ["", sample["name"], " ".join([genus, species]).strip() or "Unknown",
                          "; ".join(sorted({r.get("scheme", "") for r in self.available_profiles(sample)})) or "Not called"]
                for col, value in enumerate(values):
                    item = cell(value, sample["id"])
                    if col == 0:
                        item.setFlags(item.flags() | Qt.ItemFlag.ItemIsUserCheckable)
                        item.setCheckState(Qt.CheckState.Checked if sample["id"] in chosen else Qt.CheckState.Unchecked)
                    self.cohort_table.setItem(row, col, item)
            self.cohort_table.setSortingEnabled(True)
            self.cohort_table.blockSignals(False)
            self.cohort_count.setText(f"{len(chosen)} in cohort · {len(visible)} visible")
        finally:
            self._cohort_filling = False

    def cohort_item_changed(self, item):
        if self._cohort_filling or item.column() != 0:
            return
        if self.cohort_ids is None:
            self.cohort_ids = {s["id"] for s in self.project.samples()}
        sid = item.data(Qt.ItemDataRole.UserRole)
        if item.checkState() == Qt.CheckState.Checked:
            self.cohort_ids.add(sid)
        else:
            self.cohort_ids.discard(sid)
        self.project.set_setting("comparison_cohort", sorted(self.cohort_ids))
        self.cohort_count.setText(f"{len(self.cohort_ids)} in cohort · rebuild to apply")

    def set_cohort_visible(self, selected):
        if self.cohort_ids is None:
            self.cohort_ids = {s["id"] for s in self.project.samples()}
        visible = {self.cohort_table.item(row, 0).data(Qt.ItemDataRole.UserRole)
                   for row in range(self.cohort_table.rowCount())}
        if selected:
            self.cohort_ids.update(visible)
        else:
            self.cohort_ids.difference_update(visible)
        self.project.set_setting("comparison_cohort", sorted(self.cohort_ids))
        self.refresh_cohort_table()

    def compare_selected(self):
        if not self.selection_ids:
            self.notify("Select samples in the library first, or use the cohort checkboxes in Compare.")
            self.navigate(2)
            return
        self.cohort_ids = set(self.selection_ids)
        self.project.set_setting("comparison_cohort", sorted(self.cohort_ids))
        self.refresh_cohort_table()
        self.navigate(2)

    def comparison_scheme_changed(self):
        if not self._cohort_filling:
            self.tree_status.setText("Scheme selection changed. Build from saved evidence, or choose Call this scheme for missing profiles.")

    def comparison_results(self):
        return _requested_results(self.project, self._comparison_request())[0]

    def _comparison_request(self):
        return {'chosen': frozenset(self.cohort_ids) if self.cohort_ids is not None else None,
                'choice': self.compare_scheme.currentData(), 'genus': self.compare_genus.currentData(),
                'species': self.compare_species.currentData(), 'overlap': self.overlap.value()}

    def _comparison_key(self, request, revision):
        return (str(self.project.path), revision, request['chosen'], request['choice'],
                request['genus'], request['species'], request['overlap'])

    def _clear_comparison(self):
        self.distance_rows, self._last_comparison = [], []
        if self._pending_graph_state is None and getattr(self.tree, 'nodes', None):
            self._pending_graph_state = self.tree.export_state()
        self.tree.blockSignals(True)
        self.tree.draw_results([], [], self.cluster_threshold.value())
        self.tree.blockSignals(False)
        self.graph_legend.clear()
        self.refresh_profile_matrix([])
        self.refresh_statistics([], [])

    def refresh_comparison(self):
        if not hasattr(self, "tree") or self._comparison_closing:
            return
        self._comparison_generation += 1
        self._comparison_pending = None
        self._comparison_timer.stop()
        if self.comparison_worker and self.comparison_worker.isRunning():
            self.comparison_worker.cancel()
        request = self._comparison_request()
        revision = self.project.comparison_revision()
        request['revision'] = revision
        key = self._comparison_key(request, revision)
        if self._comparison_cache and self._comparison_cache[0] == key:
            self._apply_comparison(self._comparison_cache[1])
            return
        selected_count = len(request['chosen']) if request['chosen'] is not None else len(revision)
        if selected_count > 30:
            self._clear_comparison()
            self.tree_status.setText('Calculating cohort distances in the background… You can continue navigating. New selections replace pending work.')
            self._comparison_pending = (self._comparison_generation, request, key)
            self._comparison_timer.start(150)
            return
        try:
            payload = _calculate_comparison(self.project, request)
            self._comparison_cache = (key, payload)
            self._apply_comparison(payload)
        except Exception as exc:
            self._clear_comparison()
            self.tree_status.setText(str(exc))

    def _start_pending_comparison(self):
        if self._comparison_closing or self._comparison_pending is None:
            return
        if self.comparison_worker and self.comparison_worker.isRunning():
            self.comparison_worker.cancel()
            return
        generation, request, key = self._comparison_pending
        self._comparison_pending = None
        worker = ComparisonWorker(self.project, request, generation, self)
        self.comparison_worker = worker
        worker.succeeded.connect(lambda received, payload: self._comparison_ready(received, key, payload))
        worker.failed.connect(self._comparison_failed)
        worker.finished.connect(lambda: self._comparison_finished(worker))
        worker.start()

    def _comparison_ready(self, generation, key, payload):
        if (self._comparison_closing or generation != self._comparison_generation
                or key[0] != str(self.project.path)):
            return
        self._comparison_cache = (key, payload)
        self._apply_comparison(payload)

    def _comparison_failed(self, generation, message):
        if generation == self._comparison_generation and not self._comparison_closing:
            self.tree_status.setText('Comparison could not complete: ' + message)

    def _comparison_finished(self, worker):
        if self.comparison_worker is worker:
            self.comparison_worker = None
        worker.deleteLater()
        if self._comparison_closing:
            QTimer.singleShot(0, self.close)
        elif self._comparison_pending is not None:
            self._comparison_timer.start(0)

    def cancel_comparison_for_close(self):
        """Return True once no comparison work can read the closing project."""
        self._comparison_closing = True
        self._comparison_generation += 1
        self._comparison_timer.stop()
        self._comparison_pending = None
        if self.comparison_worker and self.comparison_worker.isRunning():
            self.comparison_worker.cancel()
            return False
        return True

    def _apply_comparison(self, payload):
        results, edges = payload['results'], payload['edges']
        try:
            self.distance_rows = payload['rows']
            self.tree.draw_results(results, edges, self.cluster_threshold.value())
            if self._pending_graph_state is not None and hasattr(self.tree, "restore_state"):
                self.tree.blockSignals(True)
                self.tree.restore_state(self._pending_graph_state)
                self.tree.blockSignals(False)
                self._pending_graph_state = None
                if hasattr(self.tree, 'legend'):
                    self.update_graph_legend(self.tree.legend())
            self._last_comparison = results
            self.refresh_profile_matrix(results)
            self.refresh_statistics(results, edges)
            excluded = sum(not row["comparable"] for row in self.distance_rows)
            total = payload['total']
            name = results[0].get("scheme", "") if results else "No compatible profiles"
            self.tree_status.setText(f"{name} · {len(results)} / {total} cohort members have this profile · {len(edges)} edges · {excluded} pairs excluded. Missing/incompatible profiles are not silently compared.")
            if payload['too_large']:
                self.tree_status.setText('Select at most 500 profiles for this interactive view. Export larger cohorts for offline comparison.')
            if hasattr(self.tree, "available_color_fields"):
                current = self.color_by.currentData()
                self.color_by.blockSignals(True)
                self.color_by.clear()
                self.color_by.addItem("Colour: cluster", "cluster")
                self.color_by.addItem("Colour: ST", "st")
                for field in self.tree.available_color_fields():
                    self.color_by.addItem(str(field), "metadata:" + str(field).removeprefix("metadata:"))
                self.color_by.setCurrentIndex(max(0, self.color_by.findData(current)))
                self.color_by.blockSignals(False)
        except Exception as exc:
            self.tree_status.setText(str(exc))

    def type_comparison_scheme(self):
        if self.busy():
            return
        choice = self.compare_scheme.currentData()
        if not choice or choice.startswith("digest:"):
            self.notify("Choose a 'Call:' entry in the scheme selector. Saved snapshots can be compared immediately.")
            return
        chosen = self.cohort_ids if self.cohort_ids is not None else {s["id"] for s in self.project.samples()}
        samples = [s for s in self.project.samples() if s["id"] in chosen and s.get("input_path") and not s.get("missing_input")]
        if not samples:
            self.error("This cohort has no available sequence inputs. Import a profile table or reconnect the assemblies.")
            return
        if QMessageBox.question(self, "Call comparison scheme", f"Call {Path(choice).name} for {len(samples)} input samples?\n\nPrimary MLST/ST results are retained. New profiles are stored separately and reused on subsequent comparisons. FASTQ files cannot be typed without assembly.") != QMessageBox.StandardButton.Yes:
            return
        samples = [dict(s, metadata={**s.get("metadata", {}), "workflow": {"typing_mode": "manual", "scheme_path": choice}}) for s in samples]
        self.worker_role = "secondary"
        self.worker = AnalysisWorker(samples, choice, parent=self)

        def saved(sid, result):
            if result.get("alleles"):
                result["scheme_path"] = choice
                self.project.set_analysis(sid, result)
                self.project.set_setting("last_comparison_digest", result["scheme_digest"])
            else:
                self.notify("A selected read file received QC only; no allelic profile was invented.")

        self.worker.sample_finished.connect(saved)
        self.worker.sample_failed.connect(lambda sid, message: self.notify(f"Comparison calling failed: {message}"))
        self.worker.progress.connect(self.job_progress)
        self.worker.finished.connect(self.secondary_finished)
        self.set_running(True)
        self.worker.start()

    def secondary_finished(self):
        self.analysis_finished()
        if self.closing_after_cancel:
            return
        digest = self.project.get_setting("last_comparison_digest", "")
        index = self.compare_scheme.findData("digest:" + digest)
        if index >= 0:
            self.compare_scheme.setCurrentIndex(index)
        self.refresh_comparison()

    def refresh_profile_matrix(self, results):
        loci = sorted({locus for result in results for locus in result.get("alleles", {})})
        self.profile_model.replace(["Sample", "ST", "Scheme", *loci],
                                  [{"Sample": r["sample_name"], "ST": r.get("st"), "Scheme": r.get("scheme"), **r.get("alleles", {})} for r in results])

    def refresh_statistics(self, results, edges):
        comparable = [p for p in self.distance_rows if p["comparable"]]
        distances = sorted(p["distance"] for p in comparable)
        total_loci = len(results[0].get("alleles", {})) if results else 0
        parts = ["<h2>Cohort evidence</h2>", f"<p>{len(results)} profiles · {total_loci} loci · {len(comparable)} comparable pairs · {len(self.distance_rows) - len(comparable)} excluded pairs</p>"]
        if distances:
            median = (distances[(len(distances) - 1) // 2] + distances[len(distances) // 2]) / 2
            parts.append(f"<p><b>Allele differences:</b> minimum {distances[0]}, median {median:g}, maximum {distances[-1]}.</p>")
            distribution = Counter(distances)
            parts.append("<h3>Distance distribution</h3><table cellpadding='5'><tr><th>Differences</th><th>Pairs</th></tr>" + "".join(f"<tr><td>{distance}</td><td>{count}</td></tr>" for distance, count in sorted(distribution.items())) + "</table>")
        parts.append("<h3>Per-sample completeness</h3><table cellpadding='5'><tr><th>Sample</th><th>Called loci</th><th>Missing</th></tr>")
        for result in results:
            alleles = result.get("alleles", {})
            called = sum(value is not None for value in alleles.values())
            parts.append(f"<tr><td>{html.escape(result['sample_name'])}</td><td>{called}/{len(alleles)}</td><td>{len(alleles) - called}</td></tr>")
        parts.append("</table><p>Threshold groups are local views of allele similarity, not curated global cluster numbers or a transmission diagnosis. MST edge lengths are not a phylogenetic clock.</p>")
        self.statistics_view.setHtml("".join(parts))

    def graph_option(self, name, value):
        if hasattr(self.tree, name):
            getattr(self.tree, name)(value)
            self.persist_graph_style()

    def set_graph_color_by(self):
        if hasattr(self.tree, "set_color_by"):
            self.tree.set_color_by(self.color_by.currentData() or "cluster")
            self.persist_graph_style()

    def color_graph_selection(self):
        color = QColorDialog.getColor(parent=self, title="Colour selected graph nodes")
        if color.isValid() and hasattr(self.tree, "set_selected_color"):
            self.tree.set_selected_color(color.name())
            self.persist_graph_style()

    def highlight_graph_selection(self):
        if hasattr(self.tree, "selected_ids"):
            ids = set(self.tree.selected_ids())
            if ids:
                self.selection_ids = ids
        self.highlight_selected()

    def persist_graph_style(self, *_args):
        if self._pending_graph_state is None and hasattr(self.tree, "export_state"):
            self.project.set_setting("graph_style", self.tree.export_state())

    def update_graph_legend(self, values):
        entries = []
        for name, value in values.items():
            color = QColor(str(value))
            swatch = color.name() if color.isValid() else '#7b8496'
            entries.append(f'<span style="color:{swatch}">●</span> {html.escape(str(name))}')
        self.graph_legend.setText(' &nbsp; '.join(entries))

    def inspect_graph_sample(self, sid):
        self.selection_ids = {sid}
        self.library_filter = ("ids", [sid])
        self.refresh_tables()
        for row in range(self.sample_table.rowCount()):
            if self.sample_table.item(row, 0).data(Qt.ItemDataRole.UserRole) == sid:
                self.sample_table.setCurrentCell(row, 0)
        self.navigate(1)

    def export_graph_action(self, index, combo):
        combo.setCurrentIndex(0)
        if index == 0:
            return
        if index == 1:
            return self.save_tree()
        if index == 5:
            return self.export_distances()
        formats = {2: ("svg", "save_svg"), 3: ("graphml", "save_graphml"), 4: ("nwk", "save_newick"), 6: ("tsv", None)}
        extension, method = formats[index]
        path, _ = QFileDialog.getSaveFileName(self, "Export comparison", f"comparison.{extension}", f"{extension.upper()} (*.{extension})")
        if not path:
            return
        try:
            self.check_output(path)
            if method:
                getattr(self.tree, method)(path)
            else:
                self.write_distance_matrix(path)
            self.notify("Comparison exported with its current cohort and evidence.")
        except Exception as exc:
            self.error(exc)

    def write_distance_matrix(self, path):
        import csv

        from wmlstudio.export import _atomic_text, _csv_cell
        results = self._last_comparison
        distances = {frozenset((p["source"], p["target"])): p["distance"] if p["comparable"] else "NA" for p in self.distance_rows}
        with _atomic_text(path) as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["sample", *[_csv_cell(r["sample_name"]) for r in results]])
            for result in results:
                sid = result["sample_id"]
                writer.writerow([_csv_cell(result["sample_name"]), *[0 if other["sample_id"] == sid
                                else distances.get(frozenset((sid, other["sample_id"])), "NA") for other in results]])

    def report_comparison(self):
        self.report_ids = {r["sample_id"] for r in self._last_comparison}
        self.report_scope.setCurrentIndex(0)
        self.refresh_report_table()
        self.navigate(5)
