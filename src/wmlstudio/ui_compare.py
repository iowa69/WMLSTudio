"""Explicit cohort/scheme selection, reusable profile views and graph presentation."""

import html
import json
import threading
from collections import Counter
from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace

from PySide6.QtCore import QAbstractTableModel, QModelIndex, Qt, QThread, QTimer, Signal
from PySide6.QtGui import QColor
from PySide6.QtWidgets import (
    QApplication,
    QCheckBox,
    QColorDialog,
    QComboBox,
    QDialog,
    QDialogButtonBox,
    QDoubleSpinBox,
    QFileDialog,
    QFormLayout,
    QHBoxLayout,
    QInputDialog,
    QLabel,
    QLineEdit,
    QListWidget,
    QListWidgetItem,
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
from wmlstudio.context_menus import TableWidgetAdapter
from wmlstudio.identification import cached_scheme, scheme_organism
from wmlstudio.investigation import (
    DIFF_ROW_FIELDS,
    SCALE_SEPARATION,
    TYPING_SCALES,
    InvestigationStore,
    build_snapshot,
    incremental_distances,
    proximity_rows,
    scheme_key_from,
    scheme_threshold_suggestion,
    snapshot_diff,
    snapshot_diff_caption,
    snapshot_diff_html,
    snapshot_diff_rows,
    snapshot_graph,
    typing_scale,
)
from wmlstudio.jobs import AnalysisWorker
from wmlstudio.project import CGMLST_LOCUS_FLOOR
from wmlstudio.project import typing_kind as stored_typing_kind
from wmlstudio.sample_workflow import current_input_sha256, hydra_evidence_status
from wmlstudio.sequence import AnalysisCancelled, SequenceReader, check_cancelled, file_sha256
from wmlstudio.ui_common import FlowLayout, cell, gene_names, make_table, organism_for
from wmlstudio.widgets import TreeView, button, label, render_side_by_side

# Laying out a forest is measurably slow well before it stops being readable:
# roughly 0.4 s at 120 nodes and 6.6 s at 500 on a development machine. Above
# this many stored profiles the baseline is drawn only when the user asks, with
# the wait stated on the button rather than spent silently.
BASELINE_AUTODRAW_LIMIT = 120
# The two comparisons this page can show. They are two answers to two different
# questions over two different locus sets, so each keeps its own reference, its
# own target count, its own threshold, its own legend and its own exports, and
# no view ever mixes their numbers.
TYPING_VIEWS = ("mlst", "cgmlst")
COUNTERPART = {"mlst": "cgmlst", "cgmlst": "mlst"}


def typing_title(kind):
    return TYPING_SCALES.get(str(kind), TYPING_SCALES["unclassified"])["title"]


def _summary_kind(row):
    """The typing kind of one stored profile, from the kind column where it exists.

    The fallback repeats ``project.classify_typing`` over exactly the fields a
    lightweight summary row carries, so a project written before the kind column
    existed is classified the same way rather than dropped out of both views.
    """
    kind = str(row.get('typing_kind') or '').strip().casefold()
    if kind in TYPING_VIEWS:
        return kind
    declared = str(row.get('analysis_kind') or '').strip().casefold()
    if declared in TYPING_VIEWS:
        return declared
    count = int(row.get('locus_count') or 0)
    if not count or not row.get('scheme_digest'):
        return 'unclassified'
    return 'cgmlst' if count > CGMLST_LOCUS_FLOOR else 'mlst'


def _profile_kinds(project, sample_id):
    """{scheme_digest: typing kind} for one sample, read from the kind column."""
    if not hasattr(project, 'analysis_summaries'):
        return {}
    try:
        rows = project.analysis_summaries(sample_id)
    except KeyError:
        return {}
    return {row['scheme_digest']: _summary_kind(row) for row in rows if row.get('scheme_digest')}


def _available_profiles(project, sample, kind=None):
    """Current profiles for one isolate, optionally only those of one typing kind.

    ``kind`` filters on the stored typing kind, never on the scheme's name: a
    classical ST and a core-genome profile are separate evidence, and neither
    stands in for the other when the requested one is absent.
    """
    if sample.get('status') != 'completed':
        return []  # Archived evidence is retained, but failed/pending work is not current.
    results = project.analysis_results(sample['id']) if hasattr(project, 'analysis_results') else []
    primary = sample.get('result') or {}
    if primary.get('alleles') and sample['status'] == 'completed':
        results = [primary, *[r for r in results if r.get('scheme_digest') != primary.get('scheme_digest')]]
    if sample['status'] in {'failed', 'interrupted', 'running'}:
        results = [r for r in results if r.get('scheme_digest') != primary.get('scheme_digest')]
    current_hash = current_input_sha256(sample) or primary.get('input_sha256')
    if current_hash:
        results = [r for r in results if not r.get('input_sha256') or r['input_sha256'] == current_hash]
    results = [r for r in results if r.get('alleles') and r.get('scheme_digest')]
    if kind is None:
        return results
    recorded = _profile_kinds(project, sample['id'])
    return [r for r in results
            if recorded.get(r.get('scheme_digest'), stored_typing_kind(r)) == kind]


def _requested_results(project, request, cancelled=None, kind=None):
    """Read profile evidence without touching Qt widgets; safe in worker threads."""
    samples = project.samples()
    chosen = request['chosen'] if request['chosen'] is not None else {s['id'] for s in samples}
    kind = kind or request.get('typing_kind')
    available = []
    for sample in samples:
        check_cancelled(cancelled)
        genus, species, _ = organism_for(sample)
        if (sample['id'] not in chosen or (request['genus'] and genus != request['genus'])
                or (request['species'] and species != request['species'])):
            continue
        for result in _available_profiles(project, sample, kind):
            available.append(dict(result, sample_id=sample['id'], sample_name=sample['name'],
                                  metadata=sample.get('metadata', {}),
                                  primary_st=(sample.get('result') or {}).get('st'),
                                  primary_scheme=(sample.get('result') or {}).get('scheme'),
                                  amr_genes=gene_names(sample), amr_evidence_status=hydra_evidence_status(sample)['status'],
                                  organism=' '.join([genus, species]).strip()))
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
    """One comparison of one typing kind, plus the other kind's tree when asked for.

    The two kinds are computed separately and returned separately. Nothing is
    merged, averaged or carried across: the counterpart is a second answer to a
    second question about the same isolates, built without an investigation of
    its own so that showing it can never write a snapshot for the wrong scale.
    """
    revision = project.comparison_revision()
    if request.get('revision', revision) != revision:
        raise ValueError('Stored profiles changed before comparison began. Build the comparison again.')
    kind = request.get('typing_kind') or 'mlst'
    payload = _kind_comparison(project, request, kind, revision,
                               investigation_id=request.get('investigation_id'), cancelled=cancelled)
    counterpart = request.get('counterpart_kind')
    if counterpart and counterpart != kind:
        # The other tree is grouped by its own link threshold. Reusing this
        # tree's number would put a core-genome cohort under a seven-locus
        # cutoff, or the reverse, which is the mistake this page exists to avoid.
        other = dict(request, threshold=request.get('counterpart_threshold',
                                                    request.get('threshold', 1)))
        payload['counterpart'] = _kind_comparison(project, other, counterpart, revision,
                                                  investigation_id=None, cancelled=cancelled)
    return payload


def _kind_comparison(project, request, kind, revision, *, investigation_id=None, cancelled=None):
    results, total = _requested_results(project, request, cancelled, kind=kind)
    available_count = len(results)
    if len(results) > 500:
        return {'results': [], 'rows': [], 'edges': [], 'total': total, 'too_large': True,
                'typing_kind': kind, 'scale': typing_scale([], kind)}
    store = InvestigationStore(project)
    prior = store.snapshot(investigation_id) if investigation_id else None
    cache = store.cache(investigation_id) if investigation_id else None
    if investigation_id:
        plan = store.get(investigation_id)
        scheme = results[0].get('scheme') if results else plan.get('scheme')
        digest = results[0].get('scheme_digest') if results else plan.get('scheme_digest')
        present = {r['sample_id'] for r in results}
        chosen = request['chosen']
        for sample in project.samples():
            check_cancelled(cancelled)
            if sample['id'] in present or (chosen is not None and sample['id'] not in chosen):
                continue
            genus, species, _ = organism_for(sample)
            if (request['genus'] and genus != request['genus']) or (request['species'] and species != request['species']):
                continue
            if not scheme or not digest or sample.get('metadata', {}).get('workflow', {}).get('source_kind') == 'read_mate':
                continue
            primary = sample.get('result') or {}
            results.append({'sample_id': sample['id'], 'sample_name': sample['name'], 'scheme': scheme,
                            'scheme_digest': digest, 'alleles': {}, 'calls': [], 'status': 'not_profiled',
                            'comparison_unavailable': True, 'metadata': sample.get('metadata', {}),
                            'input_sha256': primary.get('input_sha256'), 'primary_st': primary.get('st'),
                            'primary_scheme': primary.get('scheme'), 'amr_genes': gene_names(sample),
                            'organism': ' '.join([genus, species]).strip()})
        if len(results) > 500:
            return {'results': [], 'rows': [], 'edges': [], 'total': total, 'too_large': True,
                    'typing_kind': kind, 'scale': typing_scale([], kind)}
    if investigation_id:
        rows, cache, reuse = incremental_distances(results, request['overlap'], cache, cancelled=cancelled)
    else:
        rows = pairwise_distances(results, request['overlap'], cancelled=cancelled)
        reuse = {'reused_pairs': 0, 'computed_pairs': len(rows)}
    edges = forest_from_distances(rows, cancelled=cancelled)
    snapshot = build_snapshot(results, rows, request.get('threshold', 1), request['overlap'],
                              previous=prior, reuse=reuse, kind=kind, cancelled=cancelled)
    check_cancelled(cancelled)
    if project.comparison_revision() != revision:
        raise ValueError('Stored profiles changed during comparison. No mixed-version comparison was displayed; build again.')
    if investigation_id:
        plan = store.get(investigation_id)
        if plan.get('updated_at') != request.get('investigation_revision'):
            raise ValueError('Investigation settings changed during comparison. Build again.')
        same_cohort = request['chosen'] is not None and set(request['chosen']) == set(plan['sample_ids'])
        check_cancelled(cancelled)
        if same_cohort and plan['threshold'] == request.get('threshold', 1) and plan['min_overlap'] == request['overlap']:
            snapshot = store.save_snapshot(investigation_id, snapshot, cache,
                expected_revision=request.get('investigation_revision'), expected_project_revision=revision, cancelled=cancelled)
        else:
            snapshot.update(preview=True, investigation_id=investigation_id, investigation_name=plan['name'],
                            protocol=plan['protocol'], saved_threshold=plan['threshold'],
                            saved_min_overlap=plan['min_overlap'], threshold_evidence=plan.get('threshold_evidence', {}))
    return {'results': results, 'rows': rows, 'edges': edges, 'total': total, 'too_large': False,
            'snapshot': snapshot, 'cache': cache, 'reuse': reuse, 'available_count': available_count,
            'typing_kind': kind,
            'scale': typing_scale(results, kind, scheme=snapshot.get('scheme'),
                                  scheme_digest=snapshot.get('scheme_digest'))}


class ComparisonCallingWorker(AnalysisWorker):
    """Skip valid same-input/same-reference evidence without retyping old isolates."""

    prepared = Signal(dict)

    def __init__(self, samples, scheme_path, saved_profiles, parent=None, resource_plan=None):
        super().__init__(samples, scheme_path, parent=parent, resource_plan=resource_plan)
        self.saved_profiles = saved_profiles

    def run(self):
        original = self.samples
        pending, reused, reads = [], [], []
        try:
            scheme = cached_scheme(self.scheme_path, self.cancel_event.is_set)
            for index, sample in enumerate(original):
                check_cancelled(self.cancel_event.is_set)
                self.progress.emit(int(index / max(1, len(original)) * 15), 'Checking stored input and reference fingerprints…')
                with SequenceReader(sample['input_path'], cancelled=self.cancel_event.is_set) as reader:
                    if reader.kind == 'fastq':
                        reads.append(sample['id'])
                        continue
                candidates = [p for p in self.saved_profiles.get(sample['id'], [])
                              if p.get('scheme_digest') == scheme.digest and p.get('input_sha256')]
                current = file_sha256(sample['input_path'], self.cancel_event.is_set) if candidates else None
                if (sample.get('status') == 'completed' and candidates
                        and any(p['input_sha256'] == current for p in candidates)):
                    reused.append(sample['id'])
                else:
                    pending.append(sample)
            self.samples = pending
            self.prepared.emit({'scheme_digest': scheme.digest, 'scheme': scheme.name,
                                'reused': reused, 'pending': [s['id'] for s in pending], 'reads': reads})
            super().run()
        except AnalysisCancelled:
            self.progress.emit(100, 'Comparison calling cancelled')
        except Exception as error:
            for sample in original:
                self.sample_failed.emit(sample['id'], str(error))


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


class _ClusterTableAdapter(TableWidgetAdapter):
    """Group rows carry a group id, but the actions act on isolates.

    Reading the row's stored identifier as a sample id would send group ids to
    handlers that archive or export isolates, so the members are resolved from
    the snapshot the table was filled from.
    """

    def __init__(self, view_id, table, window):
        super().__init__(view_id, table)
        self._window = window

    def group_members(self, group_ids):
        snapshot = getattr(self._window, '_current_snapshot', None) or {}
        return tuple(dict.fromkeys(str(sid) for group in snapshot.get('groups', [])
                                   if group['id'] in group_ids for sid in group['members']))

    def resolve(self, point):
        base = super().resolve(point)
        return replace(base, sample_ids=self.group_members(set(base.sample_ids)),
                       clicked_id=None, group_id=base.clicked_id)


class ComparisonWorkspaceMixin:
    def investigation_summary(self):
        """Small, JSON-safe status for the guided overview; never loads a matrix."""
        active = getattr(self, 'active_investigation_id', None)
        kind = getattr(self, 'typing_kind', 'mlst')
        if not active:
            return {'id': None, 'name': 'Unsaved investigation', 'saved': False,
                    'typing_kind': kind, 'typing_title': typing_title(kind),
                    'cohort_count': len(self.cohort_ids) if self.cohort_ids is not None else len(self.project.samples()),
                    'profiles': len(getattr(self, '_last_comparison', []))}
        try:
            plan = InvestigationStore(self.project).get(active)
        except KeyError:
            return {'id': None, 'name': 'Investigation unavailable', 'saved': False}
        return {'id': active, 'name': plan['name'], 'saved': True, 'cohort_count': len(plan['sample_ids']),
                'typing_kind': kind, 'typing_title': typing_title(kind),
                'target_loci': (self._current_snapshot or {}).get('target_loci'),
                'scheme': plan['scheme'], 'scheme_digest': plan['scheme_digest'],
                'threshold': plan['threshold'], 'min_overlap': plan['min_overlap'], 'protocol': plan['protocol'],
                'snapshots': len(plan['snapshots']), 'review_groups': len(plan['review_groups']),
                **plan.get('summary', {})}

    def investigation_actions(self):
        return {'new': self.create_investigation_dialog, 'update': self.update_investigation,
                'compare': self.refresh_comparison, 'report_selection': self.report_graph_selection,
                'report_cohort': self.report_comparison}

    def restore_investigations(self):
        if not hasattr(self, 'investigation_combo') or getattr(self, '_comparison_closing', False):
            return
        store = InvestigationStore(self.project)
        self._report_investigation_snapshot = None
        # Switching project brings a different baseline, arrangement and dual-view
        # preference; none of the old project's state may survive the swap.
        self._reset_baseline_state()
        self._legends = {'current': {}, 'baseline': {}}
        self._pending_baseline_graph_state = self.project.get_setting('graph_style.baseline', {'version': 1})
        self.dual_toggle.blockSignals(True)
        self.dual_toggle.setChecked(bool(self.project.get_setting('compare.dual_graph', False)))
        self.dual_toggle.blockSignals(False)
        stored_kind = self.project.get_setting('compare.typing_kind', None)
        self._typing_kind_chosen = stored_kind in TYPING_VIEWS
        self.typing_kind = self.comparison_mode = stored_kind if stored_kind in TYPING_VIEWS else 'mlst'
        self._load_kind_state()
        self.advanced_toggle.blockSignals(True)
        self.advanced_toggle.setChecked(bool(self.project.get_setting('compare.advanced_open', False)))
        self.advanced_toggle.blockSignals(False)
        self.set_advanced_visible(self.advanced_toggle.isChecked())
        self.counterpart_toggle.blockSignals(True)
        self.counterpart_toggle.setChecked(bool(self.project.get_setting('compare.counterpart_graph', False)))
        self.counterpart_toggle.blockSignals(False)
        active = self.project.get_setting('investigations.active', None)
        self.active_investigation_id = active if active in {p['id'] for p in store.list()} else None
        if not self.active_investigation_id and self.cohort_ids is None:
            self.cohort_ids = set()
        self.refresh_cohort_table()
        self._restore_kind_scheme()
        self._populate_investigation_combo()
        if self.active_investigation_id:
            self.select_investigation(self.active_investigation_id)
        self.set_dual_graph(self.dual_toggle.isChecked())
        self.set_counterpart_graph(self.counterpart_toggle.isChecked())

    def _populate_investigation_combo(self):
        """List the investigations of the typing view on screen, plus the active one.

        An investigation pins one reference, so it belongs to one typing kind.
        Offering a cgMLST investigation inside the MLST tree would silently swap
        the quantity being compared, so those entries are simply not listed here.
        """
        store = InvestigationStore(self.project)
        self.investigation_combo.blockSignals(True)
        self.investigation_combo.clear()
        self.investigation_combo.addItem('Unsaved investigation', None)
        for plan in store.list():
            kind = self._digest_kinds.get(plan.get('scheme_digest'))
            if kind in TYPING_VIEWS and kind != self.typing_kind and plan['id'] != self.active_investigation_id:
                continue
            suffix = '' if kind in TYPING_VIEWS else ' · reference not in this project'
            self.investigation_combo.addItem(plan['name'] + suffix, plan['id'])
        self.investigation_combo.setCurrentIndex(max(0, self.investigation_combo.findData(self.active_investigation_id)))
        self.investigation_combo.blockSignals(False)

    def _investigation_chosen(self, *_args):
        self.select_investigation(self.investigation_combo.currentData())

    def _follow_investigation_typing_kind(self, investigation_id):
        """An investigation pins one reference, so it belongs to one typing view.

        Opening it opens that view instead of drawing its cohort as an empty tree
        under the other kind's threshold.
        """
        try:
            plan = InvestigationStore(self.project).get(investigation_id)
        except KeyError:
            return
        kind = self._digest_kinds.get(plan.get('scheme_digest'))
        if kind not in TYPING_VIEWS or kind == self.typing_kind:
            return
        self._store_kind_state()
        self.typing_kind = self.comparison_mode = kind
        self._typing_kind_chosen = True
        self.project.set_setting('compare.typing_kind', kind)
        self._load_kind_state()
        self._refresh_typing_menu()

    def select_investigation(self, investigation_id):
        if investigation_id:
            self._follow_investigation_typing_kind(investigation_id)
        self.active_investigation_id = investigation_id
        self.project.set_setting('investigations.active', investigation_id)
        self._comparison_cache = None
        self._current_snapshot = None
        self._reset_baseline_state()
        if not investigation_id:
            self.refresh_baseline_graph()
            return
        plan = InvestigationStore(self.project).get(investigation_id)
        self.threshold_evidence = plan.get('threshold_evidence', {})
        known = {s['id'] for s in self.project.samples()}
        self.cohort_ids = set(plan['sample_ids']) & known
        self.compare_genus.blockSignals(True)
        self.compare_species.blockSignals(True)
        self.compare_genus.setCurrentIndex(max(0, self.compare_genus.findData(plan.get('filters', {}).get('genus', ''))))
        self.compare_species.setCurrentIndex(max(0, self.compare_species.findData(plan.get('filters', {}).get('species', ''))))
        self.compare_genus.blockSignals(False)
        self.compare_species.blockSignals(False)
        self.refresh_cohort_table()
        choice = 'digest:' + plan['scheme_digest'] if plan['scheme_digest'] else plan.get('scheme_path')
        index = self.compare_scheme.findData(choice)
        if index < 0 and choice:
            self.compare_scheme.addItem(f"Pinned: {plan['scheme'] or 'Reference'} · {plan['scheme_digest'][:8]}", choice)
            index = self.compare_scheme.findData(choice)
        if index >= 0:
            self.compare_scheme.blockSignals(True)
            self.compare_scheme.setCurrentIndex(index)
            self.compare_scheme.blockSignals(False)
        for control, value in [(self.overlap, plan['min_overlap']), (self.cluster_threshold, plan['threshold'])]:
            control.blockSignals(True)
            control.setValue(value)
            control.blockSignals(False)
        self.tree.set_label_fields(plan.get('label_fields') or ['sample_name', 'primary_st'])
        self.cohort_count.setText(f"{plan['name']} · {len(self.cohort_ids)} isolates")
        self.refresh_comparison()

    def _investigation_dialog(self, existing=None):
        current = self._last_comparison[0] if self._last_comparison else {}
        plan = existing or {}
        dialog = QDialog(self)
        dialog.setWindowTitle('Investigation protocol')
        layout = QVBoxLayout(dialog)
        layout.addWidget(label('Pin the reference and a locally justified threshold. Similarity is not a transmission diagnosis.', 'small', True))
        form = QFormLayout()
        name = QLineEdit(plan.get('name', 'New investigation'))
        protocol = QLineEdit(plan.get('protocol', 'Exploratory local threshold; not clinically validated'))
        threshold = QSpinBox()
        threshold.setRange(0, 100000)
        threshold.setValue(self.cluster_threshold.value())
        overlap = QDoubleSpinBox()
        overlap.setRange(0.01, 1.0)
        overlap.setSingleStep(0.05)
        overlap.setValue(self.overlap.value())
        include_new = QCheckBox('Offer new matching isolates when updating')
        include_new.setChecked(plan.get('include_new', True))
        for title, widget in [('Name', name), ('Protocol / threshold rationale', protocol),
                              ('Maximum link distance', threshold), ('Minimum shared / total loci', overlap)]:
            form.addRow(title, widget)
        layout.addLayout(form)
        layout.addWidget(include_new)
        choice = self.compare_scheme.currentData()
        digest = choice.removeprefix('digest:') if choice and choice.startswith('digest:') else ''
        matches = not choice or digest == current.get('scheme_digest') or choice == current.get('scheme_path')
        digest = digest or (current.get('scheme_digest', '') if matches else '')
        scheme = current.get('scheme', '') if matches else ''
        if digest == plan.get('scheme_digest'):
            scheme = scheme or plan.get('scheme', '')
        scheme_path = choice if choice and not choice.startswith('digest:') else plan.get('scheme_path', '')
        layout.addWidget(label(f"Reference: {scheme or 'Pending selected reference'}\nSHA-256: {digest or 'Not yet called'}", 'small', True))
        buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Save | QDialogButtonBox.StandardButton.Cancel)
        buttons.accepted.connect(dialog.accept)
        buttons.rejected.connect(dialog.reject)
        layout.addWidget(buttons)
        if dialog.exec() != QDialog.DialogCode.Accepted:
            return
        ids = self.cohort_ids if self.cohort_ids is not None else {s['id'] for s in self.project.samples()}
        try:
            saved = InvestigationStore(self.project).save(name.text(), ids, investigation_id=plan.get('id'),
                scheme_digest=digest, scheme=scheme, scheme_path=scheme_path, threshold=threshold.value(),
                min_overlap=overlap.value(), protocol=protocol.text(), include_new=include_new.isChecked(),
                filters={'genus': self.compare_genus.currentData(), 'species': self.compare_species.currentData()},
                label_fields=self.tree.label_fields, threshold_evidence=getattr(self, 'threshold_evidence', {}))
            self.project.set_setting('investigations.active', saved['id'])
            self.restore_investigations()
        except Exception as error:
            self.error(error)

    def create_investigation_dialog(self):
        self._investigation_dialog()

    def edit_investigation_dialog(self):
        if not self.active_investigation_id:
            return self.create_investigation_dialog()
        self._investigation_dialog(InvestigationStore(self.project).get(self.active_investigation_id))

    def update_investigation(self):
        if not self.active_investigation_id:
            self.notify('Save an investigation first to reuse its reference, protocol and previous cohort.')
            return self.create_investigation_dialog()
        store = InvestigationStore(self.project)
        plan = store.get(self.active_investigation_id)
        candidates = []
        if plan.get('include_new'):
            for sample in self.project.samples():
                if sample['id'] in plan['sample_ids'] or sample.get('metadata', {}).get('workflow', {}).get('source_kind') == 'read_mate':
                    continue
                genus, species, _ = organism_for(sample)
                filters = plan.get('filters', {})
                if filters.get('genus') and genus != filters['genus'] or filters.get('species') and species != filters['species']:
                    continue
                candidates.append(sample['id'])
        if candidates:
            if QMessageBox.question(self, 'Review investigation update',
                    f"Add {len(candidates)} new matching isolates to {plan['name']}?\n"
                    'Existing profiles and unchanged pairwise distances are reused. Missing profiles remain visibly unassessed until called. '
                    'Frozen review groups keep their original membership.') != QMessageBox.StandardButton.Yes:
                return
            plan = store.save(plan['name'], [*plan['sample_ids'], *candidates], investigation_id=plan['id'],
                **{key: plan[key] for key in ('scheme_digest', 'scheme', 'scheme_path', 'threshold', 'min_overlap',
                                             'protocol', 'include_new', 'filters', 'label_fields')})
        self.select_investigation(plan['id'])

    def refresh_cluster_table(self):
        if not hasattr(self, 'cluster_table'):
            return
        snapshot = self._current_snapshot or {}
        groups = snapshot.get('groups', [])
        # Group rows carry distances, so they carry the target set those distances
        # were counted over; a "2" here is not a "2" in the other typing view.
        units = (f"{snapshot.get('target_loci')} "
                 f"{TYPING_SCALES.get(snapshot.get('typing_kind') or '', TYPING_SCALES['unclassified'])['target_word']}"
                 if snapshot.get('target_loci') else 'loci')
        scale = snapshot.get('scale_caption') or self.typing_view_title()
        self._cluster_filling = True
        self.cluster_table.setSortingEnabled(False)
        self.cluster_table.setRowCount(len(groups))
        for row, group in enumerate(groups):
            title = group['name'] + (' · chaining' if group.get('chained') else '')
            for col, value in enumerate((title, len(group['members']), group.get('change', ''))):
                item = cell(value, group['id'])
                item.setToolTip(f"{scale}\n"
                                f"{group['status']} · {group.get('reason', '')}\n"
                                f"Maximum direct distance: {group.get('max_direct_distance')} "
                                f"allele differences over {units}\n"
                                f"Unassessed within-group pairs: {group.get('unassessed_within_pairs', 0)}\n"
                                f"Membership change: {group.get('change')}\nClick to select all group members.")
                self.cluster_table.setItem(row, col, item)
        self.cluster_table.setSortingEnabled(True)
        self._cluster_filling = False

    def select_cluster_row(self):
        if self._cluster_filling or not self._current_snapshot:
            return
        ids = {item.data(Qt.ItemDataRole.UserRole) for item in self.cluster_table.selectedItems()}
        members = {sid for group in self._current_snapshot['groups'] if group['id'] in ids for sid in group['members']}
        if members:
            self.tree.select_ids(members)

    def select_graph_cluster(self):
        ids = self.tree.selected_ids()
        if ids:
            self.tree.select_cluster(ids[0])

    def graph_selection_changed(self, ids, role='current'):
        if not hasattr(self, 'graph_selection_label'):
            return
        self.sync_graph_selection(ids, role)
        if role == 'current':
            self.graph_report_button.setEnabled(bool(ids))
        snapshot = {'baseline': self._baseline_snapshot,
                    'counterpart': self._counterpart_snapshot}.get(role, self._current_snapshot)
        where = {'baseline': ' in the baseline tree',
                 'counterpart': f' in the {typing_title(self.counterpart_kind())} tree'}.get(role, '')
        text = f"{len(ids)} selected" + where
        if len(ids) == 1 and snapshot:
            rows = proximity_rows(snapshot, ids)
            if rows:
                nearest = rows[0]['nearest_distance']
                text += f" · nearest: {nearest if nearest is not None else 'not comparable'} allele differences"
        text += self._selection_overlap(ids, role)
        self.graph_selection_label.setText(text + ' · Ctrl-click adds; click a cluster row selects its members.')

    def _selection_overlap(self, ids, role):
        """Say plainly how much of a selection each other open tree contains."""
        if not ids:
            return ''
        keys = set(map(str, ids))
        parts = []
        if self.dual_toggle.isChecked() and role in {'current', 'baseline'}:
            other = self._graph_view('current' if role == 'baseline' else 'baseline')
            if getattr(other, '_results', None):
                shared = len(keys & set(other._results))
                parts.append(f" · {shared} still in the current comparison · {len(ids) - shared} no longer in it"
                             if role == 'baseline' else
                             f" · {shared} also in the baseline tree · {len(ids) - shared} added since the baseline")
        if self.counterpart_toggle.isChecked() and role in {'current', 'counterpart'}:
            other = self._graph_view('current' if role == 'counterpart' else 'counterpart')
            named = self.typing_view_title() if role == 'counterpart' else typing_title(self.counterpart_kind())
            shared = len(keys & set(getattr(other, '_results', {}) or {}))
            # A missing profile of the other kind is unknown evidence about that
            # isolate; it is never reported as closeness or as a zero distance.
            parts.append(f" · {shared} also have a {named} profile · "
                         f"{len(ids) - shared} have no {named} profile")
        return ''.join(parts)

    def sync_graph_selection(self, ids, role):
        """Mirror a selection into every other open tree; it is the same isolate."""
        if self._syncing_graph_selection or not self.link_selection.isChecked():
            return
        keys = set(map(str, ids))
        self._syncing_graph_selection = True
        try:
            for other_role, view in self._graph_views().items():
                if other_role == role or view is None or not self._view_active(other_role):
                    continue
                if not getattr(view, '_results', None):
                    continue
                view.select_ids(keys & set(view._results))
        finally:
            self._syncing_graph_selection = False

    def highlight_graphs(self, text):
        """Highlight matching isolates in every open tree; a view aid, never a filter."""
        counts = {}
        for role, view in self._graph_views().items():
            if view is not None and hasattr(view, 'highlight'):
                counts[role] = len(view.highlight(text))
        if not str(text).strip():
            return counts
        message = f"{counts.get('current', 0)} matched here"
        if self.dual_toggle.isChecked():
            message += f" · {counts.get('baseline', 0)} in the baseline tree"
        if self.counterpart_toggle.isChecked():
            message += (f" · {counts.get('counterpart', 0)} in the "
                        f"{typing_title(self.counterpart_kind())} tree")
        self.graph_selection_label.setText(message + ' · highlighting changes no stored evidence.')
        return counts

    def report_graph_selection(self, ids=None, role='current'):
        view = self._graph_view(role)
        ids = set(ids if isinstance(ids, (list, tuple, set)) else view.selected_ids())
        if not ids:
            self.notify('Select graph nodes, labels or a cluster row first.')
            return
        self.report_ids = ids
        self.report_scope.setCurrentIndex(0)
        # A report sourced from another tree must carry that tree's snapshot, or
        # its distances would not be the ones in the picture beside them — and a
        # cgMLST distance would be printed under a classical MLST heading.
        self._report_investigation_snapshot = self.export_snapshot(role)
        self.refresh_report_table()
        self.navigate(5)

    def report_isolate_proximity(self, sid, role='current'):
        self._graph_view(role).select_ids([sid])
        if hasattr(self, 'report_preset'):
            self.report_preset.setCurrentIndex(self.report_preset.findData('proximity'))
        self.report_graph_selection([sid], role=role)

    def freeze_investigation_selection(self):
        ids = self.tree.selected_ids()
        if not self.active_investigation_id or not ids:
            self.notify('Open a saved investigation and select graph isolates before freezing a review group.')
            return
        name, accepted = QInputDialog.getText(self, 'Freeze review membership', 'Review group name:')
        if accepted and name.strip():
            InvestigationStore(self.project).freeze_group(self.active_investigation_id, name.strip(), ids)
            if self._current_snapshot:
                self._current_snapshot['review_groups'] = InvestigationStore(self.project).get(self.active_investigation_id)['review_groups']
            self.notify(f'Frozen {len(ids)} isolate IDs. Future threshold changes will not alter this review group.')

    def choose_frozen_review_group(self):
        if not self.active_investigation_id:
            self.notify('Open a saved investigation first.')
            return
        groups = InvestigationStore(self.project).get(self.active_investigation_id)['review_groups']
        if not groups:
            self.notify('No frozen review groups yet. Select isolates, then Freeze review group.')
            return
        names = [f"{g['name']} · {len(g['sample_ids'])} isolates · {g['id'][:8]}" for g in groups]
        selected, accepted = QInputDialog.getItem(self, 'Open frozen review group', 'Fixed membership (unchanged by thresholds):', names, 0, False)
        if accepted:
            group = groups[names.index(selected)]
            self.tree.select_ids([sid for sid in group['sample_ids'] if sid in self.tree._results])
            missing = set(group['sample_ids']) - self.tree._results.keys()
            self.notify(f"Selected {len(group['sample_ids']) - len(missing)} frozen members; {len(missing)} are outside the current graph.")

    def open_investigation_history(self, snapshot_id=None):
        """The frozen evidence behind any stored snapshot, optionally opened at one."""
        if not self.active_investigation_id:
            self.notify('Save an investigation before opening its snapshot history.')
            return
        from wmlstudio.export import _atomic_text, investigation_html
        store = InvestigationStore(self.project)
        plan = store.get(self.active_investigation_id)
        if not plan['snapshots']:
            self.notify('Build the saved investigation once to create its first snapshot.')
            return
        dialog = QDialog(self)
        dialog.setWindowTitle('Immutable investigation snapshots')
        dialog.resize(900, 700)
        layout = QVBoxLayout(dialog)
        layout.addWidget(label('Past membership, thresholds and pairwise evidence are frozen. Viewing or exporting history does not replace the current cohort.', 'small', True))
        chooser = QComboBox()
        baseline = store.baseline_snapshot_id(plan['id'])
        for number, sid in reversed(list(enumerate(plan['snapshots'], 1))):
            chooser.addItem(f'Snapshot {number} · {sid[:12]}'
                            + (' · baseline tree' if sid == baseline else ''), sid)
        if isinstance(snapshot_id, str) and chooser.findData(snapshot_id) >= 0:
            chooser.setCurrentIndex(chooser.findData(snapshot_id))
        layout.addWidget(chooser)
        view = QTextBrowser()
        layout.addWidget(view, 1)
        def refresh():
            snapshot = store.snapshot(plan['id'], chooser.currentData())
            view.setHtml(investigation_html(snapshot))
        chooser.currentIndexChanged.connect(refresh)
        def export():
            path, _ = QFileDialog.getSaveFileName(dialog, 'Export immutable snapshot', 'investigation-snapshot.json', 'JSON (*.json);;HTML (*.html)')
            if not path:
                return
            try:
                self.check_output(path)
                snapshot = store.snapshot(plan['id'], chooser.currentData())
                with _atomic_text(path) as handle:
                    if Path(path).suffix.lower() == '.html':
                        handle.write('<!doctype html><html><head><meta charset="utf-8"></head><body>' + investigation_html(snapshot) + '</body></html>')
                    else:
                        json.dump(snapshot, handle, ensure_ascii=False, indent=2, allow_nan=False)
                self.notify('Exported the chosen historical snapshot without changing the current investigation.')
            except Exception as error:
                self.error(error)
        layout.addWidget(button('Export this snapshot…', export))
        refresh()
        dialog.exec()

    def edit_graph_label_fields(self):
        dialog = QDialog(self)
        dialog.setWindowTitle('Graph label fields')
        layout = QVBoxLayout(dialog)
        layout.addWidget(label('Choose up to eight fields. Primary MLST/ST is retained when the graph uses cgMLST. AMR is annotation only.', 'small', True))
        fields = QListWidget()
        for name in self.tree.available_label_fields():
            item = QListWidgetItem(name.removeprefix('metadata:'))
            item.setData(Qt.ItemDataRole.UserRole, name)
            item.setFlags(item.flags() | Qt.ItemFlag.ItemIsUserCheckable)
            item.setCheckState(Qt.CheckState.Checked if name in self.tree.label_fields else Qt.CheckState.Unchecked)
            fields.addItem(item)
        layout.addWidget(fields)
        buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel)
        buttons.accepted.connect(dialog.accept)
        buttons.rejected.connect(dialog.reject)
        layout.addWidget(buttons)
        if dialog.exec() == QDialog.DialogCode.Accepted:
            try:
                self.tree.set_label_fields([fields.item(i).data(Qt.ItemDataRole.UserRole) for i in range(fields.count())
                                            if fields.item(i).checkState() == Qt.CheckState.Checked])
                self.persist_graph_style()
                if self.active_investigation_id:
                    store = InvestigationStore(self.project)
                    plan = store.get(self.active_investigation_id)
                    store.save(plan['name'], plan['sample_ids'], investigation_id=plan['id'],
                        **{key: plan[key] for key in ('scheme_digest', 'scheme', 'scheme_path', 'threshold', 'min_overlap', 'protocol', 'include_new', 'filters')},
                        label_fields=self.tree.label_fields)
            except ValueError as error:
                self.error(error)

    def choose_comparison_cohort(self):
        from wmlstudio.cohort_picker import CohortPickerDialog
        dialog = CohortPickerDialog(self.project.samples(), self.project, 'Choose comparison isolates',
                                    selected_ids=self.cohort_ids or set(), parent=self)
        if dialog.exec() != QDialog.DialogCode.Accepted:
            return
        self.cohort_ids = set(dialog.selected_ids)
        self.compare_genus.setCurrentIndex(0)
        self.compare_species.setCurrentIndex(0)
        self.project.set_setting('comparison_cohort', sorted(self.cohort_ids))
        # A changed explicit cohort is a preview until the saved plan is edited.
        if self.active_investigation_id:
            self.notify('Cohort changed. Save the investigation protocol to record this membership.')
        self.refresh_cohort_table()
        self.refresh_comparison()

    def suggested_scheme_key(self):
        """The curated publication key this view's reference declares, or ''.

        Only an explicit declaration is accepted — from the stored profiles, from
        the installed reference's own metadata, or from a binding a person
        has already attested in the guidance dialog. Nothing is matched by name.
        """
        return self._reference_binding()['scheme_key']

    def _reference_binding(self):
        """The catalogue pin and publication key of the reference now on screen."""
        results = getattr(self, '_last_comparison', []) or []
        binding = {'catalog_key': '',
                   'scheme_key': scheme_key_from(*results[:1],
                                                 getattr(self, 'threshold_evidence', {}) or {})}
        paths = [self.compare_scheme.currentData()]
        paths.extend(result.get('scheme_path') for result in results[:1])
        for path in paths:
            if not path or str(path).startswith('digest:'):
                continue
            facts = self.comparison_scheme_facts(path)
            binding['catalog_key'] = binding['catalog_key'] or facts.get('catalog_key', '')
            binding['scheme_key'] = binding['scheme_key'] or facts.get('scheme_key', '')
        return binding

    def threshold_suggestion(self):
        """The published cutoff bound to this exact reference; never applied for you."""
        snapshot = self._current_snapshot or {}
        binding = self._reference_binding()
        return scheme_threshold_suggestion(
            binding['scheme_key'], method=self.typing_kind, targets=snapshot.get('target_loci'),
            scheme=snapshot.get('scheme', ''), catalog_key=binding['catalog_key'])

    def refresh_threshold_suggestion(self):
        if not hasattr(self, 'guidance_row'):
            return
        suggestion = self.threshold_suggestion()
        self._threshold_suggestion = suggestion
        visible = bool(self._current_snapshot) and suggestion['status'] in {'available', 'target_count_mismatch'}
        if visible:
            self.guidance_banner.setText(suggestion['headline'] + ' ' + suggestion['reason'])
            self.guidance_banner.setToolTip('\n'.join(filter(None, [
                suggestion['headline'], suggestion['reason'],
                'Catalogue ' + suggestion['catalog_version'] + ', reviewed ' + suggestion['reviewed_on'] + '.',
                *[f"{entry['scope']} — {entry['source']['citation']} ({entry['source']['url']})"
                  for entry in suggestion['entries']],
                suggestion['interpretation'], SCALE_SEPARATION])))
            self.guidance_button.setEnabled(suggestion['status'] == 'available')
        self.guidance_row.setVisible(visible)
        self._refresh_threshold_units()

    def apply_threshold_guidance(self, evidence):
        """Called after explicit reference/protocol review in the guidance dialog."""
        from copy import deepcopy
        if evidence.get('approved_threshold') is not None:
            context = evidence.get('context') or {}
            snapshot = self._current_snapshot or {}
            method = str(context.get('method') or '').strip().casefold()
            if method and method != self.typing_kind:
                raise ValueError(
                    f'That cutoff was reviewed for a {typing_title(method)} comparison and this tree is a '
                    f'{self.typing_view_title()} comparison. ' + SCALE_SEPARATION)
            if context.get('locus_count') is not None and snapshot.get('target_loci') is not None and (
                    context['locus_count'] != snapshot['target_loci']):
                raise ValueError('The reviewed target count no longer matches the displayed comparison. '
                                 'A published cutoff applies to its full target set, with no scaling.')
            choice = self.compare_scheme.currentData()
            digest = context.get('scheme_digest')
            matches = digest and digest == snapshot.get('scheme_digest')
            if choice and choice.startswith('digest:'):
                matches = matches and choice[7:] == digest
            elif choice:
                matches = matches and any(r.get('scheme_path') == choice and r.get('scheme_digest') == digest for r in self._last_comparison)
            if not matches or context.get('min_overlap') != self.overlap.value():
                raise ValueError('The reviewed reference/overlap no longer matches the displayed comparison. Rebuild and review the publication context again.')
            value = evidence['approved_threshold']
            if isinstance(value, bool) or not isinstance(value, int) or value < 0:
                raise ValueError('Reviewed threshold must be a nonnegative integer.')
            self.cluster_threshold.setValue(value)
        self.threshold_evidence = deepcopy(evidence)
        if self._current_snapshot:
            self._current_snapshot = deepcopy(self._current_snapshot)
            self._current_snapshot.update(preview=True, threshold_evidence=deepcopy(evidence))
        self.notify('Guidance recorded for this comparison preview. Save the investigation protocol to retain its citation and justification.')

    def show_threshold_guidance(self):
        callback = getattr(self, 'open_threshold_guidance', None)
        if callback:
            return callback()
        self.notify('Published guidance requires an exact scheme/protocol review. No cutoff was changed.')

    def build_compare(self):
        page, layout = self.page()
        # Dense workbench controls need breathing room at 1080×720 with native
        # Windows font metrics. Keep text size intact; reduce padding rather
        # than letting action strips consume the graph's usable height.
        page.setStyleSheet('QPushButton, QToolButton { padding: 6px 10px; } '
                           'QComboBox, QLineEdit, QSpinBox, QDoubleSpinBox '
                           '{ padding-top: 5px; padding-bottom: 5px; }')
        layout.setSpacing(8)
        # The tab's orientation strip already carries the purpose sentence and the
        # "Choose cohort…" step, so only the title is repeated here. The ~40 px this
        # saves is what lets two trees each keep a usable height at 1080x720.
        page_title = label("Investigate related isolates", "title")
        page_title.setToolTip('Core allele evidence and accessory evidence such as AMR genes are '
                              'reported separately; accessory findings are never added into the '
                              'core allele distance.')
        layout.addWidget(page_title)
        self.comparison_settings = QDialog(self)
        self.comparison_settings.setWindowTitle('Comparison settings')
        self.comparison_settings.resize(760, 340)
        settings_layout = QVBoxLayout(self.comparison_settings)
        settings_layout.addWidget(label('Reference and locally justified comparison policy. Existing compatible profiles are reused; Call missing profiles reads only inputs needing this reference.', 'small', True))
        filters = QHBoxLayout()
        self.compare_genus = QComboBox()
        self.compare_genus.setProperty('compactCharacters', 6)
        self.compare_genus.addItem("All genera", "")
        self.compare_species = QComboBox()
        self.compare_species.setProperty('compactCharacters', 6)
        self.compare_species.addItem("All species", "")
        self.compare_scheme = QComboBox()
        self.compare_scheme.setProperty('compactCharacters', 12)
        self.compare_scheme.addItem("Existing profile snapshot", None)
        self.compare_scheme.setMinimumWidth(160)
        self.compare_scheme.setSizeAdjustPolicy(QComboBox.SizeAdjustPolicy.AdjustToMinimumContentsLengthWithIcon)
        self.compare_scheme.setMinimumContentsLength(12)
        self.compare_genus.currentIndexChanged.connect(self.refresh_cohort_table)
        self.compare_species.currentIndexChanged.connect(self.refresh_cohort_table)
        self.compare_scheme.currentIndexChanged.connect(self.comparison_scheme_changed)
        filters.addWidget(self.compare_genus, 1)
        filters.addWidget(self.compare_species, 1)
        filters.addWidget(self.compare_scheme, 2)
        filters.addWidget(button("Call missing profiles…", self.type_comparison_scheme))
        settings_layout.addLayout(filters)
        controls = FlowLayout()
        self.overlap = QDoubleSpinBox()
        self.overlap.setRange(0.01, 1.0)
        self.overlap.setDecimals(2)
        self.overlap.setSingleStep(0.05)
        self.overlap.setValue(0.95)
        self.overlap.setToolTip("Minimum shared called loci divided by all loci. Missing data never becomes a matching allele.")
        # The link threshold is the one number a novice must see, so it lives on
        # the page beside the tree rather than inside this dialog; its suffix
        # always names the target set the number is counted over.
        self.cluster_threshold = QSpinBox()
        self.cluster_threshold.setRange(0, 100000)
        self.cluster_threshold.setValue(1)
        self.cluster_threshold.valueChanged.connect(self._remember_threshold)
        self.cluster_threshold.valueChanged.connect(self.refresh_comparison)
        self.overlap.valueChanged.connect(self._remember_overlap)
        overlap_row = QWidget()
        overlap_layout = QHBoxLayout(overlap_row)
        overlap_layout.setContentsMargins(0, 0, 0, 0)
        overlap_layout.addWidget(label('Shared ≥', 'small'))
        overlap_layout.addWidget(self.overlap)
        controls.addWidget(overlap_row)
        controls.addWidget(button("Apply comparison", lambda: (self.comparison_settings.hide(), self.refresh_comparison()), True))
        settings_layout.addLayout(controls)
        chooser = QWidget(self)
        chooser.hide()  # Legacy selection model; actual per-instance picker is a dialog.
        left = QVBoxLayout(chooser)
        left.setContentsMargins(0, 0, 8, 0)
        self.cohort_count = label("Choose a cohort", "cardTitle")
        left.addWidget(self.cohort_count)
        investigations = QHBoxLayout()
        self.investigation_combo = QComboBox()
        self.investigation_combo.setProperty('compactCharacters', 14)
        self.investigation_combo.addItem('Unsaved investigation', None)
        self.investigation_combo.currentIndexChanged.connect(self._investigation_chosen)
        investigations.addWidget(self.investigation_combo, 1)
        options = QToolButton()
        options.setText('…')
        options.setToolTip('Save, reopen, update or freeze an investigation review group')
        options.setPopupMode(QToolButton.ToolButtonPopupMode.InstantPopup)
        investigation_menu = QMenu(options)
        investigation_menu.addAction('Save as new investigation…', self.create_investigation_dialog)
        investigation_menu.addAction('Edit investigation protocol…', self.edit_investigation_dialog)
        investigation_menu.addAction('Update with new project isolates…', self.update_investigation)
        investigation_menu.addAction('Review snapshot history…', self.open_investigation_history)
        investigation_menu.addAction('Freeze selected review group…', self.freeze_investigation_selection)
        investigation_menu.addAction('Open frozen review group…', self.choose_frozen_review_group)
        investigation_menu.addAction('Label fields…', self.edit_graph_label_fields)
        options.setMenu(investigation_menu)
        investigations.addWidget(options)
        investigations.addWidget(button('Choose cohort…', self.choose_comparison_cohort, True))
        # The separate MLST and cgMLST trees are chosen here, and the button says
        # which one is on screen. They are never merged into one "analysis": a
        # classical ST distance and a core-genome distance are different
        # quantities and must not share a view, a threshold or a scale.
        self.typing_menu = QToolButton()
        self.typing_menu.setPopupMode(QToolButton.ToolButtonPopupMode.InstantPopup)
        typing_menu = QMenu(self.typing_menu)
        typing_menu.addAction('Classical MLST tree (sequence type)…', lambda: self.show_typing_view('mlst'))
        typing_menu.addAction('cgMLST tree (core-genome targets)…', lambda: self.show_typing_view('cgmlst'))
        typing_menu.addSeparator()
        typing_menu.addAction('SNP follow-up…', lambda: self.set_comparison_mode('snp'))
        self.typing_menu.setMenu(typing_menu)
        investigations.addWidget(self.typing_menu)
        investigations.addWidget(button('Compare', self.refresh_comparison))
        layout.addLayout(investigations)
        self.cohort_search = QLineEdit()
        self.cohort_search.setPlaceholderText("Find isolates…")
        self.cohort_search.textChanged.connect(self.refresh_cohort_table)
        left.addWidget(self.cohort_search)
        self.show_all_comparison_schemes = QCheckBox('Show all schemes')
        self.show_all_comparison_schemes.setToolTip(
            'Include schemes for other organisms. A scheme whose organism is unrecorded stays listed '
            'either way. References of the other typing kind are never listed here: they belong to '
            'the other tree.')
        self.show_all_comparison_schemes.toggled.connect(self.refresh_cohort_table)
        settings_layout.addWidget(self.show_all_comparison_schemes)
        settings_layout.addWidget(button('Graph label fields…', self.edit_graph_label_fields))
        settings_layout.addWidget(button('Published threshold guidance…', self.show_threshold_guidance))
        settings_layout.addWidget(label('Single linkage may chain A–B–C even when A–C exceeds the link threshold. Review maximum direct distances and epidemiological evidence; no universal clinical cutoff is assumed.', 'small', True))
        cohort_buttons = FlowLayout()
        cohort_buttons.addWidget(button('All visible', lambda: self.set_cohort_visible(True)))
        cohort_buttons.addWidget(button('Clear', lambda: self.set_cohort_visible(False)))
        left.addLayout(cohort_buttons)
        self.cohort_table = make_table(["Use", "Sample", "Organism", "Stored profile"])
        self.cohort_table.setColumnWidth(0, 48)
        self.cohort_table.itemChanged.connect(self.cohort_item_changed)
        self.cohort_tabs = QTabWidget()
        self.cohort_tabs.addTab(self.cohort_table, 'Isolates')
        self.cluster_table = make_table(['Group', 'N', 'Change'])
        self.cluster_table.itemSelectionChanged.connect(self.select_cluster_row)
        self.cohort_tabs.addTab(self.cluster_table, 'Clusters')
        left.addWidget(self.cohort_tabs, 1)
        left.addWidget(button("Report this cohort", self.report_comparison))
        right = QWidget()
        right.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Ignored)
        content = QVBoxLayout(right)
        content.setContentsMargins(0, 0, 0, 0)
        content.setSpacing(6)
        # Default strip: the one threshold that defines the groups, and a door to
        # everything else. This page was the most crowded in the application, so
        # only what a first-time user needs is shown before Advanced is opened.
        toolbar = FlowLayout()
        threshold_row = QWidget()
        threshold_layout = QHBoxLayout(threshold_row)
        threshold_layout.setContentsMargins(0, 0, 0, 0)
        threshold_layout.addWidget(label('Group ≤', 'small'))
        threshold_layout.addWidget(self.cluster_threshold)
        toolbar.addWidget(threshold_row)
        self.advanced_toggle = QToolButton()
        self.advanced_toggle.setCheckable(True)
        self.advanced_toggle.setText('Advanced ▸')
        self.advanced_toggle.setToolTip('Baseline tree, the other typing view, colours, highlighting and exports.')
        self.advanced_toggle.toggled.connect(self.set_advanced_visible)
        toolbar.addWidget(self.advanced_toggle)
        self.advanced_panel = QWidget()
        advanced = FlowLayout(self.advanced_panel)
        self.color_by = QComboBox()
        self.color_by.setProperty('compactCharacters', 9)
        self.color_by.addItem("Colour: cluster", "cluster")
        self.color_by.addItem("Colour: ST", "st")
        self.color_by.currentIndexChanged.connect(self.set_graph_color_by)
        graph_options = QToolButton()
        graph_options.setText('Graph options')
        graph_options.setPopupMode(QToolButton.ToolButtonPopupMode.InstantPopup)
        menu = QMenu(graph_options)
        menu.addAction('Fit graph', lambda: self.tree.fit_tree())
        menu.addAction('Reference / threshold settings…', self.comparison_settings.show)
        menu.addAction('Colour selected nodes…', self.color_graph_selection)
        menu.addAction('Highlight selected nodes…', self.highlight_graph_selection)
        menu.addAction('Report selected nodes / cluster…', self.report_graph_selection)
        menu.addAction('Select threshold group of selected node', self.select_graph_cluster)
        menu.addAction('Freeze selected review group…', self.freeze_investigation_selection)
        menu.addAction('Open frozen review group…', self.choose_frozen_review_group)
        menu.addAction('Label fields…', self.edit_graph_label_fields)
        menu.addSeparator()
        menu.addAction('Choose baseline snapshot…', self.choose_baseline_snapshot)
        menu.addAction('Pin the current comparison as the baseline', self.pin_current_as_baseline)
        menu.addAction('What changed since the baseline…', self.show_changes_tab)
        menu.addAction('View full baseline evidence…',
                       lambda: self.open_investigation_history(self.baseline_reference()[1]))
        # A checkable menu entry instead of a toolbar checkbox: the graph toolbar
        # wraps at 1080 px, and a wrapped row is taken straight out of the graph.
        self.link_selection = menu.addAction('Link selection across every open tree')
        self.link_selection.setToolTip('Selecting an isolate in one tree selects the same isolate in the '
                                       'others. It is the same isolate; the distances around it are not '
                                       'the same measurement.')
        self.link_selection.setCheckable(True)
        self.link_selection.setChecked(True)
        menu.addSeparator()
        graph_options.setMenu(menu)
        advanced.addWidget(graph_options)
        self.dual_toggle = QCheckBox('Compare with baseline tree')
        self.dual_toggle.setToolTip('Show the tree as it was first built beside the tree as it is now.\n'
                                    'Both are layouts of allele differences, not family trees, and two trees '
                                    'can differ because the cohort changed rather than because evidence did.')
        self.dual_toggle.toggled.connect(self.set_dual_graph)
        advanced.addWidget(self.dual_toggle)
        self.counterpart_toggle = QCheckBox('Show the other typing view')
        self.counterpart_toggle.setToolTip(
            'Show the same isolates typed the other way, in their own tree with their own reference, '
            'target count, threshold and legend.\n' + SCALE_SEPARATION)
        self.counterpart_toggle.toggled.connect(self.set_counterpart_graph)
        advanced.addWidget(self.counterpart_toggle)
        self.graph_search = QLineEdit()
        self.graph_search.setPlaceholderText('Highlight isolates…')
        self.graph_search.setMinimumWidth(120)
        self.graph_search.setToolTip('Highlight matching isolates in every open tree. Highlighting is a view '
                                     'aid; it changes no stored evidence and no group membership.')
        self.graph_search.textChanged.connect(self.highlight_graphs)
        advanced.addWidget(self.graph_search)
        advanced.addWidget(self.color_by)
        self.export_tree_choice = QComboBox()
        self.export_tree_choice.setProperty('compactCharacters', 9)
        self.export_tree_choice.addItem('Export: current tree', 'current')
        self.export_tree_choice.addItem('Export: baseline tree', 'baseline')
        self.export_tree_choice.addItem('Export: other typing view', 'counterpart')
        self.export_tree_choice.setToolTip('Which tree the export entries act on. Each export carries its own '
                                           'typing kind, reference and target count.')
        self.export_tree_choice.hide()
        advanced.addWidget(self.export_tree_choice)
        export = QComboBox()
        export.setProperty('compactCharacters', 9)
        export.addItems(["Export…", "PNG image", "SVG vector", "GraphML", "Newick (MST topology)",
                         "Distance JSON", "Distance matrix TSV", "JPEG image", "Profile matrix TSV",
                         "Group table TSV", "Baseline + current image (PNG)",
                         "Baseline + current image (JPEG)", "Change summary (TSV)",
                         "Change summary (JSON)"])
        export.activated.connect(lambda index: self.export_graph_action(index, export))
        advanced.addWidget(export)
        advanced.addWidget(button('Reference & overlap settings…', self.comparison_settings.show))
        content.addLayout(toolbar)
        self.advanced_panel.hide()
        content.addWidget(self.advanced_panel)
        content.addWidget(self._build_guidance_row())
        for text, method, default in [("Labels", "set_labels_visible", True), ("ST", "set_show_st", True),
                                      ("Edge distances", "set_edge_labels_visible", True),
                                      ("Merge identical", "set_merge_identical", False), ("Cluster halos", "set_halos_visible", True)]:
            check = menu.addAction(text)
            check.setCheckable(True)
            check.setChecked(default)
            check.toggled.connect(lambda enabled, name=method: self.graph_option(name, enabled))
        tabs = QTabWidget()
        tabs.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        self.tree = ComparisonTreeView()
        self.tree.setMinimumHeight(200)
        self.baseline_tree = ComparisonTreeView()
        self.baseline_tree.setMinimumHeight(200)
        self.counterpart_tree = ComparisonTreeView()
        self.counterpart_tree.setMinimumHeight(200)
        # The baseline half is hidden until it is asked for, so the single-tree
        # page is laid out exactly as it was before the second tree existed.
        self.graph_split = QSplitter(Qt.Orientation.Horizontal)
        self.graph_split.setChildrenCollapsible(False)
        self.counterpart_pane, counterpart_column = self._graph_pane()
        self.counterpart_caption = self._graph_caption()
        counterpart_column.addWidget(self.counterpart_caption)
        counterpart_column.addWidget(self.counterpart_tree, 1)
        self.counterpart_notice = label('', 'small', True)
        self.counterpart_notice.hide()
        counterpart_column.addWidget(self.counterpart_notice)
        # Its own legend, never merged into the other tree's: two kinds can both
        # have a "Cluster 001" and they are not the same group of isolates.
        self.counterpart_legend = QLabel()
        self.counterpart_legend.setWordWrap(True)
        self.counterpart_legend.setTextFormat(Qt.TextFormat.RichText)
        counterpart_column.addWidget(self.counterpart_legend)
        self.baseline_pane, baseline_column = self._graph_pane()
        self.baseline_caption = self._graph_caption()
        baseline_column.addWidget(self.baseline_caption)
        baseline_column.addWidget(self.baseline_tree, 1)
        self.baseline_notice = label('', 'small', True)
        self.baseline_notice.hide()
        baseline_column.addWidget(self.baseline_notice)
        self.baseline_action = button('', self._baseline_action_clicked)
        self.baseline_action.hide()
        baseline_column.addWidget(self.baseline_action)
        self.current_pane, current_column = self._graph_pane()
        self.current_caption = self._graph_caption()
        current_column.addWidget(self.current_caption)
        current_column.addWidget(self.tree, 1)
        self.graph_split.addWidget(self.counterpart_pane)
        self.graph_split.addWidget(self.baseline_pane)
        self.graph_split.addWidget(self.current_pane)
        self.counterpart_pane.hide()
        self.baseline_pane.hide()
        self.current_caption.hide()
        tabs.addTab(self.graph_split, "Graph")
        tabs.addTab(self.cluster_table, 'Groups')
        self.profile_table = QTableView()
        self.profile_model = EvidenceMatrixModel(self.profile_table)
        self.profile_table.setModel(self.profile_model)
        self.profile_table.setSortingEnabled(True)
        self.profile_table.setAlternatingRowColors(True)
        tabs.addTab(self.profile_table, "Profiles")
        self.statistics_view = QTextBrowser()
        tabs.addTab(self.statistics_view, "Quality / statistics")
        self.changes_view = QTextBrowser()
        tabs.addTab(self.changes_view, "Changes")
        self.graph_tabs = tabs
        if hasattr(self.pages, 'register_subtabs'):
            self.pages.register_subtabs('compare', tabs)
        content.addWidget(tabs, 1)
        self.graph_legend = QLabel()
        self.graph_legend.setWordWrap(True)
        self.graph_legend.setTextFormat(Qt.TextFormat.RichText)
        content.addWidget(self.graph_legend)
        selection_row = QHBoxLayout()
        self.graph_selection_label = label('Click nodes or their labels to select isolates.', 'small', True)
        selection_row.addWidget(self.graph_selection_label, 1)
        self.graph_report_button = button('Report selection', self.report_graph_selection)
        self.graph_report_button.setEnabled(False)
        selection_row.addWidget(self.graph_report_button)
        content.addLayout(selection_row)
        # QSplitter supplies a non-height-for-width viewport boundary. Without
        # it nested wrapped labels inflate QScrollArea's preferred graph height
        # even while the graph has ample room to shrink safely.
        graph_host = QSplitter(Qt.Orientation.Horizontal)
        graph_host.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Ignored)
        graph_host.setChildrenCollapsible(False)
        graph_host.addWidget(right)
        layout.addWidget(graph_host, 1)
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
        self.active_investigation_id = None
        self._current_snapshot = None
        self._cluster_filling = False
        # Which of the two comparisons this page is showing. comparison_mode is
        # the same value under the name the guidance dialog already reads, so a
        # published cgMLST cutoff can never be reviewed against an MLST tree.
        self.typing_kind = 'mlst'
        self.comparison_mode = 'mlst'
        self._typing_kind_chosen = False
        self._digest_kinds = {}
        self._kind_counts = dict.fromkeys(TYPING_VIEWS, 0)
        self._counterpart_snapshot = None
        self._counterpart_results = []
        self._counterpart_payload = None
        self._threshold_suggestion = {}
        self._legends = {'current': {}, 'baseline': {}, 'counterpart': {}}
        self._syncing_graph_selection = False
        self._pending_baseline_graph_state = None
        self._reset_baseline_state()
        self._comparison_timer = QTimer(self)
        self._comparison_timer.setSingleShot(True)
        self._comparison_timer.timeout.connect(self._start_pending_comparison)
        for view, role in ((self.tree, 'current'), (self.baseline_tree, 'baseline'),
                           (self.counterpart_tree, 'counterpart')):
            self._bind_graph_view(view, role)
        self.install_graph_menus()
        self._refresh_typing_menu()
        QTimer.singleShot(0, self.restore_investigations)

    def _graph_pane(self):
        """One half of the dual graph: a caption row above a tree, no extra margin."""
        pane = QWidget()
        column = QVBoxLayout(pane)
        column.setContentsMargins(0, 0, 0, 0)
        column.setSpacing(4)
        return pane, column

    def _graph_caption(self):
        caption = label('', 'small')
        # Ignored width: a long pane caption must never widen the comparison page
        # into a horizontal scroll. The full sentence stays in the tooltip.
        caption.setSizePolicy(QSizePolicy.Policy.Ignored, QSizePolicy.Policy.Preferred)
        return caption

    def _bind_graph_view(self, view, role):
        """Wire one tree, remembering which snapshot its signals speak for."""
        for name, callback in [("colorsChanged", self.persist_graph_style),
                               ("layoutChanged", self.persist_graph_style),
                               ("legendChanged", self.update_graph_legend),
                               ("selectionChanged", self.graph_selection_changed),
                               ("reportRequested", self.report_graph_selection),
                               ("proximityRequested", self.report_isolate_proximity)]:
            if hasattr(view, name):
                getattr(view, name).connect(
                    lambda *args, callback=callback, role=role: callback(*args, role=role))
        if hasattr(view, 'nodeActivated'):
            view.nodeActivated.connect(self.inspect_graph_sample)

    def _graph_views(self):
        return {'current': self.tree, 'baseline': self.baseline_tree,
                'counterpart': self.counterpart_tree}

    def _graph_view(self, role):
        return self._graph_views().get(role, self.tree)

    def _view_active(self, role):
        """Whether a pane is actually on screen; a hidden tree is not a second opinion."""
        if role == 'baseline':
            return self.dual_toggle.isChecked()
        if role == 'counterpart':
            return self.counterpart_toggle.isChecked()
        return True

    def _reset_baseline_state(self):
        self._baseline_snapshot = None
        self._baseline_key = None
        self._baseline_drawn = False
        self._baseline_diff = None
        self._baseline_draw_requested = False

    # --- the two typing views ------------------------------------------------
    def _build_guidance_row(self):
        """A published cutoff for this exact reference, offered and never applied."""
        row = QWidget()
        layout = QHBoxLayout(row)
        layout.setContentsMargins(0, 0, 0, 0)
        self.guidance_banner = label('', 'small', True)
        layout.addWidget(self.guidance_banner, 1)
        self.guidance_button = button('Review this cutoff…', self.show_threshold_guidance)
        layout.addWidget(self.guidance_button)
        self.guidance_row = row
        row.hide()
        return row

    def set_advanced_visible(self, shown):
        """Open or close the rarely-used graph controls; nothing else changes."""
        shown = bool(shown)
        self.advanced_panel.setVisible(shown)
        self.advanced_toggle.setText('Advanced ▾' if shown else 'Advanced ▸')
        if hasattr(self, 'project'):
            self.project.set_setting('compare.advanced_open', shown)

    def typing_view_title(self, kind=None):
        return typing_title(kind or self.typing_kind)

    def counterpart_kind(self):
        return COUNTERPART.get(self.typing_kind, 'cgmlst')

    def _kind_setting(self, name, kind=None):
        return f'compare.{name}.{kind or self.typing_kind}'

    def _store_kind_state(self):
        """Keep this view's threshold, overlap, reference and investigation to itself."""
        if not hasattr(self, 'cluster_threshold'):
            return
        kind = self.typing_kind
        self.project.set_setting(self._kind_setting('threshold', kind), self.cluster_threshold.value())
        self.project.set_setting(self._kind_setting('overlap', kind), self.overlap.value())
        self.project.set_setting(self._kind_setting('scheme', kind), self.compare_scheme.currentData() or '')
        self.project.set_setting(self._kind_setting('investigation', kind), self.active_investigation_id)

    def _load_kind_state(self):
        """Restore this view's own settings. Defaults are never inherited across kinds."""
        for control, name, default in ((self.cluster_threshold, 'threshold', 1),
                                       (self.overlap, 'overlap', 0.95)):
            value = self.project.get_setting(self._kind_setting(name), default)
            control.blockSignals(True)
            try:
                control.setValue(type(default)(value))
            except (TypeError, ValueError):
                control.setValue(default)
            control.blockSignals(False)
        self.threshold_evidence = {}
        self._refresh_threshold_units()

    def _restore_kind_scheme(self):
        """Re-select this view's own reference once the selector has been rebuilt."""
        stored = self.project.get_setting(self._kind_setting('scheme'), '')
        index = self.compare_scheme.findData(stored) if stored else -1
        if index >= 0:
            self.compare_scheme.blockSignals(True)
            self.compare_scheme.setCurrentIndex(index)
            self.compare_scheme.blockSignals(False)

    def _remember_threshold(self, value):
        """Each view's threshold is remembered under its own key, never shared."""
        self.project.set_setting(self._kind_setting('threshold'), int(value))

    def _remember_overlap(self, value):
        self.project.set_setting(self._kind_setting('overlap'), float(value))

    def _refresh_threshold_units(self, targets=None):
        """Say what the link threshold counts, so 3 of 7 is never read as 3 of 2 358."""
        scale = (self._current_snapshot or {}).get('scale_caption')
        if targets is None:
            targets = (self._current_snapshot or {}).get('target_loci')
        word = TYPING_SCALES.get(self.typing_kind, TYPING_SCALES['unclassified'])['target_word']
        self.cluster_threshold.setSuffix(f' of {targets} {word}' if targets else ' differences')
        suggestion = self._threshold_suggestion or {}
        self.cluster_threshold.setToolTip('\n'.join(filter(None, [
            f'Single-linkage link threshold for the {self.typing_view_title()} tree'
            + (f' ({scale}).' if scale else '.'),
            'A local protocol parameter, not a universal clinical cutoff.',
            suggestion.get('headline', ''), suggestion.get('reason', ''), SCALE_SEPARATION])))

    def _refresh_typing_menu(self):
        if not hasattr(self, 'typing_menu'):
            return
        counts = getattr(self, '_kind_counts', {}) or {}
        self.typing_menu.setText(f'Typing: {self.typing_view_title()} ▾')
        self.typing_menu.setToolTip('\n'.join([
            'Which comparison this page is showing.',
            *(f"{typing_title(kind)}: {counts.get(kind, 0)} isolate(s) with a stored profile"
              for kind in TYPING_VIEWS),
            SCALE_SEPARATION]))

    def show_typing_view(self, kind, *, chosen=True):
        """Switch the page to the classical MLST tree or to the cgMLST tree.

        The two views share no threshold, no minimum overlap, no reference, no
        investigation and no snapshot. Switching stores what this view had and
        restores what the other view had, so neither can inherit a cutoff that was
        justified for the other's target set.
        """
        if kind not in TYPING_VIEWS:
            raise ValueError(f'Unknown typing view {kind!r}; expected one of {", ".join(TYPING_VIEWS)}.')
        if chosen:
            self._typing_kind_chosen = True
            self.project.set_setting('compare.typing_kind', kind)
        if kind == self.typing_kind:
            self._refresh_typing_menu()
            return
        self._store_kind_state()
        self.typing_kind = self.comparison_mode = kind
        self._comparison_cache = None
        self._current_snapshot = None
        self._counterpart_snapshot, self._counterpart_results = None, []
        self._reset_baseline_state()
        self._load_kind_state()
        self._refresh_typing_menu()
        self.refresh_cohort_table()
        self._restore_kind_scheme()
        self._populate_investigation_combo()
        stored = self.project.get_setting(self._kind_setting('investigation'), None)
        known = {plan['id'] for plan in InvestigationStore(self.project).list()}
        self.select_investigation(stored if stored in known else None)
        if not self.active_investigation_id:
            self.refresh_comparison()

    def set_comparison_mode(self, mode):
        """Menu entry point: 'st'/'mlst' and 'cgmlst' open their own separate trees."""
        if mode == 'snp':
            callback = getattr(self, 'run_ska_selected', None)
            if callback:
                return callback()
            self.notify('SNP comparison requires a reviewed ST/cgMLST cohort and the native SNP backend. No allele distances were relabelled as SNPs.')
            return
        self.show_typing_view('mlst' if mode in {'st', 'mlst'} else 'cgmlst')

    def available_profiles(self, sample):
        return _available_profiles(self.project, sample)

    def available_profile_summaries(self, sample):
        if sample.get('status') != 'completed':
            return []
        rows = self.project.analysis_summaries(sample['id'])
        current = current_input_sha256(sample) or (sample.get('result') or {}).get('input_sha256')
        return [row for row in rows if row['locus_count'] and row['scheme_digest']
                and (not current or not row.get('input_sha256') or row['input_sha256'] == current)]

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

    def comparison_scheme_facts(self, path):
        """Everything the reference selector needs about one installed scheme folder.

        Read once per path per reference revision and cached: name, organism,
        which typing kind the folder is, how many targets it holds, and the
        curated publication key it declares for itself (empty when it declares
        none — an unbound reference gets no published cutoff offered).
        """
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
        entry = {}
        try:
            from wmlstudio.reference_index import scheme_entries
            entry = next(iter(scheme_entries((path,))), {})
        except (ImportError, OSError, ValueError):
            entry = {}          # An unreadable folder has an unknown kind, not a guessed one.
        # A scheme installed into the cgMLST library carries the catalogue pin it
        # was staged for. That pin, not the folder's name, is what binds a
        # published cutoff to this reference.
        slot = {}
        try:
            from wmlstudio.cgmlst_schemes import SLOT_FILENAME
            candidate = path / SLOT_FILENAME
            if candidate.is_file() and candidate.stat().st_size <= 1024 * 1024:
                value = json.loads(candidate.read_text(encoding='utf-8'))
                slot = value if isinstance(value, dict) else {}
        except (ImportError, OSError, UnicodeDecodeError, ValueError):
            slot = {}           # An unreadable pin binds no published cutoff at all.
        facts = {'name': name, 'genus': organism.get('genus', ''), 'species': organism.get('species', ''),
                 'kind': slot.get('kind') or entry.get('kind', 'unknown'),
                 'locus_count': entry.get('locus_count', 0),
                 'catalog_key': str(slot.get('key') or ''),
                 'scheme_key': scheme_key_from(slot, metadata)}
        self._scheme_label_cache[key] = facts
        return facts

    def comparison_scheme_label(self, path):
        """Read bounded metadata once per path/reference revision, never alleles."""
        facts = self.comparison_scheme_facts(path)
        return facts['name'], facts['genus'], facts['species']

    def refresh_cohort_table(self):
        if not hasattr(self, "cohort_table") or self._cohort_filling:
            return
        self._cohort_filling = True
        try:
            samples = self.project.samples()
            summaries = {sample['id']: self.available_profile_summaries(sample) for sample in samples}
            # Which isolates have which kind of profile, read from the stored kind
            # of each analysis rather than from a scheme's name. A project with no
            # profile of the shown kind opens on the kind it actually has, unless
            # a person asked for the empty one on purpose.
            self._digest_kinds = {row['scheme_digest']: _summary_kind(row)
                                  for rows in summaries.values() for row in rows if row.get('scheme_digest')}
            self._kind_counts = {kind: sum(any(_summary_kind(row) == kind for row in rows)
                                           for rows in summaries.values()) for kind in TYPING_VIEWS}
            if (not self._typing_kind_chosen and not self._kind_counts[self.typing_kind]
                    and any(self._kind_counts.values())):
                self.typing_kind = self.comparison_mode = next(
                    kind for kind in TYPING_VIEWS if self._kind_counts[kind])
                self._load_kind_state()
            self._refresh_typing_menu()
            word = TYPING_SCALES.get(self.typing_kind, TYPING_SCALES['unclassified'])['target_word']
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
                for result in summaries[sample['id']]:
                    # Only references of the shown typing kind are offered here.
                    # A cgMLST reference in the MLST selector would put a 2,000
                    # target distance under a 7-locus threshold.
                    if _summary_kind(result) != self.typing_kind:
                        continue
                    fingerprints[result["scheme_digest"]] = (result.get("scheme") or "Unnamed scheme",
                                                             result.get("locus_count") or 0)
            for digest, (name, loci) in sorted(fingerprints.items(), key=lambda p: (p[1][0], p[0])):
                self.compare_scheme.addItem(
                    f"Saved: {name} · {loci} {word} · {digest[:8]}" if loci else f"Saved: {name} · {digest[:8]}",
                    "digest:" + digest)
            if selected_scheme and selected_scheme.startswith('digest:') and self.compare_scheme.findData(selected_scheme) < 0:
                self.compare_scheme.addItem('Pinned reference · ' + selected_scheme[7:15], selected_scheme)
            for path in self.scheme_paths:
                facts = self.comparison_scheme_facts(path)
                name, genus, species = facts['name'], facts['genus'], facts['species']
                if facts['kind'] in TYPING_VIEWS and facts['kind'] != self.typing_kind:
                    continue
                if not self.show_all_comparison_schemes.isChecked() and (
                    (selected_genus and genus and genus != selected_genus) or
                    (self.compare_species.currentData() and species and species != self.compare_species.currentData())
                ):
                    continue
                suffix = ' · organism unknown' if not genus else ''
                if facts['kind'] not in TYPING_VIEWS:
                    suffix += ' · typing kind unknown'
                self.compare_scheme.addItem('Call: ' + name + suffix, str(path))
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
                # Name the typing kind of every stored profile: "not called" for
                # this view must never look like "not called at all".
                stored = sorted({f"{typing_title(_summary_kind(row))}: {row.get('scheme') or 'Unnamed scheme'}"
                                 for row in summaries[sample['id']]})
                values = ["", sample["name"], " ".join([genus, species]).strip() or "Unknown",
                          "; ".join(stored) or "Not called"]
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
            self.project.set_setting(self._kind_setting('scheme'), self.compare_scheme.currentData() or '')
            self.tree_status.setText("Scheme selection changed. Build from saved evidence, or choose Call this scheme for missing profiles.")

    def comparison_results(self):
        return _requested_results(self.project, self._comparison_request())[0]

    def _comparison_request(self):
        plan = InvestigationStore(self.project).get(self.active_investigation_id) if self.active_investigation_id else {}
        return {'chosen': frozenset(self.cohort_ids) if self.cohort_ids is not None else None,
                'choice': self.compare_scheme.currentData(), 'genus': self.compare_genus.currentData(),
                'species': self.compare_species.currentData(), 'overlap': self.overlap.value(),
                'threshold': self.cluster_threshold.value(), 'investigation_id': self.active_investigation_id,
                'investigation_revision': plan.get('updated_at'), 'typing_kind': self.typing_kind,
                'counterpart_kind': self.counterpart_kind() if self.counterpart_toggle.isChecked() else None,
                'counterpart_threshold': self.counterpart_threshold()}

    def counterpart_threshold(self):
        """The other typing view's own stored link threshold; never this view's."""
        try:
            return max(0, int(self.project.get_setting(
                self._kind_setting('threshold', self.counterpart_kind()), 1)))
        except (TypeError, ValueError):
            return 1

    def _comparison_key(self, request, revision):
        return (str(self.project.path), revision, request['chosen'], request['choice'],
                request['genus'], request['species'], request['overlap'], request.get('investigation_id'),
                request.get('typing_kind'), request.get('counterpart_kind'),
                request.get('counterpart_threshold'))

    def _clear_comparison(self):
        self.distance_rows, self._last_comparison = [], []
        self._current_snapshot = None
        if self._pending_graph_state is None and getattr(self.tree, 'nodes', None):
            self._pending_graph_state = self.tree.export_state()
        self.tree.blockSignals(True)
        self.tree.draw_results([], [], self.cluster_threshold.value())
        self.tree.blockSignals(False)
        self._legends['current'] = {}
        self.graph_legend.clear()
        # The drawn baseline stays: it is a stored snapshot, independent of the
        # current cohort. Only the comparison between the two is now unknown.
        self._baseline_diff = None
        self.refresh_changes_view()
        self.refresh_graph_captions()
        self.refresh_profile_matrix([])
        self.refresh_statistics([], [])
        self.refresh_cluster_table()
        self.refresh_threshold_suggestion()

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
        kind = payload.get('typing_kind') or self.typing_kind
        try:
            self.distance_rows = payload['rows']
            snapshot = payload.get('snapshot')
            if snapshot is None or snapshot['threshold'] != self.cluster_threshold.value():
                snapshot = build_snapshot(results, payload['rows'], self.cluster_threshold.value(), self.overlap.value(),
                                          previous=snapshot, reuse=payload.get('reuse'), kind=kind)
                if self.active_investigation_id:
                    plan = InvestigationStore(self.project).get(self.active_investigation_id)
                    snapshot.update(preview=True, investigation_id=plan['id'], investigation_name=plan['name'],
                                    protocol=plan['protocol'], saved_threshold=plan['threshold'],
                                    saved_min_overlap=plan['min_overlap'], threshold_evidence=plan.get('threshold_evidence', {}))
            self._current_snapshot = snapshot
            # The view states its own quantity before it draws anything, so an
            # edge label of "3" always arrives with the target set it counts over.
            self.tree.set_scale(payload.get('scale') or typing_scale(results, kind))
            self.tree.draw_results(results, edges, self.cluster_threshold.value(), groups=snapshot['groups'])
            if self._pending_graph_state is not None and hasattr(self.tree, "restore_state"):
                self.tree.blockSignals(True)
                self.tree.restore_state(self._pending_graph_state)
                self.tree.blockSignals(False)
                self._pending_graph_state = None
                if hasattr(self.tree, 'legend'):
                    self.update_graph_legend(self.tree.legend())
            self._last_comparison = results
            self.refresh_cluster_table()
            self.refresh_profile_matrix(results)
            self.refresh_statistics(results, edges)
            excluded = sum(not row["comparable"] for row in self.distance_rows)
            total = payload['total']
            scale = payload.get('scale') or typing_scale(results, kind)
            name = scale['caption'] if results else f"{scale['title']} · no compatible profiles"
            # With no profiles there is no target set to count over, and saying
            # "of loci" would imply one. State the bare threshold instead.
            units = (f" of {scale['targets']} {scale['target_word']}" if scale['targets']
                     else ' allele differences')
            self.tree_status.setText(f"{name} · {len(results)} / {total} cohort members have a "
                                     f"{scale['title']} profile · {len(edges)} edges · {excluded} pairs excluded. "
                                     'Missing/incompatible profiles are not silently compared.')
            if 'available_count' in payload:
                self.tree_status.setText(f"{name} · {payload['available_count']} / {total} stored "
                                        f"{scale['title']} profiles · link ≤ {self.cluster_threshold.value()}"
                                        f"{units} · shared ≥ {self.overlap.value():.0%} · {excluded} pairs excluded.")
            self.tree_status.setToolTip('Shared-locus fraction = shared callable loci / union of profile loci. '
                                       'Unprofiled, mixed and insufficient-overlap pairs have no assigned distance. '
                                       'Single-link components can chain distant isolates. Thresholds are exploratory unless a pinned local protocol is documented; similarity is not proof of transmission.\n'
                                       + SCALE_SEPARATION)
            if payload.get('reuse'):
                reuse = payload['reuse']
                self.tree_status.setText(self.tree_status.text() + f" Reused {reuse['reused_pairs']} pairs; computed {reuse['computed_pairs']}.")
            if snapshot.get('preview'):
                self.tree_status.setText(self.tree_status.text() + ' Unsaved threshold preview — Edit investigation protocol to retain it.')
            if payload['too_large']:
                self.tree_status.setText('Select at most 500 profiles for this interactive view. Export larger cohorts for offline comparison.')
            if hasattr(self.tree, "available_color_fields"):
                current = self.color_by.currentData()
                self.color_by.blockSignals(True)
                self.color_by.clear()
                self.color_by.addItem("Colour: cluster", "cluster")
                self.color_by.addItem("Colour: ST", "st")
                for field in self.tree.available_color_fields():
                    if field.startswith('metadata:'):
                        self.color_by.addItem(str(field), field)
                self.color_by.setCurrentIndex(max(0, self.color_by.findData(current)))
                self.color_by.blockSignals(False)
            self._apply_counterpart(payload.get('counterpart'))
            self._apply_shared_legend()
            self.refresh_baseline_graph()
            self.refresh_graph_captions()
            self.refresh_threshold_suggestion()
        except Exception as exc:
            self.tree_status.setText(str(exc))

    # --- the other typing view ----------------------------------------------
    def set_counterpart_graph(self, enabled):
        """Show or hide the other typing kind's tree for the same isolates."""
        enabled = bool(enabled)
        self.counterpart_pane.setVisible(enabled)
        self.current_caption.setVisible(enabled or self.dual_toggle.isChecked())
        self._update_export_choice()
        self.project.set_setting('compare.counterpart_graph', enabled)
        if enabled:
            self.graph_split.setSizes([1] * self.graph_split.count())
            self.refresh_comparison()
        else:
            self.clear_counterpart_graph()
            self.refresh_graph_captions()

    def clear_counterpart_graph(self):
        self._counterpart_snapshot, self._counterpart_results, self._counterpart_payload = None, [], None
        self.counterpart_tree.blockSignals(True)
        self.counterpart_tree.draw_results([], [], 1)
        self.counterpart_tree.blockSignals(False)
        self._legends['counterpart'] = {}
        self._render_graph_legend()

    def _apply_counterpart(self, payload):
        """Draw the other kind's tree from its own results, threshold and legend.

        Nothing is shared with the tree beside it: its own reference, its own
        target count and its own link threshold, because the two numbers are not
        the same measurement and must never be read off one scale.
        """
        if not self.counterpart_toggle.isChecked():
            return
        kind = self.counterpart_kind()
        self._counterpart_payload = payload
        if not payload or payload.get('too_large') or not payload.get('results'):
            self._counterpart_snapshot, self._counterpart_results = None, []
            self.counterpart_tree.blockSignals(True)
            self.counterpart_tree.set_scale(typing_scale([], kind))
            self.counterpart_tree.draw_results([], [], 1)
            self.counterpart_tree.blockSignals(False)
            self._legends['counterpart'] = {}
            self.counterpart_tree.hide()
            self.counterpart_notice.setText(
                f'This cohort has no stored {typing_title(kind)} profiles, so there is no second tree to '
                f'draw. A missing {typing_title(kind)} profile is unknown evidence, not a distance of zero.'
                if not (payload or {}).get('too_large') else
                'Too many profiles for an interactive second tree; export the larger cohort instead.')
            self.counterpart_notice.show()
            self.counterpart_caption.setText(f'{typing_title(kind)} · no stored profiles')
            self.counterpart_caption.setToolTip(self.counterpart_notice.text())
            return
        self.counterpart_notice.hide()
        self.counterpart_tree.show()
        threshold = self.counterpart_threshold()
        snapshot = payload['snapshot']
        if snapshot.get('threshold') != threshold:
            snapshot = build_snapshot(payload['results'], payload['rows'], threshold,
                                      payload['snapshot'].get('min_overlap', self.overlap.value()), kind=kind)
        self._counterpart_snapshot = snapshot
        self._counterpart_results = payload['results']
        self.counterpart_tree.blockSignals(True)
        self.counterpart_tree.set_scale(payload.get('scale') or typing_scale(payload['results'], kind))
        self.counterpart_tree.draw_results(payload['results'], payload['edges'], threshold,
                                           groups=snapshot['groups'])
        self.counterpart_tree.blockSignals(False)
        self._legends['counterpart'] = self.counterpart_tree.legend()
        self._render_graph_legend()

    def _update_export_choice(self):
        """Offer an export target only for the trees that are actually on screen."""
        if not hasattr(self, 'export_tree_choice'):
            return
        active = self.dual_toggle.isChecked() or self.counterpart_toggle.isChecked()
        self.export_tree_choice.setVisible(active)
        for index in range(self.export_tree_choice.count()):
            role = self.export_tree_choice.itemData(index)
            self.export_tree_choice.model().item(index).setEnabled(self._view_active(role))
        if not self._view_active(self.export_tree_choice.currentData() or 'current'):
            self.export_tree_choice.setCurrentIndex(0)

    def type_comparison_scheme(self):
        if self.busy():
            return
        choice = self.compare_scheme.currentData()
        if choice and choice.startswith('digest:'):
            plan = InvestigationStore(self.project).get(self.active_investigation_id) if self.active_investigation_id else {}
            choice = plan.get('scheme_path')
            if not choice:
                candidates = [r.get('scheme_path') for s in self.project.samples() for r in self.available_profiles(s)
                              if r.get('scheme_digest') == self.compare_scheme.currentData().removeprefix('digest:')]
                choice = next((path for path in candidates if path), None)
        if not choice:
            self.notify("Choose a 'Call:' reference to call missing profiles. Saved profile-only investigations can be compared without sequence inputs.")
            return
        chosen = self.cohort_ids if self.cohort_ids is not None else {s["id"] for s in self.project.samples()}
        samples = [s for s in self.project.samples() if s["id"] in chosen and s.get("input_path") and not s.get("missing_input")]
        if not samples:
            self.error("This cohort has no available sequence inputs. Import a profile table or reconnect the assemblies.")
            return
        if QMessageBox.question(self, "Call comparison scheme", f"Call {Path(choice).name} for {len(samples)} input samples?\n\nPrimary MLST/ST results are retained. New profiles are stored separately and reused on subsequent comparisons. FASTQ files cannot be typed without assembly.") != QMessageBox.StandardButton.Yes:
            return
        saved_profiles = {s['id']: self.available_profiles(s) for s in samples}
        samples = [dict(s, metadata={**s.get("metadata", {}), "workflow": {
            **s.get('metadata', {}).get('workflow', {}), "typing_mode": "manual", "scheme_path": choice}}) for s in samples]
        self.worker_role = "secondary"
        from wmlstudio.scheduler import resources_for_run
        allocation = resources_for_run(self.project.get_setting('analysis_plan', {}))
        self.worker = ComparisonCallingWorker(samples, choice, saved_profiles, parent=self, resource_plan=allocation)

        def prepared(plan):
            self.project.set_setting('last_comparison_digest', plan['scheme_digest'])
            self.notify(f"Reuse {len(plan['reused'])} profiles; call {len(plan['pending'])} assemblies; "
                        f"{len(plan['reads'])} read inputs still need assembly.")

        self.worker.prepared.connect(prepared)

        def saved(sid, result):
            if result.get("alleles"):
                result["scheme_path"] = choice
                # The typing kind is classified from the profile that was
                # produced, never stamped from whichever view happened to be
                # open: a 2,000-target result filed as classical MLST would put
                # a core-genome distance under a seven-locus threshold.
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
        self.refresh_cohort_table()
        # A newly called reference decides which tree it belongs in. If that is
        # not the tree on screen, the page follows the evidence and says so,
        # rather than showing the new profiles under the other kind's threshold.
        kind = self._digest_kinds.get(digest)
        if kind in TYPING_VIEWS and kind != self.typing_kind:
            self.notify(f'That reference produces {typing_title(kind)} profiles, so the '
                        f'{typing_title(kind)} tree is now shown. The {self.typing_view_title()} tree '
                        'keeps its own threshold and reference.')
            self.show_typing_view(kind)
        index = self.compare_scheme.findData("digest:" + digest)
        if index >= 0:
            self.compare_scheme.setCurrentIndex(index)
        self.refresh_comparison()

    def refresh_profile_matrix(self, results):
        loci = sorted({locus for result in results for locus in result.get("alleles", {})})
        if hasattr(self, 'graph_tabs'):
            index = self.graph_tabs.indexOf(self.profile_table)
            if index >= 0:
                self.graph_tabs.setTabToolTip(index, f'{self.typing_view_title()} profiles only · '
                                                     f'{len(loci)} loci in this reference.\n' + SCALE_SEPARATION)
        # "_sample_id" is not a column; it is how a right-click on a sorted matrix
        # resolves the isolate under the cursor instead of trusting the row number.
        self.profile_model.replace(["Sample", "ST", "Scheme", *loci],
                                  [{"_sample_id": r["sample_id"], "Sample": r["sample_name"], "ST": r.get("st"),
                                    "Scheme": r.get("scheme"), **r.get("alleles", {})} for r in results])

    def refresh_statistics(self, results, edges):
        comparable = [p for p in self.distance_rows if p["comparable"]]
        distances = sorted(p["distance"] for p in comparable)
        scale = typing_scale(results, (self._current_snapshot or {}).get('typing_kind') or self.typing_kind)
        total_loci = scale['targets']
        parts = [f"<h2>{html.escape(scale['title'])} cohort evidence</h2>",
                 f"<p>{html.escape(scale['caption'])}</p>",
                 f"<p class='notice'>{html.escape(SCALE_SEPARATION)}</p>",
                 f"<p>{len(results)} profiles · {total_loci} {scale['target_word']} · {len(comparable)} comparable pairs · {len(self.distance_rows) - len(comparable)} excluded pairs</p>"]
        if distances:
            median = (distances[(len(distances) - 1) // 2] + distances[len(distances) // 2]) / 2
            parts.append(f"<p><b>Allele differences over {total_loci} {scale['target_word']}:</b> minimum "
                         f"{distances[0]}, median {median:g}, maximum {distances[-1]}.</p>")
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
        # Two trees with different label or halo settings are much harder to read
        # against each other, so a presentation toggle applies to both.
        for role, view in self._graph_views().items():
            if hasattr(view, name):
                getattr(view, name)(value)
                self.persist_graph_style(role=role)
        self.refresh_graph_captions()

    def set_graph_color_by(self):
        field = self.color_by.currentData() or "cluster"
        for role, view in self._graph_views().items():
            if hasattr(view, "set_color_by") and field in view.available_color_fields():
                view.set_color_by(field)
                self.persist_graph_style(role=role)
        self._apply_shared_legend()

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

    def persist_graph_style(self, *_args, role='current'):
        # Separate keys: the baseline tree's own layout must never overwrite the
        # arrangement the user made of the current tree.
        if role == 'baseline':
            if self._pending_baseline_graph_state is None and hasattr(self.baseline_tree, "export_state"):
                self.project.set_setting("graph_style.baseline", self.baseline_tree.export_state())
            return
        if role == 'counterpart':
            if hasattr(self.counterpart_tree, 'export_state'):
                self.project.set_setting('graph_style.counterpart', self.counterpart_tree.export_state())
            return
        if self._pending_graph_state is None and hasattr(self.tree, "export_state"):
            self.project.set_setting("graph_style", self.tree.export_state())

    def update_graph_legend(self, values, role='current'):
        self._legends[role] = dict(values)
        self._render_graph_legend()

    @staticmethod
    def _legend_markup(values):
        entries = []
        for name, value in list(values.items())[:8]:
            color = QColor(str(value))
            swatch = color.name() if color.isValid() else '#7b8496'
            entries.append(f'<span style="color:{swatch}">●</span> {html.escape(str(name))}')
        if len(values) > 8:
            entries.append(f'+ {len(values) - 8} categories · all groups in the Groups tab')
        return ' &nbsp; '.join(entries)

    def _render_graph_legend(self):
        """One legend row for the current and baseline trees, which share a scale.

        The other typing view keeps its own key beneath its own tree: two kinds
        can each hold a "Cluster 001", and they are not the same isolates.
        """
        values = {**self._legends.get('baseline', {}), **self._legends.get('current', {})}
        self.graph_legend.setText(self._legend_markup(values))
        self.graph_legend.setToolTip('\n'.join(str(key) for key in values))
        counterpart = self._legends.get('counterpart', {})
        self.counterpart_legend.setText(self._legend_markup(counterpart))
        self.counterpart_legend.setToolTip(
            f'{typing_title(self.counterpart_kind())} groups only.\n' + SCALE_SEPARATION)
        self.counterpart_legend.setVisible(bool(counterpart))

    # --- right-click on the comparison views --------------------------------
    def install_graph_menus(self):
        """Give every comparison view the shared right-click actions.

        An action whose handler this window does not implement is left out of the
        menu rather than offered and then failing, so the views can be wired
        before every handler exists.
        """
        if not callable(getattr(self, 'install_view_menu', None)):
            return
        from wmlstudio.context_menus import install_context_menu
        self.install_view_menu('compare.cohort', self.cohort_table)
        self.install_view_menu('compare.profiles', self.profile_table)
        groups = _ClusterTableAdapter('compare.groups', self.cluster_table, self)
        install_context_menu(self.cluster_table, groups, self.show_context_menu)
        getattr(self, '_context_adapters', {})['compare.groups'] = groups
        for view in self._graph_views().values():
            view.context_extension = self.graph_context_entries

    def graph_context_entries(self, view, menu, node_key):
        """Append the shared isolate actions under the graph's own presentation entries."""
        from wmlstudio.context_menus import SEPARATOR, Selection, action_state, submenu_entries
        ids = tuple(view.selected_ids())
        members = tuple(view._members.get(node_key, ())) if node_key else ()
        selection = Selection('compare.graph', ids, members[0] if len(members) == 1 else None,
                              group_id=node_key, on_blank=node_key is None and not ids)
        entries = self.context_menu_plan(selection)
        if not entries:
            return {}
        menu.addSeparator()
        busy = bool(self.busy()) if callable(getattr(self, 'busy', None)) else False
        samples = self.context_samples(selection)
        dispatch = {}
        for entry in entries:
            if entry is SEPARATOR:
                menu.addSeparator()
                continue
            enabled, reason = action_state(entry, selection, window=self, busy=busy, samples=samples)
            if entry.submenu:
                options = submenu_entries(entry.submenu, self, selection)
                if not options:
                    continue
                sub = menu.addMenu(entry.format_title(selection))
                sub.setEnabled(enabled)
                for title, value in options:
                    dispatch[sub.addAction(title)] = (
                        lambda handler=entry.handler, value=value: getattr(self, handler)(selection, value))
                continue
            action = menu.addAction(entry.format_title(selection))
            action.setEnabled(enabled)
            if reason:
                action.setToolTip(reason)
            dispatch[action] = lambda handler=entry.handler: getattr(self, handler)(selection)
        return dispatch

    # context_add_to_comparison / context_remove_from_comparison are implemented on
    # WorkbenchMixin (ui_workbench.update_comparison_cohort), which also records
    # where the cohort change came from. One implementation, one wording.

    def context_select_group(self, selection):
        """Select the isolates of the group that was right-clicked, not the group id."""
        snapshot = self._current_snapshot or {}
        members = {sid for group in snapshot.get('groups', [])
                   if group['id'] == selection.group_id for sid in group['members']}
        if not members and selection.sample_ids:
            self.tree.select_cluster(selection.sample_ids[0])
            return
        self.tree.select_ids(members)

    # --- the baseline ("first/original") tree -------------------------------
    def set_dual_graph(self, enabled):
        """Show or hide the baseline tree; hidden it costs no layout and no scene."""
        enabled = bool(enabled)
        self.baseline_pane.setVisible(enabled)
        self.current_caption.setVisible(enabled or self.counterpart_toggle.isChecked())
        self._update_export_choice()
        self.project.set_setting('compare.dual_graph', enabled)
        if enabled:
            self.graph_split.setSizes([1] * self.graph_split.count())
            self.refresh_baseline_graph()
        else:
            self.clear_baseline_graph()
        self.refresh_graph_captions()

    def clear_baseline_graph(self):
        """Free the second scene but keep the pointer: the baseline itself is stored."""
        self._baseline_drawn = False
        self._baseline_draw_requested = False
        self.baseline_tree.blockSignals(True)
        self.baseline_tree.draw_results([], [], self.cluster_threshold.value())
        self.baseline_tree.blockSignals(False)
        self._legends['baseline'] = {}
        self._render_graph_legend()

    def baseline_reference(self):
        """(investigation_id, snapshot_id) for the pinned baseline, else (None, None)."""
        investigation_id = self.active_investigation_id
        if not investigation_id:
            return None, None
        try:
            return investigation_id, InvestigationStore(self.project).baseline_snapshot_id(investigation_id)
        except KeyError:
            return None, None

    def refresh_baseline_graph(self):
        """Load, diff and draw the baseline. Never raises into the current tree's path."""
        if not self.dual_toggle.isChecked():
            return
        try:
            investigation_id, snapshot_id = self.baseline_reference()
            if not snapshot_id:
                self._baseline_snapshot = self._baseline_diff = self._baseline_key = None
                self._show_baseline_empty_state()
                self.refresh_changes_view()
                return
            key = (str(self.project.path), investigation_id, snapshot_id)
            if key != self._baseline_key:
                self._baseline_snapshot = InvestigationStore(self.project).snapshot(investigation_id, snapshot_id)
                self._baseline_key, self._baseline_drawn = key, False
                self._baseline_draw_requested = False
            self._baseline_diff = (snapshot_diff(self._baseline_snapshot, self._current_snapshot)
                                   if self._baseline_snapshot and self._current_snapshot else None)
            self.refresh_changes_view()
            if not self._baseline_drawn:
                count = len(self._baseline_snapshot.get('profiles', []))
                if count > BASELINE_AUTODRAW_LIMIT and not self._baseline_draw_requested:
                    self._show_baseline_draw_button(count)
                    self.refresh_graph_captions()
                    return
                self._draw_baseline()
            self._apply_shared_legend()
            self.refresh_graph_captions()
        except Exception as error:
            self._baseline_drawn = False
            self._show_baseline_empty_state('The baseline tree could not be opened: ' + str(error))

    def _draw_baseline(self):
        results, edges, groups = snapshot_graph(self._baseline_snapshot)
        QApplication.setOverrideCursor(Qt.CursorShape.WaitCursor)
        try:
            self.baseline_tree.blockSignals(True)
            self.baseline_tree.draw_results(results, edges, self._baseline_snapshot.get('threshold', 1),
                                            groups=groups)
            if self._pending_baseline_graph_state is not None and hasattr(self.baseline_tree, 'restore_state'):
                try:
                    self.baseline_tree.restore_state(self._pending_baseline_graph_state)
                except ValueError:
                    pass  # A stored arrangement from another version is presentation only.
                self._pending_baseline_graph_state = None
            self.baseline_tree.blockSignals(False)
        finally:
            QApplication.restoreOverrideCursor()
        self._baseline_drawn = True
        self._legends['baseline'] = self.baseline_tree.legend()
        self._render_graph_legend()
        self.baseline_notice.hide()
        self.baseline_action.hide()
        self.baseline_tree.show()

    def _show_baseline_empty_state(self, message=''):
        self.baseline_tree.hide()
        self.baseline_notice.setText(message or
            'No baseline tree yet. A baseline is a frozen copy of an earlier comparison: save this '
            'comparison as an investigation and its first build is pinned automatically.')
        self.baseline_notice.show()
        self.baseline_action.setText('Save this comparison as an investigation…')
        self.baseline_action.setProperty('baselineAction', 'save')
        self.baseline_action.show()
        self.baseline_caption.setText('Baseline · none pinned')
        self.baseline_caption.setToolTip(self.baseline_notice.text())

    def _show_baseline_draw_button(self, count):
        self.baseline_tree.hide()
        self.baseline_notice.setText(
            f'This baseline holds {count} isolates. Arranging that many nodes takes several seconds, '
            'so it is drawn only when you ask for it.')
        self.baseline_notice.show()
        self.baseline_action.setText(f'Draw the baseline tree ({count} isolates; this can take several seconds)')
        self.baseline_action.setProperty('baselineAction', 'draw')
        self.baseline_action.show()

    def _baseline_action_clicked(self):
        if self.baseline_action.property('baselineAction') == 'draw':
            self._baseline_draw_requested = True
            self.refresh_baseline_graph()
            return
        self.create_investigation_dialog()

    @staticmethod
    def _scale_words(snapshot):
        """'cgMLST · 2358 targets' for a snapshot, or a plain 'not recorded'."""
        kind = (snapshot or {}).get('typing_kind') or ''
        targets = (snapshot or {}).get('target_loci')
        word = TYPING_SCALES.get(kind, TYPING_SCALES['unclassified'])['target_word']
        return ' · '.join(filter(None, [typing_title(kind) if kind else '',
                                        f'{targets} {word}' if targets else '']))

    def refresh_graph_captions(self):
        """Say which quantity, which cohort and which moment each tree is."""
        current = self._current_snapshot or {}
        threshold = self.cluster_threshold.value()
        words = self._scale_words(current)
        self.current_caption.setText(
            ' · '.join(filter(None, ['Current', words, f"{len(current.get('profiles', []))} isolates",
                                     f'link ≤ {threshold}']))
            if current else f'Current · {self.typing_view_title()} · no comparison built yet')
        self.current_caption.setToolTip(
            f"The {self.typing_view_title()} comparison as it is now · reference "
            f"{current.get('scheme') or 'not recorded'} "
            f"({str(current.get('scheme_digest') or 'no fingerprint')[:12]}) · link ≤ {threshold} "
            f"allele differences · shared ≥ {self.overlap.value():.0%} of loci.\n"
            'This is a layout of allele differences, not a phylogeny and not a transmission chain.\n'
            + SCALE_SEPARATION)
        self._refresh_counterpart_caption()
        if not self._baseline_snapshot:
            return
        snapshot = self._baseline_snapshot
        when = str(snapshot.get('created_at') or '')[:16].replace('T', ' ')
        summary = snapshot_diff_caption(self._baseline_diff) if self._baseline_diff else ''
        self.baseline_caption.setText(' · '.join(filter(None, [
            'Baseline', when or 'date not recorded', self._scale_words(snapshot),
            f"{len(snapshot.get('profiles', []))} isolates", f"link ≤ {snapshot.get('threshold')}", summary])))
        tooltip = [f"Baseline snapshot {str(snapshot.get('snapshot_id') or '')[:12]} of "
                   f"{snapshot.get('investigation_name') or 'this investigation'}, frozen "
                   f"{snapshot.get('created_at') or 'at an unrecorded time'}.",
                   f"Reference {snapshot.get('scheme') or 'not recorded'} "
                   f"({str(snapshot.get('scheme_digest') or 'no fingerprint')[:12]}) · "
                   f"link ≤ {snapshot.get('threshold')} · shared ≥ "
                   f"{(snapshot.get('min_overlap') or 0):.0%}.",
                   'Stored evidence: nothing here is recalculated, and nothing here can be edited.']
        if self.tree.merge_identical:
            tooltip.append('Identical-genotype merging is unavailable for a stored snapshot, which keeps '
                           'per-locus call evidence out of the frozen record.')
        if self._baseline_diff and not self._baseline_diff['policy']['comparable']:
            tooltip.append(self._baseline_diff['policy']['reason'])
        self.baseline_caption.setToolTip('\n'.join(tooltip))

    def _refresh_counterpart_caption(self):
        if not self.counterpart_toggle.isChecked():
            return
        kind = self.counterpart_kind()
        snapshot = self._counterpart_snapshot
        if not snapshot:
            return
        self.counterpart_caption.setText(' · '.join(filter(None, [
            self._scale_words(snapshot) or typing_title(kind),
            f"{len(snapshot.get('profiles', []))} isolates",
            f"link ≤ {snapshot.get('threshold')}"])))
        self.counterpart_caption.setToolTip('\n'.join([
            f"The same isolates, typed the other way: {typing_title(kind)} · reference "
            f"{snapshot.get('scheme') or 'not recorded'} "
            f"({str(snapshot.get('scheme_digest') or 'no fingerprint')[:12]}) · "
            f"link ≤ {snapshot.get('threshold')} over {snapshot.get('target_loci')} targets.",
            'This tree keeps its own threshold and its own legend. Its distances are never '
            'comparable with the tree beside it.', SCALE_SEPARATION]))

    def refresh_changes_view(self):
        if not hasattr(self, 'changes_view'):
            return
        if self._baseline_diff is None:
            reason = ('Build a comparison to see what changed since the baseline.'
                      if self._baseline_snapshot else
                      'No baseline tree is pinned. Save this comparison as an investigation; its first '
                      'build becomes the baseline, and later builds are compared against it.')
            self.changes_view.setHtml('<h2>What changed since the baseline</h2><p>' + html.escape(reason)
                                      + '</p>')
            return
        self.changes_view.setHtml(snapshot_diff_html(self._baseline_diff))

    def show_changes_tab(self):
        index = self.graph_tabs.indexOf(self.changes_view)
        if index >= 0:
            self.graph_tabs.setCurrentIndex(index)

    def _apply_shared_legend(self):
        """Pin one colour per category across both trees, so the legend is shared."""
        views = [view for role, view in self._graph_views().items()
                 if view is not None and self._view_active(role) and getattr(view, '_results', None)]
        if len(views) < 2 or (self.color_by.currentData() or 'cluster') == 'cluster':
            # Cluster colour already derives from the lineage-stable group number,
            # so the two trees agree by construction and pinning would only lie.
            # Across typing kinds it would lie outright: the groups are different
            # isolates measured over different loci.
            for view in views:
                view.set_pinned_legend({})
            return
        categories = sorted({category for view in views for category in view.color_categories()})
        palette = views[0]._palette
        pinned = {category: palette[index % len(palette)] for index, category in enumerate(categories)}
        for view in views:
            view.set_pinned_legend(pinned)

    def choose_baseline_snapshot(self):
        if not self.active_investigation_id:
            self.notify('Save this comparison as an investigation first. Its first build becomes the baseline.')
            return
        store = InvestigationStore(self.project)
        catalog = store.snapshot_catalog(self.active_investigation_id)
        if not catalog:
            self.notify('Build the saved investigation once to create its first snapshot.')
            return
        labels = [f"#{entry['number']} · {str(entry['created_at'] or 'date not recorded')[:16].replace('T', ' ')}"
                  f" · {entry['cohort_size'] or 'unknown'} isolates · link ≤ {entry['threshold']}"
                  + (' · current baseline' if entry['is_baseline'] else '') for entry in catalog]
        current = next((index for index, entry in enumerate(catalog) if entry['is_baseline']), 0)
        chosen, accepted = QInputDialog.getItem(self, 'Choose baseline snapshot',
            'The left tree shows this frozen snapshot. Stored snapshots are never modified:', labels, current, False)
        if not accepted:
            return
        entry = catalog[labels.index(chosen)]
        store.set_baseline(self.active_investigation_id, entry['snapshot_id'])
        self._reset_baseline_state()
        self.refresh_baseline_graph()

    def pin_current_as_baseline(self):
        if not self.active_investigation_id or not self._current_snapshot:
            self.notify('Save this comparison as an investigation and build it once before pinning a baseline.')
            return
        snapshot_id = self._current_snapshot.get('snapshot_id')
        try:
            InvestigationStore(self.project).set_baseline(self.active_investigation_id, snapshot_id)
        except KeyError:
            self.notify('This comparison is an unsaved preview, so it is not a stored snapshot yet. '
                        'Build the saved investigation to store it, then pin it.')
            return
        self._reset_baseline_state()
        self.refresh_baseline_graph()
        self.notify('Baseline pinned. The stored snapshot is unchanged; only which one is shown changed.')

    def inspect_graph_sample(self, sid):
        callback = getattr(self, 'open_isolate_record', None)
        if callback:
            return callback(sid)
        self.selection_ids = {sid}
        self.library_filter = ("ids", [sid])
        self.refresh_tables()
        for row in range(self.sample_table.rowCount()):
            if self.sample_table.item(row, 0).data(Qt.ItemDataRole.UserRole) == sid:
                self.sample_table.setCurrentCell(row, 0)
        self.navigate(1)

    def export_role(self):
        """Which tree the export entries act on; 'current' unless the user chose otherwise."""
        role = self.export_tree_choice.currentData() or 'current'
        return role if self._view_active(role) else 'current'

    def export_snapshot(self, role=None):
        role = role or self.export_role()
        return {'baseline': self._baseline_snapshot,
                'counterpart': self._counterpart_snapshot}.get(role, self._current_snapshot)

    def export_graph_action(self, index, combo):
        combo.setCurrentIndex(0)
        if index == 0:
            return
        role = self.export_role()
        if index in {10, 11, 12, 13}:
            return self.export_baseline_comparison(index)
        if role == 'baseline' and not self._baseline_drawn:
            self.notify('Show the baseline tree first: there is nothing drawn to export.')
            return
        if role == 'counterpart' and not self._counterpart_results:
            self.notify(f'There is no stored {typing_title(self.counterpart_kind())} tree to export.')
            return
        if index == 1 and role == 'current':
            return self.save_tree()
        if index == 5 and role == 'current':
            return self.export_distances()
        formats = {1: ('png', 'save_image'), 2: ("svg", "save_svg"), 3: ("graphml", "save_graphml"),
                   4: ("nwk", "save_newick"), 5: ('json', 'distances'), 6: ("tsv", None),
                   7: ('jpg', 'save_image'), 8: ('tsv', 'profiles'), 9: ('tsv', 'groups')}
        extension, method = formats[index]
        # The typing kind is in the suggested file name as well as in the file,
        # so an MLST export and a cgMLST export cannot be confused on disk.
        kind = self.counterpart_kind() if role == 'counterpart' else self.typing_kind
        try:
            snapshot = self.export_snapshot(role)
            path, _ = QFileDialog.getSaveFileName(self, "Export comparison", f"{kind}-comparison.{extension}", f"{extension.upper()} (*.{extension})")
            if not path:
                return
            self.check_output(path)
            if method in {'profiles', 'groups'}:
                self.write_comparison_table(path, method, snapshot=snapshot if role != 'current' else None)
            elif method == 'distances':
                from wmlstudio.export import write_distances
                write_distances(snapshot.get('pairs', []), path, snapshot.get('min_overlap', 0.95))
            elif method:
                view = self._graph_view(role)
                title = {'baseline': 'Baseline allele-distance minimum spanning forest',
                         'counterpart': f"{typing_title(kind)} allele-distance minimum spanning forest"}.get(role)
                if method in {'save_graphml', 'save_newick'}:
                    getattr(view, method)(path)
                else:
                    getattr(view, method)(path, title=title, subtitle=self._export_subtitle(role))
            else:
                self.write_distance_matrix(path, snapshot=snapshot if role != 'current' else None)
            self.notify({'baseline': 'Baseline snapshot exported exactly as it was stored.',
                         'counterpart': f'{typing_title(kind)} tree exported with its own reference and '
                                        'target count.'}.get(
                             role, 'Comparison exported with its current cohort and evidence.'))
        except Exception as exc:
            self.error(exc)

    def _export_subtitle(self, role):
        snapshot = self.export_snapshot(role)
        if role == 'baseline' and snapshot:
            return (f"{self._scale_words(snapshot)} · baseline snapshot frozen "
                    f"{snapshot.get('created_at') or 'at an unrecorded time'} · "
                    f"link ≤ {snapshot.get('threshold')} · not a phylogeny or transmission tree")
        if role == 'counterpart' and snapshot:
            return (f"{self._scale_words(snapshot)} · link ≤ {snapshot.get('threshold')} · "
                    'a separate measurement from the other tree; not a phylogeny or transmission tree')
        return None

    def export_baseline_comparison(self, index):
        """The two trees together, and the change summary that explains them."""
        baseline, current = self._baseline_snapshot, self._current_snapshot
        if not baseline or not current or self._baseline_diff is None:
            missing = 'baseline tree' if not baseline else 'current comparison'
            self.notify(f'This export needs both trees; there is no {missing} yet.')
            return
        if index in {10, 11} and not self._baseline_drawn:
            self.notify('Show the baseline tree first: there is nothing drawn to export.')
            return
        extension = {10: 'png', 11: 'jpg', 12: 'tsv', 13: 'json'}[index]
        path, _ = QFileDialog.getSaveFileName(self, 'Export baseline comparison',
                                              f'baseline-comparison.{extension}',
                                              f'{extension.upper()} (*.{extension})')
        if not path:
            return
        try:
            self.check_output(path)
            if index in {10, 11}:
                image = self.render_dual_graph()
                if not image.save(str(path), 'JPEG' if index == 11 else 'PNG'):
                    raise OSError(f'Could not save image: {path}')
            elif index == 12:
                self.write_change_summary(path)
            else:
                self.write_change_summary_json(path)
            self.notify('Exported both trees with the policy that produced each of them.')
        except Exception as exc:
            self.error(exc)

    def render_dual_graph(self, width=1800, height=1200):
        """One image holding both trees, headlined with whether they are comparable."""
        policy = self._baseline_diff['policy']
        headline = policy['reason']
        if policy['identical']:
            headline = ('Same reference, same link threshold, same minimum shared loci — the two '
                        'trees can be read against each other.')
        baseline, current = self._baseline_snapshot, self._current_snapshot
        return render_side_by_side(
            self.baseline_tree, self.tree,
            left_title=f"Baseline · {str(baseline.get('created_at') or '')[:16].replace('T', ' ')} · "
                       f"{self._scale_words(baseline)} · "
                       f"{len(baseline.get('profiles', []))} isolates · link ≤ {baseline.get('threshold')}",
            right_title=f"Current · {self._scale_words(current)} · "
                        f"{len(current.get('profiles', []))} isolates · "
                        f"link ≤ {current.get('threshold')}",
            left_subtitle=self._export_subtitle('baseline'),
            right_subtitle=None, headline=headline, width=width, height=height)

    def write_change_summary(self, path):
        import csv

        from wmlstudio.export import _atomic_text, _csv_cell
        rows = snapshot_diff_rows(self._baseline_diff)
        with _atomic_text(path) as handle:
            writer = csv.writer(handle, delimiter='\t')
            writer.writerow(DIFF_ROW_FIELDS)
            for row in rows:
                writer.writerow([_csv_cell(row[key]) for key in DIFF_ROW_FIELDS])
            writer.writerow([])
            for caveat in self._baseline_diff['caveats']:
                writer.writerow([_csv_cell('how_to_read'), _csv_cell(caveat)])

    def write_change_summary_json(self, path):
        from copy import deepcopy
        from datetime import datetime, timezone

        from wmlstudio import __version__
        from wmlstudio.export import _atomic_text
        document = {'application': f'WMLSTudio {__version__}',
                    'exported_at': datetime.now(timezone.utc).isoformat(),
                    **deepcopy(self._baseline_diff)}
        with _atomic_text(path) as handle:
            json.dump(document, handle, indent=2, ensure_ascii=False, allow_nan=False)

    def write_comparison_table(self, path, table, snapshot=None):
        import csv

        from wmlstudio.export import _atomic_text, _csv_cell
        if snapshot is not None and table == 'profiles':
            # known_alleles omits every uncallable locus, so a blank cell here is
            # "not called", never a matching allele.
            profiles = snapshot.get('profiles', [])
            loci = sorted({locus for p in profiles for locus in (p.get('known_alleles') or {})})
            headers = ['Sample', 'ST', 'Scheme', *loci]
            rows = [{'Sample': p.get('sample_name'), 'ST': p.get('st'), 'Scheme': snapshot.get('scheme'),
                     **(p.get('known_alleles') or {})} for p in profiles]
        elif snapshot is not None:
            rows, headers = snapshot.get('groups', []), None
        else:
            rows = self.profile_model.rows if table == 'profiles' else (self._current_snapshot or {}).get('groups', [])
            headers = self.profile_model.headers if table == 'profiles' else None
        if headers is None:
            headers = ['id', 'name', 'status', 'members', 'change', 'parents', 'max_direct_distance', 'chained', 'unassessed_within_pairs']
        with _atomic_text(path) as handle:
            writer = csv.writer(handle, delimiter='\t')
            writer.writerow(headers)
            for row in rows:
                writer.writerow([_csv_cell(row.get(key)) for key in headers])

    def write_distance_matrix(self, path, snapshot=None):
        import csv

        from wmlstudio.export import _atomic_text, _csv_cell
        if snapshot is not None:
            results = [{'sample_id': str(p['sample_id']), 'sample_name': p.get('sample_name', p['sample_id'])}
                       for p in snapshot.get('profiles', [])]
            pairs = snapshot.get('pairs', [])
        else:
            results, pairs = self._last_comparison, self.distance_rows
        distances = {frozenset((p["source"], p["target"])): p["distance"] if p["comparable"] else "NA" for p in pairs}
        with _atomic_text(path) as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["sample", *[_csv_cell(r["sample_name"]) for r in results]])
            for result in results:
                sid = result["sample_id"]
                writer.writerow([_csv_cell(result["sample_name"]), *[0 if other["sample_id"] == sid
                                else distances.get(frozenset((sid, other["sample_id"])), "NA") for other in results]])

    def report_comparison(self):
        self.report_ids = {r["sample_id"] for r in self._last_comparison}
        self._report_investigation_snapshot = self._current_snapshot
        self.report_scope.setCurrentIndex(0)
        self.refresh_report_table()
        self.navigate(5)
