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
from wmlstudio.identification import cached_scheme, scheme_organism
from wmlstudio.investigation import (
    InvestigationStore,
    build_snapshot,
    incremental_distances,
    proximity_rows,
)
from wmlstudio.jobs import AnalysisWorker
from wmlstudio.sample_workflow import current_input_sha256, hydra_evidence_status
from wmlstudio.sequence import AnalysisCancelled, SequenceReader, check_cancelled, file_sha256
from wmlstudio.ui_common import FlowLayout, cell, gene_names, make_table, organism_for
from wmlstudio.widgets import TreeView, button, label


def _available_profiles(project, sample):
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
    revision = project.comparison_revision()
    if request.get('revision', revision) != revision:
        raise ValueError('Stored profiles changed before comparison began. Build the comparison again.')
    results, total = _requested_results(project, request, cancelled)
    available_count = len(results)
    if len(results) > 500:
        return {'results': [], 'rows': [], 'edges': [], 'total': total, 'too_large': True}
    store = InvestigationStore(project)
    investigation_id = request.get('investigation_id')
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
            return {'results': [], 'rows': [], 'edges': [], 'total': total, 'too_large': True}
    if investigation_id:
        rows, cache, reuse = incremental_distances(results, request['overlap'], cache, cancelled=cancelled)
    else:
        rows = pairwise_distances(results, request['overlap'], cancelled=cancelled)
        reuse = {'reused_pairs': 0, 'computed_pairs': len(rows)}
    edges = forest_from_distances(rows, cancelled=cancelled)
    snapshot = build_snapshot(results, rows, request.get('threshold', 1), request['overlap'],
                              previous=prior, reuse=reuse, cancelled=cancelled)
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
            'snapshot': snapshot, 'cache': cache, 'reuse': reuse, 'available_count': available_count}


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


class ComparisonWorkspaceMixin:
    def investigation_summary(self):
        """Small, JSON-safe status for the guided overview; never loads a matrix."""
        active = getattr(self, 'active_investigation_id', None)
        if not active:
            return {'id': None, 'name': 'Unsaved investigation', 'saved': False,
                    'cohort_count': len(self.cohort_ids) if self.cohort_ids is not None else len(self.project.samples()),
                    'profiles': len(getattr(self, '_last_comparison', []))}
        try:
            plan = InvestigationStore(self.project).get(active)
        except KeyError:
            return {'id': None, 'name': 'Investigation unavailable', 'saved': False}
        return {'id': active, 'name': plan['name'], 'saved': True, 'cohort_count': len(plan['sample_ids']),
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
        active = self.project.get_setting('investigations.active', None)
        self.active_investigation_id = active if active in {p['id'] for p in store.list()} else None
        if not self.active_investigation_id and self.cohort_ids is None:
            self.cohort_ids = set()
        self.investigation_combo.blockSignals(True)
        self.investigation_combo.clear()
        self.investigation_combo.addItem('Unsaved investigation', None)
        for plan in store.list():
            self.investigation_combo.addItem(plan['name'], plan['id'])
        self.investigation_combo.setCurrentIndex(max(0, self.investigation_combo.findData(self.active_investigation_id)))
        self.investigation_combo.blockSignals(False)
        if self.active_investigation_id:
            self.select_investigation(self.active_investigation_id)

    def _investigation_chosen(self, *_args):
        self.select_investigation(self.investigation_combo.currentData())

    def select_investigation(self, investigation_id):
        self.active_investigation_id = investigation_id
        self.project.set_setting('investigations.active', investigation_id)
        self._comparison_cache = None
        self._current_snapshot = None
        if not investigation_id:
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
        groups = (self._current_snapshot or {}).get('groups', [])
        self._cluster_filling = True
        self.cluster_table.setSortingEnabled(False)
        self.cluster_table.setRowCount(len(groups))
        for row, group in enumerate(groups):
            title = group['name'] + (' · chaining' if group.get('chained') else '')
            for col, value in enumerate((title, len(group['members']), group.get('change', ''))):
                item = cell(value, group['id'])
                item.setToolTip(f"{group['status']} · {group.get('reason', '')}\n"
                                f"Maximum direct distance: {group.get('max_direct_distance')}\n"
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

    def graph_selection_changed(self, ids):
        if not hasattr(self, 'graph_selection_label'):
            return
        self.graph_report_button.setEnabled(bool(ids))
        text = f"{len(ids)} selected"
        if len(ids) == 1 and self._current_snapshot:
            rows = proximity_rows(self._current_snapshot, ids)
            if rows:
                nearest = rows[0]['nearest_distance']
                text += f" · nearest: {nearest if nearest is not None else 'not comparable'} allele differences"
        self.graph_selection_label.setText(text + ' · Ctrl-click adds; click a cluster row selects its members.')

    def report_graph_selection(self, ids=None):
        ids = set(ids if isinstance(ids, (list, tuple, set)) else self.tree.selected_ids())
        if not ids:
            self.notify('Select graph nodes, labels or a cluster row first.')
            return
        self.report_ids = ids
        self.report_scope.setCurrentIndex(0)
        self._report_investigation_snapshot = self._current_snapshot
        self.refresh_report_table()
        self.navigate(5)

    def report_isolate_proximity(self, sid):
        self.tree.select_ids([sid])
        if hasattr(self, 'report_preset'):
            self.report_preset.setCurrentIndex(self.report_preset.findData('proximity'))
        self.report_graph_selection([sid])

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

    def open_investigation_history(self):
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
        for number, sid in reversed(list(enumerate(plan['snapshots'], 1))):
            chooser.addItem(f'Snapshot {number} · {sid[:12]}', sid)
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

    def set_comparison_mode(self, mode):
        if mode == 'snp':
            callback = getattr(self, 'run_ska_selected', None)
            if callback:
                return callback()
            self.notify('SNP comparison requires a reviewed ST/cgMLST cohort and the native SNP backend. No allele distances were relabelled as SNPs.')
            return
        self.comparison_mode = mode
        self.comparison_settings.setWindowTitle('MLST / ST comparison settings' if mode == 'st' else 'cgMLST comparison settings')
        self.comparison_settings.show()
        self.comparison_settings.raise_()

    def apply_threshold_guidance(self, evidence):
        """Called after explicit reference/protocol review in the guidance dialog."""
        from copy import deepcopy
        if evidence.get('approved_threshold') is not None:
            context = evidence.get('context') or {}
            snapshot = self._current_snapshot or {}
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
        self.heading(layout, "Investigate related isolates", "Choose a cohort, review its groups, and report the isolates that matter. Core and accessory evidence stay separate.")
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
        self.cluster_threshold = QSpinBox()
        self.cluster_threshold.setRange(0, 100000)
        self.cluster_threshold.setValue(1)
        self.cluster_threshold.setToolTip('User-defined single-linkage threshold, not a universal clinical cutoff. Pin an organism, scheme and local protocol in a saved investigation.')
        self.cluster_threshold.valueChanged.connect(self.refresh_comparison)
        for text, control in [('Shared ≥', self.overlap), ('Group ≤', self.cluster_threshold)]:
            group = QWidget()
            row = QHBoxLayout(group)
            row.setContentsMargins(0, 0, 0, 0)
            row.addWidget(label(text, 'small'))
            row.addWidget(control)
            controls.addWidget(group)
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
        modes = QToolButton()
        modes.setText('Analysis ▾')
        modes.setPopupMode(QToolButton.ToolButtonPopupMode.InstantPopup)
        mode_menu = QMenu(modes)
        mode_menu.addAction('MLST / ST…', lambda: self.set_comparison_mode('st'))
        mode_menu.addAction('cgMLST…', lambda: self.set_comparison_mode('cgmlst'))
        mode_menu.addAction('SNP follow-up…', lambda: self.set_comparison_mode('snp'))
        modes.setMenu(mode_menu)
        investigations.addWidget(modes)
        investigations.addWidget(button('Compare', self.refresh_comparison))
        layout.addLayout(investigations)
        self.cohort_search = QLineEdit()
        self.cohort_search.setPlaceholderText("Find isolates…")
        self.cohort_search.textChanged.connect(self.refresh_cohort_table)
        left.addWidget(self.cohort_search)
        self.show_all_comparison_schemes = QCheckBox('Show all schemes')
        self.show_all_comparison_schemes.setToolTip('Include schemes for other organisms. Schemes with unknown organism metadata remain available in either mode.')
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
        toolbar = FlowLayout()
        self.color_by = QComboBox()
        self.color_by.setProperty('compactCharacters', 9)
        self.color_by.addItem("Colour: cluster", "cluster")
        self.color_by.addItem("Colour: ST", "st")
        self.color_by.currentIndexChanged.connect(self.set_graph_color_by)
        settings_layout.addWidget(self.color_by)
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
        graph_options.setMenu(menu)
        toolbar.addWidget(graph_options)
        export = QComboBox()
        export.setProperty('compactCharacters', 9)
        export.addItems(["Export…", "PNG image", "SVG vector", "GraphML", "Newick (MST topology)", "Distance JSON", "Distance matrix TSV", "JPEG image", "Profile matrix TSV", "Group table TSV"])
        export.activated.connect(lambda index: self.export_graph_action(index, export))
        toolbar.addWidget(export)
        content.addLayout(toolbar)
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
        tabs.addTab(self.tree, "Graph")
        tabs.addTab(self.cluster_table, 'Groups')
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
        self.comparison_mode = 'unspecified'
        self._comparison_timer = QTimer(self)
        self._comparison_timer.setSingleShot(True)
        self._comparison_timer.timeout.connect(self._start_pending_comparison)
        for name, callback in [("colorsChanged", self.persist_graph_style), ("layoutChanged", self.persist_graph_style),
                               ("legendChanged", self.update_graph_legend), ("nodeActivated", self.inspect_graph_sample),
                               ("selectionChanged", self.graph_selection_changed),
                               ("reportRequested", self.report_graph_selection),
                               ("proximityRequested", self.report_isolate_proximity)]:
            if hasattr(self.tree, name):
                getattr(self.tree, name).connect(callback)
        QTimer.singleShot(0, self.restore_investigations)

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
            summaries = {sample['id']: self.available_profile_summaries(sample) for sample in samples}
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
                    fingerprints[result["scheme_digest"]] = result.get("scheme") or "Unnamed scheme"
            for digest, name in sorted(fingerprints.items(), key=lambda p: (p[1], p[0])):
                self.compare_scheme.addItem(f"Saved: {name} · {digest[:8]}", "digest:" + digest)
            if selected_scheme and selected_scheme.startswith('digest:') and self.compare_scheme.findData(selected_scheme) < 0:
                self.compare_scheme.addItem('Pinned reference · ' + selected_scheme[7:15], selected_scheme)
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
                          "; ".join(sorted({r.get("scheme", "") for r in summaries[sample['id']]})) or "Not called"]
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
        plan = InvestigationStore(self.project).get(self.active_investigation_id) if self.active_investigation_id else {}
        return {'chosen': frozenset(self.cohort_ids) if self.cohort_ids is not None else None,
                'choice': self.compare_scheme.currentData(), 'genus': self.compare_genus.currentData(),
                'species': self.compare_species.currentData(), 'overlap': self.overlap.value(),
                'threshold': self.cluster_threshold.value(), 'investigation_id': self.active_investigation_id,
                'investigation_revision': plan.get('updated_at')}

    def _comparison_key(self, request, revision):
        return (str(self.project.path), revision, request['chosen'], request['choice'],
                request['genus'], request['species'], request['overlap'], request.get('investigation_id'))

    def _clear_comparison(self):
        self.distance_rows, self._last_comparison = [], []
        self._current_snapshot = None
        if self._pending_graph_state is None and getattr(self.tree, 'nodes', None):
            self._pending_graph_state = self.tree.export_state()
        self.tree.blockSignals(True)
        self.tree.draw_results([], [], self.cluster_threshold.value())
        self.tree.blockSignals(False)
        self.graph_legend.clear()
        self.refresh_profile_matrix([])
        self.refresh_statistics([], [])
        self.refresh_cluster_table()

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
            snapshot = payload.get('snapshot')
            if snapshot is None or snapshot['threshold'] != self.cluster_threshold.value():
                snapshot = build_snapshot(results, payload['rows'], self.cluster_threshold.value(), self.overlap.value(),
                                          previous=snapshot, reuse=payload.get('reuse'))
                if self.active_investigation_id:
                    plan = InvestigationStore(self.project).get(self.active_investigation_id)
                    snapshot.update(preview=True, investigation_id=plan['id'], investigation_name=plan['name'],
                                    protocol=plan['protocol'], saved_threshold=plan['threshold'],
                                    saved_min_overlap=plan['min_overlap'], threshold_evidence=plan.get('threshold_evidence', {}))
            self._current_snapshot = snapshot
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
            name = results[0].get("scheme", "") if results else "No compatible profiles"
            self.tree_status.setText(f"{name} · {len(results)} / {total} cohort members have this profile · {len(edges)} edges · {excluded} pairs excluded. Missing/incompatible profiles are not silently compared.")
            if 'available_count' in payload:
                self.tree_status.setText(f"{name} · {payload['available_count']} / {total} stored profiles · "
                                        f"link ≤ {self.cluster_threshold.value()} alleles · shared ≥ {self.overlap.value():.0%} · {excluded} pairs excluded.")
            self.tree_status.setToolTip('Shared-locus fraction = shared callable loci / union of profile loci. '
                                       'Unprofiled, mixed and insufficient-overlap pairs have no assigned distance. '
                                       'Single-link components can chain distant isolates. Thresholds are exploratory unless a pinned local protocol is documented; similarity is not proof of transmission.')
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
        except Exception as exc:
            self.tree_status.setText(str(exc))

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
        analysis_kind = self.comparison_mode
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
                if analysis_kind in {'st', 'cgmlst'}:
                    result['analysis_kind'] = 'mlst' if analysis_kind == 'st' else 'cgmlst'
                    result['analysis_kind_source'] = 'user-selected comparison mode; not independent reference validation'
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
        for name, value in list(values.items())[:8]:
            color = QColor(str(value))
            swatch = color.name() if color.isValid() else '#7b8496'
            entries.append(f'<span style="color:{swatch}">●</span> {html.escape(str(name))}')
        if len(values) > 8:
            entries.append(f'+ {len(values) - 8} categories · all groups in the Groups tab')
        self.graph_legend.setText(' &nbsp; '.join(entries))
        self.graph_legend.setToolTip('\n'.join(str(key) for key in values))

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

    def export_graph_action(self, index, combo):
        combo.setCurrentIndex(0)
        if index == 0:
            return
        if index == 1:
            return self.save_tree()
        if index == 5:
            return self.export_distances()
        formats = {2: ("svg", "save_svg"), 3: ("graphml", "save_graphml"), 4: ("nwk", "save_newick"), 6: ("tsv", None),
                   7: ('jpg', 'save_image'), 8: ('tsv', 'profiles'), 9: ('tsv', 'groups')}
        extension, method = formats[index]
        path, _ = QFileDialog.getSaveFileName(self, "Export comparison", f"comparison.{extension}", f"{extension.upper()} (*.{extension})")
        if not path:
            return
        try:
            self.check_output(path)
            if method in {'profiles', 'groups'}:
                self.write_comparison_table(path, method)
            elif method:
                getattr(self.tree, method)(path)
            else:
                self.write_distance_matrix(path)
            self.notify("Comparison exported with its current cohort and evidence.")
        except Exception as exc:
            self.error(exc)

    def write_comparison_table(self, path, table):
        import csv

        from wmlstudio.export import _atomic_text, _csv_cell
        rows = self.profile_model.rows if table == 'profiles' else (self._current_snapshot or {}).get('groups', [])
        headers = self.profile_model.headers if table == 'profiles' else ['id', 'name', 'status', 'members', 'change', 'parents', 'max_direct_distance', 'chained', 'unassessed_within_pairs']
        with _atomic_text(path) as handle:
            writer = csv.writer(handle, delimiter='\t')
            writer.writerow(headers)
            for row in rows:
                writer.writerow([_csv_cell(row.get(key)) for key in headers])

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
        self._report_investigation_snapshot = self._current_snapshot
        self.report_scope.setCurrentIndex(0)
        self.refresh_report_table()
        self.navigate(5)
