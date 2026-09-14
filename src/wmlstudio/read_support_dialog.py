"""Native one-isolate read-support review, background analysis and drill-down."""

from __future__ import annotations

import copy
import threading
from pathlib import Path

from PySide6.QtCore import Qt, QThread, Signal
from PySide6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QDialog,
    QDialogButtonBox,
    QFileDialog,
    QFormLayout,
    QHBoxLayout,
    QInputDialog,
    QSpinBox,
    QVBoxLayout,
)

from .identification import cached_scheme
from .read_support import investigate_read_support, persist_read_support, preflight_read_support
from .sample_workflow import _sha256, current_input_sha256
from .scheduler import plan_resources
from .sequence import AnalysisCancelled
from .ui_common import cell, make_table
from .widgets import button, label


class _SchemeWorker(QThread):
    loaded = Signal(object)
    failed = Signal(str)

    def __init__(self, path, parent=None):
        super().__init__(parent)
        self.path, self.cancel_event = path, threading.Event()

    def run(self):
        try:
            self.loaded.emit(cached_scheme(self.path, cancelled=self.cancel_event.is_set))
        except AnalysisCancelled:
            pass
        except Exception as error:
            self.failed.emit(str(error))


class ReadSupportPlanDialog(QDialog):
    """Explicit locus selection; native reference loading never blocks the GUI."""

    def __init__(self, sample, scheme_paths, initial_scheme=None, profiles=None, parent=None):
        super().__init__(parent)
        self.sample, self.profiles = copy.deepcopy(sample), copy.deepcopy(profiles or [])
        if sample.get('result'):
            self.profiles.append(copy.deepcopy(sample['result']))
        self.plan, self.scheme, self.worker, self._closing = None, None, None, False
        self.setWindowTitle('Investigate unresolved loci using attached reads')
        self.resize(900, 660)
        layout = QVBoxLayout(self)
        layout.addWidget(label(f'Read support · {sample["name"]}', 'title', True))
        layout.addWidget(label('Separate investigative evidence. This does not change assembly alleles, STs or cgMLST distances. '
                               'Inspect matched bases, coverage and ambiguous alternatives; no unsupported locus absence is inferred.', 'muted', True))
        reads = (sample.get('metadata') or {}).get('reads', {}).get('reads', [])
        self.reads = sorted(reads, key=lambda row: row.get('mate', 0))
        self.read_identity_ok = (len(self.reads) == 2 and [row.get('mate') for row in self.reads] == [1, 2]
                                 and all(_sha256(row.get('sha256')) and row.get('path') for row in self.reads))
        layout.addWidget(label('Attached originals: ' + (' · '.join(Path(row['path']).name for row in self.reads)
                                                       if self.read_identity_ok else 'No verified paired FASTQs attached. Attach a pair first.'), 'small', True))
        row = QHBoxLayout()
        self.scheme_choice = QComboBox()
        self.scheme_choice.setMinimumContentsLength(20)
        paths = list(dict.fromkeys(str(Path(path).resolve()) for path in scheme_paths))
        if initial_scheme and str(Path(initial_scheme).resolve()) not in paths:
            paths.insert(0, str(Path(initial_scheme).resolve()))
        for path in paths:
            self.scheme_choice.addItem(Path(path).name, path)
        if initial_scheme:
            self.scheme_choice.setCurrentIndex(max(0, self.scheme_choice.findData(str(Path(initial_scheme).resolve()))))
        row.addWidget(label('Reference scheme', 'small'))
        row.addWidget(self.scheme_choice, 1)
        self.browse_button = button('Browse…', self.browse_scheme)
        row.addWidget(self.browse_button)
        layout.addLayout(row)
        self.loci = make_table(['Investigate', 'Locus', 'Current assembly call', 'Call state', 'Alleles to test'])
        self.loci.setSortingEnabled(False)
        self.loci.setMinimumHeight(170)
        layout.addWidget(self.loci, 1)
        actions = QHBoxLayout()
        actions.addWidget(button('Select unresolved', self.select_unresolved))
        actions.addWidget(button('Clear selection', self.clear_selection))
        actions.addStretch()
        layout.addLayout(actions)
        form = QFormLayout()
        self.max_pairs = QSpinBox()
        self.max_pairs.setRange(1, 1_000_000)
        self.max_pairs.setSingleStep(10_000)
        self.max_pairs.setValue(100_000)
        form.addRow('Maximum read pairs (first N, not random)', self.max_pairs)
        layout.addLayout(form)
        layout.addWidget(label('Limits: 20 loci, 5,000 alleles, 10 million reference bases, 60 million read bases and 128 MiB alignment output. '
                               'Exceeding a bound fails explicitly; no alleles are silently discarded. Default screening uses Phred+33 ≥20, '
                               '≥90% read identity and ≥80% read alignment; 95% locus breadth at ≥3-read depth is a support heuristic, not an allele call.', 'small', True))
        self.feedback = label('Choose a scheme and explicit loci to review the plan.', 'muted', True)
        layout.addWidget(self.feedback)
        self.buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel)
        self.run_button = self.buttons.button(QDialogButtonBox.StandardButton.Ok)
        self.run_button.setText('Run reviewed read investigation')
        self.run_button.setEnabled(False)
        self.buttons.accepted.connect(self.accept)
        self.buttons.rejected.connect(self.reject)
        layout.addWidget(self.buttons)
        self.loci.itemChanged.connect(self.update_plan)
        self.max_pairs.valueChanged.connect(self.update_plan)
        self.scheme_choice.currentIndexChanged.connect(self.load_scheme)
        if self.scheme_choice.count():
            self.load_scheme()

    def browse_scheme(self):
        path = QFileDialog.getExistingDirectory(self, 'Choose a local typing scheme')
        if path:
            index = self.scheme_choice.findData(str(Path(path).resolve()))
            if index < 0:
                self.scheme_choice.addItem(Path(path).name, str(Path(path).resolve()))
                index = self.scheme_choice.count() - 1
            self.scheme_choice.setCurrentIndex(index)

    def load_scheme(self, *_):
        path = self.scheme_choice.currentData()
        if not path or self.worker is not None and self.worker.isRunning():
            return
        self.scheme = None
        self.loci.setRowCount(0)
        self.run_button.setEnabled(False)
        self.scheme_choice.setEnabled(False)
        self.browse_button.setEnabled(False)
        self.feedback.setText('Loading and validating allele references in the background…')
        self.worker = _SchemeWorker(path, self)
        self.worker.loaded.connect(self.scheme_loaded)
        self.worker.failed.connect(self.feedback.setText)
        self.worker.finished.connect(self.load_finished)
        self.worker.start()

    def load_finished(self):
        self.scheme_choice.setEnabled(True)
        self.browse_button.setEnabled(True)
        if self._closing:
            super().reject()

    def scheme_loaded(self, scheme):
        if self._closing:
            return
        self.scheme = scheme
        current_sha = current_input_sha256(self.sample)
        matched = next((row for row in self.profiles if row.get('scheme_digest') == scheme.digest
                        and current_sha and row.get('input_sha256') == current_sha), {})
        calls = {row['locus']: row for row in matched.get('calls', [])}
        self.loci.blockSignals(True)
        self.loci.setRowCount(len(scheme.loci))
        for index, locus in enumerate(scheme.loci):
            call = calls.get(locus, {})
            item = cell('', locus)
            item.setFlags(item.flags() | Qt.ItemFlag.ItemIsUserCheckable)
            item.setCheckState(Qt.CheckState.Unchecked)
            self.loci.setItem(index, 0, item)
            self.loci.setItem(index, 1, cell(locus))
            self.loci.setItem(index, 2, cell((matched.get('alleles') or {}).get(locus) or '—'))
            self.loci.setItem(index, 3, cell(call.get('status') or 'No current profile under this scheme'))
            self.loci.setItem(index, 4, cell(str(len(scheme.alleles[locus]))))
        self.loci.blockSignals(False)
        self.update_plan()

    def select_unresolved(self):
        rows = [index for index in range(self.loci.rowCount())
                if self.loci.item(index, 3).text() not in {'exact', 'novel'}
                and not self.loci.item(index, 3).text().startswith('No current')]
        if not rows:
            self.feedback.setText('No current unresolved assembly calls are available under this scheme. Select specific loci manually for a clearly labelled investigation or positive control.')
            return
        if len(rows) > 20:
            self.feedback.setText(f'{len(rows)} unresolved loci exceed the 20-locus bound. Select a smaller explicit set; none was silently selected or omitted.')
            return
        self.loci.blockSignals(True)
        for index in rows:
            self.loci.item(index, 0).setCheckState(Qt.CheckState.Checked)
        self.loci.blockSignals(False)
        self.update_plan()

    def clear_selection(self):
        self.loci.blockSignals(True)
        for index in range(self.loci.rowCount()):
            self.loci.item(index, 0).setCheckState(Qt.CheckState.Unchecked)
        self.loci.blockSignals(False)
        self.update_plan()

    def update_plan(self, *_):
        self.plan = None
        self.run_button.setEnabled(False)
        if self.scheme is None or not self.read_identity_ok:
            return
        loci = [self.loci.item(index, 0).data(Qt.ItemDataRole.UserRole) for index in range(self.loci.rowCount())
                if self.loci.item(index, 0).checkState() == Qt.CheckState.Checked]
        try:
            panel = preflight_read_support(self.scheme, loci)
            resources = plan_resources(threads_per_sample=2, memory_gb=2)
        except ValueError as error:
            self.feedback.setText(str(error))
            return
        self.plan = {'sample_id': self.sample['id'], 'scheme': self.scheme, 'scheme_path': str(self.scheme.path),
                     'loci': loci, 'max_pairs': self.max_pairs.value(), 'threads': resources.threads_per_sample,
                     'assembly_path': self.sample['input_path'], 'read_paths': [row['path'] for row in self.reads],
                     'expected_read_sha256': [row['sha256'] for row in self.reads], 'panel': panel,
                     'memory_reservation_gb': resources.memory_gb}
        self.feedback.setText(f'Reviewed plan: {len(loci)} loci · all {panel["allele_count"]:,} allele references · '
                              f'first ≤{self.max_pairs.value():,} read pairs · {resources.threads_per_sample} threads. '
                              'Full input hashes are checked. No typing result will be changed.')
        self.run_button.setEnabled(True)

    def accept(self):
        if self.worker is not None and self.worker.isRunning():
            return
        self.update_plan()
        if self.plan is not None:
            super().accept()

    def reject(self):
        if self.worker is not None and self.worker.isRunning():
            self._closing = True
            self.worker.cancel_event.set()
            self.feedback.setText('Cancelling reference loading safely…')
            return
        super().reject()

    def closeEvent(self, event):
        if self.worker is not None and self.worker.isRunning():
            event.ignore()
            self.reject()
        else:
            super().closeEvent(event)


class ReadSupportResultDialog(QDialog):
    def __init__(self, sample_name, result, parent=None):
        super().__init__(parent)
        self.result_data = result
        self.setWindowTitle(f'Read-support evidence · {sample_name}')
        self.resize(1080, 680)
        layout = QVBoxLayout(self)
        sampling = result['sampling']
        layout.addWidget(label(f'{sampling["reads_aligned"]:,} reads aligned · '
                               + ('sampled first-pair prefix' if sampling['sampled'] else 'complete paired files'), 'title', True))
        layout.addWidget(label('Investigative support, not an allele call. Ambiguous and incomplete findings remain saved. '
                               'No-support does not establish absence; candidate depth includes non-unique reads and overlapping mates.', 'muted', True))
        controls = QHBoxLayout()
        self.locus_choice = QComboBox()
        for locus in result['loci']:
            self.locus_choice.addItem(f'{locus["locus"]} · {locus["status"]}', locus)
        controls.addWidget(self.locus_choice, 1)
        self.include_zero = QCheckBox('Show no-support candidates (all are saved)')
        controls.addWidget(self.include_zero)
        layout.addLayout(controls)
        self.table = make_table(['Candidate allele', 'Support state', 'Breadth at depth threshold', 'Mean read depth',
                                 'Observed identity %', 'Unique-best reads', 'Tied-best reads', 'Mixed positions', 'Indel positions'])
        layout.addWidget(self.table, 1)
        self.summary = label('', 'small', True)
        layout.addWidget(self.summary)
        layout.addWidget(label('Scheme SHA-256: ' + result['scheme_digest'] + '\nAssembly SHA-256: ' + str(result.get('input_sha256')), 'small', True))
        buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Close)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)
        self.locus_choice.currentIndexChanged.connect(self.refresh_candidates)
        self.include_zero.toggled.connect(self.refresh_candidates)
        self.refresh_candidates()

    def refresh_candidates(self, *_):
        locus = self.locus_choice.currentData() or {}
        all_rows = locus.get('candidates', [])
        rows = [row for row in all_rows if self.include_zero.isChecked() or row['status'] != 'no_support']
        self.table.setSortingEnabled(False)
        self.table.setRowCount(len(rows))
        for index, row in enumerate(rows):
            values = [row['allele'], row['status'], f'{100 * row["breadth_min_depth"]:.1f}%', f'{row["mean_depth"]:.2f}',
                      f'{row["observed_identity_pct"]:.2f}' if row['observed_identity_pct'] is not None else '—',
                      row['uniquely_best_reads'], row['tied_best_reads'], row['mixed_positions'],
                      row['deletion_positions'] + row['insertion_anchors']]
            for column, value in enumerate(values):
                item = cell(str(value))
                item.setToolTip('Review positions: ' + str(row.get('review_positions') or 'none'))
                self.table.setItem(index, column, item)
        self.table.setSortingEnabled(True)
        self.summary.setText(f'{len(rows):,} of {len(all_rows):,} tested candidate alleles displayed. '
                             f'Coverage threshold: ≥{self.result_data["parameters"]["min_depth"]} reads at Phred+33 ≥'
                             f'{self.result_data["parameters"]["min_base_quality_phred33"]}. No allele was assigned.')


def launch_read_support(window, sample_id=None):
    """Native entry point; uses the application's cancellable background task."""
    if window.busy():
        return False
    eligible = [sample for sample in window.project.samples() if sample.get('input_path')
                and len((sample.get('metadata') or {}).get('reads', {}).get('reads', [])) == 2]
    if sample_id is None:
        if not eligible:
            window.error('Attach a validated original FASTQ pair to an assembly before investigating read support.')
            return False
        names = [f'{sample["name"]} · {sample["id"][:8]}' for sample in eligible]
        name, accepted = QInputDialog.getItem(window, 'Choose one isolate', 'Assembly with attached paired reads', names, 0, False)
        if not accepted:
            return False
        sample_id = eligible[names.index(name)]['id']
    sample = next((sample for sample in eligible if sample['id'] == sample_id), None)
    if sample is None:
        window.error('The selected isolate has no validated attached FASTQ pair.')
        return False
    workflow = (sample.get('metadata') or {}).get('workflow') or {}
    initial = workflow.get('scheme_path') or (sample.get('result') or {}).get('scheme_path') or window.project.get_setting('scheme_path', '')
    profiles = window.available_profiles(sample) if hasattr(window, 'available_profiles') else []
    review = ReadSupportPlanDialog(sample, window.scheme_paths, initial_scheme=initial, profiles=profiles, parent=window)
    if review.exec() != QDialog.DialogCode.Accepted:
        return False
    plan = review.plan

    def operation(cancelled, progress):
        return investigate_read_support(*plan['read_paths'], plan['scheme'], plan['loci'],
            cancelled=cancelled, progress=progress, assembly_path=plan['assembly_path'],
            expected_read_sha256=plan['expected_read_sha256'], max_pairs=plan['max_pairs'], threads=plan['threads'])

    def finish(result):
        try:
            persist_read_support(window.project, sample_id, result)
        except (ValueError, OSError) as error:
            window.error(str(error))
            return
        window.notify('Read-support investigation saved separately; assembly alleles and ST are unchanged.')
        ReadSupportResultDialog(sample['name'], result, window).exec()

    return window.launch_task(operation, 'read_support', finish)
