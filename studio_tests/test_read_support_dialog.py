import threading

from PySide6.QtCore import Qt
from PySide6.QtWidgets import QDialog

from wmlstudio import read_support_dialog
from wmlstudio.read_support_dialog import ReadSupportPlanDialog, ReadSupportResultDialog
from wmlstudio.sequence import AnalysisCancelled
from wmlstudio.typing import load_scheme


def sample_scheme(tmp_path):
    folder = tmp_path / 'scheme'
    folder.mkdir()
    (folder / 'a.tfa').write_text('>a_1\nACGTACGT\n')
    (folder / 'b.tfa').write_text('>b_1\nATGCATGC\n')
    scheme = load_scheme(folder)
    sample = {'id': 'sample1', 'name': 'One isolate', 'input_path': str(tmp_path / 'assembly.fasta'),
              'metadata': {'reads': {'reads': [{'mate': 1, 'path': 'R1.fastq', 'sha256': '1' * 64},
                                               {'mate': 2, 'path': 'R2.fastq', 'sha256': '2' * 64}]}},
              'result': {'scheme_digest': scheme.digest, 'input_sha256': 'a' * 64,
                         'alleles': {'a': '1', 'b': None}, 'calls': [{'locus': 'a', 'status': 'exact'}, {'locus': 'b', 'status': 'missing'}]}}
    return sample, scheme


def test_dialog_loads_in_background_requires_explicit_loci_and_preserves_calls(qtbot, tmp_path):
    sample, scheme = sample_scheme(tmp_path)
    dialog = ReadSupportPlanDialog(sample, [scheme.path])
    qtbot.addWidget(dialog)
    dialog.show()
    qtbot.waitUntil(lambda: dialog.scheme is not None and not dialog.worker.isRunning())
    assert not dialog.run_button.isEnabled() and dialog.plan is None
    assert dialog.loci.rowCount() == 2
    dialog.select_unresolved()
    assert dialog.plan['loci'] == ['b'] and dialog.plan['panel']['allele_count'] == 1
    assert dialog.plan['expected_read_sha256'] == ['1' * 64, '2' * 64]
    assert sample['result']['alleles'] == {'a': '1', 'b': None}
    dialog.max_pairs.setValue(123)
    dialog.accept()
    assert dialog.result() == QDialog.DialogCode.Accepted
    assert dialog.plan['max_pairs'] == 123 and dialog.plan['assembly_path'] == sample['input_path']


def test_no_current_profile_is_not_silently_selected_as_missing(qtbot, tmp_path):
    sample, scheme = sample_scheme(tmp_path)
    sample['metadata']['profiles'] = []
    sample['result']['scheme_digest'] = 'b' * 64
    dialog = ReadSupportPlanDialog(sample, [scheme.path])
    qtbot.addWidget(dialog)
    qtbot.waitUntil(lambda: dialog.scheme is not None and not dialog.worker.isRunning())
    dialog.select_unresolved()
    assert dialog.plan is None and 'No current unresolved' in dialog.feedback.text()
    dialog.loci.item(0, 0).setCheckState(Qt.CheckState.Checked)
    assert dialog.plan['loci'] == ['a']  # An explicit positive/control investigation remains possible.
    dialog.reject()


def test_dialog_cancel_during_scheme_load_waits_for_thread_without_blocking_ui(qtbot, tmp_path, monkeypatch):
    sample, scheme = sample_scheme(tmp_path)
    started = threading.Event()
    def cancellable_load(path, cancelled=None):
        started.set()
        event = threading.Event()
        while not cancelled():
            event.wait(.005)
        raise AnalysisCancelled()
    monkeypatch.setattr(read_support_dialog, 'cached_scheme', cancellable_load)
    dialog = ReadSupportPlanDialog(sample, [scheme.path])
    qtbot.addWidget(dialog)
    dialog.show()
    qtbot.waitUntil(started.is_set)
    dialog.reject()
    qtbot.waitUntil(lambda: not dialog.worker.isRunning() and not dialog.isVisible())
    assert dialog.plan is None


def test_result_dialog_discloses_zero_candidate_filter_and_sampling(qtbot):
    base = {'allele': '1', 'status': 'supported', 'breadth_min_depth': 1, 'mean_depth': 5,
            'observed_identity_pct': 100, 'uniquely_best_reads': 0, 'tied_best_reads': 10,
            'mixed_positions': 0, 'deletion_positions': 0, 'insertion_anchors': 0}
    result = {'sampling': {'reads_aligned': 100, 'sampled': True},
              'loci': [{'locus': 'a', 'status': 'supported', 'candidates': [base, dict(base, allele='2', status='no_support')]}],
              'scheme_digest': 'a' * 64, 'input_sha256': 'b' * 64,
              'parameters': {'min_depth': 3, 'min_base_quality_phred33': 20}}
    dialog = ReadSupportResultDialog('One isolate', result)
    qtbot.addWidget(dialog)
    assert dialog.table.rowCount() == 1 and '1 of 2 tested' in dialog.summary.text()
    assert 'all are saved' in dialog.include_zero.text()
    dialog.include_zero.setChecked(True)
    assert dialog.table.rowCount() == 2
    assert 'No allele was assigned' in dialog.summary.text()
