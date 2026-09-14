import csv

import pytest
from PySide6.QtWidgets import QDialog

from wmlstudio.app import MainWindow
from wmlstudio.cohort_picker import CohortPickerDialog


@pytest.fixture
def window(qtbot, tmp_path):
    window = MainWindow(storage_root=tmp_path / 'app')
    qtbot.addWidget(window)
    window.show()
    qtbot.wait(30)
    yield window
    window.close()


def sample(window, name, gene):
    digest = 'a' * 64
    return window.project.add_profile(name, {'sample_name': name, 'st': '20', 'input_sha256': digest,
        'scheme': 'MLST', 'scheme_digest': 'reference', 'status': 'profile_imported', 'alleles': {'a': '1'}},
        {'organism': {'genus': 'Staphylococcus', 'species': 'aureus'},
         'annotations': {'ward': 'ICU'}, 'hydra': {'evidence_input_sha256': digest,
         'source_sample': name, 'hits': [{'gene': gene, 'element_type': 'AMR', 'primary': True}]}})


def test_evidence_cohort_is_independent_and_exports_only_visible_columns(window, tmp_path, monkeypatch):
    first, second = sample(window, 'A', 'blaA'), sample(window, 'B', 'vanB')
    window.selection_ids = {second}
    window.report_ids = {second}
    window.refresh_features()
    assert window.feature_model.rows == [] and window.amr_model.rows == []
    def choose(dialog):
        assert dialog.selected_ids == set()
        dialog.selected_ids = {first}
        return QDialog.DialogCode.Accepted
    monkeypatch.setattr(CohortPickerDialog, 'exec', choose)
    window.choose_feature_cohort()
    assert window.feature_ids == {first}
    assert window.selection_ids == window.report_ids == {second}
    assert [row['_sample_id'] for row in window.feature_model.rows] == [first]
    assert 'blaA' in window.amr_model.headers and 'vanB' not in window.amr_model.headers
    opened = []
    monkeypatch.setattr(window, 'open_isolate_record', opened.append)
    window.feature_table.doubleClicked.emit(window.feature_model.index(0, 0))
    assert opened == [first]
    path = tmp_path / 'evidence.tsv'
    window.write_feature_table(path, 'amr')
    with path.open(encoding='utf-8-sig') as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    assert len(rows) == 1 and rows[0]['sample_id'] == first
    assert second not in path.read_text() and 'vanB' not in path.read_text()
    window.gene_filter.setText('not-present')
    window.write_feature_table(path, 'amr')
    assert 'blaA' not in path.read_text()


def test_advanced_source_controls_are_not_in_main_cohort_view(window):
    tabs = window.evidence_tabs
    assert any('Advanced' in tabs.tabText(index) for index in range(tabs.count()))
    assert not window.hydra_table.isVisible()
    window.navigate(4)
    assert not window.hydra_table.isVisible()
    assert window.feature_scope_label.isVisible()
