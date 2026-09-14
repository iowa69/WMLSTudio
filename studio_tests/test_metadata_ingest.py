from copy import deepcopy

import pytest
from PySide6.QtCore import Qt
from PySide6.QtWidgets import QDialog

from wmlstudio.metadata_ingest import (
    CLEAR_TOKEN,
    annotation_values,
    apply_metadata_preview,
    export_metadata,
    preview_metadata,
    read_metadata_table,
)
from wmlstudio.project import Project


@pytest.fixture
def project(tmp_path):
    with Project(tmp_path / 'metadata.wmlstudio') as project:
        yield project


def isolate(project, name, metadata=None):
    return project.add_profile(name, {'sample_name': name, 'scheme': 'test', 'scheme_digest': 'digest',
                                    'st': '20', 'alleles': {'a': '1'}, 'status': 'profile_imported'}, metadata)


def test_identity_mapping_is_explicit_and_ambiguous_names_block(project):
    one, two = isolate(project, 'duplicate'), isolate(project, 'duplicate')
    unique = isolate(project, 'Unique')
    preview = preview_metadata(project, [{'sample_name': 'Unique', 'ward': 'ICU'}])
    assert preview['errors']
    preview = preview_metadata(project, [{'sample_name': 'duplicate', 'ward': 'ICU'}], allow_name_fallback=True)
    assert 'ambiguous' in str(preview['errors'])
    preview = preview_metadata(project, [{'sample_name': 'Unique', 'ward': 'ICU'}], allow_name_fallback=True)
    assert preview['rows'][0]['sample_id'] == unique and not preview['errors']
    assert apply_metadata_preview(project, preview) == [unique]
    supplied_bad_id = preview_metadata(project, [{'sample_id': 'unknown', 'sample_name': 'Unique', 'ward': 'ICU'}], allow_name_fallback=True)
    assert supplied_bad_id['errors'] and not supplied_bad_id['rows'][0]['changes']
    assert annotation_values(project.get_sample(one)) == annotation_values(project.get_sample(two)) == {}


def test_duplicate_rows_and_outside_cohort_do_not_apply_partially(project):
    one, two = isolate(project, 'A'), isolate(project, 'B')
    preview = preview_metadata(project, [{'sample_id': one, 'ward': 'ICU'}, {'sample_id': one, 'ward': 'ER'}])
    with pytest.raises(ValueError, match='Resolve every'):
        apply_metadata_preview(project, preview)
    assert annotation_values(project.get_sample(one)) == {}
    restricted = preview_metadata(project, [{'sample_id': two, 'ward': 'ICU'}], sample_ids={one})
    assert 'outside' in str(restricted['errors'])


@pytest.mark.parametrize('bad_date', ['2025-02-29', '12/09/2026', '2026-9-12', '2026-13-01', '2026-09-12T12:00:00'])
def test_dates_are_typed_unambiguous_iso_calendar_values(project, bad_date):
    sid = isolate(project, 'A')
    preview = preview_metadata(project, [{'sample_id': sid, 'collection_date': bad_date}])
    assert preview['errors']
    with pytest.raises(ValueError):
        apply_metadata_preview(project, preview)


def test_blank_preserves_and_clear_requires_explicit_review(project):
    sid = isolate(project, 'A', {'annotations': {'ward': 'ICU', 'collection_date': '2024-02-29'}, 'hydra': {'keep': True}})
    before = project.get_sample(sid)
    blank = preview_metadata(project, [{'sample_id': sid, 'ward': '', 'collection_date': ' '}])
    assert apply_metadata_preview(project, blank) == []
    clearing = preview_metadata(project, [{'sample_id': sid, 'ward': CLEAR_TOKEN}], clear_token=CLEAR_TOKEN)
    with pytest.raises(ValueError, match='approve overwrites'):
        apply_metadata_preview(project, clearing)
    apply_metadata_preview(project, clearing, accept_conflicts=True)
    after = project.get_sample(sid)
    assert after['metadata']['annotations']['ward'] is None
    assert after['metadata']['annotations']['collection_date'] == '2024-02-29'
    assert after['metadata']['hydra'] == before['metadata']['hydra']
    assert after['result'] == before['result'] and after['input_path'] == before['input_path']


def test_protected_columns_and_preview_tampering_cannot_overwrite_evidence(project):
    sid = isolate(project, 'A')
    for field in ['ST', 'input_path', 'result', 'workflow', 'hydra', 'organism', 'characterization']:
        preview = preview_metadata(project, [{'sample_id': sid, field: 'replacement'}])
        assert preview['errors']
        with pytest.raises(ValueError):
            apply_metadata_preview(project, preview, accept_conflicts=True)
    preview = preview_metadata(project, [{'sample_id': sid, 'ward': 'ICU'}])
    preview['rows'][0]['changes'] = {'input_path': '/different/input'}
    with pytest.raises(ValueError, match='Protected'):
        apply_metadata_preview(project, preview)
    assert project.get_sample(sid)['result']['st'] == '20'


def test_stale_preview_and_second_write_failure_roll_back_whole_batch(project, monkeypatch):
    first, second = isolate(project, 'A'), isolate(project, 'B')
    rows = [{'sample_id': first, 'ward': 'ICU'}, {'sample_id': second, 'ward': 'ER'}]
    preview = preview_metadata(project, rows)
    project.set_metadata(second, {'annotations': {'specimen': 'blood'}})
    with pytest.raises(ValueError, match='changed after'):
        apply_metadata_preview(project, preview)
    assert annotation_values(project.get_sample(first)) == {}
    preview = preview_metadata(project, rows)
    original = project.set_metadata
    def fail_second(sid, metadata):
        if sid == second:
            raise OSError('Injected storage failure')
        return original(sid, metadata)
    history_before = deepcopy(project.history())
    monkeypatch.setattr(project, 'set_metadata', fail_second)
    with pytest.raises(OSError):
        apply_metadata_preview(project, preview)
    assert annotation_values(project.get_sample(first)) == {}
    assert annotation_values(project.get_sample(second)) == {'specimen': 'blood'}
    assert project.history() == history_before


def test_csv_tsv_read_export_and_protected_destination(project, tmp_path):
    sid = isolate(project, '=A', {'annotations': {'ward': '=2+3', 'collection_date': '2026-09-12'}})
    other = isolate(project, 'Excluded')
    for suffix in ['csv', 'tsv']:
        path = tmp_path / ('metadata.' + suffix)
        export_metadata(project, path, sample_ids={sid})
        rows = read_metadata_table(path)
        assert len(rows) == 1 and rows[0]['sample_id'] == sid
        assert rows[0]['sample_name'] == "'=A" and rows[0]['ward'] == "'=2+3"
        assert other not in path.read_text()
    with pytest.raises(ValueError, match='different output'):
        export_metadata(project, project.path)
    for text in ['sample_id,ward,ward\nx,ICU,ER\n', 'sample_id,ward\nx,ICU,extra\n', 'sample_id,ward\nx\n']:
        path = tmp_path / 'bad.csv'
        path.write_text(text)
        with pytest.raises(ValueError):
            read_metadata_table(path)


def test_native_preview_requires_conflict_approval_and_maps_unique_name(qtbot, project):
    from wmlstudio.metadata_grid import MetadataPreviewDialog
    sid = isolate(project, 'A', {'annotations': {'ward': 'ICU'}})
    dialog = MetadataPreviewDialog(project, [{'sample_name': 'A', 'ward': 'ER'}])
    qtbot.addWidget(dialog)
    assert not dialog.apply_button.isEnabled()
    dialog.allow_names.setChecked(True)
    assert not dialog.apply_button.isEnabled()
    dialog.approve_conflicts.setChecked(True)
    assert dialog.apply_button.isEnabled()
    dialog.apply_changes()
    assert dialog.result() == QDialog.DialogCode.Accepted
    assert dialog.changed_ids == [sid]
    assert annotation_values(project.get_sample(sid))['ward'] == 'ER'


def test_native_grid_keeps_ids_readonly_and_sorting_does_not_remap(qtbot, tmp_path, monkeypatch):
    from wmlstudio.app import MainWindow
    from wmlstudio.metadata_grid import MetadataGridDialog, MetadataPreviewDialog
    window = MainWindow(storage_root=tmp_path / 'app')
    qtbot.addWidget(window)
    first, second = isolate(window.project, 'B'), isolate(window.project, 'A')
    grid = MetadataGridDialog(window, {first, second})
    qtbot.addWidget(grid)
    assert not grid.table.item(0, 0).flags() & Qt.ItemFlag.ItemIsEditable
    grid.table.sortItems(1, Qt.SortOrder.AscendingOrder)
    column = grid.fields.index('ward') + 2
    grid.table.item(0, column).setText('ER')
    def approve(dialog):
        assert not dialog.preview['errors']
        dialog.apply_changes()
        return dialog.result()
    monkeypatch.setattr(MetadataPreviewDialog, 'exec', approve)
    grid.review_edits()
    assert annotation_values(window.project.get_sample(second))['ward'] == 'ER'
    assert 'ward' not in annotation_values(window.project.get_sample(first))
    assert window.selection_ids == set()
    grid.close()
    window.close()
