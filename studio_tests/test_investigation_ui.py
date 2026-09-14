"""Outbreak-review workflows: explicit cohorts, immutable evidence, readable exports."""

import json
from copy import deepcopy

import pytest
from PySide6.QtCore import Qt
from PySide6.QtWidgets import QDialog, QMessageBox

from wmlstudio.app import MainWindow
from wmlstudio.export import REPORT_PRESETS, review_report_html, write_json
from wmlstudio.investigation import InvestigationStore, build_snapshot


@pytest.fixture
def window(qtbot, tmp_path, monkeypatch):
    widget = MainWindow(storage_root=tmp_path / 'application')
    errors = []
    monkeypatch.setattr(widget, 'error', lambda message: errors.append(str(message)))
    widget.test_errors = errors
    qtbot.addWidget(widget)
    widget.show()
    qtbot.wait(40)
    yield widget
    widget.cancel_comparison_for_close()
    qtbot.waitUntil(lambda: widget.comparison_worker is None, timeout=5000)
    widget.close()


def add(window, name, vector='1111'):
    return window.project.add_profile(name, {
        'sample_name': name, 'scheme': 'Study cgMLST', 'scheme_digest': 'study-reference',
        'status': 'profile_imported', 'alleles': dict(zip('abcd', vector)), 'calls': [],
        'st': '20', 'input_sha256': 'a' * 64,
    }, {'organism': {'genus': 'Staphylococcus', 'species': 'aureus'}, 'annotations': {'ward': 'ICU'}})


def activate(window, ids):
    window.refresh_cohort_table()
    plan = InvestigationStore(window.project).save('Week 1 <review>', ids, scheme='Study cgMLST',
        scheme_digest='study-reference', threshold=1, include_new=True,
        protocol='Synthetic test protocol; not a clinical cutoff', filters={'genus': 'Staphylococcus', 'species': 'aureus'})
    window.select_investigation(plan['id'])
    assert window._current_snapshot, window.tree_status.text()
    return plan['id']


def test_graph_labels_and_all_groups_select_isolates(window, qtbot):
    a, b, c, d = [add(window, name, vector) for name, vector in [('A', '1111'), ('B', '2111'), ('C', '2211'), ('D', '9999')]]
    activate(window, [a, b, c, d])
    window.navigate(2)
    window.tree.fit_tree()
    qtbot.wait(100)
    assert not window.cohort_table.isVisible()
    assert not window.comparison_settings.isVisible()
    assert len(window.tree._halos) == 2
    group = next(g for g in window._current_snapshot['groups'] if a in g['members'])
    assert group['chained'] and set(group['members']) == {a, b, c}
    position = window.tree.mapFromScene(window.tree.labels[a].sceneBoundingRect().center())
    qtbot.mouseClick(window.tree.viewport(), Qt.MouseButton.LeftButton, pos=position)
    assert set(window.tree.selected_ids()) == {a}
    window.select_graph_cluster()
    assert set(window.tree.selected_ids()) == {a, b, c}
    window.report_graph_selection()
    assert window.report_sample_ids() == {a, b, c}
    assert 'Chaining' in window.report_preview.toPlainText()
    assert window.test_errors == []


def test_weekly_update_reuses_pairs_and_retains_fixed_review_groups(window, monkeypatch):
    a, b = add(window, 'A'), add(window, 'B', '2111')
    iid = activate(window, [a, b])
    store = InvestigationStore(window.project)
    original = store.snapshot(iid)
    frozen = store.freeze_group(iid, 'Original review', [a, b])
    c = add(window, 'C', '2211')
    monkeypatch.setattr(QMessageBox, 'question', lambda *args: QMessageBox.StandardButton.Yes)
    monkeypatch.setattr('wmlstudio.ui_compare.ComparisonCallingWorker.start', lambda *args: pytest.fail('Update must not retype'))
    window.update_investigation()
    latest = store.snapshot(iid)
    assert latest['reuse']['reused_pairs'] == 1 and latest['reuse']['computed_pairs'] == 2
    assert latest['groups'][0]['id'] == original['groups'][0]['id']
    assert latest['groups'][0]['change'] == 'grown'
    assert latest['review_groups'][0]['sample_ids'] == frozen['sample_ids'] == sorted([a, b])
    assert c not in frozen['sample_ids']
    assert store.snapshot(iid, original['snapshot_id']) == original


def test_missing_profiles_are_not_inferred_singletons(window, tmp_path):
    a, b = add(window, 'A'), add(window, 'B', '2111')
    missing_path = tmp_path / 'pending.fasta'
    missing_path.write_text('>one\nACGT\n')
    missing = window.project.add_sample(missing_path)
    window.project.set_metadata(missing, {'organism': {'genus': 'Staphylococcus', 'species': 'aureus'}})
    activate(window, [a, b, missing])
    groups = window._current_snapshot['groups']
    unavailable = next(g for g in groups if missing in g['members'])
    assert unavailable['status'] == 'not_comparable'
    assert window.project.analysis_results(missing) == []
    assert all(p['distance'] is None for p in window.distance_rows if missing in {p['source'], p['target']})


def test_per_isolate_report_scopes_neighbours_and_actual_pdf(window, tmp_path, monkeypatch):
    a, b = add(window, 'FOCAL <A>'), add(window, 'NEIGHBOUR B', '2111')
    activate(window, [a, b])
    window.report_isolate_proximity(a)
    assert window.report_sample_ids() == {a}
    html_path = tmp_path / 'focal.html'
    monkeypatch.setattr('wmlstudio.ui_reports.QFileDialog.getSaveFileName', lambda *args: (str(html_path), ''))
    window.export_report('html')
    html = html_path.read_text()
    assert 'FOCAL &lt;A&gt;' in html and 'Comparison-cohort neighbour' in html
    assert 'NEIGHBOUR B' in html and '4/4' in html
    assert 'data:image/png;base64,' in html
    assert 'Full stored evidence below' not in html
    json_path = tmp_path / 'focal.json'
    write_json(window.report_records(), json_path, selected_ids={a}, investigation=window._current_snapshot)
    document = json.loads(json_path.read_text())
    assert [s['sample_id'] for s in document['samples']] == [a]
    assert document['investigation']['proximity'][0]['nearest'][0]['neighbour_id'] == b
    pdf = tmp_path / 'focal.pdf'
    window.write_pdf_report(pdf, selected_ids={a})
    assert pdf.read_bytes().startswith(b'%PDF') and pdf.stat().st_size > 1000
    assert window.test_errors == []


def test_saved_report_template_and_new_picker_do_not_use_global_selection(window, monkeypatch):
    from wmlstudio.cohort_picker import CohortPickerDialog
    a, b = add(window, 'A'), add(window, 'B')
    window.selection_ids = {b}
    seen = []
    def choose(dialog):
        seen.append(set(dialog.selected_ids))
        dialog.selected_ids = {a}
        return QDialog.DialogCode.Accepted
    monkeypatch.setattr(CohortPickerDialog, 'exec', choose)
    window.choose_report_cohort()
    assert seen == [set()] and window.report_sample_ids() == {a}
    window._report_options_override = {'title': 'Saved IPC review', 'amr': False, 'provenance': False}
    window.save_report_template()
    window._report_options_override = None
    assert window.report_options()['title'] == 'Saved IPC review'
    assert not window.report_options()['amr']
    assert window.selection_ids == {b}


def test_reopen_investigation_clears_old_report_context_and_restores_ids(window, tmp_path):
    a, b = add(window, 'A'), add(window, 'B')
    iid = activate(window, [a, b])
    original_path = window.project.path
    window.report_isolate_proximity(a)
    assert window.switch_project(tmp_path / 'other.wmlstudio')
    assert window._report_investigation_snapshot is None
    assert window.switch_project(original_path)
    assert window.active_investigation_id == iid
    assert window.cohort_ids == {a, b}
    assert window._current_snapshot['scheme_digest'] == 'study-reference'
    assert window.project.get_sample(a)['result']['st'] == '20'


def test_canonical_characterization_state_prevents_stale_labels_and_report(window):
    sid = add(window, 'A')
    sample = window.project.get_sample(sid)
    evidence = {'input_sha256': 'b' * 64, 'species_evidence': {'status': 'completed', 'genus': 'FAKE', 'species': 'STALE'},
                'virulence': {'status': 'detected', 'hits': [{'gene': 'STALE_VIRULENCE'}]},
                'drug_associations': {'status': 'completed', 'associations': [{'class': 'STALE_DRUG'}]}}
    sample['metadata']['characterization'] = evidence
    window.project.set_metadata(sid, sample['metadata'])
    activate(window, [sid])
    result = window._last_comparison[0]
    labels = window.tree._metadata_values(result)
    assert labels['characterization.state'] == 'stale'
    assert 'STALE_VIRULENCE' not in json.dumps(labels)
    markup = review_report_html(window.report_records(), selected_ids={sid})
    assert 'STALE_VIRULENCE' not in markup and 'Archived or unverified' in markup
    sample['metadata']['characterization']['input_sha256'] = 'a' * 64
    window.project.set_metadata(sid, sample['metadata'])
    current = review_report_html(window.report_records(), selected_ids={sid}, options={**REPORT_PRESETS['cohort'], 'virulence': False, 'drug_associations': False})
    assert 'STALE_VIRULENCE' not in current and 'STALE_DRUG' not in current


def test_threshold_preview_does_not_rewrite_saved_protocol_or_members(window):
    a, b, c = add(window, 'A'), add(window, 'B', '2111'), add(window, 'C', '9999')
    iid = activate(window, [a, b])
    store = InvestigationStore(window.project)
    before = store.get(iid)
    window.cluster_threshold.setValue(4)
    assert window._current_snapshot['preview']
    assert store.get(iid)['snapshots'] == before['snapshots']
    window.cluster_threshold.setValue(1)
    window.cohort_ids = {a, b, c}
    window.refresh_comparison()
    assert window._current_snapshot['preview']
    assert store.get(iid)['sample_ids'] == sorted([a, b])
    assert store.get(iid)['snapshots'] == before['snapshots']


def test_reference_change_starts_new_lineage():
    from wmlstudio.comparison import pairwise_distances
    results = [{'sample_id': sid, 'sample_name': sid, 'scheme': 'one', 'scheme_digest': 'one',
                'alleles': {'a': '1'}, 'status': 'complete'} for sid in ['a', 'b']]
    first = build_snapshot(results, pairwise_distances(results), 1, .95)
    changed = [dict(result, scheme='two', scheme_digest='two') for result in deepcopy(results)]
    second = build_snapshot(changed, pairwise_distances(changed), 1, .95, previous=first)
    assert second['groups'][0]['change'] == 'new'
    assert second['groups'][0]['parents'] == []


def test_threshold_citation_persists_in_snapshot_and_printable_report(window):
    a, b = add(window, 'A'), add(window, 'B')
    iid = activate(window, [a, b])
    store = InvestigationStore(window.project)
    plan = store.get(iid)
    evidence = {'citation': {'title': 'Example <study>', 'doi': '10.test/example'},
                'approved_threshold': 1, 'application_justification': 'Exact synthetic test reference only'}
    store.save(plan['name'], plan['sample_ids'], investigation_id=iid,
               **{key: plan[key] for key in ('scheme', 'scheme_digest', 'threshold', 'min_overlap', 'protocol')},
               threshold_evidence=evidence)
    window.select_investigation(iid)
    snapshot = store.snapshot(iid)
    assert snapshot['threshold_evidence'] == evidence
    markup = review_report_html(window.report_records(), selected_ids={a}, investigation=snapshot)
    assert 'Example &lt;study&gt;' in markup and '10.test/example' in markup


def test_secondary_calling_reuses_only_full_matching_input_hash(qtbot, tmp_path):
    from wmlstudio.typing import call_assembly, load_scheme
    from wmlstudio.ui_compare import ComparisonCallingWorker
    scheme_path = tmp_path / 'scheme'
    scheme_path.mkdir()
    (scheme_path / 'arcA.tfa').write_text('>arcA_1\nAACCGTACGTTAG\n')
    (scheme_path / 'profiles.tsv').write_text('ST\tarcA\n17\t1\n')
    source = tmp_path / 'original.fasta'
    source.write_text('>contig\nAACCGTACGTTAG\n')
    previous = call_assembly(source, load_scheme(scheme_path))
    sample = {'id': 'one', 'name': 'one', 'input_path': str(source), 'status': 'completed'}
    def execute(status='completed'):
        worker = ComparisonCallingWorker([dict(sample, status=status)], str(scheme_path), {'one': [previous]})
        plans, results, failures = [], [], []
        worker.prepared.connect(plans.append)
        worker.sample_finished.connect(lambda sid, result: results.append(result))
        worker.sample_failed.connect(lambda sid, error: failures.append(error))
        with qtbot.waitSignal(worker.finished, timeout=15000):
            worker.start()
        assert not failures
        return plans[0], results
    plan, called = execute()
    assert plan['reused'] == ['one'] and called == []
    source.write_text('>contig\nAACCGTACGTTAG\n>new_contig\nGCGCGC\n')
    plan, called = execute()
    assert plan['pending'] == ['one'] and called[0]['st'] == '17'
    assert called[0]['input_sha256'] != previous['input_sha256']
    previous = called[0]
    plan, called = execute(status='failed')
    assert plan['pending'] == ['one'] and called[0]['st'] == '17'


def test_failed_or_queued_sample_does_not_reuse_archived_secondary(window):
    sid = add(window, 'A')
    result = window.project.get_sample(sid)['result']
    window.project.set_analysis(sid, dict(result, scheme='Secondary', scheme_digest='secondary'))
    for state in ['failed', 'queued', 'running', 'interrupted']:
        window.project.set_status(sid, state)
        assert window.available_profiles(window.project.get_sample(sid)) == []
        assert len(window.project.analysis_results(sid)) == 2
