import json
from collections import Counter
from pathlib import Path

import pytest
from PySide6.QtCore import Qt
from PySide6.QtWidgets import QDialog

from wmlstudio import organism_modules
from wmlstudio.characterization import characterize_assembly, persist_characterization
from wmlstudio.characterization_refs import reference_digest, validate_characterization_references
from wmlstudio.context_menus import SEPARATOR, Selection
from wmlstudio.organism_modules import (
    BOUNDARY,
    CORE_SECTIONS,
    OrganismMatch,
    OrganismModule,
    cohort_applicability,
    cohort_sentence,
    modules_for,
    organism_of,
    record_skipped,
    registered_modules,
    selection_for_record,
    summarize_record,
)
from wmlstudio.project import Project
from wmlstudio.sequence import AnalysisCancelled, file_sha256
from wmlstudio.ui_characterization import CharacterizationPlanDialog, characterization_html


def module_panel(tmp_path, *, sccmec=None, name='references'):
    """A minimal format-2 snapshot; organism-module sections are declared, not implied."""
    root = tmp_path / name
    root.mkdir()
    (root / 'marker.fasta').write_text('>mecA__x__I\nATGAAATAA\n')
    manifest = {'format_version': 2, 'source_revision': 'synthetic-truth',
                'source_repository': 'synthetic-test-fixture',
                'sources': {'synthetic': {'repository': 'synthetic-test-fixture', 'revision': 'synthetic-truth',
                                          'license': 'synthetic', 'license_file': 'source-LICENSE-synthetic'}},
                'species': [], 'virulence': {}, 'locus_profiles': {},
                'sccmec': sccmec if sccmec is not None else {}, 'capsule': {}, 'files': []}
    for path in sorted(root.iterdir()):
        manifest['files'].append({'path': path.name, 'bytes': path.stat().st_size, 'sha256': file_sha256(path)})
    manifest['reference_digest'] = reference_digest(manifest)
    (root / 'manifest.json').write_text(json.dumps(manifest))
    return root


@pytest.fixture
def registered():
    """Register throw-away modules and always remove them, however the test ends."""
    added = []

    def make(key='probe', *, runner=None, manifest_sections=(), match=None, option_keys=('threads',),
             locus_st_provider=None):
        def default_runner(path, reference_root, cancelled=None, progress=None, **options):
            return {'status': 'completed', 'observed_options': sorted(options)}
        module = OrganismModule(
            key=key, title=f'{key} assay', column_title=key.title(),
            match=match or OrganismMatch(genera=frozenset({'Testgenus'}), species=frozenset({'target'})),
            runner=runner or default_runner, option_keys=option_keys, manifest_sections=manifest_sections,
            summary=lambda evidence: f'{key}:{evidence.get("status")}',
            detail_html=lambda evidence: f'<p>{key}</p>', locus_st_provider=locus_st_provider)
        organism_modules._load()
        organism_modules.register(module)
        added.append(key)
        return module

    yield make
    for key in added:
        organism_modules.REGISTRY.pop(key, None)


def assembly_at(tmp_path, name='a.fasta'):
    path = tmp_path / name
    path.write_text('>c\n' + 'ACGT' * 100 + '\n')
    return path


def test_match_rule_separates_recommended_possible_and_off_panel_organisms():
    match = OrganismMatch(genera=frozenset({'Staphylococcus'}), species=frozenset({'aureus'}))
    assert match.evaluate('Staphylococcus', 'aureus')[0] == 'recommended'
    assert match.evaluate('staphylococcus', 'AUREUS')[0] == 'recommended'
    state, reason = match.evaluate('Staphylococcus', 'epidermidis')
    assert state == 'possible' and 'outside the taxa this panel was curated against' in reason
    state, reason = match.evaluate('Klebsiella', 'pneumoniae')
    assert state == 'off_panel' and 'Klebsiella' in reason
    excluded = OrganismMatch(genera=frozenset({'Staphylococcus'}), exclude_species=frozenset({'lugdunensis'}))
    assert excluded.evaluate('Staphylococcus', 'lugdunensis')[0] == 'off_panel'
    assert excluded.evaluate('Staphylococcus', 'epidermidis')[0] == 'recommended'


def test_unknown_organism_never_blocks_a_module_but_is_recorded_as_unknown(tmp_path, registered):
    registered()
    state, reason = OrganismMatch(genera=frozenset({'Testgenus'})).evaluate('', '')
    assert state == 'unknown_organism' and 'No genus or species' in reason
    path = assembly_at(tmp_path)
    root = module_panel(tmp_path)
    result = characterize_assembly(path, reference_root=root, species=False, virulence=False,
                                   modules={'probe': True}, organism=None)
    assert result['probe']['status'] == 'completed'
    assert result['probe']['applicability'] == 'unknown_organism'
    assert 'No genus or species' in result['probe']['applicability_reason']


def test_off_panel_organisms_still_run_the_module_and_are_stamped(tmp_path, registered):
    registered()
    result = characterize_assembly(assembly_at(tmp_path), reference_root=module_panel(tmp_path),
                                   species=False, virulence=False, modules={'probe': True},
                                   organism=('Klebsiella', 'pneumoniae'))
    assert result['probe']['status'] == 'completed' and result['probe']['applicability'] == 'off_panel'


def test_registry_refuses_a_key_colliding_with_a_core_characterization_section(registered):
    for key in sorted(CORE_SECTIONS):
        with pytest.raises(ValueError, match='collides with a core characterization section'):
            registered(key=key)
    registered(key='probe')
    with pytest.raises(ValueError, match='already registered'):
        registered(key='probe')
    with pytest.raises(ValueError, match='unsupported options'):
        registered(key='probe_two', option_keys=('threads', 'reference_root'))


def test_selected_module_missing_its_reference_section_is_not_run_with_a_named_reason(tmp_path, registered):
    registered(manifest_sections=('sccmec.targets', 'sccmec.rules'))
    result = characterize_assembly(assembly_at(tmp_path), reference_root=module_panel(tmp_path),
                                   species=False, virulence=False, modules={'probe': True})
    assert result['probe']['status'] == 'not_run'
    assert 'sccmec.targets' in result['probe']['reason'] and 'sccmec.rules' in result['probe']['reason']
    assert 'format 2' in result['probe']['reason']
    # A missing section is never presented as a negative result.
    assert 'not_detected' not in json.dumps(result['probe'])


def test_an_unselected_module_reports_not_run_rather_than_being_absent(tmp_path, registered):
    registered()
    result = characterize_assembly(assembly_at(tmp_path), reference_root=module_panel(tmp_path),
                                   species=False, virulence=False)
    assert result['probe'] == {'status': 'not_run', 'reason': 'Assay not selected.',
                               'applicability': 'unknown_organism',
                               'applicability_reason': result['probe']['applicability_reason']}


def test_module_status_participates_in_the_overall_characterization_status(tmp_path, registered):
    def ambiguous(path, reference_root, cancelled=None, progress=None, **options):
        return {'status': 'ambiguous', 'reason': 'deliberately unresolved'}
    registered(key='probe', runner=ambiguous)
    path, root = assembly_at(tmp_path), module_panel(tmp_path)
    assert characterize_assembly(path, reference_root=root, species=False, virulence=False)['status'] == 'not_run'
    result = characterize_assembly(path, reference_root=root, species=False, virulence=False,
                                   modules={'probe': True})
    assert result['probe']['status'] == 'ambiguous' and result['status'] == 'ambiguous'


def test_module_failure_is_retained_and_cancellation_is_never_swallowed(tmp_path, registered):
    def fail(path, reference_root, cancelled=None, progress=None, **options):
        raise ValueError('Reference corruption')

    def cancel(path, reference_root, cancelled=None, progress=None, **options):
        raise AnalysisCancelled('cancelled')

    registered(key='probe', runner=fail)
    path, root = assembly_at(tmp_path), module_panel(tmp_path)
    result = characterize_assembly(path, reference_root=root, species=False, virulence=False, modules={'probe': True})
    assert result['probe'] == {'status': 'failed', 'reason': 'Reference corruption',
                               'applicability': 'unknown_organism',
                               'applicability_reason': result['probe']['applicability_reason']}
    assert result['status'] == 'failed'
    organism_modules.REGISTRY.pop('probe')
    registered(key='probe', runner=cancel)
    with pytest.raises(AnalysisCancelled):
        characterize_assembly(path, reference_root=root, species=False, virulence=False, modules={'probe': True})


def test_module_receives_only_the_options_it_declared(tmp_path, registered):
    registered(key='probe', option_keys=('threads', 'blastn_path'))
    result = characterize_assembly(assembly_at(tmp_path), reference_root=module_panel(tmp_path),
                                   species=False, virulence=False, modules={'probe': True}, threads=3)
    assert result['probe']['observed_options'] == ['blastn_path', 'threads']


def test_no_module_writes_metadata_organism(tmp_path, registered):
    def opinionated(path, reference_root, cancelled=None, progress=None, **options):
        return {'status': 'completed', 'proposed_organism': {'genus': 'Testgenus', 'species': 'target'}}

    registered(key='probe', runner=opinionated)
    path = assembly_at(tmp_path)
    with Project(tmp_path / 'project.sqlite') as project:
        sample_id = project.add_sample(path)
        result = characterize_assembly(path, reference_root=module_panel(tmp_path), species=False,
                                       virulence=False, modules={'probe': True},
                                       organism=('Testgenus', 'target'))
        record = persist_characterization(project, sample_id, result)
    assert result['probe']['applicability'] == 'recommended'
    # Applicability is stamped, but organism assignment stays an explicit user action.
    assert not (record.get('metadata') or {}).get('organism')
    assert organism_of(record)[2] == 'Unknown'


def test_modules_for_orders_recommended_before_possible_and_off_panel(registered):
    registered(key='probe')
    registered(key='probe_genus', match=OrganismMatch(genera=frozenset({'Testgenus'})))
    ordered = [(state, module.key) for state, _, module in modules_for('Testgenus', 'other')]
    assert ordered.index(('recommended', 'probe_genus')) < ordered.index(('possible', 'probe'))
    assert all(state == 'off_panel' for state, key in ordered if key not in {'probe', 'probe_genus'})


def test_cohort_applicability_counts_every_selected_isolate(registered):
    module = registered(key='probe')
    records = [{'metadata': {'organism': {'genus': 'Testgenus', 'species': 'target'}}},
               {'metadata': {'organism': 'Testgenus other'}},
               {'result': {'identification': {'organism': {'genus': 'Klebsiella', 'species': 'pneumoniae'}}}},
               {}]
    counts = cohort_applicability(module, records)
    assert counts == {'recommended': 1, 'possible': 1, 'off_panel': 1, 'unknown_organism': 1}
    assert sum(counts.values()) == len(records)


def test_summarize_record_reports_only_modules_that_actually_ran(registered):
    registered(key='probe')
    assert summarize_record({'probe': {'status': 'not_run', 'reason': 'Assay not selected.'}}) == ''
    assert 'probe:completed' in summarize_record({'probe': {'status': 'completed'}})


def test_module_panel_fixture_is_a_valid_installed_snapshot(tmp_path):
    manifest = validate_characterization_references(module_panel(tmp_path))
    assert manifest['format_version'] == 2 and manifest['sccmec'] == {}
    assert Path(manifest['files'][0]['path']).name == 'marker.fasta'


def test_cohort_sentence_names_the_isolates_a_module_will_and_will_not_run_on(registered):
    module = registered(key='probe')
    counts = Counter({'recommended': 3, 'possible': 1, 'unknown_organism': 1, 'off_panel': 5})
    skipping = cohort_sentence(module, counts)
    assert skipping.startswith('Will run on 5 isolates of 10 selected')
    assert '3 within the curated taxa' in skipping
    assert '1 same genus, other species' in skipping and '1 no organism assigned yet' in skipping
    assert 'Will not run on the 5 isolates outside these taxa' in skipping
    assert 'not a negative result' in skipping
    including = cohort_sentence(module, counts, off_panel_included=True)
    assert including.startswith('Will run on 10 isolates of 10 selected')
    assert '5 outside these taxa' in including and 'weak evidence in either direction' in including
    single = cohort_sentence(module, Counter({'recommended': 1}))
    assert single == 'Will run on 1 isolate of 1 selected (1 within the curated taxa).'


def test_off_panel_isolates_are_skipped_unless_the_user_asks_for_them(registered):
    registered(key='probe')
    klebsiella = {'metadata': {'organism': {'genus': 'Klebsiella', 'species': 'pneumoniae'}}}
    running, skipped = selection_for_record(klebsiella, {'probe': True})
    assert running == {} and 'outside the taxa' in skipped['probe']
    running, skipped = selection_for_record(klebsiella, {'probe': True}, include_off_panel=True)
    assert running == {'probe': True} and skipped == {}
    # An unassigned isolate is never guessed into a skip: it runs and is stamped unknown.
    running, skipped = selection_for_record({}, {'probe': True})
    assert running == {'probe': True} and skipped == {}
    assert selection_for_record(klebsiella, {'probe': False}) == ({}, {})


def test_a_skipped_isolate_records_why_rather_than_reading_as_not_selected(tmp_path, registered):
    registered(key='probe')
    klebsiella = {'metadata': {'organism': {'genus': 'Klebsiella', 'species': 'pneumoniae'}}}
    running, skipped = selection_for_record(klebsiella, {'probe': True})
    result = characterize_assembly(assembly_at(tmp_path), reference_root=module_panel(tmp_path),
                                   species=False, virulence=False, modules=running,
                                   organism=('Klebsiella', 'pneumoniae'))
    assert result['probe']['reason'] == 'Assay not selected.'
    record_skipped(result, skipped)
    assert result['probe']['status'] == 'not_run'
    assert result['probe']['reason'].startswith('Not run for this isolate.')
    assert 'not a negative result' in result['probe']['reason']
    assert 'Assay not selected' not in result['probe']['reason']


def test_every_registered_module_says_in_plain_words_what_it_does_not_establish():
    modules = registered_modules()
    assert set(modules) >= {'sccmec', 'klebsiella_locus_st', 'klebsiella_capsule'}
    for key, module in modules.items():
        assert 'does not establish' in module.purpose, key
        assert module.limitations and module.column_title
    for project in ('Kleborate', 'Kaptive', 'AMRFinderPlus', 'SCCmecFinder'):
        assert project in BOUNDARY
    assert 'none is equivalent to those tools' in BOUNDARY
    assert 'no susceptibility, phenotype, serotype or transmission conclusion' in BOUNDARY


# --- the workspace these modules surface in ----------------------------------

@pytest.fixture
def window(qtbot, tmp_path):
    from wmlstudio.app import MainWindow
    widget = MainWindow(storage_root=tmp_path / 'workspace')
    qtbot.addWidget(widget)
    widget.show()
    yield widget
    widget.close()


def plan_sample(identifier, name, genus='', species=''):
    organism = {'genus': genus, 'species': species} if genus or species else {}
    return {'id': identifier, 'name': name, 'input_path': f'/tmp/{name}.fasta',
            'metadata': {'organism': organism} if organism else {},
            'result': {'input_sha256': 'a' * 64}}


def installed_snapshot(tmp_path):
    """Only `manifest.json` needs to exist for the plan dialog's own gate."""
    root = tmp_path / 'refs'
    root.mkdir()
    (root / 'manifest.json').write_text('{}')
    return root


def test_the_plan_offers_the_tools_the_cohort_recommends_and_names_who_is_skipped(qtbot, tmp_path):
    samples = [plan_sample('s1', 'SA-1', 'Staphylococcus', 'aureus'),
               plan_sample('s2', 'SA-2', 'Staphylococcus', 'aureus'),
               plan_sample('k1', 'KP-1', 'Klebsiella', 'pneumoniae')]
    dialog = CharacterizationPlanDialog(samples, installed_snapshot(tmp_path))
    qtbot.addWidget(dialog)
    assert dialog.module_boxes['sccmec'].isChecked()
    assert dialog.module_boxes['klebsiella_capsule'].isChecked()
    note = dialog.module_notes['sccmec'].text()
    assert note.startswith('Will run on 2 isolates of 3 selected')
    assert 'Will not run on the 1 isolate outside these taxa' in note
    assert 'not a negative result' in note
    # Why a tool is or is not offered has to be visible, with its source.
    assert dialog.table.item(2, 2).text() == 'Klebsiella pneumoniae (Assigned)'
    assert 'does not establish' in dialog.module_boxes['sccmec'].toolTip()
    dialog.include_off_panel.setChecked(True)
    assert dialog.module_notes['sccmec'].text().startswith('Will run on 3 isolates of 3 selected')
    assert 'weak evidence in either direction' in dialog.module_notes['sccmec'].text()


def test_excluding_an_isolate_recounts_the_tools_and_a_wholly_off_panel_choice_warns(qtbot, tmp_path):
    samples = [plan_sample('s1', 'SA-1', 'Staphylococcus', 'aureus'),
               plan_sample('k1', 'KP-1', 'Klebsiella', 'pneumoniae')]
    dialog = CharacterizationPlanDialog(samples, installed_snapshot(tmp_path))
    qtbot.addWidget(dialog)
    assert dialog.module_warning.text() == ''
    dialog.set_included(['k1'], False)
    assert dialog.module_notes['sccmec'].text().startswith('Will run on 1 isolate of 1 selected')
    # A Klebsiella tool ticked for a Staphylococcus-only cohort must say so.
    dialog.module_boxes['klebsiella_capsule'].setChecked(True)
    dialog.module_toggled('klebsiella_capsule')
    assert 'Klebsiella wzi / wzc capsule markers' in dialog.module_warning.text()
    assert 'weak evidence in either direction' in dialog.module_warning.text()
    dialog.accept()
    assert dialog.plan['modules'] == {'sccmec': True, 'klebsiella_locus_st': False,
                                      'klebsiella_capsule': True}
    assert dialog.plan['modules_off_panel'] is False
    assert dialog.plan['sample_ids'] == ['s1']


def test_organism_tools_refuse_to_run_without_an_installed_snapshot(qtbot, tmp_path):
    dialog = CharacterizationPlanDialog([plan_sample('s1', 'SA-1', 'Staphylococcus', 'aureus')],
                                        tmp_path / 'absent')
    qtbot.addWidget(dialog)
    dialog.species.setChecked(False)
    dialog.virulence.setChecked(False)
    dialog.accept()
    assert dialog.plan is None
    assert 'Organism-specific tools need a characterization snapshot' in dialog.feedback.text()


def test_correcting_the_organism_changes_which_tools_are_offered(qtbot, tmp_path, monkeypatch):
    from wmlstudio.workflow_dialogs import BatchAssignmentDialog
    with Project(tmp_path / 'project.sqlite') as project:
        identifier = project.add_profile('KP-014', {'sample_name': 'KP-014', 'st': '258',
            'input_sha256': 'a' * 64, 'status': 'profile_imported', 'alleles': {}}, {})
        dialog = CharacterizationPlanDialog(list(project.samples()), installed_snapshot(tmp_path),
                                            project=project)
        qtbot.addWidget(dialog)
        assert dialog.table.item(0, 2).text() == 'Unknown organism'
        assert not dialog.module_boxes['klebsiella_capsule'].isChecked()

        def assign(self):
            self.assignments[0].update(genus='Klebsiella', species='pneumoniae', typing_mode='manual')
            return QDialog.DialogCode.Accepted

        monkeypatch.setattr(BatchAssignmentDialog, 'exec', assign)
        dialog.assign_organisms([identifier])
        assert dialog.table.item(0, 2).text() == 'Klebsiella pneumoniae (Assigned)'
        assert dialog.module_boxes['klebsiella_capsule'].isChecked()
        assert not dialog.module_boxes['sccmec'].isChecked()
        assert project.get_sample(identifier)['metadata']['organism']['genus'] == 'Klebsiella'


def characterized(window, name, evidence, *, stored='a' * 64):
    return window.project.add_profile(name, {'sample_name': name, 'st': '20', 'input_sha256': 'a' * 64,
        'scheme': 'MLST', 'scheme_digest': 'reference', 'status': 'profile_imported', 'alleles': {'a': '1'}},
        {'organism': {'genus': 'Staphylococcus', 'species': 'aureus'},
         'characterization': dict(evidence, input_sha256=stored)})


def test_module_columns_read_not_run_and_stale_module_evidence_is_never_shown_as_current(window):
    call = {'status': 'completed', 'type': 'IV', 'candidate_types': ['IV'],
            'ccr_complexes': [{'name': 'ccr Type 2', 'state': 'present', 'kind': 'ccr', 'targets': {}}],
            'mec_classes': [{'name': 'mec Class B', 'state': 'present', 'kind': 'mec', 'targets': {}}]}
    fragmented = {'status': 'ambiguous', 'type': None, 'candidate_types': ['IV'],
                  'contig_fragmented': True, 'ccr_complexes': [], 'mec_classes': []}
    current = characterized(window, 'SA-1', {'status': 'completed', 'sccmec': call})
    stale = characterized(window, 'SA-2', {'status': 'completed', 'sccmec': call}, stored='b' * 64)
    broken = characterized(window, 'SA-3', {'status': 'ambiguous', 'sccmec': fragmented})
    untouched = characterized(window, 'SA-4', {'status': 'not_run'})
    window.feature_ids = {current, stale, broken, untouched}
    window.refresh_features()
    table = window.characterization_table
    headers = [table.horizontalHeaderItem(index).text() for index in range(table.columnCount())]
    assert headers[8:] == ['SCCmec', 'Locus STs', 'Capsule markers']
    column = headers.index('SCCmec')
    values = {table.item(row, 0).text(): table.item(row, column).text() for row in range(table.rowCount())}
    assert values['SA-1'] == 'IV (candidate) - ccr2 + mec B'
    assert values['SA-2'] == 'not_run'
    # A cassette split across contigs reads as withheld, never as the nearest type.
    assert values['SA-3'] == 'ambiguous - IV withheld'
    assert values['SA-4'] == 'not_run'
    assert table.horizontalHeaderItem(column).toolTip().startswith('Looks for the ccr and mec marker genes')


def test_the_results_table_offers_a_right_click_menu_for_the_selected_isolates(window):
    table = window.characterization_table
    assert table.contextMenuPolicy() == Qt.ContextMenuPolicy.CustomContextMenu
    assert 'characterization' in window._context_adapters
    plan = window.context_menu_plan(Selection('characterization', ('s1', 's2')))
    titles = [entry.format_title(Selection('characterization', ('s1', 's2')))
              for entry in plan if entry is not SEPARATOR]
    assert titles, 'the characterization table offered no runnable actions'
    assert any('Copy' in title for title in titles)
    assert any('Export' in title for title in titles)


def test_the_drill_down_shows_the_matches_behind_a_call_and_repeats_the_module_limits():
    hit = {'gene': 'mecA<script>', 'contig': 'contig_1', 'identity_pct': 99.87, 'coverage_pct': 100.0,
           'start': 10, 'end': 2016, 'strand': '+',
           'cds_qc': {'valid': False, 'reasons': ['internal in-frame stop codon']}}
    evidence = {'status': 'ambiguous', 'input_sha256': 'a' * 64,
                'sccmec': {'status': 'ambiguous', 'candidate_types': ['IV', 'IVa'], 'hits': [hit],
                           'mecA': 'detected', 'mecC': 'not_assayed', 'contig_fragmented': True,
                           'limitations': ['mecC is not in this reference panel and is not assayed here.']}}
    sample = {'name': 'SA-1', 'result': {'input_sha256': 'a' * 64},
              'metadata': {'organism': {'genus': 'Staphylococcus', 'species': 'aureus'},
                           'characterization': evidence}}
    body = characterization_html(sample)
    assert 'Organism-specific tools' in body
    assert 'Staphylococcus aureus (Assigned)' in body and BOUNDARY in body
    assert 'mecA&lt;script&gt;' in body and '<script>' not in body
    assert '99.9' in body and '100.0' in body
    assert 'no: internal in-frame stop codon' in body
    assert 'mecC is not in this reference panel and is not assayed here.' in body
    assert 'ambiguous - 2 candidate types' in body and 'IV, IVa' in body
    # An ambiguous cassette is never rendered as a call.
    assert 'IV (candidate)' not in body
    assert 'not proof that the gene is intact' in body


def test_the_drill_down_truncates_a_long_hit_list_and_states_the_total():
    hits = [{'gene': f'ccrB{index}', 'contig': 'c1', 'identity_pct': 90 + index / 100,
             'coverage_pct': 95.0, 'start': index, 'end': index + 10, 'strand': '+',
             'cds_qc': {'valid': True, 'reasons': []}} for index in range(60)]
    evidence = {'status': 'completed', 'input_sha256': 'a' * 64,
                'sccmec': {'status': 'completed', 'type': 'IV', 'candidate_types': ['IV'],
                           'ccr_complexes': [], 'mec_classes': [], 'hits': hits, 'limitations': []}}
    sample = {'name': 'SA-1', 'result': {'input_sha256': 'a' * 64},
              'metadata': {'characterization': evidence}}
    body = characterization_html(sample)
    assert 'Showing the 40 highest-identity matches of 60.' in body
    assert body.count('<tr>') <= 60
    assert 'Unknown organism' in body


def test_a_right_click_run_characterizes_only_the_clicked_isolates(window, tmp_path, monkeypatch):
    identifiers = {}
    for name in ('SA-1', 'SA-2'):
        path = tmp_path / f'{name}.fasta'
        path.write_text('>c\n' + 'ACGT' * 100 + '\n')
        identifiers[name] = window.project.add_sample(path)
    profile_only = characterized(window, 'PR-1', {'status': 'not_run'})
    seen = []

    def decline(self):
        seen.append(sorted(sample['name'] for sample in self.samples))
        return QDialog.DialogCode.Rejected

    monkeypatch.setattr(CharacterizationPlanDialog, 'exec', decline)
    window.context_run_characterization(Selection('characterization', (identifiers['SA-2'],)))
    assert seen == [['SA-2.fasta']]
    # A profile row carries no assembly, so the plan is never opened for it.
    window.context_run_characterization(Selection('characterization', (profile_only,)))
    window.context_run_characterization(Selection('characterization', ()))
    assert seen == [['SA-2.fasta']]


def test_the_run_gives_each_isolate_only_the_tools_its_organism_supports(window, qtbot, tmp_path, monkeypatch):
    from wmlstudio.cohort_picker import CohortPickerDialog
    from wmlstudio.scheduler import plan_resources
    from wmlstudio.storage import assign_organism
    root = module_panel(tmp_path)
    identifiers = {}
    for name, genus, species in (('SA-1', 'Staphylococcus', 'aureus'), ('KP-1', 'Klebsiella', 'pneumoniae')):
        path = tmp_path / f'{name}.fasta'
        path.write_text('>c\n' + 'ACGT' * 100 + '\n')
        identifiers[genus] = window.project.add_sample(path)
        assign_organism(window.project, [identifiers[genus]], genus, species, typing_mode='manual')
    resources = plan_resources(threads_per_sample=1, memory_gb=1, policy='low_memory').to_dict()

    def accept_picker(self):
        self.selected_ids = set(identifiers.values())
        return QDialog.DialogCode.Accepted

    def accept_plan(self):
        self.plan = {'sample_ids': sorted(identifiers.values()), 'reference_path': str(root),
                     'species': False, 'virulence': False, 'hydra': False, 'resources': resources,
                     'modules': {'sccmec': True, 'klebsiella_capsule': True, 'klebsiella_locus_st': False},
                     'modules_off_panel': False}
        return QDialog.DialogCode.Accepted

    monkeypatch.setattr(CohortPickerDialog, 'exec', accept_picker)
    monkeypatch.setattr(CharacterizationPlanDialog, 'exec', accept_plan)
    window.run_characterization_selected()
    qtbot.waitUntil(lambda: window.worker is not None and not window.worker.isRunning(), timeout=30000)
    stored = {genus: window.project.get_sample(identifier)['metadata']['characterization']
              for genus, identifier in identifiers.items()}
    # The Staphylococcus isolate got the Staphylococcus tool, and the snapshot gap is named.
    assert stored['Staphylococcus']['sccmec']['applicability'] == 'recommended'
    assert 'sccmec.targets' in stored['Staphylococcus']['sccmec']['reason']
    # The Klebsiella tool was chosen but is not run on a Staphylococcus isolate, and says why.
    skipped = stored['Staphylococcus']['klebsiella_capsule']
    assert skipped['applicability'] == 'off_panel'
    assert skipped['reason'].startswith('Not run for this isolate.')
    assert 'not a negative result' in skipped['reason']
    assert stored['Klebsiella']['sccmec']['reason'].startswith('Not run for this isolate.')
    assert 'capsule.loci' in stored['Klebsiella']['klebsiella_capsule']['reason']
    # A module the user never chose is still distinguishable from one that was skipped.
    assert stored['Klebsiella']['klebsiella_locus_st']['reason'] == 'Assay not selected.'
