import json
from pathlib import Path

import pytest

from wmlstudio import organism_modules
from wmlstudio.characterization import characterize_assembly, persist_characterization
from wmlstudio.characterization_refs import reference_digest, validate_characterization_references
from wmlstudio.organism_modules import (
    CORE_SECTIONS,
    OrganismMatch,
    OrganismModule,
    cohort_applicability,
    modules_for,
    organism_of,
    summarize_record,
)
from wmlstudio.project import Project
from wmlstudio.sequence import AnalysisCancelled, file_sha256


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
