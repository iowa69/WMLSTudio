import json
import os
import random
from pathlib import Path

import pytest

from wmlstudio import characterization
from wmlstudio.characterization import (
    characterize_assembly,
    current_characterization,
    hydra_report_for_record,
    persist_characterization,
    synthesize_accessory_evidence,
)
from wmlstudio.characterization_refs import reference_digest, validate_characterization_references
from wmlstudio.project import Project
from wmlstudio.sample_workflow import current_input_sha256, link_hydra
from wmlstudio.sequence import AnalysisCancelled, file_sha256
from wmlstudio.species_evidence import identify_species, interpret_species_hits
from wmlstudio.virulence_evidence import (
    _unique_locations,
    screen_virulence,
    summarize_virulence_hits,
)


def reference_panel(tmp_path, genomes=None, genes=None, format_version=1):
    root = tmp_path / 'references'
    root.mkdir()
    manifest = {'format_version': format_version, 'source_revision': 'synthetic-truth',
                'source_repository': 'synthetic-test-fixture', 'species': [], 'virulence': {}, 'files': []}
    if format_version >= 2:
        manifest.update(sources={'synthetic': {'repository': 'synthetic-test-fixture', 'revision': 'synthetic-truth',
                                               'license': 'synthetic', 'license_file': 'source-LICENSE-synthetic'}},
                        locus_profiles={}, sccmec={}, capsule={})
    for index, (taxon, sequence) in enumerate(genomes or []):
        relative = f'ref{index}.fasta'
        (root / relative).write_text(f'>ref{index}\n{sequence}\n')
        manifest['species'].append({'id': f'ref{index}', 'path': relative, 'genus': 'Testgenus',
                                    'species': taxon, 'subspecies': '', 'outgroup': index != 0})
    for gene, sequence in (genes or {}).items():
        relative = f'{gene}.fasta'
        (root / relative).write_text(f'>{gene}_1\n{sequence}\n')
        manifest['virulence'].setdefault('test-locus', []).append({'gene': gene, 'path': relative})
    for path in sorted(root.iterdir()):
        manifest['files'].append({'path': path.name, 'bytes': path.stat().st_size, 'sha256': file_sha256(path)})
    manifest['reference_digest'] = reference_digest(manifest)
    (root / 'manifest.json').write_text(json.dumps(manifest))
    return root


def species_refs():
    return [{'id': 'kp', 'genus': 'Klebsiella', 'species': 'pneumoniae'},
            {'id': 'kq', 'genus': 'Klebsiella', 'species': 'quasipneumoniae'},
            {'id': 'ec', 'genus': 'Escherichia', 'species': 'coli', 'outgroup': True}]


def species_hit(reference_id, ani, qf=.9, rf=.9):
    return {'reference_id': reference_id, 'ani': ani, 'query_fraction': qf, 'reference_fraction': rf}


def test_species_kernel_requires_independent_ani_fraction_and_gap_evidence():
    result = interpret_species_hits([species_hit('kq', 93), species_hit('kp', 99)], species_refs())
    assert result['status'] == 'completed' and result['species'] == 'pneumoniae'
    assert result['gap_ani'] == 6 and result['subspecies'] == ''
    result = interpret_species_hits([species_hit('kp', 99), species_hit('kq', 98.8)], species_refs())
    assert result['status'] == 'ambiguous' and result['species'] == ''
    assert 'Competing' in result['reason']
    result = interpret_species_hits([species_hit('kp', 100, .49, 1), species_hit('kq', 100, .49, 1)], species_refs())
    assert result['status'] == 'ambiguous' and result['species'] == ''
    result = interpret_species_hits([species_hit('kp', 94.9)], species_refs())
    assert result['status'] == 'ambiguous'
    result = interpret_species_hits([species_hit('kp', 99, .75, .95), species_hit('kq', 99.5, .2, .95)], species_refs())
    assert result['status'] == 'ambiguous' and 'possible mixture' in result['reason']
    assert result['discordant_low_fraction_hits'][0]['reference_id'] == 'kq'


def test_subspecies_is_only_provisional_and_requires_competing_subspecies():
    references = [{'id': 'a', 'genus': 'Klebsiella', 'species': 'quasipneumoniae', 'subspecies': 'quasipneumoniae'},
                  {'id': 'b', 'genus': 'Klebsiella', 'species': 'quasipneumoniae', 'subspecies': 'similipneumoniae'}]
    result = interpret_species_hits([species_hit('a', 99.5)], references)
    assert result['subspecies'] == ''
    result = interpret_species_hits([species_hit('a', 99.5), species_hit('b', 96.5)], references)
    assert result['subspecies'] == 'quasipneumoniae'
    assert result['subspecies_status'] == 'provisional_reference_match'
    result = interpret_species_hits([species_hit('a', 99.5), species_hit('b', 99.4)], references)
    assert result['subspecies'] == ''


@pytest.mark.parametrize('hit', [species_hit('kp', float('nan')), species_hit('kp', 99, 1.01), species_hit('unknown', 99)])
def test_species_rejects_invalid_output(hit):
    with pytest.raises(ValueError):
        interpret_species_hits([hit], species_refs())


def test_ecoli_does_not_silently_exclude_shigella():
    result = interpret_species_hits([species_hit('ec', 100)], species_refs())
    assert result['status'] == 'ambiguous' and result['confidence'] == 'complex_only'
    assert 'Shigella' in result['reason']


def test_real_pyskani_computation_and_changed_reference_rejection(tmp_path):
    random_source = random.Random(922)
    first = ''.join(random_source.choices('ACGT', k=200_000))
    second = ''.join(random_source.choices('ACGT', k=200_000))
    references = reference_panel(tmp_path, [('target', first), ('outgroup', second)])
    assembly = tmp_path / 'assembly.fasta'
    assembly.write_text('>one\n' + first + '\n')
    result = identify_species(assembly, references)
    assert result['status'] == 'completed' and result['species'] == 'target'
    assert result['nearest']['ani'] > 99.9
    assert result['nearest']['query_fraction'] > .99
    assert result['provenance']['version'] == '0.2.0'
    assert result['provenance']['skani_version'] == '0.3.0'
    assert len(result['provenance']['binary_sha256']) == 64
    # Assembly compression/filename changes do not become biological evidence.
    second_result = identify_species(assembly, references)
    assert second_result['hits'] == result['hits']
    (references / 'ref0.fasta').write_text('>changed\nACGT\n')
    with pytest.raises(ValueError, match='changed or is missing'):
        identify_species(assembly, references)


def test_species_cancellation_and_no_outgroup_are_not_valid_calls(tmp_path):
    with pytest.raises(AnalysisCancelled):
        identify_species(tmp_path / 'missing', tmp_path / 'missing', cancelled=lambda: True)
    refs = reference_panel(tmp_path, [('target', 'ACGT' * 50)])
    assembly = tmp_path / 'a.fasta'
    assembly.write_text('>a\nACGT\n')
    with pytest.raises(ValueError, match='outgroup'):
        identify_species(assembly, refs)


def test_virulence_no_hit_is_not_absence_when_assay_qc_fails():
    loci = {'iuc': [{'gene': 'iucA'}, {'gene': 'iucB'}]}
    assert summarize_virulence_hits([], loci, adequate_negative_assay=False)[0]['status'] == 'ambiguous'
    adequate = summarize_virulence_hits([], loci, adequate_negative_assay=True)[0]
    assert adequate['status'] == 'not_detected' and adequate['phenotype'] == 'not_inferred'
    hit = {'gene': 'iucA', 'identity_pct': 100, 'coverage_pct': 90, 'cds_qc': {'valid': False}}
    result = summarize_virulence_hits([hit], loci, adequate_negative_assay=True)[0]
    assert result['status'] == 'detected' and result['completeness'] == 'incomplete_or_unresolved'
    assert result['intact_genes_detected'] == 0 and result['official_locus_st'] is None


def test_virulence_dedup_preserves_distinct_physical_copies():
    base = {'gene': 'iucA', 'contig': 'c', 'strand': '+', 'start': 10, 'end': 109,
            'bitscore': 100, 'identity_pct': 100, 'coverage_pct': 100, 'reference_allele': 'iucA_1'}
    alternate = dict(base, reference_allele='iucA_2', identity_pct=99, bitscore=98)
    other_copy = dict(base, start=200, end=299)
    assert len(_unique_locations([alternate, other_copy, base])) == 2


def amr_hit(gene='blaKPC-2', kind='AMR', **changes):
    return dict({'gene': gene, 'element_type': kind, 'database': 'ncbi', 'primary': True,
                 'class': 'BETA-LACTAM', 'subclass': 'CARBAPENEM', 'method': 'BLASTN', 'resolution': 'COMPLETE',
                 'sequence': 'c', 'start': 101, 'end': 200, 'coverage_pct': 100, 'identity_pct': 100}, **changes)


def test_drug_associations_preserve_curated_classes_without_inventing_ast():
    sample = {'hits': [amr_hit(), amr_hit('partial', resolution='PARTIAL'), amr_hit('uncurated', **{'class': '', 'subclass': ''})]}
    result = synthesize_accessory_evidence(sample, {'status': 'completed', 'databases': ['ncbi']}, {'c': 2000})
    drugs = result['drug_associations']
    assert [row['gene'] for row in drugs['associations']] == ['blaKPC-2']
    assert drugs['associations'][0]['affected_drugs'] == []
    assert drugs['associations'][0]['subclass'] == 'CARBAPENEM'
    assert len(drugs['requires_review']) == 2 and drugs['phenotype'] == 'not_inferred'
    assert result['plasmid_hypotheses']['status'] == 'not_run'


def test_specific_drug_labels_require_explicit_reference_tokens_not_class_or_gene_guesses():
    sample = {'hits': [amr_hit('aac-truth', **{'class': 'AMINOGLYCOSIDE', 'subclass': 'AMIKACIN/GENTAMICIN'}),
                       amr_hit('blaKPC-2'), amr_hit('unknown', **{'class': 'MEROPENEM', 'subclass': 'UNKNOWN'}),
                       amr_hit('porin', resolution='POINT', method='POINT_DISRUPT')]}
    result = synthesize_accessory_evidence(sample, {'databases': ['ncbi']}, {'contig': 10000})['drug_associations']
    assert result['associations'][0]['affected_drugs'] == ['AMIKACIN', 'GENTAMICIN']
    assert result['associations'][1]['affected_drugs'] == []  # CARBAPENEM is not a specific drug.
    assert result['associations'][2]['affected_drugs'] == []  # No class or gene-name inference.
    assert result['requires_review'][0]['gene'] == 'porin'
    assert result['phenotype'] == 'not_inferred'


def test_plasmid_hypothesis_requires_actual_same_contig_coordinates():
    replicon = amr_hit('IncFIB', 'PLASMID', database='plasmidfinder', start=301, end=400)
    sample = {'hits': [amr_hit(), replicon, amr_hit('other_contig', sequence='other')]}
    result = synthesize_accessory_evidence(sample, {'status': 'completed', 'databases': ['ncbi', 'plasmidfinder']}, {'c': 2000, 'other': 2000})
    links = result['plasmid_hypotheses']['contig_associations']
    assert len(links) == 1 and links[0]['marker'] == 'blaKPC-2' and links[0]['gap_bp'] == 100
    assert links[0]['status'] == 'hypothesis' and result['plasmid_hypotheses']['reconstruction'] == 'not_run'
    sample['hits'][1]['start'] = 3000
    assert not synthesize_accessory_evidence(sample, {'databases': ['plasmidfinder']}, {'c': 2000})['plasmid_hypotheses']['contig_associations']


def test_characterization_hashless_or_stale_hydra_is_not_promoted(tmp_path):
    path = tmp_path / 'a.fasta'
    path.write_text('>c\n' + 'ACGT' * 100 + '\n')
    report = {'hydra_version': '1.4.0', 'samples': [{'sample': 'a', 'hits': [amr_hit()]}], 'databases': ['ncbi']}
    result = characterize_assembly(path, hydra_report=report)
    assert result['drug_associations']['status'] == 'ambiguous'
    assert result['drug_associations']['associations'] == []
    assert result['virulence']['status'] == 'not_run'
    report['execution_provenance'] = {'inputs': [{'sha256': '0' * 64}]}
    result = characterize_assembly(path, hydra_report=report)
    assert result['drug_associations']['status'] == 'failed'
    report['execution_provenance']['inputs'][0]['sha256'] = file_sha256(path)
    result = characterize_assembly(path, hydra_report=report)
    assert result['drug_associations']['status'] == 'completed'
    assert result['drug_associations']['associations'][0]['gene'] == 'blaKPC-2'


def test_module_failure_retained_and_cancellation_never_swallowed(tmp_path, monkeypatch):
    path = tmp_path / 'a.fasta'
    path.write_text('>c\nACGT\n')
    def fail(*args, **kwargs):
        raise ValueError('Reference corruption')
    monkeypatch.setattr(characterization, 'identify_species', fail)
    result = characterize_assembly(path, reference_root=tmp_path, virulence=False)
    assert result['status'] == 'failed' and result['species_evidence']['reason'] == 'Reference corruption'
    def cancel(*args, **kwargs):
        raise AnalysisCancelled('cancelled')
    monkeypatch.setattr(characterization, 'identify_species', cancel)
    with pytest.raises(AnalysisCancelled):
        characterize_assembly(path, reference_root=tmp_path, virulence=False)


def test_characterization_persistence_preserves_typing_state_and_history(tmp_path):
    path = tmp_path / 'a.fasta'
    path.write_text('>c\nACGT\n')
    with Project(tmp_path / 'project.sqlite') as project:
        sid = project.add_sample(path)
        before = project.get_sample(sid)['status']
        result = characterize_assembly(path)
        first = persist_characterization(project, sid, result)
        assert first['status'] == before and first['result'] is None
        assert current_input_sha256(first) == file_sha256(path)
        assert current_characterization(first)['status'] == 'current'
        persist_characterization(project, sid, result)
        assert any(row['action'] == 'characterization_superseded' for row in project.history(sid))
        project.update_metadata(sid, {'assembly': {'provenance': {'assembly_sha256': 'f' * 64}}})
        stale = current_characterization(project.get_sample(sid))
        assert stale['status'] == 'stale' and stale['evidence'] is None
        with pytest.raises(ValueError, match='recorded current'):
            persist_characterization(project, sid, result)


def test_characterization_rejects_input_changed_before_attach(tmp_path):
    path = tmp_path / 'a.fasta'
    path.write_text('>c\nACGT\n')
    with Project(tmp_path / 'project.sqlite') as project:
        sid = project.add_sample(path)
        result = characterize_assembly(path)
        path.write_text('>c\nTGCA\n')
        with pytest.raises(ValueError, match='bytes changed'):
            persist_characterization(project, sid, result)
        assert 'characterization' not in project.get_sample(sid)['metadata']


def test_per_record_hydra_bridge_keeps_source_and_never_hashes_user_mapping(tmp_path):
    path = tmp_path / 'a.fasta'
    path.write_text('>c\nACGT\n')
    digest = file_sha256(path)
    report = {'hydra_version': '1.4.0', 'samples': [{'sample': 'custom', 'hits': [amr_hit()]}],
              'databases': ['ncbi'], 'import_provenance': {'sha256': 'r' * 64},
              'execution_provenance': {'version': '1.4.0', 'inputs': [{'sha256': digest}]}}
    with Project(tmp_path / 'project.sqlite') as project:
        sid = project.add_sample(path)
        persist_characterization(project, sid, characterize_assembly(path))
        link_hydra(project, report, {'custom': sid})
        reconstructed = hydra_report_for_record(project.get_sample(sid))
        assert reconstructed['databases'] == ['ncbi']
        assert reconstructed['samples'][0]['input_sha256'] == digest
        assert reconstructed['execution_provenance'] == report['execution_provenance']
        record = project.get_sample(sid)
        record['metadata']['hydra']['evidence_input_sha256'] = None
        record['metadata']['hydra']['execution_provenance'] = {}
        assert hydra_report_for_record(record) is None


@pytest.mark.parametrize('format_version', [1, 2])
def test_reference_manifest_does_not_trust_modified_sequence_bytes(tmp_path, format_version):
    root = reference_panel(tmp_path, genes={'a': 'ATGAAATAA'}, format_version=format_version)
    assert validate_characterization_references(root)['format_version'] == format_version
    (root / 'a.fasta').write_text('>a_1\nATGCCCTAA\n')
    with pytest.raises(ValueError, match='changed or is missing'):
        validate_characterization_references(root)


def test_format_version_one_snapshot_still_validates_and_new_modules_report_not_run(tmp_path):
    root = reference_panel(tmp_path, genes={'a': 'ATGAAATAA'})
    manifest = validate_characterization_references(root)
    assert manifest['format_version'] == 1 and 'sccmec' not in manifest
    assembly = tmp_path / 'a.fasta'
    assembly.write_text('>c\nACGT\n')
    result = characterize_assembly(assembly, reference_root=root, species=False, virulence=False,
                                   modules={'sccmec': True, 'klebsiella_locus_st': True})
    for key in ('sccmec', 'klebsiella_locus_st'):
        assert result[key]['status'] == 'not_run'
        assert 'format 1' in result[key]['reason'] and key.split('_')[0] in result[key]['reason'].lower()


def test_v2_manifest_missing_a_declared_section_is_rejected(tmp_path):
    root = reference_panel(tmp_path, genes={'a': 'ATGAAATAA'}, format_version=2)
    manifest = json.loads((root / 'manifest.json').read_text())
    del manifest['capsule']
    (root / 'manifest.json').write_text(json.dumps(manifest))
    with pytest.raises(ValueError, match='omits a declared organism-module section'):
        validate_characterization_references(root)


def test_reference_manifest_rejects_a_stripped_organism_module_section(tmp_path):
    root = reference_panel(tmp_path, genes={'a': 'ATGAAATAA'}, format_version=2)
    manifest = json.loads((root / 'manifest.json').read_text())
    manifest['sccmec'] = {'targets': [{'gene': 'mecA', 'path': 'a.fasta'}]}
    (root / 'manifest.json').write_text(json.dumps(manifest))
    with pytest.raises(ValueError, match='fingerprint is invalid'):
        validate_characterization_references(root)


def test_v2_manifest_requires_a_source_hash_for_every_organism_module_reference(tmp_path):
    from wmlstudio.characterization_refs import reference_digest as digest
    root = reference_panel(tmp_path, genes={'a': 'ATGAAATAA'}, format_version=2)
    manifest = json.loads((root / 'manifest.json').read_text())
    manifest['capsule'] = {'loci': [{'gene': 'wzi', 'path': 'capsule/wzi.fasta'}]}
    manifest['reference_digest'] = digest(manifest)
    (root / 'manifest.json').write_text(json.dumps(manifest))
    with pytest.raises(ValueError, match='no source hash'):
        validate_characterization_references(root)


def test_reference_provision_failure_or_cancel_never_publishes_partial_snapshot(tmp_path, monkeypatch):
    from wmlstudio import characterization_refs
    root = tmp_path / 'installed'
    calls = []
    def broken_fetch(source_key, relative, target, cancelled=None, **options):
        calls.append((source_key, relative))
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(b'partial network response')
        raise OSError('connection lost')
    monkeypatch.setattr(characterization_refs, '_fetch', broken_fetch)
    with pytest.raises(OSError, match='connection lost'):
        characterization_refs.provision_characterization_references(root)
    assert calls and list(root.iterdir()) == []
    with pytest.raises(AnalysisCancelled):
        characterization_refs.provision_characterization_references(root, cancelled=lambda: True)
    assert list(root.iterdir()) == []


def test_pyskani_frozen_collection_keeps_actual_native_module():
    from PyInstaller.utils.hooks import collect_submodules
    assert 'pyskani._skani' in collect_submodules('pyskani')
    recipe = (Path(__file__).resolve().parents[1] / 'studio_packaging/wmlstudio.spec').read_text()
    assert 'collect_submodules("pyskani")' in recipe
    assert '"pyskani"' in recipe and 'characterization_database' in recipe


def test_native_blast_virulence_complete_and_disrupted_truth(tmp_path):
    tools = os.environ.get('WMLSTUDIO_TEST_BLAST_DIR')
    if not tools:
        pytest.skip('Set WMLSTUDIO_TEST_BLAST_DIR for the real native BLAST assay.')
    random_source = random.Random(421)
    codons = ['AAA', 'GCT', 'TTC', 'GGC', 'CCA', 'GTT', 'CAT']
    gene = 'ATG' + ''.join(random_source.choices(codons, k=150)) + 'TAA'
    refs = reference_panel(tmp_path, genes={'iucA': gene})
    assembly = tmp_path / 'isolate.fasta'
    assembly.write_text('>c\n' + 'N' * 100 + gene + 'N' * 100 + '\n')
    suffix = '.exe' if os.name == 'nt' else ''
    options = {'blastn_path': Path(tools) / ('blastn' + suffix), 'makeblastdb_path': Path(tools) / ('makeblastdb' + suffix)}
    result = screen_virulence(assembly, refs, **options)
    gene_result = result['loci'][0]['genes'][0]
    assert gene_result['status'] == 'detected' and gene_result['intact_cds_copies'] == 1
    assert gene_result['hits'][0]['start'] == 101
    assert result['assay']['adequate_negative_assay'] is False
    damaged = gene[:90] + 'TAA' + gene[93:]
    assembly.write_text('>c\n' + 'N' * 100 + damaged + 'N' * 100 + '\n')
    result = screen_virulence(assembly, refs, **options)
    gene_result = result['loci'][0]['genes'][0]
    assert gene_result['status'] == 'detected' and gene_result['intact_cds_copies'] == 0
    assert 'internal in-frame stop codon' in gene_result['hits'][0]['cds_qc']['reasons']


def isolate_record(identifier='iso-1', name='Isolate 1', genus='Staphylococcus', species='aureus'):
    """A project sample row in the shape the characterization plan dialog reads."""
    return {'id': identifier, 'name': name, 'input_path': f'/inputs/{identifier}.fasta',
            'metadata': {'organism': {'genus': genus, 'species': species}},
            'result': {'kind': 'fasta'}}


def test_a_characterization_run_is_refused_when_hydra_has_no_reference_data(qtbot, tmp_path):
    """Refuse up front instead of failing each isolate after its other assays ran."""
    from wmlstudio.ui_characterization import CharacterizationPlanDialog

    empty = tmp_path / 'amr-store'
    empty.mkdir()
    dialog = CharacterizationPlanDialog([isolate_record()], database_root=str(empty))
    qtbot.addWidget(dialog)
    # Whatever else this machine has installed, an empty store always blocks a run.
    assert dialog.hydra_check is not None and dialog.hydra_check['ready'] is False
    assert 'every gene and mutation this screen can name' in dialog.hydra_state.text()
    # Not offered as a choice that would fail: switched off, disabled, reason shown.
    assert dialog.hydra.isEnabled() is False and dialog.hydra.isChecked() is False
    assert 'HYDRA cannot start' in dialog.hydra.toolTip()
    dialog.species.setChecked(False)
    dialog.virulence.setChecked(False)
    for control in dialog.module_boxes.values():
        control.setChecked(False)
    dialog.hydra.setChecked(True)
    dialog.accept()
    assert dialog.plan is None
    assert 'HYDRA cannot start' in dialog.feedback.text()
    assert 'Nothing is ever downloaded during a run.' in dialog.feedback.text()


def test_the_other_assays_still_run_when_hydra_is_not_part_of_the_plan(qtbot, tmp_path):
    from wmlstudio.ui_characterization import CharacterizationPlanDialog

    empty = tmp_path / 'amr-store'
    empty.mkdir()
    dialog = CharacterizationPlanDialog([isolate_record()], database_root=str(empty))
    qtbot.addWidget(dialog)
    dialog.hydra.setChecked(False)
    dialog.species.setChecked(False)
    dialog.virulence.setChecked(False)
    for control in dialog.module_boxes.values():
        control.setChecked(False)
    dialog.accept()
    assert dialog.plan is not None and dialog.plan['hydra'] is False
    assert dialog.plan['sample_ids'] == ['iso-1']


def test_an_isolate_outside_every_installed_catalogue_is_named_before_the_run(tmp_path):
    """Gene screening still runs for it; that is not evidence it carries no mutation."""
    from wmlstudio.ui_characterization import CharacterizationPlanDialog

    accepted = ['Escherichia', 'Staphylococcus_aureus']
    covered = CharacterizationPlanDialog.organism_covered(
        isolate_record(genus='Staphylococcus', species='aureus'), accepted)
    genus_level = CharacterizationPlanDialog.organism_covered(
        isolate_record(genus='Escherichia', species='coli'), accepted)
    outside = CharacterizationPlanDialog.organism_covered(
        isolate_record(genus='Listeria', species='monocytogenes'), accepted)
    unknown = CharacterizationPlanDialog.organism_covered(
        isolate_record(genus='', species=''), accepted)
    assert covered is True and genus_level is True
    assert outside is False and unknown is False


def test_a_transient_download_failure_does_not_lose_the_whole_snapshot(tmp_path, monkeypatch):
    """Staging fetches hundreds of files, so one dropped connection cost everything.

    A rate-limit reply or a short read is worth retrying; a 404 is not, and
    retrying it would only delay an answer the user needs. A retried attempt must
    also clear the partial file, because the writer creates its target exclusively.
    """
    import http.client
    import io
    import urllib.error

    from wmlstudio import characterization_refs as refs

    monkeypatch.setattr(refs, "_wait_before_retry", lambda *args, **kwargs: None)
    attempts = {"count": 0}

    def rate_limited_once(request, timeout=None):
        attempts["count"] += 1
        if attempts["count"] == 1:
            raise urllib.error.HTTPError(request.full_url, 429, "Too Many Requests",
                                         {"Retry-After": "1"}, None)
        return io.BytesIO(b">seq\nACGT\n")

    monkeypatch.setattr(refs.urllib.request, "urlopen", rate_limited_once)
    record = refs._fetch("kleborate", "x.fasta", tmp_path / "a.fasta")
    assert attempts["count"] == 2 and record["bytes"] == 10
    assert (tmp_path / "a.fasta").is_file()

    class Truncated(io.BytesIO):
        def read(self, size=-1):
            raise http.client.IncompleteRead(b"partial")

    attempts["count"] = 0

    def cut_then_complete(request, timeout=None):
        attempts["count"] += 1
        return Truncated(b"") if attempts["count"] == 1 else io.BytesIO(b">seq\nTTTT\n")

    monkeypatch.setattr(refs.urllib.request, "urlopen", cut_then_complete)
    assert refs._fetch("kleborate", "y.fasta", tmp_path / "b.fasta")["bytes"] == 10

    attempts["count"] = 0

    def missing(request, timeout=None):
        attempts["count"] += 1
        raise urllib.error.HTTPError(request.full_url, 404, "Not Found", {}, None)

    monkeypatch.setattr(refs.urllib.request, "urlopen", missing)
    with pytest.raises(urllib.error.HTTPError):
        refs._fetch("kleborate", "z.fasta", tmp_path / "c.fasta")
    assert attempts["count"] == 1, "a permanent refusal must not be retried"
