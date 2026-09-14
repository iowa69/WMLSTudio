import html
import json
import random

import pytest

from wmlstudio.characterization_refs import _capsule_split_name, _split_fasta, reference_digest
from wmlstudio.klebsiella_evidence import (
    CAPSULE_LIMITATIONS,
    capsule_html,
    locus_st_assignments,
    locus_st_html,
    summarize_capsule,
    summarize_locus_sts,
    type_capsule_markers,
    type_virulence_loci,
)
from wmlstudio.sequence import file_sha256
from wmlstudio.virulence_evidence import apply_locus_sts, summarize_virulence_hits

SOURCE = random.Random(5150)
ALLELES = {'geneA': {'1': ''.join(SOURCE.choices('ACGT', k=400)), '2': ''.join(SOURCE.choices('ACGT', k=400))},
           'geneB': {'1': ''.join(SOURCE.choices('ACGT', k=400)), '2': ''.join(SOURCE.choices('ACGT', k=400))}}
CAPSULE = {'wzi': {'1': ''.join(SOURCE.choices('ACGT', k=420)), '2': ''.join(SOURCE.choices('ACGT', k=420))},
           'wzc': {'5': ''.join(SOURCE.choices('ACGT', k=130)), '6': ''.join(SOURCE.choices('ACGT', k=130))}}
PROFILES = [('1', ('1', '1'), 'tst 1; ICEKp0'), ('2', ('2', '2'), 'tst 2; plasmid')]


def klebsiella_panel(tmp_path, *, alleles=None, profiles=None, capsule=None, name='references'):
    """A format-2 snapshot carrying a locus profile table and a capsule marker panel."""
    root = tmp_path / name
    root.mkdir()
    manifest = {'format_version': 2, 'source_revision': 'synthetic-truth',
                'source_repository': 'synthetic-test-fixture',
                'sources': {'kleborate': {'repository': 'synthetic', 'revision': 'synthetic-truth',
                                          'license': 'GPL-3.0-or-later', 'license_file': 'source-LICENSE-kleborate'}},
                'species': [], 'virulence': {}, 'locus_profiles': {}, 'sccmec': {}, 'capsule': {}, 'files': []}
    if alleles is not None:
        (root / 'virulence/tst').mkdir(parents=True)
        manifest['virulence']['tst'] = []
        for gene, values in alleles.items():
            relative = f'virulence/tst/{gene}.fasta'
            (root / relative).write_text(''.join(f'>{gene}_{allele}\n{sequence}\n' for allele, sequence in values.items()))
            manifest['virulence']['tst'].append({'gene': gene, 'path': relative})
        if profiles is not None:
            header = ['ST', *alleles, 'tst_lineage']
            rows = ['\t'.join(header)]
            rows += ['\t'.join([st, *vector, lineage]) for st, vector, lineage in profiles]
            (root / 'virulence/tst/profiles.tsv').write_text('\n'.join(rows) + '\n')
            manifest['locus_profiles']['tst'] = {'module': 'klebsiella__tst', 'path': 'virulence/tst/profiles.tsv',
                                                 'st_field': 'ST', 'lineage_field': 'tst_lineage',
                                                 'source': 'kleborate'}
    if capsule is not None:
        (root / 'capsule').mkdir()
        entries = []
        for gene, values in capsule.items():
            relative = f'capsule/{gene}.fasta'
            (root / relative).write_text(''.join(f'>{gene}_{allele}\n{sequence}\n' for allele, sequence in values.items()))
            entries.append({'gene': gene, 'path': relative, 'allele_count': len(values)})
        manifest['capsule'] = {'source': 'kaptive', 'loci': entries, 'k_locus_typing': 'not_staged'}
    for path in sorted(root.rglob('*')):
        if path.is_file():
            relative = path.relative_to(root).as_posix()
            manifest['files'].append({'path': relative, 'bytes': path.stat().st_size, 'sha256': file_sha256(path)})
    manifest['reference_digest'] = reference_digest(manifest)
    (root / 'manifest.json').write_text(json.dumps(manifest))
    return root


def assembly_of(tmp_path, *sequences, name='isolate.fasta', spacer='N' * 60):
    path = tmp_path / name
    path.write_text('>contig_1\n' + spacer.join(sequences) + '\n')
    return path


def test_locus_st_is_assigned_only_from_a_complete_exact_allele_vector(tmp_path):
    root = klebsiella_panel(tmp_path, alleles=ALLELES, profiles=PROFILES)
    assembly = assembly_of(tmp_path, ALLELES['geneA']['1'], ALLELES['geneB']['1'])
    result = type_virulence_loci(assembly, root)
    assert result['status'] == 'completed'
    locus = result['loci']['tst']
    assert locus['locus_st'] == '1' and locus['status'] == 'complete'
    assert locus['alleles'] == {'geneA': '1', 'geneB': '1'}
    assert result['official_type'] is None and locus['phenotype'] == 'not_inferred'
    assert summarize_locus_sts(result) == 'tst 1'


def test_incomplete_or_novel_profile_yields_no_locus_st_and_retains_the_engine_status(tmp_path):
    root = klebsiella_panel(tmp_path, alleles=ALLELES, profiles=PROFILES)
    incomplete = type_virulence_loci(assembly_of(tmp_path, ALLELES['geneA']['1']), root)
    assert incomplete['loci']['tst']['status'] == 'incomplete'
    assert incomplete['loci']['tst']['locus_st'] is None and incomplete['loci']['tst']['lineage'] is None
    novel = type_virulence_loci(assembly_of(tmp_path, ALLELES['geneA']['1'], ALLELES['geneB']['2'],
                                            name='novel.fasta'), root)
    assert novel['loci']['tst']['status'] == 'novel_profile'
    assert novel['loci']['tst']['locus_st'] is None
    assert any('absent from this local profile table' in note for note in novel['loci']['tst']['notes'])
    assert summarize_locus_sts(novel) == 'no locus ST assigned'


def test_ambiguous_and_mixed_engine_statuses_propagate_unchanged(tmp_path):
    root = klebsiella_panel(tmp_path, alleles=ALLELES, profiles=PROFILES)
    assembly = assembly_of(tmp_path, ALLELES['geneA']['1'], ALLELES['geneA']['2'], ALLELES['geneB']['1'])
    result = type_virulence_loci(assembly, root)
    assert result['loci']['tst']['status'] == 'mixed' and result['loci']['tst']['locus_st'] is None
    assert result['status'] == 'ambiguous'


def test_lineage_is_reported_verbatim_from_the_profile_table_and_labelled_as_such(tmp_path):
    root = klebsiella_panel(tmp_path, alleles=ALLELES, profiles=PROFILES)
    result = type_virulence_loci(assembly_of(tmp_path, ALLELES['geneA']['2'], ALLELES['geneB']['2']), root)
    locus = result['loci']['tst']
    assert locus['locus_st'] == '2' and locus['lineage'] == 'tst 2; plasmid'
    assert locus['lineage_field'] == 'tst_lineage'
    assert locus['lineage_source'] == 'Kleborate profile table, reported verbatim; not inferred by WMLSTudio'
    assert 'not an inference by this software' in ' '.join(result['limitations'])
    assert 'tst 2; plasmid' in locus_st_html(locus and result)


def test_a_snapshot_without_profile_tables_reports_not_run_not_an_absent_locus_st(tmp_path):
    root = klebsiella_panel(tmp_path, alleles=ALLELES)
    result = type_virulence_loci(assembly_of(tmp_path, ALLELES['geneA']['1']), root)
    assert result['status'] == 'not_run' and 'no Kleborate locus-ST profile tables' in result['reason']
    assert summarize_locus_sts(result) == 'not_run'


def test_official_locus_st_stays_none_when_the_locus_st_module_did_not_run():
    loci = {'iuc': [{'gene': 'iucA'}, {'gene': 'iucB'}]}
    groups = summarize_virulence_hits([], loci, adequate_negative_assay=True)
    assert groups[0]['official_locus_st'] is None and groups[0]['official_locus_st_source'] is None
    assigned = summarize_virulence_hits([], loci, adequate_negative_assay=True,
                                        locus_sts={'iuc': {'locus_st': '3', 'source': 'exact-allele call'}})
    assert assigned[0]['official_locus_st'] == '3' and assigned[0]['official_locus_st_source'] == 'exact-allele call'


def test_a_locus_st_reaches_the_virulence_screen_only_when_one_was_actually_assigned(tmp_path):
    root = klebsiella_panel(tmp_path, alleles=ALLELES, profiles=PROFILES)
    complete = type_virulence_loci(assembly_of(tmp_path, ALLELES['geneA']['1'], ALLELES['geneB']['1']), root)
    assignments = locus_st_assignments(complete)
    assert assignments['tst']['locus_st'] == '1' and 'exact-allele' in assignments['tst']['source']
    incomplete = type_virulence_loci(assembly_of(tmp_path, ALLELES['geneA']['1'], name='partial.fasta'), root)
    assert locus_st_assignments(incomplete) == {}
    evidence = {'status': 'completed',
                'loci': [{'locus': 'tst', 'official_locus_st': None, 'official_locus_st_source': None},
                         {'locus': 'other', 'official_locus_st': None, 'official_locus_st_source': None}]}
    apply_locus_sts(evidence, assignments)
    assert evidence['loci'][0]['official_locus_st'] == '1'
    assert evidence['loci'][1]['official_locus_st'] is None
    # A screen that did not complete is never retro-fitted with a locus ST.
    unfinished = {'status': 'failed', 'loci': [{'locus': 'tst', 'official_locus_st': None}]}
    apply_locus_sts(unfinished, assignments)
    assert unfinished['loci'][0]['official_locus_st'] is None


def test_wzi_allele_is_reported_without_any_k_locus_inference(tmp_path):
    root = klebsiella_panel(tmp_path, capsule=CAPSULE)
    result = type_capsule_markers(assembly_of(tmp_path, CAPSULE['wzi']['2'], CAPSULE['wzc']['5']), root)
    assert result['wzi'] == {'allele': '2', 'status': 'exact', 'candidates': ['2'],
                             'reason': result['wzi']['reason']}
    assert result['wzc']['allele'] == '5'
    assert result['k_locus'] is None and result['k_locus_status'] == 'not_assayed'
    assert 'neither shipped nor applied' in ' '.join(result['limitations'])
    assert summarize_capsule(result) == 'wzi 2; wzc 5'
    rendered = html.unescape(capsule_html(result))
    assert 'not assayed' in rendered and 'K locus' in rendered
    for limitation in CAPSULE_LIMITATIONS:
        assert limitation in rendered


def test_capsule_scheme_never_acquires_a_sequence_type(tmp_path):
    root = klebsiella_panel(tmp_path, capsule=CAPSULE)
    result = type_capsule_markers(assembly_of(tmp_path, CAPSULE['wzi']['1'], CAPSULE['wzc']['6']), root)
    assert result['st'] is None and 'st' not in result['markers']
    assert result['scheme_status'] == 'profile_unavailable'
    assert any('without ST assignment' in note for note in result['notes'])
    assert json.dumps(result).count('"k_locus": null') == 1


def test_ambiguous_wzc_match_is_surfaced_not_tie_broken(tmp_path):
    shared = CAPSULE['wzc']['5']
    capsule = {'wzi': CAPSULE['wzi'], 'wzc': {'5': shared, '6': shared, '7': CAPSULE['wzc']['6']}}
    root = klebsiella_panel(tmp_path, capsule=capsule)
    result = type_capsule_markers(assembly_of(tmp_path, CAPSULE['wzi']['1'], shared), root)
    assert result['wzc']['allele'] is None and result['wzc']['status'] == 'ambiguous'
    assert result['wzc']['candidates'] == ['5', '6']
    assert result['status'] == 'ambiguous'
    assert summarize_capsule(result) == 'wzi 1; wzc ambiguous'
    assert 'more prone to several exact matches' in ' '.join(result['limitations'])


def test_a_snapshot_without_capsule_markers_reports_not_run(tmp_path):
    root = klebsiella_panel(tmp_path, alleles=ALLELES)
    result = type_capsule_markers(assembly_of(tmp_path, ALLELES['geneA']['1']), root)
    assert result['status'] == 'not_run' and 'wzi/wzc capsule marker panel' in result['reason']


def test_wzi_header_rewrite_preserves_sequence_bytes(tmp_path):
    from wmlstudio.sequence import SequenceReader
    stage = tmp_path / 'stage'
    (stage / 'sources/kaptive').mkdir(parents=True)
    records = [('1__wzi__1__1', 'ACGTACGTAC'), ('1__wzi__14__14', 'TTGGCCAATT'),
               ('2__wzc__941__603', 'GGGGCCCCAA')]
    (stage / 'sources/kaptive/wzi_wzc_db.fasta').write_text(
        ''.join(f'>{header}\n{sequence}\n' for header, sequence in records))
    written = _split_fasta(stage, 'sources/kaptive/wzi_wzc_db.fasta', 'capsule', _capsule_split_name, None)
    by_group = {group: relative for group, relative, _ in written}
    assert sorted(by_group) == ['wzc', 'wzi']
    with SequenceReader(stage / by_group['wzi'], None) as reader:
        staged = {record.identifier: record.sequence for record in reader}
    assert staged == {'wzi_1': 'ACGTACGTAC', 'wzi_14': 'TTGGCCAATT'}
    with pytest.raises(ValueError, match='Unexpected capsule marker header'):
        _capsule_split_name('1__wzi__1')


def test_capsule_markers_require_one_panel_directory(tmp_path):
    root = klebsiella_panel(tmp_path, capsule=CAPSULE)
    manifest = json.loads((root / 'manifest.json').read_text())
    manifest['capsule']['loci'][0]['path'] = 'virulence/tst/wzi.fasta'
    manifest['reference_digest'] = reference_digest(manifest)
    (root / 'manifest.json').write_text(json.dumps(manifest))
    with pytest.raises(ValueError, match='no source hash'):
        type_capsule_markers(assembly_of(tmp_path, CAPSULE['wzi']['1']), root)
