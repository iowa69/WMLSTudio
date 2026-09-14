import copy

import pytest

from wmlstudio import snp_annotation
from wmlstudio.sequence import AnalysisCancelled, file_sha256
from wmlstudio.snp_annotation import annotate_snps
from wmlstudio.typing import reverse_complement


def reference_annotation(tmp_path, *, embedded=True, partial=False):
    coding = 'ATGGAATGGTAA'  # M E W STOP
    sequence = coding + 'CCCCCC' + reverse_complement(coding) + 'AAAAAA'
    reference = tmp_path / 'reference.fasta'
    reference.write_text('>chromosome\n' + sequence + '\n')
    annotation = tmp_path / 'reference.gff3'
    annotation.write_text('##gff-version 3\n'
                          'chromosome\ttruth\tCDS\t1\t12\t.\t+\t0\tID=plus;gene=positive;locus_tag=P1;product=Test%20protein' +
                          (';partial=true' if partial else '') + '\n'
                          'chromosome\ttruth\tCDS\t19\t30\t.\t-\t0\tID=minus;gene=negative\n' +
                          ('##FASTA\n>chromosome\n' + sequence + '\n' if embedded else ''))
    return reference, annotation


def variant(position, ref, alt):
    return {'contig': 'chromosome', 'pos': position, 'ref': ref, 'alt': alt}


def effect(result, index=0, alternate=0):
    return result['sites'][index]['alternatives'][alternate]['consequences'][0]


def test_exact_reference_plus_minus_synonymous_missense_and_stop_controls(tmp_path):
    reference, annotation = reference_annotation(tmp_path)
    inputs = [variant(6, 'A', 'G'), variant(5, 'A', 'C'), variant(8, 'G', 'A'), variant(10, 'T', 'C'),
              variant(1, 'A', 'C'), variant(25, 'T', 'C'), variant(26, 'T', 'G')]
    result = annotate_snps(reference, inputs, annotation_path=annotation)
    assert [effect(result, i)['effect'] for i in range(len(inputs))] == [
        'synonymous', 'missense', 'stop_gained', 'stop_lost', 'synonymous', 'synonymous', 'missense']
    assert effect(result, 1)['ref_codon'] == 'GAA' and effect(result, 1)['alt_codon'] == 'GCA'
    assert effect(result, 6)['ref_codon'] == 'GAA' and effect(result, 6)['alt_codon'] == 'GCA'
    assert effect(result, 6)['strand'] == '-' and effect(result, 6)['cds_position'] == 5
    assert effect(result, 0)['product'] == 'Test protein'
    assert result['reference_sha256'] == file_sha256(reference)
    assert result['protocol']['annotation']['sha256'] == file_sha256(annotation)
    assert result['protocol']['annotation']['binding'] == 'embedded_reference_sequence_verified'
    assert result['protocol']['recombination_masking'] == 'not_run'
    assert len(result['protocol_sha256']) == 64
    assert result['summary'] == {'sites': 7, 'included': 7, 'excluded': 0}
    assert inputs[0] == variant(6, 'A', 'G')  # Never mutate caller variants.


def test_start_loss_multiallelic_indel_and_reference_allele_are_distinct(tmp_path):
    reference, annotation = reference_annotation(tmp_path)
    result = annotate_snps(reference, [variant(2, 'T', ['C', 'A', 'T', '<DEL>']), variant(3, 'G', 'GA')], annotation_path=annotation)
    assert effect(result)['effect'] == 'start_lost'
    assert effect(result, alternate=1)['effect'] == 'start_lost'
    assert result['sites'][0]['alternatives'][2]['status'] == 'reference_allele'
    assert result['sites'][0]['alternatives'][3]['status'] == 'unsupported_non_SNP'
    assert result['sites'][0]['included'] is True  # Supported alternative retained separately.
    assert result['sites'][1]['included'] is False and result['sites'][1]['exclusion_reasons'] == ['not_a_supported_SNP']


def test_bed_zero_based_half_open_coding_selection_and_excluded_context_retained(tmp_path):
    reference, annotation = reference_annotation(tmp_path)
    mask = tmp_path / 'repeat.bed'
    mask.write_text('chromosome\t5\t6\tIS-associated-repeat\n')  # Exactly one-based base 6.
    result = annotate_snps(reference, [variant(5, 'A', 'C'), variant(6, 'A', 'G'), variant(13, 'C', 'T')],
                           annotation_path=annotation, mask_bed=mask, coding_only=True)
    assert [row['included'] for row in result['sites']] == [True, False, False]
    assert result['sites'][1]['exclusion_reasons'] == ['explicit_BED_mask']
    assert result['sites'][1]['coding_context'] is True  # Repeats and coding regions overlap.
    assert effect(result, 1)['effect'] == 'synonymous'  # Exclusion does not erase annotation.
    assert result['sites'][2]['exclusion_reasons'] == ['outside_selected_CDS_annotation']
    assert result['sites'][2]['alternatives'][0]['status'] == 'noncoding'
    assert result['protocol']['mask']['sha256'] == file_sha256(mask)
    assert result['protocol']['eligible_reference_intervals'] == [
        {'contig': 'chromosome', 'start': 0, 'end': 5},
        {'contig': 'chromosome', 'start': 6, 'end': 12},
        {'contig': 'chromosome', 'start': 18, 'end': 30}]
    assert result['protocol']['eligible_reference_bases'] == 23
    assert any('regulatory' in note for note in result['limitations'])


def test_partial_and_compound_cds_do_not_invent_complete_frame(tmp_path):
    reference, annotation = reference_annotation(tmp_path, partial=True)
    result = annotate_snps(reference, [variant(6, 'A', 'G')], annotation_path=annotation)
    assert effect(result)['effect'] == 'unresolved' and 'Partial CDS' in effect(result)['reason']
    annotation.write_text(annotation.read_text().replace('ID=minus', 'ID=plus'))
    result = annotate_snps(reference, [variant(25, 'T', 'C')], annotation_path=annotation)
    assert effect(result)['effect'] == 'unresolved' and 'Multipart' in effect(result)['reason']


def test_wrong_reference_binding_coordinates_and_format_fail_explicitly(tmp_path):
    reference, annotation = reference_annotation(tmp_path)
    with pytest.raises(ValueError, match='REF does not exactly match'):
        annotate_snps(reference, [variant(6, 'T', 'G')], annotation_path=annotation)
    annotation.write_text(annotation.read_text().replace('ATGGAATGGTAA', 'ATGGAATGGTAC'))
    with pytest.raises(ValueError, match='embedded reference sequences'):
        annotate_snps(reference, [], annotation_path=annotation)
    bad_format = tmp_path / 'annotation.gbff'
    bad_format.write_text('not-a-valid-GenBank-file')
    with pytest.raises(ValueError, match='GFF3 only'):
        annotate_snps(reference, [], annotation_path=bad_format)
    mask = tmp_path / 'bad.bed'
    mask.write_text('chromosome\t-1\t20\n')
    with pytest.raises(ValueError, match='outside the exact reference'):
        annotate_snps(reference, [], predict_coding=False, mask_bed=mask)


def test_unverified_plain_gff_cannot_silently_exclude_noncoding_sites(tmp_path):
    reference, annotation = reference_annotation(tmp_path, embedded=False)
    result = annotate_snps(reference, [variant(6, 'A', 'G')], annotation_path=annotation)
    assert effect(result)['confidence'] == 'coordinate_only_unverified'
    with pytest.raises(ValueError, match='Coding-only filtering requires'):
        annotate_snps(reference, [], annotation_path=annotation, coding_only=True)
    with pytest.raises(ValueError, match='Declared annotation reference SHA'):
        annotate_snps(reference, [], annotation_path=annotation, annotation_reference_sha256='0' * 64)
    result = annotate_snps(reference, [], annotation_path=annotation, coding_only=True,
                           annotation_reference_sha256=file_sha256(reference))
    assert result['protocol']['annotation']['binding'] == 'explicit_user_reference_hash_binding'


def test_no_annotation_is_not_noncoding_and_bounds_and_cancellation_are_explicit(tmp_path):
    reference, _ = reference_annotation(tmp_path)
    result = annotate_snps(reference, [variant(6, 'A', 'G')], predict_coding=False)
    assert result['sites'][0]['alternatives'][0]['status'] == 'annotation_not_run'
    with pytest.raises(ValueError, match='site bound'):
        annotate_snps(reference, [variant(6, 'A', 'G')] * 2, predict_coding=False, max_sites=1)
    with pytest.raises(AnalysisCancelled):
        annotate_snps('missing', [], cancelled=lambda: True)


def test_reference_change_during_prediction_does_not_create_accepted_evidence(tmp_path, monkeypatch):
    reference, _ = reference_annotation(tmp_path)
    def changing_prediction(path, **kwargs):
        path.write_text(path.read_text() + '\n')
        return [], {'tool': 'test'}
    monkeypatch.setattr(snp_annotation, 'predict_cds', changing_prediction)
    with pytest.raises(ValueError, match='changed during annotation'):
        annotate_snps(reference, [])


def test_predicted_cds_provenance_is_distinct_and_reference_bound(tmp_path, monkeypatch):
    reference, _ = reference_annotation(tmp_path)
    genes = [{'id': 'cds0', 'contig': 'chromosome', 'start': 1, 'end': 12, 'strand': '+',
              'sequence': 'ATGGAATGGTAA', 'cds_qc': {'genetic_code': 11, 'partial_begin': False, 'partial_end': False}}]
    monkeypatch.setattr(snp_annotation, 'predict_cds', lambda *args, **kwargs: (copy.deepcopy(genes), {'tool': 'pyrodigal', 'version': 'test'}))
    result = annotate_snps(reference, [variant(6, 'A', 'G')], coding_only=True)
    assert effect(result)['confidence'] == 'predicted_CDS'
    assert result['protocol']['annotation']['reference_sha256'] == file_sha256(reference)
    assert result['protocol']['annotation']['kind'] == 'predicted_Prodigal'


def test_actual_native_pyrodigal_prediction_runs_on_bounded_reference(tmp_path):
    reference, _ = reference_annotation(tmp_path)
    result = annotate_snps(reference, [], predict_coding=True)
    assert result['status'] == 'completed'
    assert result['protocol']['annotation']['provenance']['tool'] == 'pyrodigal'
    assert result['protocol']['annotation']['binding'] == 'exact_reference_prediction'
