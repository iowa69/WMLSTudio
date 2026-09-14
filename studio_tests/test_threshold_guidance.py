"""Cutoffs are scoped evidence, not species-wide predictions."""

import copy
from datetime import date

import pytest

from wmlstudio.threshold_guidance import (
    ORGANISMS,
    SOURCES,
    catalog_entries,
    guidance_for,
    record_decision,
    review_age_days,
)


def context():
    return {'organism': 'Klebsiella pneumoniae', 'method': 'cgmlst',
            'scheme_key': 'cgmlst.org:kpneumoniae-2358', 'scheme_digest': 'a' * 64,
            'locus_count': 2358, 'caller': 'WMLSTudio',
            'missing_policy': 'Shared callable loci; 95% overlap'}


def kp():
    return next(e for e in guidance_for('Klebsiella pneumoniae')['entries'] if e['source_id'] == 'glasgow2025')


def test_catalog_has_explicit_coverage_gaps_and_sources():
    assert 20 <= len(ORGANISMS) <= 30
    entries = catalog_entries()
    assert len({e['id'] for e in entries}) == len(entries)
    assert all(e['source_id'] in SOURCES and not e['auto_apply'] for e in entries)
    assert guidance_for('Klebsiella variicola')['entries'] == []
    assert guidance_for('Klebsiella')['entries'] == []
    assert guidance_for('unknown species')['status'] == 'no_curated_transferable_cutoff'
    assert review_age_days(date(2026, 9, 20)) == 8


def test_citation_only_never_changes_cutoff_even_without_profiles():
    result = record_decision(kp()['id'], {})
    assert result['approved_threshold'] is None
    assert result['published_threshold'] == 15
    assert result['doi'] == '10.1128/jcm.00646-25'
    assert result['status'] == 'citation_only_no_threshold_change'
    assert len(result['evidence_digest']) == 64


def test_local_adaptation_is_not_validated_and_input_is_immutable():
    before = context()
    saved = copy.deepcopy(before)
    result = record_decision(kp()['id'], before, selected_threshold=12,
        justification='Exploratory sensitivity analysis; confirm epidemiology independently.',
        protocol_reviewed=True, schema_reviewed=True, epi_reviewed=True)
    assert before == saved
    assert result['approved_threshold'] == 12
    assert result['departs_from_published_number']
    assert result['status'] == 'local_adaptation_requires_validation'


@pytest.mark.parametrize('change', [
    {'method': 'snp'}, {'method': 'wgmlst'}, {'method': 'st'},
    {'organism': 'Klebsiella quasipneumoniae'}, {'locus_count': 694},
    {'locus_count': 4891}, {'scheme_key': None}, {'scheme_digest': None},
    {'missing_policy': None}, {'caller': None},
])
def test_wrong_protocol_cannot_apply_published_value(change):
    with pytest.raises(ValueError):
        record_decision(kp()['id'], {**context(), **change}, selected_threshold=15,
            justification='These details were reviewed for this research analysis.',
            protocol_reviewed=True, schema_reviewed=True, epi_reviewed=True)


def test_incomplete_review_or_study_only_entry_cannot_apply():
    with pytest.raises(ValueError, match='Review'):
        record_decision(kp()['id'], context(), selected_threshold=15)
    local = next(e for e in guidance_for('Klebsiella pneumoniae')['entries'] if e['source_id'] == 'siddall2025')
    with pytest.raises(ValueError, match='exact published'):
        record_decision(local['id'], context(), selected_threshold=15)


def test_faecium_ska_and_cg_contexts_never_conflated():
    cg = guidance_for('Enterococcus faecium', 'cgmlst')['entries']
    assert {e['published_threshold'] for e in cg} >= {7, 20, 25}
    snp = guidance_for('Enterococcus faecium', 'snp')['entries']
    assert len(snp) == 1 and snp[0]['published_threshold'] == 7
    assert 'short-read' in snp[0]['scope']
    assert 'SKA2' in snp[0]['limitations']


def test_native_guidance_view_defaults_to_citation_only(qtbot, tmp_path):
    from wmlstudio.app import MainWindow
    from wmlstudio.threshold_dialog import ThresholdGuideDialog
    window = MainWindow(storage_root=tmp_path)
    qtbot.addWidget(window)
    dialog = ThresholdGuideDialog(window)
    qtbot.addWidget(dialog)
    dialog.organism.setCurrentText('Klebsiella pneumoniae')
    assert not dialog.apply_value.isChecked()
    assert '2358' in dialog.text.toPlainText()
    dialog.accept()
    assert dialog.evidence['approved_threshold'] is None
    window.close()
