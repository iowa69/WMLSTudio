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
    short_citation,
    suggested_threshold,
    typing_scale,
)


def context():
    return {'organism': 'Klebsiella pneumoniae', 'method': 'cgmlst',
            'scheme_key': 'cgmlst.org:kpneumoniae-2358', 'scheme_digest': 'a' * 64,
            'locus_count': 2358, 'caller': 'WMLSTudio',
            'missing_policy': 'Shared callable loci; 95% overlap'}


def kp():
    return next(e for e in guidance_for('Klebsiella pneumoniae')['entries'] if e['source_id'] == 'glasgow2025')


def test_catalog_has_explicit_coverage_gaps_and_sources():
    # A sanity bound on the surveillance work-list, not a scientific claim: it
    # must stay populated and hand-curated rather than becoming a species dump.
    assert 20 <= len(ORGANISMS) <= 40
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


def test_a_curated_caveat_is_a_slice_of_the_curated_limitations_not_new_wording():
    """A short caveat must be the authors' curated words, never a paraphrase."""
    for entry in catalog_entries():
        assert entry["caveat"], entry["id"]
        assert entry["limitations"].startswith(entry["caveat"].rstrip(".")), entry["id"]
        assert entry["bindable"] == bool(entry["scheme_key"] and entry["published_threshold"] is not None)


def test_most_recent_is_decided_inside_one_method_and_never_across_two():
    faecium = guidance_for("Enterococcus faecium")["entries"]
    published = [entry["published"] for entry in faecium]
    assert published == sorted(published, reverse=True), "entries are offered newest first"
    for method in ("cgmlst", "snp"):
        group = [entry for entry in faecium if entry["method"] == method]
        assert sum(entry["most_recent"] for entry in group) >= 1
    # The 2022 SKA SNP entry is the most recent of its own method even though a
    # 2025 cgMLST entry exists: the two measure different quantities.
    snp = next(entry for entry in faecium if entry["method"] == "snp")
    assert snp["most_recent"] and snp["published"] == "2022-01-26"
    newest_cg = next(entry for entry in faecium if entry["method"] == "cgmlst")
    assert newest_cg["source_id"] == "glasgow2025" and newest_cg["most_recent"]


def test_a_suggestion_is_never_an_applied_setting_and_prefers_the_most_recent_source():
    suggestion = suggested_threshold("Klebsiella pneumoniae", locus_count=2358)
    assert suggestion["applied"] is False and suggestion["auto_apply"] is False
    assert suggestion["status"] == "suggestion_requires_review"
    assert suggestion["suggestion"]["source_id"] == "glasgow2025"
    assert suggestion["suggestion"]["published_threshold"] == 15
    assert suggestion["suggestion"]["source"]["doi"] == "10.1128/jcm.00646-25"
    assert suggestion["scheme_match"]["matches"] is True
    assert "WMLSTudio has not applied it" in suggestion["notice"]


def test_a_published_number_is_not_rescaled_onto_a_different_target_set():
    suggestion = suggested_threshold("Klebsiella pneumoniae", locus_count=1200)
    assert suggestion["suggestion"]["published_threshold"] == 15, "the number itself is never adjusted"
    assert suggestion["scheme_match"]["matches"] is False
    assert "measured over 2358 targets and this comparison used 1200" in suggestion["scheme_match"]["reason"]
    assert "No scaling for a different target set exists." in suggestion["scheme_match"]["reason"]


def test_no_core_genome_cutoff_is_suggested_for_a_classical_scheme():
    suggestion = suggested_threshold("Klebsiella pneumoniae", locus_count=7)
    assert suggestion["status"] == "scale_mismatch"
    assert suggestion["suggestion"] is None and suggestion["alternatives"] == []
    assert "classical MLST over 7 loci" in suggestion["headline"]
    assert typing_scale(7)["kind"] == "mlst" and typing_scale(2358)["kind"] == "cgmlst"
    assert typing_scale(None)["kind"] == "unknown"


@pytest.mark.parametrize("organism,expected", [
    ("Enterococcus faecium", "do not agree on one number"),
    ("Listeria monocytogenes", "The same integer on a different scheme is not the same cutoff."),
])
def test_disagreeing_sources_are_shown_rather_than_silently_dropped(organism, expected):
    loci = {"Enterococcus faecium": 1423, "Listeria monocytogenes": 1701}[organism]
    suggestion = suggested_threshold(organism, locus_count=loci)
    assert expected in suggestion["disagreement"]
    assert suggestion["alternatives"], "the sources that were not chosen stay visible"


def test_the_newest_source_is_named_even_when_it_cannot_supply_a_number():
    difficile = suggested_threshold("Clostridioides difficile", locus_count=2270)
    assert difficile["suggestion"]["source_id"] == "bletz2018"
    assert difficile["newest"]["source_id"] == "siddall2025"
    assert "is not bound to a scheme this catalog can bind" in difficile["headline"]
    # A scheme with a deliberately uncurated number suggests nothing at all.
    aeruginosa = suggested_threshold("Pseudomonas aeruginosa", locus_count=3867)
    assert aeruginosa["status"] == "no_bindable_cutoff" and aeruginosa["suggestion"] is None
    absent = suggested_threshold("Nocardia farcinica", locus_count=2000)
    assert absent["status"] == "no_curated_transferable_cutoff"
    assert "evidence gap" in absent["headline"]


def test_a_recorded_decision_carries_the_citation_caveat_and_published_target_set():
    citation_only = record_decision(kp()["id"], {})
    assert citation_only["caveat"].startswith("Selection was based on prior clustering")
    assert citation_only["short_citation"] == "Glasgow et al. (2025)"
    assert citation_only["published_locus_count"] == 2358
    assert citation_only["published_unit"] == "allele differences"
    assert citation_only["approved_threshold"] is None
    assert short_citation(SOURCES["debeen2015"]) == "de Been et al. (2015)"
    assert short_citation(SOURCES["vanwalle2018"]) == "Van Walle et al. (2018)"


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
