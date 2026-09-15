"""A laboratory's own cutoff is its own cutoff, never a publication."""

import pytest

from wmlstudio.threshold_guidance import (
    catalog_entries,
    operational_cutoffs,
    record_decision,
    suggested_threshold,
)


def klebsiella():
    return operational_cutoffs('Klebsiella pneumoniae', locus_count=2358)[0]


def stub_window():
    """The little the guidance dialog asks of its window, and nothing more.

    The dialog only reads the current cohort and comparison key, so the test
    builds those instead of a whole main window: a change anywhere else in the
    application must not decide whether this claim is still true.
    """
    from PySide6.QtWidgets import QWidget

    class Project:
        @staticmethod
        def samples():
            return []

        @staticmethod
        def comparison_revision():
            return 0

    window = QWidget()
    window.project, window.cohort_ids, window.comparison_mode = Project(), set(), 'cgmlst'
    window._comparison_request = lambda: {}
    window._comparison_key = lambda request, revision: (repr(request), revision)
    return window


def test_the_klebsiella_operational_cutoff_is_five_on_the_full_2358_target_scheme():
    """The number the laboratory actually runs must be reachable, and bound.

    Without it the application can only offer 15, so a user whose rule is 5 has
    no way to see their own cutoff stated with its scheme and target count.
    """
    row = klebsiella()
    assert row['operational_threshold'] == 5
    assert row['unit'] == 'allele differences' and row['operator'] == '<='
    assert row['scheme_key'] == 'cgmlst.org:kpneumoniae-2358'
    assert row['locus_count'] == 2358 and row['method'] == 'cgmlst'


def test_the_operational_cutoff_claims_no_publication_and_carries_none():
    """Prevents an unsourced number from acquiring the authority of a citation.

    A fabricated or borrowed citation beside a local rule is the worst failure
    this catalogue can produce, so the row must have no source, no DOI and no
    published number at all, and must say in words where it came from.
    """
    row = klebsiella()
    assert row['published_threshold'] is None and row['bindable'] is False
    assert 'source_id' not in row and 'doi' not in row and 'citation' not in row
    assert row['declared_by'] and row['declared_on']
    assert 'found no publication establishing 5 allele differences' in row['provenance']
    assert 'not a published one' in row['notice']


def test_an_operational_cutoff_is_never_a_catalogue_publication_row():
    """Every consumer of a catalogue row prints it as evidence with a DOI.

    If a local number ever entered catalog_entries(), the report and the scheme
    binding would render it as "Published cutoff for this reference", which would
    be false.
    """
    ids = {entry['id'] for entry in catalog_entries()}
    assert klebsiella()['id'] not in ids
    for entry in catalog_entries():
        assert entry.get('evidence_class') != 'user_supplied_operational'


def test_an_operational_cutoff_cannot_be_recorded_as_a_citation():
    """A saved decision is read downstream as an adopted publication.

    Recording the local number that way would make the report say the threshold
    in force is a published cutoff, when no publication backs it.
    """
    with pytest.raises(ValueError, match='not publication evidence'):
        record_decision(klebsiella()['id'], {}, selected_threshold=5,
                        justification='Our laboratory applies five allele differences.',
                        protocol_reviewed=True, schema_reviewed=True, epi_reviewed=True)
    with pytest.raises(ValueError, match='not publication evidence'):
        record_decision(klebsiella()['id'], {})


@pytest.mark.parametrize('binding', [
    {'locus_count': 2357}, {'locus_count': 2390}, {'locus_count': 7}, {'locus_count': 'many'},
    {'scheme_key': 'pasteur:kpneumoniae-629'}, {'scheme_key': 'cgmlst.org:kpneumoniae-2365'},
])
def test_a_local_cutoff_declared_for_one_target_set_is_not_offered_for_another(binding):
    """Five differences over 2,358 targets says nothing about any other panel.

    A distance measured over a different target set is a different quantity, so
    the row must disappear rather than be rescaled or quietly reused.
    """
    assert operational_cutoffs('Klebsiella pneumoniae', **binding) == []


def test_a_local_cgmlst_cutoff_is_not_offered_beside_a_seven_locus_distance():
    """A classical MLST distance and a cgMLST cutoff share no scale.

    The scale guard must cover local cutoffs exactly as it covers published ones.
    """
    payload = suggested_threshold('Klebsiella pneumoniae', locus_count=7)
    assert payload['status'] == 'scale_mismatch'
    assert payload['operational'] == [] and payload['operational_notice'] == ''
    assert operational_cutoffs('Klebsiella pneumoniae', 'snp', locus_count=2358) == []
    assert operational_cutoffs('Escherichia coli') == [], 'no local cutoff was declared for this organism'


def test_the_suggestion_keeps_the_local_cutoff_apart_from_the_published_one():
    """Both numbers are shown, and neither is allowed to stand in for the other.

    Merging them would let 5 be read as published, or let 15 be read as the rule
    the laboratory actually runs.
    """
    payload = suggested_threshold('Klebsiella pneumoniae', locus_count=2358,
                                  scheme_key='cgmlst.org:kpneumoniae-2358')
    assert payload['suggestion']['published_threshold'] == 15
    assert payload['suggestion']['source_id'] == 'glasgow2025'
    assert [row['operational_threshold'] for row in payload['operational']] == [5]
    assert all(row.get('evidence_class') != 'user_supplied_operational'
               for row in [payload['suggestion'], *payload['alternatives']])
    assert payload['applied'] is False and payload['auto_apply'] is False


def test_the_local_cutoff_is_printed_beside_the_published_numbers_it_departs_from():
    """A local rule three times stricter than the published one must show that.

    Reading the notice on its own has to be enough to see that 5 is nobody's
    published number and that the reviewed source for the same scheme says 15.
    """
    notice = suggested_threshold('Klebsiella pneumoniae', locus_count=2358)['operational_notice']
    assert 'at most 5 allele differences on cgmlst.org:kpneumoniae-2358 over 2358 targets' in notice
    assert 'at most 15 allele differences from Glasgow et al. (2025)' in notice
    assert 'This local cutoff is not one of them.' in notice
    assert 'A stricter cutoff is not automatically the safer one.' in notice


def test_every_offerable_cgmlst_cutoff_names_both_a_scheme_and_its_target_count():
    """The binding rule is only enforceable when both halves are recorded.

    A cgMLST row with a number and a scheme key but no target count could be
    offered against any panel of that scheme's name, which is the exact transfer
    this catalogue exists to refuse.
    """
    for entry in catalog_entries():
        if entry['method'] == 'cgmlst' and entry['bindable']:
            assert entry['scheme_key'] and entry['locus_count'], entry['id']
    for row in operational_cutoffs(method=None):
        assert row['scheme_key'] and row['locus_count'], row['id']


def test_the_guidance_dialog_offers_the_local_cutoff_without_offering_to_cite_it(qtbot):
    """The dialog's Save path produces evidence; a local rule must not enter it.

    A user picking their own cutoff here must be told to set it, and must not
    walk away with a citation record that reports would print as published.
    """
    from wmlstudio.threshold_dialog import ThresholdGuideDialog
    window = stub_window()
    qtbot.addWidget(window)
    dialog = ThresholdGuideDialog(window)
    qtbot.addWidget(dialog)
    dialog.organism.setCurrentText('Klebsiella pneumoniae')
    titles = [dialog.entries.item(row).text() for row in range(dialog.entries.count())]
    local = next(row for row, title in enumerate(titles) if 'local cutoff' in title)
    assert 'not published' in titles[local]
    dialog.entries.setCurrentRow(local)
    shown = dialog.text.toPlainText()
    assert 'Your own operational cutoff: ≤ 5 allele differences' in shown
    assert 'cgmlst.org:kpneumoniae-2358' in shown and 'doi' not in shown.casefold()
    assert not dialog.apply_value.isEnabled()
    dialog.accept()
    assert dialog.evidence is None
    assert "your laboratory's own cutoff" in dialog.feedback.text()
    window.close()
