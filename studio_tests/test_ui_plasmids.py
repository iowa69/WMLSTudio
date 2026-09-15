"""The dedicated plasmid page: every table, and the question it says was never asked.

Nothing here runs a screen. Each isolate carries a stored plasmid block of the
exact shape :func:`wmlstudio.characterization.synthesize_accessory_evidence`
writes, and the reference store is a real manifest on disk, so the page is
exercised end to end against the evidence a real project holds.

The failure these tests exist to prevent is one picture: a reference store with no
plasmid set installed, an empty replicon table, and a reader concluding that a
cohort of carbapenemase producers carries no plasmid.
"""

import json

import pytest
from PySide6.QtWidgets import QLabel, QWidget

from wmlstudio import plasmid_evidence
from wmlstudio.ui_plasmids import TAB_TITLES, PlasmidPanel

DIGEST = "a" * 64


class FakeProject:
    """Just enough project for a page that only ever reads stored metadata."""

    def __init__(self, samples):
        self.rows = list(samples)

    def samples(self):
        return [dict(row) for row in self.rows]


class Host(QWidget):
    """A workspace window offering exactly the three things the contract names."""

    def __init__(self, samples, database=None):
        super().__init__()
        self.project = FakeProject(samples)
        self.messages = []
        self.roles = []
        self._database = database

    def launch_task(self, operation, role, completed=None, *, caption=None):
        self.roles.append(role)
        result = operation(lambda: False, lambda *args: None)
        if completed:
            completed(result)
        return True

    def notify(self, message):
        self.messages.append(str(message))

    def active_amr_database(self):
        return self._database


def store(tmp_path, *names):
    """A reference store holding exactly the named sets, as its manifest records them."""
    root = tmp_path / "hydra"
    root.mkdir(exist_ok=True)
    entries = {}
    for name in names:
        (root / name).mkdir(exist_ok=True)
        entries[name] = {"path": name, "kind": "nucl", "version": "2026-08-07.1"}
    (root / "manifest.json").write_text(json.dumps({"databases": entries}), encoding="utf-8")
    return str(root)


def isolate(sample_id, *, replicons=(), contigs=(), placements=(), associations=(),
            databases=("ncbi", "plasmidfinder"), status="completed", reason=""):
    """One isolate whose stored characterization holds a plasmid block."""
    block = {"status": status, "reason": reason,
             "replicons": [{"gene": gene} for gene in replicons],
             "source": {"databases": list(databases)},
             "contig_evidence": list(contigs),
             "determinant_placement": list(placements),
             "contig_associations": [{"replicon": pair[0], "marker": pair[1], "marker_type": "AMR"}
                                     for pair in associations]}
    return {"id": sample_id, "name": sample_id.upper(), "result": {"input_sha256": DIGEST},
            "metadata": {"characterization": {"input_sha256": DIGEST,
                                              "plasmid_hypotheses": block}}}


def contig(name="Contig_2_201.434_Circ", *, replicons=("IncFIB",), determinants=("blaKPC-2",)):
    return {"contig": name, "length_bp": 4372, "closure": "declared_circular",
            "declared_coverage": 201.434, "assembler_convention": "skesa",
            "declared_length_disagrees": False, "depth_ratio": 4.1,
            "signals": ["closure_claim", "depth_departure"],
            "support_state": "replicon_marker_plus_closure_and_depth",
            "replicons": [{"gene": gene} for gene in replicons],
            "determinants": [{"gene": gene} for gene in determinants]}


def placement(gene="blaKPC-2", *, where="Contig_2_201.434_Circ", replicons=("IncFIB",)):
    return {"gene": gene, "marker_type": "AMR", "contig": where,
            "placement": "co_located_with_replicon" if replicons else "no_replicon_on_this_contig",
            "replicons_on_contig": list(replicons),
            "interpretation": "Shares one assembled contig with a replicon marker: a co-location "
                              "hypothesis, not a plasmid-borne determinant."}


@pytest.fixture
def build(qtbot):
    """Build a panel on a bare host, the way a page-less caller would."""
    def make(samples, database=None):
        host = Host(samples, database)
        qtbot.addWidget(host)
        panel = PlasmidPanel(host)
        qtbot.addWidget(panel)
        return panel
    return make


def column(panel, key, heading):
    table = panel.tables[key]
    headings = [table.horizontalHeaderItem(index).text() for index in range(table.columnCount())]
    return headings.index(heading)


def texts(panel, key, heading):
    index = column(panel, key, heading)
    table = panel.tables[key]
    return [table.item(row, index).text() for row in range(table.rowCount())]


def test_a_panel_built_on_a_host_missing_the_contract_says_which_part_is_missing(qtbot):
    """A page mounted on the wrong window must name what it needs, not fail obscurely."""
    host = QWidget()
    qtbot.addWidget(host)
    with pytest.raises(TypeError) as error:
        PlasmidPanel(host)
    for name in PlasmidPanel.REQUIRED:
        assert name in str(error.value)


def test_the_page_offers_one_tab_for_every_table_the_payload_builds(build):
    """A table with no tab is a table nobody finds: that is how the plasmid work went missing."""
    panel = build([isolate("a", replicons=["IncFIB"])])
    assert list(TAB_TITLES) == list(plasmid_evidence.TABLE_NOTES)
    assert panel.tabs.count() == len(TAB_TITLES)
    assert [panel.tabs.tabText(index) for index in range(panel.tabs.count())] \
        == list(TAB_TITLES.values())


def test_no_plasmid_database_installed_is_never_shown_as_no_replicons_found(build, tmp_path):
    """The exact false negative: an unasked question rendering as a clean cohort."""
    panel = build([isolate("a", databases=["ncbi"])], database=store(tmp_path, "ncbi"))
    assert "No plasmid database is installed" in panel.database_state.text()
    assert texts(panel, "replicons_by_isolate", "Replicon screen") == [
        plasmid_evidence.SCREEN_WORDS["plasmid_database_not_installed"]]
    reason = texts(panel, "replicons_by_isolate", "What this row means")[0]
    assert "not evidence that these isolates carry no plasmid" in reason
    # The empty grid says the same thing where the missing rows would have been.
    assert "Install plasmidfinder" in panel.empties["replicon_matrix"].text()


def test_a_search_that_found_nothing_is_kept_apart_from_a_search_that_never_ran(build, tmp_path):
    """"Searched and found none" and "never looked" must not share a cell."""
    panel = build([isolate("a", databases=["ncbi", "plasmidfinder"]),
                   isolate("b", databases=["ncbi"])],
                  database=store(tmp_path, "ncbi", "plasmidfinder"))
    assert "plasmidfinder" in panel.database_state.text()
    assert texts(panel, "replicons_by_isolate", "Replicon screen") == [
        plasmid_evidence.SCREEN_WORDS["none_detected"],
        plasmid_evidence.SCREEN_WORDS["not_searched"]]


def test_an_unreadable_reference_store_is_unknown_rather_than_empty(build):
    """Only a store can say a set is absent; a page that never reached one says so."""
    panel = build([isolate("a", replicons=["IncFIB"])])
    assert "unknown" in panel.database_state.text()
    assert panel.payload["database"]["status"] == "unknown"


def test_an_isolate_with_no_plasmid_result_is_listed_rather_than_counted_as_carrying_nothing(build):
    panel = build([isolate("a", replicons=["IncFIB"]),
                   isolate("b", status="not_run", reason="No plasmid-reference assay was run.")])
    assert texts(panel, "isolates_not_assayed", "Isolate") == ["B"]
    assert texts(panel, "isolates_not_assayed", "Why") == ["No plasmid-reference assay was run."]
    assert panel.payload["denominator"] == 1 and panel.payload["isolate_count"] == 2
    assert "1 of 2 isolate(s) had a plasmid-marker assay" in panel.headline.text()


def test_the_replicon_grid_never_reads_an_unassayed_isolate_as_a_replicon_absence(build):
    """A blank in a presence grid is the same false negative in a different shape."""
    panel = build([isolate("a", replicons=["IncFIB"]),
                   isolate("b", status="not_run", reason="No plasmid-reference assay was run.")])
    assert texts(panel, "replicon_matrix", "IncFIB") == ["present", "not assayed"]


def test_the_per_isolate_contig_evidence_is_a_table_now_and_not_only_html(build):
    """It existed only inside a detail pane, where it could not be sorted or copied."""
    panel = build([isolate("a", replicons=["IncFIB"], contigs=[contig()],
                           placements=[placement()])])
    assert texts(panel, "contigs", "Contig") == ["Contig_2_201.434_Circ"]
    assert texts(panel, "contigs", "Closure (assembler's claim)") == ["declared circular"]
    assert texts(panel, "contigs", "Coverage vs backbone") == ["4.1×"]
    assert texts(panel, "contigs", "Determinants on the same contig") == ["blaKPC-2"]


def test_a_determinant_with_no_usable_coordinates_appears_rather_than_vanishing(build):
    """An omitted row reads as a determinant that was not found, which it is not."""
    unplaced = {"gene": "blaOXA-48", "marker_type": "AMR", "contig": None, "placement": "unplaced",
                "replicons_on_contig": [],
                "interpretation": "No usable contig coordinates were reported for this "
                                  "determinant, so nothing is claimed about where it sits."}
    panel = build([isolate("a", replicons=["IncFIB"], placements=[placement(), unplaced])])
    assert texts(panel, "determinant_colocation", "Determinant") == ["blaKPC-2", "blaOXA-48"]
    assert texts(panel, "determinant_colocation", "Contig")[1] == "no usable contig coordinates"


def test_every_table_prints_its_own_limit_beside_the_numbers(build):
    """These limits are severe; a reader who must open a help page will over-read the table."""
    panel = build([isolate("a", replicons=["IncFIB"], contigs=[contig()])])
    for index, key in enumerate(TAB_TITLES):
        page = panel.tabs.widget(index)
        printed = [item.text() for item in page.findChildren(QLabel)]
        assert plasmid_evidence.TABLE_NOTES[key] in printed
    assert "not a plasmid" in panel.limitations.text()
    assert "MOB-suite" in panel.gap.text()


def test_the_page_states_what_a_mob_suite_run_would_have_added(build):
    panel = build([isolate("a", replicons=["IncFIB"])])
    for missing in ("relaxase", "MPF", "oriT", "reconstruction", "cluster code"):
        assert missing in panel.gap.text()


def test_the_report_payload_is_the_same_tables_the_page_is_showing(build):
    """One payload for the page and the report, so neither can drift from the other."""
    panel = build([isolate("a", replicons=["IncFIB"], contigs=[contig()],
                           placements=[placement()], associations=[("IncFIB", "blaKPC-2")])])
    payload = panel.report_payload()
    assert set(payload["tables"]) == set(TAB_TITLES)
    for key, spec in payload["tables"].items():
        table = panel.tables[key]
        assert table.columnCount() == len(spec["columns"])
        assert table.rowCount() == len(spec["rows"])
    assert payload["tables"]["cohort_cooccurrence"]["rows"][0][:3] == ["IncFIB", "blaKPC-2", "AMR"]


def test_a_cohort_too_large_to_read_inline_is_built_on_the_background_runner(build):
    """The same function on both paths: a slower cohort must not become a different one."""
    panel = build([isolate(f"s{index:03d}", replicons=["IncFIB"]) for index in range(60)])
    assert panel.host.roles == ["plasmids"]
    assert panel.payload["isolate_count"] == 60 and panel.payload["denominator"] == 60


def test_an_empty_cohort_says_nothing_was_read_rather_than_nothing_was_found(build):
    panel = build([])
    assert "Choose isolates" in panel.headline.text()
    assert all(panel.tables[key].rowCount() == 0 for key in TAB_TITLES)
