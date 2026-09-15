import json
import random
from pathlib import Path

import pytest

from wmlstudio import identification
from wmlstudio.paths import scheme_locations
from wmlstudio.reference_index import (
    CGMLST_LOCUS_THRESHOLD,
    INDEX_DIRNAME,
    UNRESOLVED_NODE,
    clear_index_cache,
    clear_organism_index,
    filter_scheme_locations,
    organism_tree,
    panel_entries,
    scheme_entries,
    scheme_library,
    write_organism_index,
)
from wmlstudio.sequence import file_signature

BUNDLED = Path(identification.__file__).resolve().parent / "resources" / "schemes"


def scheme(root, name, *, loci=7, metadata=None):
    path = root / name
    path.mkdir(parents=True, exist_ok=True)
    rng = random.Random(len(name))
    for index in range(loci):
        (path / f"locus{index}.tfa").write_text(
            f">locus{index}_1\n{''.join(rng.choice('ACGT') for _ in range(60))}\n")
    if metadata is not None:
        (path / "scheme.json").write_text(json.dumps(metadata))
    return path


@pytest.fixture(autouse=True)
def fresh_cache():
    clear_index_cache()
    yield
    clear_index_cache()


def test_a_scheme_reports_the_organism_it_records_and_how_it_was_read(tmp_path):
    declared = scheme(tmp_path, "declared", metadata={"name": "Kp", "organism": "Klebsiella pneumoniae"})
    slug = scheme(tmp_path, "slug", metadata={
        "name": "ecoli", "API": "https://rest.pubmlst.org/db/pubmlst_ecoli_seqdef/schemes/1"})
    bare = scheme(tmp_path, "bare")
    entries = {entry["id"]: entry for entry in scheme_entries([declared, slug, bare])}
    assert entries["declared"]["genus"] == "Klebsiella" and entries["declared"]["organism_basis"] == "metadata"
    assert entries["slug"] == {**entries["slug"], "genus": "Escherichia", "species": "coli",
                               "organism_basis": "api_slug"}
    assert entries["bare"]["organism_basis"] == "unresolved" and entries["bare"]["genus"] == ""
    assert all(entry["kind"] == "mlst" and entry["locus_count"] == 7 for entry in entries.values())


def test_a_large_or_declared_panel_is_classified_as_cgmlst(tmp_path):
    big = scheme(tmp_path, "big", loci=31, metadata={"name": "cg"})
    declared = scheme(tmp_path, "declared", loci=5, metadata={"name": "cg", "type": "cgMLST"})
    kinds = {entry["id"]: entry["kind"] for entry in scheme_entries([big, declared])}
    assert kinds == {"big": "cgmlst", "declared": "cgmlst"}


def test_entries_are_recomputed_when_a_reference_folder_changes(tmp_path):
    path = scheme(tmp_path, "panel", metadata={"name": "a", "organism": "Klebsiella pneumoniae"})
    first = scheme_entries([path])[0]
    (path / "scheme.json").write_text(json.dumps({"name": "a", "organism": "Listeria monocytogenes"}))
    second = scheme_entries([path])[0]
    assert first["genus"] == "Klebsiella" and second["genus"] == "Listeria"
    assert first["directory_digest"] != second["directory_digest"]


def test_the_derived_index_holds_pointers_and_never_a_single_sequence_byte(tmp_path):
    root = tmp_path / "schemes"
    paths = [scheme(root, "kp", metadata={"name": "kp", "organism": "Klebsiella pneumoniae"}),
             scheme(root, "ec", metadata={"name": "ec", "organism": "Escherichia coli"}),
             scheme(root, "kleb", metadata={"name": "kleb", "organism": "Klebsiella spp."}),
             scheme(root, "mystery")]
    before = {path: file_signature(path) for scheme_path in paths
              for path in scheme_path.rglob("*") if path.is_file()}
    index = write_organism_index(root, scheme_entries(paths))
    written = sorted(path for path in index.rglob("*") if path.is_file())
    assert {path.suffix for path in written} == {".json", ".txt"}
    assert sum(path.stat().st_size for path in written) < 1024 * 1024
    assert (index / "Klebsiella" / "pneumoniae" / "kp.json").is_file()
    assert (index / "Klebsiella" / "Unknown_species" / "kleb.json").is_file()
    assert (index / "Escherichia" / "coli" / "ec.json").is_file()
    assert (index / UNRESOLVED_NODE / "mystery.json").is_file()
    pointer = json.loads((index / "Klebsiella" / "pneumoniae" / "kp.json").read_text())
    assert pointer["scheme_path"] == str(paths[0]) and pointer["locus_count"] == 7
    assert pointer["organism"] == {"genus": "Klebsiella", "species": "pneumoniae"}
    assert not any(path.is_symlink() for path in index.rglob("*"))
    assert {path: file_signature(path) for scheme_path in paths
            for path in scheme_path.rglob("*") if path.is_file()} == before
    readme = (index / "README.txt").read_text(encoding="utf-8")
    assert "Derived index — safe to delete." in readme
    assert "not an independent verification" in readme


def test_rebuilding_the_index_is_idempotent_and_clearing_removes_only_the_index(tmp_path):
    root = tmp_path / "schemes"
    paths = [scheme(root, "kp", metadata={"name": "kp", "organism": "Klebsiella pneumoniae"})]
    entries = scheme_entries(paths)
    index = write_organism_index(root, entries)
    first = sorted(path.relative_to(index).as_posix() for path in index.rglob("*"))
    write_organism_index(root, entries)
    assert sorted(path.relative_to(index).as_posix() for path in index.rglob("*")) == first
    assert sorted(path.name for path in root.iterdir()) == [INDEX_DIRNAME, "kp"]
    clear_organism_index(root)
    assert not index.exists()
    assert sorted(path.name for path in root.iterdir()) == ["kp"]
    assert (root / "kp" / "locus0.tfa").is_file()


def test_a_derived_index_is_never_offered_as_an_installable_scheme(tmp_path):
    root = tmp_path / "schemes"
    scheme(root, "kp", metadata={"name": "kp", "organism": "Klebsiella pneumoniae"})
    write_organism_index(root, scheme_entries([root / "kp"]))
    (root / ".hidden").mkdir()
    locations = [path for path in scheme_locations(tmp_path) if path.parent == root]
    # That one line has now changed: scheme_locations filters at the source, so a
    # derived index cannot reach the typing code even through a caller that forgets
    # to filter. The filter stays idempotent for callers that already applied it.
    assert (root / INDEX_DIRNAME) not in locations
    assert [path.name for path in locations] == ["kp"]
    filtered = filter_scheme_locations(locations)
    assert [path.name for path in filtered] == ["kp"]
    assert all(not path.name.startswith(("_", ".")) for path in filter_scheme_locations(
        scheme_locations(tmp_path)))


def test_an_unfiltered_derived_index_is_still_harmless_to_the_discovery_engine(tmp_path):
    root = tmp_path / "schemes"
    paths = [scheme(root, "kp", metadata={"name": "kp", "organism": "Klebsiella pneumoniae"})]
    index = write_organism_index(root, scheme_entries(paths))
    assembly = tmp_path / "isolate.fasta"
    assembly.write_text(">c\n" + "ACGT" * 100 + "\n")
    identification.clear_identification_cache()
    result = identification.identify_assembly(assembly, [*paths, index])
    identification.clear_identification_cache()
    excluded = {Path(item["scheme_path"]).name for item in result["excluded_schemes"]}
    assert INDEX_DIRNAME in excluded
    assert result["identification_status"] in {"insufficient", "no_schemes", "ambiguous"}


def test_installed_ani_panels_contribute_one_row_per_covered_taxon(tmp_path):
    panel = tmp_path / "panels" / "characterization-abc"
    panel.mkdir(parents=True)
    (panel / "manifest.json").write_text(json.dumps({
        "format_version": 1, "reference_digest": "abc123", "panel_kind": "species_panel",
        "species": [{"id": "GCF_1", "genus": "Klebsiella", "species": "pneumoniae",
                     "taxonomy_basis": "NCBI RefSeq assembly GCF_1 (Klebsiella pneumoniae)"},
                    {"id": "GCF_2", "genus": "Citrobacter", "species": ""}], "virulence": {},
        "files": []}))
    entries = panel_entries(tmp_path / "panels")
    assert [entry["id"] for entry in entries] == ["GCF_2", "GCF_1"]
    assert all(entry["kind"] == "ani_panel" for entry in entries)
    assert all(entry["organism_basis"] == "panel_manifest" for entry in entries)
    assert panel_entries(panel) == entries
    assert panel_entries(tmp_path / "missing") == []


def test_unlabelled_references_keep_their_own_node_rather_than_being_dropped(tmp_path):
    paths = [scheme(tmp_path, "mystery"), scheme(tmp_path, "kp", metadata={"organism": "Klebsiella pneumoniae"})]
    tree = organism_tree(scheme_entries(paths))
    assert [entry["id"] for entry in tree[""][""]] == ["mystery"]
    assert [entry["id"] for entry in tree["Klebsiella"]["pneumoniae"]] == ["kp"]
    assert sum(len(rows) for children in tree.values() for rows in children.values()) == len(paths)


def test_indexing_the_bundled_snapshot_stays_cheap_and_complete(tmp_path):
    paths = sorted(path for path in BUNDLED.iterdir() if path.is_dir())
    entries = scheme_entries(paths)
    assert len(entries) == len(paths) >= 160
    index = write_organism_index(tmp_path / "schemes", entries)
    pointers = [path for path in index.rglob("*.json")]
    assert len(pointers) == len(paths)
    assert sum(path.stat().st_size for path in index.rglob("*") if path.is_file()) < 1024 * 1024
    unresolved = list((index / UNRESOLVED_NODE).glob("*.json"))
    assert 0 < len(unresolved) < 10
    assert (index / "Klebsiella").is_dir() and (index / "Escherichia" / "coli").is_dir()


def test_the_kind_of_a_scheme_is_reported_with_the_evidence_it_rests_on(tmp_path):
    classical = scheme(tmp_path, "classical", metadata={"name": "c", "type": "MLST"})
    counted = scheme(tmp_path, "counted", loci=40, metadata={"name": "b"})
    partial = scheme(tmp_path, "partial", loci=5, metadata={"name": "p", "type": "cgMLST"})
    rows = {entry["id"]: entry for entry in scheme_entries([classical, counted, partial])}
    assert rows["classical"]["kind"] == "mlst"
    assert rows["classical"]["kind_basis"] == "scheme metadata declares type 'mlst'"
    assert rows["counted"]["kind"] == "cgmlst"
    assert "40 targets" in rows["counted"]["kind_basis"]
    assert f"{CGMLST_LOCUS_THRESHOLD}-locus" in rows["counted"]["kind_basis"]
    assert rows["partial"]["kind"] == "cgmlst"
    assert "declares type 'cgmlst'" in rows["partial"]["kind_basis"]
    assert "partial snapshot" in rows["partial"]["kind_conflict"]


def test_a_folder_whose_kind_cannot_be_read_is_unknown_and_belongs_to_neither_tab(tmp_path):
    empty = tmp_path / "empty"
    empty.mkdir()
    (empty / "README.txt").write_text("this slot is not filled yet")
    row = scheme_entries([empty])[0]
    assert row["kind"] == "unknown"
    assert row["kind_basis"] == "the folder holds no allele FASTA files"
    library = scheme_library([empty])
    assert library["mlst"] == [] and library["cgmlst"] == []
    assert [entry["id"] for entry in library["unknown"]] == ["empty"]


def test_a_contradiction_between_size_and_declared_type_is_shown_not_resolved_silently(tmp_path):
    # A two-thousand-target scheme that calls itself MLST is not offered where a
    # seven-locus scheme is expected: the measured target count decides, and the
    # contradiction is reported so the metadata can be fixed.
    mislabelled = scheme(tmp_path, "mislabelled", loci=45, metadata={"name": "m", "type": "MLST"})
    row = scheme_entries([mislabelled])[0]
    assert row["kind"] == "cgmlst"
    assert "declares type 'MLST'" in row["kind_conflict"]
    assert "the installed target count decides" in row["kind_conflict"]


def test_each_library_can_be_asked_for_on_its_own(tmp_path):
    small = scheme(tmp_path, "small", metadata={"organism": "Klebsiella pneumoniae"})
    big = scheme(tmp_path, "big", loci=35, metadata={"organism": "Listeria monocytogenes"})
    paths = [small, big]
    assert [row["id"] for row in scheme_entries(paths, kind="mlst")] == ["small"]
    assert [row["id"] for row in scheme_entries(paths, kind="cgmlst")] == ["big"]
    assert set(scheme_library(paths)) == {"mlst", "cgmlst", "unknown"}
    with pytest.raises(ValueError, match="scheme kind"):
        scheme_entries(paths, kind="cgmlst_maybe")


def test_a_row_carries_a_readable_title_instead_of_a_folder_name(tmp_path):
    # The reported bug rendered a downloaded scheme as "cgmlst org kpneumoniae
    # abcdef0123456789". A row now says what a microbiologist would say.
    folder = scheme(tmp_path, "cgmlst_org_Kpneumoniae_complex_abcdef0123456789", loci=40,
                    metadata={"name": "Klebsiella pneumoniae sensu lato cgMLST",
                              "organism": "Klebsiella pneumoniae sensu lato", "type": "cgMLST",
                              "source": "cgMLST.org", "last_updated": "2026-09-14",
                              "API": "https://www.cgmlst.org/ncs/schema/Kpneumoniae_complex/"})
    row = scheme_entries([folder])[0]
    assert row["title"] == ("Klebsiella pneumoniae sensu lato · cgMLST · 40 targets · "
                            "cgMLST.org · updated 2026-09-14")
    assert row["organism_label"] == "Klebsiella pneumoniae sensu lato"
    assert row["provider"] == "cgMLST.org"
    assert row["version"] == "2026-09-14"
    assert "abcdef0123456789" not in row["title"]
    # The unit is named, because 40 targets and 7 loci are different quantities.
    classical = scheme(tmp_path, "kp_mlst", metadata={"organism": "Klebsiella pneumoniae",
                                                      "description": "MLST", "source": "pubmlst"})
    assert scheme_entries([classical])[0]["title"] == \
        "Klebsiella pneumoniae · MLST · 7 loci · PubMLST"


def test_a_scheme_that_records_nothing_still_gets_a_title_a_person_can_read(tmp_path):
    bare = scheme(tmp_path, "practice_7")
    row = scheme_entries([bare])[0]
    assert row["title"] == "practice_7 · 7 loci"
    assert row["kind"] == "mlst"
    assert row["provider"] == ""


def test_a_scheme_says_whether_its_target_set_is_core_or_accessory(tmp_path):
    core = scheme(tmp_path, "core", loci=40, metadata={"name": "c", "type": "cgMLST"})
    accessory = scheme(tmp_path, "accessory", loci=60, metadata={"name": "a", "type": "wgMLST"})
    unstated = scheme(tmp_path, "unstated", loci=40, metadata={"name": "u"})
    classical = scheme(tmp_path, "classical", metadata={"name": "m", "type": "MLST"})
    rows = {entry["id"]: entry["target_set"]
            for entry in scheme_entries([core, accessory, unstated, classical])}
    assert rows == {"core": "core", "accessory": "accessory", "unstated": "", "classical": ""}


def test_the_derived_index_pointer_records_the_kind_and_the_basis_for_it(tmp_path):
    root = tmp_path / "schemes"
    path = scheme(root, "kp", loci=40, metadata={"organism": "Klebsiella pneumoniae",
                                                 "type": "cgMLST", "name": "kp"})
    index = write_organism_index(root, scheme_entries([path]))
    pointer = json.loads((index / "Klebsiella" / "pneumoniae" / "kp.json").read_text())
    assert pointer["kind"] == "cgmlst"
    assert "40 targets" in pointer["kind_basis"]
    assert pointer["target_set"] == "core"
    assert "Klebsiella pneumoniae" in pointer["title"]
