import json
import random
from pathlib import Path

import pytest

from wmlstudio import identification
from wmlstudio.paths import scheme_locations
from wmlstudio.reference_index import (
    INDEX_DIRNAME,
    UNRESOLVED_NODE,
    clear_index_cache,
    clear_organism_index,
    filter_scheme_locations,
    organism_tree,
    panel_entries,
    scheme_entries,
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
    # Documents today's unfiltered behaviour: paths.scheme_locations returns every
    # subdirectory, so a caller must apply this filter until that one line changes.
    assert (root / INDEX_DIRNAME) in locations
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
