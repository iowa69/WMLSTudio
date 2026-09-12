import json
import random
import shutil

import pytest

from wmlstudio.identification import (
    clear_identification_cache,
    identify_assembly,
    scheme_organism,
)
from wmlstudio.sequence import AnalysisCancelled, SequenceError
from wmlstudio.typing import load_scheme, reverse_complement


@pytest.fixture
def panels(tmp_path):
    clear_identification_cache()
    paths, sequences = [], []
    for seed, organism in ((1, "Staphylococcus epidermidis"), (2, "Pseudomonas aeruginosa")):
        path = tmp_path / f"panel_{seed}"
        path.mkdir()
        rng = random.Random(seed)
        alleles = ["".join(rng.choice("ACGT") for _ in range(120)) for _ in range(7)]
        for index, sequence in enumerate(alleles):
            variant = ("A" if sequence[0] != "A" else "C") + sequence[1:]
            (path / f"gene{index}.tfa").write_text(f">gene{index}_1\n{sequence}\n>gene{index}_2\n{variant}\n")
        (path / "profiles.tsv").write_text("ST\t" + "\t".join(f"gene{i}" for i in range(7)) + "\n91\t" + "\t".join(["1"] * 7) + "\n")
        (path / "scheme.json").write_text(json.dumps({"name": f"Scheme {seed}", "organism": organism}))
        paths.append(path)
        sequences.append(alleles)
    yield paths, sequences
    clear_identification_cache()


def write_assembly(path, sequences):
    path.write_text("".join(f">contig{i}\n{sequence}\n" for i, sequence in enumerate(sequences)))
    return path


def test_strong_unique_match_selects_panel_and_retains_provisional_label(tmp_path, panels):
    paths, sequences = panels
    assembly = write_assembly(tmp_path / "unknown.fa", [reverse_complement(s) for s in sequences[0]])
    result = identify_assembly(assembly, paths)
    assert result["identification_status"] == "assigned"
    assert result["organism"] == {"genus": "Staphylococcus", "species": "epidermidis"}
    assert result["best_scheme_path"] == str(paths[0])
    assert result["typing_result"]["st"] == "91"
    assert result["candidates"][0]["matched_loci"] == 7
    assert any("not species confirmation" in note for note in result["notes"])


def test_shared_seeds_without_complete_alleles_do_not_assign(tmp_path, panels):
    paths, sequences = panels
    changed = [sequence[:-1] + ("A" if sequence[-1] != "A" else "C") for sequence in sequences[0]]
    assembly = write_assembly(tmp_path / "unknown.fa", changed)
    result = identify_assembly(assembly, paths)
    assert result["candidates"][0]["seed_loci"] == 7
    assert result["candidates"][0]["matched_loci"] == 0
    assert result["identification_status"] == "insufficient"
    assert result["best_scheme_path"] is None
    assert "typing_result" not in result


def test_tied_distinct_panels_never_pick_alphabetically(tmp_path, panels):
    paths, sequences = panels
    duplicate = tmp_path / "another_panel"
    shutil.copytree(paths[0], duplicate)
    (duplicate / "scheme.json").write_text('{"name":"Another panel","organism":"Othergenus species"}')
    result = identify_assembly(write_assembly(tmp_path / "unknown.fa", sequences[0]), [paths[0], duplicate])
    assert result["identification_status"] == "ambiguous"
    assert result["best_scheme_path"] is None
    assert result["provisional_label"] == "Unknown"
    assert len(result["candidates"]) == 2


def test_identical_snapshot_copies_do_not_create_false_ties(tmp_path, panels):
    paths, sequences = panels
    duplicate = tmp_path / "copy"
    shutil.copytree(paths[0], duplicate)
    result = identify_assembly(write_assembly(tmp_path / "unknown.fa", sequences[0]), [paths[0], duplicate])
    assert result["identification_status"] == "assigned"
    assert len(result["candidates"]) == 1


def test_weak_and_mixed_evidence_remains_unassigned(tmp_path, panels):
    paths, sequences = panels
    weak = identify_assembly(write_assembly(tmp_path / "weak.fa", sequences[0][:3]), paths)
    assert weak["identification_status"] == "insufficient"
    variant = ("A" if sequences[0][0][0] != "A" else "C") + sequences[0][0][1:]
    mixed = identify_assembly(write_assembly(tmp_path / "mixed.fa", sequences[0] + [variant]), paths)
    assert mixed["identification_status"] == "ambiguous"
    assert mixed["best_scheme_path"] is None
    assert mixed["candidates"][0]["status"] == "mixed"


def test_partial_strong_assignment_does_not_fabricate_st(tmp_path, panels):
    paths, sequences = panels
    result = identify_assembly(write_assembly(tmp_path / "partial.fa", sequences[0][:6]), paths)
    assert result["identification_status"] == "assigned"
    assert result["typing_result"]["status"] == "incomplete"
    assert result["typing_result"]["st"] is None


def test_cache_invalidates_after_reference_change(tmp_path, panels):
    paths, sequences = panels
    assembly = write_assembly(tmp_path / "isolate.fa", sequences[0])
    first = identify_assembly(assembly, paths)
    old = paths[0] / "gene0.tfa"
    old.write_text(old.read_text().replace(sequences[0][0], "T" * 120))
    second = identify_assembly(assembly, paths)
    assert second["typing_result"]["scheme_digest"] != first["typing_result"]["scheme_digest"]
    assert second["typing_result"]["alleles"]["gene0"] is None


def test_group_catalog_label_does_not_invent_species(panels):
    paths, _ = panels
    (paths[0] / "scheme.json").write_text('{"name":"Escherichia panel","organism":"Escherichia spp."}')
    _, parts = scheme_organism(load_scheme(paths[0]))
    assert parts == {"genus": "Escherichia", "species": ""}


@pytest.mark.parametrize('label,expected', [
    ('Escherichia coli complex', {'genus': 'Escherichia', 'species': ''}),
    ('Klebsiella pneumoniae/variicola/quasipneumoniae', {'genus': 'Klebsiella', 'species': ''}),
    ('Campylobacter jejuni/coli', {'genus': 'Campylobacter', 'species': ''}),
    ('Candidatus Liberibacter asiaticus', {'genus': '', 'species': ''}),
    ('Unknown panel', {'genus': '', 'species': ''}),
])
def test_taxonomic_label_parser_does_not_split_groups_into_single_species(label, expected):
    from wmlstudio.reference_catalog import organism_parts
    assert organism_parts(label) == expected


def test_plasmid_panel_cannot_be_used_as_host_identification(tmp_path, panels):
    paths, sequences = panels
    (paths[0] / 'scheme.json').write_text('{"name":"Plasmid MLST","organism":"Escherichia coli"}')
    result = identify_assembly(write_assembly(tmp_path / 'plasmid.fa', sequences[0]), [paths[0]])
    assert result['identification_status'] == 'no_schemes'
    assert result['organism'] == {'genus': '', 'species': ''}
    assert 'Plasmid panels' in result['excluded_schemes'][0]['reason']


def test_discovery_rejects_reads_and_observes_cancel(tmp_path, panels):
    paths, _ = panels
    reads = tmp_path / "reads.fa"
    reads.write_text("@r\nACGT\n+\nIIII\n")
    with pytest.raises(SequenceError, match="FASTA assembly"):
        identify_assembly(reads, paths)
    with pytest.raises(AnalysisCancelled):
        identify_assembly(reads, paths, cancelled=lambda: True)


def test_ineligible_and_invalid_schemes_are_disclosed(tmp_path, panels):
    paths, sequences = panels
    large = tmp_path / "large"
    large.mkdir()
    for index in range(31):
        (large / f"locus{index}.fa").write_text(f">locus{index}_1\nACGT\n")
    result = identify_assembly(write_assembly(tmp_path / "a.fa", sequences[0]), [large, tmp_path / "missing"])
    assert result["identification_status"] == "no_schemes"
    assert len(result["excluded_schemes"]) == 2
    assert any("31" in entry["reason"] for entry in result["excluded_schemes"])
