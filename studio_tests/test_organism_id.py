import json
import random
from pathlib import Path
from types import SimpleNamespace

import pytest

from wmlstudio import identification, species_evidence
from wmlstudio.characterization_refs import reference_digest
from wmlstudio.identification import clear_identification_cache, scheme_organism
from wmlstudio.organism_id import (
    AUTO_CONFIRM_FLOOR,
    AUTO_FILE_MIN,
    CONFIDENCE_ORDER,
    apply_policy,
    confidence_rank,
    evidence_sentences,
    identify_batch,
    identify_input,
    meets_policy,
    proposed_destination,
)
from wmlstudio.sequence import AnalysisCancelled, file_sha256
from wmlstudio.storage import QUARANTINE_ROOT, preview_target

SCHEME_ROOT = Path(identification.__file__).resolve().parent / "resources" / "schemes"
# Measured by replaying scheme_organism over the bundled snapshot. These five are
# out of clinical scope, and the last two carry labels organism_parts refuses to
# split on purpose ("Candidatus …", "Oral Streptococcus spp.").
UNLABELLED_SCHEMES = {"halobacteria", "liberibacter", "llactis_phage", "mamphoriforme",
                      "oralstrep", "streptothermophilus", "tvaginalis"}


def bundled_scheme_organisms():
    organisms = {}
    for directory in sorted(path for path in SCHEME_ROOT.iterdir() if path.is_dir()):
        files = sorted(path for path in directory.iterdir()
                       if path.name == "scheme.json" or path.name.endswith("_info.json"))
        metadata = json.loads(files[0].read_text(encoding="utf-8")) if files else {}
        organisms[directory.name] = scheme_organism(
            SimpleNamespace(metadata=metadata, name=directory.name))[1]
    return organisms


def reference_panel(root, genomes):
    """A minimal panel in the shape validate_characterization_references accepts."""
    root.mkdir(parents=True, exist_ok=True)
    manifest = {"format_version": 1, "source_revision": "synthetic-truth",
                "source_repository": "synthetic-test-fixture", "species": [], "virulence": {},
                "files": [], "limitations": ["Synthetic panel limitation sentence."]}
    for identifier, genus, species, sequence in genomes:
        relative = f"{identifier}.fasta"
        (root / relative).write_text(f">{identifier}\n{sequence}\n")
        manifest["species"].append({"id": identifier, "path": relative, "genus": genus,
                                    "species": species, "subspecies": "", "outgroup": True})
    for path in sorted(root.iterdir()):
        manifest["files"].append({"path": path.name, "bytes": path.stat().st_size,
                                  "sha256": file_sha256(path)})
    manifest["reference_digest"] = reference_digest(manifest)
    (root / "manifest.json").write_text(json.dumps(manifest))
    return root


def genome(seed, length=200_000):
    return "".join(random.Random(seed).choices("ACGT", k=length))


def assembly(path, sequence):
    path.write_text(f">contig\n{sequence}\n")
    return path


@pytest.fixture
def mlst_panel(tmp_path):
    clear_identification_cache()
    path = tmp_path / "panel_klebsiella"
    path.mkdir()
    rng = random.Random(11)
    alleles = ["".join(rng.choice("ACGT") for _ in range(120)) for _ in range(7)]
    for index, sequence in enumerate(alleles):
        variant = ("A" if sequence[0] != "A" else "C") + sequence[1:]
        (path / f"gene{index}.tfa").write_text(
            f">gene{index}_1\n{sequence}\n>gene{index}_2\n{variant}\n")
    (path / "profiles.tsv").write_text("ST\t" + "\t".join(f"gene{i}" for i in range(7))
                                       + "\n258\t" + "\t".join(["1"] * 7) + "\n")
    (path / "scheme.json").write_text(json.dumps({"name": "Synthetic Kp",
                                                  "organism": "Klebsiella pneumoniae"}))
    yield path, alleles
    clear_identification_cache()


def stub_species(monkeypatch, result):
    calls = []

    def fake(path, root, cancelled=None, progress=None, **options):
        calls.append((Path(path), Path(root)))
        return dict(result)

    monkeypatch.setattr(species_evidence, "identify_species", fake)
    return calls


def test_bundled_scheme_labels_never_regress_and_fewer_schemes_are_unlabelled():
    organisms = bundled_scheme_organisms()
    unlabelled = {name for name, parts in organisms.items() if not parts["genus"]}
    assert unlabelled == UNLABELLED_SCHEMES
    # The previously measured baseline was nineteen unlabelled schemes.
    assert len(unlabelled) < 19
    assert all(parts["genus"] or name in UNLABELLED_SCHEMES for name, parts in organisms.items())
    for name in ("ecoli", "ecoli_achtman_4", "senterica_achtman_2", "listeria_2",
                 "cdiphtheriae", "diphtheria_3", "kingella", "staphlugdunensis",
                 "mcatarrhalis_achtman_6"):
        assert organisms[name]["genus"] and organisms[name]["species"], name


def test_whole_genus_databases_stay_genus_only_rather_than_naming_a_species():
    organisms = bundled_scheme_organisms()
    for name, genus in (("klebsiella", "Klebsiella"), ("yersinia", "Yersinia"),
                        ("bordetella_3", "Bordetella")):
        assert organisms[name] == {"genus": genus, "species": ""}


def test_reads_are_routed_to_assembly_without_invoking_either_engine(tmp_path, monkeypatch):
    reads = tmp_path / "sample_R1.fastq"
    reads.write_text("@read1\nACGTACGT\n+\nIIIIIIII\n")

    def refuse(*args, **options):
        raise AssertionError("Reads must never reach an identification engine.")

    monkeypatch.setattr(species_evidence, "identify_species", refuse)
    monkeypatch.setattr(identification, "identify_assembly", refuse)
    verdict = identify_input(reads, species_panel_root=tmp_path, scheme_paths=[tmp_path])
    assert verdict["kind"] == "fastq" and verdict["basis"] == "none"
    assert verdict["confidence"] == "unresolved"
    assert verdict["quarantine_reason"] == "reads_not_assembled"
    assert verdict["proposed"] == {"genus": "", "species": ""}
    assert proposed_destination(verdict) == ("reads_not_assembled", "", "")


def test_an_mlst_panel_match_never_exceeds_panel_compatibility(tmp_path, mlst_panel):
    path, alleles = mlst_panel
    target = assembly(tmp_path / "isolate.fasta", "".join(alleles))
    verdict = identify_input(target, scheme_paths=[path])
    assert verdict["basis"] == "mlst_panel"
    assert verdict["confidence"] == "panel_compatibility"
    assert verdict["proposed"] == {"genus": "Klebsiella", "species": "pneumoniae"}
    assert verdict["detail"]["mlst_top"][0]["matched_loci"] == 7
    assert any("not species confirmation" in note for note in verdict["notes"])
    assert verdict["reason"] == verdict["notes"][0]
    assert confidence_rank(verdict["confidence"]) < confidence_rank(AUTO_FILE_MIN)
    assert confidence_rank(verdict["confidence"]) < confidence_rank(AUTO_CONFIRM_FLOOR)


def test_a_complete_panel_match_still_cannot_be_confirmed_by_any_policy(tmp_path, mlst_panel):
    path, alleles = mlst_panel
    verdict = identify_input(assembly(tmp_path / "isolate.fasta", "".join(alleles)),
                             scheme_paths=[path])
    for floor in CONFIDENCE_ORDER:
        policy = {"auto_confirm": True, "min_confidence": floor}
        assert meets_policy(verdict, policy) is False
        assert apply_policy(verdict, policy)["status"] == "proposed"
        assert proposed_destination(verdict, policy=policy)[0] == "awaiting_identification"


def test_a_reference_supported_call_is_offered_with_the_engines_own_sentences(tmp_path):
    target, other = genome(51), genome(52)
    panel = reference_panel(tmp_path / "panel", [("kp", "Klebsiella", "pneumoniae", target),
                                                 ("pa", "Pseudomonas", "aeruginosa", other)])
    query = assembly(tmp_path / "isolate.fasta", target)
    engine = species_evidence.identify_species(query, panel)
    verdict = identify_input(query, species_panel_root=panel)
    assert verdict["basis"] == "genomic_ani"
    assert verdict["confidence"] == "genomic_reference_supported"
    assert verdict["proposed"] == {"genus": "Klebsiella", "species": "pneumoniae"}
    assert verdict["reason"] == engine["reason"]
    assert verdict["limitations"] == engine["limitations"]
    assert "Synthetic panel limitation sentence." in verdict["limitations"]
    assert verdict["panel"]["reference_digest"] == engine["reference_digest"]
    assert len(verdict["detail"]["ani_top"]) <= 5
    assert verdict["reason"] in evidence_sentences(verdict)


def test_a_genus_level_reference_never_invents_a_species(tmp_path):
    target = genome(61)
    panel = reference_panel(tmp_path / "panel", [("cf", "Citrobacter", "", target),
                                                 ("pa", "Pseudomonas", "aeruginosa", genome(62))])
    verdict = identify_input(assembly(tmp_path / "isolate.fasta", target), species_panel_root=panel)
    assert verdict["confidence"] == "genus_only"
    assert verdict["proposed"] == {"genus": "Citrobacter", "species": ""}
    assert proposed_destination(dict(verdict, status="confirmed",
                                     accepted=verdict["proposed"])) == (None, "Citrobacter", "")


def test_a_complex_only_call_keeps_its_species_but_stays_below_the_filing_floor(tmp_path):
    target = genome(71)
    panel = reference_panel(tmp_path / "panel", [("ec", "Escherichia", "coli", target),
                                                 ("pa", "Pseudomonas", "aeruginosa", genome(72))])
    verdict = identify_input(assembly(tmp_path / "isolate.fasta", target), species_panel_root=panel)
    assert verdict["confidence"] == "complex_only"
    assert verdict["proposed"] == {"genus": "Escherichia", "species": "coli"}
    assert confidence_rank(verdict["confidence"]) < confidence_rank(AUTO_FILE_MIN)
    assert meets_policy(verdict, {"auto_confirm": True, "min_confidence": AUTO_FILE_MIN}) is False
    assert any("Shigella" in sentence for sentence in evidence_sentences(verdict))


def test_a_klebsiella_call_is_refined_against_the_focused_complex_panel(tmp_path):
    target = genome(81)
    broad = reference_panel(tmp_path / "broad", [("kleb", "Klebsiella", "", target),
                                                 ("pa", "Pseudomonas", "aeruginosa", genome(82))])
    focused = reference_panel(tmp_path / "focused", [("kp", "Klebsiella", "pneumoniae", target),
                                                     ("ec", "Escherichia", "coli", genome(83))])
    query = assembly(tmp_path / "isolate.fasta", target)
    verdict = identify_input(query, species_panel_root=broad, kpsc_panel_root=focused)
    assert verdict["basis"] == "genomic_ani_kpsc"
    assert verdict["confidence"] == "genomic_reference_supported"
    assert verdict["proposed"] == {"genus": "Klebsiella", "species": "pneumoniae"}
    assert verdict["panel"]["kind"] == "characterization_starter"


def test_a_non_klebsiella_call_is_not_re_queried_against_the_focused_panel(tmp_path):
    target = genome(91)
    broad = reference_panel(tmp_path / "broad", [("pa", "Pseudomonas", "aeruginosa", target),
                                                 ("kp", "Klebsiella", "pneumoniae", genome(92))])
    focused = reference_panel(tmp_path / "focused", [("kp", "Klebsiella", "pneumoniae", genome(93)),
                                                     ("ec", "Escherichia", "coli", genome(94))])
    queried = []
    original = species_evidence.identify_species

    def counted(path, root, *args, **options):
        queried.append(Path(root))
        return original(path, root, *args, **options)

    species_evidence.identify_species = counted
    try:
        verdict = identify_input(assembly(tmp_path / "isolate.fasta", target),
                                 species_panel_root=broad, kpsc_panel_root=focused)
    finally:
        species_evidence.identify_species = original
    assert queried == [broad]
    assert verdict["proposed"] == {"genus": "Pseudomonas", "species": "aeruginosa"}


def test_nothing_close_enough_is_reported_as_a_panel_limit_not_a_new_species(tmp_path, monkeypatch):
    stub_species(monkeypatch, {
        "status": "ambiguous", "genus": "", "species": "", "confidence": "unresolved",
        "nearest": None, "runner_up": None, "gap_ani": None, "hits": [],
        "thresholds": {"min_ani_pct": 95.0, "min_aligned_fraction": .65, "min_species_gap_ani_pct": .5},
        "reason": "Nearest adequate reference is below the configured species ANI threshold.",
        "limitations": ["Panel limitation."]})
    verdict = identify_input(assembly(tmp_path / "isolate.fasta", "ACGT" * 100),
                             species_panel_root=tmp_path)
    assert verdict["quarantine_reason"] == "not_in_reference_panel"
    assert verdict["confidence"] == "unresolved" and verdict["proposed"]["genus"] == ""
    assert verdict["reason"] == "Nearest adequate reference is below the configured species ANI threshold."
    assert verdict["limitations"] == ["Panel limitation."]


def test_competing_references_are_named_as_a_conflict_rather_than_a_weak_call(tmp_path, monkeypatch):
    stub_species(monkeypatch, {
        "status": "ambiguous", "genus": "", "species": "", "confidence": "unresolved",
        "nearest": {"reference_id": "a", "genus": "Klebsiella", "species": "pneumoniae", "ani": 99.0},
        "runner_up": {"reference_id": "b", "genus": "Klebsiella", "species": "variicola", "ani": 98.8},
        "gap_ani": 0.2, "hits": [],
        "thresholds": {"min_ani_pct": 95.0, "min_aligned_fraction": .65, "min_species_gap_ani_pct": .5},
        "reason": "Competing species references are too close; review possible hybrid, mixture or reference-panel limitations.",
        "limitations": []})
    verdict = identify_input(assembly(tmp_path / "isolate.fasta", "ACGT" * 100),
                             species_panel_root=tmp_path)
    assert verdict["quarantine_reason"] == "conflicting_evidence"
    assert verdict["margin_ani"] == 0.2
    assert verdict["runner_up"]["species"] == "variicola"
    assert "Competing species references" in verdict["reason"]


def test_an_unreadable_input_becomes_one_unresolved_verdict_not_a_failed_batch(tmp_path):
    broken = tmp_path / "broken.fasta"
    broken.write_text("this is not a sequence file\n")
    good = assembly(tmp_path / "good.fasta", "ACGTACGTACGT")
    verdicts = identify_batch([broken, good])
    assert len(verdicts) == 2
    assert verdicts[0]["confidence"] == "unresolved" and verdicts[0]["errors"]
    assert verdicts[0]["quarantine_reason"] == "awaiting_identification"
    assert [Path(verdict["input_path"]).name for verdict in verdicts] == ["broken.fasta", "good.fasta"]


def test_cancellation_is_never_swallowed(tmp_path, monkeypatch):
    target = assembly(tmp_path / "isolate.fasta", genome(101, 5_000))
    with pytest.raises(AnalysisCancelled):
        identify_input(target, species_panel_root=tmp_path, cancelled=lambda: True)
    with pytest.raises(AnalysisCancelled):
        identify_batch([target], cancelled=lambda: True)

    def cancel_inside(*args, **options):
        raise AnalysisCancelled("Analysis cancelled.")

    monkeypatch.setattr(species_evidence, "identify_species", cancel_inside)
    monkeypatch.setattr(identification, "identify_assembly", cancel_inside)
    with pytest.raises(AnalysisCancelled):
        identify_input(target, species_panel_root=tmp_path)
    with pytest.raises(AnalysisCancelled):
        identify_input(target, scheme_paths=[tmp_path])


def test_a_batch_returns_one_verdict_per_input_in_the_order_given(tmp_path, mlst_panel):
    path, alleles = mlst_panel
    inputs = [assembly(tmp_path / f"isolate{index}.fasta", "".join(alleles) + "A" * index)
              for index in range(3)]
    verdicts = identify_batch(inputs, scheme_paths=[path])
    assert [Path(verdict["input_path"]) for verdict in verdicts] == inputs
    assert all(verdict["input_sha256"] == file_sha256(path) for verdict, path in zip(verdicts, inputs))


def test_the_goes_to_column_and_the_storage_layer_agree_for_every_tier(tmp_path):
    root = tmp_path / "managed"
    original = tmp_path / "isolate.fasta"
    original.write_text(">c\nACGT\n")
    verdicts = [
        {"proposed": {"genus": "Klebsiella", "species": "pneumoniae"}, "confidence": "genomic_reference_supported"},
        {"proposed": {"genus": "Escherichia", "species": "coli"}, "confidence": "complex_only"},
        {"proposed": {"genus": "Enterobacter", "species": ""}, "confidence": "genus_only"},
        {"proposed": {"genus": "Klebsiella", "species": "pneumoniae"}, "confidence": "panel_compatibility"},
        {"proposed": {"genus": "", "species": ""}, "confidence": "unresolved"},
    ]
    for index, partial in enumerate(verdicts):
        verdict = {"status": "proposed", "quarantine_reason": None, "input_path": str(original),
                   "accepted": {"genus": "", "species": ""}, **partial}
        bucket, genus, species = proposed_destination(verdict)
        assert bucket == "awaiting_identification" and not genus
        accepted = dict(verdict, status="confirmed", accepted=dict(verdict["proposed"]))
        bucket, genus, species = proposed_destination(accepted)
        target = preview_target(root, f"sample{index}", original, genus, species, quarantine=bucket)
        parts = target.relative_to(root).parts
        if verdict["proposed"]["genus"]:
            assert bucket is None and parts[0] == verdict["proposed"]["genus"]
            assert parts[1] == (verdict["proposed"]["species"] or "Unknown_species")
        else:
            assert parts[0] == QUARANTINE_ROOT and parts[1] == "User_deferred"
        assert f"sample{index}" in parts


def test_auto_confirming_records_the_policy_as_the_confirming_party(tmp_path):
    verdict = {"status": "proposed", "confidence": "genomic_reference_supported",
               "proposed": {"genus": "Klebsiella", "species": "pneumoniae"},
               "accepted": {"genus": "", "species": ""}, "quarantine_reason": None,
               "input_path": str(tmp_path / "isolate.fasta")}
    assert apply_policy(verdict, {"auto_confirm": False})["status"] == "proposed"
    confirmed = apply_policy(verdict, {"auto_confirm": True, "min_confidence": AUTO_FILE_MIN})
    assert confirmed["status"] == "confirmed" and confirmed["confirmed_by"] == "auto_policy"
    assert proposed_destination(confirmed) == (None, "Klebsiella", "pneumoniae")


def test_an_unknown_confidence_token_is_ranked_as_the_weakest_never_the_strongest():
    assert confidence_rank("probably_novel") == 0
    assert confidence_rank(None) == 0
    assert meets_policy({"proposed": {"genus": "Klebsiella"}, "confidence": "probably_novel"},
                        {"auto_confirm": True, "min_confidence": "unresolved"}) is False
