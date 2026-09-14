import gzip
import hashlib
import json
import os
import random
import re
from pathlib import Path

import pytest

from wmlstudio import organism_panel, species_evidence
from wmlstudio.characterization_refs import validate_characterization_references
from wmlstudio.organism_panel import (
    LICENSE_NOTICE,
    PANEL_LIMITATIONS,
    SPECIES_PANEL,
    assembly_source,
    bundled_species_panel,
    installed_species_panel,
    panel_digest,
    panel_manifest,
    provision_species_panel,
    verify_pins,
)
from wmlstudio.sequence import AnalysisCancelled


def synthetic_source(root):
    """A local mirror of the NCBI tree, so provisioning is exercised without network."""
    pinned = {}
    for index, row in enumerate(SPECIES_PANEL):
        accession, assembly_dir = row[0], row[1]
        sequence = "".join(random.Random(index).choices("ACGT", k=4000))
        payload = gzip.compress(f">{accession}\n{sequence}\n".encode(), mtime=0)
        directory = root / assembly_source(assembly_dir).rstrip("/")
        directory.mkdir(parents=True, exist_ok=True)
        name = f"{assembly_dir}_genomic.fna.gz"
        (directory / name).write_bytes(payload)
        (directory / "md5checksums.txt").write_text(
            f"{hashlib.md5(payload).hexdigest()}  ./{name}\n"
            "0123456789abcdef0123456789abcdef  ./unrelated.txt\n")
        pinned[accession] = (len(payload), hashlib.sha256(payload).hexdigest(),
                             hashlib.sha256(gzip.decompress(payload)).hexdigest())
    return root, pinned


@pytest.fixture
def mirrored(tmp_path, monkeypatch):
    source, pinned = synthetic_source(tmp_path / "ncbi")
    rows = tuple((accession, assembly_dir, genus, species, strain, *pinned[accession])
                 for accession, assembly_dir, genus, species, strain, *_ in SPECIES_PANEL)
    monkeypatch.setattr(organism_panel, "SPECIES_PANEL", rows)
    return source


def test_every_pinned_row_is_a_complete_verifiable_record():
    # A literal count, so adding or losing a reference is a deliberate edit.
    assert len(SPECIES_PANEL) == 18
    assert len({row[0] for row in SPECIES_PANEL}) == len(SPECIES_PANEL)
    assert len({row[1] for row in SPECIES_PANEL}) == len(SPECIES_PANEL)
    for accession, assembly_dir, genus, species, strain, size, gz_sha, fasta_sha in SPECIES_PANEL:
        assert re.fullmatch(r"GCF_\d{9}\.\d+", accession), accession
        assert assembly_dir.startswith(accession + "_")
        assert genus[:1].isupper() and species and species.islower() and strain
        assert 0 < size < organism_panel.MAX_FILE_BYTES
        assert re.fullmatch(r"[0-9a-f]{64}", gz_sha) and re.fullmatch(r"[0-9a-f]{64}", fasta_sha)
        assert gz_sha != fasta_sha
    assert sum(row[5] for row in SPECIES_PANEL) < organism_panel.MAX_TOTAL_BYTES
    # A taxon may carry more than one reference where its lineages straddle the
    # species line, but never the same strain twice, and the exceptions are named
    # so a second reference elsewhere has to be a deliberate edit.
    assert len({(row[2], row[3], row[4]) for row in SPECIES_PANEL}) == len(SPECIES_PANEL)
    repeated = {(row[2], row[3]) for row in SPECIES_PANEL
                if sum(1 for other in SPECIES_PANEL if other[2:4] == row[2:4]) > 1}
    assert repeated == {("Listeria", "monocytogenes")}
    assert len({row[2] for row in SPECIES_PANEL}) >= 13


def test_the_ftp_layout_matches_the_published_accession_directory_shape():
    assert assembly_source("GCF_000016305.1_ASM1630v1", "md5checksums.txt") == (
        "GCF/000/016/305/GCF_000016305.1_ASM1630v1/md5checksums.txt")
    assert assembly_source("GCF_015732555.1_ASM1573255v1") == (
        "GCF/015/732/555/GCF_015732555.1_ASM1573255v1/")
    with pytest.raises(ValueError, match="Unexpected assembly directory"):
        assembly_source("not_an_assembly")


def test_the_manifest_states_the_panels_limits_and_never_claims_a_species_database():
    manifest = panel_manifest("pinned_https")
    assert manifest["license_notice"] == LICENSE_NOTICE and "RefSeq" in LICENSE_NOTICE
    assert list(PANEL_LIMITATIONS) == manifest["limitations"]
    assert ("A small reference set is a triage panel, not a representation of "
            "within-species diversity.") in manifest["limitations"]
    assert "This panel does not distinguish Escherichia coli from Shigella." in manifest["limitations"]
    assert all(entry["outgroup"] for entry in manifest["species"])
    assert all(entry["taxonomy_basis"].startswith("NCBI RefSeq assembly") for entry in manifest["species"])
    assert manifest["virulence"] == {}
    assert bundled_species_panel() is None


def test_a_provisioned_panel_validates_and_is_accepted_as_an_ani_database(tmp_path, mirrored):
    report = provision_species_panel(tmp_path / "installed", source_root=mirrored)
    path = Path(report["path"])
    manifest = validate_characterization_references(path)
    assert manifest["virulence"] == {} and len(manifest["species"]) == len(SPECIES_PANEL)
    assert manifest["reference_digest"] == panel_digest(manifest) == report["reference_digest"]
    assert path.name == "species-panel-" + manifest["reference_digest"][:20]
    assert installed_species_panel(tmp_path) is None
    assert installed_species_panel(tmp_path / "installed" / "..") is None
    database, loaded = species_evidence._database(path)
    assert loaded["reference_digest"] == manifest["reference_digest"]
    assert database is not None
    for entry in manifest["files"]:
        assert entry["uncompressed_sha256"] and entry["md5"]
        assert entry["source_url"].startswith(organism_panel.NCBI_BASE)


def test_installed_panels_are_discovered_under_the_conventional_data_folder(tmp_path, mirrored):
    data_root = tmp_path / "Data"
    assert installed_species_panel(data_root) is None
    report = provision_species_panel(data_root / organism_panel.PANEL_DIRECTORY, source_root=mirrored)
    assert installed_species_panel(data_root) == Path(report["path"])


def test_a_single_changed_byte_is_refused_and_nothing_is_published(tmp_path, mirrored):
    row = organism_panel.SPECIES_PANEL[0]
    target = mirrored / assembly_source(row[1], f"{row[1]}_genomic.fna.gz")
    target.write_bytes(gzip.compress(b">tampered\nACGT\n", mtime=0))
    root = tmp_path / "installed"
    with pytest.raises(ValueError, match="expected .* bytes"):
        provision_species_panel(root, source_root=mirrored)
    assert list(root.iterdir()) == []


def test_bytes_that_disagree_with_the_published_checksum_are_refused(tmp_path, mirrored):
    row = organism_panel.SPECIES_PANEL[0]
    directory = mirrored / assembly_source(row[1]).rstrip("/")
    (directory / "md5checksums.txt").write_text(
        f"0123456789abcdef0123456789abcdef  ./{row[1]}_genomic.fna.gz\n")
    root = tmp_path / "installed"
    with pytest.raises(ValueError, match="checksum NCBI publishes"):
        provision_species_panel(root, source_root=mirrored)
    assert list(root.iterdir()) == []


def test_a_sequence_that_does_not_match_its_pinned_hash_is_refused(tmp_path, mirrored, monkeypatch):
    rows = list(organism_panel.SPECIES_PANEL)
    rows[0] = (*rows[0][:7], "f" * 64)
    monkeypatch.setattr(organism_panel, "SPECIES_PANEL", tuple(rows))
    root = tmp_path / "installed"
    with pytest.raises(ValueError, match="decompressed sequence does not match"):
        provision_species_panel(root, source_root=mirrored)
    assert list(root.iterdir()) == []


def test_a_failed_or_cancelled_provision_publishes_nothing(tmp_path, mirrored, monkeypatch):
    root = tmp_path / "installed"
    calls = []

    def broken(row, target, cancelled=None, **options):
        calls.append(row[0])
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(b"partial network response")
        raise OSError("connection lost")

    monkeypatch.setattr(organism_panel, "_fetch", broken)
    with pytest.raises(OSError, match="connection lost"):
        provision_species_panel(root, source_root=mirrored)
    assert calls and list(root.iterdir()) == []
    with pytest.raises(AnalysisCancelled):
        provision_species_panel(root, source_root=mirrored, cancelled=lambda: True)
    assert list(root.iterdir()) == []


def test_the_panel_fingerprint_is_stable_and_a_different_one_is_never_overwritten(tmp_path, mirrored):
    root = tmp_path / "installed"
    first = provision_species_panel(root, source_root=mirrored)
    second = provision_species_panel(root, source_root=mirrored)
    assert first["reference_digest"] == second["reference_digest"]
    assert first["path"] == second["path"]
    assert [path.name for path in root.iterdir()] == [Path(first["path"]).name]
    manifest = json.loads((Path(first["path"]) / "manifest.json").read_text())
    manifest["species"][0]["species"] = "relabelled"
    manifest["reference_digest"] = panel_digest(manifest)
    (Path(first["path"]) / "manifest.json").write_text(json.dumps(manifest))
    with pytest.raises(ValueError, match="different fingerprint"):
        provision_species_panel(root, source_root=mirrored)


def test_verify_pins_reads_only_the_published_checksums(tmp_path, mirrored):
    report = verify_pins(source_root=mirrored)
    assert report["checked"] == len(SPECIES_PANEL) and not report["missing"]
    assert len(report["available"]) == len(SPECIES_PANEL)
    row = organism_panel.SPECIES_PANEL[0]
    (mirrored / assembly_source(row[1], "md5checksums.txt")).write_text("\n")
    degraded = verify_pins(source_root=mirrored)
    assert degraded["missing"] and degraded["missing"][0]["accession"] == row[0]


def test_two_panels_stay_cached_together_and_are_still_revalidated(tmp_path, mirrored, monkeypatch):
    species_evidence._CACHE.clear()
    first = Path(provision_species_panel(tmp_path / "a", source_root=mirrored)["path"])
    second = Path(provision_species_panel(tmp_path / "b", source_root=mirrored)["path"])
    # Two installs of the same pins share a digest, so give the second its own identity.
    manifest = json.loads((second / "manifest.json").read_text())
    manifest["source_revision"] = "second-panel"
    manifest["reference_digest"] = panel_digest(manifest)
    (second / "manifest.json").write_text(json.dumps(manifest))
    sketched = []
    original = species_evidence._sequences

    def counted(path, cancelled=None):
        sketched.append(Path(path).parent.parent)
        return original(path, cancelled)

    monkeypatch.setattr(species_evidence, "_sequences", counted)
    for _ in range(2):
        species_evidence._database(first)
        species_evidence._database(second)
    assert sorted({str(path) for path in sketched}) == sorted({str(first), str(second)})
    assert len(sketched) == 2 * len(SPECIES_PANEL)
    (second / "manifest.json").write_text(json.dumps({**manifest, "files": []}))
    with pytest.raises(ValueError):
        species_evidence._database(second)
    species_evidence._CACHE.clear()


@pytest.mark.realdata
@pytest.mark.skipif(os.environ.get("WMLSTUDIO_REALDATA") != "1",
                    reason="Set WMLSTUDIO_REALDATA=1 to download the pinned panel from NCBI.")
def test_the_real_pinned_panel_downloads_and_matches_every_hash(tmp_path):
    report = provision_species_panel(tmp_path / "installed")
    manifest = validate_characterization_references(Path(report["path"]))
    pinned = {row[0]: row for row in SPECIES_PANEL}
    for entry in manifest["files"]:
        accession = Path(entry["path"]).name.removesuffix(".fna.gz")
        assert entry["sha256"] == pinned[accession][6]
        assert entry["uncompressed_sha256"] == pinned[accession][7]
        assert entry["bytes"] == pinned[accession][5]


def test_listeria_carries_both_lineages_because_they_straddle_the_species_line():
    """One reference per taxon silently fails where a species is internally diverse.

    L. monocytogenes lineages I and II sit either side of the 95% ANI species
    cutoff: a lineage II genome measured 94.83% against the lineage I reference
    and was left unresolved, which is honest but useless to the user. Two
    references fix that without weakening the ANI gate, which must not move.
    """
    from wmlstudio.organism_panel import SPECIES_PANEL
    listeria = [row for row in SPECIES_PANEL if row[2] == "Listeria" and row[3] == "monocytogenes"]
    assert len(listeria) == 2, "both L. monocytogenes lineages must be represented"
    strains = {row[4] for row in listeria}
    assert any("4b" in s for s in strains) and any("1/2a" in s for s in strains)
    # Never a self-match against the practice cohort: that would prove nothing.
    from wmlstudio.practice_cohorts import COHORTS
    cohort = {row[0] for entry in COHORTS.values() for row in entry["genomes"]}
    assert not cohort & {row[0] for row in SPECIES_PANEL}
