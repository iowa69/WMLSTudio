import ast
import importlib.util
import io
import json
import tarfile
from pathlib import Path

import pytest

from wmlstudio import organism_panel, practice_cohorts, threshold_guidance
from wmlstudio.characterization_refs import SOURCES, reference_digest
from wmlstudio.sequence import file_sha256

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "studio_scripts/stage_schemes.py"
SPEC = importlib.util.spec_from_file_location("stage_schemes", SCRIPT)
stage = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(stage)
RECIPE = ROOT / "studio_packaging/wmlstudio.spec"
NOTICES = ROOT / "studio_packaging/THIRD_PARTY_NOTICES.md"


def panels_module():
    spec = importlib.util.spec_from_file_location(
        "stage_reference_panels", ROOT / "studio_packaging/stage_reference_panels.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def staged_panel(tmp_path, *, format_version=1, with_modules=False):
    """A minimal valid characterization snapshot, shaped like the real staging output."""
    root = tmp_path / "references"
    root.mkdir()
    manifest = {"format_version": format_version, "source_revision": "synthetic-truth",
                "species": [], "virulence": {}, "files": []}
    if format_version >= 2:
        manifest.update(sources={"synthetic": {"repository": "synthetic-test-fixture",
                                               "revision": "synthetic-truth", "license": "synthetic",
                                               "license_file": "source-LICENSE-synthetic"}},
                        locus_profiles={}, sccmec={}, capsule={})
    if with_modules:
        (root / "IVa.fasta").write_text(">IVa__AB063172__IV\nACGTACGTAC\n", encoding="utf-8")
        (root / "wzi.fasta").write_text(">wzi_1\nACGT\n", encoding="utf-8")
        (root / "ybt.tsv").write_text("ST\tybtS\tlineage_ICE\n1\t1\tybt 1\n", encoding="utf-8")
        manifest["locus_profiles"] = {"ybt": {"path": "ybt.tsv", "st_field": "ST"}}
        manifest["sccmec"] = {"regions": [{"gene": "IVa", "path": "IVa.fasta"}], "targets": [],
                              "rules": {"types": [{"name": "IV"}]}, "not_assayed": ["mecC"]}
        manifest["capsule"] = {"loci": [{"gene": "wzi", "path": "wzi.fasta", "allele_count": 1}]}
    for path in sorted(root.iterdir()):
        manifest["files"].append({"path": path.name, "bytes": path.stat().st_size,
                                  "sha256": file_sha256(path)})
    manifest["reference_digest"] = reference_digest(manifest)
    (root / "manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
    return root


def statically_imported_names():
    """Every module name reachable by an import statement anywhere in the package."""
    names = set()
    for path in (ROOT / "src/wmlstudio").rglob("*.py"):
        for node in ast.walk(ast.parse(path.read_text(encoding="utf-8"))):
            if isinstance(node, ast.Import):
                names.update(alias.name.rsplit(".", 1)[-1] for alias in node.names)
            elif isinstance(node, ast.ImportFrom):
                names.update(alias.name for alias in node.names)
                if node.module:
                    names.add(node.module.rsplit(".", 1)[-1])
    return names


def test_frozen_pyrodigal_explicitly_collects_namespace_implementations():
    from PyInstaller.utils.hooks import collect_submodules
    modules = collect_submodules("pyrodigal.impl")
    assert "pyrodigal.impl" in modules and "pyrodigal.impl.generic" in modules
    recipe = RECIPE.read_text()
    assert 'collect_submodules("pyrodigal.impl")' in recipe


def test_frozen_bundle_declares_the_assay_modules_the_registry_imports_by_name():
    """The registry reaches its assays through a computed __import__.

    That compiles to no IMPORT_NAME opcode and nothing else in the package
    references them, so the packaging scan cannot find them. Every
    characterization run calls module_tasks, so a bundle without them would fail
    on its first characterization rather than only when a module is selected.
    """
    modules = panels_module().assay_module_imports()
    assert modules, "the registry reported no assay modules to declare"
    reachable = statically_imported_names()
    assert not [name for name in modules if name.rsplit(".", 1)[-1] in reachable], (
        "an assay module became statically importable; confirm the packaging scan "
        "now finds it before relaxing this gate")
    recipe = RECIPE.read_text()
    assert "*assay_module_imports()" in recipe and "from stage_reference_panels import" in recipe


def test_packaged_documentation_includes_the_new_reference_guides():
    recipe = RECIPE.read_text()
    for name in ("ORGANISM_MODULES.md", "THRESHOLDS.md", "TEST_DATASETS.md"):
        assert (ROOT / "docs" / name).is_file()
        assert f'docs/{name}"), "docs"' in recipe


def test_notices_name_every_pinned_reference_source_its_revision_and_its_licence():
    notices = NOTICES.read_text(encoding="utf-8")
    for source in SOURCES.values():
        assert source.repository in notices, f"{source.key} repository is unattributed"
        assert source.revision in notices, f"{source.key} pinned revision is unrecorded"
        assert source.license_id in notices, f"{source.key} licence identifier is unrecorded"


def test_notices_carry_the_ncbi_statement_for_data_downloaded_on_request():
    notices = NOTICES.read_text(encoding="utf-8")
    quoted = " ".join(practice_cohorts.LICENSE_NOTICE.split())
    assert quoted in " ".join(notices.replace(">", " ").split())
    assert organism_panel.PANEL_REVISION in notices
    for name in practice_cohorts.cohort_names():
        assert practice_cohorts.cohort_digest(name) in notices, f"{name} pin is unrecorded"


def test_a_format_one_snapshot_is_reported_as_absent_with_the_reason_not_as_a_negative(tmp_path):
    summary = panels_module().panel_summary(staged_panel(tmp_path, format_version=1))
    assert summary["organism_modules"] == "absent"
    assert "locus_profiles, sccmec, capsule" in summary["organism_modules_reason"]
    assert "never a negative finding" in summary["organism_modules_reason"]
    assert summary["sccmec_targets"] == 0 and summary["capsule_loci"] == {}


def test_a_declared_but_empty_module_section_is_absent_not_staged(tmp_path):
    summary = panels_module().panel_summary(staged_panel(tmp_path, format_version=2))
    assert summary["organism_modules"] == "absent"
    assert summary["sources"] == {"synthetic": "synthetic-truth (synthetic)"}


def test_a_staged_panel_reports_its_sections_and_locates_its_cassette_reference(tmp_path):
    panels = panels_module()
    root = staged_panel(tmp_path, format_version=2, with_modules=True)
    summary = panels.panel_summary(root)
    assert summary["organism_modules"] == "staged" and summary["organism_modules_reason"] == ""
    assert summary["sccmec_regions"] == 1 and summary["sccmec_types"] == ["IV"]
    assert summary["capsule_loci"] == {"wzi": 1} and summary["locus_profiles"] == ["ybt"]
    assert summary["sccmec_not_assayed"] == ["mecC"]
    assert panels.sccmec_region_reference(root) == root / "IVa.fasta"
    assert panels.sccmec_region_reference(root, "nothing-like-this") is None


def test_a_changed_reference_byte_fails_the_packaging_gate(tmp_path):
    root = staged_panel(tmp_path, format_version=2, with_modules=True)
    (root / "wzi.fasta").write_text(">wzi_1\nTGCA\n", encoding="utf-8")
    with pytest.raises(ValueError, match="changed or is missing"):
        panels_module().panel_summary(root)


@pytest.mark.parametrize("kind", ["cohort-directory", "cohort-genome", "species-panel"])
def test_packaging_refuses_sequence_data_the_user_is_meant_to_download(tmp_path, kind):
    panels = panels_module()
    bundle = tmp_path / "resources"
    (bundle / "schemes").mkdir(parents=True)
    (bundle / "schemes/keep.txt").write_text("a legitimately bundled file", encoding="utf-8")
    assert panels.verify_bundle_payload([bundle])["unbundled_payload"] == []
    if kind == "cohort-directory":
        (bundle / practice_cohorts.default_destination(Path("."), "kpneumoniae-10").parent.name).mkdir()
    elif kind == "cohort-genome":
        name = practice_cohorts.genome_filename(practice_cohorts.cohort_genomes("kpneumoniae-10")[0])
        (bundle / "schemes" / name).write_bytes(b"not really a genome")
    else:
        (bundle / f"{organism_panel.PANEL_PREFIX}0123456789ab").mkdir()
    with pytest.raises(SystemExit, match="never travel in the portable build"):
        panels.verify_bundle_payload([bundle])
    assert (bundle / "schemes/keep.txt").is_file()


def frozen_module():
    spec = importlib.util.spec_from_file_location(
        "check_frozen", ROOT / "studio_packaging/check_frozen.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def frozen_bundle(tmp_path, **panel):
    bundle = tmp_path / "WMLSTudio"
    destination = bundle / "_internal/wmlstudio/resources/characterization/starter"
    destination.parent.mkdir(parents=True)
    staged_panel(tmp_path, **panel).rename(destination)
    return bundle


def test_the_frozen_check_reverifies_the_panel_that_actually_reached_the_bundle(tmp_path):
    bundle = frozen_bundle(tmp_path, format_version=2, with_modules=True)
    report = frozen_module().check_reference_panels(bundle)
    assert report["status"] == "passed" and report["unbundled_payload"] == []
    assert report["characterization"]["organism_modules"] == "staged"
    (bundle / "_internal/wmlstudio/resources/characterization/starter/wzi.fasta").unlink()
    with pytest.raises(ValueError, match="changed or is missing"):
        frozen_module().check_reference_panels(bundle)


def test_the_frozen_check_rejects_a_bundle_carrying_the_users_own_downloads(tmp_path):
    bundle = frozen_bundle(tmp_path, format_version=2, with_modules=True)
    (bundle / "_internal" / f"{organism_panel.PANEL_PREFIX}0123456789ab").mkdir()
    with pytest.raises(ValueError, match="never travel in the portable build"):
        frozen_module().check_reference_panels(bundle)


def test_an_unstaged_panel_withholds_the_sccmec_control_instead_of_reporting_a_negative(tmp_path):
    """A missing panel must read as 'not run, and here is why', never as 'nothing found'."""
    from wmlstudio.organism_modules import registered_modules
    frozen = frozen_module()
    panel = frozen.check_reference_panels(frozen_bundle(tmp_path, format_version=1))["characterization"]
    report = frozen.check_organism_modules(
        tmp_path / "absent-cli", tmp_path, panel, sorted(registered_modules()))
    assert report["status"] == "passed" and report["registered_modules"] == sorted(registered_modules())
    assert report["sccmec"]["status"] == "not_run"
    assert "format 1" in report["sccmec"]["reason"]


def test_a_frozen_report_without_module_evidence_fails_rather_than_passing_quietly(tmp_path):
    frozen = frozen_module()
    panel = frozen.check_reference_panels(frozen_bundle(tmp_path, format_version=1))["characterization"]
    with pytest.raises(ValueError, match="registry did not load"):
        frozen.check_organism_modules(
            tmp_path / "absent-cli", tmp_path, panel, ["species_evidence", "virulence"])


def test_threshold_documentation_reports_the_catalogue_as_it_actually_is():
    """The published audit must move when the catalogue moves, not drift behind it."""
    document = (ROOT / "docs/THRESHOLDS.md").read_text(encoding="utf-8")
    entries = threshold_guidance.catalog_entries()
    bindable = [entry for entry in entries
                if entry["scheme_key"] and entry["published_threshold"] is not None]
    covered = {entry["organism"] for entry in entries}
    uncurated = [organism for organism in threshold_guidance.ORGANISMS if organism not in covered]
    assert f"**{len(threshold_guidance.ORGANISMS)} organisms**" in document
    assert f"**{len(entries)} entries**" in document
    assert f"**{len(threshold_guidance.SOURCES)} cited sources**" in document
    assert f"| A — a published cutoff bound to a named scheme | **{len({e['organism'] for e in bindable})}**" in document
    assert f"| C — listed for surveillance, no curated publication | **{len(uncurated)}**" in document
    for organism in uncurated:
        assert organism in document, f"{organism} has no curated cutoff and is not named in the audit"
    assert threshold_guidance.CATALOG_VERSION in document


def test_threshold_documentation_never_promises_a_cutoff_the_catalogue_refuses():
    document = (ROOT / "docs/THRESHOLDS.md").read_text(encoding="utf-8")
    for entry in threshold_guidance.catalog_entries():
        if entry["scheme_key"] and entry["published_threshold"] is not None:
            assert entry["scheme_key"] in document, f"{entry['id']} is bindable but undocumented"
    # Tonnies curates a scheme and deliberately no number; the audit must say so.
    assert "**none**" in document and "no universal numeric rule" in document


def source_scheme(tmp_path):
    source = tmp_path / "source"
    directory = source / "tiny"
    directory.mkdir(parents=True)
    (directory / "abc.tfa").write_text(">abc_1\nACGT\n", encoding="utf-8")
    (directory / "tiny.txt").write_text("ST\tabc\n1\t1\n", encoding="utf-8")
    (directory / "tiny_info.json").write_text('{"source":"synthetic test"}', encoding="utf-8")
    (directory / "ignore.py").write_text("raise Exception('not data')", encoding="utf-8")
    return source


def test_staging_hashes_and_preserves_metadata(tmp_path):
    source = source_scheme(tmp_path)
    destination = tmp_path / "staged"
    first = stage.stage_schemes(source, destination)
    second = stage.stage_schemes(source, destination)
    assert first["scheme_count"] == 1
    assert first["file_count"] == 3
    assert first["snapshot_sha256"] == second["snapshot_sha256"]
    assert first["schemes"][0]["source_metadata"] == {"source": "synthetic test"}
    assert not (destination / "tiny/ignore.py").exists()
    for entry in first["files"]:
        assert stage.sha256(destination / entry["path"]) == entry["sha256"]
    assert json.loads((destination / "manifest.json").read_text())["files"] == first["files"]


def test_staging_refuses_stale_or_overlapping_destination(tmp_path):
    source = source_scheme(tmp_path)
    with pytest.raises(ValueError, match="separate"):
        stage.stage_schemes(source, source / "nested")
    destination = tmp_path / "staged"
    destination.mkdir()
    (destination / "unrelated.txt").write_text("keep me")
    with pytest.raises(ValueError, match="stale"):
        stage.stage_schemes(source, destination)
    assert (destination / "unrelated.txt").read_text() == "keep me"


@pytest.mark.parametrize("member_name,kind", [
    ("repo/db/pubmlst/../../escape.tfa", "file"),
    ("repo/db/pubmlst/tiny/link.tfa", "symlink"),
    ("repo/db/pubmlst/tiny/C:escape.tfa", "file"),
])
def test_download_archive_rejects_paths_and_links(tmp_path, member_name, kind):
    archive = tmp_path / "bad.tar.gz"
    with tarfile.open(archive, "w:gz") as handle:
        info = tarfile.TarInfo(member_name)
        if kind == "symlink":
            info.type = tarfile.SYMTYPE
            info.linkname = "../../other"
            handle.addfile(info)
        else:
            info.size = 4
            handle.addfile(info, io.BytesIO(b"ACGT"))
    with pytest.raises(ValueError):
        stage.unpack_schemes(archive, tmp_path / "out")


def test_download_archive_only_copies_scheme_data(tmp_path):
    archive = tmp_path / "good.tar.gz"
    with tarfile.open(archive, "w:gz") as handle:
        for name in ("repo/db/pubmlst/tiny/a.tfa", "repo/src/legacy.py"):
            info = tarfile.TarInfo(name)
            info.size = 4
            handle.addfile(info, io.BytesIO(b"ACGT"))
    result = stage.unpack_schemes(archive, tmp_path / "out")
    assert (result / "tiny/a.tfa").read_bytes() == b"ACGT"
    assert len(list(result.rglob("*.*"))) == 1


def test_validation_resume_requires_unchanged_bytes_and_limits(tmp_path):
    script = SCRIPT.with_name("validate_real_data.py")
    spec = importlib.util.spec_from_file_location("validate_real_data", script)
    validate = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(validate)
    reads = tmp_path / "reads.fastq"
    reads.write_text("@r1\nACGT\n+\nIIII\n@r2\nAAAA\n+\nIIII\n")
    report = tmp_path / "validation.json"
    arguments = ["--fastq", str(reads), "--max-reads", "1", "--output", str(report)]
    assert validate.main(arguments) == 0
    assert validate.main([*arguments, "--resume"]) == 0
    assert json.loads(report.read_text())["records"][0]["resumed"]
    reads.write_text("@r1\nTGCA\n+\nIIII\n@r2\nAAAA\n+\nIIII\n")
    assert validate.main([*arguments, "--resume"]) == 0
    assert not json.loads(report.read_text())["records"][0]["resumed"]
    arguments[3] = "2"
    assert validate.main([*arguments, "--resume"]) == 0
    assert not json.loads(report.read_text())["records"][0]["resumed"]
