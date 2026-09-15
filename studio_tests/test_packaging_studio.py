import ast
import importlib.util
import io
import json
import tarfile
import zipfile
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


def budget_module():
    spec = importlib.util.spec_from_file_location(
        "check_size_budget", ROOT / "studio_packaging/check_size_budget.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def sized_bundle(tmp_path):
    """A frozen-bundle shape whose parts are the components the gate reports."""
    root = tmp_path / "WMLSTudio"
    payload = {
        "_internal/Tools/blast/bin/blastn.exe": b"native search " * 400,
        "_internal/Tools/fastp/fastp": b"read trimming " * 300,
        "_internal/Tools/fastqc/jre/bin/java.dll": b"private java runtime " * 500,
        "_internal/wmlstudio/resources/schemes/manifest.json": b'{"scheme_count": 1}',
        "_internal/wmlstudio/resources/hydra/starter/manifest.json": b'{"databases": {}}',
        "_internal/PySide6/Qt6Core.dll": b"qt " * 900,
        "WMLSTudio.exe": b"MZ application",
    }
    for relative, content in payload.items():
        path = root / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(content)
    return root


def test_the_portable_download_is_measured_per_component_and_gated(tmp_path):
    """A size promise nobody measures is a size promise that quietly stops being true."""
    budget = budget_module()
    report = budget.check(budget.measure(sized_bundle(tmp_path)))
    named = {entry["component"]: entry for entry in report["components"]}
    assert "fastp — read trimming" in named
    assert named["BLAST+ — the nucleotide and translated search"]["files"] == 1
    assert named["FastQC and its private Java runtime"]["files"] == 1
    assert sum(entry["files"] for entry in report["components"]) == report["files"] == 7
    assert report["status"] == "within_budget" and report["budget_bytes"] == 1_000_000_000
    # The same measurement, against a budget this bundle cannot meet, must fail.
    over = budget.check(budget.measure(sized_bundle(tmp_path)), budget=10)
    assert over["status"] == "over_budget" and over["headroom_bytes"] < 0
    warned = budget.check(budget.measure(sized_bundle(tmp_path)),
                          budget=int(report["total_bytes"] / 0.9))
    assert warned["status"] == "approaching_budget"


def test_every_byte_of_the_package_is_attributed_to_some_component(tmp_path):
    """An unattributed byte is exactly how a bundle grows without anyone noticing."""
    budget = budget_module()
    report = budget.measure_tree(sized_bundle(tmp_path), compress=False)
    assert sum(entry["bytes"] for entry in report["components"]) == report["total_bytes"]
    assert budget.component_for("nothing/known/here.dll") == budget.OTHER
    assert budget.component_for("wmlstudio/resources/tools/skesa/skesa.exe").startswith("SKESA")


def test_a_directory_prediction_is_not_smaller_than_the_archive_it_predicts(tmp_path):
    """The gate on the built folder must not pass a package the real ZIP would fail."""
    budget = budget_module()
    bundle = sized_bundle(tmp_path)
    archive = tmp_path / "portable.zip"
    with zipfile.ZipFile(archive, "w", zipfile.ZIP_DEFLATED) as handle:
        for path in sorted(bundle.rglob("*")):
            if path.is_file():
                handle.write(path, Path("WMLSTudio", path.relative_to(bundle)).as_posix())
    predicted = budget.measure_tree(bundle)
    measured = budget.measure_archive(archive)
    assert predicted["total_bytes"] >= measured["total_bytes"]
    assert {entry["component"] for entry in predicted["components"]} == \
           {entry["component"] for entry in measured["components"]}


def test_bytes_on_disk_are_never_compared_with_the_download_budget(tmp_path):
    """Extracted size is a different quantity; reporting it as headroom would mislead."""
    budget = budget_module()
    report = budget.check(budget.measure_tree(sized_bundle(tmp_path), compress=False), budget=None)
    assert report["status"] == "reported" and report["gated"] is False
    assert "budget_bytes" not in report and "headroom_bytes" not in report
    assert "not the download size" in report["measured"]


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


def test_the_package_bundles_verified_fastp_where_it_exists_and_builds_without_it(tmp_path):
    """Upstream ships no Windows fastp, so an absent tool is a build state, not a failure."""
    import wmlstudio.read_tools as read_tools
    recipe = RECIPE.read_text(encoding="utf-8")
    assert "from stage_read_tools import verify as verify_fastp" in recipe
    assert 'datas.append((str(fastp), "Tools/fastp"))' in recipe
    # The staged tree must land exactly where the application looks for it.
    assert 'Tools/fastp' in Path(read_tools.__file__).read_text(encoding="utf-8")
    # Nothing about fastp may stop a build, and its libraries must not be hoisted
    # to the application root any more than SKA2's or the JRE's are.
    fastp_block = recipe.split("fastp = root /")[1].split("\nif json.loads(")[0]
    assert 'if (fastp / "manifest.json").is_file():' in fastp_block
    assert "raise SystemExit" not in fastp_block and "ships without read trimming" in fastp_block
    assert "(fastqc, ska, tools, fastp)" in recipe
    # With nothing staged the application states the absence rather than trimming.
    capabilities = read_tools.runtime_capabilities(root=tmp_path / "nothing-staged")
    assert capabilities["available"] is False
    assert "unavailable" in capabilities["reason"] and "download" in capabilities["reason"]


def core_reference_bundle(tmp_path, *, dna=True):
    """A frozen bundle carrying only the AMR core store, in its packaged location."""
    store = tmp_path / "WMLSTudio/_internal/wmlstudio/resources/hydra/starter"
    (store / "nucl/ncbi").mkdir(parents=True)
    (store / "prot/protein").mkdir(parents=True)
    (store / "prot/protein/AMRProt-mutation.tsv").write_text(
        "organism\tgene\nEscherichia\tgyrA\n", encoding="utf-8")
    (store / "prot/protein/taxgroup.tsv").write_text(
        "organism\tgroup\nEscherichia\tEscherichia\nBurkholderia_mallei\tBurkholderia\n",
        encoding="utf-8")
    (store / "prot/protein/meta.tsv").write_text(
        "gene\telement_type\nblaZ\tAMR\nfimH\tVIRULENCE\n", encoding="utf-8")
    if dna:
        (store / "mutation/dna").mkdir(parents=True)
        (store / "mutation/dna/Escherichia.fna").write_text(">gyrA\nACGT\n", encoding="utf-8")
    (store / "manifest.json").write_text(json.dumps({"databases": {
        "ncbi": {"path": "nucl/ncbi"}, "protein": {"path": "prot/protein"}}}), encoding="utf-8")
    return tmp_path / "WMLSTudio"


def test_the_frozen_bundle_must_carry_the_point_mutations_it_promises(tmp_path):
    """An unperformed mutation search and a clean mutation result look identical."""
    report = frozen_module().check_core_reference(core_reference_bundle(tmp_path))
    assert report["status"] == "passed"
    assert report["dna_point_mutation_catalogues"] == ["Escherichia"]
    # An accepted organism with no catalogue is named, not quietly counted as covered.
    assert report["accepted_without_any_catalogue"] == ["Burkholderia_mallei"]
    assert "never a negative mutation result" in report["boundary"]
    assert report["elements"]["VIRULENCE"] == 1 and report["elements"]["read"] is True


def test_a_frozen_bundle_that_lost_its_mutation_catalogues_fails_the_build(tmp_path):
    with pytest.raises(ValueError, match="point-mutation"):
        frozen_module().check_core_reference(core_reference_bundle(tmp_path, dna=False))


def test_the_frozen_check_names_an_absent_fastp_instead_of_reporting_trimming_done(tmp_path):
    bundle = tmp_path / "WMLSTudio"
    (bundle / "_internal").mkdir(parents=True)
    report = frozen_module().check_fastp(bundle, tmp_path, "")
    assert report["status"] == "not_bundled"
    assert "unavailable" in report["reason"] and "untrimmed" in report["reason"]


def test_the_capability_audit_stops_calling_fastp_preprocessing_absent():
    """The audit is the page a reader checks before trusting a claim; it must move."""
    document = (ROOT / "docs/STUDIO_CAPABILITIES.md").read_text(encoding="utf-8")
    assert "no fastp/SPAdes/long-read pipeline" not in document
    assert "fastp preprocessing, validated" not in document
    assert "fastp 1.3.7" in document
    guide = (ROOT / "docs/WORKFLOW_GUIDE.md").read_text(encoding="utf-8")
    assert "fastp" in guide and "does not validate an isolate" in guide
    # A tool that is present on one platform and absent on another must say which.
    assert "trimming is unavailable" in guide


def test_the_audit_never_claims_interface_reach_for_a_tab_that_says_it_is_planned():
    """The Reach column is the audit's whole point: it must track the built tabs.

    A station listed in `ui_tabs.PLANNED` prints "nothing runs on this tab yet".
    Publishing "Interface" for the same capability would make the audit the more
    optimistic of the two documents a reader can check, which is exactly backwards.
    """
    from wmlstudio.ui_tabs import PLANNED, TAB_LABELS
    document = (ROOT / "docs/STUDIO_CAPABILITIES.md").read_text(encoding="utf-8")
    revision = document.split("## What the portable package costs")[0]
    rows = [line for line in revision.splitlines() if line.startswith("| ") and line.endswith(" |")]
    for key in PLANNED:
        label = TAB_LABELS[key]
        for row in rows:
            if label.casefold() in row.casefold():
                assert not row.rstrip(" |").endswith("Interface"), (
                    f"the audit claims Interface reach for {label}, which ui_tabs still "
                    f"lists as planned: {row}")


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
