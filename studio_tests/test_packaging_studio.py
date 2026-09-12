import importlib.util
import io
import json
import tarfile
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[1] / "studio_scripts/stage_schemes.py"
SPEC = importlib.util.spec_from_file_location("stage_schemes", SCRIPT)
stage = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(stage)


def test_frozen_pyrodigal_explicitly_collects_namespace_implementations():
    from PyInstaller.utils.hooks import collect_submodules
    modules = collect_submodules("pyrodigal.impl")
    assert "pyrodigal.impl" in modules and "pyrodigal.impl.generic" in modules
    recipe = (SCRIPT.parent.parent / "studio_packaging/wmlstudio.spec").read_text()
    assert 'collect_submodules("pyrodigal.impl")' in recipe


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
