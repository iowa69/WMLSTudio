import importlib.util
import json
from pathlib import Path

import pytest

SPEC = importlib.util.spec_from_file_location("stage_ska", Path(__file__).resolve().parents[1] / "studio_packaging/stage_ska.py")
stage = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(stage)


def fixture_bundle(root):
    root.mkdir()
    files = {"ska.exe": b"fixture", "LICENSE": b"Apache", "NOTICE": b"Attribution",
             "Cargo.lock": b"locked", "ska-0.5.1-corresponding-source.tar.gz": b"source archive",
             "smoke-test.json": b'{"passed":true}'}
    for name, content in files.items():
        (root / name).write_bytes(content)
    manifest = {"tool": "SKA2", "version": stage.VERSION, "source_commit": stage.REVISION,
                "platform": "windows-x64", "files": [{"path": name, "sha256": stage.digest(root / name)} for name in files]}
    (root / "manifest.json").write_text(json.dumps(manifest))
    return manifest


def test_stage_only_manifest_listed_source_and_notice_files(tmp_path):
    source, target = tmp_path / "source", tmp_path / "target"
    manifest = fixture_bundle(source)
    (source / "unrelated.txt").write_text("do not distribute")
    assert stage.stage(target, "windows-x64", source) == manifest
    assert stage.stage(target, "windows-x64") == manifest
    assert not (target / "unrelated.txt").exists()
    assert (source / "unrelated.txt").read_text() == "do not distribute"
    assert (target / "ska-0.5.1-corresponding-source.tar.gz").is_file()


def test_stage_rejects_modified_binary_or_missing_source_and_keeps_original(tmp_path):
    source = tmp_path / "source"
    fixture_bundle(source)
    (source / "ska.exe").write_bytes(b"different")
    with pytest.raises(ValueError, match="checksum"):
        stage.stage(tmp_path / "target", "windows-x64", source)
    assert not (tmp_path / "target").exists()
    assert (source / "ska.exe").read_bytes() == b"different"


def test_ska_platform_and_unsafe_manifest_paths_are_rejected(tmp_path):
    manifest = fixture_bundle(tmp_path / "source")
    with pytest.raises(ValueError, match="platform"):
        stage.verify(tmp_path / "source", "linux-x64")
    manifest["files"][0]["path"] = "../escape"
    (tmp_path / "source/manifest.json").write_text(json.dumps(manifest))
    with pytest.raises(ValueError, match="Unsafe"):
        stage.verify(tmp_path / "source", "windows-x64")


def test_windows_stage_cannot_silently_omit_native_engine(tmp_path):
    with pytest.raises(ValueError, match="artifact"):
        stage.stage(tmp_path / "out", "windows-x64")
