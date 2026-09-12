"""Pinned native tool and starter-database staging safeguards."""

import hashlib
import importlib.util
import io
import json
import tarfile
from pathlib import Path

import pytest

_spec = importlib.util.spec_from_file_location("stage_bio_tools", Path(__file__).resolve().parents[1]
                                             / "studio_packaging/stage_bio_tools.py")
staging = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(staging)


def test_windows_tls_filter_preserves_pinned_python_and_native_schannel():
    retain = [(name, "source", "BINARY") for name in (
        "PySide6/plugins/tls/qschannelbackend.dll", "PySide6/plugins/tls/qcertonlybackend.dll",
        "libcrypto-3.dll", "libssl-3.dll", "_ssl.pyd", "Tools/blast/bin/vcruntime140.dll")]
    ambient = [(name, "runner-PATH", "BINARY") for name in (
        "PySide6\\plugins\\tls\\qopensslbackend.dll", "libcrypto-3-x64.dll", "LIBSSL-3-X64.DLL")]
    original = retain + ambient
    assert staging.filter_windows_qt_tls(original) == retain
    assert len(original) == 9  # No mutation of installed files or original TOC.


def archive(tmp_path, monkeypatch, extras=None):
    path = tmp_path / "official-test.tar.gz"
    entries = {f"bin/{name}.exe": b"MZ-test-fixture" for name in staging.PROGRAMS}
    entries.update({"LICENSE": b"test license", "README": b"test tools"})
    entries.update(extras or {})
    with tarfile.open(path, "w:gz") as bundle:
        for relative, content in entries.items():
            info = tarfile.TarInfo(f"ncbi-blast-{staging.VERSION}+/{relative}")
            info.size = len(content)
            bundle.addfile(info, io.BytesIO(content))
    monkeypatch.setitem(staging.PACKAGES, "windows-x64", {"archive": path.name, "suffix": ".exe",
        "sha256": hashlib.sha256(path.read_bytes()).hexdigest()})
    return path


def test_stage_requires_matching_archive_hash_and_never_publishes_partial(tmp_path):
    source = tmp_path / "corrupt.tar.gz"
    source.write_bytes(b"not the pinned archive")
    destination = tmp_path / "tools"
    with pytest.raises(ValueError, match="SHA-256 mismatch"):
        staging.stage("windows-x64", destination, source)
    assert not destination.exists()
    assert not list(tmp_path.glob("wmlstudio-blast-stage-*"))


def test_verified_programs_license_hashes_and_idempotence(tmp_path, monkeypatch):
    source = archive(tmp_path, monkeypatch, {"bin/unused.exe": b"not-needed"})
    destination = tmp_path / "tools"
    manifest = staging.stage("windows-x64", destination, source)
    assert (destination / "bin/blastn.exe").read_bytes() == b"MZ-test-fixture"
    assert not (destination / "bin/unused.exe").exists()
    assert (destination / "LICENSE").is_file()
    assert manifest == staging.stage("windows-x64", destination, source)
    (destination / "bin/blastn.exe").write_bytes(b"changed")
    with pytest.raises(ValueError, match="already exists"):
        staging.stage("windows-x64", destination, source)


def test_empty_manifest_cannot_fake_a_complete_tool_install(tmp_path):
    destination = tmp_path / "tools"
    destination.mkdir()
    (destination / "manifest.json").write_text(json.dumps({
        "archive_sha256": staging.PACKAGES["windows-x64"]["sha256"], "files": []}))
    with pytest.raises(ValueError, match="already exists"):
        staging.stage("windows-x64", destination)


def test_archive_path_traversal_is_rejected_before_publication(tmp_path, monkeypatch):
    source = archive(tmp_path, monkeypatch, {"../outside.exe": b"not allowed"})
    destination = tmp_path / "tools"
    with pytest.raises(ValueError, match="Unsafe"):
        staging.stage("windows-x64", destination, source)
    assert not destination.exists()
    assert not (tmp_path / "outside.exe").exists()


def test_hydra_starter_never_silently_redistributes_other_providers(tmp_path):
    root = tmp_path / "source"
    (root / "nucl/card").mkdir(parents=True)
    (root / "manifest.json").write_text(json.dumps({"databases": {"card": {"path": "nucl/card"}}}))
    with pytest.raises(ValueError, match="exactly the NCBI"):
        staging.stage_hydra_database(root, tmp_path / "staged")
    assert not (tmp_path / "staged").exists()


def skesa_fixture(tmp_path):
    source = tmp_path / "skesa-artifact"
    source.mkdir()
    records = []
    for name in ("skesa.exe", "SKESA-LICENSE.txt", "AGPL-3.0.txt",
                 "skesa-2.4.0-wmlstudio-source.tar.gz", "dependency-licenses/runtime/LICENSE"):
        path = source / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(b"synthetic staging test; not a real binary")
        records.append({"path": name, "sha256": staging.digest(path)})
    (source / "manifest.json").write_text(json.dumps({"tool": "SKESA", "platform": "windows-x64-ucrt",
        "source_commit": "c1413581e4f37211892d3c4310d01f3d9a9b3490", "files": records}))
    (source / "smoke-test.json").write_text(json.dumps({"all_contigs_match_reference": True,
                                                     "native_adapter": {"passed": True}}))
    return source


def test_skesa_requires_sources_licenses_integrity_and_native_smoke(tmp_path):
    source = skesa_fixture(tmp_path)
    destination = tmp_path / "staged"
    first = staging.stage_skesa_bundle(source, destination)
    assert first == staging.stage_skesa_bundle(source, destination)
    (destination / "skesa.exe").write_bytes(b"tampered")
    with pytest.raises(ValueError, match="integrity"):
        staging.verify_skesa_bundle(destination)
    (source / "smoke-test.json").write_text('{}')
    with pytest.raises(ValueError, match="smoke evidence"):
        staging.verify_skesa_bundle(source)


def test_skesa_manifest_cannot_omit_corresponding_source(tmp_path):
    source = skesa_fixture(tmp_path)
    path = source / "manifest.json"
    manifest = json.loads(path.read_text())
    manifest["files"] = [item for item in manifest["files"] if not item["path"].endswith(".tar.gz")]
    path.write_text(json.dumps(manifest))
    with pytest.raises(ValueError, match="missing required"):
        staging.verify_skesa_bundle(source)
