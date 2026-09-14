import importlib.util
import io
import json
import tarfile
import zipfile
from pathlib import Path

import pytest

SPEC = importlib.util.spec_from_file_location("stage_fastqc", Path(__file__).resolve().parents[1] / "studio_packaging/stage_fastqc.py")
stage = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(stage)


@pytest.mark.parametrize("name", ["../escape", "/absolute/path", "root/../escape", "root/C:escape", "root\\escape"])
def test_fastqc_archive_rejects_unsafe_names(tmp_path, name):
    archive = tmp_path / "bad.zip"
    with zipfile.ZipFile(archive, "w") as handle:
        entry = zipfile.ZipInfo("placeholder")
        # Windows ZipInfo normally normalizes backslashes at construction;
        # retain the malicious bytes to test the archive reader, not the writer.
        entry.filename = name
        handle.writestr(entry, "payload")
    with pytest.raises(ValueError, match="Unsafe"):
        stage.extract(archive, tmp_path / "output")


def test_fastqc_archive_rejects_casefold_duplicate_paths(tmp_path):
    archive = tmp_path / "duplicate.zip"
    with zipfile.ZipFile(archive, "w") as handle:
        handle.writestr("root/File", "one")
        handle.writestr("root/file", "two")
    with pytest.raises(ValueError, match="Duplicate"):
        stage.extract(archive, tmp_path / "output")


def test_jre_license_links_are_preserved_as_safe_regular_copies(tmp_path):
    archive = tmp_path / "licenses.tar"
    with tarfile.open(archive, "w") as handle:
        license = tarfile.TarInfo("jre/legal/java.base/LICENSE")
        license.size = 3
        handle.addfile(license, io.BytesIO(b"GPL"))
        link = tarfile.TarInfo("jre/legal/java.xml/LICENSE")
        link.type, link.linkname = tarfile.SYMTYPE, "../java.base/LICENSE"
        handle.addfile(link)
    stage.extract(archive, tmp_path / "output")
    target = tmp_path / "output/legal/java.xml/LICENSE"
    assert target.read_text() == "GPL" and not target.is_symlink()


def test_jre_license_link_may_not_escape_archive(tmp_path):
    archive = tmp_path / "bad.tar"
    with tarfile.open(archive, "w") as handle:
        link = tarfile.TarInfo("jre/link")
        link.type, link.linkname = tarfile.SYMTYPE, "../../outside"
        handle.addfile(link)
    with pytest.raises(ValueError):
        stage.extract(archive, tmp_path / "output")


def test_cached_upstream_archive_is_still_hash_checked(tmp_path):
    name = stage.PACKAGES["fastqc"][0]
    (tmp_path / name).write_bytes(b"not the original binary")
    with pytest.raises(ValueError, match="checksum"):
        stage.acquire("fastqc", tmp_path, offline=True)
    assert (tmp_path / name).read_bytes() == b"not the original binary"
    with pytest.raises(ValueError, match="not cached"):
        stage.acquire("java_source", tmp_path, offline=True)


def test_staged_manifest_rejects_modified_bytes_and_wrong_target(tmp_path):
    data = tmp_path / "LICENSE"
    data.write_text("GPL")
    manifest = {"fastqc_version": stage.FASTQC_VERSION, "java_version": stage.JRE_VERSION,
                "platform": "windows-x64", "files": [{"path": "LICENSE", "sha256": stage.digest(data)}]}
    (tmp_path / "manifest.json").write_text(json.dumps(manifest))
    assert stage.verify(tmp_path, "windows-x64") == manifest
    with pytest.raises(ValueError, match="platform"):
        stage.verify(tmp_path, "linux-x64")
    data.write_text("changed")
    with pytest.raises(ValueError, match="checksum"):
        stage.verify(tmp_path)
