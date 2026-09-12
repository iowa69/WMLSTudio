import importlib.util
import io
from pathlib import Path
from types import SimpleNamespace

import pytest


def notices_module():
    spec = importlib.util.spec_from_file_location("stage_notices", Path(__file__).resolve().parents[1]
                                                / "studio_packaging/stage_notices.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_incidental_distribution_notices_are_preserved(tmp_path, monkeypatch):
    module = notices_module()
    installed = tmp_path / "installed"
    (installed / "helper/licenses").mkdir(parents=True)
    license_path = installed / "helper/licenses/LICENSE"
    license_path.write_text("An incidental helper's attribution", encoding="utf-8")
    relative = Path("helper/licenses/LICENSE")
    package = SimpleNamespace(metadata={"Name": "Unexpected_Helper"}, version="1.2.3",
                              files=[relative], locate_file=lambda path: installed / path)
    monkeypatch.setattr(module, "distributions", lambda: [package])
    target, manifest = tmp_path / "notices", {"packages": {}, "files": []}
    module.preserve_installed_notices(target, manifest)
    assert manifest["packages"] == {"unexpected-helper": "1.2.3"}
    assert (target / manifest["files"][0]["path"]).read_bytes() == license_path.read_bytes()
    assert len(manifest["files"][0]["sha256"]) == 64


@pytest.mark.parametrize("url,authorized", [
    ("https://api.github.com/repos/qt/qtbase", True),
    ("https://raw.githubusercontent.com/qt/qtbase/main/LICENSE", False),
    ("https://ftp.ncbi.nlm.nih.gov/file", False),
    ("https://api.github.com.example.com/file", False),
    ("http://api.github.com/file", False),
])
def test_github_credential_is_limited_to_https_api_host(monkeypatch, url, authorized):
    module = notices_module()
    observed = []
    monkeypatch.setenv("GH_TOKEN", "unit-test-secret")
    monkeypatch.setattr(module.urllib.request, "urlopen", lambda request, **kwargs: observed.append(request) or io.BytesIO(b"license"))
    assert module.fetch(url) == b"license"
    assert observed[0].get_header("Authorization") == ("Bearer unit-test-secret" if authorized else None)
    if authorized:
        redirected = module.urllib.request.HTTPRedirectHandler().redirect_request(
            observed[0], None, 302, "redirect", {}, "https://raw.githubusercontent.com/file")
        assert redirected.get_header("Authorization") is None


@pytest.mark.parametrize("use_fallback", [False, True])
def test_windows_interpreter_full_notices_are_preserved(tmp_path, monkeypatch, use_fallback):
    module = notices_module()
    base, executable = tmp_path / "base", tmp_path / "runtime/python.exe"
    runtime = executable.parent if use_fallback else base
    runtime.mkdir(parents=True)
    original = runtime / "LICENSE.txt"
    original.write_bytes(b"Combined target interpreter license; OpenSSL and native runtime notices")
    monkeypatch.setattr(module.sys, "base_prefix", str(base))
    monkeypatch.setattr(module.sys, "executable", str(executable))
    destination = tmp_path / "notices"
    destination.mkdir()
    manifest = {"files": []}
    module.preserve_windows_python_notices(destination, manifest)
    assert (destination / manifest["files"][0]["path"]).read_bytes() == original.read_bytes()
    assert len(manifest["files"][0]["sha256"]) == 64


def test_windows_runtime_notices_are_required(tmp_path, monkeypatch):
    module = notices_module()
    monkeypatch.setattr(module.sys, "base_prefix", str(tmp_path / "no-base"))
    monkeypatch.setattr(module.sys, "executable", str(tmp_path / "no-runtime/python.exe"))
    with pytest.raises(ValueError, match="missing its bundled LICENSE"):
        module.preserve_windows_python_notices(tmp_path, {"files": []})
