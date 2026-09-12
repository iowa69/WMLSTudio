"""Snapshot provenance must be identical under Windows and POSIX path ordering."""

import importlib.util
import json
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / 'studio_scripts' / 'stage_schemes.py'
SPEC = importlib.util.spec_from_file_location('stage_schemes_order_test', SCRIPT)
stage = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(stage)


def test_manifest_uses_case_sensitive_posix_path_order_on_every_platform(tmp_path, monkeypatch):
    source = tmp_path / 'source'
    data = {
        'Zebra/Oxf_cpn60.tfa': b'>Oxf_cpn60_1\nACGT\n',
        'Zebra/abaumannii.txt': b'ST\tOxf_cpn60\n1\t1\n',
        'alpha/Beta.tfa': b'>Beta_1\nACGT\n',
        'alpha/alpha_info.json': b'{"name":"Alpha"}\n',
        'alpha-2/Zed.tfa': b'>Zed_1\nTGCA\n',
    }
    for relative, content in data.items():
        path = source / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(content)

    def windows_path_order(left, right):
        return str(left).casefold() < str(right).casefold()

    with monkeypatch.context() as patch:
        # Exercise the Windows failure on Linux as well, without emulating I/O.
        patch.setattr(type(tmp_path), '__lt__', windows_path_order)
        manifest = stage.stage_schemes(source, tmp_path / 'staged', provenance={'method': 'test'})

    paths = [entry['path'] for entry in manifest['files']]
    assert paths == [
        'Zebra/Oxf_cpn60.tfa', 'Zebra/abaumannii.txt', 'alpha-2/Zed.tfa',
        'alpha/Beta.tfa', 'alpha/alpha_info.json',
    ]
    assert paths == sorted(paths)
    assert [entry['id'] for entry in manifest['schemes']] == ['Zebra', 'alpha', 'alpha-2']
    assert manifest['file_count'] == 5
    assert json.loads((tmp_path / 'staged' / 'manifest.json').read_text())['files'] == manifest['files']
    for relative, content in data.items():
        assert (tmp_path / 'staged' / relative).read_bytes() == content
    # A build using the host's ordinary Path comparisons gets the exact same
    # content fingerprint, not merely matching individual file hashes.
    repeated = stage.stage_schemes(source, tmp_path / 'ordinary', provenance={'method': 'test'})
    assert repeated['snapshot_sha256'] == manifest['snapshot_sha256']
