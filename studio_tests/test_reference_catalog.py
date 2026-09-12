import io
import json
import urllib.error
import zipfile

import pytest

from wmlstudio.reference_catalog import (
    API_ROOT,
    CatalogError,
    CGMLSTOrgCatalog,
    PasteurCatalog,
    PubMLSTCatalog,
)
from wmlstudio.sequence import AnalysisCancelled
from wmlstudio.typing import call_assembly, load_scheme

DB = API_ROOT + "/db/pubmlst_example_seqdef"
SCHEME = DB + "/schemes/1"
ARC = "AACCGTACGTTAG"
GYR = "TTGGCATACCTGA"


class Response(io.BytesIO):
    def __init__(self, url, body, length=None):
        super().__init__(body)
        self.url = url
        self.headers = {"Content-Length": str(len(body) if length is None else length)}

    def geturl(self):
        return self.url


class Opener:
    def __init__(self, routes):
        self.routes = routes
        self.requests = []

    def open(self, request, timeout):
        assert request.method == "GET"
        assert request.data is None
        assert "Authorization" not in request.headers
        self.requests.append(request.full_url)
        payload = self.routes[request.full_url]
        if callable(payload):
            return payload(request.full_url)
        if isinstance(payload, Exception):
            raise payload
        raw = payload if isinstance(payload, bytes) else json.dumps(payload).encode()
        return Response(request.full_url, raw)


@pytest.fixture
def remote():
    details = {"id": 1, "description": "MLST", "locus_count": 2,
               "loci": [DB + "/loci/arc", DB + "/loci/gyr"],
               "last_updated": "2026-08-01", "records": 1,
               "profiles_csv": SCHEME + "/profiles_csv",
               "primary_key_field": SCHEME + "/fields/ST",
               "message": "Public data are restricted to records submitted before 2025."}
    routes = {
        API_ROOT + "/db": [{"name": "example", "description": "Examplegenus species",
                             "databases": [{"name": "pubmlst_example_seqdef",
                                            "description": "Examplegenus species sequence/profile definitions", "href": DB}]}],
        DB + "/schemes": {"schemes": [{"scheme": SCHEME, "description": "MLST"}], "records": 1},
        SCHEME: details,
        DB + "/loci/arc/alleles_fasta": f">arc_1\n{ARC}\n".encode(),
        DB + "/loci/gyr/alleles_fasta": f">gyr_1\n{GYR}\n".encode(),
        SCHEME + "/profiles_csv": b"ST\tarc\tgyr\n42\t1\t1\n",
    }
    opener = Opener(routes)
    return PubMLSTCatalog(opener=opener, retries=0), opener, routes


def entry(client):
    return client.list_schemes(client.list_organisms()[0])[0]


def test_catalog_discovers_organisms_filters_and_access_notice(remote):
    client, opener, routes = remote
    schemes = client.search_schemes("Examplegenus", scheme_type="MLST", min_loci=2, max_loci=7)
    assert schemes["errors"] == []
    assert schemes["organisms_searched"] == 1
    found = schemes["schemes"][0]
    assert found["locus_count"] == 2
    assert found["type"] == "MLST"
    assert found["organism"] == "Examplegenus species"
    assert "restricted" in found["access_notice"]
    assert client.search_schemes("Examplegenus", min_loci=3)["schemes"] == []
    assert all(url.startswith(API_ROOT) for url in opener.requests)
    routes[SCHEME]["description"] = "cgMLST"
    updated = client.list_schemes(client.list_organisms()[0], refresh=True)
    assert updated[0]["type"] == "cgMLST"


def test_valid_download_is_versioned_hashed_and_callable(tmp_path, remote):
    client, _, _ = remote
    result = client.download_scheme(entry(client), tmp_path / "library")
    assert result["created"] is True
    scheme = load_scheme(result["path"])
    assert scheme.digest == result["scheme_digest"]
    assert any("restricted" in note for note in scheme.notes)
    assembly = tmp_path / "assembly.fa"
    assembly.write_text(f">a\n{ARC}\n>b\n{GYR}\n")
    assert call_assembly(assembly, scheme)["st"] == "42"
    manifest = json.loads((scheme.path / "reference_manifest.json").read_text())
    assert manifest["authenticated"] is False
    assert len(manifest["sources"]) == 3
    assert all(len(source["sha256"]) == 64 for source in manifest["sources"])
    repeated = client.download_scheme(entry(client), tmp_path / "library")
    assert repeated["created"] is False
    assert repeated["path"] == result["path"]


def test_update_creates_new_snapshot_without_changing_previous(tmp_path, remote):
    client, _, routes = remote
    original = client.download_scheme(entry(client), tmp_path / "library")
    original_path = load_scheme(original["path"]).path
    old = (original_path / "profiles.tsv").read_bytes()
    routes[SCHEME + "/profiles_csv"] = b"ST\tarc\tgyr\n88\t1\t1\n"
    updated = client.update_scheme(original_path)
    assert updated["created"] is True
    assert updated["path"] != original["path"]
    assert (original_path / "profiles.tsv").read_bytes() == old
    assert len(list((tmp_path / "library").iterdir())) == 2


def test_failed_locus_download_never_installs_partial_reference(tmp_path, remote):
    client, _, routes = remote
    selected = entry(client)
    routes[DB + "/loci/gyr/alleles_fasta"] = urllib.error.URLError("network interrupted")
    with pytest.raises(CatalogError, match="reach PubMLST"):
        client.download_scheme(selected, tmp_path / "library")
    assert list((tmp_path / "library").iterdir()) == []


def test_truncation_and_download_quota_fail_without_publication(tmp_path, remote):
    client, _, routes = remote
    selected = entry(client)
    routes[DB + "/loci/arc/alleles_fasta"] = lambda url: Response(url, b">arc_1\nACGT\n", length=99)
    with pytest.raises(CatalogError, match="truncated"):
        client.download_scheme(selected, tmp_path / "library")
    assert list((tmp_path / "library").iterdir()) == []
    with pytest.raises(CatalogError, match="size limit"):
        client.download_scheme(selected, tmp_path / "library", max_bytes=2)


def test_cancellation_cleans_download_staging(tmp_path, remote):
    client, _, _ = remote
    cancelled = False

    def progress(current, total, message):
        nonlocal cancelled
        cancelled = True

    with pytest.raises(AnalysisCancelled):
        client.download_scheme(entry(client), tmp_path / "library",
                               cancelled=lambda: cancelled, progress=progress)
    assert list((tmp_path / "library").iterdir()) == []


def test_remote_change_during_download_is_rejected(tmp_path, remote):
    client, _, routes = remote
    selected = entry(client)
    initial = dict(routes[SCHEME])
    counter = 0

    def changing_details(url):
        nonlocal counter
        counter += 1
        details = dict(initial, records=1 if counter == 1 else 2)
        return Response(url, json.dumps(details).encode())

    routes[SCHEME] = changing_details
    with pytest.raises(CatalogError, match="changed during"):
        client.download_scheme(selected, tmp_path / "library")
    assert list((tmp_path / "library").iterdir()) == []


def test_external_resource_url_and_case_collisions_are_rejected(tmp_path, remote):
    client, opener, routes = remote
    selected = entry(client)
    routes[SCHEME]["loci"][0] = "https://example.com/private"
    with pytest.raises(CatalogError, match="official HTTPS"):
        client.download_scheme(selected, tmp_path / "library")
    assert all("example.com" not in url for url in opener.requests)
    routes[SCHEME]["loci"] = [DB + "/loci/arc", DB + "/loci/ARC"]
    with pytest.raises(CatalogError, match="Windows case"):
        client.download_scheme(selected, tmp_path / "library")
    assert list((tmp_path / "library").iterdir()) == []


def test_tampered_existing_snapshot_is_not_overwritten(tmp_path, remote):
    client, _, _ = remote
    selected = entry(client)
    installed = client.download_scheme(selected, tmp_path / "library")
    path = load_scheme(installed["path"]).path / "profiles.tsv"
    changed = "ST\tarc\tgyr\n999\t1\t1\n"
    path.write_text(changed)
    with pytest.raises(CatalogError, match="modified"):
        client.download_scheme(selected, tmp_path / "library")
    assert path.read_text() == changed


def test_nonstandard_cg_profile_identifier_is_retained(tmp_path, remote):
    client, _, routes = remote
    routes[SCHEME]["primary_key_field"] = SCHEME + "/fields/cgST"
    routes[SCHEME + "/profiles_csv"] = b"cgST\tarc\tgyr\n123\t1\t1\n"
    installed = client.download_scheme(entry(client), tmp_path / "library")
    scheme = load_scheme(installed["path"])
    assert scheme.metadata["profile_field"] == "cgST"
    assert scheme.profiles[("1", "1")] == ("123",)


def test_pasteur_uses_official_api_and_has_provider_namespaced_identity(remote):
    _, _, routes = remote
    base = 'https://bigsdb.pasteur.fr/api'
    converted = {key.replace(API_ROOT, base): json.loads(json.dumps(value).replace(API_ROOT, base))
                 for key, value in routes.items() if not isinstance(value, bytes)}
    client = PasteurCatalog(opener=Opener(converted), retries=0)
    result = client.search_schemes('Examplegenus')['schemes'][0]
    assert result['provider'] == 'BIGSdb-Pasteur'
    assert result['url'].startswith(base + '/db/')
    assert result['id'].startswith('BIGSdb-Pasteur:')


def cg_remote(members=None):
    base = 'https://www.cgmlst.org/ncs/schema/Test/'
    archive = io.BytesIO()
    with zipfile.ZipFile(archive, 'w') as handle:
        for name, value in (members or {'bundle/arc.fasta': f'>arc_1\n{ARC}\n',
                                         'bundle/gyr.fasta': f'>gyr_1\n{GYR}\n'}).items():
            handle.writestr(name, value)
    routes = {
        'https://www.cgmlst.org/ncs': b'<table><tr><td><a href="/ncs/schema/Test/">Examplegenus species</a></td><td>2</td></tr></table>',
        base + 'locus/?content-type=csv': b'Locus\tName\n\tarc\tArc\n\tgyr\tGyr\n',
        base + 'alleles/': archive.getvalue(),
    }
    return CGMLSTOrgCatalog(opener=Opener(routes), retries=0)


def test_cgmlst_org_html_archive_and_manifest_are_real_local_scheme(tmp_path):
    client = cg_remote()
    selected = client.search_schemes('Examplegenus', scheme_type='cgMLST')['schemes'][0]
    assert selected['locus_count'] == 2
    with pytest.raises(CatalogError, match='explicitly acknowledge'):
        client.download_scheme(selected, tmp_path)
    installed = client.download_scheme(selected, tmp_path, terms_acknowledged=True)
    scheme = load_scheme(installed['path'])
    assert scheme.loci == ('arc', 'gyr')
    assert scheme.profile_count == 0
    assert not (scheme.path / 'source_alleles.zip').exists()
    manifest = json.loads((scheme.path / 'reference_manifest.json').read_text())
    assert len(manifest['sources']) == 4
    assert len(manifest['sources'][1]['sha256']) == 64
    assert manifest['user_acknowledged_permitted_use'] is True
    assert 'serverpolicy' in manifest['terms_url']
    assert client.download_scheme(selected, tmp_path, terms_acknowledged=True)['created'] is False


@pytest.mark.parametrize('members,reason', [
    ({'../arc.fasta': '>arc_1\nACGT\n'}, 'Unsafe path'),
    ({'C:/arc.fasta': '>arc_1\nACGT\n'}, 'Unsafe path'),
    ({'arc.fasta': '>arc_1\nACGT\n'}, 'missing loci'),
    ({'unexpected.fasta': '>unexpected_1\nACGT\n'}, 'not in'),
    ({'a/arc.fasta': '>arc_1\nACGT\n', 'b/arc.fasta': '>arc_1\nACGT\n'}, 'colliding'),
])
def test_cg_archive_failures_never_publish_partial_snapshot(tmp_path, members, reason):
    client = cg_remote(members)
    selected = client.search_schemes('Examplegenus')['schemes'][0]
    with pytest.raises(CatalogError, match=reason):
        client.download_scheme(selected, tmp_path, terms_acknowledged=True)
    assert list(tmp_path.iterdir()) == []
