import http.client
import io
import json
import urllib.error
import zipfile
from pathlib import Path

import pytest

from wmlstudio import cgmlst_schemes
from wmlstudio.reference_catalog import (
    API_ROOT,
    CatalogError,
    CGMLSTOrgCatalog,
    PasteurCatalog,
    PubMLSTCatalog,
    install_destination,
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


def published(root):
    """Scheme folders a user would see, in either library.

    A kept partial download is hidden and is not one: reference_index.
    filter_scheme_locations drops '.'-prefixed names. A labelled but empty cgMLST
    slot is not one either: it holds a README and no allele file.
    """
    bases = [Path(root), Path(root) / cgmlst_schemes.LIBRARY_DIRNAME]
    return sorted(path.name for base in bases if base.is_dir() for path in base.iterdir()
                  if path.is_dir() and not path.name.startswith(".")
                  and cgmlst_schemes.has_alleles(path))


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
    assert published(tmp_path / "library") == []


def test_truncation_and_download_quota_fail_without_publication(tmp_path, remote):
    client, _, routes = remote
    selected = entry(client)
    routes[DB + "/loci/arc/alleles_fasta"] = lambda url: Response(url, b">arc_1\nACGT\n", length=99)
    with pytest.raises(CatalogError, match="truncated"):
        client.download_scheme(selected, tmp_path / "library")
    assert published(tmp_path / "library") == []
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
    assert published(tmp_path / "library") == []


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
    assert published(tmp_path / "library") == []


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
    assert published(tmp_path / "library") == []


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
    assert published(tmp_path) == []


def failing_once(error, then):
    """A route that raises once and then answers, to exercise retry behaviour."""
    state = {"raised": False}

    def handler(url):
        if not state["raised"]:
            state["raised"] = True
            raise error
        return Response(url, then)

    return handler


def test_non_iupac_allele_is_excluded_and_named_instead_of_failing_the_scheme(tmp_path):
    # cgMLST.org publishes occasional alleles containing a literal 'X'. One such
    # record used to abort the whole multi-gigabyte download at validation.
    client = cg_remote({'bundle/arc.fasta': f'>arc_1\n{ARC}\n>arc_2\nACGXTACG\n',
                        'bundle/gyr.fasta': f'>gyr_1\n{GYR}\n'})
    selected = client.search_schemes('Examplegenus')['schemes'][0]
    installed = client.download_scheme(selected, tmp_path, terms_acknowledged=True)
    scheme = load_scheme(installed['path'])
    assert scheme.loci == ('arc', 'gyr')
    assert set(scheme.alleles['arc']) == {'1'}, 'the unreadable record is gone, the locus stays'
    assert installed['excluded_alleles'] == ["arc:arc_2 (non-IUPAC X)"]
    manifest = json.loads((scheme.path / 'reference_manifest.json').read_text())
    assert manifest['excluded_alleles'] == ["arc:arc_2 (non-IUPAC X)"]
    note = next(text for text in installed['notes'] if 'excluded' in text)
    assert 'not rewritten' in note and 'unmatched sequence' in note
    # The received bytes stay auditable; the stored bytes are recorded separately.
    source = next(item for item in manifest['sources'] if item['file'] == 'arc.fasta')
    assert source['sha256'] != source['stored_sha256']


def test_a_locus_whose_every_allele_is_unreadable_installs_nothing(tmp_path):
    client = cg_remote({'bundle/arc.fasta': '>arc_1\nACGXT\n',
                        'bundle/gyr.fasta': f'>gyr_1\n{GYR}\n'})
    selected = client.search_schemes('Examplegenus')['schemes'][0]
    with pytest.raises(CatalogError, match='Every allele record'):
        client.download_scheme(selected, tmp_path, terms_acknowledged=True)
    assert published(tmp_path) == []


def test_pubmlst_non_iupac_allele_is_excluded_and_the_snapshot_still_installs(tmp_path, remote):
    client, _, routes = remote
    routes[DB + "/loci/arc/alleles_fasta"] = f">arc_1\n{ARC}\n>arc_9\nACG?T\n".encode()
    installed = client.download_scheme(entry(client), tmp_path / "library")
    assert installed['excluded_alleles'] == ["arc:arc_9 (non-IUPAC ?)"]
    scheme = load_scheme(installed['path'])
    assert set(scheme.alleles['arc']) == {'1'}
    assert call_assembly(_assembly(tmp_path), scheme)["st"] == "42"


def _assembly(tmp_path):
    path = tmp_path / "assembly.fa"
    path.write_text(f">a\n{ARC}\n>b\n{GYR}\n")
    return path


def test_a_scheme_error_is_reported_as_a_catalog_error_not_leaked(tmp_path, remote):
    # A scheme failure during a download is a download failure. Callers documented
    # to catch CatalogError must not also have to catch SchemeError from typing.
    client, _, routes = remote
    routes[DB + "/loci/arc/alleles_fasta"] = f">arc_1\n{ARC}\n>arc_1\n{GYR}\n".encode()
    with pytest.raises(CatalogError, match='could not be read as a scheme'):
        client.download_scheme(entry(client), tmp_path / "library")
    assert published(tmp_path / "library") == []


def test_an_entirely_unreadable_locus_stops_the_download(tmp_path, remote):
    client, _, routes = remote
    routes[DB + "/loci/arc/alleles_fasta"] = b">arc_1\n\n"
    with pytest.raises(CatalogError, match='Every allele record'):
        client.download_scheme(entry(client), tmp_path / "library")
    assert published(tmp_path / "library") == []


def test_an_interrupted_download_resumes_instead_of_refetching_every_locus(tmp_path, remote):
    client, opener, routes = remote
    selected = entry(client)
    routes[DB + "/loci/gyr/alleles_fasta"] = urllib.error.URLError("network interrupted")
    with pytest.raises(CatalogError, match="reach PubMLST"):
        client.download_scheme(selected, tmp_path / "library")
    state = client.resume_state(tmp_path / "library", selected)
    assert state["files_held"] == 1 and state["expected"] == 2
    assert Path(state["path"]).name.startswith(".resume-")
    # A partial download is never a scheme, and it is kept OUTSIDE the library that
    # paths.scheme_locations hands to the typing code.
    assert list((tmp_path / "library").iterdir()) == []
    assert Path(state["path"]).parent.name == ".wmlstudio-partial-downloads"
    routes[DB + "/loci/gyr/alleles_fasta"] = f">gyr_1\n{GYR}\n".encode()
    opener.requests.clear()
    installed = client.download_scheme(selected, tmp_path / "library")
    assert installed["created"] is True
    assert DB + "/loci/arc/alleles_fasta" not in opener.requests, "already held, not refetched"
    assert DB + "/loci/gyr/alleles_fasta" in opener.requests
    assert call_assembly(_assembly(tmp_path), load_scheme(installed["path"]))["st"] == "42"
    assert client.resume_state(tmp_path / "library", selected) is None


def test_a_changed_remote_scheme_discards_the_partial_rather_than_blending_it(tmp_path, remote):
    client, opener, routes = remote
    selected = entry(client)
    routes[DB + "/loci/gyr/alleles_fasta"] = urllib.error.URLError("network interrupted")
    with pytest.raises(CatalogError):
        client.download_scheme(selected, tmp_path / "library")
    assert client.resume_state(tmp_path / "library", selected)["files_held"] == 1
    routes[SCHEME]["last_updated"] = "2026-09-14"
    routes[DB + "/loci/gyr/alleles_fasta"] = f">gyr_1\n{GYR}\n".encode()
    routes[DB + "/loci/arc/alleles_fasta"] = f">arc_1\n{ARC}\n".encode()
    opener.requests.clear()
    client.download_scheme(selected, tmp_path / "library")
    assert DB + "/loci/arc/alleles_fasta" in opener.requests, "the stale copy was not trusted"


def test_a_partial_download_can_be_discarded_without_touching_an_installed_scheme(tmp_path, remote):
    client, _, routes = remote
    selected = entry(client)
    installed = client.download_scheme(selected, tmp_path / "library")
    routes[SCHEME]["records"] = 2
    routes[DB + "/loci/gyr/alleles_fasta"] = urllib.error.URLError("network interrupted")
    with pytest.raises(CatalogError):
        client.download_scheme(selected, tmp_path / "library")
    assert client.clear_resume(tmp_path / "library", selected) is True
    assert client.clear_resume(tmp_path / "library", selected) is False
    assert Path(installed["path"]).is_dir()
    assert published(tmp_path / "library") == [Path(installed["path"]).name]


def test_resume_can_be_switched_off_and_then_leaves_nothing_behind(tmp_path, remote):
    client, _, routes = remote
    selected = entry(client)
    routes[DB + "/loci/gyr/alleles_fasta"] = urllib.error.URLError("network interrupted")
    with pytest.raises(CatalogError):
        client.download_scheme(selected, tmp_path / "library", resume=False)
    assert list((tmp_path / "library").iterdir()) == []
    assert not (tmp_path / ".wmlstudio-partial-downloads").exists()
    assert client.resume_state(tmp_path / "library", selected) is None


def test_rate_limiting_is_waited_out_and_then_reported_in_words(tmp_path, remote):
    client, _, routes = remote
    client.retries = 1
    routes[DB + "/loci/arc/alleles_fasta"] = failing_once(
        urllib.error.HTTPError(DB, 429, "Too Many Requests", {"Retry-After": "0"}, None),
        f">arc_1\n{ARC}\n".encode())
    installed = client.download_scheme(entry(client), tmp_path / "library")
    assert installed["created"] is True
    routes[DB + "/loci/gyr/alleles_fasta"] = urllib.error.HTTPError(
        DB, 429, "Too Many Requests", {"Retry-After": "0"}, None)
    with pytest.raises(CatalogError, match="rate-limiting this computer"):
        client.download_scheme(entry(client), tmp_path / "library2")


def test_a_truncated_chunked_response_is_retried_not_raised_as_an_http_error(tmp_path, remote):
    client, _, routes = remote
    client.retries = 1
    routes[DB + "/loci/arc/alleles_fasta"] = failing_once(
        http.client.IncompleteRead(b"", 10), f">arc_1\n{ARC}\n".encode())
    installed = client.download_scheme(entry(client), tmp_path / "library")
    assert load_scheme(installed["path"]).locus_count == 2


def test_cgmlst_org_error_page_is_named_instead_of_read_as_a_locus_table(tmp_path):
    client = cg_remote()
    selected = client.search_schemes('Examplegenus')['schemes'][0]
    client.opener.routes['https://www.cgmlst.org/ncs/schema/Test/locus/?content-type=csv'] = (
        b'<html><body><b>ERROR Occured #20260914</b><br>Illegal Typing ID!')
    with pytest.raises(CatalogError, match='does not recognise the scheme identifier'):
        client.download_scheme(selected, tmp_path, terms_acknowledged=True)
    assert published(tmp_path) == []


def test_a_retargeted_scheme_is_refused_with_both_counts_named(tmp_path):
    client = cg_remote()
    selected = client.search_schemes('Examplegenus')['schemes'][0]
    selected['locus_count'] = 3
    with pytest.raises(CatalogError, match='2 targets but the catalogue entry says 3'):
        client.download_scheme(selected, tmp_path, terms_acknowledged=True)
    assert published(tmp_path) == []


def test_download_estimates_warn_before_a_thousand_request_transfer():
    per_locus = PubMLSTCatalog(opener=Opener({})).download_estimate(
        {"locus_count": 1748, "has_profiles": True})
    assert per_locus["requests"] == 1749
    assert per_locus["resumable"] is True
    assert "1748 separate downloads" in per_locus["notice"]
    assert "cancelled and resumed" in per_locus["notice"]
    small = PubMLSTCatalog(opener=Opener({})).download_estimate({"locus_count": 7})
    assert small["requests"] == 7 and "quick" in small["notice"]
    archive = CGMLSTOrgCatalog(opener=Opener({})).download_estimate({"locus_count": 1861})
    assert archive["requests"] == 2 and archive["resumable"] is False
    assert "gigabytes" in archive["notice"]


def test_download_progress_reports_bytes_while_a_length_free_archive_streams(tmp_path):
    payload = b'>arc_1\n' + b'ACGT' * 2_000_000 + b'\n'
    client = cg_remote({'bundle/arc.fasta': payload.decode(), 'bundle/gyr.fasta': f'>gyr_1\n{GYR}\n'})
    selected = client.search_schemes('Examplegenus')['schemes'][0]
    seen = []
    client.download_scheme(selected, tmp_path, terms_acknowledged=True,
                           progress=lambda current, total, text: seen.append(text))
    assert any('MB received' in text for text in seen), 'a long transfer must not look frozen'


def cg_scheme(client, tmp_path, *, root=None):
    selected = client.search_schemes('Examplegenus')['schemes'][0]
    return client.download_scheme(selected, root or tmp_path / "schemes",
                                  terms_acknowledged=True)


def test_a_downloaded_cgmlst_scheme_lands_in_the_cgmlst_library_under_a_readable_name(tmp_path):
    # The reported bug in one assertion: a downloaded cgMLST scheme used to land in
    # <data root>/schemes as cgmlst_org_<slug>_<digest16>, sorted alphabetically
    # among 162 seven-locus schemes under a name nobody could recognise.
    installed = cg_scheme(cg_remote(), tmp_path)
    path = Path(installed['path'])
    assert path.parent == tmp_path / cgmlst_schemes.LIBRARY_DIRNAME
    assert path.name == 'Examplegenus_species__cgmlst_org_2'
    assert 'cgmlst_org_Test_' not in path.name
    assert load_scheme(path).locus_count == 2
    # Nothing was left in the classical library the caller named.
    assert not (tmp_path / "schemes").exists() or published(tmp_path / "schemes") == []


def test_a_download_fills_the_labelled_slot_the_library_already_describes(tmp_path):
    # cgmlst_schemes.prepare_library pre-creates a README-bearing folder for every
    # catalogued scheme. A download must install INTO it, not beside it.
    cgmlst_schemes.prepare_library(tmp_path, keys=["cgmlst.org:saureus-1861"])
    slot = cgmlst_schemes.library_root(tmp_path) / \
        cgmlst_schemes.entry_for("cgmlst.org:saureus-1861")["slot"]
    assert slot.is_dir() and not cgmlst_schemes.has_alleles(slot)
    client = cg_remote()
    selected = client.search_schemes('Examplegenus')['schemes'][0]
    selected.update({'slug': 'Saureus', 'organism': 'Staphylococcus aureus',
                     'name': 'Staphylococcus aureus cgMLST'})
    client.opener.routes.update({
        'https://www.cgmlst.org/ncs/schema/Saureus/locus/?content-type=csv':
            client.opener.routes['https://www.cgmlst.org/ncs/schema/Test/locus/?content-type=csv'],
        'https://www.cgmlst.org/ncs/schema/Saureus/alleles/':
            client.opener.routes['https://www.cgmlst.org/ncs/schema/Test/alleles/']})
    installed = client.download_scheme(selected, tmp_path / "schemes", terms_acknowledged=True)
    assert Path(installed['path']) == slot
    assert (slot / 'README.txt').is_file(), 'the licence note stayed with the scheme'
    assert 'Ridom' in (slot / 'README.txt').read_text(encoding='utf-8')
    assert (slot / 'scheme_slot.json').is_file()
    assert load_scheme(slot).locus_count == 2


def test_a_classical_scheme_keeps_its_content_addressed_folder(tmp_path, remote):
    # Saved projects store these exact paths, and a seven-locus scheme is small:
    # only the gene-by-gene library is reorganised.
    client, _, _ = remote
    installed = client.download_scheme(entry(client), tmp_path / "schemes")
    path = Path(installed['path'])
    assert path.parent == tmp_path / "schemes"
    assert path.name.startswith('pubmlst_example_seqdef_1_')
    assert not (tmp_path / cgmlst_schemes.LIBRARY_DIRNAME).exists()


def test_a_pubmlst_cgmlst_scheme_is_routed_by_its_target_count_not_its_folder(tmp_path, remote):
    client, _, routes = remote
    routes[SCHEME]['description'] = 'cgMLST'
    selected = client.list_schemes(client.list_organisms()[0], refresh=True)[0]
    assert selected['type'] == 'cgMLST'
    installed = client.download_scheme(selected, tmp_path / "schemes")
    assert Path(installed['path']).parent == tmp_path / cgmlst_schemes.LIBRARY_DIRNAME
    assert Path(installed['path']).name == 'Examplegenus_species__pubmlst_2'
    assert load_scheme(installed['path']).locus_count == 2
    # A kept partial download still lives outside BOTH libraries, so it is never
    # handed to the typing code as a scheme.
    assert published(tmp_path / "schemes") == []
    assert not (tmp_path / cgmlst_schemes.LIBRARY_DIRNAME / ".wmlstudio-partial-downloads").exists()


def test_install_destination_never_sends_two_kinds_to_one_library(tmp_path):
    classical = {'database': 'pubmlst_example_seqdef', 'scheme_id': '1', 'locus_count': 7,
                 'type': 'MLST'}
    gene_by_gene = {'database': 'pubmlst_example_seqdef', 'scheme_id': '6', 'locus_count': 2513,
                    'type': 'cgMLST', 'organism': 'Escherichia coli', 'provider': 'PubMLST'}
    library = tmp_path / "schemes"
    assert install_destination(library, classical, 'f' * 64).parent == library
    assert install_destination(library, gene_by_gene, 'f' * 64).parent == \
        tmp_path / cgmlst_schemes.LIBRARY_DIRNAME
    # A small but explicitly declared cgMLST snapshot is still gene-by-gene.
    partial = {**gene_by_gene, 'locus_count': 12}
    assert install_destination(library, partial, 'f' * 64).parent == \
        tmp_path / cgmlst_schemes.LIBRARY_DIRNAME
