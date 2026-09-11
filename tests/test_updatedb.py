# SPDX-License-Identifier: GPL-2.0-only
# Copyright (C) 2025-2026 IOWA-Tech - Giovanni Lorenzin
# Copyright (C) Torsten Seemann (upstream `mlst`, from which WMLST is ported)
"""Tests for wmlst/updatedb.py — docs/ARCHITECTURE.md sections 4.9, 6, 7.

Runnable two ways::

    python3 -m pytest tests/test_updatedb.py
    python3 tests/test_updatedb.py

The offline tests drive the updater through a fake HTTP opener, so the whole
update algorithm (check -> apply -> hash oracle -> rollback -> bundle) is
exercised without a network. The tests marked ``live`` really talk to
rest.pubmlst.org and bigsdb.pasteur.fr and skip themselves when the hosts are
unreachable or when WMLST_SKIP_NET is set.
"""

from __future__ import annotations

import glob
import io
import json
import os
import shutil
import sys
import tempfile
import urllib.error

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from wmlst import updatedb as U
from wmlst.engine import Cancelled, UpdateError

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BUNDLED_DB = os.path.join(REPO, "db")
BUNDLED_SCHEMES = os.path.join(BUNDLED_DB, "pubmlst")


# ---------------------------------------------------------------------------
# fake upstream
# ---------------------------------------------------------------------------

class _FakeHandle:
    def __init__(self, status, headers, body):
        self.status = status
        self.headers = headers
        self._body = body

    def read(self):
        return self._body

    def getcode(self):
        return self.status

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False


class FakeOpener:
    """Minimal stand-in for a urllib opener; records every URL it serves."""

    def __init__(self, routes):
        self.routes = routes
        self.calls = []

    def open(self, request, timeout=None):
        url = request.full_url if hasattr(request, "full_url") else str(request)
        self.calls.append(url)
        route = self.routes.get(url)
        if route is None:
            raise urllib.error.HTTPError(
                url, 404, "Not Found", {"Content-Type": "application/json"},
                io.BytesIO(b'{"message":"Scheme has not been defined"}'))
        if callable(route):
            route = route(url, len([c for c in self.calls if c == url]))
        status, content_type, body = route
        headers = {"Content-Type": content_type}
        if status != 200:
            raise urllib.error.HTTPError(url, status, "err", headers, io.BytesIO(body))
        return _FakeHandle(status, headers, body)


ROOT = U.REST_ROOTS["pubmlst"]
DB = "pubmlst_tiny_seqdef"
API = "%s/db/%s/schemes/1" % (ROOT, DB)
LOCI = ("abc", "def")

PROFILES_V1 = b"ST\tabc\tdef\n1\t1\t1\n2\t1\t2\n"
PROFILES_V2 = b"ST\tabc\tdef\n1\t1\t1\n2\t1\t2\n3\t2\t2\n"
ALLELES_V1 = {"abc": b">abc_1\nACGTACGTAA\n", "def": b">def_1\nTTGGCCAATT\n>def_2\nTTGGCCAATG\n"}
ALLELES_V2 = {"abc": b">abc_1\nACGTACGTAA\n>abc_2\nACGTACGTAC\n",
              "def": b">def_1\nTTGGCCAATT\n>def_2\nTTGGCCAATG\n"}


def _scheme_doc(last_updated="2026-01-01", records=2, locus_count=2):
    return json.dumps({
        "id": 1, "description": "MLST", "locus_count": locus_count,
        "last_updated": last_updated, "records": records,
        "loci": ["%s/db/%s/loci/%s" % (ROOT, DB, line) for line in LOCI],
        "fields": ["%s/fields/ST" % API],
    }).encode("utf-8")


def _routes(profiles, alleles, **kwargs):
    text = "text/plain; charset=UTF-8"
    routes = {API: (200, "application/json", _scheme_doc(**kwargs)),
              API + "/profiles_csv": (200, text, profiles)}
    for locus, body in alleles.items():
        routes["%s/db/%s/loci/%s/alleles_fasta" % (ROOT, DB, locus)] = (200, text, body)
    return routes


def _fetcher(routes):
    return U._Fetcher(opener=FakeOpener(routes), delay=0.0, sleep=lambda s: None)


def _tiny_db(tmpdir, profiles=PROFILES_V1, alleles=ALLELES_V1, info_extra=None):
    dbdir = os.path.join(tmpdir, "db")
    scheme_dir = os.path.join(dbdir, "pubmlst", "tiny")
    os.makedirs(scheme_dir)
    with open(os.path.join(scheme_dir, "tiny.txt"), "wb") as fh:
        fh.write(profiles)
    for locus, body in alleles.items():
        with open(os.path.join(scheme_dir, locus + ".tfa"), "wb") as fh:
            fh.write(body)
    info = {"name": "tiny", "description": "MLST", "locus": 2,
            "download_date": "2026-01-01", "last_updated": "2026-01-01",
            "source": "pubmlst", "API": API, "authenticated": False}
    info.update(info_extra or {})
    U.write_scheme_info(os.path.join(scheme_dir, "tiny_info.json"), info)
    U.stamp_db_version(dbdir, "2026-01-01")
    U.write_manifest(dbdir, U.load_manifest(dbdir))
    return dbdir


# ---------------------------------------------------------------------------
# 4.9 content_sha256 — the change oracle
# ---------------------------------------------------------------------------

def test_content_sha256_matches_the_specified_recipe():
    """4.9: sha256 of a fixed prefix, the profiles, then byte-sorted loci."""
    import hashlib
    expected = hashlib.sha256()
    expected.update(b"WMLST-DBv1\n")
    expected.update(PROFILES_V1)
    for locus in ("abc", "def"):
        expected.update(locus.encode("utf-8") + b"\x00")
        expected.update(ALLELES_V1[locus])
    assert U.content_sha256(PROFILES_V1, ALLELES_V1) == expected.hexdigest()


def test_content_sha256_is_insensitive_to_dict_order():
    reordered = {k: ALLELES_V1[k] for k in reversed(list(ALLELES_V1))}
    assert U.content_sha256(PROFILES_V1, reordered) == \
        U.content_sha256(PROFILES_V1, ALLELES_V1)


def test_content_sha256_notices_a_one_byte_edit():
    """7.2: a row-count check would miss Pasteur's listeria_2 curation edit."""
    edited = PROFILES_V1.replace(b"2\t1\t2", b"2\t1\t3")
    assert len(edited) == len(PROFILES_V1)
    assert U.content_sha256(edited, ALLELES_V1) != U.content_sha256(PROFILES_V1, ALLELES_V1)


def test_local_content_sha256_reads_the_files_on_disk():
    with tempfile.TemporaryDirectory() as tmp:
        dbdir = _tiny_db(tmp)
        assert U.local_content_sha256(dbdir, "tiny") == \
            U.content_sha256(PROFILES_V1, ALLELES_V1)
        assert U.local_content_sha256(dbdir, "absent") is None


# ---------------------------------------------------------------------------
# 6.3 the manifest
# ---------------------------------------------------------------------------

def test_bundled_manifest_covers_every_scheme_directory():
    refs = U.load_manifest(BUNDLED_DB)
    names = {r.name for r in refs}
    on_disk = {n for n in os.listdir(BUNDLED_SCHEMES)
               if os.path.isdir(os.path.join(BUNDLED_SCHEMES, n))
               and not n.startswith(".")}
    assert names == on_disk
    assert len(refs) == 162


def test_bundled_manifest_has_no_unresolved_rows():
    unresolved = [r.name for r in U.load_manifest(BUNDLED_DB) if not r.resolved]
    assert unresolved == []


def test_manifest_triples_are_unique_except_for_the_alias_pair():
    refs = U.load_manifest(BUNDLED_DB)
    primaries = [r for r in refs if r.alias_of is None]
    triples = [r.triple for r in primaries]
    assert len(triples) == len(set(triples))
    aliases = {r.name: r.alias_of for r in refs if r.alias_of}
    assert aliases == {"cdiphtheriae": "diphtheria_3"}


def test_manifest_records_the_three_documented_dirname_violations():
    """6.3: salmonella -> 2, mgenitalium -> 2, cdiphtheriae -> diphtheria scheme 3."""
    refs = {r.name: r for r in U.load_manifest(BUNDLED_DB)}
    assert refs["salmonella"].scheme_id == "2"
    assert refs["mgenitalium"].scheme_id == "2"
    assert refs["cdiphtheriae"].db == "pubmlst_diphtheria_seqdef"
    assert refs["cdiphtheriae"].scheme_id == "3"
    assert refs["cdiphtheriae"].source == "pasteur"


def test_dirname_derivation_holds_for_161_of_162():
    """6.3: a validation assertion only, never the resolver."""
    mismatched = [r.name for r in U.load_manifest(BUNDLED_DB)
                  if U._derive_dirname(r.db, r.scheme_id) != r.name]
    assert sorted(mismatched) == ["cdiphtheriae", "mgenitalium", "salmonella"]


def test_manifest_lists_exactly_nine_pasteur_schemes():
    refs = U.load_manifest(BUNDLED_DB)
    pasteur = sorted(r.name for r in refs if r.source == "pasteur")
    assert pasteur == ["bordetella_3", "cdiphtheriae", "diphtheria_3", "ecoli",
                       "kingella", "klebsiella", "listeria_2", "staphlugdunensis",
                       "streptothermophilus"]


def test_manifest_nloci_matches_the_tfa_count():
    for ref in U.load_manifest(BUNDLED_DB):
        tfas = glob.glob(os.path.join(BUNDLED_SCHEMES, ref.name, "*.tfa"))
        assert ref.nloci == len(tfas), ref.name


def test_manifest_round_trips_through_disk():
    refs = U.load_manifest(BUNDLED_DB)
    with tempfile.TemporaryDirectory() as tmp:
        dbdir = os.path.join(tmp, "db")
        os.makedirs(os.path.join(dbdir, "pubmlst"))
        for ref in refs:
            os.makedirs(os.path.join(dbdir, "pubmlst", ref.name))
        U.write_manifest(dbdir, refs)
        assert U.load_manifest(dbdir) == refs


def test_unresolved_rows_survive_a_round_trip():
    with tempfile.TemporaryDirectory() as tmp:
        dbdir = os.path.join(tmp, "db")
        os.makedirs(os.path.join(dbdir, "pubmlst", "mystery"))
        refs = U.load_manifest(dbdir)
        assert refs[0].source == U.UNRESOLVED and not refs[0].resolved
        U.write_manifest(dbdir, refs)
        assert U.load_manifest(dbdir) == refs
        assert "unresolved" in open(U.manifest_path(dbdir)).read()


def test_api_url_is_rebuilt_from_the_triple():
    refs = {r.name: r for r in U.load_manifest(BUNDLED_DB)}
    for name in refs:
        info = U.scheme_info(BUNDLED_DB, name)
        assert refs[name].api == info["API"], name


# ---------------------------------------------------------------------------
# 6.2 _info.json
# ---------------------------------------------------------------------------

def test_every_bundled_info_json_round_trips_byte_for_byte():
    """6.2: 8 canonical keys in order, 2-space indent, trailing newline."""
    for path in sorted(glob.glob(os.path.join(BUNDLED_SCHEMES, "*", "*_info.json"))):
        with open(path, "rb") as fh:
            raw = fh.read()
        assert U._info_bytes(json.loads(raw)) == raw, path


def test_info_json_puts_wmlst_keys_after_the_canonical_eight():
    info = json.loads(U._info_bytes({
        "wmlst_content_sha256": "x", "authenticated": False, "name": "tiny",
        "description": "MLST", "locus": 2, "download_date": "d",
        "last_updated": "u", "source": "pubmlst", "API": API}))
    assert list(info)[:8] == ["name", "description", "locus", "download_date",
                              "last_updated", "source", "API", "authenticated"]
    assert list(info)[8] == "wmlst_content_sha256"


def test_no_version_information_is_written_verbatim():
    """6.2: mgenitalium legitimately has no upstream version string."""
    info = U.scheme_info(BUNDLED_DB, "mgenitalium")
    assert info["last_updated"] == U.NO_VERSION
    assert U._remote_date({"last_updated": None}) == U.NO_VERSION
    assert U._remote_date({}) == U.NO_VERSION


# ---------------------------------------------------------------------------
# 7.6 validation gates
# ---------------------------------------------------------------------------

def _response(body, content_type="text/plain; charset=UTF-8", status=200):
    return U._Response(status, {"Content-Type": content_type}, body, "http://x")


def test_validate_tfa_accepts_a_bundled_file():
    with open(os.path.join(BUNDLED_SCHEMES, "kingella", "abcZ.tfa"), "rb") as fh:
        U.validate_tfa("abcZ", _response(fh.read()))


@pytest.mark.parametrize("body,why", [
    (b"", "empty"),
    (b"ACGT\n", "does not start with '>'"),
    (b">abcZ_1\r\nACGT\n", "carriage return"),
    (b">wrong_1\nACGT\n", "malformed header"),
    (b">abcZ_1\nACGT\n\n>abcZ_2\nACGT\n", "blank line"),
    (b"<html><body>login</body></html>\n", "does not start with '>'"),
])
def test_validate_tfa_rejects_bad_payloads(body, why):
    with pytest.raises(UpdateError):
        U.validate_tfa("abcZ", _response(body))


def test_validate_tfa_rejects_a_json_message_body():
    """7.6: never accept a {"message": ...} body as allele data."""
    with pytest.raises(UpdateError):
        U.validate_tfa("abcZ", _response(b'{"message":"Please authenticate"}',
                                         content_type="application/json"))


def test_validate_tfa_rejects_a_non_text_content_type():
    with pytest.raises(UpdateError):
        U.validate_tfa("abcZ", _response(b">abcZ_1\nACGT\n", "text/html"))


def test_validate_profiles_accepts_the_bundled_table():
    with open(os.path.join(BUNDLED_SCHEMES, "kingella", "kingella.txt"), "rb") as fh:
        body = fh.read()
    loci = [os.path.basename(p)[:-4] for p in
            glob.glob(os.path.join(BUNDLED_SCHEMES, "kingella", "*.tfa"))]
    rows = U.validate_profiles("kingella", _response(body), loci, ("ST",))
    assert rows == 76


def test_validate_profiles_allows_zero_data_rows():
    """7.6: mgenitalium ships a 59-byte header-only profile table."""
    body = b"ST\tabc\tdef\n"
    assert U.validate_profiles("tiny", _response(body), ["abc", "def"], ("ST",)) == 0


def test_validate_profiles_rejects_a_locus_set_mismatch():
    with pytest.raises(UpdateError):
        U.validate_profiles("tiny", _response(b"ST\tabc\n1\t1\n"),
                            ["abc", "def"], ("ST",))


def test_validate_profiles_rejects_a_ragged_row():
    with pytest.raises(UpdateError):
        U.validate_profiles("tiny", _response(b"ST\tabc\tdef\n1\t1\n"),
                            ["abc", "def"], ("ST",))


def test_validate_profiles_ignores_the_six_non_locus_columns():
    body = b"ST\tabc\tdef\tclonal_complex\n1\t1\t1\tcc1\n"
    assert U.validate_profiles("tiny", _response(body), ["abc", "def"], ("ST",)) == 1


def test_locus_names_that_would_break_windows_are_rejected():
    """7.6: a ':' would silently create an NTFS alternate data stream."""
    for bad in ("ad:k", "ad/k", "ad\\k", "ad*k", 'ad"k', "ad|k", ""):
        with pytest.raises(UpdateError):
            U._validate_locus_name(bad)
    for good in ("abcZ", "EF-2_N", "rpoB'_N", "Pas_cpn60"):
        U._validate_locus_name(good)


def test_unparseable_loci_are_recorded():
    """5.19(b): the three loci the sseqid regex cannot parse stay visible."""
    assert U._unparseable_loci(["abcZ", "EF-2_N", "rpoB'_N"]) == ("EF-2_N", "rpoB'_N")


# ---------------------------------------------------------------------------
# 7.4 etiquette: retry, backoff, circuit breaker
# ---------------------------------------------------------------------------

def test_a_503_is_retried_and_then_succeeds():
    state = {"n": 0}

    def flaky(url, count):
        state["n"] += 1
        if state["n"] < 3:
            return (503, "text/plain", b"busy")
        return (200, "text/plain", b"ok")

    fetcher = _fetcher({"https://rest.pubmlst.org/x": flaky})
    assert fetcher.get("https://rest.pubmlst.org/x").body == b"ok"
    assert state["n"] == 3


def test_a_404_is_never_retried():
    opener = FakeOpener({})
    fetcher = U._Fetcher(opener=opener, delay=0.0, sleep=lambda s: None)
    with pytest.raises(U.HttpStatusError) as excinfo:
        fetcher.get("https://rest.pubmlst.org/db/nope/schemes/1")
    assert excinfo.value.code == 404
    assert excinfo.value.is_undefined_scheme
    assert len(opener.calls) == 1


def test_exhausted_retries_raise_an_actionable_update_error():
    fetcher = _fetcher({"https://rest.pubmlst.org/x": (503, "text/plain", b"busy")})
    with pytest.raises(UpdateError) as excinfo:
        fetcher.get("https://rest.pubmlst.org/x")
    assert "attempts" in str(excinfo.value)


def test_the_user_agent_names_the_vendor_and_the_author():
    assert "IOWA-Tech" in U.USER_AGENT and "Giovanni Lorenzin" in U.USER_AGENT
    assert U.USER_AGENT.startswith("WMLST/")


def test_cancellation_is_honoured_before_a_request():
    class Flag:
        def is_set(self):
            return True
    fetcher = U._Fetcher(opener=FakeOpener({}), delay=0.0, cancel=Flag(),
                         sleep=lambda s: None)
    with pytest.raises(Cancelled):
        fetcher.get("https://rest.pubmlst.org/x")


# ---------------------------------------------------------------------------
# 7.2 / 7.3 check and apply, driven by the fake upstream
# ---------------------------------------------------------------------------

def test_check_writes_nothing_and_reports_changed():
    with tempfile.TemporaryDirectory() as tmp:
        dbdir = _tiny_db(tmp)
        before = _snapshot(dbdir)
        plan = U.check(dbdir, fetcher=_fetcher(
            _routes(PROFILES_V2, ALLELES_V2, records=3)))
        assert [u.status for u in plan.updates] == ["changed"]
        assert plan.updates[0].added_types == 1     # 3 upstream rows vs 2 local
        assert plan.total_bytes > 0
        assert _snapshot(dbdir) == before


def test_check_uses_the_ledger_date_to_skip_a_scheme():
    """7.2: the cheap fast path costs zero payload bytes."""
    with tempfile.TemporaryDirectory() as tmp:
        dbdir = _tiny_db(tmp, info_extra={
            "wmlst_upstream_last_updated_at_verification": "2026-01-01"})
        fetcher = _fetcher(_routes(PROFILES_V2, ALLELES_V2))
        plan = U.check(dbdir, fetcher=fetcher)
        assert [u.status for u in plan.updates] == ["up-to-date"]
        assert fetcher._opener.calls == [API]


def test_a_retired_scheme_keeps_its_local_data():
    """7.2: only a 404 'has not been defined' justifies retired-upstream."""
    with tempfile.TemporaryDirectory() as tmp:
        dbdir = _tiny_db(tmp)
        plan = U.check(dbdir, fetcher=_fetcher({}))
        assert [u.status for u in plan.updates] == ["retired-upstream"]
        assert os.path.isfile(os.path.join(dbdir, "pubmlst", "tiny", "tiny.txt"))


def test_a_locus_count_change_quarantines_the_scheme():
    with tempfile.TemporaryDirectory() as tmp:
        dbdir = _tiny_db(tmp)
        routes = _routes(PROFILES_V2, ALLELES_V2, locus_count=9)
        plan = U.check(dbdir, fetcher=_fetcher(routes))
        assert [u.status for u in plan.updates] == ["quarantined"]


def test_apply_commits_new_content_and_updates_the_ledger():
    with tempfile.TemporaryDirectory() as tmp:
        dbdir = _tiny_db(tmp)
        fetcher = _fetcher(_routes(PROFILES_V2, ALLELES_V2, last_updated="2026-02-02"))
        plan = U.check(dbdir, fetcher=fetcher)
        version = U.apply(dbdir, plan, None, fetcher=fetcher, tools=_FakeTools(),
                          force_blast=False)
        scheme = os.path.join(dbdir, "pubmlst", "tiny")
        with open(os.path.join(scheme, "tiny.txt"), "rb") as fh:
            assert fh.read() == PROFILES_V2
        with open(os.path.join(scheme, "abc.tfa"), "rb") as fh:
            assert fh.read() == ALLELES_V2["abc"]
        info = U.scheme_info(dbdir, "tiny")
        assert info["wmlst_content_sha256"] == U.content_sha256(PROFILES_V2, ALLELES_V2)
        assert info["wmlst_upstream_last_updated_at_verification"] == "2026-02-02"
        assert info["last_updated"] == "2026-02-02"
        assert U.db_version(dbdir) == version == U._today()


def test_a_moved_date_with_identical_content_rewrites_nothing():
    """7.2: upstream's last_updated advances daily; the payload does not."""
    with tempfile.TemporaryDirectory() as tmp:
        dbdir = _tiny_db(tmp)
        scheme = os.path.join(dbdir, "pubmlst", "tiny")
        before = {n: os.path.getmtime(os.path.join(scheme, n))
                  for n in os.listdir(scheme) if not n.endswith("_info.json")}
        fetcher = _fetcher(_routes(PROFILES_V1, ALLELES_V1, last_updated="2026-09-10"))
        plan = U.check(dbdir, fetcher=fetcher)
        assert plan.updates[0].status == "changed"      # the date moved
        U.apply(dbdir, plan, None, fetcher=fetcher, tools=_FakeTools())
        after = {n: os.path.getmtime(os.path.join(scheme, n))
                 for n in os.listdir(scheme) if not n.endswith("_info.json")}
        assert after == before                          # the content did not
        commits = [r for r in U.read_journal(dbdir) if r["event"] == "scheme_commit"]
        assert commits[-1]["status"] == "date-moved-content-identical"
        assert U.scheme_info(dbdir, "tiny")[
            "wmlst_upstream_last_updated_at_verification"] == "2026-09-10"


def test_a_failed_validation_leaves_the_live_directory_untouched():
    """7.6: any gate failure deletes the staging directory and commits nothing."""
    with tempfile.TemporaryDirectory() as tmp:
        dbdir = _tiny_db(tmp)
        before = _snapshot(dbdir)
        routes = _routes(PROFILES_V2, {"abc": b">WRONG_1\nACGT\n",
                                       "def": ALLELES_V1["def"]})
        fetcher = _fetcher(routes)
        plan = U.check(dbdir, fetcher=fetcher)
        with pytest.raises(UpdateError):
            U.apply(dbdir, plan, None, fetcher=fetcher, tools=_FakeTools())
        assert _snapshot(dbdir) == before
        assert not os.path.isdir(os.path.join(dbdir, "pubmlst", U.STAGING_DIRNAME,
                                              "tiny.%d" % os.getpid()))
        errors = [r for r in U.read_journal(dbdir) if r["event"] == "scheme_error"]
        assert errors


def test_a_shrunken_download_is_refused():
    with tempfile.TemporaryDirectory() as tmp:
        big = {"abc": b"".join(b">abc_%d\n%s\n" % (i, b"ACGT" * 40)
                               for i in range(1, 60)),
               "def": ALLELES_V1["def"]}
        dbdir = _tiny_db(tmp, alleles=big)
        fetcher = _fetcher(_routes(PROFILES_V1, ALLELES_V1))
        plan = U.check(dbdir, fetcher=fetcher)
        with pytest.raises(UpdateError) as excinfo:
            U.apply(dbdir, plan, None, fetcher=fetcher, tools=_FakeTools())
        assert "half the size" in str(excinfo.value)


def test_rollback_restores_the_previous_content():
    with tempfile.TemporaryDirectory() as tmp:
        dbdir = _tiny_db(tmp)
        original = _snapshot(dbdir)["pubmlst/tiny/tiny.txt"]
        fetcher = _fetcher(_routes(PROFILES_V2, ALLELES_V2))
        plan = U.check(dbdir, fetcher=fetcher)
        U.apply(dbdir, plan, None, fetcher=fetcher, tools=_FakeTools(), backup=True)
        assert _snapshot(dbdir)["pubmlst/tiny/tiny.txt"] == PROFILES_V2
        restored = U.rollback(dbdir)
        assert _snapshot(dbdir)["pubmlst/tiny/tiny.txt"] == original
        assert restored == "2026-01-01"
        assert U.db_version(dbdir) == "2026-01-01"


def test_resume_skips_an_already_committed_scheme():
    with tempfile.TemporaryDirectory() as tmp:
        dbdir = _tiny_db(tmp)
        fetcher = _fetcher(_routes(PROFILES_V2, ALLELES_V2))
        plan = U.check(dbdir, fetcher=fetcher)
        U.apply(dbdir, plan, None, fetcher=fetcher, tools=_FakeTools())
        calls = len(fetcher._opener.calls)
        U.apply(dbdir, plan, None, fetcher=fetcher, tools=_FakeTools(), resume=True)
        assert len(fetcher._opener.calls) == calls      # nothing re-fetched
        skips = [r for r in U.read_journal(dbdir) if r["event"] == "scheme_skip"]
        assert skips and skips[-1]["scheme"] == "tiny"


def test_cancellation_between_commits_raises_and_journals():
    class Flag:
        def __init__(self):
            self.armed = False

        def is_set(self):
            return self.armed

    with tempfile.TemporaryDirectory() as tmp:
        dbdir = _tiny_db(tmp)
        flag = Flag()
        fetcher = _fetcher(_routes(PROFILES_V2, ALLELES_V2))
        plan = U.check(dbdir, fetcher=fetcher)
        flag.armed = True
        with pytest.raises(Cancelled):
            U.apply(dbdir, plan, None, fetcher=fetcher, tools=_FakeTools(), cancel=flag)
        ends = [r for r in U.read_journal(dbdir) if r["event"] == "run_end"]
        assert ends[-1]["status"] == "cancelled"


def test_progress_reports_a_fraction_and_a_message():
    with tempfile.TemporaryDirectory() as tmp:
        dbdir = _tiny_db(tmp)
        seen = []
        fetcher = _fetcher(_routes(PROFILES_V2, ALLELES_V2))
        U.check(dbdir, fetcher=fetcher, progress=lambda f, m: seen.append((f, m)))
        assert seen and all(0.0 <= f <= 1.0 and isinstance(m, str) for f, m in seen)
        assert seen[-1][0] == 1.0


def test_progress_also_supports_the_three_argument_gui_hook():
    """The GUI task runner's hook is progress(done, total, text) (4.4, 10.5)."""
    seen = []

    def gui_progress(done=0, total=0, text=""):
        seen.append((done, total, text))

    with tempfile.TemporaryDirectory() as tmp:
        dbdir = _tiny_db(tmp)
        U.check(dbdir, fetcher=_fetcher(_routes(PROFILES_V1, ALLELES_V1)),
                progress=gui_progress)
    assert seen and all(total == 1.0 and isinstance(text, str)
                        for _done, total, text in seen)
    assert seen[-1][0] == 1.0 and seen[-1][2]


def test_a_broken_progress_callback_never_breaks_an_update():
    def explode(fraction, message):
        raise RuntimeError("the GUI went away")
    with tempfile.TemporaryDirectory() as tmp:
        dbdir = _tiny_db(tmp)
        plan = U.check(dbdir, fetcher=_fetcher(_routes(PROFILES_V1, ALLELES_V1)),
                       progress=explode)
        assert plan.updates


def test_alias_pairs_are_fetched_once_and_materialised_twice():
    """7.3: cdiphtheriae and diphtheria_3 hold identical payloads."""
    with tempfile.TemporaryDirectory() as tmp:
        dbdir = _tiny_db(tmp)
        root = os.path.join(dbdir, "pubmlst")
        shutil.copytree(os.path.join(root, "tiny"), os.path.join(root, "tinyalias"))
        os.rename(os.path.join(root, "tinyalias", "tiny.txt"),
                  os.path.join(root, "tinyalias", "tinyalias.txt"))
        os.rename(os.path.join(root, "tinyalias", "tiny_info.json"),
                  os.path.join(root, "tinyalias", "tinyalias_info.json"))
        info = U.scheme_info(dbdir, "tinyalias")
        info["name"] = "tinyalias"
        U.write_scheme_info(os.path.join(root, "tinyalias", "tinyalias_info.json"), info)
        refs = U.load_manifest(dbdir)
        assert {r.name: r.alias_of for r in refs} == \
            {"tiny": None, "tinyalias": "tiny"}
        fetcher = _fetcher(_routes(PROFILES_V2, ALLELES_V2))
        plan = U.check(dbdir, fetcher=fetcher)
        U.apply(dbdir, plan, ["tiny", "tinyalias"], fetcher=fetcher, tools=_FakeTools())
        assert fetcher._opener.calls.count(API + "/profiles_csv") == 1
        for name in ("tiny", "tinyalias"):
            with open(os.path.join(root, name, name + ".txt"), "rb") as fh:
                assert fh.read() == PROFILES_V2
            assert U.scheme_info(dbdir, name)["name"] == name


def test_payload_files_are_written_in_binary_mode():
    """6.1: the text layer must never be able to insert a \\r on Windows."""
    with tempfile.TemporaryDirectory() as tmp:
        dbdir = _tiny_db(tmp)
        fetcher = _fetcher(_routes(PROFILES_V2, ALLELES_V2))
        plan = U.check(dbdir, fetcher=fetcher)
        U.apply(dbdir, plan, None, fetcher=fetcher, tools=_FakeTools())
        for path in glob.glob(os.path.join(dbdir, "pubmlst", "tiny", "*")):
            with open(path, "rb") as fh:
                assert b"\r" not in fh.read(), path


def test_write_version_files_is_opt_in():
    """D12: --info's DATE column stays 'Unknown' by default."""
    with tempfile.TemporaryDirectory() as tmp:
        dbdir = _tiny_db(tmp)
        fetcher = _fetcher(_routes(PROFILES_V2, ALLELES_V2, last_updated="2026-02-02"))
        plan = U.check(dbdir, fetcher=fetcher)
        U.apply(dbdir, plan, None, fetcher=fetcher, tools=_FakeTools())
        stamp = os.path.join(dbdir, "pubmlst", "tiny", "database_version.txt")
        assert not os.path.exists(stamp)
    with tempfile.TemporaryDirectory() as tmp:
        dbdir = _tiny_db(tmp)
        fetcher = _fetcher(_routes(PROFILES_V2, ALLELES_V2, last_updated="2026-02-02"))
        plan = U.check(dbdir, fetcher=fetcher)
        U.apply(dbdir, plan, None, fetcher=fetcher, tools=_FakeTools(),
                write_version_files=True)
        with open(os.path.join(dbdir, "pubmlst", "tiny",
                               "database_version.txt")) as fh:
            assert fh.read() == "2026-02-02\n"


# ---------------------------------------------------------------------------
# 7.3 crash repair
# ---------------------------------------------------------------------------

def test_sweep_repairs_a_half_finished_swap():
    with tempfile.TemporaryDirectory() as tmp:
        dbdir = _tiny_db(tmp)
        root = os.path.join(dbdir, "pubmlst")
        os.rename(os.path.join(root, "tiny"), os.path.join(root, "tiny.old.4242"))
        assert not os.path.isdir(os.path.join(root, "tiny"))
        U.sweep_staging(dbdir)
        assert os.path.isdir(os.path.join(root, "tiny"))
        assert not os.path.isdir(os.path.join(root, "tiny.old.4242"))


def test_sweep_removes_an_orphaned_old_copy_and_staging_dir():
    with tempfile.TemporaryDirectory() as tmp:
        dbdir = _tiny_db(tmp)
        root = os.path.join(dbdir, "pubmlst")
        shutil.copytree(os.path.join(root, "tiny"), os.path.join(root, "tiny.old.99"))
        os.makedirs(os.path.join(root, U.STAGING_DIRNAME, "tiny.1234"))
        U.sweep_staging(dbdir)
        assert os.path.isdir(os.path.join(root, "tiny"))
        assert not os.path.isdir(os.path.join(root, "tiny.old.99"))
        assert os.listdir(os.path.join(root, U.STAGING_DIRNAME)) == []


def test_a_staging_directory_is_invisible_to_the_scheme_list():
    """6.1: dot-prefixed entries are skipped, so .staging never becomes a scheme."""
    with tempfile.TemporaryDirectory() as tmp:
        dbdir = _tiny_db(tmp)
        os.makedirs(os.path.join(dbdir, "pubmlst", U.STAGING_DIRNAME, "x"))
        assert U._scheme_names_on_disk(dbdir) == ("tiny",)


# ---------------------------------------------------------------------------
# 4.9 build_blast_db / index_is_stale
# ---------------------------------------------------------------------------

class _FakeTools:
    makeblastdb = shutil.which("makeblastdb") or "makeblastdb"


def test_concat_fasta_rewrites_headers_and_orders_deterministically():
    with tempfile.TemporaryDirectory() as tmp:
        dbdir = _tiny_db(tmp)
        out = os.path.join(tmp, "mlst.fa")
        assert U.concat_fasta(dbdir, out) == 3
        with open(out, "rb") as fh:
            data = fh.read()
        assert data.startswith(b">tiny.abc_1\n")
        assert [line for line in data.split(b"\n") if line.startswith(b">")] == [
            b">tiny.abc_1", b">tiny.def_1", b">tiny.def_2"]


def test_concat_fasta_drops_the_whole_not_a_locus_record():
    """D4: upstream's grep -v orphans the sequence lines onto the previous record."""
    with tempfile.TemporaryDirectory() as tmp:
        dbdir = _tiny_db(tmp)
        with open(os.path.join(dbdir, "pubmlst", "tiny", "abc.tfa"), "wb") as fh:
            fh.write(b">abc_1\nACGT\n>abc not a locus\nNNNNNNNN\n>abc_2\nACGA\n")
        out = os.path.join(tmp, "mlst.fa")
        assert U.concat_fasta(dbdir, out) == 4
        with open(out, "rb") as fh:
            data = fh.read()
        assert b"not a locus" not in data
        assert b"NNNNNNNN" not in data
        assert b">tiny.abc_2\nACGA\n" in data


def test_concat_fasta_matches_the_bundled_record_set():
    """The shipped mlst.fa was built by the shell script; the record SET must match.

    Order differs by design: the shell uses readdir + LC_COLLATE, WMLST uses byte
    order (4.9). Verified out-of-band to give byte-identical blastn output.
    """
    shipped = os.path.join(BUNDLED_DB, "blast", "mlst.fa")
    if not os.path.isfile(shipped):
        pytest.skip("db/blast/mlst.fa has not been built")
    with tempfile.TemporaryDirectory() as tmp:
        out = os.path.join(tmp, "mlst.fa")
        records = U.concat_fasta(BUNDLED_DB, out)
        assert records == 225777
        assert os.path.getsize(out) == os.path.getsize(shipped)
        mine = sorted(line for line in open(out, "rb") if line.startswith(b">"))
        theirs = sorted(line for line in open(shipped, "rb") if line.startswith(b">"))
        assert mine == theirs
        assert len(set(mine)) == len(mine)          # 0 duplicate seqids


def test_index_is_stale_notices_a_missing_and_an_outdated_index():
    with tempfile.TemporaryDirectory() as tmp:
        dbdir = _tiny_db(tmp)
        assert U.index_is_stale(dbdir) is True
        nsq = os.path.join(dbdir, "blast", "mlst.fa.nsq")
        os.makedirs(os.path.dirname(nsq))
        with open(nsq, "wb") as fh:
            fh.write(b"x")
        assert U.index_is_stale(dbdir) is False
        tfa = os.path.join(dbdir, "pubmlst", "tiny", "abc.tfa")
        os.utime(tfa, (os.path.getmtime(nsq) + 10, os.path.getmtime(nsq) + 10))
        assert U.index_is_stale(dbdir) is True


@pytest.mark.skipif(shutil.which("makeblastdb") is None,
                    reason="makeblastdb is not on PATH")
def test_build_blast_db_produces_a_complete_v5_index():
    with tempfile.TemporaryDirectory() as tmp:
        dbdir = _tiny_db(tmp)
        fasta = U.build_blast_db(dbdir, _FakeTools())
        assert fasta == os.path.join(dbdir, "blast", "mlst.fa")
        for ext in U.BLAST_INDEX_EXTENSIONS:
            assert os.path.isfile(fasta + "." + ext), ext
        assert not os.path.isdir(os.path.join(dbdir, "blast", U.STAGING_DIRNAME))
        assert U.index_is_stale(dbdir) is False


@pytest.mark.skipif(shutil.which("makeblastdb") is None,
                    reason="makeblastdb is not on PATH")
def test_a_failed_makeblastdb_leaves_no_partial_index():
    class Broken:
        makeblastdb = shutil.which("makeblastdb")
    with tempfile.TemporaryDirectory() as tmp:
        dbdir = _tiny_db(tmp)
        with open(os.path.join(dbdir, "pubmlst", "tiny", "abc.tfa"), "wb") as fh:
            fh.write(b">abc_1\n" + b"!" * 10 + b"\n")   # not nucleotide data
        blastbin = sys.modules["wmlst.blastbin"]
        original = blastbin.run_tool

        class Result:
            returncode = 1
            stdout = ""
            stderr = "synthetic failure"

        blastbin.run_tool = lambda *a, **k: Result()
        try:
            with pytest.raises(UpdateError):
                U.build_blast_db(dbdir, Broken())
        finally:
            blastbin.run_tool = original
        assert not os.path.isfile(os.path.join(dbdir, "blast", "mlst.fa.nsq"))


# ---------------------------------------------------------------------------
# 4.9 export_bundle / import_bundle
# ---------------------------------------------------------------------------

def test_bundle_round_trip_restores_every_byte():
    with tempfile.TemporaryDirectory() as tmp:
        dbdir = _tiny_db(tmp)
        before = _snapshot(dbdir)
        bundle = os.path.join(tmp, "wmlst-db.zip")
        U.export_bundle(dbdir, bundle)
        target = os.path.join(tmp, "target", "db")
        os.makedirs(os.path.join(target, "pubmlst"))
        version = U.import_bundle(target, bundle)
        assert version == U.db_version(dbdir)
        assert _snapshot(target) == before


def test_import_rejects_a_zip_that_is_not_a_bundle():
    import zipfile
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "x.zip")
        with zipfile.ZipFile(path, "w") as zf:
            zf.writestr("hello.txt", "hi")
        with pytest.raises(UpdateError):
            U.import_bundle(os.path.join(tmp, "db"), path)


def test_import_rejects_a_traversal_path():
    import zipfile
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "x.zip")
        with zipfile.ZipFile(path, "w") as zf:
            zf.writestr(U._BUNDLE_MANIFEST, json.dumps({"db_version": "x"}))
            zf.writestr("../escape.txt", "nope")
        with pytest.raises(UpdateError):
            U.import_bundle(os.path.join(tmp, "db"), path)


# ---------------------------------------------------------------------------
# live tests (skipped when the APIs are unreachable)
# ---------------------------------------------------------------------------

_ONLINE = None


def _online():
    global _ONLINE
    if _ONLINE is None:
        if os.environ.get("WMLST_SKIP_NET"):
            _ONLINE = False
            return _ONLINE
        import urllib.request
        try:
            request = urllib.request.Request(
                U.REST_ROOTS["pubmlst"] + "/db/pubmlst_mgenitalium_seqdef/schemes/2",
                headers={"User-Agent": U.USER_AGENT})
            with urllib.request.urlopen(request, timeout=10):
                _ONLINE = True
        except Exception:
            _ONLINE = False
    return _ONLINE


live = pytest.mark.skipif(not _online(), reason="PubMLST/Pasteur are unreachable")


@live
def test_live_both_hosts_answer_a_scheme_document():
    fetcher = U._Fetcher()
    for url in ("%s/db/pubmlst_mgenitalium_seqdef/schemes/2" % U.REST_ROOTS["pubmlst"],
                "%s/db/pubmlst_kingella_seqdef/schemes/1" % U.REST_ROOTS["pasteur"]):
        meta = fetcher.get_json(url)
        assert isinstance(meta.get("locus_count"), int)
        assert meta.get("loci")


@live
@pytest.mark.skipif(shutil.which("makeblastdb") is None,
                    reason="makeblastdb is not on PATH")
def test_live_download_of_a_small_scheme_matches_the_bundled_bytes():
    """7.1/6.1: the anonymous view IS the bundled database, byte for byte."""
    with tempfile.TemporaryDirectory() as tmp:
        dbdir = os.path.join(tmp, "db")
        os.makedirs(os.path.join(dbdir, "pubmlst"))
        shutil.copy(os.path.join(BUNDLED_DB, "VERSION.txt"),
                    os.path.join(dbdir, "VERSION.txt"))
        refs = tuple(r for r in U.load_manifest(BUNDLED_DB) if r.name == "kingella")
        plan = U.check(dbdir, refs)
        assert [u.status for u in plan.updates] == ["changed"]
        U.apply(dbdir, plan, ["kingella"], refs=refs, tools=_FakeTools(),
                force_blast=False)
        for name in os.listdir(os.path.join(BUNDLED_SCHEMES, "kingella")):
            if name.endswith("_info.json"):
                continue
            with open(os.path.join(BUNDLED_SCHEMES, "kingella", name), "rb") as fh:
                expected = fh.read()
            with open(os.path.join(dbdir, "pubmlst", "kingella", name), "rb") as fh:
                assert fh.read() == expected, name
        assert U.local_content_sha256(dbdir, "kingella") == \
            U.local_content_sha256(BUNDLED_DB, "kingella")


@live
@pytest.mark.skipif(shutil.which("makeblastdb") is None,
                    reason="makeblastdb is not on PATH")
def test_live_a_second_pass_finds_nothing_to_do():
    """7.2: the hash oracle, not the date, decides. A no-op costs no payload."""
    with tempfile.TemporaryDirectory() as tmp:
        dbdir = os.path.join(tmp, "db")
        os.makedirs(os.path.join(dbdir, "pubmlst"))
        shutil.copytree(os.path.join(BUNDLED_SCHEMES, "kingella"),
                        os.path.join(dbdir, "pubmlst", "kingella"))
        refs = tuple(r for r in U.load_manifest(BUNDLED_DB) if r.name == "kingella")
        plan = U.check(dbdir, refs)
        U.apply(dbdir, plan, ["kingella"], refs=refs, tools=_FakeTools(),
                force_blast=False)
        commits = [r for r in U.read_journal(dbdir) if r["event"] == "scheme_commit"]
        assert commits[-1]["status"] == "date-moved-content-identical"
        plan2 = U.check(dbdir, refs)
        assert [u.status for u in plan2.updates] == ["up-to-date"]


@live
def test_live_every_manifest_row_resolves_upstream():
    """6.3: absence from /db is never a reason to prune; only a 404 is."""
    if not os.environ.get("WMLST_FULL_NET_TESTS"):
        pytest.skip("set WMLST_FULL_NET_TESTS=1 for the full 161-request pass")
    refs = U.load_manifest(BUNDLED_DB)
    verified = U.build_manifest(BUNDLED_DB, verify=True)
    assert len(verified) == len(refs)
    by_name = {r.name: r for r in refs}
    for ref in verified:
        assert ref.nloci == by_name[ref.name].nloci, ref.name


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def _snapshot(dbdir):
    """-> ``{relative path: bytes}`` for everything under ``pubmlst/``."""
    out = {}
    root = os.path.join(dbdir, "pubmlst")
    for base, _dirs, files in os.walk(root):
        for name in files:
            path = os.path.join(base, name)
            rel = os.path.relpath(path, dbdir).replace(os.sep, "/")
            with open(path, "rb") as fh:
                out[rel] = fh.read()
    return out


def _install_fake_blastbin():
    """Provide wmlst.blastbin when the real one is not in the tree yet.

    The updater is not allowed to import subprocess; a test harness is.
    """
    try:
        import wmlst.blastbin  # noqa: F401
        return
    except ImportError:
        pass
    import subprocess
    import types
    module = types.ModuleType("wmlst.blastbin")

    def run_tool(argv, *, cwd=None, env=None, timeout=None, cancel=None):
        return subprocess.run(argv, cwd=cwd, stdin=subprocess.DEVNULL,
                              capture_output=True, text=True)

    module.run_tool = run_tool
    module.find_blast = lambda explicit=None: _FakeTools()
    sys.modules["wmlst.blastbin"] = module


_install_fake_blastbin()


if __name__ == "__main__":
    raise SystemExit(pytest.main([os.path.abspath(__file__), "-q"]))
