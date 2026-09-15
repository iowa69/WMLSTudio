"""A catalogued scheme is an identity plus a licence verdict, never a free number.

Nothing here contacts a network service. The pins are compared against themselves
and against threshold_guidance; studio_scripts/stage_cgmlst_schemes.py --verify is
the tool that re-reads the live services.
"""

import importlib.util
import json
from pathlib import Path

import pytest

from wmlstudio import cgmlst_schemes
from wmlstudio.cgmlst_schemes import (
    CGMLST_TARGET_FLOOR,
    LIBRARY_DIRNAME,
    MIGRATIONS_FILENAME,
    PROVIDERS,
    TARGET_SET_CORE,
    SchemeCatalogError,
    bundled_entries,
    catalog_digest,
    catalog_entries,
    classical_root,
    download_only_entries,
    download_plan,
    entry_for,
    entry_for_source,
    install_folder,
    install_root,
    installed_entries,
    installed_scheme,
    library_root,
    library_status,
    migrate_downloads,
    prepare_library,
    resolve_migrated_path,
    scheme_variants,
    slot_for,
    threshold_citations,
    threshold_for,
)
from wmlstudio.threshold_guidance import catalog_entries as publication_entries
from wmlstudio.threshold_guidance import suggested_threshold


def install(root, key, *, loci=None, api=None):
    """A minimal but real installed scheme folder in the slot the catalogue names."""
    entry = entry_for(key)
    folder = library_root(root) / entry["slot"]
    folder.mkdir(parents=True, exist_ok=True)
    for index in range(loci if loci is not None else entry["locus_count"]):
        (folder / f"locus{index:05d}.fasta").write_text(">1\nACGT\n")
    (folder / "scheme.json").write_text(json.dumps(
        {"name": entry["scheme_name"], "type": "cgMLST",
         "API": entry["source_url"] if api is None else api}))
    return folder


def legacy_download(root, key, *, loci=None, digest="abcdef0123456789"):
    """A scheme installed the way releases up to 0.3.0 installed one.

    Every download landed in <data root>/schemes under a digest-shaped folder name,
    which is exactly the reported bug: the scheme is on disk and nobody can find it.
    """
    entry = entry_for(key)
    slug = entry["scheme_id"] if entry["provider"] == "cgmlst.org" else entry["database"]
    folder = classical_root(root) / f"cgmlst_org_{slug}_{digest}"
    folder.mkdir(parents=True, exist_ok=True)
    count = loci if loci is not None else entry["locus_count"]
    for index in range(count):
        (folder / f"t{index:05d}.fasta").write_text(">1\nACGT\n")
    (folder / "scheme.json").write_text(json.dumps(
        {"name": f"{entry['organism']} cgMLST", "organism": entry["organism"],
         "type": "cgMLST", "source": "cgMLST.org", "API": entry["source_url"],
         "locus_count": count}))
    (folder / "reference_manifest.json").write_text(json.dumps(
        {"format_version": 1, "provider": "cgMLST.org", "scheme_digest": digest * 4}))
    return folder


def test_every_pinned_scheme_carries_a_complete_verifiable_identity():
    entries = catalog_entries()
    assert len(entries) >= 20
    assert len({entry["key"] for entry in entries}) == len(entries)
    assert len({entry["slot"] for entry in entries}) == len(entries)
    for entry in entries:
        assert entry["provider"] in PROVIDERS
        assert len(entry["target_list_sha256"]) == 64
        assert entry["locus_count"] > CGMLST_TARGET_FLOOR
        assert entry["kind"] == "cgmlst"
        assert entry["source_url"].startswith("https://")
        assert entry["binding_basis"].strip()
    # The catalogue digest is what a staging manifest pins; it must be stable and
    # must move when any pinned field moves.
    assert catalog_digest() == catalog_digest()
    assert len(catalog_digest()) == 64


def test_only_a_provider_that_grants_redistribution_may_be_packed():
    packable = {entry["provider"] for entry in bundled_entries()}
    assert packable == {"pubmlst"}
    assert PROVIDERS["pubmlst"]["may_bundle"] is True
    # Every refusal is recorded with the sentence it rests on, so the verdict can
    # be re-checked rather than taken on trust.
    for name in ("pasteur", "cgmlst.org", "enterobase", "chewie-ns"):
        provider = PROVIDERS[name]
        assert provider["may_bundle"] is False
        assert provider["restriction"].strip()
        assert provider["terms_url"].startswith("https://")
    assert {entry["provider"] for entry in download_only_entries()} == {"pasteur", "cgmlst.org"}
    assert all(entry["requires_terms_acknowledgement"] for entry in download_only_entries())


def test_threshold_is_offerable_only_when_scheme_and_full_target_count_match():
    salmonella = threshold_for("pubmlst:senterica-cgmlst-3002", locus_count=3002)
    assert salmonella["status"] == "threshold_offerable"
    assert salmonella["threshold"] == 10
    assert salmonella["entries"][0]["source"]["doi"] == "10.3389/fmicb.2023.1254777"
    # One missing target is a different target set, so the cutoff is withdrawn -
    # never scaled, never rounded.
    partial = threshold_for("pubmlst:senterica-cgmlst-3002", locus_count=3001)
    assert partial["status"] == "target_count_mismatch"
    assert partial["threshold"] is None
    assert "3001" in partial["reason"] and "3002" in partial["reason"]
    assert partial["entries"], "the citation stays visible even when the number is withheld"


def test_equal_target_counts_from_different_providers_never_share_a_threshold():
    # PubMLST and cgMLST.org both publish a 2,692-target S. marcescens scheme and
    # they share no target name at all. Only the one the publication names is bound.
    ridom = threshold_for("cgmlst.org:smarcescens-2692", locus_count=2692)
    pubmlst = threshold_for("pubmlst:smarcescens-cgmlst-2692", locus_count=2692)
    assert entry_for("cgmlst.org:smarcescens-2692")["locus_count"] == \
        entry_for("pubmlst:smarcescens-cgmlst-2692")["locus_count"]
    assert entry_for("cgmlst.org:smarcescens-2692")["target_list_sha256"] != \
        entry_for("pubmlst:smarcescens-cgmlst-2692")["target_list_sha256"]
    assert ridom["status"] == "threshold_offerable"
    assert ridom["threshold"] == 12
    assert pubmlst["status"] == "no_bound_cutoff"
    assert pubmlst["threshold"] is None
    assert "coincidence" in pubmlst["reason"]


def test_a_bound_scheme_with_no_curated_number_returns_the_citation_and_no_number():
    guidance = threshold_for("cgmlst.org:paeruginosa-3867", locus_count=3867)
    assert guidance["status"] == "citation_only"
    assert guidance["threshold"] is None
    assert guidance["entries"] and guidance["entries"][0]["published_threshold"] is None
    citations = threshold_citations("cgmlst.org:paeruginosa-3867")
    assert citations[0]["doi"] == "10.1128/jcm.01987-20"
    assert citations[0]["published_threshold"] is None
    assert citations[0]["scheme_key"] == "cgmlst.org:paeruginosa-3867"


def test_a_superseded_scheme_version_does_not_inherit_the_old_cutoff():
    # Bletz 2018 published 6 alleles on the 2,270-target scheme; the provider now
    # serves 2,147 targets. The number must not follow the organism name across.
    guidance = threshold_for("cgmlst.org:cdifficile-2147", locus_count=2147)
    assert guidance["status"] == "no_bound_cutoff"
    assert guidance["threshold"] is None
    assert "2,270" in guidance["reason"]
    assert any(entry["scheme_key"] == "cgmlst.org:cdifficile-2270"
               for entry in publication_entries())


def test_every_bound_scheme_key_exists_in_the_publication_catalogue():
    published = {entry["scheme_key"] for entry in publication_entries()}
    bound = {entry["threshold_scheme_key"] for entry in catalog_entries()
             if entry["threshold_scheme_key"]}
    assert bound
    assert bound <= published, sorted(bound - published)
    for entry in catalog_entries():
        if not entry["threshold_scheme_key"]:
            continue
        rows = [row for row in publication_entries()
                if row["scheme_key"] == entry["threshold_scheme_key"]
                and row["locus_count"] is not None]
        # Where the publication pins a target count it must be this scheme's count.
        assert all(row["locus_count"] == entry["locus_count"] for row in rows), entry["key"]


def test_an_mlst_question_is_refused_rather_than_answered_with_a_cgmlst_number():
    with pytest.raises(SchemeCatalogError, match="different quantity"):
        threshold_for("cgmlst.org:saureus-1861", method="mlst")
    with pytest.raises(SchemeCatalogError, match="Unknown cgMLST scheme"):
        threshold_for("cgmlst.org:not-a-scheme")


def test_library_layout_is_created_with_a_labelled_folder_for_every_scheme(tmp_path):
    report = prepare_library(tmp_path)
    base = library_root(tmp_path)
    assert Path(report["root"]) == base
    assert len(report["created"]) == len(catalog_entries())
    assert (base / "README.txt").is_file()
    for entry in catalog_entries():
        folder = base / entry["slot"]
        assert folder.is_dir()
        readme = (folder / "README.txt").read_text(encoding="utf-8")
        assert entry["organism"] in readme
        assert str(entry["locus_count"]) in readme
        assert entry["terms_url"] in readme
        slot = json.loads((folder / "scheme_slot.json").read_text(encoding="utf-8"))
        assert slot["key"] == entry["key"]
        assert slot["bundled"] is entry["bundled"]
        assert slot["target_list_sha256"] == entry["target_list_sha256"]
    # A folder that cannot be packed says so, in words, before any download starts.
    ridom = (base / entry_for("cgmlst.org:saureus-1861")["slot"] / "README.txt").read_text()
    assert "CANNOT be packed" in ridom
    assert "does not grant" in ridom


def test_preparing_twice_is_idempotent_and_never_touches_installed_data(tmp_path):
    prepare_library(tmp_path)
    folder = install(tmp_path, "cgmlst.org:efaecium-1423", loci=3)
    (folder / "locus00000.fasta").write_text(">1\nACGTACGT\n")
    before = (folder / "locus00000.fasta").read_bytes()
    second = prepare_library(tmp_path)
    assert second["created"] == []
    # The only thing a second run rewrites is the slot's own description, and only
    # because the slot is no longer empty: no allele file and no manifest is read,
    # written or moved.
    assert second["refreshed"] == ["Enterococcus_faecium__cgmlst_org_1423/README.txt"]
    assert "IS installed in this folder" in (folder / "README.txt").read_text(encoding="utf-8")
    assert (folder / "locus00000.fasta").read_bytes() == before
    assert prepare_library(tmp_path)["refreshed"] == []


def test_an_installed_scheme_is_recognised_only_when_its_identity_matches(tmp_path):
    prepare_library(tmp_path)
    assert installed_scheme(tmp_path, "cgmlst.org:efaecalis-1972") is None
    install(tmp_path, "cgmlst.org:efaecalis-1972", loci=1972)
    found = installed_scheme(tmp_path, "cgmlst.org:efaecalis-1972")
    assert found["count_matched"] is True and found["identity_matched"] is True
    # The right number of files with somebody else's identity is not this scheme.
    other = install(tmp_path, "cgmlst.org:abaumannii-2390", loci=2390,
                    api="https://www.cgmlst.org/ncs/schema/Efaecalis/")
    mismatch = installed_scheme(tmp_path, "cgmlst.org:abaumannii-2390")
    assert mismatch["count_matched"] is True
    assert mismatch["identity_matched"] is False
    assert Path(mismatch["path"]) == other


def test_library_status_withholds_a_cutoff_from_a_short_installed_scheme(tmp_path):
    prepare_library(tmp_path)
    install(tmp_path, "cgmlst.org:kpneumoniae-2358", loci=2358)
    install(tmp_path, "cgmlst.org:efaecium-1423", loci=1400)
    rows = {row["key"]: row for row in library_status(tmp_path)}
    complete = rows["cgmlst.org:kpneumoniae-2358"]
    assert complete["ready"] is True
    assert complete["threshold"]["status"] == "threshold_offerable"
    assert complete["threshold"]["threshold"] == 15
    short = rows["cgmlst.org:efaecium-1423"]
    assert short["ready"] is False
    assert short["threshold"]["status"] == "target_count_mismatch"
    assert short["threshold"]["threshold"] is None
    absent = rows["cgmlst.org:cfreundii-3250"]
    assert absent["installed"] is None
    assert absent["ready"] is False
    # Nothing is installed, so nothing is offerable - but the folder still exists.
    assert Path(absent["folder"]).is_dir()


def test_download_plan_states_the_cost_and_the_terms_before_any_transfer():
    ridom = download_plan("cgmlst.org:saureus-1861")
    assert ridom["requires_acknowledgement"] is True
    assert ridom["may_be_bundled"] is False
    assert "archive" in ridom["method"]
    assert "gigabytes" in ridom["method"]
    assert ridom["terms_url"] == PROVIDERS["cgmlst.org"]["terms_url"]
    pasteur = download_plan("pasteur:lmonocytogenes-1748")
    assert "1748 separate requests" in pasteur["method"]
    assert pasteur["requires_acknowledgement"] is True
    oxford = download_plan("pubmlst:senterica-cgmlst-3002")
    assert oxford["may_be_bundled"] is True
    assert oxford["requires_acknowledgement"] is False


def test_the_threshold_organisms_the_catalogue_cannot_pack_are_named(tmp_path):
    # The organisms the threshold catalogue binds are overwhelmingly served by a
    # provider that forbids redistribution. That is a fact about the licences, and
    # the catalogue must keep stating it rather than quietly shipping fewer schemes.
    bound = [entry for entry in catalog_entries() if entry["threshold_scheme_key"]]
    assert bound
    unpackable = [entry["organism"] for entry in bound if not entry["bundled"]]
    assert "Klebsiella pneumoniae" in unpackable
    assert "Listeria monocytogenes" in unpackable
    assert "Staphylococcus aureus" in unpackable
    packable = [entry["organism"] for entry in bound if entry["bundled"]]
    assert "Salmonella enterica" in packable
    assert "Escherichia coli" in packable


def test_slot_descriptions_keep_mlst_and_cgmlst_apart(tmp_path):
    prepare_library(tmp_path)
    library = (library_root(tmp_path) / "README.txt").read_text(encoding="utf-8")
    assert "seven-locus" in library
    assert "different quantities" in library
    slot = json.loads((library_root(tmp_path) / entry_for("cgmlst.org:efaecium-1423")["slot"]
                       / "scheme_slot.json").read_text(encoding="utf-8"))
    assert slot["kind"] == "cgmlst"
    assert "never shares a scale" in slot["interpretation"]
    assert cgmlst_schemes.CGMLST_TARGET_FLOOR == 30


def staging_module():
    spec = importlib.util.spec_from_file_location(
        "stage_cgmlst_schemes",
        Path(__file__).resolve().parents[1] / "studio_scripts/stage_cgmlst_schemes.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_staging_writes_a_definition_pack_that_contains_no_sequence_data(tmp_path):
    staging = staging_module()
    report = staging.stage_definitions(tmp_path, fetch=False)
    manifest = json.loads(Path(report["root"], "manifest.json").read_text(encoding="utf-8"))
    assert manifest["contains_sequence_data"] is False
    assert manifest["catalog_digest"] == catalog_digest()
    assert len(manifest["slots"]) == len(catalog_entries())
    assert {slot["key"] for slot in manifest["slots"]} == {e["key"] for e in catalog_entries()}
    written = sorted(path.name for path in Path(report["root"]).rglob("*") if path.is_file())
    assert set(written) <= {"README.txt", "scheme_slot.json", "manifest.json"}


def test_staging_refuses_to_pack_a_scheme_its_provider_does_not_license(tmp_path):
    staging = staging_module()
    with pytest.raises(SystemExit) as refused:
        staging.stage_alleles("cgmlst.org:saureus-1861", tmp_path)
    message = str(refused.value)
    assert "may not be staged into a release" in message
    assert "Ridom" in message
    assert list(tmp_path.iterdir()) == []
    with pytest.raises(SystemExit, match="may not be staged"):
        staging.stage_alleles("pasteur:lmonocytogenes-1748", tmp_path)


def test_staged_notices_record_each_provider_verdict_and_its_source(tmp_path):
    text = staging_module().notices()
    for provider in PROVIDERS.values():
        assert provider["name"] in text
        assert provider["terms_url"] in text
    assert "May be packed into a WMLSTudio release: NO" in text
    assert "May be packed into a WMLSTudio release: yes" in text
    assert cgmlst_schemes.CATALOG_VERSION in text


def test_a_scheme_folder_on_disk_resolves_to_its_cutoff_and_citation(tmp_path):
    prepare_library(tmp_path)
    folder = install(tmp_path, "cgmlst.org:kpneumoniae-2358", loci=2358)
    identified = cgmlst_schemes.identify_installed(folder)
    assert identified["key"] == "cgmlst.org:kpneumoniae-2358"
    assert identified["locus_count_on_disk"] == 2358
    guidance = cgmlst_schemes.threshold_for_installed(folder)
    assert guidance["status"] == "threshold_offerable"
    assert guidance["threshold"] == 15
    assert guidance["path"] == str(folder)
    assert guidance["entries"][0]["source"]["doi"] == "10.1128/jcm.00646-25"


def test_an_incomplete_snapshot_on_disk_does_not_inherit_the_published_cutoff(tmp_path):
    prepare_library(tmp_path)
    folder = install(tmp_path, "cgmlst.org:kpneumoniae-2358", loci=2300)
    guidance = cgmlst_schemes.threshold_for_installed(folder)
    assert guidance["status"] == "target_count_mismatch"
    assert guidance["threshold"] is None
    assert "2300" in guidance["reason"]


def test_an_unrecognised_scheme_folder_says_so_instead_of_guessing(tmp_path):
    folder = tmp_path / "somebody_elses_scheme"
    folder.mkdir()
    (folder / "locus.fasta").write_text(">1\nACGT\n")
    (folder / "scheme.json").write_text(json.dumps({"name": "Local scheme", "type": "cgMLST"}))
    assert cgmlst_schemes.identify_installed(folder) is None
    guidance = cgmlst_schemes.threshold_for_installed(folder)
    assert guidance["status"] == "scheme_not_catalogued"
    assert guidance["threshold"] is None
    assert "Distances can still be computed" in guidance["reason"]
    assert cgmlst_schemes.identify_installed(tmp_path / "absent") is None


def test_several_cutoffs_on_one_scheme_are_ranked_the_same_way_the_organism_view_ranks_them():
    # Three publications cite the 1,423-target E. faecium scheme for different
    # questions. This lookup and threshold_guidance.suggested_threshold must never
    # show a person two different numbers for the same comparison.
    guidance = threshold_for("cgmlst.org:efaecium-1423", locus_count=1423)
    assert guidance["status"] == "threshold_offerable"
    organism = suggested_threshold("Enterococcus faecium", "cgmlst",
                                   locus_count=1423, scheme_key="cgmlst.org:efaecium-1423")
    assert guidance["suggestion"]["id"] == organism["suggestion"]["id"]
    assert guidance["threshold"] == organism["suggestion"]["published_threshold"]
    assert guidance["alternatives"], "the disagreeing entries stay visible"
    assert {row["published_threshold"] for row in guidance["entries"]} == {20, 25}
    assert all(row["method"] == "cgmlst" for row in guidance["entries"]), \
        "the SKA SNP cutoff for the same organism is a different quantity"


def test_a_scheme_downloaded_into_the_existing_schemes_folder_is_still_recognised(tmp_path):
    # The online reference dialog installs into <data root>/schemes, not into the
    # catalogue slot. Passing those paths must still bind the scheme and its cutoff.
    prepare_library(tmp_path)
    elsewhere = tmp_path / "schemes" / "cgmlst_org_Saureus_0123456789abcdef"
    elsewhere.mkdir(parents=True)
    for index in range(1861):
        (elsewhere / f"t{index:05d}.fasta").write_text(">1\nACGT\n")
    (elsewhere / "scheme.json").write_text(json.dumps(
        {"name": "Staphylococcus aureus cgMLST", "type": "cgMLST",
         "API": "https://www.cgmlst.org/ncs/schema/Saureus/"}))
    assert installed_scheme(tmp_path, "cgmlst.org:saureus-1861") is None
    found = installed_scheme(tmp_path, "cgmlst.org:saureus-1861", extra_paths=[elsewhere])
    assert found["identity_matched"] is True and found["count_matched"] is True
    row = next(row for row in library_status(tmp_path, extra_paths=[elsewhere])
               if row["key"] == "cgmlst.org:saureus-1861")
    assert row["ready"] is True
    assert row["threshold"]["threshold"] == 24


# ----------------------------------------------------------- one library, not two

def test_a_scheme_downloaded_into_the_old_place_is_moved_to_where_it_is_looked_for(tmp_path):
    # The reported bug: "I can search and download the cgMLST scheme but then I
    # cannot find it and it does not work." It landed in <data root>/schemes under
    # a digest-shaped name, sorted alphabetically between two seven-locus schemes.
    old = legacy_download(tmp_path, "cgmlst.org:kpneumoniae-2358", loci=2358)
    report = prepare_library(tmp_path)
    moved = report["migrated"]
    assert [row["from"] for row in moved] == [str(old)]
    new = library_root(tmp_path) / entry_for("cgmlst.org:kpneumoniae-2358")["slot"]
    assert Path(moved[0]["to"]) == new
    assert not old.exists()
    assert len(list(new.glob("*.fasta"))) == 2358
    # The slot's own description travelled with it rather than being overwritten.
    assert (new / "README.txt").is_file() and (new / "scheme_slot.json").is_file()
    found = installed_scheme(tmp_path, "cgmlst.org:kpneumoniae-2358")
    assert found["identity_matched"] is True and found["count_matched"] is True
    assert Path(found["path"]) == new


def test_a_moved_scheme_keeps_a_trail_back_to_the_path_a_project_stored(tmp_path):
    old = legacy_download(tmp_path, "cgmlst.org:efaecalis-1972", loci=1972)
    prepare_library(tmp_path)
    new = library_root(tmp_path) / entry_for("cgmlst.org:efaecalis-1972")["slot"]
    ledger = json.loads((library_root(tmp_path) / MIGRATIONS_FILENAME).read_text(encoding="utf-8"))
    assert [row["from"] for row in ledger["moved"]] == [str(old)]
    assert resolve_migrated_path(tmp_path, old) == new
    # A path that was never moved is handed back unchanged, so a caller can route
    # every stored scheme_path through this without deciding first.
    assert resolve_migrated_path(tmp_path, tmp_path / "elsewhere") == tmp_path / "elsewhere"


def test_a_classical_scheme_is_never_moved_out_of_the_mlst_library(tmp_path):
    classical = classical_root(tmp_path) / "pubmlst_kpneumoniae_seqdef_1_0123456789abcdef"
    classical.mkdir(parents=True)
    for locus in ("gapA", "infB", "mdh", "pgi", "phoE", "rpoB", "tonB"):
        (classical / f"{locus}.tfa").write_text(">1\nACGT\n")
    (classical / "scheme.json").write_text(json.dumps({"name": "MLST", "type": "MLST"}))
    assert prepare_library(tmp_path)["migrated"] == []
    assert classical.is_dir()
    assert not (library_root(tmp_path) / classical.name).exists()


def test_a_migration_never_overwrites_a_scheme_already_in_the_slot(tmp_path):
    prepare_library(tmp_path)
    kept = install(tmp_path, "cgmlst.org:abaumannii-2390", loci=2390)
    before = sorted(path.name for path in kept.iterdir())
    stale = legacy_download(tmp_path, "cgmlst.org:abaumannii-2390", loci=2390, digest="0" * 16)
    moved = migrate_downloads(tmp_path)
    # The second copy is kept, under its own digest-suffixed folder. Two snapshots
    # of one scheme can differ in their target set, so neither is discarded and
    # neither is merged into the other.
    assert [row["from"] for row in moved] == [str(stale)]
    assert Path(moved[0]["to"]).name == kept.name + "__000000000000"
    assert not stale.exists()
    assert sorted(path.name for path in kept.iterdir()) == before


def test_an_unreadable_folder_in_the_mlst_library_is_left_alone(tmp_path):
    # No allele files and nothing declared: the kind cannot be read, so the folder
    # is not moved on a guess.
    mystery = classical_root(tmp_path) / "mystery"
    mystery.mkdir(parents=True)
    (mystery / "notes.txt").write_text("nothing here")
    assert prepare_library(tmp_path)["migrated"] == []
    assert mystery.is_dir()


# --------------------------------------------------------------- readable naming

def test_a_library_folder_is_named_for_the_organism_and_never_for_a_digest():
    entry = entry_for("cgmlst.org:kpneumoniae-2358")
    assert entry["slot"] == "Klebsiella_pneumoniae__cgmlst_org_2358"
    assert "Klebsiella pneumoniae" in entry["title"]
    assert "2358 targets" in entry["title"]
    assert PROVIDERS["cgmlst.org"]["name"] in entry["title"]
    for row in catalog_entries():
        assert row["organism"] in row["title"]
        assert f"{row['locus_count']} targets" in row["title"]
        assert not any(part in row["slot"] for part in ("  ", "__cgmlst_org_0"))


def test_a_download_is_matched_to_its_pinned_slot_by_the_providers_own_identifier():
    # cgMLST.org names the scheme by slug, PubMLST by database and scheme id. The
    # organism label a catalogue page shows is never what resolves the identity.
    ridom = {"organism": "Klebsiella pneumoniae sensu lato", "provider": "cgMLST.org",
             "slug": "Kpneumoniae_complex", "locus_count": 2358,
             "url": "https://www.cgmlst.org/ncs/schema/Kpneumoniae_complex/"}
    assert entry_for_source(ridom)["key"] == "cgmlst.org:kpneumoniae-2358"
    assert slot_for(ridom) == entry_for("cgmlst.org:kpneumoniae-2358")["slot"]
    oxford = {"organism": "Salmonella spp.", "provider": "PubMLST",
              "database": "pubmlst_salmonella_seqdef", "scheme_id": "4", "locus_count": 3002}
    assert entry_for_source(oxford)["key"] == "pubmlst:senterica-cgmlst-3002"
    pasteur = {"organism": "Listeria monocytogenes", "provider": "BIGSdb-Pasteur",
               "url": "https://bigsdb.pasteur.fr/api/db/pubmlst_listeria_seqdef/schemes/3"}
    assert entry_for_source(pasteur)["key"] == "pasteur:lmonocytogenes-1748"
    # An uncatalogued scheme still gets a readable folder rather than a digest.
    unknown = {"organism": "Vibrio cholerae", "provider": "cgMLST.org", "slug": "Vcholerae",
               "locus_count": 1234, "url": "https://www.cgmlst.org/ncs/schema/Vcholerae/"}
    assert entry_for_source(unknown) is None
    assert slot_for(unknown) == "Vibrio_cholerae__cgmlst_org_1234"


def test_equal_target_counts_from_two_providers_never_resolve_to_one_identity():
    ridom = {"provider": "cgMLST.org", "slug": "Smarcescens", "locus_count": 2692,
             "url": "https://www.cgmlst.org/ncs/schema/Smarcescens/"}
    oxford = {"provider": "PubMLST", "database": "pubmlst_serratia_seqdef", "scheme_id": "2",
              "locus_count": 2692}
    assert entry_for_source(ridom)["key"] == "cgmlst.org:smarcescens-2692"
    assert entry_for_source(oxford)["key"] == "pubmlst:smarcescens-cgmlst-2692"
    assert slot_for(ridom) != slot_for(oxford)


def test_a_redefined_scheme_of_the_same_size_gets_its_own_folder(tmp_path):
    prepare_library(tmp_path)
    entry = entry_for("cgmlst.org:efaecium-1423")
    slot = library_root(tmp_path) / entry["slot"]
    descriptor = {"organism": entry["organism"], "provider": "cgMLST.org",
                  "slug": entry["scheme_id"], "locus_count": entry["locus_count"],
                  "url": entry["source_url"]}
    # An empty labelled slot IS the destination: that is what the folder is for.
    assert install_folder(library_root(tmp_path), descriptor, digest="a" * 64) == slot
    install(tmp_path, "cgmlst.org:efaecium-1423", loci=1423)
    (slot / "reference_manifest.json").write_text(json.dumps({"scheme_digest": "a" * 64}))
    assert install_folder(library_root(tmp_path), descriptor, digest="a" * 64) == slot
    second = install_folder(library_root(tmp_path), descriptor, digest="b" * 64)
    assert second != slot and second.name.startswith(entry["slot"])
    assert "bbbbbbbbbbbb" in second.name


def test_a_gene_by_gene_scheme_installs_into_the_one_cgmlst_library(tmp_path):
    assert install_root(tmp_path / "schemes") == tmp_path / LIBRARY_DIRNAME
    assert install_root(tmp_path / LIBRARY_DIRNAME) == tmp_path / LIBRARY_DIRNAME
    assert install_root(tmp_path / "anywhere") == tmp_path / "anywhere" / LIBRARY_DIRNAME
    # A classical scheme keeps the library the caller named; existing projects
    # store those exact folder paths.
    assert install_root(tmp_path / "schemes", kind="mlst") == tmp_path / "schemes"


# ------------------------------------------------------------- core vs accessory

def test_every_catalogued_scheme_says_which_target_set_it_is():
    """A row that named no target set would be a scheme of unknown meaning.

    Nothing may be inferred from a target count: a distance is a core-genome
    distance, an accessory one or a pan-genome one because the pin says so.
    """
    rows = catalog_entries()
    assert {row["target_set"] for row in rows} <= set(cgmlst_schemes.TARGET_SETS)
    assert {row["target_set_detail"] for row in rows} <= set(cgmlst_schemes.TARGET_SET_DETAILS)
    assert all(row["target_set"] == TARGET_SET_CORE
               for row in rows if row["target_set_detail"] == TARGET_SET_CORE)
    # Every set that is not the core set falls in the one bucket the rest of the
    # application renders, so no row can be shown as a core set by omission.
    assert all(row["target_set"] == cgmlst_schemes.TARGET_SET_ACCESSORY
               for row in rows if row["target_set_detail"] != TARGET_SET_CORE)
    assert all(row["scheme_group"].startswith(("pubmlst:", "pasteur:", "cgmlst_org:"))
               for row in rows)


def test_a_target_set_the_catalogue_does_not_know_is_refused_rather_than_read_as_core(monkeypatch):
    """Prevents a typo in a future pin from silently becoming a core-genome scheme.

    A row whose target set defaulted to "core" on an unrecognised value would let a
    cutoff published for a core set be offered for something that is not one.
    """
    key = "pubmlst:banthracis-accessory-1263"
    monkeypatch.setattr(cgmlst_schemes, "_SCHEMES", tuple(
        row if row.get("key") != key else {**row, "target_set": "acessory"}
        for row in cgmlst_schemes._SCHEMES))
    with pytest.raises(SchemeCatalogError) as refused:
        catalog_entries()
    assert "acessory" in str(refused.value)
    assert "guessing" in str(refused.value)


def test_an_organism_with_only_a_core_set_says_so_instead_of_offering_an_empty_choice():
    variants = scheme_variants("Klebsiella pneumoniae")
    assert [row["key"] for row in variants["core"]] == ["cgmlst.org:kpneumoniae-2358"]
    assert variants["accessory"] == []
    assert variants["has_core"] is True and variants["has_accessory"] is False
    assert "Only a core target set is catalogued" in variants["message"]
    absent = scheme_variants("Vibrio cholerae")
    assert absent["core"] == [] and absent["accessory"] == []
    assert "gap in this catalogue" in absent["message"]
    # Two providers publish a core Staphylococcus aureus scheme; both are offered,
    # and neither is an accessory set.
    aureus = scheme_variants("Staphylococcus aureus")
    assert {row["provider"] for row in aureus["core"]} == {"cgmlst.org", "pubmlst"}
    assert aureus["has_accessory"] is False


def test_the_anthracis_core_and_accessory_sets_are_both_offered_as_an_explicit_choice():
    """The reported bug: the target-set chooser had nothing to choose between.

    PubMLST publishes 3,803 core targets and 1,263 accessory targets for
    B. anthracis. Without both pinned, a user who wants core-plus-accessory has no
    way to ask for it and no way to see that the two are different runs.
    """
    variants = scheme_variants("Bacillus anthracis")
    assert [row["key"] for row in variants["core"]] == ["pubmlst:banthracis-cgmlst-3803"]
    assert [row["key"] for row in variants["accessory"]] == ["pubmlst:banthracis-accessory-1263"]
    assert variants["has_accessory"] is True and variants["has_whole_genome"] is False
    assert "1 core and 1 accessory target set(s)" in variants["message"]
    assert "never shares a scale, an axis or a threshold" in variants["message"]
    core, accessory = variants["core"][0], variants["accessory"][0]
    # One provider, one database, two schemes: they are told apart by the provider's
    # own scheme id and by the digest of the target list, never by organism or size.
    assert core["database"] == accessory["database"] == "pubmlst_bcereus_seqdef"
    assert (core["scheme_id"], accessory["scheme_id"]) == ("2", "3")
    assert core["target_list_sha256"] != accessory["target_list_sha256"]
    assert core["slot"] != accessory["slot"]
    assert core["scheme_group"] == accessory["scheme_group"] == "pubmlst:bacillus_anthracis"


def test_the_gonococcal_pan_genome_set_is_named_a_whole_genome_set_and_not_an_accessory_one():
    """A 251-target accessory set and a 1,907-target pan-genome set are not one thing.

    Both are offered as the alternative to a core run, but a report that called the
    pan-genome set "accessory" would misdescribe what its distances were measured on.
    """
    variants = scheme_variants("Neisseria gonorrhoeae")
    assert [row["key"] for row in variants["core"]] == ["pubmlst:ngonorrhoeae-cgmlst-1430",
                                                        "pubmlst:ngonorrhoeae-cgmlst-1649"]
    assert [row["key"] for row in variants["accessory_only"]] == \
        ["pubmlst:ngonorrhoeae-accessory-251"]
    assert [row["key"] for row in variants["whole_genome"]] == ["pubmlst:ngonorrhoeae-pgmlst-1907"]
    # The menu of everything that is not the core set holds both of them.
    assert len(variants["accessory"]) == 2
    assert "2 core, 1 accessory and 1 whole-genome target set(s)" in variants["message"]
    pan = entry_for("pubmlst:ngonorrhoeae-pgmlst-1907")
    assert pan["target_set"] == "accessory" and pan["target_set_detail"] == "whole_genome"
    assert cgmlst_schemes.target_set_label(pan) == "Whole genome"
    assert cgmlst_schemes.target_set_label(entry_for("pubmlst:ngonorrhoeae-accessory-251")) == \
        "Accessory"
    assert "whole-genome (pan-genome) target set" in pan["title"]
    assert "1907 targets" in pan["title"]


def test_a_non_core_target_set_never_inherits_a_cutoff_published_for_a_core_set():
    """The core cutoff must not follow the organism name onto a bigger target set.

    Abdel-Glil et al. published five allele differences on the 3,803-target
    B. anthracis core set. Offering that number for a 1,263-target accessory run, or
    for a 1,907-target pan-genome run, would present a threshold as applying to a
    quantity it was never measured on.
    """
    for key in ("pubmlst:banthracis-accessory-1263", "pubmlst:ngonorrhoeae-accessory-251",
                "pubmlst:ngonorrhoeae-pgmlst-1907"):
        entry = entry_for(key)
        assert entry["threshold_scheme_key"] is None
        guidance = threshold_for(key, locus_count=entry["locus_count"])
        assert guidance["status"] == "no_bound_cutoff"
        assert guidance["threshold"] is None
        assert guidance["entries"] == []
        assert guidance["target_set"] == "accessory"
    assert "5,066" in " ".join(entry_for("pubmlst:banthracis-accessory-1263")["notes"])
    assert "different quantity" in entry_for("pubmlst:banthracis-accessory-1263")["binding_basis"]
    # The interpretation sentence a report carries says the same thing in words.
    assert "core-plus-accessory" in cgmlst_schemes.INTERPRETATION


def test_a_publication_with_no_stated_target_count_binds_only_to_a_core_set(monkeypatch):
    """A number published without a target count cannot be claimed for a non-core set.

    It is read as a cutoff on the pinned core set, which is what an unqualified
    cgMLST cutoff means; read the same way for an accessory or pan-genome set it
    would silently move a core threshold onto a different quantity.
    """
    key = "cgmlst.org:efaecalis-1972"
    template = next(row for row in publication_entries()
                    if row["scheme_key"] == key and row["method"] == "cgmlst")
    # Every curated publication states its target count today, so the case is built
    # here rather than waiting for a future entry that does not.
    monkeypatch.setattr(cgmlst_schemes, "_threshold_guidance_entries",
                        lambda: [{**template, "locus_count": None}])
    assert threshold_for(key, locus_count=1972)["status"] == "threshold_offerable"
    declared_accessory = tuple(
        row if row.get("key") != key else {**row, "target_set": cgmlst_schemes.TARGET_SET_ACCESSORY}
        for row in cgmlst_schemes._SCHEMES)
    monkeypatch.setattr(cgmlst_schemes, "_SCHEMES", declared_accessory)
    guidance = threshold_for(key, locus_count=1972)
    assert guidance["status"] == "citation_only"
    assert guidance["threshold"] is None
    assert "not a core genome scheme" in guidance["reason"]


def test_a_slot_for_a_non_core_set_says_in_words_that_it_is_not_a_cgmlst_distance(tmp_path):
    """A folder labelled only "cgMLST" would mislabel an accessory download.

    The slot description is what a person reads months later to find out what the
    numbers in it were measured on.
    """
    prepare_library(tmp_path, keys=["pubmlst:banthracis-accessory-1263",
                                    "pubmlst:ngonorrhoeae-pgmlst-1907",
                                    "pubmlst:banthracis-cgmlst-3803"])
    base = library_root(tmp_path)
    accessory = (base / entry_for("pubmlst:banthracis-accessory-1263")["slot"]
                 / "README.txt").read_text(encoding="utf-8")
    assert "Target set : Accessory" in accessory
    assert "ACCESSORY target set, not a cgMLST scheme" in accessory
    assert "reported as not assayed" in accessory
    pan = (base / entry_for("pubmlst:ngonorrhoeae-pgmlst-1907")["slot"]
           / "README.txt").read_text(encoding="utf-8")
    assert "Target set : Whole genome" in pan
    assert "Run it INSTEAD of a core scheme" in pan
    core = (base / entry_for("pubmlst:banthracis-cgmlst-3803")["slot"]
            / "README.txt").read_text(encoding="utf-8")
    assert "Target set : Core" in core
    assert "ACCESSORY target set" not in core
    slot = json.loads((base / entry_for("pubmlst:banthracis-accessory-1263")["slot"]
                       / "scheme_slot.json").read_text(encoding="utf-8"))
    assert slot["target_set_detail"] == "accessory" and slot["target_set_label"] == "Accessory"


def test_a_download_of_a_non_core_set_says_what_it_is_before_any_byte_moves():
    """The download button is the last moment before the run's meaning is fixed."""
    plan = download_plan("pubmlst:ngonorrhoeae-pgmlst-1907")
    assert plan["target_set_label"] == "Whole genome"
    assert "Run it INSTEAD of a core scheme" in plan["target_set_notice"]
    assert plan["destination_slot"] == "Neisseria_gonorrhoeae__pubmlst_1907"
    assert "1907 separate requests" in plan["method"]
    core = download_plan("pubmlst:ngonorrhoeae-cgmlst-1430")
    assert core["target_set_label"] == "Core"
    assert "not a seven-locus MLST distance" in core["target_set_notice"]
    assert "Run it INSTEAD" not in core["target_set_notice"]
    assert core["destination_slot"] != plan["destination_slot"]


def test_each_pinned_target_set_carries_the_identity_a_download_is_verified_against():
    """A wrong scheme id, database or target count makes a download fail verification.

    These six rows were read from the live PubMLST API on 2026-09-15; the test pins
    what was read so a later edit cannot quietly change one of them.
    """
    expected = {
        "pubmlst:banthracis-cgmlst-3803": ("pubmlst_bcereus_seqdef", "2", 3803, True),
        "pubmlst:banthracis-accessory-1263": ("pubmlst_bcereus_seqdef", "3", 1263, False),
        "pubmlst:ngonorrhoeae-cgmlst-1430": ("pubmlst_neisseria_seqdef", "89", 1430, True),
        "pubmlst:ngonorrhoeae-cgmlst-1649": ("pubmlst_neisseria_seqdef", "62", 1649, True),
        "pubmlst:ngonorrhoeae-accessory-251": ("pubmlst_neisseria_seqdef", "80", 251, False),
        "pubmlst:ngonorrhoeae-pgmlst-1907": ("pubmlst_neisseria_seqdef", "81", 1907, False),
    }
    for key, (database, scheme_id, count, profiles) in expected.items():
        entry = entry_for(key)
        assert (entry["database"], entry["scheme_id"]) == (database, scheme_id)
        assert entry["locus_count"] == count
        assert entry["has_profiles"] is profiles
        assert entry["profile_field"] == ("cgST" if profiles else "")
        assert entry["provider"] == "pubmlst"
        assert entry["source_url"] == \
            f"https://rest.pubmlst.org/db/{database}/schemes/{scheme_id}"
        assert len(entry["target_list_sha256"]) == 64
        assert entry["locus_count"] > CGMLST_TARGET_FLOOR
    # Six new target sets, six distinct pins: nothing is shared by accident.
    digests = {entry_for(key)["target_list_sha256"] for key in expected}
    slots = {entry_for(key)["slot"] for key in expected}
    assert len(digests) == len(slots) == len(expected)


# -------------------------------------------------------- what the cgMLST tab lists

def test_the_cgmlst_tab_lists_what_is_installed_with_a_readable_title(tmp_path):
    prepare_library(tmp_path)
    install(tmp_path, "cgmlst.org:efaecalis-1972", loci=1972)
    rows = installed_entries(tmp_path)
    assert [row["catalog_key"] for row in rows] == ["cgmlst.org:efaecalis-1972"]
    row = rows[0]
    assert row["kind"] == "cgmlst"
    assert "Enterococcus faecalis" in row["title"]
    assert "1972 targets" in row["title"]
    assert row["count_matched"] is True
    assert row["threshold"]["status"] == "threshold_offerable"
    assert row["threshold"]["threshold"] == \
        threshold_for("cgmlst.org:efaecalis-1972", locus_count=1972)["threshold"]
    # An empty labelled slot is a description, not an installed scheme.
    assert all("Klebsiella" not in entry["title"] for entry in rows)


def test_an_uncatalogued_installed_scheme_is_listed_without_a_borrowed_cutoff(tmp_path):
    prepare_library(tmp_path)
    folder = library_root(tmp_path) / "Vibrio_cholerae__cgmlst_org_40"
    folder.mkdir(parents=True)
    for index in range(40):
        (folder / f"t{index:03d}.fasta").write_text(">1\nACGT\n")
    (folder / "scheme.json").write_text(json.dumps(
        {"name": "Vibrio cholerae cgMLST", "organism": "Vibrio cholerae", "type": "cgMLST",
         "source": "cgMLST.org", "API": "https://www.cgmlst.org/ncs/schema/Vcholerae/"}))
    row = next(entry for entry in installed_entries(tmp_path) if "Vibrio" in entry["title"])
    assert row["catalogued"] is False and row["catalog_key"] is None
    assert row["threshold"]["status"] == "scheme_not_catalogued"
    assert row["threshold"]["threshold"] is None
    assert "Distances can still be computed" in row["threshold"]["reason"]
