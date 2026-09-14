"""Pinned cgMLST scheme catalogue, licence verdicts and the on-disk library layout.

A cgMLST scheme is a *target set*, not a file. This module pins the identity of
each scheme WMLSTudio knows about -- provider, scheme identifier, revision, target
count and the SHA-256 of its sorted target list -- so that an installed folder can
be recognised, a published threshold can be bound to it, and an upstream
redefinition is detected instead of silently adopted. No allele sequence is stored
here and none is committed to this repository; studio_scripts/stage_cgmlst_schemes.py
stages the definition pack, and allele data is downloaded on explicit request.

REDISTRIBUTION. Every provider below was read directly on the review date and the
verdict recorded with the sentence it rests on. Only PubMLST (University of Oxford)
grants redistribution, and only for records submitted on or before 2024-12-31 --
which is exactly the subset its unauthenticated API serves. Institut Pasteur,
cgMLST.org/Ridom and EnteroBase each forbid or do not grant redistribution, and
chewie-NS states no terms for its schema data at all; those schemes are never
packed, only offered as a download into a pre-created, clearly-labelled folder.

A 7-locus MLST distance and a 2,000-target cgMLST distance are different
quantities. Nothing in this module will return a cgMLST threshold for an MLST
scheme, and a cutoff is offerable only when the bound scheme key AND the full
published target count both match the installed scheme.
"""

from __future__ import annotations

import hashlib
import json
import re
from copy import deepcopy
from datetime import UTC, datetime
from pathlib import Path

from .threshold_guidance import SOURCES
from .threshold_guidance import catalog_entries as _threshold_guidance_entries

CATALOG_VERSION = "2026-09-14.1"
LICENCE_REVIEWED_ON = "2026-09-14"
LIBRARY_DIRNAME = "cgmlst"
SLOT_FILENAME = "scheme_slot.json"
TARGETS_FILENAME = "targets.txt"
README_FILENAME = "README.txt"
SLOT_FORMAT_VERSION = 1
# The same floor project.classify_typing and reference_index use to separate a
# classical seven-locus scheme from a gene-by-gene one. Repeated, not imported,
# to keep this module free of the storage/typing import chain.
CGMLST_TARGET_FLOOR = 30

INTERPRETATION = (
    "An installed cgMLST scheme is a target set and an allele nomenclature. Two profiles are "
    "comparable only when both were called against the same installed scheme; a distance never "
    "crosses schemes, providers or revisions, and never shares a scale with a seven-locus MLST "
    "distance."
)

PROVIDERS = {
    "pubmlst": {
        "name": "PubMLST (University of Oxford)",
        "terms_url": "https://pubmlst.org/terms-conditions",
        "api_root": "https://rest.pubmlst.org",
        "redistribution": "permitted_for_submissions_up_to_2024-12-31",
        "quote": ("All data submitted on or before 31 December 2024 may continue to be "
                  "downloaded, used, and redistributed without restriction, subject only to "
                  "proper scientific citation and acknowledgement."),
        "restriction": ("Data submitted on or after 1 January 2025 may not be redistributed or "
                        "incorporated into a software product without written authorisation from "
                        "the University of Oxford. Unauthenticated API access is itself limited to "
                        "records submitted on or before 2024-12-31, so an anonymous snapshot "
                        "contains only the redistributable subset -- and is therefore NOT the "
                        "current allele set."),
        "may_bundle": True,
        "requires_acknowledgement": False,
        "reviewed_on": LICENCE_REVIEWED_ON,
    },
    "pasteur": {
        "name": "Institut Pasteur BIGSdb",
        "terms_url": "https://bigsdb.pasteur.fr/policy/",
        "api_root": "https://bigsdb.pasteur.fr/api",
        "redistribution": "prohibited",
        "quote": ("You acknowledge and agree that you may not share or publicly post significant "
                  "parts of the Data (including alleles and profiles definitions)."),
        "restriction": ("All SeqDef data are for academic use only; commercial services and "
                        "software tools that enable use or downloading of the data require a "
                        "specific licence contract from Institut Pasteur."),
        "may_bundle": False,
        "requires_acknowledgement": True,
        "reviewed_on": LICENCE_REVIEWED_ON,
    },
    "cgmlst.org": {
        "name": "cgMLST.org Nomenclature Server (Ridom GmbH)",
        "terms_url": "https://www.cgmlst.org/serverpolicy.html",
        "api_root": "https://www.cgmlst.org/ncs",
        "redistribution": "prohibited_without_permission",
        "quote": ("Reuse of database copies in a product or service requires permission."),
        "restriction": ("cgMLST.org nomenclature belongs to Ridom GmbH. Individual downloads are "
                        "limited to non-commercial use and publication acknowledgement is "
                        "required. Confirm that your intended use is permitted before "
                        "downloading; this software grants no redistribution rights."),
        "may_bundle": False,
        "requires_acknowledgement": True,
        "reviewed_on": LICENCE_REVIEWED_ON,
    },
    "enterobase": {
        "name": "EnteroBase (University of Warwick)",
        "terms_url": "https://enterobase.readthedocs.io/en/latest/enterobase-terms-of-use.html",
        "api_root": "https://enterobase.warwick.ac.uk",
        "redistribution": "not_granted",
        "quote": ("Reverse engineering and Reproducing EnteroBase-like database copies without "
                  "explicit licencing [is not permitted]."),
        "restriction": ("Use is academic; usage outside academic use requires explicit licensing "
                        "from the University of Warwick. Packing the allele database into a "
                        "product is reproducing a database copy, so it is not done here. The "
                        "EnteroBase-derived Salmonella core target set is instead taken from "
                        "PubMLST, which hosts it under terms that do permit redistribution."),
        "may_bundle": False,
        "requires_acknowledgement": True,
        "reviewed_on": LICENCE_REVIEWED_ON,
    },
    "chewie-ns": {
        "name": "Chewie Nomenclature Server (chewbbaca.online)",
        "terms_url": "https://chewbbaca.online/",
        "api_root": "https://chewbbaca.online/NS/api",
        "redistribution": "not_stated",
        "quote": "",
        "restriction": ("The chewie-NS software is GPLv3, but no licence or terms of use is "
                        "stated for the hosted schema data itself. Silence is not permission, so "
                        "no chewie-NS schema is packed and none is catalogued here. A curator who "
                        "obtains written permission for a specific schema can add it as a pinned "
                        "entry."),
        "may_bundle": False,
        "requires_acknowledgement": True,
        "reviewed_on": LICENCE_REVIEWED_ON,
    },
}

# Every row was read from the live service on the review date. locus_count and
# target_list_sha256 are the pins: target_list_sha256 is the SHA-256 of the sorted,
# newline-joined target names, so it survives a change in listing order but not a
# change in the target set. threshold_scheme_key is the threshold_guidance key this
# scheme IS, never the key of a scheme that merely has the same number of targets.
_SCHEMES = (
    # --- Redistributable: PubMLST, unauthenticated (submissions up to 2024-12-31) ---
    {"key": "pubmlst:senterica-cgmlst-3002", "genus": "Salmonella", "species": "enterica",
     "provider": "pubmlst", "database": "pubmlst_salmonella_seqdef", "scheme_id": "4",
     "scheme_name": "cgMLST v2 (Enterobase)", "revision": "read 2026-09-14; PubMLST lists no "
     "last_updated date for this definition-only scheme", "locus_count": 3002,
     "target_list_sha256": "e819029ed5808af4f98d0358ab2da335958f61111cdfa80c3ee6a123629a7c44",
     "has_profiles": False, "profile_field": "",
     "threshold_scheme_key": "enterobase:senterica-cgmlst-3002",
     "binding_basis": ("The sorted target-name set of this PubMLST scheme is byte-identical to the "
                       "cgMLST.org Senterica 3,002-target set, verified by SHA-256 on 2026-09-14, "
                       "and PubMLST itself names it the EnteroBase v2 scheme. Leeper et al. 2023 "
                       "evaluated the EnteroBase-derived 3,002-locus core scheme."),
     "notes": ("PubMLST supplies no cgST profile table for this scheme, so no cgMLST sequence type "
               "is assigned; allelic distances between isolates called against this snapshot are "
               "still computed.",)},
    {"key": "pubmlst:ecoli-cgmlst-2513", "genus": "Escherichia", "species": "coli",
     "provider": "pubmlst", "database": "pubmlst_escherichia_seqdef", "scheme_id": "6",
     "scheme_name": "cgMLST", "revision": "read 2026-09-14; PubMLST lists no last_updated date "
     "for this definition-only scheme", "locus_count": 2513,
     "target_list_sha256": "cadaab19f8de3bfc6b75a34725cb560a3fd38e4d3231030b4ba8bbdb013d9971",
     "has_profiles": False, "profile_field": "",
     "threshold_scheme_key": "cgmlst.org:ecoli-2513",
     "binding_basis": ("The sorted target-name set is byte-identical to the cgMLST.org Ecoli "
                       "2,513-target set, verified by SHA-256 on 2026-09-14, so the target set "
                       "Glasgow et al. 2025 compared is the target set installed here."),
     "notes": ("The target set matches, the allele nomenclature and the caller do not. The "
               "published cutoff was derived with SeqSphere+ on the Ridom allele database; it is "
               "context, not a validated WMLSTudio setting.",
               "This scheme does not separate Escherichia coli from Shigella.")},
    {"key": "pubmlst:saureus-cgmlst-1716", "genus": "Staphylococcus", "species": "aureus",
     "provider": "pubmlst", "database": "pubmlst_saureus_seqdef", "scheme_id": "20",
     "scheme_name": "cgMLST", "revision": "last_updated 2026-09-02", "locus_count": 1716,
     "target_list_sha256": "af49fbade41f5b711444f6e3de34ced882c8ad0c2a0b80ce8a41e49a802ba609",
     "has_profiles": True, "profile_field": "cgST",
     "threshold_scheme_key": None,
     "binding_basis": ("No threshold is bound. The curated S. aureus cutoffs cite the cgMLST.org "
                       "1,861-target scheme; this PubMLST scheme has 1,716 targets, so it is a "
                       "different quantity and no cutoff is offered for it."),
     "notes": ()},
    {"key": "pubmlst:abaumannii-cgmlst-2133", "genus": "Acinetobacter", "species": "baumannii",
     "provider": "pubmlst", "database": "pubmlst_abaumannii_seqdef", "scheme_id": "3",
     "scheme_name": "cgMLST v1", "revision": "last_updated 2026-09-12", "locus_count": 2133,
     "target_list_sha256": "3ecb4e140409efb49dcd03ce610227d6ffae9d67ea3b514a3ee1b6f0a5bae31d",
     "has_profiles": True, "profile_field": "cgST",
     "threshold_scheme_key": None,
     "binding_basis": ("No threshold is bound. The curated A. baumannii cutoff cites the "
                       "cgMLST.org 2,390-target scheme; this scheme has 2,133 targets."),
     "notes": ()},
    {"key": "pubmlst:smarcescens-cgmlst-2692", "genus": "Serratia", "species": "marcescens",
     "provider": "pubmlst", "database": "pubmlst_serratia_seqdef", "scheme_id": "2",
     "scheme_name": "S. marcescens cgMLST", "revision": "last_updated 2026-09-13",
     "locus_count": 2692,
     "target_list_sha256": "a745b4ea68a2d22518d952194b807f9945c3012e307ae44898c714dec3cd006e",
     "has_profiles": True, "profile_field": "cgST",
     "threshold_scheme_key": None,
     "binding_basis": ("No threshold is bound. This scheme has the same NUMBER of targets (2,692) "
                       "as the Kampmeier 2022 cgMLST.org scheme but a completely different target "
                       "set: none of the 2,692 target names is shared (SERR* here, SMDB11_RS* "
                       "there), verified on 2026-09-14. Equal target counts are a coincidence, "
                       "not a scheme identity."),
     "notes": ()},
    {"key": "pubmlst:mtbc-cgmlst-1561", "genus": "Mycobacterium", "species": "",
     "organism": "Mycobacterium tuberculosis complex",
     "provider": "pubmlst", "database": "pubmlst_mycobacteria_seqdef", "scheme_id": "3",
     "scheme_name": "MTBC cgMLST", "revision": "last_updated 2026-09-10", "locus_count": 1561,
     "target_list_sha256": "40cd85297b4b9028771cf5df6f2af15b4a03b9f5a5dfea9313f6c891d09abba8",
     "has_profiles": True, "profile_field": "cgST",
     "threshold_scheme_key": None,
     "binding_basis": ("No threshold is bound. The curated MTBC cutoff cites the cgMLST.org "
                       "2,891-target scheme; this scheme has 1,561 targets. Kohl et al. state "
                       "explicitly that a reduced target set voids the cutoff."),
     "notes": ()},
    {"key": "pubmlst:campylobacter-cgmlst-1142", "genus": "Campylobacter", "species": "jejuni",
     "provider": "pubmlst", "database": "pubmlst_campylobacter_seqdef", "scheme_id": "8",
     "scheme_name": "C. jejuni / C. coli cgMLST v2", "revision": "last_updated 2026-09-05",
     "locus_count": 1142,
     "target_list_sha256": "0391699d0fcf5b416b4b2eb711c76994b14a2491a0327fcb57f5619dd21533ae",
     "has_profiles": True, "profile_field": "cgST",
     "threshold_scheme_key": None,
     "binding_basis": "No cutoff for Campylobacter is curated in the threshold catalogue.",
     "notes": ("The scheme spans C. jejuni and C. coli; the organism label on an isolate is a "
               "separate determination.",)},
    {"key": "pubmlst:nmeningitidis-cgmlst-1329", "genus": "Neisseria", "species": "meningitidis",
     "provider": "pubmlst", "database": "pubmlst_neisseria_seqdef", "scheme_id": "88",
     "scheme_name": "N. meningitidis cgMLST v3", "revision": "last_updated 2026-09-12",
     "locus_count": 1329,
     "target_list_sha256": "b7330a96db12f2525d93dfd02519c1ce5342ed9dc6100a20660e71fb82dfd2f0",
     "has_profiles": True, "profile_field": "cgST",
     "threshold_scheme_key": None,
     "binding_basis": "No cutoff for N. meningitidis is curated in the threshold catalogue.",
     "notes": ()},
    {"key": "pubmlst:hinfluenzae-cgmlst-1037", "genus": "Haemophilus", "species": "influenzae",
     "provider": "pubmlst", "database": "pubmlst_hinfluenzae_seqdef", "scheme_id": "56",
     "scheme_name": "cgMLST v1", "revision": "last_updated 2026-09-10", "locus_count": 1037,
     "target_list_sha256": "977aaee3a66ffe6d491d4af017bd5c9ef31029d755b00ca33afdcdba1e28dc19",
     "has_profiles": True, "profile_field": "cgST",
     "threshold_scheme_key": None,
     "binding_basis": "No cutoff for H. influenzae is curated in the threshold catalogue.",
     "notes": ()},
    {"key": "pubmlst:bcc-cgmlst-1925", "genus": "Burkholderia", "species": "",
     "organism": "Burkholderia cepacia complex",
     "provider": "pubmlst", "database": "pubmlst_bcc_seqdef", "scheme_id": "2",
     "scheme_name": "cgMLST", "revision": "last_updated 2026-08-31", "locus_count": 1925,
     "target_list_sha256": "7153c3561f64ef5eafef127bfce592f46e8688917fed44634b198c7c0b868749",
     "has_profiles": True, "profile_field": "cgST",
     "threshold_scheme_key": None,
     "binding_basis": "No cutoff for the B. cepacia complex is curated in the threshold catalogue.",
     "notes": ("A complex-level scheme. Species within the complex are not separated by it.",)},
    {"key": "pubmlst:senterica-cgmlst-2750", "genus": "Salmonella", "species": "enterica",
     "provider": "pubmlst", "database": "pubmlst_salmonella_seqdef", "scheme_id": "3",
     "scheme_name": "SalmcgMLST v1.0", "revision": "last_updated 2026-04-06", "locus_count": 2750,
     "target_list_sha256": "b1649c4935c8ee19133240f389efe4faafd436faea7594cb57c66796cc4d66f5",
     "has_profiles": True, "profile_field": "cgST",
     "threshold_scheme_key": None,
     "binding_basis": ("No threshold is bound. The PulseNet criterion cites the 3,002-locus "
                       "EnteroBase-derived scheme, not this 2,750-target v1.0 scheme."),
     "notes": ("Unlike the 3,002-target v2 scheme this one carries cgST profiles, so a cgMLST "
               "sequence type can be assigned -- but no curated cutoff applies to it.",)},

    # --- Download-only: licence does not permit packing, directory is pre-created ---
    {"key": "cgmlst.org:lmonocytogenes-1701", "genus": "Listeria", "species": "monocytogenes",
     "provider": "cgmlst.org", "database": "", "scheme_id": "Lmonocytogenes",
     "scheme_name": "Listeria monocytogenes cgMLST", "revision": "catalogue read 2026-09-14",
     "locus_count": 1701,
     "target_list_sha256": "a48b9e9b77d7eeb3ca55709847e1cc132294a931bc931261c74c3540f86ab944",
     "has_profiles": False, "profile_field": "",
     "threshold_scheme_key": "cgmlst.org:lmonocytogenes-1701",
     "binding_basis": ("The scheme the ECDC multi-country evaluation (Van Walle et al. 2018) used, "
                       "identified by provider and its full 1,701-target set."),
     "notes": ("Ten alleles in this scheme carry a literal 'X' where a base is unknown, which is "
               "not an IUPAC nucleotide code. Those allele records are excluded on download and "
               "named in the snapshot manifest; a genome carrying one of them is reported as an "
               "unmatched sequence, never as that allele number.",)},
    {"key": "cgmlst.org:kpneumoniae-2358", "genus": "Klebsiella", "species": "pneumoniae",
     "provider": "cgmlst.org", "database": "", "scheme_id": "Kpneumoniae_complex",
     "scheme_name": "Klebsiella pneumoniae/variicola/quasipneumoniae cgMLST",
     "revision": "catalogue read 2026-09-14", "locus_count": 2358,
     "target_list_sha256": "1f6a89f2ef8b62fd40a8dfc69f0af257a5922da9b9f5ffe9c49fdd75145034c9",
     "has_profiles": False, "profile_field": "",
     "threshold_scheme_key": "cgmlst.org:kpneumoniae-2358",
     "binding_basis": "Provider and full 2,358-target set match the curated Glasgow 2025 entry.",
     "notes": ("A species-complex scheme: K. variicola and K. quasipneumoniae share it with "
               "K. pneumoniae, and the curated cutoff was reported for K. pneumoniae.",)},
    {"key": "cgmlst.org:efaecium-1423", "genus": "Enterococcus", "species": "faecium",
     "provider": "cgmlst.org", "database": "", "scheme_id": "Efaecium",
     "scheme_name": "Enterococcus faecium cgMLST", "revision": "catalogue read 2026-09-14",
     "locus_count": 1423,
     "target_list_sha256": "516ed3fa4c55cfc39ac86385f8957bb53fe94d602c214b4eca744f50dc043d3c",
     "has_profiles": False, "profile_field": "",
     "threshold_scheme_key": "cgmlst.org:efaecium-1423",
     "binding_basis": "Provider and full 1,423-target de Been scheme match the curated entries.",
     "notes": ("Three curated cutoffs cite this target set for different questions (7, 20 and 25). "
               "They are not interchangeable; read each entry's scope before choosing.",)},
    {"key": "cgmlst.org:efaecalis-1972", "genus": "Enterococcus", "species": "faecalis",
     "provider": "cgmlst.org", "database": "", "scheme_id": "Efaecalis",
     "scheme_name": "Enterococcus faecalis cgMLST", "revision": "catalogue read 2026-09-14",
     "locus_count": 1972,
     "target_list_sha256": "300c79fd7c6de8d97ecb00317253617c7d7d75613b51b839770cd49402ef2080",
     "has_profiles": False, "profile_field": "",
     "threshold_scheme_key": "cgmlst.org:efaecalis-1972",
     "binding_basis": "Provider and full 1,972-target set match the curated Glasgow 2025 entry.",
     "notes": ()},
    {"key": "cgmlst.org:abaumannii-2390", "genus": "Acinetobacter", "species": "baumannii",
     "provider": "cgmlst.org", "database": "", "scheme_id": "Abaumannii",
     "scheme_name": "Acinetobacter baumannii cgMLST", "revision": "catalogue read 2026-09-14",
     "locus_count": 2390,
     "target_list_sha256": "2e62dc87375d23c471a2d0cdcc4e5e6fee5464e39d710602ea07119299c83e1a",
     "has_profiles": False, "profile_field": "",
     "threshold_scheme_key": "cgmlst.org:abaumannii-2390",
     "binding_basis": "Provider and full 2,390-target set match the curated Glasgow 2025 entry.",
     "notes": ()},
    {"key": "cgmlst.org:saureus-1861", "genus": "Staphylococcus", "species": "aureus",
     "provider": "cgmlst.org", "database": "", "scheme_id": "Saureus",
     "scheme_name": "Staphylococcus aureus cgMLST", "revision": "catalogue read 2026-09-14",
     "locus_count": 1861,
     "target_list_sha256": "7958a5fe362488e9b30dc1b6dcb20b8488d3f824171c9f69bd33b750a334aa26",
     "has_profiles": False, "profile_field": "",
     "threshold_scheme_key": "cgmlst.org:saureus-1861",
     "binding_basis": "Provider and full 1,861-target set match the curated Glasgow 2025 entry.",
     "notes": ("Measured on 2026-09-14: this scheme's allele files occupy about 4.6 GB once "
               "installed, and the download took roughly 3.5 minutes on a fast connection. Check "
               "free disk space before starting.",)},
    {"key": "cgmlst.org:cdifficile-2147", "genus": "Clostridioides", "species": "difficile",
     "provider": "cgmlst.org", "database": "", "scheme_id": "Cdifficile",
     "scheme_name": "Clostridioides difficile cgMLST", "revision": "catalogue read 2026-09-14",
     "locus_count": 2147,
     "target_list_sha256": "e2c3dbdd72e3e60ac87ad0c8728d357143052a28658164023d1d395a0caa0ee7",
     "has_profiles": False, "profile_field": "",
     "threshold_scheme_key": None,
     "binding_basis": ("No threshold is bound. The curated Bletz 2018 cutoff of 6 cites the "
                       "2,270-target v1 scheme, which the provider now lists as deprecated. The "
                       "current scheme has 2,147 targets, so the published cutoff does not carry "
                       "over and no number is offered."),
     "notes": ("If you must use the published 6-allele criterion, install the deprecated 2,270-"
               "target v1 scheme and record which version you ran.",)},
    {"key": "cgmlst.org:smarcescens-2692", "genus": "Serratia", "species": "marcescens",
     "provider": "cgmlst.org", "database": "", "scheme_id": "Smarcescens",
     "scheme_name": "Serratia marcescens cgMLST", "revision": "catalogue read 2026-09-14",
     "locus_count": 2692,
     "target_list_sha256": "0aee3a6a8a4fed9138565a797defa34f0e4db23c8ee7dbc3d439971d1a9a56e1",
     "has_profiles": False, "profile_field": "",
     "threshold_scheme_key": "cgmlst.org:smarcescens-2692",
     "binding_basis": "Provider and full 2,692-target Db11-derived set match the Kampmeier 2022 entry.",
     "notes": ()},
    {"key": "cgmlst.org:cfreundii-3250", "genus": "Citrobacter", "species": "freundii",
     "provider": "cgmlst.org", "database": "", "scheme_id": "Cfreundii",
     "scheme_name": "Citrobacter freundii cgMLST", "revision": "catalogue read 2026-09-14",
     "locus_count": 3250,
     "target_list_sha256": "986350d4057ba7d8104586edf46d89b16cea72963e34ee422696dba108c7568a",
     "has_profiles": False, "profile_field": "",
     "threshold_scheme_key": "kieninger2025:cfreundii-3250",
     "binding_basis": ("Provider and full 3,250-target species-specific set match the Kieninger "
                       "2025 entry, whose scheme key names the publication rather than a slug."),
     "notes": ("Not the paper's separate combined-species Citrobacter scheme; that is the 2,307-"
               "target Cfreundii_complex scheme and no cutoff is curated for it.",)},
    {"key": "cgmlst.org:paeruginosa-3867", "genus": "Pseudomonas", "species": "aeruginosa",
     "provider": "cgmlst.org", "database": "", "scheme_id": "Paeruginosa",
     "scheme_name": "Pseudomonas aeruginosa cgMLST", "revision": "catalogue read 2026-09-14",
     "locus_count": 3867,
     "target_list_sha256": "5cbf019a4ac19126f96a5b3cdc9123a2d74230a032e2ca279f819a6a90aac1b2",
     "has_profiles": False, "profile_field": "",
     "threshold_scheme_key": "cgmlst.org:paeruginosa-3867",
     "binding_basis": ("Provider and full 3,867-target set match the Tonnies 2021 entry, which "
                       "deliberately curates NO numeric cutoff."),
     "notes": ("The bound entry supplies a citation and scope only. No transferable number exists "
               "for this organism in the catalogue.",)},
    {"key": "pasteur:lmonocytogenes-1748", "genus": "Listeria", "species": "monocytogenes",
     "provider": "pasteur", "database": "pubmlst_listeria_seqdef", "scheme_id": "3",
     "scheme_name": "cgMLST1748", "revision": "last_updated 2023-09-22", "locus_count": 1748,
     "has_profiles": True, "profile_field": "profile_id",
     "target_list_sha256": "5ee68c96379be84827987bc0977d630a329b02e7e5b330f2be748b19d9ad6288",
     "threshold_scheme_key": "pasteur:lmonocytogenes-1748",
     "binding_basis": ("The Institut Pasteur BIGSdb-Lm 1,748-locus scheme in which Moura et al. "
                       "2016 defined the cgMLST type."),
     "notes": ("Measured on 2026-09-14: downloading this scheme from the BIGSdb REST API means "
               "1,748 separate requests and roughly 2.5 GB over about 17 minutes, because BIGSdb "
               "serves every allele ever submitted. Expect a long download.",
               "Anonymous access is limited to submissions up to 2024-12-31, so the snapshot is "
               "not the current allele set. Institut Pasteur forbids redistribution of it.",
               "The provider also publishes cgMLST1748_v2, flagged experimental, on the same 1,748 "
               "targets with a separate profile numbering. They are not interchangeable.")},
)


class SchemeCatalogError(ValueError):
    """A catalogue entry or library slot could not be resolved without a guess."""


def _organism(entry: dict) -> str:
    if entry.get("organism"):
        return str(entry["organism"])
    return " ".join(part for part in (entry.get("genus"), entry.get("species")) if part)


def catalog_entries() -> list[dict]:
    """Every pinned scheme, with its provider licence verdict resolved onto the row."""
    rows = []
    for raw in _SCHEMES:
        entry = deepcopy(dict(raw))
        provider = PROVIDERS[entry["provider"]]
        entry["organism"] = _organism(entry)
        entry["provider_name"] = provider["name"]
        entry["terms_url"] = provider["terms_url"]
        entry["redistribution"] = provider["redistribution"]
        entry["bundled"] = bool(provider["may_bundle"])
        entry["requires_terms_acknowledgement"] = bool(provider["requires_acknowledgement"])
        entry["licence_quote"] = provider["quote"]
        entry["licence_restriction"] = provider["restriction"]
        entry["licence_reviewed_on"] = provider["reviewed_on"]
        entry["kind"] = "cgmlst"
        entry["slot"] = slot_name(entry)
        entry["source_url"] = source_url(entry)
        entry["notes"] = list(entry.get("notes", ()))
        rows.append(entry)
    return rows


def entry_for(key) -> dict:
    """One pinned scheme by key; a dict entry is accepted and returned resolved."""
    wanted = key.get("key") if isinstance(key, dict) else key
    for entry in catalog_entries():
        if entry["key"] == wanted:
            return entry
    raise SchemeCatalogError(f"Unknown cgMLST scheme key: {wanted!r}")


def bundled_entries() -> list[dict]:
    """Schemes whose provider grants redistribution, so a build may stage them."""
    return [entry for entry in catalog_entries() if entry["bundled"]]


def download_only_entries() -> list[dict]:
    """Schemes that may not be packed; the folder is created and the user downloads."""
    return [entry for entry in catalog_entries() if not entry["bundled"]]


def source_url(entry: dict) -> str:
    provider = PROVIDERS[entry["provider"]]
    if entry["provider"] == "cgmlst.org":
        return f"{provider['api_root']}/schema/{entry['scheme_id']}/"
    return f"{provider['api_root']}/db/{entry['database']}/schemes/{entry['scheme_id']}"


def slot_name(entry: dict) -> str:
    """The library folder for one scheme: organism first, then provider and target count."""
    organism = re.sub(r"[^A-Za-z0-9]+", "_", _organism(entry) or "Unknown_organism").strip("_")
    provider = re.sub(r"[^A-Za-z0-9]+", "_", entry["provider"]).strip("_")
    return f"{organism}__{provider}_{entry['locus_count']}"


def catalog_digest() -> str:
    """A fingerprint of the pinned catalogue, for staging manifests and tests."""
    pinned = [{key: entry[key] for key in
               ("key", "provider", "database", "scheme_id", "locus_count",
                "target_list_sha256", "threshold_scheme_key")}
              for entry in catalog_entries()]
    return hashlib.sha256(json.dumps(pinned, sort_keys=True).encode()).hexdigest()


# ---------------------------------------------------------------- threshold binding

def _bound_guidance(scheme_key: str) -> list[dict]:
    """Publication rows whose own scheme_key is this scheme, newest reviewed source first.

    The ordering deliberately matches threshold_guidance.suggested_threshold, so the
    per-organism suggestion and this per-scheme lookup can never show a person two
    different numbers for the same comparison.
    """
    rows = []
    for entry in _threshold_guidance_entries():
        if entry["scheme_key"] == scheme_key and entry["method"] == "cgmlst":
            row = deepcopy(entry)
            row["source"] = deepcopy(SOURCES[entry["source_id"]])
            row["published"] = row["source"]["published"]
            rows.append(row)
    rows.sort(key=lambda row: (row["published"], row["id"]), reverse=True)
    return rows


def threshold_for(key, *, locus_count=None, method="cgmlst") -> dict:
    """The suggested cutoff and citation for an installed scheme, or why there is none.

    A number is offerable only when the pinned scheme carries a threshold_scheme_key,
    the publication entry pins the same key, and the installed target count equals
    BOTH the catalogue's pinned count and the published count. `locus_count` is the
    count actually found on disk; pass it whenever a scheme is installed, so a
    truncated or partial snapshot cannot inherit a cutoff derived on the full set.
    """
    if method != "cgmlst":
        raise SchemeCatalogError(
            "This catalogue binds cgMLST target-set cutoffs only. A seven-locus MLST distance is a "
            "different quantity and must not be compared against a cgMLST threshold.")
    entry = entry_for(key)
    result = {"key": entry["key"], "organism": entry["organism"], "method": "cgmlst",
              "scheme_key": entry["threshold_scheme_key"],
              "catalogue_locus_count": entry["locus_count"], "installed_locus_count": locus_count,
              "binding_basis": entry["binding_basis"], "catalog_version": CATALOG_VERSION,
              "interpretation": INTERPRETATION, "entries": [], "threshold": None,
              "suggestion": None, "alternatives": [], "status": "", "reason": ""}
    if entry["locus_count"] <= CGMLST_TARGET_FLOOR:
        result["status"] = "not_a_cgmlst_scheme"
        result["reason"] = (f"{entry['organism']} {entry['scheme_name']} declares only "
                            f"{entry['locus_count']} targets, at or below the {CGMLST_TARGET_FLOOR}-"
                            "locus floor that separates classical MLST from gene-by-gene typing.")
        return result
    if not entry["threshold_scheme_key"]:
        result["status"] = "no_bound_cutoff"
        result["reason"] = entry["binding_basis"]
        return result
    rows = _bound_guidance(entry["threshold_scheme_key"])
    if not rows:
        result["status"] = "no_curated_cutoff"
        result["reason"] = (f"No reviewed cutoff in the publication catalogue names "
                            f"{entry['threshold_scheme_key']}. This is an evidence gap, not proof "
                            "that no publication exists.")
        return result
    result["entries"] = rows
    if locus_count is not None and locus_count != entry["locus_count"]:
        result["status"] = "target_count_mismatch"
        result["reason"] = (f"The installed scheme reports {locus_count} targets but "
                            f"{entry['organism']} {entry['scheme_name']} is pinned at "
                            f"{entry['locus_count']}. A cutoff derived on the full published target "
                            "set is not offered for a different one; no scaling is applied.")
        return result
    numeric = [row for row in rows if row["published_threshold"] is not None
               and (row["locus_count"] is None or row["locus_count"] == entry["locus_count"])]
    if not numeric:
        result["status"] = "citation_only"
        result["reason"] = ("The bound publication supplies scope and citation but no transferable "
                            "number for this target set.")
        return result
    result["status"] = "threshold_offerable"
    result["suggestion"] = numeric[0]
    result["threshold"] = numeric[0]["published_threshold"]
    result["alternatives"] = [row for row in rows if row is not numeric[0]]
    result["reason"] = (f"{len(numeric)} reviewed cutoff(s) are bound to this exact scheme and its "
                        f"full {entry['locus_count']}-target set. The most recently reviewed source "
                        "is suggested and the others are returned as alternatives, because picking "
                        "whichever number is most convenient is the mistake this catalogue exists "
                        "to prevent. A published cutoff is evidence to review, never a setting.")
    return result


def threshold_citations(key) -> list[dict]:
    """Citation rows for a scheme, in the shape a PDF report footnote needs."""
    rows = []
    for row in threshold_for(key)["entries"]:
        source = row["source"]
        rows.append({"entry_id": row["id"], "organism": row["organism"],
                     "published_threshold": row["published_threshold"],
                     "operator": row["operator"], "unit": row["unit"],
                     "scheme_key": row["scheme_key"], "locus_count": row["locus_count"],
                     "citation": source["citation"], "doi": source["doi"], "url": source["url"],
                     "published": source["published"], "locator": source["locator"],
                     "scope": row["scope"], "missing_policy": row["missing_policy"],
                     "limitations": row["limitations"], "reviewed_on": row["reviewed_on"]})
    return rows


# ---------------------------------------------------------------- library layout

SLOT_README = """{organism} -- {scheme_name}

Put the {locus_count}-target {provider_name} cgMLST scheme in THIS folder.

{state_line}

Provider : {provider_name}
Scheme   : {source_url}
Targets  : {locus_count}
Terms    : {terms_url}

{restriction}

This folder is created by WMLSTudio so that you never have to make one. It is
safe to delete: the app recreates it, and deleting it removes no scheme you have
already installed. Nothing in this folder is sequence data until you install a
scheme into it.

A cgMLST distance from this scheme is a number of differing targets out of the
targets called in BOTH isolates. It is not a seven-locus MLST distance, it is not
a SNP count, and it is never proof of transmission.
"""

_BUNDLED_STATE = ("This scheme's provider permits redistribution, so the release may already "
                  "carry it. If the folder is empty, use Reference data > Install cgMLST scheme.")
_DOWNLOAD_STATE = ("This scheme CANNOT be packed into WMLSTudio: its provider does not grant "
                   "redistribution. Use Reference data > Download cgMLST scheme. You will be "
                   "shown the provider's terms and asked to confirm that your use is permitted "
                   "before anything is downloaded.")


def library_root(root) -> Path:
    """<data root>/cgmlst -- the one place installed cgMLST schemes live."""
    return Path(root).expanduser() / LIBRARY_DIRNAME


def slot_payload(entry: dict) -> dict:
    """The machine-readable description written into an empty slot."""
    return {"format_version": SLOT_FORMAT_VERSION, "catalog_version": CATALOG_VERSION,
            "key": entry["key"], "organism": entry["organism"], "genus": entry["genus"],
            "species": entry["species"], "kind": "cgmlst",
            "provider": entry["provider"], "provider_name": entry["provider_name"],
            "database": entry["database"], "scheme_id": entry["scheme_id"],
            "scheme_name": entry["scheme_name"], "revision": entry["revision"],
            "locus_count": entry["locus_count"],
            "target_list_sha256": entry["target_list_sha256"],
            "has_profiles": entry["has_profiles"], "profile_field": entry["profile_field"],
            "source_url": entry["source_url"], "terms_url": entry["terms_url"],
            "redistribution": entry["redistribution"], "bundled": entry["bundled"],
            "requires_terms_acknowledgement": entry["requires_terms_acknowledgement"],
            "licence_quote": entry["licence_quote"],
            "licence_restriction": entry["licence_restriction"],
            "licence_reviewed_on": entry["licence_reviewed_on"],
            "threshold_scheme_key": entry["threshold_scheme_key"],
            "binding_basis": entry["binding_basis"], "notes": list(entry["notes"]),
            "interpretation": INTERPRETATION}


def prepare_library(root, *, keys=None, targets=None) -> dict:
    """Create the cgMLST library layout so a first run never asks for a folder.

    Every catalogued scheme gets a named, empty folder carrying a README and a
    scheme_slot.json. Existing folders and any scheme already installed in them are
    left untouched: only the two description files are refreshed, and only when
    their content changed. `targets` optionally maps a key to its target-name list,
    which is written as targets.txt beside the slot.
    """
    base = library_root(root)
    base.mkdir(parents=True, exist_ok=True)
    wanted = set(keys) if keys is not None else None
    report = {"root": str(base), "created": [], "existing": [], "refreshed": [],
              "catalog_version": CATALOG_VERSION, "catalog_digest": catalog_digest(),
              "generated_utc": datetime.now(UTC).isoformat()}
    for entry in catalog_entries():
        if wanted is not None and entry["key"] not in wanted:
            continue
        folder = base / entry["slot"]
        (report["existing"] if folder.is_dir() else report["created"]).append(entry["slot"])
        folder.mkdir(parents=True, exist_ok=True)
        readme = SLOT_README.format(
            organism=entry["organism"], scheme_name=entry["scheme_name"],
            locus_count=entry["locus_count"], provider_name=entry["provider_name"],
            source_url=entry["source_url"], terms_url=entry["terms_url"],
            restriction=entry["licence_restriction"],
            state_line=_BUNDLED_STATE if entry["bundled"] else _DOWNLOAD_STATE)
        payload = json.dumps(slot_payload(entry), indent=2, sort_keys=True) + "\n"
        for name, text in ((README_FILENAME, readme), (SLOT_FILENAME, payload)):
            path = folder / name
            if not path.is_file() or path.read_text(encoding="utf-8") != text:
                path.write_text(text, encoding="utf-8")
                report["refreshed"].append(f"{entry['slot']}/{name}")
        names = (targets or {}).get(entry["key"])
        if names:
            text = "\n".join(names) + "\n"
            path = folder / TARGETS_FILENAME
            if not path.is_file() or path.read_text(encoding="utf-8") != text:
                path.write_text(text, encoding="utf-8")
                report["refreshed"].append(f"{entry['slot']}/{TARGETS_FILENAME}")
    (base / README_FILENAME).write_text(LIBRARY_README, encoding="utf-8")
    return report


LIBRARY_README = """WMLSTudio cgMLST scheme library

One folder per catalogued cgMLST scheme. Each folder explains, in its own
README.txt, which scheme belongs there, who publishes it, under what terms, and
how to install it. An empty folder means that scheme is not installed yet -- it
does not mean anything is broken.

These are cgMLST schemes: hundreds to thousands of targets. Classical seven-locus
MLST schemes live in the separate schemes folder. Distances from the two are
different quantities and are never mixed, compared or thresholded together.

Deleting a folder here deletes only its description files unless you have
installed a scheme into it. WMLSTudio recreates the empty layout on the next run.
"""


def _installed_locus_count(folder: Path) -> int:
    suffixes = {".tfa", ".fasta", ".fa", ".fna"}
    count = 0
    for path in folder.iterdir():
        if not path.is_file():
            continue
        name = path.with_suffix("") if path.suffix.casefold() in {".gz", ".bz2"} else path
        if name.suffix.casefold() in suffixes:
            count += 1
    return count


def installed_scheme(root, key, *, extra_paths=()) -> dict | None:
    """The scheme installed for one catalogue key, or None. Counts files, parses none.

    A folder is accepted as this scheme only when its recorded provider and scheme
    identifier match the pin. A folder with the right number of allele files but a
    different or absent identity is reported with `identity_matched` False, so a
    caller can show it without ever binding a threshold to it.
    """
    entry = entry_for(key)
    candidates = [library_root(root) / entry["slot"]]
    candidates.extend(Path(path) for path in extra_paths)
    for folder in candidates:
        if not folder.is_dir():
            continue
        count = _installed_locus_count(folder)
        if not count:
            continue
        metadata = {}
        for name in ("scheme.json", "reference_manifest.json"):
            path = folder / name
            if path.is_file() and path.stat().st_size <= 4 * 1024 * 1024:
                try:
                    value = json.loads(path.read_text(encoding="utf-8"))
                except (OSError, UnicodeDecodeError, ValueError):
                    continue
                if isinstance(value, dict):
                    metadata[name] = value
        matched = _identity_matches(entry, metadata)
        return {"key": entry["key"], "path": str(folder), "locus_count": count,
                "expected_locus_count": entry["locus_count"],
                "count_matched": count == entry["locus_count"],
                "identity_matched": matched,
                "organism": entry["organism"], "scheme_name": entry["scheme_name"],
                "provider_name": entry["provider_name"]}
    return None


def _identity_matches(entry: dict, metadata: dict) -> bool:
    scheme = metadata.get("scheme.json") or {}
    manifest = metadata.get("reference_manifest.json") or {}
    recorded = manifest.get("entry") if isinstance(manifest.get("entry"), dict) else {}
    api = str(scheme.get("API") or recorded.get("url") or "")
    if entry["provider"] == "cgmlst.org":
        return f"/schema/{entry['scheme_id']}/" in api
    if entry["database"] and entry["scheme_id"]:
        return api.endswith(f"/db/{entry['database']}/schemes/{entry['scheme_id']}")
    return False


def library_status(root, *, extra_paths=()) -> list[dict]:
    """One row per catalogued scheme: where it belongs, whether it is there, what is
    offerable. This is the model a "cgMLST schemes" page renders directly."""
    rows = []
    for entry in catalog_entries():
        installed = installed_scheme(root, entry["key"], extra_paths=extra_paths)
        count = installed["locus_count"] if installed else None
        guidance = threshold_for(entry["key"], locus_count=count)
        rows.append({**entry, "folder": str(library_root(root) / entry["slot"]),
                     "installed": installed, "ready": bool(installed and installed["count_matched"]
                                                           and installed["identity_matched"]),
                     "threshold": guidance})
    return rows


def identify_installed(path) -> dict | None:
    """Which catalogued scheme a scheme folder is, read from what it records itself.

    Returns None when the folder records no identity this catalogue recognises.
    That is the honest answer: an unrecognised scheme is still perfectly usable for
    comparison, it simply has no bound cutoff.
    """
    folder = Path(path).expanduser()
    if not folder.is_dir():
        return None
    metadata = {}
    for name in ("scheme.json", "reference_manifest.json"):
        candidate = folder / name
        if candidate.is_file() and candidate.stat().st_size <= 4 * 1024 * 1024:
            try:
                value = json.loads(candidate.read_text(encoding="utf-8"))
            except (OSError, UnicodeDecodeError, ValueError):
                continue
            if isinstance(value, dict):
                metadata[name] = value
    for entry in catalog_entries():
        if _identity_matches(entry, metadata):
            return {**entry, "path": str(folder), "locus_count_on_disk": _installed_locus_count(folder)}
    return None


def threshold_for_installed(path) -> dict:
    """The cutoff lookup for a scheme folder on disk, counted rather than assumed.

    This is the entry point for a report or a comparison view that holds a
    scheme_path. The target count is read from the folder, so a snapshot that is
    missing targets cannot inherit a cutoff published for the complete set.
    """
    identified = identify_installed(path)
    if identified is None:
        return {"key": None, "path": str(path), "status": "scheme_not_catalogued",
                "threshold": None, "entries": [], "catalog_version": CATALOG_VERSION,
                "interpretation": INTERPRETATION,
                "reason": ("This scheme records no identity the cgMLST catalogue recognises, so no "
                           "published cutoff is bound to it. Distances can still be computed and "
                           "compared within this scheme.")}
    guidance = threshold_for(identified["key"], locus_count=identified["locus_count_on_disk"])
    guidance["path"] = str(path)
    return guidance


def download_plan(key) -> dict:
    """Everything a download button must show BEFORE any byte is transferred."""
    entry = entry_for(key)
    provider = PROVIDERS[entry["provider"]]
    if entry["provider"] == "cgmlst.org":
        method = ("One archive of all allele FASTA files, then a locus table to verify it. "
                  "Typically tens of megabytes compressed; several gigabytes once installed.")
    elif entry["provider"] in {"pubmlst", "pasteur"}:
        method = (f"{entry['locus_count']} separate requests, one per target, because the service "
                  "offers no whole-scheme archive. Expect a long download and gigabytes on disk; "
                  "it can be cancelled and resumed.")
    else:
        method = "No supported download route is catalogued for this provider."
    return {"key": entry["key"], "organism": entry["organism"],
            "scheme_name": entry["scheme_name"], "locus_count": entry["locus_count"],
            "provider": entry["provider"], "provider_name": entry["provider_name"],
            "source_url": entry["source_url"], "terms_url": provider["terms_url"],
            "terms_notice": provider["restriction"], "licence_quote": provider["quote"],
            "requires_acknowledgement": bool(provider["requires_acknowledgement"]),
            "may_be_bundled": bool(provider["may_bundle"]), "method": method,
            "destination_slot": entry["slot"], "notes": list(entry["notes"])}
