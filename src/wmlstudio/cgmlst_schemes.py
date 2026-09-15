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

TARGET SETS. A cgMLST scheme is a CORE set. PubMLST publishes an accessory set
beside the core one for Bacillus anthracis, and an accessory and a pan-genome set
beside the two gonococcal core schemes; those are pinned here as separate rows so
"core" and "core plus accessory" are an explicit choice rather than one silently
standing in for the other. A core-set distance and a core-plus-accessory distance
are different quantities for the same reason MLST and cgMLST are: they count
differences over different numbers of different targets. They are separate rows,
separate keys, separate library folders and separate SHA-256 pins, and a cutoff
published on a core set is never offered for a run on any other set.
"""

from __future__ import annotations

import hashlib
import json
import re
import shutil
import urllib.parse
from copy import deepcopy
from datetime import UTC, datetime
from pathlib import Path

from .threshold_guidance import SOURCES
from .threshold_guidance import catalog_entries as _threshold_guidance_entries

CATALOG_VERSION = "2026-09-15.1"
LICENCE_REVIEWED_ON = "2026-09-14"
LIBRARY_DIRNAME = "cgmlst"
SCHEME_DIRNAME = "schemes"
SLOT_FILENAME = "scheme_slot.json"
TARGETS_FILENAME = "targets.txt"
README_FILENAME = "README.txt"
MIGRATIONS_FILENAME = "migrations.json"
SLOT_FORMAT_VERSION = 1
# The same floor project.classify_typing and reference_index use to separate a
# classical seven-locus scheme from a gene-by-gene one. Repeated, not imported,
# to keep this module free of the storage/typing import chain.
CGMLST_TARGET_FLOOR = 30
# A cgMLST scheme is a CORE target set: every target is expected in every isolate
# of the organism. Some providers publish an accessory or whole-genome set beside
# it. The two are different quantities, never share a cutoff and never share an
# axis, so which one a scheme is is recorded on the row rather than inferred.
TARGET_SET_CORE = "core"
TARGET_SET_ACCESSORY = "accessory"
# A pan-genome (wgMLST) set is core and accessory targets in ONE set, which is not
# the same thing as an accessory set: 251 accessory targets and 1,907 pan-genome
# targets are themselves two different quantities. The rest of the application
# stores and renders the two-valued bucket in TARGET_SETS -- a set either IS the
# core set or is not -- so the finer name is carried beside it, never instead of it.
TARGET_SET_WHOLE_GENOME = "whole_genome"
TARGET_SETS = (TARGET_SET_CORE, TARGET_SET_ACCESSORY)
TARGET_SET_DETAILS = (TARGET_SET_CORE, TARGET_SET_ACCESSORY, TARGET_SET_WHOLE_GENOME)
_TARGET_SET_LABELS = {TARGET_SET_CORE: "Core", TARGET_SET_ACCESSORY: "Accessory",
                      TARGET_SET_WHOLE_GENOME: "Whole genome"}
_TARGET_SET_PHRASES = {
    TARGET_SET_CORE: "cgMLST core target set",
    TARGET_SET_ACCESSORY: "accessory target set",
    TARGET_SET_WHOLE_GENOME: "whole-genome (pan-genome) target set",
}
_ALLELE_SUFFIXES = {".tfa", ".fasta", ".fa", ".fna"}
# Provider spellings differ between the pinned catalogue ("pasteur") and the
# download clients ("BIGSdb-Pasteur"); both name one provider and must resolve to
# one library folder.
_PROVIDER_TOKENS = {"bigsdb_pasteur": "pasteur", "bigsdb": "pasteur", "ridom": "cgmlst_org",
                    "cgmlst_org_nomenclature_server": "cgmlst_org"}
_PROVIDER_HOSTS = {"rest.pubmlst.org": "pubmlst", "pubmlst.org": "pubmlst",
                   "bigsdb.pasteur.fr": "pasteur", "www.cgmlst.org": "cgmlst_org",
                   "cgmlst.org": "cgmlst_org"}
_KNOWN_PROVIDERS = frozenset({"pubmlst", "pasteur", "cgmlst_org", "enterobase", "chewie_ns"})

INTERPRETATION = (
    "An installed cgMLST scheme is a target set and an allele nomenclature. Two profiles are "
    "comparable only when both were called against the same installed scheme; a distance never "
    "crosses schemes, providers or revisions, and never shares a scale with a seven-locus MLST "
    "distance. A core target set and a core-plus-accessory set are two such schemes: their "
    "distances never share a scale, an axis or a threshold either."
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

    # --- Core and accessory sets published side by side, so the choice is explicit ---
    # These are the only two organisms for which any provider this catalogue supports
    # publishes a non-core target set. Every PubMLST and every Institut Pasteur
    # seqdef database was listed on 2026-09-15 and searched for a scheme describing
    # itself as accessory, pan-genome or whole-genome; cgMLST.org publishes core
    # schemes only. C. chauvoei has such a pair too and is left out as veterinary.
    {"key": "pubmlst:banthracis-cgmlst-3803", "genus": "Bacillus", "species": "anthracis",
     "provider": "pubmlst", "database": "pubmlst_bcereus_seqdef", "scheme_id": "2",
     "scheme_name": "B. anthracis cgMLST", "revision": "last_updated 2026-09-07",
     "locus_count": 3803,
     "target_list_sha256": "e4a911b675778ede04da21b3aad4adde951d4f6e1ccb1df5cdb55440aa6da5f2",
     "has_profiles": True, "profile_field": "cgST",
     "threshold_scheme_key": None,
     "binding_basis": ("No threshold is bound. Abdel-Glil et al. 2021 propose five differing "
                       "alleles on this 3,803-target core set to trace epidemiologically linked "
                       "strains, but that number is not curated in this application's publication "
                       "catalogue, so no number is offered. Read the paper before using it."),
     "notes": ("Abdel-Glil et al. 2021 (J Clin Microbiol 59:e02889-20) defined this scheme on 57 "
               "B. anthracis genomes spanning the phylogeny and evaluated it on 584 genomes from "
               "50 countries.",
               "The same PubMLST database also hosts a 1,568-target B. cereus cgMLST scheme. It "
               "shares 1,225 targets with this one, verified on 2026-09-15, and is a different "
               "scheme: the two are not comparable and are not two versions of one thing.")},
    {"key": "pubmlst:banthracis-accessory-1263", "genus": "Bacillus", "species": "anthracis",
     "provider": "pubmlst", "database": "pubmlst_bcereus_seqdef", "scheme_id": "3",
     "scheme_name": "B. anthracis accessory genes", "target_set": TARGET_SET_ACCESSORY,
     "revision": "read 2026-09-15; PubMLST lists no last_updated date for this definition-only "
     "scheme", "locus_count": 1263,
     "target_list_sha256": "9001693fde7ce0eb6095349244877fe3fe2ea6a00ad27eea237ad39a941ccc9a",
     "has_profiles": False, "profile_field": "",
     "threshold_scheme_key": None,
     "binding_basis": ("No threshold is bound and none can be. The only published B. anthracis "
                       "cutoff -- five alleles, Abdel-Glil et al. 2021 -- was derived on the "
                       "3,803-target CORE set. A distance over these 1,263 accessory targets is a "
                       "different quantity: the core cutoff is not offered for it, and the two "
                       "distances are never added together into one number."),
     "notes": ("These are the 1,263 accessory targets Abdel-Glil et al. 2021 host at PubMLST "
               "beside the core scheme; their wgMLST is the two sets run together, 3,803 + 1,263 "
               "= 5,066 targets. Verified on 2026-09-15: the two sets share no target name.",
               "An accessory target is absent from some B. anthracis isolates by design, so an "
               "uncalled target here is biology and not a failed call. It is reported as not "
               "assayed, never as a difference.",
               "PubMLST assigns no cgST for this set, so no sequence type is produced; allelic "
               "distances within this set are still computed.")},
    {"key": "pubmlst:ngonorrhoeae-cgmlst-1430", "genus": "Neisseria", "species": "gonorrhoeae",
     "provider": "pubmlst", "database": "pubmlst_neisseria_seqdef", "scheme_id": "89",
     "scheme_name": "N. gonorrhoeae cgMLST v2", "revision": "last_updated 2026-09-15",
     "locus_count": 1430,
     "target_list_sha256": "9f50a56e0c826d398b98ed9cba20225a0563cf80499a0fb9aa9e8e52da340582",
     "has_profiles": True, "profile_field": "cgST",
     "threshold_scheme_key": None,
     "binding_basis": "No cutoff for N. gonorrhoeae is curated in the threshold catalogue.",
     "notes": ("Unitt et al. 2025 (eLife 14) describe this refined scheme and the LIN code "
               "nomenclature PubMLST publishes on it. LIN code bin thresholds are a nomenclature "
               "for naming lineages, not an outbreak cutoff, and none is applied here.",
               "Verified on 2026-09-15: 12 of these 1,430 targets are in neither the 1,649-target "
               "cgMLST v1.0 set nor the 1,907-target pgMLST v1.0 set, so this scheme is NOT the "
               "core half of that pan-genome set.")},
    {"key": "pubmlst:ngonorrhoeae-cgmlst-1649", "genus": "Neisseria", "species": "gonorrhoeae",
     "provider": "pubmlst", "database": "pubmlst_neisseria_seqdef", "scheme_id": "62",
     "scheme_name": "N. gonorrhoeae cgMLST v1.0", "revision": "last_updated 2026-09-15",
     "locus_count": 1649,
     "target_list_sha256": "e99d48e926952cb8c6a00ad9ab60348c18e1320ad2d414d27140fb2d7799dbd2",
     "has_profiles": True, "profile_field": "cgST",
     "threshold_scheme_key": None,
     "binding_basis": "No cutoff for N. gonorrhoeae is curated in the threshold catalogue.",
     "notes": ("Harrison et al. 2020 (J Infect Dis 222:1816-1825) defined this gonococcal core "
               "genome. It is the core set the accessory and pan-genome v1.0 schemes are built "
               "around: verified on 2026-09-15, these 1,649 targets and the 251 accessory targets "
               "are disjoint and together account for 1,900 of the 1,907 pgMLST targets.",
               "PubMLST curates cgMLST v2 as the typing scheme now. v1.0 and v2 are different "
               "target sets, not two revisions of one, and their distances never share a scale.")},
    {"key": "pubmlst:ngonorrhoeae-accessory-251", "genus": "Neisseria", "species": "gonorrhoeae",
     "provider": "pubmlst", "database": "pubmlst_neisseria_seqdef", "scheme_id": "80",
     "scheme_name": "N. gonorrhoeae agMLST v1.0", "target_set": TARGET_SET_ACCESSORY,
     "revision": "read 2026-09-15; PubMLST lists no last_updated date and flags this scheme as in "
     "development", "locus_count": 251,
     "target_list_sha256": "d6ee3dc3344878a4953def93ae358b56f68c3246eefe97d3a86018a8bb31cb6a",
     "has_profiles": False, "profile_field": "",
     "threshold_scheme_key": None,
     "binding_basis": ("No threshold is bound and none can be. No cutoff for N. gonorrhoeae is "
                       "curated at all, and a cutoff published on a core set would not carry over "
                       "to a 251-target accessory set in any case."),
     "notes": ("PubMLST describes this scheme as the gonococcal accessory genome and labels it "
               "'in development' and 'unpublished'. It downloads and it is usable within itself, "
               "but the target set may be redefined without notice; the pinned SHA-256 is what "
               "detects that instead of silently adopting it.",
               "Verified on 2026-09-15: none of these 251 targets appears in either gonococcal "
               "core scheme. An accessory target is absent from some isolates by design, so an "
               "uncalled target is biology and not a failed call.")},
    {"key": "pubmlst:ngonorrhoeae-pgmlst-1907", "genus": "Neisseria", "species": "gonorrhoeae",
     "provider": "pubmlst", "database": "pubmlst_neisseria_seqdef", "scheme_id": "81",
     "scheme_name": "N. gonorrhoeae pgMLST v1.0", "target_set": TARGET_SET_WHOLE_GENOME,
     "revision": "read 2026-09-15; PubMLST lists no last_updated date and flags this scheme as in "
     "development", "locus_count": 1907,
     "target_list_sha256": "22f8f05e2c648b530c7f12ce96128ba6c6b67fa867ccaa9cf2970ae6f2d5e049",
     "has_profiles": False, "profile_field": "",
     "threshold_scheme_key": None,
     "binding_basis": ("No threshold is bound. A pan-genome distance over 1,907 targets is not a "
                       "cgMLST distance; no cutoff is published for it, and none from a core "
                       "scheme is offered in its place."),
     "notes": ("This is the 'cgMLST plus accessory genes' set for N. gonorrhoeae. Run it INSTEAD "
               "of a core scheme, never beside one as something added to a core distance.",
               "Verified on 2026-09-15: it contains all 1,649 cgMLST v1.0 targets and all 251 "
               "agMLST v1.0 targets, plus 7 further targets in neither (NEIS1391, NEIS3185, "
               "NEIS3191, NEIS3199, NEIS3201, NEIS3235, NEIS3236), so it is not the arithmetic "
               "union of those two schemes.",
               "PubMLST labels this scheme 'in development' and 'unpublished'; the target set may "
               "be redefined without notice, and the pinned SHA-256 is what detects that.")},

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


def _version_of(entry: dict) -> str:
    """A short revision date pulled out of the pinned revision sentence, or ''."""
    match = re.search(r"\d{4}-\d{2}-\d{2}", str(entry.get("revision") or ""))
    return match.group(0) if match else ""


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
        # A row that says nothing is a core set, because that is what a cgMLST
        # scheme is; a row pinning a provider's accessory or pan-genome set names
        # it, so the two appear as an explicit choice instead of one silently
        # standing in for the other. Two fields carry one fact on purpose:
        # target_set is the two-valued bucket the rest of the application stores
        # and renders -- a set either IS the core set or is not -- and
        # target_set_detail says which kind of not-core set it is, because an
        # accessory distance and a pan-genome distance are themselves different
        # quantities and must not be read as one.
        detail = str(entry.get("target_set") or TARGET_SET_CORE)
        if detail not in TARGET_SET_DETAILS:
            raise SchemeCatalogError(
                f"{entry['key']} declares the target set {detail!r}, which this catalogue does not "
                f"know. It must be one of {', '.join(TARGET_SET_DETAILS)}: guessing would be "
                "guessing what a distance measured on it means.")
        entry["target_set_detail"] = detail
        entry["target_set"] = (TARGET_SET_CORE if detail == TARGET_SET_CORE
                               else TARGET_SET_ACCESSORY)
        entry["scheme_group"] = scheme_group(entry)
        entry["slot"] = slot_name(entry)
        entry["source_url"] = source_url(entry)
        entry["version"] = _version_of(entry)
        entry["title"] = scheme_title(entry)
        entry["notes"] = list(entry.get("notes", ()))
        rows.append(entry)
    return rows


def target_set_detail(entry) -> str:
    """Which of the three target sets a row is, defaulting to core for a row that
    says nothing -- because a cgMLST scheme with nothing said about it is the core
    set. An unrecognised value is returned as-is so a caller renders the words the
    row actually carries rather than quietly calling an unknown set a core one."""
    if not isinstance(entry, dict):
        return str(entry or TARGET_SET_CORE)
    return str(entry.get("target_set_detail") or entry.get("target_set") or TARGET_SET_CORE)


def target_set_label(entry) -> str:
    """'Core', 'Accessory' or 'Whole genome' -- one column's worth of the truth.

    A row whose target set this catalogue does not recognise reads 'Not recorded',
    never 'Core': an unlabelled set is an unknown quantity, not a core genome.
    """
    return _TARGET_SET_LABELS.get(target_set_detail(entry), "Not recorded")


def scheme_title(entry: dict) -> str:
    """One readable row title for a catalogued scheme, never a folder name.

    The unit is spelled out because it is the point: 'targets' here and 'loci' on a
    classical scheme are different quantities that must never share a scale. A set
    that is not the core set says so in the title as well, so a person choosing
    between two rows of a menu cannot mistake one for the other.
    """
    detail = target_set_detail(entry)
    parts = [_organism(entry), str(entry.get("scheme_name") or ""),
             _TARGET_SET_PHRASES.get(detail, "") if detail != TARGET_SET_CORE else "",
             f"{entry['locus_count']} targets",
             PROVIDERS[entry["provider"]]["name"] if entry.get("provider") in PROVIDERS else "",
             f"updated {_version_of(entry)}" if _version_of(entry) else ""]
    return " · ".join(part for part in parts if part)


def provider_token(value) -> str:
    """One canonical token per provider, whatever spelling a caller arrived with."""
    token = re.sub(r"[^A-Za-z0-9]+", "_", str(value or "")).strip("_").casefold()
    return _PROVIDER_TOKENS.get(token, token)


def scheme_group(entry: dict) -> str:
    """The provider-and-organism family whose core and accessory sets belong together."""
    organism = re.sub(r"[^A-Za-z0-9]+", "_", _organism(entry) or "unknown_organism").strip("_")
    return f"{provider_token(entry.get('provider') or entry.get('source'))}:{organism.casefold()}"


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


def _provider_of(descriptor) -> str:
    """One canonical provider token, preferring the service the URL actually names.

    A scheme records its provider as anything from "pubmlst" to "cgMLST.org
    Nomenclature Server (Ridom GmbH)". A free-text name that resolves to nothing
    this catalogue knows is not allowed to become a folder name: the host of the
    scheme's own API answers the same question and answers it the same way twice.
    """
    token = provider_token(descriptor.get("provider") or descriptor.get("source"))
    if token in _KNOWN_PROVIDERS:
        return token
    url = str(descriptor.get("url") or descriptor.get("API") or descriptor.get("source_url") or "")
    host = (urllib.parse.urlsplit(url).hostname or "").casefold()
    return _PROVIDER_HOSTS.get(host, token)


def slot_name(entry: dict) -> str:
    """The library folder for one scheme: organism first, then provider and target count.

    Named for what a microbiologist recognises. A folder called
    ``cgmlst_org_Kpneumoniae_abcdef0123456789`` tells a person nothing and sorts
    between two unrelated schemes; ``Klebsiella_pneumoniae__cgmlst_org_2358`` says
    the organism, who defined the targets and how many there are.
    """
    organism = re.sub(r"[^A-Za-z0-9]+", "_",
                      _organism(entry) or "Unknown_organism").strip("_") or "Unknown_organism"
    count = entry.get("locus_count") or 0
    return f"{organism}__{_provider_of(entry) or 'unknown_provider'}_{count}"


def _source_identity(descriptor) -> tuple[str, str, str]:
    """(provider token, database, scheme id) as the provider itself names them.

    Read from an online catalogue entry, a downloaded scheme.json or a snapshot
    manifest alike, because those are the three shapes an installed scheme's
    identity arrives in.
    """
    if not isinstance(descriptor, dict):
        return ("", "", "")
    url = str(descriptor.get("url") or descriptor.get("API") or descriptor.get("source_url") or "")
    path = urllib.parse.unquote(urllib.parse.urlsplit(url).path)
    provider = _provider_of(descriptor)
    database = str(descriptor.get("database") or "")
    scheme_id = str(descriptor.get("slug") or descriptor.get("scheme_id") or "")
    schema = re.search(r"/schema/([^/]+)/?$", path)
    bigsdb = re.search(r"/db/([^/]+)/schemes/([^/]+)/?$", path)
    if schema:
        provider = provider or "cgmlst_org"
        scheme_id = scheme_id or schema.group(1)
    elif bigsdb:
        database = database or bigsdb.group(1)
        scheme_id = scheme_id or bigsdb.group(2)
    return (provider, database, scheme_id)


def entry_for_source(descriptor) -> dict | None:
    """The pinned catalogue row a descriptor names, or None when it names none.

    Matched on the provider and on the provider's OWN scheme identifier -- never on
    the organism name and never on the target count. PubMLST and cgMLST.org both
    publish a 2,692-target Serratia marcescens scheme that shares no target name
    at all, so equal counts are a coincidence and must not resolve an identity.
    """
    provider, database, scheme_id = _source_identity(descriptor)
    if not provider or not scheme_id:
        return None
    for entry in catalog_entries():
        if provider_token(entry["provider"]) != provider or entry["scheme_id"] != scheme_id:
            continue
        if entry["database"] and database and entry["database"] != database:
            continue
        return entry
    return None


def slot_for(descriptor) -> str:
    """The library folder a scheme belongs in, catalogued or not.

    A scheme this catalogue pins lands in the slot the library already describes,
    so a download fills the labelled, licence-annotated folder a user has been
    looking at instead of creating a second one beside it.
    """
    pinned = entry_for_source(descriptor)
    return pinned["slot"] if pinned else slot_name(descriptor)


def scheme_variants(organism=None, *, group=None) -> dict:
    """The core, accessory and whole-genome target sets catalogued for one organism.

    A cgMLST scheme is a core set. Where a provider also publishes an accessory or
    pan-genome set, all of them are returned so the choice can be offered
    explicitly. Where only the core set is catalogued this says so in words,
    because an empty "accessory" list on its own reads as a broken menu rather than
    as an answer.

    "accessory" is every set that is not the core set, which is what a menu offering
    the alternative to a core run needs. "accessory_only" and "whole_genome" split
    it, because 251 accessory targets and 1,907 pan-genome targets are not one
    quantity either.
    """
    wanted = str(organism or "").strip().casefold()
    rows = [entry for entry in catalog_entries()
            if (group is None or entry["scheme_group"] == group)
            and (not wanted or entry["organism"].casefold() == wanted
                 or entry["genus"].casefold() == wanted)]
    core = [entry for entry in rows if entry["target_set"] == TARGET_SET_CORE]
    accessory = [entry for entry in rows if entry["target_set"] == TARGET_SET_ACCESSORY]
    accessory_only = [entry for entry in accessory
                      if entry["target_set_detail"] == TARGET_SET_ACCESSORY]
    whole_genome = [entry for entry in accessory
                    if entry["target_set_detail"] == TARGET_SET_WHOLE_GENOME]
    if not rows:
        message = (f"No cgMLST scheme is catalogued for {organism or group}. That is a gap in this "
                   "catalogue, not proof that no scheme exists.")
    elif not accessory:
        message = ("Only a core target set is catalogued for this organism. No accessory or "
                   "whole-genome set is pinned, so there is nothing to choose between: the core "
                   "set is the scheme.")
    elif not core:
        message = ("Only an accessory or whole-genome target set is catalogued for this organism. "
                   "Neither is a core genome scheme and neither measures the same quantity one "
                   "would.")
    else:
        counts = [f"{len(core)} core"]
        if accessory_only:
            counts.append(f"{len(accessory_only)} accessory")
        if whole_genome:
            counts.append(f"{len(whole_genome)} whole-genome")
        counted = (", ".join(counts[:-1]) + " and " + counts[-1]) if len(counts) > 1 else counts[0]
        message = (f"{counted} target set(s) are catalogued for this organism. They are different "
                   "quantities: a distance from one never shares a scale, an axis or a threshold "
                   "with a distance from another, and a cutoff published on the core set is not "
                   "offered for a core-plus-accessory run. Choose one and record which you chose.")
    return {"organism": str(organism or ""), "group": group, "core": core,
            "accessory": accessory, "accessory_only": accessory_only,
            "whole_genome": whole_genome, "has_core": bool(core),
            "has_accessory": bool(accessory), "has_whole_genome": bool(whole_genome),
            "message": message}


def catalog_digest() -> str:
    """A fingerprint of the pinned catalogue, for staging manifests and tests."""
    pinned = [{key: entry[key] for key in
               ("key", "provider", "database", "scheme_id", "locus_count",
                "target_list_sha256", "threshold_scheme_key", "target_set_detail")}
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
              "target_set": entry["target_set"], "target_set_detail": entry["target_set_detail"],
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
    # A publication that states no target count is read as having been derived on
    # the pinned CORE set, because that is what an unqualified cgMLST cutoff means.
    # It is never read that way for an accessory or pan-genome set: a number with no
    # stated target count cannot be shown to have been measured on that set, and a
    # core cutoff quietly offered for a core-plus-accessory run is precisely the
    # mistake this catalogue exists to prevent.
    numeric = [row for row in rows if row["published_threshold"] is not None
               and (row["locus_count"] == entry["locus_count"]
                    or (row["locus_count"] is None
                        and entry["target_set"] == TARGET_SET_CORE))]
    if not numeric:
        result["status"] = "citation_only"
        result["reason"] = ("The bound publication supplies scope and citation but no transferable "
                            "number for this target set.")
        if entry["target_set"] != TARGET_SET_CORE:
            result["reason"] += (
                " This is not a core genome scheme: it is the "
                f"{_TARGET_SET_PHRASES.get(entry['target_set_detail'], 'non-core target set')} of "
                f"{entry['locus_count']} targets, so a cutoff is offered only by a publication "
                "that states this exact target count.")
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

Put the {locus_count}-target {provider_name} {set_phrase} in THIS folder.

{state_line}

Provider   : {provider_name}
Scheme     : {source_url}
Targets    : {locus_count}
Target set : {set_label}
Terms      : {terms_url}

{restriction}

This folder is created by WMLSTudio so that you never have to make one, and a
download of this scheme installs into THIS folder rather than somewhere else.
Deleting it removes these description files and whatever scheme you installed
here; the app recreates the empty, labelled folder on the next run.

{set_note}"""

# What a distance from this folder means depends on which target set is in it, so
# each set gets the closing paragraph it actually needs rather than one written for
# a core scheme and then contradicted.
_SET_NOTES = {
    TARGET_SET_CORE: """A cgMLST distance from this scheme is a number of differing targets out of the
targets called in BOTH isolates. It is not a seven-locus MLST distance, it is not
a SNP count, and it is never proof of transmission.
""",
    TARGET_SET_ACCESSORY: """This is an ACCESSORY target set, not a cgMLST scheme. Its targets are absent from
some isolates of this organism by design, so a target that is not called here is
biology and not a failed call: it is reported as not assayed, never as a
difference. A distance measured on this set never shares a scale, an axis or a
published cutoff with a core-genome distance, the two are never added together
into one number, and neither is proof of transmission.
""",
    TARGET_SET_WHOLE_GENOME: """This is a WHOLE-GENOME (pan-genome) target set: core and accessory targets in one
set. Run it INSTEAD of a core scheme, never beside one as something added to a
core distance. A distance measured on it is not a cgMLST distance, never shares a
scale, an axis or a published cutoff with one, and is never proof of transmission.
""",
}

_BUNDLED_STATE = ("This scheme's provider permits redistribution, so the release may already "
                  "carry it. If the folder is empty, use Reference data > Install cgMLST scheme.")
_DOWNLOAD_STATE = ("This scheme CANNOT be packed into WMLSTudio: its provider does not grant "
                   "redistribution. Use Reference data > Download cgMLST scheme. You will be "
                   "shown the provider's terms and asked to confirm that your use is permitted "
                   "before anything is downloaded.")
_INSTALLED_STATE = ("This scheme IS installed in this folder. The allele files beside this README "
                    "are the installed target set; reference_manifest.json records where every "
                    "byte came from and which allele records, if any, were excluded.")


def library_root(root) -> Path:
    """<data root>/cgmlst -- the one place installed cgMLST schemes live."""
    return Path(root).expanduser() / LIBRARY_DIRNAME


def classical_root(root) -> Path:
    """<data root>/schemes -- where the classical seven-locus schemes live instead."""
    return Path(root).expanduser() / SCHEME_DIRNAME


def install_root(library_root_path, kind: str = "cgmlst") -> Path:
    """Where a scheme of this kind installs, given whichever library a caller named.

    There is ONE cgMLST library. A caller that hands over the classical scheme
    folder is not installing a 2,000-target scheme into it: the sibling cgMLST
    library is returned instead, so a downloaded cgMLST scheme lands in the library
    a user browses for cgMLST schemes rather than between two seven-locus ones.
    """
    base = Path(library_root_path).expanduser()
    if kind != "cgmlst":
        return base
    if base.name == LIBRARY_DIRNAME:
        return base
    if base.name == SCHEME_DIRNAME:
        return base.parent / LIBRARY_DIRNAME
    return base / LIBRARY_DIRNAME


def is_gene_by_gene(descriptor) -> bool:
    """Whether a scheme descriptor names a gene-by-gene scheme, not a classical one.

    One definition for the whole application: the declared target count decides
    once it clears the floor, and below the floor a declared cgMLST/wgMLST type
    decides. Nothing else is consulted, because the organism name and the provider
    say nothing about how many targets a scheme has.
    """
    if not isinstance(descriptor, dict):
        return False
    declared = str(descriptor.get("type") or "").strip().casefold().replace(" ", "")
    try:
        count = int(descriptor.get("locus_count") or 0)
    except (TypeError, ValueError):
        count = 0
    return count > CGMLST_TARGET_FLOOR or declared in {"cgmlst", "wgmlst", "core", "accessory"}


def has_alleles(folder) -> bool:
    """Whether a folder holds an installed scheme rather than a labelled empty slot."""
    folder = Path(folder)
    try:
        for path in folder.iterdir():
            if not path.is_file():
                continue
            name = path.with_suffix("") if path.suffix.casefold() in {".gz", ".bz2"} else path
            if name.suffix.casefold() in _ALLELE_SUFFIXES:
                return True
    except OSError:
        return False
    return False


def _folder_digest(folder: Path) -> str:
    """The scheme digest a snapshot manifest records for an installed folder, or ''."""
    path = folder / "reference_manifest.json"
    if not path.is_file() or path.stat().st_size > 4 * 1024 * 1024:
        return ""
    try:
        manifest = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, ValueError):
        return ""
    return str(manifest.get("scheme_digest") or "") if isinstance(manifest, dict) else ""


def install_folder(base, descriptor, *, digest="") -> Path:
    """The folder this scheme installs into: its readable slot, never a digest name.

    An empty, pre-created slot IS the destination -- that is what the labelled
    folder is for. A slot already holding a DIFFERENT scheme is never overwritten:
    a second, digest-suffixed folder is returned instead, so an upstream
    redefinition that keeps the same target count cannot silently replace the
    target set a result was called against.
    """
    base = Path(base).expanduser()
    slot = slot_for(descriptor)
    candidate = base / slot
    if not candidate.is_dir() or not has_alleles(candidate):
        return candidate
    if not digest or _folder_digest(candidate) == digest:
        # The same scheme, or a caller with no digest to tell them apart; the
        # caller checks whether the folder is occupied before writing to it.
        return candidate
    return base / f"{slot}__{digest[:12]}"


def slot_payload(entry: dict) -> dict:
    """The machine-readable description written into an empty slot."""
    return {"format_version": SLOT_FORMAT_VERSION, "catalog_version": CATALOG_VERSION,
            "key": entry["key"], "organism": entry["organism"], "genus": entry["genus"],
            "species": entry["species"], "kind": "cgmlst", "target_set": entry["target_set"],
            "target_set_detail": entry["target_set_detail"],
            "target_set_label": target_set_label(entry),
            "scheme_group": entry["scheme_group"], "title": entry["title"],
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


def _descriptor(folder: Path) -> dict:
    """What an installed folder records about itself, from the files it carries."""
    descriptor: dict = {}
    for name in ("scheme.json", "reference_manifest.json", SLOT_FILENAME):
        path = folder / name
        if not path.is_file() or path.stat().st_size > 4 * 1024 * 1024:
            continue
        try:
            value = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, UnicodeDecodeError, ValueError):
            continue
        if not isinstance(value, dict):
            continue
        if name == "reference_manifest.json":
            value = value.get("entry") if isinstance(value.get("entry"), dict) else {}
        for key, item in (value or {}).items():
            descriptor.setdefault(key, item)
    return descriptor


def _is_cgmlst_folder(folder: Path) -> bool:
    """Whether an installed folder holds a gene-by-gene scheme, on its own evidence.

    The installed target count decides once it clears the floor; below the floor a
    declared cgMLST/wgMLST type decides, because a partial snapshot of a
    gene-by-gene scheme is still one. A folder that records neither is left where
    it is rather than moved on a guess.
    """
    count = _installed_locus_count(folder)
    if not count:
        return False
    return is_gene_by_gene({**_descriptor(folder), "locus_count": count})


def _read_migrations(base: Path) -> dict:
    path = base / MIGRATIONS_FILENAME
    if not path.is_file() or path.stat().st_size > 4 * 1024 * 1024:
        return {"format_version": 1, "moved": []}
    try:
        stored = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, ValueError):
        return {"format_version": 1, "moved": []}
    if not isinstance(stored, dict) or not isinstance(stored.get("moved"), list):
        return {"format_version": 1, "moved": []}
    return stored


def migrate_downloads(root, *, sources=None) -> list[dict]:
    """Move cgMLST schemes out of <data root>/schemes and into the cgMLST library.

    Earlier releases installed every downloaded scheme into <data root>/schemes,
    where a 2,358-target cgMLST scheme sat between two seven-locus MLST schemes
    under a digest-shaped folder name -- found by nothing a person would look for
    and offered where a classical scheme was expected. The bytes are MOVED, never
    copied and never deleted, into the labelled slot the library already describes.

    Each move is recorded in migrations.json so a project that stored the old
    scheme_path can still be pointed at the scheme it was actually called against.
    A folder is left exactly where it is unless it clearly records itself as a
    gene-by-gene scheme, and an occupied destination is never overwritten.
    """
    base = library_root(root)
    bases = [Path(path).expanduser() for path in sources] if sources is not None \
        else [classical_root(root)]
    moved = []
    for source in bases:
        if not source.is_dir() or source.resolve() == base.resolve():
            continue
        for folder in sorted((p for p in source.iterdir() if p.is_dir()),
                             key=lambda p: p.name.casefold()):
            if folder.name.startswith(("_", ".")) or not _is_cgmlst_folder(folder):
                continue
            descriptor = _descriptor(folder)
            descriptor.setdefault("locus_count", _installed_locus_count(folder))
            destination = install_folder(base, descriptor, digest=_folder_digest(folder))
            if destination.exists() and has_alleles(destination):
                continue
            base.mkdir(parents=True, exist_ok=True)
            try:
                install_into(folder, destination)
            except OSError:
                continue  # A locked or in-use folder stays put; nothing is lost.
            moved.append({"from": str(folder), "to": str(destination),
                          "locus_count": int(descriptor.get("locus_count") or 0),
                          "moved_utc": datetime.now(UTC).isoformat()})
    if moved:
        ledger = _read_migrations(base)
        ledger["moved"] = [row for row in ledger["moved"] if isinstance(row, dict)] + moved
        ledger["format_version"] = 1
        (base / MIGRATIONS_FILENAME).write_text(
            json.dumps(ledger, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return moved


def install_into(source: Path, destination: Path, *, progress=None) -> None:
    """Move a scheme folder into place, filling a pre-created empty slot if there is one.

    A labelled slot already carries its README and scheme_slot.json, so the slot is
    filled file by file rather than replaced; those two descriptions are what tell
    a person whose scheme this is and under what terms it was obtained.

    Filling a slot is thousands of moves, and a copying move across a filesystem
    boundary carries gigabytes, so it reports each file rather than going quiet
    at the very last step of a long download.
    """
    if not destination.exists():
        if progress:
            progress(0, 1, f"Moving the scheme into {destination.name}")
        try:
            source.rename(destination)
        except OSError:
            shutil.move(str(source), str(destination))
        if progress:
            progress(1, 1, f"Moved the scheme into {destination.name}")
        return
    children = sorted(source.iterdir())
    for number, child in enumerate(children, 1):
        target = destination / child.name
        if progress:
            progress(number, len(children), f"Filing {child.name}")
        if target.exists():
            continue
        try:
            child.rename(target)
        except OSError:
            shutil.move(str(child), str(target))
    remaining = sorted(child.name for child in source.iterdir())
    if remaining:
        # Nothing already in the destination is replaced. A caller that cannot
        # complete the move keeps its staging folder and publishes nothing.
        raise OSError(f"{destination} already holds {', '.join(remaining)}; nothing was replaced.")
    source.rmdir()


def resolve_migrated_path(root, path) -> Path:
    """Where a scheme folder went, for a stored path that names its old location.

    Returns the path unchanged when it was never moved, so a caller can route every
    stored scheme_path through this without deciding first.
    """
    wanted = str(Path(path).expanduser())
    for row in reversed(_read_migrations(library_root(root))["moved"]):
        if isinstance(row, dict) and str(row.get("from") or "") == wanted:
            return Path(str(row.get("to") or wanted))
    return Path(path)


def prepare_library(root, *, keys=None, targets=None, migrate=True) -> dict:
    """Create the cgMLST library layout so a first run never asks for a folder.

    Every catalogued scheme gets a named folder carrying a README and a
    scheme_slot.json; a download installs INTO that folder, so the layout a user
    browses and the place a scheme actually lands are one and the same. Existing
    folders and any scheme already installed in them are left untouched: only the
    description files are refreshed, and only when their content changed.
    `targets` optionally maps a key to its target-name list, written as targets.txt
    beside the slot.

    With `migrate` (the default) any cgMLST scheme still sitting in the classical
    <data root>/schemes folder from an earlier release is moved into the library
    first, so a scheme downloaded yesterday is where it is now looked for.
    """
    base = library_root(root)
    base.mkdir(parents=True, exist_ok=True)
    wanted = set(keys) if keys is not None else None
    report = {"root": str(base), "created": [], "existing": [], "refreshed": [], "migrated": [],
              "catalog_version": CATALOG_VERSION, "catalog_digest": catalog_digest(),
              "generated_utc": datetime.now(UTC).isoformat()}
    if migrate:
        report["migrated"] = migrate_downloads(root)
    for entry in catalog_entries():
        if wanted is not None and entry["key"] not in wanted:
            continue
        folder = base / entry["slot"]
        (report["existing"] if folder.is_dir() else report["created"]).append(entry["slot"])
        folder.mkdir(parents=True, exist_ok=True)
        if has_alleles(folder):
            state_line = _INSTALLED_STATE
        else:
            state_line = _BUNDLED_STATE if entry["bundled"] else _DOWNLOAD_STATE
        detail = entry["target_set_detail"]
        readme = SLOT_README.format(
            organism=entry["organism"], scheme_name=entry["scheme_name"],
            locus_count=entry["locus_count"], provider_name=entry["provider_name"],
            source_url=entry["source_url"], terms_url=entry["terms_url"],
            restriction=entry["licence_restriction"], state_line=state_line,
            set_phrase=_TARGET_SET_PHRASES.get(detail, "target set"),
            set_label=target_set_label(entry), set_note=_SET_NOTES.get(detail, ""))
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

Every cgMLST scheme WMLSTudio installs lives here, in a folder named for its
organism, its provider and its target count -- never a folder named after a
digest. One folder is pre-created for each catalogued scheme and explains, in its
own README.txt, who publishes it, under what terms, and how to install it. An
empty folder means that scheme is not installed yet; it does not mean anything is
broken. A scheme this catalogue does not pin gets its own folder here too, named
the same way.

These are cgMLST schemes: hundreds to thousands of targets. Classical seven-locus
MLST schemes live in the separate schemes folder. Distances from the two are
different quantities and are never mixed, compared or thresholded together.

Each folder's own README names its target set: Core, Accessory or Whole genome.
A core set is expected in every isolate; an accessory set is not, and a
whole-genome set is the two together. They too are different quantities. A
distance from one never shares a scale, an axis or a published cutoff with a
distance from another, and they are never added together into one number.

migrations.json, if present, records schemes moved here from the classical
schemes folder by an earlier release, so a saved project that stored the old
location can still be pointed at the scheme it used.

Deleting a folder here deletes its description files and any scheme installed
into it. WMLSTudio recreates the empty, labelled layout on the next run.
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


def _folder_metadata(folder: Path) -> dict:
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
    return metadata


def library_index(root, *, extra_paths=()) -> list[dict]:
    """Every folder in the cgMLST library that actually holds a scheme, read once.

    Counted and identified, never parsed. Building this once and handing it to the
    per-scheme lookups is what keeps a twenty-three row library page from walking a
    four-thousand-file folder twenty-three times.
    """
    base = library_root(root)
    folders = [path for path in sorted(base.iterdir(), key=lambda p: p.name.casefold())
               if path.is_dir()] if base.is_dir() else []
    folders.extend(Path(path).expanduser() for path in extra_paths)
    seen, rows = set(), []
    for folder in folders:
        resolved = str(folder)
        if resolved in seen or not folder.is_dir():
            continue
        seen.add(resolved)
        count = _installed_locus_count(folder)
        if not count:
            continue
        rows.append({"path": folder, "locus_count": count,
                     "metadata": _folder_metadata(folder)})
    return rows


def installed_scheme(root, key, *, extra_paths=(), index=None) -> dict | None:
    """The scheme installed for one catalogue key, or None. Counts files, parses none.

    The whole library is searched, not only the folder this scheme is meant to sit
    in: a scheme installed under any name is still that scheme, and a user who
    cannot find what they downloaded is the bug this exists to prevent. A folder is
    accepted as this scheme only when its recorded provider and scheme identifier
    match the pin, so the named slot is a preference and never a proof of identity.
    A folder with the right number of allele files but a different or absent
    identity is reported with `identity_matched` False, so a caller can show it
    without ever binding a threshold to it.
    """
    entry = entry_for(key)
    rows = library_index(root, extra_paths=extra_paths) if index is None else index
    slot = library_root(root) / entry["slot"]
    matched = [row for row in rows if _identity_matches(entry, row["metadata"])]
    if matched:
        row = next((item for item in matched if item["path"] == slot), matched[0])
        return _installed_row(entry, row, identity_matched=True)
    unnamed = next((row for row in rows if row["path"] == slot), None)
    return _installed_row(entry, unnamed, identity_matched=False) if unnamed else None


def _installed_row(entry: dict, row: dict, *, identity_matched: bool) -> dict:
    count = row["locus_count"]
    return {"key": entry["key"], "path": str(row["path"]), "locus_count": count,
            "expected_locus_count": entry["locus_count"],
            "count_matched": count == entry["locus_count"],
            "identity_matched": identity_matched,
            "organism": entry["organism"], "scheme_name": entry["scheme_name"],
            "provider_name": entry["provider_name"], "target_set": entry["target_set"],
            "target_set_detail": entry["target_set_detail"],
            "target_set_label": target_set_label(entry)}


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
    index = library_index(root, extra_paths=extra_paths)
    for entry in catalog_entries():
        installed = installed_scheme(root, entry["key"], index=index)
        count = installed["locus_count"] if installed else None
        guidance = threshold_for(entry["key"], locus_count=count)
        rows.append({**entry, "folder": str(library_root(root) / entry["slot"]),
                     "installed": installed, "ready": bool(installed and installed["count_matched"]
                                                           and installed["identity_matched"]),
                     "threshold": guidance})
    return rows


def installed_entries(root, *, extra_paths=(), cancelled=None) -> list[dict]:
    """Every cgMLST scheme actually installed, catalogued or not, with a readable row.

    This is what a "cgMLST schemes" tab lists: what the user HAS, titled by organism
    and provider rather than by folder name. A scheme this catalogue does not pin is
    listed all the same, with `catalog_key` None and no bound cutoff -- an
    unrecognised scheme is still perfectly usable within itself.
    """
    from .reference_index import scheme_entries as _scheme_entries
    index = library_index(root, extra_paths=extra_paths)
    rows = _scheme_entries([row["path"] for row in index], cancelled=cancelled, kind="cgmlst")
    listed = []
    for row in rows:
        identified = identify_installed(row["path"])
        guidance = threshold_for_installed(row["path"])
        detail = (identified["target_set_detail"] if identified else "") or row["target_set"]
        listed.append({**row, "catalog_key": identified["key"] if identified else None,
                       "catalogued": identified is not None,
                       "target_set": row["target_set"] or (
                           identified["target_set"] if identified else ""),
                       "target_set_detail": detail,
                       "target_set_label": target_set_label(detail) if detail else "Not recorded",
                       "expected_locus_count": identified["locus_count"] if identified else None,
                       "count_matched": bool(identified)
                       and row["locus_count"] == identified["locus_count"],
                       "threshold": guidance})
    return listed


def identify_installed(path) -> dict | None:
    """Which catalogued scheme a scheme folder is, read from what it records itself.

    Returns None when the folder records no identity this catalogue recognises.
    That is the honest answer: an unrecognised scheme is still perfectly usable for
    comparison, it simply has no bound cutoff.
    """
    folder = Path(path).expanduser()
    if not folder.is_dir():
        return None
    metadata = _folder_metadata(folder)
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
    return {"key": entry["key"], "organism": entry["organism"], "title": entry["title"],
            "scheme_name": entry["scheme_name"], "locus_count": entry["locus_count"],
            "provider": entry["provider"], "provider_name": entry["provider_name"],
            "target_set": entry["target_set"], "target_set_detail": entry["target_set_detail"],
            "target_set_label": target_set_label(entry), "scheme_group": entry["scheme_group"],
            # What the person is about to download decides what the numbers it
            # produces mean, so the download button says it before the first byte
            # rather than the report saying it afterwards.
            "target_set_notice": _SET_NOTES.get(entry["target_set_detail"], "").strip(),
            "source_url": entry["source_url"], "terms_url": provider["terms_url"],
            "terms_notice": provider["restriction"], "licence_quote": provider["quote"],
            "requires_acknowledgement": bool(provider["requires_acknowledgement"]),
            "may_be_bundled": bool(provider["may_bundle"]), "method": method,
            "destination_slot": entry["slot"],
            "destination_library": LIBRARY_DIRNAME, "notes": list(entry["notes"])}
