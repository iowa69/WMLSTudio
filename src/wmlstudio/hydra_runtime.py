"""Cancellable, native assembly execution of the pinned upstream HYDRA engine.

The wrapper does not reimplement any scientific calling logic. It runs HYDRA
in a dedicated child process, never auto-downloads databases during analysis,
and imports the engine's JSON with execution and reference provenance.

It also answers what a finished run could and could not report: which reference
sets supply which kind of element, what each element type found for one isolate,
and which of those blanks are unasked questions rather than clean results. That
last part lives here rather than in a view because it is a statement about the
reference data, and because every view that showed only the resistance genes was
showing an isolate's virulence, stress and plasmid evidence as nothing at all.
"""

from __future__ import annotations

import html
import importlib.metadata
import json
import math
import os
import re
import shutil
import signal
import subprocess
import sys
import tempfile
import time
import urllib.request
import uuid
from datetime import UTC, datetime
from pathlib import Path

from wmlstudio.hydra import ELEMENT_TITLES, ELEMENT_TYPES, is_point_mutation, load_hydra_report
from wmlstudio.sequence import (
    AnalysisCancelled,
    SequenceReader,
    file_sha256,
    file_signature,
    sample_name,
)

HYDRA_COMMIT = "6d36c109491c16544e8919fe6962b4b62e97d3d7"
HYDRA_VERSION = "1.4.0"
TOOLS = ("blastn", "blastx", "blastp", "tblastn", "makeblastdb")
NCBI_REFERENCE_ROOT = "https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/latest"
# The stores a run uses when the caller names none, in the order a report reads.
# These two ship with the application, so a first run works with no download.
DEFAULT_DATABASES = ("ncbi", "protein")
# What each store is for, said in the words a microbiologist reads. A run without
# one of these does not fail: it returns nothing for that whole class of evidence,
# which reads exactly like a negative result. Every refusal below names the store
# and this sentence, so "nothing was found" is never confused with "nothing ran".
# The names are the engine's own; these sentences are ours, because a registry
# entry titled "MEGARes" tells a microbiologist nothing about what it answers.
DATABASE_PURPOSE = {
    "ncbi": ("acquired resistance, stress and virulence genes, searched with blastn against the "
             "NCBI AMRFinderPlus nucleotide catalogue"),
    "protein": ("the translated protein search and every organism point-mutation catalogue, "
                "searched with blastx against AMRProt"),
    "card": ("a second, broader opinion on acquired resistance, including the efflux pumps and "
             "regulators the NCBI catalogue deliberately leaves out; a hit here that NCBI does "
             "not report is a difference between two catalogues, not a stronger finding"),
    "resfinder": ("the acquired-resistance gene set behind the CGE tools, for comparing a call "
                  "against what a ResFinder-based report would have named"),
    "argannot": ("an older acquired-resistance gene set, useful when a gene name in a paper "
                 "predates the current NCBI catalogue"),
    "megares": ("resistance together with biocide and heavy-metal determinants, for disinfectant "
                "and metal-exposure questions the AMR sets do not cover"),
    "vfdb": "the curated core virulence factors of VFDB (setA)",
    "vfdb_full": ("every predicted virulence factor VFDB publishes (setB): a much broader screen "
                  "whose hits are hypotheses rather than curated calls"),
    "ecoli_vf": ("Escherichia coli virulence factors, for pathotype questions the general "
                 "virulence sets do not answer"),
    "plasmidfinder": ("plasmid replicon types, which suggest what kind of plasmid a contig may "
                      "belong to; a replicon is not a plasmid, a location or a transfer event"),
    "ecoh": "Escherichia coli O and H antigen loci, for serotype prediction from sequence",
    "pubmlst": ("the engine's own 7-locus MLST schemes; WMLSTudio types MLST itself from its own "
                "scheme library, so this is only needed to reproduce the engine's own ST call"),
    "lineage": ("Kleborate-derived lineage and sublineage loci for the Klebsiella complex; this "
                "is a screen over that project's published reference data, not Kleborate"),
    "sccmec": ("whole SCCmec cassette references for staphylococci, read by the engine's typing "
               "step rather than by gene screening"),
    "species": ("optional Mash sketches for the engine's own species guess; WMLSTudio identifies "
                "organisms from its own species panel and does not need this"),
}
# Who publishes each set. The engine's registry records an address and a citation
# but not a name a reader would recognise, and "whose data is this" is the first
# question asked of any screen that is not the tool it screens for.
DATABASE_PROVIDER = {
    "ncbi": "NCBI, US National Library of Medicine",
    "protein": "NCBI, US National Library of Medicine",
    "card": "McMaster University (CARD)",
    "resfinder": "Center for Genomic Epidemiology, DTU",
    "argannot": "IHU Méditerranée Infection",
    "megares": "Colorado State University (MEGARes)",
    "vfdb": "Institute of Pathogen Biology, CAMS (VFDB)",
    "vfdb_full": "Institute of Pathogen Biology, CAMS (VFDB)",
    "ecoli_vf": "Public Health Agency of Canada",
    "plasmidfinder": "Center for Genomic Epidemiology, DTU",
    "ecoh": "SRST2 / Holt laboratory",
    "pubmlst": "PubMLST, University of Oxford",
    "lineage": "Kleborate, Holt laboratory",
    "sccmec": "Center for Genomic Epidemiology, DTU",
    "species": "Kleborate, Holt laboratory",
}
# Which reference sets can report which kind of element. Without this the
# application cannot tell a clean isolate from an unasked question: a store
# holding only the two core sets reports no plasmid replicon for any isolate on
# earth, and "0 plasmids" is what the reader sees.
#
# The engine's registry records one headline element type per set and that is
# where every row below comes from, with one measured exception: the two
# AMRFinderPlus catalogues are registered as AMR but their meta tables carry
# virulence and stress rows as well (8794 AMR, 1025 VIRULENCE and 259 STRESS in
# the bundled 2026-08-07.1 protein snapshot, and the same 1025/259 in the
# nucleotide one). The nucleotide engine applies no element-type filter, so
# those virulence and stress genes are reported whatever --plus is set to; only
# the translated protein search is gated by it. Nothing is credited with an
# element type nobody here has counted.
ELEMENT_SOURCES = {
    "AMR": ("ncbi", "protein", "card", "resfinder", "argannot", "megares"),
    "VIRULENCE": ("ncbi", "protein", "vfdb", "vfdb_full", "ecoli_vf"),
    "STRESS": ("ncbi", "protein"),
    "PLASMID": ("plasmidfinder",),
}
# How large a set is. Only four figures exist here, and each says where it came
# from: the two core sets are measured from the copy bundled with this release,
# and two more are the estimate the engine's own registry records. No provider
# publishes a size for the rest, so those rows say so and report a measured size
# once the set is installed, rather than showing a number nobody checked.
DATABASE_SIZE = {
    "ncbi": (14, "measured from the copy bundled with this release"),
    "protein": (10, "measured from the copy bundled with this release, point-mutation "
                    "catalogues included"),
    "vfdb_full": (70, "the estimate the engine's registry records"),
    "lineage": (70, "the engine's registry records one ~70 MB archive shared with the species "
                    "sketches"),
    "species": (70, "the engine's registry records one ~70 MB archive shared with the lineage "
                    "loci"),
}
# A licence a reader can act on without asking anyone. Anything else — including a
# licence the provider never recorded — is listed, downloadable on request, and
# never swept into an "install everything" that nobody read the terms for.
OPEN_LICENCE_TOKENS = ("public domain", "apache", "cc0", "cc-0", "mit", "bsd", "creative commons")
# What each executable is for, so a missing one is a sentence and not a filename.
TOOL_PURPOSE = {
    "blastn": "the nucleotide gene search",
    "blastx": "the translated protein search, which is where point mutations are read",
    "blastp": "protein-to-protein confirmation of a translated hit",
    "tblastn": "recovering a gene that the translated search split across a contig break",
    "makeblastdb": "building the search index a reference store needs before it can be read",
}
# A snapshot older than this still runs and is still reported as the exact release
# it is; the age is surfaced so a user can decide whether to update before
# reporting. Determinants named after this date are simply not in it.
REFERENCE_AGE_DAYS = 180


class HydraRuntimeError(ValueError):
    """A readable configuration, execution, or provenance failure."""


def tool_directories():
    roots = [Path(sys.executable).resolve().parent]
    if getattr(sys, "frozen", False):
        roots.insert(0, Path(sys._MEIPASS))
    directories = [root / relative for root in roots
                   for relative in ("Tools/blast/bin", "Tools/blast+/bin", "Tools/bin")]
    directories.append(Path(__file__).resolve().parent / "resources/tools/blast/bin")
    return directories


def bundled_database_root():
    """Return the explicitly staged starter snapshot, or None in a source-only checkout."""
    candidates = [Path(__file__).resolve().parent / "resources/hydra/starter"]
    if getattr(sys, "frozen", False):
        candidates.insert(0, Path(sys._MEIPASS) / "wmlstudio/resources/hydra/starter")
    return next((path for path in candidates if (path / "manifest.json").is_file()), None)


def resolve_tool(name, explicit=None):
    if name not in TOOLS:
        raise HydraRuntimeError(f"Unsupported native tool: {name}")
    filename = name + (".exe" if os.name == "nt" else "")
    if explicit is not None:
        path = Path(explicit).resolve()
        if not path.is_file():
            raise HydraRuntimeError(f"Tool was not found: {path}")
        return str(path)
    for directory in tool_directories():
        if (directory / filename).is_file():
            return str(directory / filename)
    return shutil.which(filename)


def installed_databases(db_root):
    """Read the upstream manifest without modifying or creating a database store."""
    root = Path(db_root).resolve()
    path = root / "manifest.json"
    if not path.is_file():
        return {}
    try:
        manifest = json.loads(path.read_text(encoding="utf-8"))
        records = manifest.get("databases", {})
        if not isinstance(records, dict):
            raise ValueError("databases must be an object")
        result = {}
        for name, entry in records.items():
            if not isinstance(entry, dict) or not isinstance(entry.get("path"), str):
                raise ValueError(f"invalid entry for {name}")
            # A store staged on Windows records "nucl\ncbi"; on POSIX that is one
            # filename, not a folder, and the whole store would read as empty.
            # Normalising here makes one staged snapshot portable, and tightens
            # the containment check rather than loosening it: "..\\outside" is
            # now recognised as an escape instead of a legal filename.
            target = (root / entry["path"].replace("\\", "/")).resolve()
            if not target.is_relative_to(root):
                raise ValueError(f"database path leaves the selected store: {name}")
            if target.exists():
                result[name] = dict(entry)
        return result
    except (OSError, ValueError, TypeError, AttributeError) as exc:
        raise HydraRuntimeError(f"Cannot read HYDRA database manifest: {exc}") from exc


def _entry_directory(root, entry):
    """The folder one manifest entry names, with the same separator normalisation."""
    return (Path(root) / str(entry.get("path", "")).replace("\\", "/")).resolve()


def _release_date(version):
    """The NCBI release date inside a version like '2026-08-07.1', or None."""
    match = re.match(r"([0-9]{4}-[0-9]{2}-[0-9]{2})", str(version or ""))
    if not match:
        return None
    try:
        return datetime.strptime(match.group(1), "%Y-%m-%d").replace(tzinfo=UTC)
    except ValueError:
        return None


def _first_column(path, *, skip_header=True):
    """The first tab-separated column of a reference table, comments excluded."""
    values = set()
    try:
        with path.open(encoding="utf-8") as handle:
            if skip_header:
                handle.readline()
            for line in handle:
                value = line.split("\t")[0].strip()
                if value and not value.startswith("#"):
                    values.add(value)
    except OSError as exc:
        raise HydraRuntimeError(f"Cannot read the organism table: {exc}") from exc
    return values


def organism_catalogue(db_root):
    """The organism names the installed store will accept, and what each one buys.

    Upstream resolves --organism against its own taxgroup table, its protein
    mutation table and its DNA mutation catalogues, so this reads the same three
    places rather than inventing a list. They are three different things and are
    reported separately:

    "accepted" is every name the engine will take. "protein_point_mutations" is
    the set with curated protein-level mutations in AMRProt-mutation.tsv;
    "dna_point_mutations" is the smaller set that also has a DNA catalogue staged
    for it. "point_mutations" is their union, kept because an isolate outside it
    is screened for genes only.

    An organism outside "accepted" stops the run upstream. One inside it but
    outside "point_mutations" is accepted and simply has no mutation catalogue to
    report from, which is not the same thing as having no mutations.
    """
    root = Path(db_root).resolve() if db_root is not None else None
    empty = {"accepted": [], "point_mutations": [], "dna_point_mutations": [],
             "protein_point_mutations": [], "suppressed": [], "root": str(root or "")}
    if root is None or not root.is_dir():
        return empty
    try:
        entries = installed_databases(root)
    except HydraRuntimeError:
        entries = {}
    dna = {name.strip() for name in (entries.get("protein", {}).get("organisms") or [])
           if isinstance(name, str) and name.strip()}
    mutation = root / "mutation" / "dna"
    if mutation.is_dir():
        dna.update(path.stem for path in mutation.glob("*.fna"))
    protein_entry = entries.get("protein")
    directory = _entry_directory(root, protein_entry) if protein_entry else None
    protein, accepted, suppressed = set(), set(), set()
    if directory is not None:
        if (directory / "AMRProt-mutation.tsv").is_file():
            protein = _first_column(directory / "AMRProt-mutation.tsv")
        if (directory / "AMRProt-suppress.tsv").is_file():
            suppressed = _first_column(directory / "AMRProt-suppress.tsv")
        if (directory / "taxgroup.tsv").is_file():
            accepted = _first_column(directory / "taxgroup.tsv")
    accepted |= dna | protein
    return {"accepted": sorted(accepted), "point_mutations": sorted(dna | protein),
            "dna_point_mutations": sorted(dna), "protein_point_mutations": sorted(protein),
            "suppressed": sorted(suppressed), "root": str(root)}


# Counting element types means reading a 1.5 MB table, and the same store is asked
# about on every preflight and every refresh of the database page. The answer only
# changes when the file does, so it is cached against the file's own identity.
_ELEMENT_COUNTS = {}


def element_counts(db_root):
    """How many acquired-gene, virulence and stress elements the protein set holds.

    This is what decides whether offering virulence would offer anything at all.
    A store whose protein reference is absent returns zeros and says so through
    the caller, rather than letting an unperformed search look like a clean result.
    """
    root = Path(db_root).resolve() if db_root is not None else None
    counts = {"AMR": 0, "VIRULENCE": 0, "STRESS": 0, "point": 0, "total": 0, "read": False}
    if root is None or not root.is_dir():
        return counts
    try:
        entries = installed_databases(root)
    except HydraRuntimeError:
        return counts
    entry = entries.get("protein")
    if entry is None:
        return counts
    path = _entry_directory(root, entry) / "meta.tsv"
    if not path.is_file():
        return counts
    try:
        stamp = path.stat()
        key = (str(path), stamp.st_mtime_ns, stamp.st_size)
        if key in _ELEMENT_COUNTS:
            return dict(_ELEMENT_COUNTS[key])
        with path.open(encoding="utf-8") as handle:
            columns = handle.readline().rstrip("\n").split("\t")
            kind = columns.index("element_type") if "element_type" in columns else None
            subtype = columns.index("element_subtype") if "element_subtype" in columns else None
            for line in handle:
                fields = line.rstrip("\n").split("\t")
                if kind is None or len(fields) <= kind:
                    continue
                counts["total"] += 1
                counts[fields[kind]] = counts.get(fields[kind], 0) + 1
                if subtype is not None and len(fields) > subtype and fields[subtype] == "POINT":
                    counts["point"] += 1
    except (OSError, ValueError) as exc:
        raise HydraRuntimeError(f"Cannot read the protein reference table: {exc}") from exc
    counts["read"] = True
    _ELEMENT_COUNTS[key] = dict(counts)
    return counts


def virulence_support(db_root=None, *, organism=None, catalogue=None, counts=None):
    """Whether this store can report virulence and stress elements, and for whom.

    This is about the protein reference only. The nucleotide catalogues carry
    virulence-typed genes of their own and report them whatever this says, so
    "virulence off" never means "no virulence gene was reported" — it means the
    translated search was limited to acquired resistance.

    AMRFinderPlus's "plus" elements are one curated catalogue that is not split by
    organism: the same virulence genes are searched for whatever the isolate is.
    What the organism does change is curation — the suppression table is keyed by
    taxgroup — so a run with an established organism is curated for that organism
    and a run without one is not. Both are reported for what they are here, so a
    generic screen is never presented as an organism's curated virulence panel.
    """
    catalogue = catalogue if catalogue is not None else organism_catalogue(db_root)
    counts = counts if counts is not None else element_counts(db_root)
    resolved = match_organism(organism, accepted=catalogue["accepted"]) if organism else None
    if organism and not catalogue["accepted"]:
        resolved = str(organism)
    support = {"available": bool(counts.get("VIRULENCE") or counts.get("STRESS")),
               "virulence_elements": int(counts.get("VIRULENCE") or 0),
               "stress_elements": int(counts.get("STRESS") or 0),
               "organism": resolved or "", "organism_curated": bool(resolved),
               "suppression_organisms": list(catalogue["suppressed"]), "reason": ""}
    if not counts.get("read"):
        support["reason"] = ("The protein reference set is not installed, so no virulence or "
                             "stress element can be read from it.")
    elif not support["available"]:
        support["reason"] = ("The installed protein reference carries no virulence or stress "
                             "elements, so there is nothing for this option to add.")
    elif resolved:
        support["reason"] = (
            f"{support['virulence_elements']} virulence and {support['stress_elements']} stress "
            f"elements in the protein reference are searched, curated for '{resolved}'. This is "
            "NCBI's published reference data, not a validated virulence prediction, and a gene is "
            "not a demonstrated phenotype.")
    else:
        support["reason"] = (
            "No organism is established for this isolate, so the protein virulence and stress "
            "search would run with no organism curation applied: the same elements are reported, "
            "but nothing an organism's curation would have suppressed is suppressed. Assign a "
            "genus and species to have it curated.")
    return support


def _one(title):
    """The singular of a plain-words element title, for counting one of them."""
    return title[:-1] if title.endswith("s") and not title.endswith("ss") else title


def _set_note(name):
    """One reference set, named with what it answers and who publishes it."""
    provider = DATABASE_PROVIDER.get(name, "")
    return {"name": name, "provider": provider,
            "purpose": DATABASE_PURPOSE.get(name, "an additional reference set this engine reads")}


def element_type_support(db_root=None, databases=None, *, status=None):
    """Which element types the chosen reference sets can report, and which they cannot.

    The engine is silent about what it was never asked: it reports acquired
    resistance, virulence, stress and plasmid elements from whichever sets are
    installed and says nothing at all about the rest. A reader then sees no
    replicon and no virulence gene and has no way to tell that from a clean
    isolate. So every element type is listed, whether or not anything can report
    it, together with the sets that are being searched for it, the sets that are
    installed but not selected for this run, and the sets that would supply it
    and are simply not here. Nothing is downloaded or created by asking.
    """
    status = status if status is not None else database_status(db_root, databases)
    installed = set(status["installed"])
    chosen = [name for name in status["requested"] if name in installed]
    rows = {}
    for element_type in ELEMENT_TYPES:
        sources = ELEMENT_SOURCES[element_type]
        title = ELEMENT_TITLES[element_type]
        searched = [name for name in sources if name in chosen]
        idle = [name for name in sources if name in installed and name not in searched]
        absent = [_set_note(name) for name in sources if name not in installed]
        if searched:
            reason = f"{title.capitalize()} are reported from {', '.join(searched)}."
            if absent:
                reason += (" " + ", ".join(entry["name"] for entry in absent) + " would widen this "
                           "screen and " + ("is" if len(absent) == 1 else "are") + " not installed.")
        elif idle:
            reason = (f"No reference set that reports {title} is among the sets this run searches, "
                      f"although {', '.join(idle)} {'is' if len(idle) == 1 else 'are'} installed. "
                      "Nothing was looked for, so a blank is not a negative result.")
        else:
            reason = (f"No reference set that reports {title} is installed, so none can be reported "
                      f"for any isolate. {', '.join(entry['name'] for entry in absent)} would supply "
                      "them; install " + ("it" if len(absent) == 1 else "one")
                      + " from the AMR databases page. Nothing was looked for, so a blank is not a "
                        "negative result.")
        rows[element_type] = {"element_type": element_type, "title": title,
                              "searched": searched, "installed_not_selected": idle,
                              "not_installed": absent, "available": bool(searched),
                              "reason": reason}
    return rows


def element_evidence(evidence, *, execution=None, databases=None, installed=None):
    """Every element type one isolate's HYDRA screen could report, and what it found.

    This is the answer to "I ran Klebsiella and cannot find the virulence
    output". The engine returns acquired resistance, virulence, stress and
    plasmid-replicon hits in one list and the application keeps all of them, but
    an interface that reads only the AMR ones shows nothing for the rest, which
    is indistinguishable from an isolate that carries none.

    Four states, and they are deliberately not interchangeable: ``detected`` is
    a call, ``none_detected`` is a search that found nothing, ``not_searched``
    is a search that never ran because no reference set that reports this type
    was read, and ``no_report`` is an isolate with no HYDRA result at all.
    Point mutations are returned separately from the genes because an organism
    with no mutation catalogue is screened for genes only, and a list that
    merged the two would present that as an isolate with no mutations.

    ``evidence`` is a stored ``metadata['hydra']`` block or an upstream report
    sample. ``installed`` is optional and names what the reference store holds,
    so a set that is present but was not selected can be distinguished from one
    that is not there at all; without it an unsearched set is reported as
    unsearched and never as absent, because only the store can answer that.
    Nothing here reads a database, a file or the network.
    """
    evidence = evidence if isinstance(evidence, dict) else {}
    execution = execution if execution is not None else (evidence.get("execution_provenance") or {})
    execution = execution if isinstance(execution, dict) else {}
    if databases is None:
        databases = evidence.get("databases")
    if databases is None:
        snapshot = execution.get("reference_snapshot") or {}
        databases = sorted((snapshot.get("databases") or {}))
    searched = [str(name) for name in (databases or [])]
    # Whether the store was asked what it holds. Without that, a set that simply
    # was not read must not be reported as one the computer does not have.
    known = installed is not None
    present = {str(name) for name in installed} if known else set(searched)
    hits = [hit for hit in (evidence.get("hits") or []) if isinstance(hit, dict)]
    primary = [hit for hit in hits if hit.get("primary") is True]
    linked = bool(hits or searched or execution)
    virulence = execution.get("virulence") if isinstance(execution.get("virulence"), dict) else {}
    organism = execution.get("organism") if isinstance(execution.get("organism"), dict) else {}
    curated_for = str(virulence.get("organism") or organism.get("resolved") or "")
    elements = []
    for element_type in ELEMENT_TYPES:
        sources = ELEMENT_SOURCES[element_type]
        title = ELEMENT_TITLES[element_type]
        reading = [name for name in sources if name in searched]
        # --plus gates the translated search alone. With it off the protein set
        # reports acquired resistance only, so it contributes nothing to these
        # two types and must not be counted as having searched for them: a run
        # against protein alone would otherwise claim it looked and found none.
        if element_type in {"VIRULENCE", "STRESS"} and virulence and not virulence.get("enabled"):
            reading = [name for name in reading if name != "protein"]
        absent = [_set_note(name) for name in sources if name not in present]
        genes = sorted({str(hit["gene"]) for hit in primary
                        if hit.get("element_type") == element_type and hit.get("gene")})
        matches = [hit for hit in hits if hit.get("element_type") == element_type]
        caveats = []
        if not linked:
            status = "no_report"
            reason = (f"No HYDRA screen is linked to this isolate, so nothing is known about its "
                      f"{title}. That is unknown, not absent.")
        elif genes:
            # A hit is decisive even when the report did not record which set
            # found it: an imported report that lists no databases must not turn
            # its own findings into "nothing was searched for".
            status = "detected"
            reason = (f"{len(genes)} {title if len(genes) != 1 else _one(title)} matched "
                      + (", ".join(reading) if reading else "this screen's reference data")
                      + ". A reference match is sequence evidence, not a measured phenotype.")
        elif not reading:
            status = "not_searched"
            reason = (f"No reference set that reports {title} was searched for this isolate"
                      + (": " + ", ".join(entry["name"] for entry in absent) + " would supply them."
                         if absent else ".")
                      + " Nothing was looked for, so this is not a negative result.")
        else:
            status = "none_detected"
            reason = (f"{', '.join(reading)} {'was' if len(reading) == 1 else 'were'} searched and "
                      f"no {_one(title)} met this run's thresholds. A screen over public reference "
                      "data is not proof of absence.")
        # The protein search is the only one --plus gates, so switching virulence
        # off narrows what could be found without meaning nothing was sought. The
        # engine's own sentence is repeated verbatim rather than paraphrased.
        if element_type in {"VIRULENCE", "STRESS"} and virulence:
            if virulence.get("enabled"):
                caveats.append(f"The translated protein search was curated for '{curated_for}'."
                               if virulence.get("organism_curated") and curated_for else
                               "The translated protein search ran with no organism curation, so "
                               "nothing an organism's curation would have suppressed is suppressed.")
            elif reading:
                caveats.append("The translated protein search was limited to acquired resistance "
                               "for this run, so these come from the nucleotide catalogue only.")
            else:
                caveats.append("The translated protein search was limited to acquired resistance "
                               "for this run and no nucleotide catalogue that reports these was "
                               "read, so none could be reported at all.")
            if virulence.get("reason"):
                caveats.append(str(virulence["reason"]))
        if absent and reading:
            caveats.append(("Not installed, so not searched: " if known else "Not searched: ")
                           + "; ".join(f"{entry['name']} ({entry['purpose']})" for entry in absent) + ".")
        elements.append({"element_type": element_type, "title": title, "status": status,
                         "reason": reason, "genes": genes, "hits": matches,
                         "searched": reading, "not_installed": absent, "caveats": caveats})
    return {"linked": linked, "databases": searched, "elements": elements,
            "point_mutations": _point_mutation_evidence(primary, execution, searched, linked)}


def _point_mutation_evidence(primary, execution, searched, linked):
    """The catalogued resistance mutations, and whose catalogue could have named one.

    An organism outside the installed release's catalogues is screened for genes
    only. Reporting that as an isolate with no mutations is the failure this
    block exists to prevent, so the organism the run resolved and the engine's
    own sentence about it travel with the list.
    """
    organism = execution.get("organism") if isinstance(execution.get("organism"), dict) else {}
    matches = [hit for hit in primary if is_point_mutation(hit)]
    genes = sorted({str(hit["gene"]) for hit in matches if hit.get("gene")})
    level = str(organism.get("point_mutation_level") or "unknown")
    resolved = str(organism.get("resolved") or "")
    requested = organism.get("point_mutations", True) if organism else True
    if not linked:
        status = "no_report"
        reason = ("No HYDRA screen is linked to this isolate, so no resistance mutation was "
                  "assessed. That is unknown, not absent.")
    elif organism and not requested:
        status = "not_searched"
        reason = ("Point mutations were switched off for this run, so none was assessed. That is "
                  "not evidence that this isolate carries none.")
    elif "protein" not in searched:
        status = "not_searched"
        reason = ("The protein reference set was not searched, and it is where point mutations are "
                  "read, so none could be reported. That is not evidence that none is present.")
    elif organism and level in {"none", "unknown"}:
        status = "not_searched"
        reason = str(organism.get("reason") or
                     "No point-mutation catalogue was selected for this isolate, so none was "
                     "assessed. Assign a genus and species to have them assessed.")
    elif matches:
        status = "detected"
        reason = (f"{len(genes)} gene(s) carry a catalogued resistance mutation"
                  + (f", read against the '{resolved}' catalogue" if resolved else "")
                  + ". A catalogued mutation is a genomic association, not a susceptibility result.")
    else:
        status = "none_detected"
        reason = ("The catalogue"
                  + (f" for '{resolved}'" if resolved else "")
                  + " was searched and no catalogued resistance mutation was found. Mutations "
                    "outside that catalogue are not assessed by this screen.")
    return {"status": status, "reason": reason, "genes": genes, "hits": matches,
            "organism": resolved, "catalogue_level": level,
            "catalogue_reason": str(organism.get("reason") or "")}


ELEMENT_BOUNDARY = (
    'This is a BLAST screen over public reference data run by the pinned HYDRA engine. It is not '
    'AMRFinderPlus, Kleborate, Kaptive or MOB-suite and is not equivalent to them. A gene is sequence '
    'evidence, never a measured susceptibility, a virulence phenotype or a plasmid.')
HIT_LIMIT = 40


def element_evidence_html(evidence, *, execution=None, databases=None, installed=None,
                          limit=HIT_LIMIT):
    """Every element type of one isolate's screen, blanks explained, as a drill-down.

    Rendered here beside the evidence it describes, as the organism-specific
    assays render theirs, so a view needs one call to show all four element
    types instead of reimplementing the distinction between a negative and an
    unasked question. Every value is escaped: none of this text is trusted.
    """
    def escape(value):
        return html.escape(str(value if value is not None else '—'))

    result = element_evidence(evidence, execution=execution, databases=databases,
                              installed=installed)
    parts = ['<h3>HYDRA screen · resistance, virulence, stress and plasmid elements</h3>']
    if not result['linked']:
        return ''.join(parts + [
            '<p>No HYDRA screen is linked to this isolate, so nothing is known about its '
            'resistance, virulence, stress or plasmid elements. That is unknown, not absent.</p>',
            f'<p>{escape(ELEMENT_BOUNDARY)}</p>'])
    parts.append('<p>Searched: ' + escape(', '.join(result['databases']) or 'no reference set')
                 + '. Every element type the engine can report is listed, including the ones '
                   'nothing searched for.</p>')
    for row in result['elements']:
        parts.append(f'<h4>{escape(row["title"].capitalize())} · {escape(row["status"].replace("_", " "))}</h4>')
        parts.append(f'<p>{escape(row["reason"])}</p>')
        if row['genes']:
            parts.append('<p><b>Detected:</b> ' + escape('; '.join(row['genes'])) + '</p>')
            parts.append(_hit_table(row['hits'], limit, escape))
        for note in row['caveats']:
            parts.append(f'<p class="muted">{escape(note)}</p>')
    mutations = result['point_mutations']
    parts.append('<h4>Resistance point mutations · '
                 + escape(mutations['status'].replace('_', ' ')) + '</h4>')
    parts.append(f'<p>{escape(mutations["reason"])}</p>')
    if mutations['genes']:
        parts.append('<p><b>Mutated genes:</b> ' + escape('; '.join(mutations['genes'])) + '</p>')
        parts.append(_hit_table(mutations['hits'], limit, escape))
    parts.append(f'<p>{escape(ELEMENT_BOUNDARY)}</p>')
    return ''.join(parts)


def _hit_table(hits, limit, escape):
    """The matches behind a call, highest identity first, with the total always stated."""
    ordered = sorted(hits, key=lambda hit: (-float(hit.get('identity_pct') or 0),
                                            -float(hit.get('coverage_pct') or 0),
                                            str(hit.get('gene') or '')))
    body = ("<table cellpadding='5'><tr><th>Gene</th><th>Reference set</th><th>Class</th>"
            "<th>Identity %</th><th>Reference covered %</th><th>Method</th><th>Resolution</th></tr>")
    for hit in ordered[:limit]:
        body += ('<tr>' + ''.join(f'<td>{escape(hit.get(key))}</td>' for key in (
            'gene', 'database', 'class', 'identity_pct', 'coverage_pct', 'method', 'resolution'))
            + '</tr>')
    body += '</table>'
    if len(ordered) > limit:
        body += f'<p>Showing the {limit} highest-identity matches of {len(ordered)}.</p>'
    return body


def catalogue_covers(organism, names):
    """The engine's own taxgroup rule: a catalogue row applies to a parent or a child.

    Mirrors hydra_amr.engines.mutations.MutationCatalog._taxgroup_matches, where
    "Escherichia" covers "Escherichia_coli" and "Campylobacter" covers every
    campylobacter. Testing exact membership instead would tell a user that an
    isolate has no catalogue when the engine is about to apply one — a false
    statement about coverage in the direction that looks safe and is not.
    """
    org = str(organism or "")
    if not org:
        return False
    return any(org == name or org.startswith(name + "_") or name.startswith(org + "_")
               for name in names)


def match_organism(name, db_root=None, *, accepted=None):
    """Map an assigned organism onto a name the installed catalogues really use.

    Upstream's taxgroups mix species ("Klebsiella_pneumoniae") with genus-level
    groups ("Escherichia", "Salmonella"), and they are underscore-joined, so an
    assignment of "Escherichia coli" belongs to "Escherichia" and one of
    "Klebsiella pneumoniae" to "Klebsiella_pneumoniae". Nothing is resolved by
    similarity: a name with no catalogue returns None so the caller can say that
    plainly instead of substituting a neighbouring organism's mutation list.
    """
    text = " ".join(str(name or "").replace("_", " ").split())
    if not text:
        return None
    choices = list(accepted) if accepted is not None else organism_catalogue(db_root)["accepted"]
    index = {choice.replace("_", " ").casefold(): choice for choice in choices}
    return index.get(text.casefold()) or index.get(text.split(" ")[0].casefold())


def database_status(db_root=None, databases=None):
    """What a store holds, which release it is, how old that release is, what is absent.

    Plain data only: nothing is downloaded, created or repaired by asking. The
    release and its age are reported separately from the day the copy was staged,
    because it is the reference release — not the download — that decides which
    determinants the run could possibly name.
    """
    root = Path(db_root).resolve() if db_root is not None else None
    wanted = [str(name) for name in (databases if databases is not None else DEFAULT_DATABASES)]
    status = {"root": str(root or ""), "bundled": False, "installed": {}, "requested": wanted,
              "missing": [], "error": "", "release": "", "released_utc": "", "staged": "",
              "age_days": None, "stale": False, "label": "No AMR reference database is installed.",
              "organisms": [], "point_mutation_organisms": [],
              "dna_point_mutation_organisms": [], "protein_point_mutation_organisms": [],
              "suppression_organisms": []}
    bundled = bundled_database_root()
    status["bundled"] = bool(root is not None and bundled is not None
                             and root == Path(bundled).resolve())
    if root is None:
        return status
    try:
        entries = installed_databases(root)
    except HydraRuntimeError as exc:
        status["error"] = str(exc)
        status["label"] = str(exc)
        return status
    for name, entry in sorted(entries.items()):
        status["installed"][name] = {
            "path": str(entry.get("path", "")), "kind": str(entry.get("kind", "")),
            "version": str(entry.get("version", "") or ""),
            "installed": str(entry.get("installed", "") or ""),
            "sequences": entry.get("sequences"),
            "title": str(entry.get("title", "") or ""),
            "purpose": DATABASE_PURPOSE.get(name, "an additional reference set this engine reads")}
    status["missing"] = [{"name": name, "purpose": DATABASE_PURPOSE.get(
        name, "an additional reference set this engine reads")} for name in wanted
        if name not in entries]
    if not entries:
        return status
    catalogue = organism_catalogue(root)
    status["organisms"] = catalogue["accepted"]
    status["point_mutation_organisms"] = catalogue["point_mutations"]
    status["dna_point_mutation_organisms"] = catalogue["dna_point_mutations"]
    status["protein_point_mutation_organisms"] = catalogue["protein_point_mutations"]
    status["suppression_organisms"] = catalogue["suppressed"]
    versions = sorted({item["version"] for item in status["installed"].values() if item["version"]})
    status["release"] = versions[0] if len(versions) == 1 else ("; ".join(versions) or "unrecorded")
    staged = sorted({item["installed"] for item in status["installed"].values()
                     if item["installed"]})
    status["staged"] = staged[0] if staged else ""
    released = _release_date(versions[0]) if len(versions) == 1 else None
    if released is not None:
        status["released_utc"] = released.isoformat()
        status["age_days"] = max(0, (datetime.now(UTC) - released).days)
        status["stale"] = status["age_days"] > REFERENCE_AGE_DAYS
    where = "bundled with this release" if status["bundled"] else "installed"
    age = (f", {status['age_days']} days old" if status["age_days"] is not None else "")
    status["label"] = (f"NCBI AMRFinderPlus reference release {status['release']} "
                       f"({', '.join(sorted(entries))}), {where}"
                       + (f" on {status['staged'][:10]}" if status["staged"] else "") + age + ".")
    return status


def _registry_specs():
    """The engine's own database registry, or an explanation of why it is absent."""
    try:
        from hydra_amr.db.fetch import can_fetch
        from hydra_amr.db.registry import DATABASES
    except ImportError as exc:
        return {}, None, (f"The HYDRA engine is not installed in this environment ({exc}), so "
                          "only the reference sets already in the store can be listed.")
    return dict(DATABASES), can_fetch, ""


def _database_bytes(root, entry, name):
    """Measured size of one installed set, the point-mutation catalogues included."""
    directories = [_entry_directory(root, entry)]
    if name == "protein" and (Path(root) / "mutation").is_dir():
        directories.append(Path(root) / "mutation")
    total = 0
    try:
        for directory in directories:
            for path in directory.rglob("*"):
                if path.is_file():
                    total += path.stat().st_size
    except OSError:
        return None
    return total


def _human_size(count):
    if not count:
        return ""
    megabytes = count / (1024 * 1024)
    return f"{megabytes:.0f} MB" if megabytes >= 1 else f"{count / 1024:.0f} KB"


def _licence_is_open(licence):
    text = str(licence or "").casefold()
    return any(token in text for token in OPEN_LICENCE_TOKENS)


def database_catalogue(db_root=None, *, measure=True):
    """Every reference set HYDRA can search, by name, and whether this computer has it.

    This exists because the engine is powerful and silent about its inputs: it
    will happily run against whatever is installed and report nothing for
    everything that is not, and until a person can see the whole list there is no
    way to tell a clean isolate from an unasked question. So every set is listed
    whether or not it is installed, with what it is for, who publishes it, under
    what licence, how large it is and which release is here.

    The names, titles, addresses, citations and licences come from the engine's
    own registry rather than a second list that could drift away from it. The
    plain-language purpose, the provider's name and the size are ours, because the
    registry does not carry them. Nothing is downloaded and no server is
    contacted: a row that says "not installed" is a statement about this computer.
    """
    root = Path(db_root).resolve() if db_root is not None else None
    error = ""
    try:
        entries = installed_databases(root) if root is not None else {}
    except HydraRuntimeError as exc:
        entries, error = {}, str(exc)
    specs, can_fetch, registry_error = _registry_specs()
    error = error or registry_error
    bundled = bundled_database_root()
    is_bundled = bool(root is not None and bundled is not None
                      and root == Path(bundled).resolve())
    rows = []
    for name in sorted(set(specs) | set(entries) | set(DEFAULT_DATABASES)):
        spec, entry = specs.get(name), entries.get(name)
        licence = getattr(spec, "licence", "") or ""
        size, basis = DATABASE_SIZE.get(name, (None, ""))
        measured = _database_bytes(root, entry, name) if (entry and measure and root) else None
        if measured:
            readable, basis = _human_size(measured), "measured in the installed store"
        elif size:
            readable = f"about {size} MB"
        else:
            readable = ""
            basis = ("no size is published for this set; it is measured once it is installed")
        rows.append({
            "name": name,
            "title": getattr(spec, "title", "") or (entry or {}).get("title") or name,
            "kind": getattr(spec, "kind", "") or (entry or {}).get("kind", ""),
            "element_type": getattr(spec, "element_type", "") or (entry or {}).get("element_type", ""),
            "purpose": DATABASE_PURPOSE.get(name, "an additional reference set this engine reads"),
            "provider": DATABASE_PROVIDER.get(name, ""),
            "url": getattr(spec, "url", ""), "citation": getattr(spec, "citation", ""),
            "notes": getattr(spec, "notes", ""),
            "licence": licence or "not recorded by the provider",
            "open_licence": _licence_is_open(licence),
            "licence_note": ("" if _licence_is_open(licence) else
                             "The provider's own terms apply; read them before installing, using "
                             "or sharing this data."),
            "download": ("automatic" if (can_fetch and can_fetch(name)) else "by hand"),
            "core": name in DEFAULT_DATABASES,
            "installed": entry is not None,
            "bundled": is_bundled and entry is not None,
            "version": str((entry or {}).get("version", "") or ""),
            "staged": str((entry or {}).get("installed", "") or ""),
            "sequences": (entry or {}).get("sequences"),
            "bytes": measured, "size": readable, "size_basis": basis,
            "state": "installed" if entry is not None else "not installed",
        })
    # Core sets first, then whatever else is already here, then the rest by name:
    # the order somebody reads the list in, not the order a dict happened to hold.
    rows.sort(key=lambda row: (not row["core"], not row["installed"], row["name"]))
    here = [row["name"] for row in rows if row["installed"]]
    elsewhere = [row["name"] for row in rows if not row["installed"]]
    summary = (f"{len(here)} of {len(rows)} reference sets are installed in "
               f"{root or 'no store'}" + (": " + ", ".join(here) if here else "")
               + f". The other {len(elsewhere)} are listed here and downloaded only when you ask "
                 "for one; an isolate is never screened against a set that is not installed, and "
                 "a result does not say so unless you read this page.")
    return {"root": str(root or ""), "error": error, "entries": rows, "installed": here,
            "available": elsewhere, "core": list(DEFAULT_DATABASES), "summary": summary,
            "automatic": [row["name"] for row in rows if row["download"] == "automatic"]}


def runtime_capabilities(db_root=None):
    try:
        distribution = importlib.metadata.distribution("hydra-amr")
        version = distribution.version
        origin = json.loads(distribution.read_text("direct_url.json") or "{}")
    except importlib.metadata.PackageNotFoundError:
        version = None
        origin = {}
    tools = {name: resolve_tool(name) for name in TOOLS}
    databases = installed_databases(db_root) if db_root is not None else {}
    available = version == HYDRA_VERSION and all(tools.values())
    missing = [name for name, path in tools.items() if not path]
    if not available:
        message = f"HYDRA {version or 'not installed'}; missing tools: {', '.join(missing) or 'none'}."
    elif db_root is None or not databases:
        # The engine being ready is not the same as the screen being runnable, and
        # reporting only the first was how "HYDRA does not work" stayed unexplained.
        message = ("The HYDRA engine and BLAST+ tools are ready, but no AMR reference database is "
                   "installed, so a run would have nothing to search.")
    else:
        message = ("Native HYDRA assembly runtime is ready with "
                   + ", ".join(sorted(databases)) + ".")
    return {"hydra_version": version, "expected_hydra_version": HYDRA_VERSION,
            "source_commit": HYDRA_COMMIT, "distribution_origin": origin,
            "tools": tools, "databases": databases,
            "assembly_available": available, "available": available, "message": message,
            "blastn": tools["blastn"], "makeblastdb": tools["makeblastdb"],
            "reads_available": False,
            "limitations": ["Assembly-only HYDRA execution; its direct-read and minority-allele pipelines are not enabled. Paired-read assembly is handled separately by native SKESA.",
                            "Gene and mutation evidence is not a validated susceptibility phenotype."]}


def preflight(db_root=None, *, databases=None, organism=None, point_mutations=True,
              virulence=None, capabilities=None):
    """Everything a run needs, checked before a single sequence file is opened.

    The failure this exists to stop is the quiet one: HYDRA is installed, BLAST is
    installed, no reference store is, and the run either refuses with a sentence
    about selecting a database or — worse — completes and reports nothing, which a
    reader cannot tell from a clean isolate. So every piece is checked first, each
    missing one is named together with what it is for and what the run loses
    without it, and `ready` is False before any work starts.

    Nothing here downloads, creates or repairs anything, and `warnings` is kept
    separate from `missing`: a warning narrows what the run can report and is
    recorded with the result, a missing item stops it.

    `virulence` is None for "use it where the organism makes it meaningful", True
    to search for virulence and stress elements whatever the organism is, and
    False to leave them alone. Whichever is chosen, the returned `virulence` block
    says what was searched for and under whose curation, because an uncurated
    plus-element screen and an organism's curated one are not the same result.
    """
    capabilities = capabilities if capabilities is not None else runtime_capabilities(db_root)
    tools = capabilities.get("tools") or {}
    version = capabilities.get("hydra_version")
    missing, warnings = [], []
    if version != HYDRA_VERSION:
        missing.append({"kind": "engine", "name": f"HYDRA {HYDRA_VERSION}",
                        "purpose": "the analysis engine itself",
                        "reason": (f"{version or 'No HYDRA engine'} is installed and this "
                                   f"application is pinned to {HYDRA_VERSION}; it will not report "
                                   "results from an engine it was not validated against.")})
    for name in TOOLS:
        if not tools.get(name):
            missing.append({"kind": "tool", "name": name, "purpose": TOOL_PURPOSE[name],
                            "reason": (f"{name} was not found beside the application or on PATH, "
                                       "so this search cannot be run at all.")})
    status = database_status(db_root, databases)
    if db_root is None:
        missing.append({"kind": "database", "name": "AMR reference database",
                        "purpose": "every gene and mutation this screen can name",
                        "reason": "No database store has been chosen for this project."})
    elif status["error"]:
        missing.append({"kind": "database", "name": "AMR reference database",
                        "purpose": "every gene and mutation this screen can name",
                        "reason": status["error"]})
    chosen = [name for name in status["requested"] if name in status["installed"]]
    if db_root is not None and not status["error"]:
        if not status["installed"]:
            missing.append({"kind": "database", "name": "AMR reference database",
                            "purpose": "every gene and mutation this screen can name",
                            "reason": (f"{status['root']} holds no reference database, so the "
                                       "screen has nothing to search and would report no "
                                       "determinants for every isolate.")})
        elif databases is not None:
            # A database the caller asked for by name is not optional: running
            # without it would answer a different question than the one asked.
            for entry in status["missing"]:
                missing.append({"kind": "database", **entry,
                                "reason": (f"The database '{entry['name']}' is not installed in "
                                           f"{status['root']}.")})
        else:
            for entry in status["missing"]:
                warnings.append(
                    f"The '{entry['name']}' reference set is not installed, so {entry['purpose']} "
                    "did not run. Nothing was reported for it; that is not a negative result.")
    resolved, organism_reason, level = None, "", "unknown"
    if organism:
        accepted = status["organisms"]
        resolved = match_organism(organism, accepted=accepted) if accepted else str(organism)
        dna = catalogue_covers(resolved, status["dna_point_mutation_organisms"])
        protein_level = catalogue_covers(resolved, status["protein_point_mutation_organisms"])
        if accepted and resolved is None:
            level = "none"
            organism_reason = (
                f"The installed reference release has no catalogue for '{organism}', so point "
                "mutations were not assessed for this isolate. Genes were still searched for. "
                "An absent catalogue is not an absence of mutations.")
            warnings.append(organism_reason)
        elif resolved and accepted and dna and protein_level:
            level = "dna_and_protein"
        elif resolved and accepted and protein_level:
            level = "protein_only"
            organism_reason = (
                f"'{resolved}' has a curated protein mutation catalogue in this release but no DNA "
                "catalogue, so mutations in genes that are read at the DNA level (23S rRNA and the "
                "other non-coding targets) could not be assessed for it.")
            warnings.append(organism_reason)
        elif resolved and accepted and dna:
            level = "dna_only"
            organism_reason = (
                f"'{resolved}' has a DNA mutation catalogue in this release but no curated protein "
                "mutation entries, so only DNA-level mutations could be reported for it.")
            warnings.append(organism_reason)
        elif resolved and accepted:
            level = "none"
            organism_reason = (
                f"'{resolved}' is accepted by the installed release but has no point-mutation "
                "catalogue in it at all, so this isolate was screened for genes only. That is not "
                "evidence that it carries no resistance mutation.")
            warnings.append(organism_reason)
    elif point_mutations:
        level = "none"
        organism_reason = ("No organism was given, so no point-mutation catalogue was selected. "
                           "Assign a genus and species to this isolate to have them assessed.")
        warnings.append(organism_reason)
    if point_mutations and "protein" not in chosen and status["installed"]:
        level = "none"
        warnings.append("Point mutations were requested but the protein reference set is not "
                        "among the databases being searched, so none can be reported.")
    # Virulence and stress elements from the PROTEIN reference: offered where the
    # isolate's organism is established, because the engine's curation is keyed to
    # that organism and a run without one is a different, uncurated search. "auto"
    # is the caller saying "use it where it applies" and is the only setting that
    # decides by itself; True and False are the user's own choice and are obeyed.
    # The nucleotide catalogues report their own virulence-typed genes either way,
    # so none of these sentences may say that nothing virulent was looked for.
    support = {"available": False, "virulence_elements": 0, "stress_elements": 0,
               "organism": resolved or "", "organism_curated": bool(resolved),
               "suppression_organisms": [], "reason": ""}
    if status["installed"] and "protein" in chosen:
        try:
            # database_status has already read the organism tables; handing them
            # over keeps one preflight to one read of each reference table.
            support = virulence_support(
                status["root"], organism=resolved or organism,
                catalogue={"accepted": status["organisms"],
                           "suppressed": status["suppression_organisms"]})
        except HydraRuntimeError as exc:
            support["reason"] = str(exc)
    elif status["installed"]:
        support["reason"] = ("The protein reference set is not among the databases being searched, "
                             "so no virulence or stress element can be read from it.")
    requested = "auto" if virulence is None else bool(virulence)
    virulence_reason = support["reason"]
    if requested is False:
        enabled = False
        virulence_reason = ("The translated protein search was limited to acquired resistance, so "
                            "no virulence or stress element was read from the protein reference. "
                            "Virulence-typed genes in the nucleotide catalogues are reported "
                            "regardless: this setting governs the protein search only.")
    elif not support["available"]:
        enabled = False
    elif requested is True:
        enabled = True
    else:
        enabled = bool(support["organism_curated"])
    if virulence_reason and requested is not False and (not enabled
                                                        or not support["organism_curated"]):
        warnings.append(virulence_reason)
    if status["stale"]:
        warnings.append(f"This reference release is {status['age_days']} days old "
                        f"({status['release']}). Determinants named after it are not in it; use "
                        "the update action to fetch the current NCBI release.")
    message = ""
    if missing:
        lines = [f"  • {item['name']} — needed for {item['purpose']}. {item['reason']}"
                 for item in missing]
        fix = ("Install the BLAST+ tools that ship beside the application"
               if any(item["kind"] == "tool" for item in missing) else
               "Download the reference databases from HYDRA > AMR databases / updates")
        message = ("HYDRA cannot start, so nothing was run and no isolate was marked screened.\n"
                   + "\n".join(lines) + f"\n{fix}. Nothing is ever downloaded during a run.")
    return {"ready": not missing, "message": message, "missing": missing, "warnings": warnings,
            "databases": chosen, "database_status": status,
            "organism": {"requested": str(organism or ""), "resolved": resolved or "",
                         "reason": organism_reason,
                         "point_mutations": bool(point_mutations),
                         "point_mutation_level": level},
            "virulence": {"requested": requested, "enabled": enabled,
                          "reason": virulence_reason,
                          "organism_curated": bool(support["organism_curated"]),
                          "available": bool(support["available"]),
                          "virulence_elements": support["virulence_elements"],
                          "stress_elements": support["stress_elements"]},
            "engine": {"installed": version or "", "expected": HYDRA_VERSION},
            "tools": dict(tools), "limitations": list(capabilities.get("limitations") or ())}


def _worker_command(arguments):
    if getattr(sys, "frozen", False):
        suffix = ".exe" if os.name == "nt" else ""
        worker = Path(sys.executable).resolve().parent / ("WMLSTudio-HYDRA" + suffix)
        if not worker.is_file():
            raise HydraRuntimeError("The portable HYDRA worker is missing. Extract the complete WMLSTudio folder.")
        return [str(worker), "--upstream", *arguments]
    return [sys.executable, "-m", "wmlstudio.hydra_runtime", "--upstream", *arguments]


def _child_environment():
    environment = dict(os.environ)
    directories = [str(directory) for directory in tool_directories() if directory.is_dir()]
    environment["PATH"] = os.pathsep.join(directories + [environment.get("PATH", "")])
    environment["PYTHONUTF8"] = "1"
    environment["HYDRA_NO_BANNER"] = "1"
    # BLAST receives its explicit per-sample --threads allocation. Prevent
    # numerical helper libraries from creating an additional machine-wide pool
    # in every concurrently running HYDRA worker.
    for key in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
        environment[key] = "1"
    # Frozen Linux applications prepend their own shared libraries. BLAST must
    # resolve its native dependencies from its own runtime, not Qt's libraries.
    if "LD_LIBRARY_PATH_ORIG" in environment:
        environment["LD_LIBRARY_PATH"] = environment["LD_LIBRARY_PATH_ORIG"]
    elif getattr(sys, "frozen", False):
        environment.pop("LD_LIBRARY_PATH", None)
    return environment


def _stop_process(process):
    if process.poll() is not None:
        return
    if os.name == "nt":
        # PID belongs to the child we launched, and /T includes its BLAST children.
        subprocess.run(["taskkill.exe", "/PID", str(process.pid), "/T", "/F"],
                       capture_output=True, timeout=15, check=False,
                       creationflags=subprocess.CREATE_NO_WINDOW)
    else:
        try:
            os.killpg(process.pid, signal.SIGTERM)
        except ProcessLookupError:
            pass
    try:
        process.wait(timeout=5)
    except subprocess.TimeoutExpired:
        if os.name != "nt":
            try:
                os.killpg(process.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
        process.kill()
        process.wait(timeout=5)


def _run_child(arguments, directory, cancelled=None, progress=None):
    if cancelled and cancelled():
        raise AnalysisCancelled("HYDRA cancelled before execution.")
    command = _worker_command(arguments)
    log_path = Path(directory) / "hydra-execution.log"
    settings = ({"creationflags": subprocess.CREATE_NO_WINDOW | subprocess.CREATE_NEW_PROCESS_GROUP}
                if os.name == "nt" else {"start_new_session": True})
    previous_size = 0
    with log_path.open("wb") as log:
        process = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT,
                                   cwd=directory, env=_child_environment(), **settings)
        try:
            while process.poll() is None:
                if cancelled and cancelled():
                    raise AnalysisCancelled("HYDRA cancelled; partial results were not imported.")
                size = log_path.stat().st_size
                if progress and size != previous_size:
                    with log_path.open("rb") as reader:
                        reader.seek(max(0, size - 2048))
                        tail = reader.read().decode("utf-8", errors="replace").splitlines()
                    if tail:
                        progress(0, 0, tail[-1][-300:])
                    previous_size = size
                time.sleep(0.1)
        except BaseException:
            _stop_process(process)
            raise
    with log_path.open("rb") as log:
        log.seek(max(0, log_path.stat().st_size - 128_000))
        log_text = log.read().decode("utf-8", errors="replace")
    if process.returncode:
        raise HydraRuntimeError(f"HYDRA exited with code {process.returncode}.\n{log_text[-5000:]}")
    return command, log_text


def _database_provenance(db_root, names, cancelled=None):
    root = Path(db_root).resolve()
    manifest = installed_databases(root)
    entries = {}
    for name in names:
        if name not in manifest:
            raise HydraRuntimeError(f"HYDRA database '{name}' is not installed. Download or choose a populated database store first.")
        entry = manifest[name]
        directory = _entry_directory(root, entry)
        files = {}
        references = list(directory.rglob("*"))
        if name == "protein" and (root / "mutation").is_dir():
            references.extend((root / "mutation").rglob("*"))
        for path in sorted(set(references)):
            if not path.is_file():
                continue
            if not path.resolve().is_relative_to(root):
                raise HydraRuntimeError(f"Database contains an external file link: {path}")
            # Index files are included: the exact reference/index snapshot is auditable.
            files[path.relative_to(root).as_posix()] = {"size": path.stat().st_size,
                                                       "sha256": file_sha256(path, cancelled)}
        entries[name] = {"manifest_entry": entry, "files": files}
    return {"root": str(root), "manifest_sha256": file_sha256(root / "manifest.json", cancelled),
            "databases": entries}


def run_assemblies(inputs, db_root, databases=None, *, sample_names=None, organism=None,
                   threads=2, protein=True, point_mutations=True, virulence=None, cancelled=None,
                   progress=None, work_root=None, min_identity=80, min_coverage=60,
                   protein_min_identity=90, protein_min_coverage=90):
    """Run the original HYDRA assembly pipeline and return its validated JSON report.

    Every run goes through `preflight` first, so a run that could not finish never
    starts and the caller gets a refusal naming each missing tool or reference set.
    What the run could not look for — an absent reference set, an organism with no
    catalogue, the age of the release — is recorded in the report's provenance and
    repeated in its import warnings, because an empty result and an unperformed
    search read identically otherwise.

    Input paths and reference files are never rewritten by this wrapper. The
    database store should be a managed writable cache because upstream may build
    missing indexes there. The caller must not concurrently update that store.
    """
    capabilities = runtime_capabilities(db_root)
    # One gate for every caller: a run that cannot finish never starts, and the
    # refusal names each missing tool or reference set and what it was for.
    checked = preflight(db_root, databases=databases, organism=organism,
                        point_mutations=point_mutations, virulence=virulence,
                        capabilities=capabilities)
    if not checked["ready"]:
        raise HydraRuntimeError(checked["message"])
    if not isinstance(threads, int) or isinstance(threads, bool) or not 1 <= threads <= 256:
        raise HydraRuntimeError("Threads must be an integer between 1 and 256.")
    thresholds = {"min-identity": min_identity, "min-coverage": min_coverage,
                  "protein-min-identity": protein_min_identity,
                  "protein-min-coverage": protein_min_coverage}
    for name, value in thresholds.items():
        if (isinstance(value, bool) or not isinstance(value, (int, float))
                or not math.isfinite(value) or not 0 <= value <= 100):
            raise HydraRuntimeError(f"{name} must be a finite percentage between 0 and 100.")
    paths = [Path(path).resolve() for path in inputs]
    if not paths or len(paths) != len(set(paths)):
        raise HydraRuntimeError("Choose one or more distinct assembly files.")
    names = list(sample_names) if sample_names is not None else None
    if names is not None and (len(names) != len(paths)
                              or any(not isinstance(name, str) or not name.strip() or any(c in name for c in "\n\r\t") for name in names)
                              or len(set(names)) != len(names)):
        raise HydraRuntimeError("Sample names must be unique nonempty labels, one per assembly.")
    signatures = {}
    input_provenance = []
    for index, path in enumerate(paths):
        if not path.is_file():
            raise HydraRuntimeError(f"Assembly was not found: {path}")
        with SequenceReader(path, cancelled) as reader:
            if reader.kind != "fasta":
                raise HydraRuntimeError(f"HYDRA execution accepts FASTA assemblies only: {path.name}")
        signatures[path] = file_signature(path)
        input_provenance.append({"path": str(path), "sha256": file_sha256(path, cancelled),
                                 "sample_name": names[index] if names is not None else sample_name(path)})
    names_db = list(databases) if databases is not None else list(checked["databases"])
    if not names_db or any(not isinstance(name, str) or not name or name.startswith("-") for name in names_db):
        raise HydraRuntimeError("Select at least one installed HYDRA database.")
    if progress:
        progress(0, 0, "Recording the selected database snapshot…")
    reference = _database_provenance(db_root, names_db, cancelled)
    arguments = ["run", "--db-dir", str(Path(db_root).resolve()), "--format", "json", "--prefix", "hydra",
                 "--threads", str(threads), "--no-banner", "--no-mlst", "--no-typing",
                 "--no-heteroresistance", "--no-reads-mlst", "--no-reads-variants"]
    for name, value in thresholds.items():
        arguments.extend(["--" + name, str(value)])
    # WMLSTudio's native typing evidence is not silently replaced by another engine.
    for name in names_db:
        arguments.extend(["--db", name])
    for path in paths:
        arguments.extend(["--assembly", str(path)])
    if names is not None:
        for name in names:
            arguments.extend(["--name", name])
    # The engine takes underscore-joined taxgroup names and stops the whole run on
    # one it does not know, so the assigned organism is resolved against the
    # installed catalogue first. A name with no catalogue is not passed and not
    # silently replaced by a neighbouring one: the run proceeds without mutation
    # catalogs and says so in the report.
    selected_organism = checked["organism"]["resolved"]
    if selected_organism:
        arguments.extend(["--organism", selected_organism])
    else:
        # Without native species sketches and an explicit organism, mutation
        # catalogs must not be chosen by an unvalidated inference.
        arguments.append("--no-auto-organism")
    if not protein:
        arguments.append("--no-protein")
    if not point_mutations:
        arguments.append("--no-point-mutations")
    # Stress and virulence elements are reported only when preflight decided they
    # apply to this isolate; the flag is always passed explicitly so the recorded
    # command says which of the two searches was performed.
    arguments.append("--plus" if checked["virulence"]["enabled"] else "--no-plus")
    if work_root is not None:
        Path(work_root).mkdir(parents=True, exist_ok=True)
    started = time.monotonic()
    with tempfile.TemporaryDirectory(prefix="wmlstudio-hydra-", dir=work_root) as work:
        arguments.extend(["--outdir", str(Path(work) / "results"), "--tmpdir", str(Path(work) / "temporary")])
        if progress:
            progress(0, 0, "Running pinned HYDRA assembly searches with native BLAST+…")
        command, log = _run_child(arguments, work, cancelled, progress)
        for path, signature in signatures.items():
            if file_signature(path) != signature:
                raise HydraRuntimeError(f"Assembly changed during HYDRA analysis: {path.name}. Run it again.")
        if file_sha256(Path(db_root) / "manifest.json", cancelled) != reference["manifest_sha256"]:
            raise HydraRuntimeError("The HYDRA database manifest changed during analysis. Do not update databases while a run is active.")
        report = load_hydra_report(Path(work) / "results/hydra.json")
    report["execution_provenance"] = {
        "engine": "HYDRA", "version": HYDRA_VERSION, "source_commit": HYDRA_COMMIT,
        "mode": "native-assembly", "command": command, "inputs": input_provenance,
        "distribution_origin": capabilities.get("distribution_origin", {}),
        "reference_snapshot": reference, "tools": capabilities["tools"],
        "completed_utc": datetime.now(UTC).isoformat(), "runtime_seconds": round(time.monotonic() - started, 3),
        "log": log, "scientific_limits": capabilities["limitations"],
        # What the run could and could not look for travels with the result, so a
        # report is never read as a complete screen when part of it never ran.
        "organism": dict(checked["organism"]),
        "virulence": dict(checked["virulence"]),
        "reference_release": {key: checked["database_status"][key]
                              for key in ("release", "released_utc", "staged", "age_days",
                                          "stale", "bundled", "label")},
        "coverage_warnings": list(checked["warnings"]),
    }
    if checked["warnings"]:
        report.setdefault("import_warnings", []).extend(checked["warnings"])
    return report


def _ncbi_version():
    request = urllib.request.Request(NCBI_REFERENCE_ROOT + "/version.txt",
                                     headers={"User-Agent": "WMLSTudio-reference-provenance"})
    with urllib.request.urlopen(request, timeout=30) as response:
        value = response.read(256).decode("ascii").strip()
    if not re.fullmatch(r"[0-9]{4}-[0-9]{2}-[0-9]{2}\.[0-9]+", value):
        raise HydraRuntimeError("NCBI returned an unrecognized reference version; snapshot was not published.")
    return value


def latest_release(*, cancelled=None):
    """Ask NCBI which reference release is current today, without downloading it.

    An explicit, user-initiated network call: nothing in startup, project open or
    analysis reaches here, and no sample data is sent — only a request for one
    version string. It exists so an update menu can say whether there is anything
    to fetch before asking a person to wait for 25 MB.
    """
    if cancelled and cancelled():
        raise AnalysisCancelled("Cancelled before the version check.")
    return _ncbi_version()


def update_databases(db_root, names, *, cancelled=None, progress=None, work_root=None):
    """Explicit user-initiated provider downloads; never called during analysis.

    Upstream owns normalization, version discovery and index construction. Its
    manifest records the selected providers. A selected-database-only snapshot is
    staged beside the requested store and atomically published under a new name.
    Existing stores and their analysis provenance remain untouched.
    """
    selected = list(names)
    if not selected or any(not isinstance(name, str) or not name or name.startswith("-") for name in selected):
        raise HydraRuntimeError("Select explicit database names for download/update.")
    # A name the engine does not know is refused here, by name, against the list a
    # person can read; without this the typo surfaces minutes later as an engine
    # error in the middle of a download, with the previous store already replaced
    # in the user's mind if not on disk.
    known = {row["name"] for row in database_catalogue(measure=False)["entries"]}
    unknown = sorted(set(selected) - known) if known else []
    if unknown:
        raise HydraRuntimeError(
            f"This engine has no reference set called {', '.join(unknown)}. The sets it can "
            f"install are: {', '.join(sorted(known))}.")
    # Nucleotide normalization consumes protein family annotations when present.
    # Fetch protein first so the original upstream importer can attach those
    # curated classes instead of leaving them absent in a brand-new store.
    selected = sorted(set(selected), key=lambda name: (name != "protein", name))
    capabilities = runtime_capabilities()
    if not capabilities["assembly_available"]:
        raise HydraRuntimeError("Install the complete native HYDRA / BLAST runtime before building reference databases.")
    root = Path(db_root).resolve()
    root.parent.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now(UTC).strftime("%Y%m%d-%H%M%S")
    destination = root.parent / f"{root.name}-snapshot-{stamp}-{uuid.uuid4().hex[:8]}"
    provider_version = _ncbi_version() if set(selected).intersection({"ncbi", "protein"}) else None
    with tempfile.TemporaryDirectory(prefix=".wmlstudio-hydra-stage-", dir=root.parent) as work:
        staged = Path(work) / "snapshot"
        try:
            command, log = _run_child(["db", "download", *selected, "--db-dir", str(staged)],
                                      work, cancelled, progress)
        except HydraRuntimeError as exc:
            # Nothing has been published at this point: the staging directory is
            # discarded with the temporary folder and the store in use is exactly
            # as it was. Saying so is the difference between a failed download and
            # a user who believes their reference data is now in an unknown state.
            raise HydraRuntimeError(
                f"Downloading {', '.join(selected)} failed, so no new snapshot was published and "
                f"the reference store you are using is unchanged.\n{exc}") from exc
        if provider_version:
            if _ncbi_version() != provider_version:
                raise HydraRuntimeError("NCBI changed its current database version during download; retry to obtain one consistent snapshot.")
            manifest_path = staged / "manifest.json"
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            for name in set(selected).intersection({"ncbi", "protein"}):
                entry = manifest.get("databases", {}).get(name)
                if entry:
                    entry["upstream_reported_version"] = entry.get("version", "unknown")
                    entry["version"] = provider_version
                    entry["version_source"] = NCBI_REFERENCE_ROOT + "/version.txt"
                    entry["source_url"] = NCBI_REFERENCE_ROOT
            manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        reference = _database_provenance(staged, selected, cancelled)
        if cancelled and cancelled():
            raise AnalysisCancelled("Database download cancelled before publishing; previous snapshots are unchanged.")
        os.replace(staged, destination)
    reference["root"] = str(destination)
    return {"command": command, "log": log, "reference_snapshot": reference,
            "database_root": str(destination), "previous_database_root": str(root),
            "snapshot_scope": "Selected databases only; previous snapshots are retained.",
            "completed_utc": datetime.now(UTC).isoformat(), "source_commit": HYDRA_COMMIT}


def main(argv=None):
    arguments = list(sys.argv[1:] if argv is None else argv)
    if not arguments or arguments.pop(0) != "--upstream":
        raise SystemExit("This worker is launched by WMLSTudio; use the desktop to select assemblies and databases.")
    from hydra_amr.cli import main as hydra_main
    sys.argv = ["hydra", *arguments]
    return hydra_main(arguments)


if __name__ == "__main__":
    raise SystemExit(main())
