"""The cohort AMR matrix: acquired determinants by isolate, and resistance point mutations.

Two Qt-free tables live here, so the evidence page can be tested without a screen.

THE DETERMINANT MATRIX is one row per acquired determinant and one column per isolate,
grouped by the antimicrobial class the reference curates the determinant under, because
that is the axis a shared pattern across an outbreak is read along. One isolate's column
is what that isolate's own report says, never what the cohort's reports say together.

THE POINT MUTATIONS are a separate table, because they are separate evidence. An
acquired gene is a gene the isolate either carries or does not; a resistance point
mutation is one catalogued substitution at one position of one gene, and the catalogue
it is read from is curated for one organism. So the mutation table names the gene, the
substitution and the organism catalogue it came from, and an isolate whose organism has
no catalogue in the installed release is said to have none rather than left blank beside
isolates that were screened.

Every cell states what was actually found, in words that survive printing and a
colour-blind reading: colour, where a screen adds it, is only a second carrier of what
the cell already says. Detected, detected as a partial or disrupted match, not detected
in a report that ran, and each separate reason the question was never asked of this
isolate -- no report, a report belonging to another input, a reference set that was not
searched, an organism with no mutation catalogue -- are distinct states and never
collapse into one another.

Nothing here is a susceptibility result. A determinant is sequence similarity to a
public reference sequence at a threshold; it is not an MIC, not a category, and no
determinant detected is not a susceptible isolate.
"""

from __future__ import annotations

import csv
import re
from collections.abc import Mapping
from datetime import UTC, datetime
from pathlib import Path

from . import __version__
from .sample_workflow import amr_search_scope, current_hydra_evidence, sample_organism

# The one sentence that has to reach every reader of this table, on screen and in an
# exported file. It is printed, not implied, and it is never softened.
NOT_A_PHENOTYPE = (
    "Not a susceptibility result. A determinant is sequence similarity to a public reference "
    "sequence; no determinant detected is not a susceptible isolate, and a determinant detected "
    "is not a resistant one.")

# Every state a determinant cell can hold, with the sentence the page prints beside it.
# `evidence` is the only thing that may be counted: 'detected' and 'absent' are answers a
# report gave, 'unknown' is a question that was never asked of this isolate and is never
# folded into either. The order is the order a legend reads in.
DETERMINANT_STATES = {
    "detected": {
        "label": "Detected", "evidence": "detected",
        "meaning": "A complete match to a reference determinant at this screen's thresholds. "
                   "Carriage of the sequence, not an expressed gene and not a phenotype."},
    "detected_partial": {
        "label": "Detected · partial", "evidence": "detected",
        "meaning": "Part of the reference was matched -- a contig edge, a truncation or a "
                   "fragment. That the determinant is intact is not established."},
    "detected_disrupted": {
        "label": "Detected · disrupted", "evidence": "detected",
        "meaning": "The match carries an internal stop or a disrupting change, so the reading "
                   "frame the reference describes is broken here. Whether anything is expressed "
                   "is not established."},
    "detected_unclassified": {
        "label": "Detected · completeness not recorded", "evidence": "detected",
        "meaning": "The report records a match but not how complete it is, so nothing is claimed "
                   "about its completeness."},
    "not_detected": {
        "label": "Not detected", "evidence": "absent",
        "meaning": "This isolate was screened against the reference sets recorded for it and this "
                   "determinant was not reported. It is not a susceptibility result and not "
                   "evidence that the isolate is susceptible to anything."},
    "not_searched": {
        "label": "Not in a set searched here", "evidence": "unknown",
        "meaning": "Every report naming this determinant in this cohort read it from a reference "
                   "set this isolate's run did not search, so it was never looked for here. "
                   "Unknown, never an absence."},
    "no_report": {
        "label": "No report", "evidence": "unknown",
        "meaning": "No AMR report is linked to this isolate, so nothing was searched for at all. "
                   "Unknown, never an absence."},
    "report_not_current": {
        "label": "Report not current", "evidence": "unknown",
        "meaning": "The stored report refers to an earlier or different input, so it says nothing "
                   "about the assembly in front of you. Rerun or review the input identity. "
                   "Unknown, never an absence."},
}
DETERMINANT_STATE_ORDER = tuple(DETERMINANT_STATES)

# Which detected state is shown when one isolate carries several matches of one
# determinant. The most complete match wins the cell; every state seen is listed in the
# cell's own reason, so a complete copy never hides a disrupted one.
_DETECTED_ORDER = ("detected", "detected_unclassified", "detected_partial", "detected_disrupted")

MUTATION_STATES = {
    "detected": {
        "label": "Detected", "evidence": "detected",
        "meaning": "The catalogued substitution was reported at this position for this isolate. "
                   "A resistance mutation in a catalogue is not a measured susceptibility."},
    "not_detected": {
        "label": "Not detected", "evidence": "absent",
        "meaning": "A point-mutation catalogue covering this isolate's organism was searched and "
                   "did not report this substitution. Not a susceptibility result."},
    "no_catalogue": {
        "label": "No catalogue for this organism", "evidence": "unknown",
        "meaning": "The installed reference release holds no point-mutation catalogue this "
                   "isolate could be screened against, so its mutations were never assessed. No "
                   "other organism's catalogue is substituted. Unknown, never an absence."},
    "not_searched": {
        "label": "Point mutations not searched", "evidence": "unknown",
        "meaning": "The run that produced this isolate's report did not search for point "
                   "mutations. Unknown, never an absence."},
    "target_not_assessed": {
        "label": "This target not assessed", "evidence": "unknown",
        "meaning": "This isolate's catalogue covers only the other level -- protein or DNA -- so "
                   "the target this row is read from was not examined for it. Unknown, never an "
                   "absence."},
    "catalogue_not_recorded": {
        "label": "Catalogue not recorded", "evidence": "unknown",
        "meaning": "This report does not record which point-mutation catalogue, if any, was "
                   "searched for this isolate, so nothing is claimed either way."},
    "no_report": {
        "label": "No report", "evidence": "unknown",
        "meaning": "No AMR report is linked to this isolate, so nothing was searched for at all. "
                   "Unknown, never an absence."},
    "report_not_current": {
        "label": "Report not current", "evidence": "unknown",
        "meaning": "The stored report refers to an earlier or different input, so it says nothing "
                   "about the assembly in front of you. Unknown, never an absence."},
}
MUTATION_STATE_ORDER = tuple(MUTATION_STATES)

MATRIX_LIMITATIONS = (
    NOT_A_PHENOTYPE,
    "A cell reads only the report linked to that isolate. An isolate without a current report is "
    "unknown evidence and is never counted as carrying nothing.",
    "Only the primary match of a locus is shown, as in the engine's own counts: the same locus "
    "reported again by a second reference set is redundant, not a second determinant.",
    "This is a screen over public reference data. It is not AMRFinderPlus, ResFinder, CARD or "
    "Kleborate, and it does not reproduce their curation, thresholds or reports.",
    "Two isolates sharing a row share a reference match. That is not transmission, not one "
    "mobile element and not a measure of relatedness.",
)

MUTATION_LIMITATIONS = (
    NOT_A_PHENOTYPE,
    "A resistance point mutation is one catalogued substitution in one gene, curated for the one "
    "organism its catalogue was written for. An isolate whose organism has no catalogue in the "
    "installed release was never screened for mutations, and no other organism's catalogue is "
    "substituted for the missing one.",
    "A substitution outside the installed release's catalogue cannot be reported by it, however "
    "well described it is elsewhere. This table is what that release holds, not what is known.",
    "The position and the substitution are the catalogue's own, read from the report. Nothing "
    "here re-reads the assembly or re-derives a coordinate.",
)

# Which catalogue level examines which kind of target, in the engine's own words. A level
# is about what could be looked in, never about what was found.
MUTATION_LEVEL_WORDS = {
    "dna_and_protein": "protein and DNA point-mutation catalogues were searched",
    "protein_only": "the protein catalogue was searched; no DNA-level target (23S rRNA and the "
                    "other non-coding targets) was assessed",
    "dna_only": "the DNA catalogue was searched; no curated protein mutation was assessed",
    "none": "no point-mutation catalogue is installed for this organism, so mutations were not "
            "assessed at all",
    "unknown": "this run recorded no point-mutation catalogue",
    "unrecorded": "this report records no point-mutation catalogue",
}
# The kinds of target each level can answer for. A hit from reads is reported by whichever
# catalogue the run held, so it is never claimed as a level gap.
_COVERED_KINDS = {
    "dna_and_protein": frozenset({"protein", "dna", "reads"}),
    "protein_only": frozenset({"protein", "reads"}),
    "dna_only": frozenset({"dna", "reads"}),
}
# The engine's own method names for a point-mutation call, and which reference each reads.
_MUTATION_KINDS = {"POINTX": "protein", "SUSCEPTIBLEX": "protein", "POINTN": "dna",
                   "POINTR": "reads", "VARIANTR": "reads"}
_KIND_WORDS = {"protein": "the protein catalogue", "dna": "the DNA catalogue",
               "reads": "reads", "unknown": "a target this report does not name"}
_RESOLUTION_STATES = {"COMPLETE": "detected", "PARTIAL": "detected_partial",
                      "INTERNAL_STOP": "detected_disrupted"}
# Upstream writes a mutation as `gene_S83L`, with the observed residue in brackets after
# it. The catalogue's own symbols are more varied than a substitution: `pbp4_T-266A`
# counts back from the end of a reference, `rplD_WR65del` is a deletion and
# `cirA_S90YfsTer15ins2` a frameshift, and all of them are catalogued changes this table
# has to name. So anything carrying a position is kept exactly as the catalogue wrote it,
# and anything else -- the divergent-from-susceptible finding, which is prose and not a
# change at a position -- is never turned into a substitution and is printed as it stands.
_OBSERVED = re.compile(r"\s*\(([^()]*)\)\s*$")
_SUBSTITUTION = re.compile(r"^[A-Za-z*]*-?\d+[A-Za-z0-9*]*$")

_NO_CLASS = "Class not recorded by the reference"


def _text(value):
    return str(value).strip() if value is not None else ""


def _hits(record):
    """The primary AMR hits of one isolate's current report, or none at all.

    ``current_hydra_evidence`` already withholds a report that belongs to another
    input, so a stale or missing one reaches here as no hits and the caller states
    which of the two it was.
    """
    evidence = current_hydra_evidence(record)
    hits = evidence.get("hits") if isinstance(evidence, Mapping) else None
    return [hit for hit in (hits or ()) if isinstance(hit, Mapping)
            and hit.get("primary") is True and _text(hit.get("element_type")).upper() == "AMR"]


def is_point_mutation(hit) -> bool:
    """Whether one AMR hit is a catalogued point mutation rather than an acquired gene."""
    return (_text(hit.get("element_subtype")).upper() == "POINT"
            or _text(hit.get("resolution")).upper() == "POINT")


def mutation_identity(hit) -> dict:
    """The gene, the substitution and the observed residue behind one point-mutation hit.

    Upstream records the catalogue's own symbol in the hit's note, optionally
    followed by what was actually observed. Nothing is invented here: a note that
    is not a substitution symbol -- the divergent-from-susceptible finding, for
    one -- keeps an empty substitution and is reported in the catalogue's own
    words instead.
    """
    note = _text(hit.get("note"))
    gene = _text(hit.get("gene"))
    observed, symbol = "", note
    match = _OBSERVED.search(note)
    if match:
        observed, symbol = match.group(1).strip(), note[:match.start()].strip()
    core = symbol
    if gene and core.startswith(gene + "_"):
        core = core[len(gene) + 1:]
    elif "_" in core:
        core = core.rpartition("_")[2]
    recorded = bool(core) and bool(_SUBSTITUTION.fullmatch(core))
    return {"gene": gene, "symbol": symbol if recorded else "",
            "substitution": core if recorded else "", "observed": observed,
            "recorded": recorded, "finding": symbol if recorded else note,
            "kind": _MUTATION_KINDS.get(_text(hit.get("method")).upper(), "unknown"),
            "catalogue_entry": _text(hit.get("product")), "accession": _text(hit.get("accession"))}


def _determinant_state(hit) -> str:
    """How complete the match is, from the report's own resolution and method."""
    method = _text(hit.get("method")).upper()
    subtype = _text(hit.get("element_subtype")).upper()
    resolution = _text(hit.get("resolution")).upper()
    if "DISRUPT" in method or "DISRUPT" in subtype:
        return "detected_disrupted"
    if resolution in _RESOLUTION_STATES:
        return _RESOLUTION_STATES[resolution]
    if "PARTIAL" in method:
        return "detected_partial"
    return "detected_unclassified"


def _percent(hit, key):
    value = hit.get(key)
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        return None
    return float(value)


def _match_sentence(hits, states):
    """One line saying how many matches there are, and what each of them was."""
    identities = [value for value in (_percent(hit, "identity_pct") for hit in hits)
                  if value is not None]
    coverages = [value for value in (_percent(hit, "coverage_pct") for hit in hits)
                 if value is not None]
    parts = [f"{len(hits)} primary match" + ("es" if len(hits) != 1 else "")]
    if identities and coverages:
        parts.append(f"best {max(identities):.1f}% identity over {max(coverages):.1f}% of the "
                     "reference")
    databases = sorted({_text(hit.get("database")) for hit in hits if _text(hit.get("database"))})
    if databases:
        parts.append("reported by " + ", ".join(databases))
    seen = [DETERMINANT_STATES[state]["label"] for state in _DETECTED_ORDER
            if state in states and (len(states) > 1 or state != "detected")]
    if seen:
        parts.append("matches read as " + "; ".join(seen))
    return ". ".join(parts) + "."


def _isolates(records):
    """One column per isolate: its identity, its organism and what its report could answer."""
    columns = []
    for record in records:
        if not isinstance(record, Mapping):
            raise ValueError("Every isolate in an AMR matrix must be a project sample record.")
        result = record.get("result") if isinstance(record.get("result"), Mapping) else {}
        sample_id = _text(record.get("id") or record.get("sample_id"))
        if not sample_id:
            raise ValueError("Every isolate in an AMR matrix needs a sample id.")
        organism = sample_organism(record)
        columns.append({
            "sample_id": sample_id,
            "sample_name": _text(record.get("name") or result.get("sample_name")) or sample_id,
            "organism": organism["organism"] or "Organism not established",
            "organism_source": organism["source"],
            "scope": amr_search_scope(record),
            "hits": _hits(record),
        })
    identifiers = [column["sample_id"] for column in columns]
    repeated = sorted({value for value in identifiers if identifiers.count(value) > 1})
    if repeated:
        raise ValueError("Each isolate appears once in an AMR matrix. Repeated identifiers: "
                         + ", ".join(repeated))
    columns.sort(key=lambda column: (column["sample_name"].casefold(), column["sample_id"]))
    names = [column["sample_name"] for column in columns]
    for column in columns:
        # Two isolates may legitimately share a name; an exported column heading may not,
        # so the ambiguous ones carry their own identifier and the rest stay readable.
        column["column_title"] = (column["sample_name"] if names.count(column["sample_name"]) == 1
                                  else f"{column['sample_name']} · {column['sample_id'][:8]}")
    return columns


def _unknown_cell(scope):
    """The one unknown a whole column is in, before any row is considered."""
    if scope["status"] == "missing":
        return {"state": "no_report", "reason": scope["reason"]}
    if scope["status"] == "stale":
        return {"state": "report_not_current", "reason": scope["reason"]}
    return None


def _cell(state, reason, states, totals, hits=0):
    states[state] = states.get(state, 0) + 1
    totals[state] = totals.get(state, 0) + 1
    return {"state": state, "display": "", "reason": reason, "hits": hits}


def _finish(row, columns, states, words):
    """Fill in each cell's printed words and the row's own numerator and denominator."""
    detected = sum(count for state, count in states.items()
                   if words[state]["evidence"] == "detected" and count)
    answered = sum(count for state, count in states.items()
                   if words[state]["evidence"] in {"detected", "absent"} and count)
    for cell in row["cells"].values():
        cell["display"] = words[cell["state"]]["label"]
    row.update(states=states, detected_in=detected, answered=answered,
               isolates=len(columns), unknown=len(columns) - answered,
               detected_summary=f"{detected} of {answered} screened"
                                + (f" · {len(columns) - answered} unknown"
                                   if answered != len(columns) else ""))
    return row


def build_determinant_matrix(records) -> dict:
    """One row per acquired determinant, one column per isolate, grouped by drug class.

    ``records`` are project sample records. Point mutations are deliberately not here:
    they are a different kind of evidence and are built by
    :func:`build_mutation_matrix` instead.
    """
    columns = _isolates(records)
    if not columns:
        raise ValueError("Choose at least one isolate to build an AMR determinant matrix.")
    grouped: dict[str, dict] = {}
    for column in columns:
        for hit in column["hits"]:
            if is_point_mutation(hit):
                continue
            gene = _text(hit.get("gene"))
            if not gene:
                continue
            row = grouped.setdefault(gene, {"gene": gene, "classes": set(), "subclasses": set(),
                                            "databases": set(), "products": set(), "found": {}})
            row["classes"].add(_text(hit.get("class")))
            row["subclasses"].add(_text(hit.get("subclass")))
            row["databases"].add(_text(hit.get("database")))
            row["products"].add(_text(hit.get("product")))
            row["found"].setdefault(column["sample_id"], []).append(hit)
    totals = dict.fromkeys(DETERMINANT_STATE_ORDER, 0)
    rows = []
    for gene, found in sorted(grouped.items(), key=lambda item: item[0].casefold()):
        classes = sorted(value for value in found["classes"] if value)
        group = "; ".join(classes) or _NO_CLASS
        cells, states = {}, dict.fromkeys(DETERMINANT_STATE_ORDER, 0)
        for column in columns:
            scope, hits = column["scope"], found["found"].get(column["sample_id"]) or []
            unknown = _unknown_cell(scope)
            if hits:
                seen = {_determinant_state(hit) for hit in hits}
                state = next(value for value in _DETECTED_ORDER if value in seen)
                cells[column["sample_id"]] = _cell(state, _match_sentence(hits, seen), states,
                                                   totals, len(hits))
            elif unknown is not None:
                cells[column["sample_id"]] = _cell(unknown["state"], unknown["reason"], states,
                                                   totals)
            elif scope["databases"] and found["databases"] - {""} and not (
                    found["databases"] & set(scope["databases"])):
                cells[column["sample_id"]] = _cell(
                    "not_searched",
                    f"{gene} is reported here from " + ", ".join(sorted(found["databases"] - {""}))
                    + "; this isolate's run searched " + ", ".join(scope["databases"]) + ".",
                    states, totals)
            else:
                cells[column["sample_id"]] = _cell(
                    "not_detected",
                    "Screened against "
                    + (", ".join(scope["databases"]) or "reference sets this report does not name")
                    + f"; {gene} was not reported for this isolate.", states, totals)
        rows.append(_finish({"key": gene, "group": group, "label": gene,
                             "classes": classes,
                             "subclass": "; ".join(sorted(value for value in found["subclasses"]
                                                          if value)),
                             "databases": sorted(value for value in found["databases"] if value),
                             "product": "; ".join(sorted(value for value in found["products"]
                                                         if value)),
                             "cells": cells}, columns, states, DETERMINANT_STATES))
    rows.sort(key=lambda row: (row["group"] == _NO_CLASS, row["group"].casefold(),
                               row["label"].casefold()))
    return _table("determinants", columns, rows, DETERMINANT_STATES, DETERMINANT_STATE_ORDER,
                  totals, MATRIX_LIMITATIONS,
                  leading=({"key": "group", "title": "Antimicrobial class"},
                           {"key": "label", "title": "Determinant"},
                           {"key": "subclass", "title": "Reference subclass"},
                           {"key": "detected_summary", "title": "Detected in"}))


def _mutation_state(scope, row, hits):
    """Why this isolate has no row here -- which is four different answers, not one."""
    if hits:
        best = max(hits, key=lambda hit: (_percent(hit, "identity_pct") or 0.0))
        identity = _percent(best, "identity_pct")
        observed = mutation_identity(best)["observed"]
        return "detected", ". ".join(filter(None, [
            f"Reported as {mutation_identity(best)['finding'] or row['label']}",
            f"observed {observed}" if observed else "",
            f"{identity:.1f}% identity to the catalogue reference" if identity is not None else "",
            _text(best.get("database")) and f"read from the {_text(best.get('database'))} set"])) + "."
    unknown = _unknown_cell(scope)
    if unknown is not None:
        return unknown["state"], unknown["reason"]
    level = scope["point_mutation_level"]
    if scope["point_mutations"] is False:
        return "not_searched", ("This isolate's run was asked not to search for point mutations, "
                                "so none was reported for it.")
    if level == "none":
        return "no_catalogue", scope["point_mutation_reason"] or (
            "The installed release holds no point-mutation catalogue for this isolate's organism, "
            "so it was screened for genes only.")
    if level not in _COVERED_KINDS:
        return "catalogue_not_recorded", (
            scope["point_mutation_reason"]
            or "This report records no point-mutation catalogue for this isolate.")
    kinds = row["kinds"] - {"unknown"}
    if kinds and not (kinds & _COVERED_KINDS[level]):
        return "target_not_assessed", (
            "This row is read from "
            + "; ".join(_KIND_WORDS[kind] for kind in sorted(kinds))
            + "; for this isolate " + MUTATION_LEVEL_WORDS[level] + ".")
    return "not_detected", (f"Screened against the point-mutation catalogue for "
                            f"{scope['organism'] or 'this isolate'}: "
                            + MUTATION_LEVEL_WORDS.get(level, "a catalogue was searched")
                            + f"; {row['label']} was not reported.")


def build_mutation_matrix(records) -> dict:
    """One row per catalogued substitution, one column per isolate, with its catalogue named.

    The catalogue an isolate was screened against is the isolate's own, recorded by the
    run that produced its report. An isolate whose organism has no catalogue in the
    installed release carries that answer in every cell, and again in
    ``catalogue_coverage``, so an empty section can never be read as "none found".
    """
    columns = _isolates(records)
    if not columns:
        raise ValueError("Choose at least one isolate to build a point-mutation table.")
    grouped: dict[tuple[str, str], dict] = {}
    for column in columns:
        for hit in column["hits"]:
            if not is_point_mutation(hit):
                continue
            identity = mutation_identity(hit)
            gene = identity["gene"] or _text(hit.get("gene"))
            finding = identity["finding"] or _text(hit.get("note"))
            if not gene and not finding:
                continue
            row = grouped.setdefault((gene, finding), {
                "gene": gene, "finding": finding, "substitution": identity["substitution"],
                "recorded": identity["recorded"], "classes": set(), "subclasses": set(),
                "kinds": set(), "entries": set(), "catalogues": set(), "found": {}})
            row["classes"].add(_text(hit.get("class")))
            row["subclasses"].add(_text(hit.get("subclass")))
            row["kinds"].add(identity["kind"])
            row["entries"].add(identity["catalogue_entry"])
            if column["scope"]["organism"]:
                row["catalogues"].add(column["scope"]["organism"])
            row["found"].setdefault(column["sample_id"], []).append(hit)
    totals = dict.fromkeys(MUTATION_STATE_ORDER, 0)
    rows = []
    for (gene, finding), found in sorted(grouped.items(),
                                         key=lambda item: (item[0][0].casefold(),
                                                           item[0][1].casefold())):
        classes = sorted(value for value in found["classes"] if value)
        label = f"{gene} {found['substitution']}".strip() if found["recorded"] else (
            f"{gene}: {finding}".strip(": ") or finding)
        row = {"key": f"{gene}|{finding}", "group": "; ".join(classes) or _NO_CLASS,
               "label": label, "gene": gene,
               "substitution": found["substitution"] if found["recorded"]
               else "no substitution recorded — see the finding",
               "finding": finding, "kinds": set(found["kinds"]),
               "target": "; ".join(_KIND_WORDS[kind] for kind in sorted(found["kinds"] - {"unknown"}))
               or "not recorded by the report",
               "catalogue": "; ".join(sorted(found["catalogues"])) or "not recorded by the report",
               "catalogue_entry": "; ".join(sorted(value for value in found["entries"] if value)),
               "subclass": "; ".join(sorted(value for value in found["subclasses"] if value)),
               "cells": {}}
        states = dict.fromkeys(MUTATION_STATE_ORDER, 0)
        for column in columns:
            state, reason = _mutation_state(column["scope"], row,
                                            found["found"].get(column["sample_id"]) or [])
            row["cells"][column["sample_id"]] = _cell(
                state, reason, states, totals,
                len(found["found"].get(column["sample_id"]) or []))
        row["kinds"] = sorted(row["kinds"])
        rows.append(_finish(row, columns, states, MUTATION_STATES))
    rows.sort(key=lambda row: (row["group"] == _NO_CLASS, row["group"].casefold(),
                               row["gene"].casefold(), row["label"].casefold()))
    table = _table("point_mutations", columns, rows, MUTATION_STATES, MUTATION_STATE_ORDER,
                   totals, MUTATION_LIMITATIONS,
                   leading=({"key": "group", "title": "Antimicrobial class"},
                            {"key": "gene", "title": "Gene"},
                            {"key": "substitution", "title": "Substitution"},
                            {"key": "catalogue", "title": "Organism catalogue"},
                            {"key": "target", "title": "Read from"},
                            {"key": "detected_summary", "title": "Detected in"}))
    table["catalogue_coverage"] = [_coverage(column) for column in columns]
    table["without_catalogue"] = [row["sample_name"] for row in table["catalogue_coverage"]
                                  if not row["searched"]]
    table["catalogue_note"] = _catalogue_note(table["catalogue_coverage"])
    return table


def _coverage(column) -> dict:
    """Whether one isolate's mutations were assessed at all, and in the run's own words."""
    scope = column["scope"]
    level = scope["point_mutation_level"]
    searched = bool(scope["has_evidence"] and scope["point_mutations"] is not False
                    and level in _COVERED_KINDS)
    if not scope["has_evidence"]:
        reason = scope["reason"]
    elif scope["point_mutations"] is False:
        reason = ("This isolate's run was asked not to search for point mutations. That is not "
                  "evidence that it carries none.")
    elif scope["point_mutation_reason"]:
        # The engine's own sentence about this isolate, printed as it wrote it.
        reason = scope["point_mutation_reason"]
    elif level in _COVERED_KINDS:
        reason = ("The catalogue for " + (scope["organism"] or "this organism")
                  + " was chosen by the run that produced this report: "
                  + MUTATION_LEVEL_WORDS[level] + ".")
    else:
        reason = ("This report records no point-mutation catalogue for this isolate, so nothing "
                  "is claimed about mutations either way.")
    return {"sample_id": column["sample_id"], "sample_name": column["sample_name"],
            "column_title": column["column_title"], "organism": column["organism"],
            "organism_source": column["organism_source"],
            "catalogue": scope["organism"] or "none chosen",
            "level": level, "level_words": MUTATION_LEVEL_WORDS.get(level, MUTATION_LEVEL_WORDS["unrecorded"]),
            "evidence_state": scope["status"], "searched": searched, "reason": reason}


def _catalogue_note(coverage) -> str:
    """The sentence above the table, naming the isolates it can say nothing about."""
    missing = [row["sample_name"] for row in coverage if not row["searched"]]
    if not missing:
        return (f"All {len(coverage)} isolate(s) here were screened against a point-mutation "
                "catalogue chosen for their own organism. A catalogued mutation is not a "
                "measured susceptibility.")
    return (f"{len(missing)} of {len(coverage)} isolate(s) were not screened for point mutations "
            "at all: " + "; ".join(sorted(missing)) + ". Their columns say so rather than reading "
            "as an absence, and no other organism's catalogue is substituted for the missing one.")


def _table(kind, columns, rows, words, order, totals, limitations, *, leading) -> dict:
    answered = sum(row["answered"] for row in rows)
    return {
        "format_version": 1, "kind": kind,
        "isolate_count": len(columns), "row_count": len(rows),
        "samples": [{key: column[key] for key in
                     ("sample_id", "sample_name", "column_title", "organism", "organism_source")}
                    | {"scope": column["scope"],
                       "detected": sum(1 for row in rows
                                       if words[row["cells"][column["sample_id"]]["state"]]
                                       ["evidence"] == "detected")}
                    for column in columns],
        "leading": [dict(entry) for entry in leading],
        "rows": rows, "states": {state: dict(words[state]) for state in order},
        "state_order": list(order), "state_totals": dict(totals),
        "cells": len(rows) * len(columns), "answered_cells": answered,
        "unknown_cells": len(rows) * len(columns) - answered,
        "headline": NOT_A_PHENOTYPE,
        "denominator_note": (
            f"{len(rows)} row(s) × {len(columns)} isolate(s). Each row's count is over the "
            "isolates whose own report could answer for it, never over the cohort: "
            f"{len(rows) * len(columns) - answered} of {len(rows) * len(columns)} cells are "
            "unknown rather than negative."),
        "limitations": list(limitations),
        "provenance": {"software": "WMLSTudio", "version": __version__,
                       "created_utc": datetime.now(UTC).isoformat()},
    }


def matrix_table(table) -> tuple[list[str], list[list[str]]]:
    """The built matrix as a header row and text rows, exactly as it is read on screen."""
    headers = [entry["title"] for entry in table["leading"]] + [column["column_title"]
                                                                for column in table["samples"]]
    rows = []
    for row in table["rows"]:
        values = [_text(row.get(entry["key"])) for entry in table["leading"]]
        values += [row["cells"][column["sample_id"]]["display"] for column in table["samples"]]
        rows.append(values)
    return headers, rows


def write_matrix(table, path) -> Path:
    """Write one built matrix as a delimited file, with every distinction still in words.

    The cells carry their own state as text, so a printed or colour-blind reading of
    the export is the reading on screen. The sentences that say what the table is not
    follow the data, after one blank row, rather than being left on a screen the file
    will be read away from.
    """
    destination = Path(path).expanduser().resolve()
    from .export import _atomic_text, _csv_cell
    headers, rows = matrix_table(table)
    with _atomic_text(destination, newline="") as handle:
        writer = csv.writer(handle,
                            delimiter="," if destination.suffix.casefold() == ".csv" else "\t")
        # The isolate names are column headings here, so the heading row is guarded the
        # same way the cells are: a sample called "=cmd" must not reach a spreadsheet as
        # a formula because it happened to be written across the top.
        writer.writerow([_csv_cell(value) for value in headers])
        writer.writerows([_csv_cell(value) for value in row] for row in rows)
        writer.writerow([])
        for note in [table["denominator_note"], *table["limitations"]]:
            writer.writerow([_csv_cell(note)])
        writer.writerow([])
        writer.writerow([_csv_cell("What each cell state means")])
        for state in table["state_order"]:
            words = table["states"][state]
            writer.writerow([_csv_cell(words["label"]), _csv_cell(words["meaning"])])
    return destination
