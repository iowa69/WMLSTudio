"""Stable-identity evidence linking and flattened features for native views."""

from __future__ import annotations

import copy
import re
from pathlib import Path

from wmlstudio.sequence import file_sha256, sample_name

CONFIGURATION_VERSION = 1
# One sample configuration, under these names and no other. A menu that keeps its
# own copy of "which scheme" or "which organism" is how two menus come to disagree
# about one isolate; every surface reads sample_configuration and writes through
# storage.set_sample_configuration.
ORGANISM_FIELDS = ("genus", "species")
WORKFLOW_FIELDS = ("typing_mode", "scheme_path", "calling_mode", "genetic_code",
                   "cgmlst_scheme_key", "cgmlst_scheme_path", "run_hydra", "run_cgmlst")
CONFIGURATION_FIELDS = (*ORGANISM_FIELDS, *WORKFLOW_FIELDS)
TYPING_MODES = ("auto", "manual", "unknown")
CALLING_MODES = ("auto", "exact", "full_cds")

# An MLST panel match is compatibility with an installed panel. The caveat travels
# with the value it produced, so no later surface can read it as taxonomy.
MLST_ORGANISM_CAVEAT = (
    "An MLST scheme match is evidence for a lineage compatible with an installed panel, "
    "not an independent species determination.")

# The four things an allele call is made against. Nothing else can change the call,
# and a rerun is asked for only when one of them differs from what was stored.
RERUN_FIELDS = ("input_sha256", "scheme_digest", "scheme", "scheme_version")
_RERUN_LABELS = {"input_sha256": "the sequence input", "scheme_digest": "the scheme's contents",
                 "scheme": "the scheme", "scheme_version": "the scheme version"}
_SCHEME_VERSION_KEYS = ("revision", "version", "last_updated", "scheme_version")


def _sha256(value):
    return value.lower() if isinstance(value, str) and re.fullmatch(r'[0-9a-fA-F]{64}', value) else None


def current_input_sha256(record):
    """Best recorded current-input identity; no filesystem reads in presentation code."""
    result = record.get('result') or (record if 'result' not in record else {})
    metadata = record.get('metadata') or {}
    assembly = metadata.get('assembly') or {}
    workflow = metadata.get('workflow') or {}
    return (_sha256(result.get('input_sha256'))
            or _sha256((assembly.get('provenance') or {}).get('assembly_sha256'))
            or _sha256(workflow.get('managed_sha256'))
            or _sha256((metadata.get('input_identity') or {}).get('sha256')))


def _execution_input_sha256(execution, source_name, *, single_sample=False):
    inputs = execution.get('inputs') or []
    if not isinstance(inputs, list):
        return None
    if single_sample and len(inputs) == 1 and isinstance(inputs[0], dict):
        return _sha256(inputs[0].get('sha256'))
    matches = []
    for entry in inputs:
        if not isinstance(entry, dict):
            continue
        names = {entry.get('sample'), entry.get('sample_id'), entry.get('sample_name')}
        if entry.get('path'):
            names.add(sample_name(Path(str(entry['path']))))
        if source_name in names and _sha256(entry.get('sha256')):
            matches.append(_sha256(entry['sha256']))
    return matches[0] if len(matches) == 1 else None


def hydra_evidence_status(record):
    """Classify stored AMR evidence without presenting old inputs as current ones.

    Link-time baseline hashes detect later changes but never establish identity
    of an imported report that did not supply its own assembly hash.
    """
    metadata = record.get('metadata') or {}
    evidence = metadata.get('hydra') or {}
    current = current_input_sha256(record)
    if not isinstance(evidence, dict) or not evidence:
        return {'status': 'missing', 'reason': 'No linked HYDRA report.',
                'input_sha256': current, 'evidence_input_sha256': None}
    source_hash = _sha256(evidence.get('evidence_input_sha256'))
    if source_hash is None and 'evidence_input_sha256' not in evidence:
        source_hash = _execution_input_sha256(evidence.get('execution_provenance') or {},
                                              evidence.get('source_sample'), single_sample=True)
    baseline = _sha256(evidence.get('linked_input_sha256'))
    reference = source_hash or baseline
    if current and reference and current != reference:
        status, reason = 'stale', 'Archived HYDRA evidence refers to an earlier or different input; rerun or explicitly relink after reviewing the input identity.'
    elif (record.get('job_status') or (record.get('status') if 'result' in record else None)) in {'failed', 'interrupted', 'running', 'queued'}:
        status, reason = 'unverified', 'Saved HYDRA evidence is retained, but the current sample input has no completed active analysis. Review the input identity before treating this as current evidence.'
    elif current and source_hash and current == source_hash:
        status, reason = 'current', 'HYDRA assembly SHA-256 agrees with the recorded current typing input.'
    else:
        status, reason = 'unverified', 'User-mapped HYDRA evidence; no matching report/input SHA-256 confirmation is available. A link-time baseline is not report identity verification.'
    return {'status': status, 'reason': reason, 'input_sha256': current,
            'evidence_input_sha256': source_hash, 'linked_input_sha256': baseline}


def current_hydra_evidence(record):
    """Visible evidence, including explicitly unverified mappings, excluding stale records."""
    if hydra_evidence_status(record)['status'] in {'stale', 'missing'}:
        return {}
    evidence = (record.get('metadata') or {}).get('hydra') or {}
    return evidence if isinstance(evidence, dict) else {}


def amr_search_scope(record) -> dict:
    """What one isolate's stored AMR report was actually able to report, from its own provenance.

    A determinant list means nothing without the question it answered. Which
    reference sets were read, which release they were, whether point mutations
    were searched at all and which organism catalogue was chosen are all recorded
    by the run itself, so they are read back rather than inferred: an older or
    imported report that recorded none of it says "unrecorded" instead of
    borrowing the settings of today's run.

    ``point_mutation_level`` is the engine's own answer -- dna_and_protein,
    protein_only, dna_only, none -- and 'none' is not 'no mutations'. It means
    the installed release holds no catalogue this isolate could be screened
    against, which is unknown evidence, not a clean result.
    """
    state = hydra_evidence_status(record)
    evidence = current_hydra_evidence(record)
    execution = evidence.get('execution_provenance') or {}
    organism = execution.get('organism') if isinstance(execution.get('organism'), dict) else {}
    virulence = execution.get('virulence') if isinstance(execution.get('virulence'), dict) else {}
    release = execution.get('reference_release') if isinstance(execution.get('reference_release'), dict) else {}
    names = evidence.get('databases')
    if not isinstance(names, (list, tuple)):
        snapshot = execution.get('reference_snapshot') or {}
        names = list((snapshot.get('databases') or {})) if isinstance(snapshot, dict) else []
    requested = organism.get('point_mutations')
    return {
        'status': state['status'], 'reason': state['reason'],
        'has_evidence': bool(evidence),
        'databases': sorted({str(name) for name in names if str(name)}),
        'release': str(release.get('release') or ''),
        'organism': str(organism.get('resolved') or ''),
        'organism_requested': str(organism.get('requested') or ''),
        'point_mutations': bool(requested) if isinstance(requested, bool) else None,
        'point_mutation_level': str(organism.get('point_mutation_level') or '') or 'unrecorded',
        'point_mutation_reason': str(organism.get('reason') or ''),
        'virulence_enabled': bool(virulence['enabled']) if 'enabled' in virulence else None,
        'virulence_curated': bool(virulence.get('organism_curated')) if virulence else None,
    }


def select_records(records, selected_ids=None):
    records = list(records)
    if selected_ids is None:
        return records
    requested = set(map(str, selected_ids))
    available = {str(row.get("id") or row.get("sample_id")) for row in records}
    missing = requested - available
    if missing:
        raise KeyError(f"Selected sample identifiers were not found: {', '.join(sorted(missing))}")
    return [row for row in records if str(row.get("id") or row.get("sample_id")) in requested]


def suggest_hydra_links(project, report) -> dict:
    """Suggest exact name/stem links; ambiguous matches never become a mapping."""
    samples = project.samples()
    mapping, ambiguous, missing = {}, {}, []
    source_names = [sample["sample"] for sample in report.get("samples", [])]
    if len(source_names) != len(set(source_names)):
        raise ValueError("HYDRA report contains duplicate sample names.")
    for name in source_names:
        normalized = sample_name(Path(name)).casefold()
        candidates = [sample["id"] for sample in samples if normalized in {
            sample_name(Path(sample["name"])).casefold(),
            sample_name(Path(sample["input_path"])).casefold() if sample.get("input_path") else "",
        }]
        if len(candidates) == 1:
            mapping[name] = candidates[0]
        elif candidates:
            ambiguous[name] = candidates
        else:
            missing.append(name)
    # Multiple report names can normalize to one isolate. Such suggestions must
    # also be resolved explicitly instead of letting the last report row win.
    for sample_id in set(mapping.values()):
        sources = [name for name, target in mapping.items() if target == sample_id]
        if len(sources) > 1:
            for name in sources:
                ambiguous[name] = [sample_id]
                del mapping[name]
    return {"mapping": mapping, "ambiguous": ambiguous, "missing": missing}


def link_hydra(project, report, mapping: dict[str, str]) -> list[str]:
    """Link only explicit upstream-name → stable-ID mappings in one transaction."""
    upstream = report.get("samples", [])
    by_name = {sample["sample"]: sample for sample in upstream}
    if len(by_name) != len(upstream):
        raise ValueError("HYDRA report contains duplicate sample names.")
    if len(set(mapping.values())) != len(mapping):
        raise ValueError("Multiple HYDRA samples cannot be linked to one project isolate.")
    provenance = report.get("import_provenance", {})
    digest = provenance.get("sha256")
    if not digest:
        raise ValueError("Import and validate the HYDRA report before linking its evidence.")
    for source_name, sample_id in mapping.items():
        if source_name not in by_name:
            raise KeyError(f"Unknown HYDRA sample: {source_name}")
        project.get_sample(sample_id)
    with project.transaction():
        for source_name, sample_id in mapping.items():
            source = copy.deepcopy(by_name[source_name])
            sample = project.get_sample(sample_id)
            baseline = current_input_sha256(sample)
            if baseline is None and sample.get('input_path'):
                try:
                    # An optional link-time baseline must not synchronously scan
                    # an enormous raw-read dataset from the report-mapping UI.
                    if Path(sample['input_path']).stat().st_size <= 32 * 1024 * 1024:
                        baseline = file_sha256(sample['input_path'])
                except OSError:
                    pass
            source_hash = (_sha256(source.get('input_sha256'))
                           or _execution_input_sha256(report.get('execution_provenance') or {},
                                                     source_name, single_sample=len(upstream) == 1))
            evidence = {
                "report_sha256": digest, "source_sample": source_name,
                "summary": source.get("summary", {}), "hits": source.get("hits", []),
                "mlst": source.get("mlst", {}), "species": source.get("species", {}),
                "provenance": copy.deepcopy(provenance), "upstream": source,
                "execution_provenance": copy.deepcopy(report.get("execution_provenance", {})),
                "databases": copy.deepcopy(report.get("databases", [])),
                'evidence_input_sha256': source_hash,
                'linked_input_sha256': baseline,
                'link_identity_basis': 'report assembly SHA-256' if source_hash else 'explicit user mapping; report input hash unavailable',
            }
            metadata = sample["metadata"]
            if metadata.get('hydra'):
                project.record_history(sample_id, 'hydra_superseded', {'evidence': metadata['hydra']})
            metadata["hydra"] = evidence
            project.set_metadata(sample_id, metadata)
            project.record_history(sample_id, "hydra_linked", {
                "report_sha256": digest, "source_sample": source_name,
                'report_input_sha256': source_hash, 'linked_input_sha256': baseline,
            })
    return list(mapping.values())


def set_cluster(project, sample_ids, label: str, color: str = "#2F8A78", highlight: bool = True):
    if not re.fullmatch(r"#[0-9A-Fa-f]{6}", color):
        raise ValueError("Choose a cluster colour in #RRGGBB format.")
    with project.transaction():
        for sample_id in dict.fromkeys(sample_ids):
            project.update_metadata(sample_id, {"cluster": {
                "label": str(label), "color": color, "highlight": bool(highlight),
                "meaning": "User-defined grouping; not a transmission inference",
            }})


def typing_profiles(record) -> dict:
    """The classical MLST and the cgMLST profile this record carries, kept apart.

    Read from what the record already holds — its headline result and any
    additional profiles attached to it — so presentation code never has to ask
    the database which of the two a scheme was. Either may be None: a sample with
    an ST and no core-genome profile reports no cgMLST, and the ST is never
    offered in its place.
    """
    from wmlstudio.project import typing_kind
    result = record.get("result") or (record if "result" not in record else {})
    profiles: dict[str, dict | None] = {"mlst": None, "cgmlst": None}
    seen = set()
    for profile in [result, *(record.get("analyses") or ()),
                    *(record.get("additional_profiles") or ())]:
        if not isinstance(profile, dict) or not profile.get("scheme_digest"):
            continue
        digest = str(profile["scheme_digest"])
        if digest in seen:
            continue
        seen.add(digest)
        kind = typing_kind(profile)
        if kind in profiles and profiles[kind] is None:
            profiles[kind] = profile
    return profiles


def _metadata(record) -> dict:
    metadata = record.get("metadata") if isinstance(record, dict) else None
    return metadata if isinstance(metadata, dict) else {}


def assigned_organism(record) -> dict:
    """The organism stored on the sample, tolerating the old free-text spelling of it.

    Only what was written down: unlike sample_organism this never falls back to what
    an analysis detected, so a write can tell a stored assignment from a detection
    that is merely being read through.
    """
    assigned = _metadata(record).get("organism") or {}
    if isinstance(assigned, str):
        parts = assigned.split(maxsplit=1)
        assigned = {"genus": parts[0] if parts else "", "species": parts[1] if len(parts) > 1 else ""}
    return assigned if isinstance(assigned, dict) else {}


def sample_organism(record) -> dict:
    """Genus and species for one sample, with where the value came from and how strong it is.

    The single read every surface uses, so nothing asks again for what the first
    analysis already established. An assignment somebody made outranks a detection,
    a detection is reported as a detection, and the basis travels with the value:
    an MLST scheme match names a lineage compatible with an installed panel, never a
    species. Genus and species always come from one source, never one from each.
    """
    metadata = _metadata(record)
    assigned = assigned_organism(record)
    evidence = metadata.get("organism_evidence")
    evidence = evidence if isinstance(evidence, dict) else {}
    result = record.get("result") or (record if "result" not in record else {})
    result = result if isinstance(result, dict) else {}
    detected = result.get("identification") or {}
    detected = detected if isinstance(detected, dict) else {}
    named = detected.get("organism") if isinstance(detected.get("organism"), dict) else {}
    genus, species = str(assigned.get("genus") or ""), str(assigned.get("species") or "")
    if genus:
        source = "assigned"
        basis = str(evidence.get("basis") or "")
        confidence = str(evidence.get("confidence") or "")
    else:
        genus = str(named.get("genus") or detected.get("genus") or "")
        species = str(named.get("species") or detected.get("species") or "")
        # The basis describes the value, so a detection never borrows the basis of a
        # decision record about some other proposal. The only engine that writes an
        # identification into a typing result is identification.identify_assembly.
        source = "scheme_match" if genus else "unknown"
        basis = str(detected.get("basis") or "mlst_panel") if genus else ""
        confidence = str(detected.get("confidence") or "panel_compatibility") if genus else ""
    authored = _authorship_record(metadata).get("genus") or {}
    return {"genus": genus, "species": species,
            "organism": " ".join(part for part in (genus, species) if part),
            "source": source, "basis": basis, "confidence": confidence,
            "status": str(evidence.get("status") or ""),
            "note": MLST_ORGANISM_CAVEAT if basis == "mlst_panel" else "",
            "set_by": str(authored.get("by") or evidence.get("confirmed_by") or ""),
            "set_utc": str(authored.get("utc") or evidence.get("confirmed_utc") or ""),
            "surface": str(authored.get("surface") or ""),
            "automatic": bool(authored["automatic"]) if "automatic" in authored
            else evidence.get("confirmed_by") in {None, "", "auto_policy"}}


def _authorship_record(metadata) -> dict:
    record = metadata.get("configuration") if isinstance(metadata, dict) else None
    fields = record.get("fields") if isinstance(record, dict) else None
    if not isinstance(fields, dict):
        return {}
    return {field: dict(entry) for field, entry in fields.items()
            if field in CONFIGURATION_FIELDS and isinstance(entry, dict)}


def configuration_authorship(record) -> dict:
    """Who set each configuration field and when, as far as the sample records it.

    A field nobody has touched is simply absent. The organism falls back to its own
    decision record, so a label a person accepted before this bookkeeping existed is
    still recognised as theirs and is not quietly replaced by an automatic proposal.
    """
    metadata = _metadata(record)
    authorship = _authorship_record(metadata)
    organism = sample_organism(record)
    if organism["source"] == "assigned" and organism["set_by"]:
        for field in ORGANISM_FIELDS:
            authorship.setdefault(field, {"by": organism["set_by"], "utc": organism["set_utc"],
                                          "surface": organism["surface"],
                                          "automatic": organism["automatic"]})
    return authorship


def sample_configuration(record) -> dict:
    """The one configuration a sample has: organism, scheme choice, typing mode, run flags.

    Every menu and submenu reads this and writes through
    storage.set_sample_configuration, so a change made anywhere is what every other
    surface shows next. Older top-level metadata keys are still read as a fallback
    so a project made before the single home keeps working; they are never written
    back to. The organism is reported with its source, and is not promoted to an
    assignment by being read here.
    """
    metadata = _metadata(record)
    workflow = metadata.get("workflow") if isinstance(metadata.get("workflow"), dict) else {}
    organism = sample_organism(record)

    def choice(field, default=None):
        return workflow[field] if field in workflow else metadata.get(field, default)

    def text(field):
        value = choice(field)
        return str(value).strip() or None if value not in (None, "") else None

    code = choice("genetic_code")
    return {"genus": organism["genus"], "species": organism["species"],
            "organism": organism["organism"], "organism_source": organism["source"],
            "organism_basis": organism["basis"], "organism_confidence": organism["confidence"],
            "organism_status": organism["status"], "organism_note": organism["note"],
            "typing_mode": str(choice("typing_mode", "") or ""),
            "scheme_path": text("scheme_path"),
            "calling_mode": str(choice("calling_mode", "auto") or "auto"),
            "genetic_code": int(code) if (isinstance(code, int) and not isinstance(code, bool))
            or (isinstance(code, str) and code.strip().isdigit()) else None,
            "cgmlst_scheme_key": text("cgmlst_scheme_key"),
            "cgmlst_scheme_path": text("cgmlst_scheme_path"),
            "run_hydra": bool(choice("run_hydra", False)),
            "run_cgmlst": bool(choice("run_cgmlst", False)),
            "set_by": configuration_authorship(record)}


def _scheme_version(state) -> str | None:
    if state.get("scheme_version"):
        return str(state["scheme_version"])
    metadata = state.get("scheme_metadata")
    if isinstance(metadata, dict):
        for key in _SCHEME_VERSION_KEYS:
            if metadata.get(key):
                return str(metadata[key])
    return None


def typing_state(value) -> dict:
    """What one typing run was made against, read out of a stored result or a request.

    The caller settings are whatever the caller recorded for itself; a request
    should state only the settings that change the answer, because every setting it
    states is one the stored result has to be able to answer for.
    """
    value = value if isinstance(value, dict) else {}
    caller = value.get("caller")
    if not isinstance(caller, dict):
        caller = value.get("parameters") if isinstance(value.get("parameters"), dict) else {}
    return {"input_sha256": _sha256(value.get("input_sha256")),
            "scheme_digest": str(value.get("scheme_digest") or "") or None,
            "scheme": str(value.get("scheme") or "") or None,
            "scheme_version": _scheme_version(value), "caller": dict(caller)}


def _rerun_compare(field, label, before, now, decision) -> None:
    if before is None or now is None:
        missing = "the request does not state it" if now is None else "the stored result does not record it"
        decision["unverified"].append({"field": field, "reason": f"{label} cannot be compared: {missing}."})
        return
    decision["compared"].append(field)
    if before != now:
        decision["changes"].append({"field": field, "was": before, "now": now,
                                    "reason": f"{label} changed."})


def rerun_decision(stored, requested) -> dict:
    """Whether typing has to run again, and exactly what changed since it last ran.

    Pure: it reads two mappings, touches no database and no file, so a Run button
    can call it to enable or explain itself. Typing is repeated only when something
    it was run against changed -- the input bytes, the scheme's identity, the
    scheme's version, or a caller setting. When nothing changed, the stored result
    stands and may be said to stand.

    "Cannot be compared" is a third answer, kept apart from both: a result whose
    input hash nobody recorded is not shown to be current, so a rerun is still asked
    for and the reason names what is missing rather than implying staleness.
    """
    before, now = typing_state(stored), typing_state(requested)
    decision = {"rerun": True, "status": "never_run", "changes": [], "unverified": [],
                "compared": [], "summary": ""}
    if not before["scheme_digest"]:
        decision["summary"] = ("Nothing is stored for this sample against a scheme, so typing has "
                               "not been run yet.")
        return decision
    for field in ("input_sha256", "scheme_digest"):
        _rerun_compare(field, _RERUN_LABELS[field], before[field], now[field], decision)
    for field in ("scheme", "scheme_version"):
        # Optional: a request that does not name the scheme version is not asking
        # about it, and is not told that it could not be checked.
        if now[field] is not None:
            _rerun_compare(field, _RERUN_LABELS[field], before[field], now[field], decision)
    for key in sorted(now["caller"]):
        _rerun_compare(f"caller.{key}", f"the caller setting {key!r}",
                       before["caller"].get(key), now["caller"][key], decision)
    if decision["changes"]:
        decision["status"] = "changed"
        decision["summary"] = ("Typing has to run again: "
                               + " ".join(change["reason"] for change in decision["changes"]))
    elif decision["unverified"]:
        decision["status"] = "unverifiable"
        decision["summary"] = ("The stored result cannot be shown to be current. "
                               + " ".join(entry["reason"] for entry in decision["unverified"]))
    else:
        decision.update(rerun=False, status="current")
        decision["summary"] = ("Nothing typing was run against has changed, so the stored result "
                               "stands. Checked: " + ", ".join(decision["compared"]) + ".")
    return decision


def _cohort_label(record) -> str:
    result = record.get("result") if isinstance(record.get("result"), dict) else {}
    return str(record.get("name") or result.get("sample_name") or record.get("sample_name")
               or record.get("id") or record.get("sample_id") or "this sample")


def _name_list(entries, limit=5) -> str:
    names = [f"{entry['name']} ({entry['organism']})" if entry["organism"] else entry["name"]
             for entry in entries[:limit]]
    remaining = len(entries) - len(names)
    return ", ".join(names) + (f", and {remaining} more" if remaining > 0 else "")


def cgmlst_scheme_applies(records, scheme) -> dict:
    """Whether one cgMLST scheme applies to the organisms in a cohort, and why it does not.

    A cgMLST target set is defined for one organism or one species complex. Calling
    a genome of another genus against it does not produce a larger distance, it
    produces absent targets, and absent targets are not differences. So a scheme is
    refused for a cohort that contradicts it, and refused for a cohort in which no
    sample's organism is known yet -- an unidentified isolate is not a match.

    A sample of the same genus but another species is allowed and named in
    ``warnings``: several catalogued schemes are complex-level on purpose, and which
    species a target set was defined on is something the person launching the run
    has to be told rather than have decided for them.
    """
    from wmlstudio.project import CGMLST_LOCUS_FLOOR
    scheme = scheme if isinstance(scheme, dict) else {}
    genus = str(scheme.get("genus") or "").strip()
    species = str(scheme.get("species") or "").strip()
    label = (str(scheme.get("organism") or "").strip()
             or " ".join(part for part in (genus, species) if part) or "no organism")
    name = str(scheme.get("scheme_name") or scheme.get("key") or "the chosen cgMLST scheme")
    verdict = {"applies": False, "status": "", "scheme_key": scheme.get("key"), "scheme_name": name,
               "scheme_organism": label, "matched": [], "mismatched": [], "unknown": [],
               "warnings": [], "message": ""}
    count = scheme.get("locus_count")
    if isinstance(count, int) and count <= CGMLST_LOCUS_FLOOR:
        verdict["status"] = "not_a_cgmlst_scheme"
        verdict["message"] = (f"{name} declares {count} targets, at or below the "
                              f"{CGMLST_LOCUS_FLOOR}-locus floor that separates classical MLST "
                              "from gene-by-gene typing. It cannot be run as a cgMLST scheme.")
        return verdict
    if not genus:
        verdict["status"] = "scheme_organism_unknown"
        verdict["message"] = (f"{name} records no organism, so nothing confirms it applies to "
                              "these samples. Choose a catalogued cgMLST scheme.")
        return verdict
    for record in records:
        organism = sample_organism(record)
        entry = {"sample_id": str(record.get("id") or record.get("sample_id") or ""),
                 "name": _cohort_label(record), "organism": organism["organism"]}
        if not organism["genus"]:
            verdict["unknown"].append(entry)
        elif organism["genus"].casefold() != genus.casefold():
            verdict["mismatched"].append(entry)
        else:
            verdict["matched"].append(entry)
            if species and organism["species"] and organism["species"].casefold() != species.casefold():
                verdict["warnings"].append(
                    f"{entry['name']} is {organism['organism']}; {name} was defined on {label}. "
                    "Check the scheme's scope before comparing these profiles.")
    if verdict["mismatched"]:
        verdict["status"] = "organism_mismatch"
        verdict["message"] = (f"{name} is defined for {label}. "
                              f"{len(verdict['mismatched'])} selected sample(s) are a different "
                              f"organism: {_name_list(verdict['mismatched'])}. A cgMLST distance "
                              "never crosses organisms; choose that organism's scheme or deselect "
                              "those samples.")
        return verdict
    if not verdict["matched"]:
        verdict["status"] = "organism_unknown"
        verdict["message"] = (f"No selected sample has an organism yet, so nothing confirms that "
                              f"{name} ({label}) applies to it. Identify or assign the organism "
                              "first, then run cgMLST.")
        return verdict
    verdict["applies"] = True
    verdict["status"] = "applies"
    verdict["message"] = (f"{name} applies to {len(verdict['matched'])} of {len(verdict['matched']) + len(verdict['unknown'])} "
                          f"selected sample(s) as {label}.")
    if verdict["unknown"]:
        verdict["warnings"].append(
            f"{len(verdict['unknown'])} selected sample(s) have no organism yet and were not "
            f"checked against {label}: {_name_list(verdict['unknown'])}.")
    return verdict


def _organism_typing(record):
    """Organism-specific calls as one line, or empty when none are current.

    Imported here because the registry imports its assay modules on first use,
    and a flat export must not pay for that until it asks for the column.
    """
    from wmlstudio.characterization import current_characterization
    from wmlstudio.organism_modules import summarize_record
    return summarize_record(current_characterization(record).get("evidence") or {})


def feature_fields(record, highlight=None) -> dict:
    metadata = record.get("metadata", {})
    metadata = metadata if isinstance(metadata, dict) else {}
    # One reader for the organism, so this row and every menu show the same answer.
    identity = sample_organism(record)
    genus, species, organism = identity["genus"], identity["species"], identity["organism"]
    evidence = current_hydra_evidence(record)
    evidence_state = hydra_evidence_status(record)
    primary = [hit for hit in evidence.get("hits", []) if hit.get("primary") is True]
    def genes(kind):
        return sorted({hit["gene"] for hit in primary if hit.get("element_type") == kind})
    cluster = highlight if highlight is not None else metadata.get("cluster", {})
    cluster = cluster if isinstance(cluster, dict) else {"label": str(cluster), "highlight": True}
    workflow = metadata.get("workflow", {})
    # The decision record says which engine produced a label and how strong it is.
    # It is reported beside organism_source rather than inside it, because the three
    # existing organism_source values are consumed verbatim by report wording.
    # Named apart from the HYDRA evidence bound above: rebinding that name here
    # silently emptied the hydra_* fields in every export.
    decision = metadata.get("organism_evidence")
    decision = decision if isinstance(decision, dict) else {}
    # Classical MLST and cgMLST are reported in their own columns: the headline
    # result is whichever analysis ran last, and an ST must not disappear from a
    # row because a core-genome run happened after it.
    profiles = typing_profiles(record)
    mlst, cgmlst = profiles["mlst"] or {}, profiles["cgmlst"] or {}
    cg_alleles = cgmlst.get("alleles") if isinstance(cgmlst.get("alleles"), dict) else {}
    return {
        "genus": genus, "species": species, "organism": organism,
        "organism_source": "assigned" if identity["source"] == "assigned"
        else "local_scheme_detection" if organism else "unknown",
        "organism_basis": str(decision.get("basis") or ""),
        "organism_confidence": str(decision.get("confidence") or ""),
        "organism_status": str(decision.get("status") or ""),
        "organism_quarantine": str(decision.get("quarantine_reason") or ""),
        # Gated on current_characterization, so typing belonging to an earlier
        # assembly leaves the column empty rather than reading as this one's.
        "organism_typing": _organism_typing(record),
        "typing_mode": workflow.get("typing_mode", metadata.get("typing_mode", "")),
        "typing_kinds": [kind for kind in ("mlst", "cgmlst") if profiles[kind]],
        "mlst_st": mlst.get("st"), "mlst_scheme": mlst.get("scheme") or "",
        "mlst_status": mlst.get("status") or "",
        "cgmlst_scheme": cgmlst.get("scheme") or "", "cgmlst_status": cgmlst.get("status") or "",
        # Scheme size and called loci are separate numbers: a missing call is not
        # an allele, and the denominator belongs beside the numerator.
        "cgmlst_loci": len(cg_alleles), "cgmlst_called": sum(1 for v in cg_alleles.values() if v),
        "amr_genes": genes("AMR"), "virulence_genes": genes("VIRULENCE"),
        "plasmid_replicons": genes("PLASMID"), "stress_genes": genes("STRESS"),
        "amr_classes": sorted({hit["class"] for hit in primary if hit.get("element_type") == "AMR" and hit.get("class")}),
        "hydra_source_sample": evidence.get("source_sample", ""),
        'hydra_evidence_status': evidence_state['status'],
        'hydra_evidence_reason': evidence_state['reason'],
        "hydra_report_sha256": evidence.get("report_sha256", ""),
        "hydra_amr_genes": evidence.get("summary", {}).get("amr_genes") if evidence else None,
        "cluster_label": cluster.get("label", ""),
        "cluster_color": cluster.get("color", "#2F8A78"),
        "cluster_highlight": bool(cluster.get("highlight", False)),
    }


def combined_feature_rows(records, selected_ids=None, highlight_clusters=None) -> list[dict]:
    rows = []
    for record in select_records(records, selected_ids):
        sample_id = str(record.get("id") or record.get("sample_id") or "")
        result = record.get("result") or (record if "result" not in record else {})
        fields = feature_fields(record, (highlight_clusters or {}).get(sample_id))
        rows.append({"sample_id": sample_id, "sample_name": record.get("name", result.get("sample_name", "")),
                     "input_path": record.get("input_path", ""), "job_status": record.get("status", ""),
                     "status": result.get("status", ""), "st": result.get("st"),
                     "scheme": result.get("scheme"), "scheme_digest": result.get("scheme_digest"),
                     **fields})
    return rows
