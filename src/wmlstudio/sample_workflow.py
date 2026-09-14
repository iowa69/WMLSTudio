"""Stable-identity evidence linking and flattened features for native views."""

from __future__ import annotations

import copy
import re
from pathlib import Path

from wmlstudio.sequence import file_sha256, sample_name


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


def feature_fields(record, highlight=None) -> dict:
    result = record.get("result") or (record if "result" not in record else {})
    metadata = record.get("metadata", {})
    metadata = metadata if isinstance(metadata, dict) else {}
    assigned = metadata.get("organism", {})
    if isinstance(assigned, str):
        parts = assigned.split(maxsplit=1)
        assigned = {"genus": parts[0] if parts else "", "species": parts[1] if len(parts) > 1 else ""}
    assigned = assigned if isinstance(assigned, dict) else {}
    detected = result.get("identification", {}) or {}
    detected_organism = detected.get("organism") or {}
    detected_organism = detected_organism if isinstance(detected_organism, dict) else {}
    genus = assigned.get("genus") or detected_organism.get("genus") or detected.get("genus") or ""
    species = assigned.get("species") or detected_organism.get("species") or detected.get("species") or ""
    organism = " ".join(filter(None, [genus, species]))
    evidence = current_hydra_evidence(record)
    evidence_state = hydra_evidence_status(record)
    primary = [hit for hit in evidence.get("hits", []) if hit.get("primary") is True]
    def genes(kind):
        return sorted({hit["gene"] for hit in primary if hit.get("element_type") == kind})
    cluster = highlight if highlight is not None else metadata.get("cluster", {})
    cluster = cluster if isinstance(cluster, dict) else {"label": str(cluster), "highlight": True}
    workflow = metadata.get("workflow", {})
    return {
        "genus": genus, "species": species, "organism": organism,
        "organism_source": "assigned" if assigned.get("genus") else "local_scheme_detection" if organism else "unknown",
        "typing_mode": workflow.get("typing_mode", metadata.get("typing_mode", "")),
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
