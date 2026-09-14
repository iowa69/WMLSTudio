"""Layered organism evidence for one sequence file, with its basis retained.

Neither installed engine is a species diagnostic. Whole-genome ANI compares an
assembly against whatever reference panel is installed, and MLST discovery
reports compatibility with an installed allele panel — its own notes say so.
This module therefore never produces "the species". It produces a verdict with a
basis, a confidence tier, the runner-up margin and the engine's own sentences,
and when nothing supports a label it says so instead of guessing.

"No suitable match" means unresolved. It is not a novel species, it is not a
contaminated isolate, and it is not a reason to file a file into the organism
folder someone expected. Unresolved inputs go to the needs-review tree.
"""

from __future__ import annotations

from pathlib import Path

from . import identification, species_evidence
from .scheduler import resources_for_run, run_bounded
from .sequence import AnalysisCancelled, SequenceReader, check_cancelled, file_sha256

EVIDENCE_FORMAT_VERSION = 1
# Weakest to strongest. Ordering is a filing policy, not a statistical scale.
CONFIDENCE_ORDER = ("unresolved", "panel_compatibility", "genus_only",
                    "complex_only", "genomic_reference_supported")
AUTO_FILE_MIN = "genomic_reference_supported"
# Nothing weaker than this may be accepted without a person, whatever the
# configured floor says: an MLST panel match is explicitly not a species call.
AUTO_CONFIRM_FLOOR = "genus_only"
SETTING_AUTO_CONFIRM = "filing.auto_confirm"
SETTING_MIN_CONFIDENCE = "filing.min_confidence"
DEFAULT_POLICY = {"auto_confirm": False, "min_confidence": AUTO_FILE_MIN}
# Plain language first; the technical token belongs in a tooltip, never the cell.
CONFIDENCE_LABELS = {
    "genomic_reference_supported": ("Strong", "Whole-genome comparison met every configured gate against the installed reference panel."),
    "complex_only": ("Review — complex only", "The nearest reference is inside a group this comparison cannot separate."),
    "genus_only": ("Review — genus only", "The nearest reference carries a genus label and no species name."),
    "panel_compatibility": ("Panel match only", "An installed typing panel matched. That is panel compatibility, not a species confirmation."),
    "unresolved": ("Not identified", "Nothing installed supported a label. This is a limit of the installed references, not evidence of a new organism."),
}
_LABELLED = {"genomic_reference_supported", "complex_only", "genus_only"}
_ANI_DETAIL_KEYS = ("reference_id", "genus", "species", "ani", "query_fraction", "reference_fraction")
_MLST_DETAIL_KEYS = ("scheme", "scheme_path", "coverage", "matched_loci", "st", "organism_label")
_MAX_DETAIL = 5


def confidence_rank(token) -> int:
    """Rank a confidence token; an unknown token ranks as the weakest, never as strong."""
    try:
        return CONFIDENCE_ORDER.index(str(token))
    except ValueError:
        return 0


def _taxon(row) -> dict | None:
    if not isinstance(row, dict):
        return None
    return {"reference_id": row.get("reference_id"), "genus": row.get("genus", ""),
            "species": row.get("species", ""), "ani": row.get("ani")}


def _ani_detail(result: dict) -> dict:
    return {"ani_top": [{key: hit.get(key) for key in _ANI_DETAIL_KEYS}
                        for hit in (result.get("hits") or [])[:_MAX_DETAIL]],
            "gap_ani": result.get("gap_ani"), "thresholds": result.get("thresholds") or {}}


def _mlst_detail(result: dict) -> dict:
    return {"mlst_top": [{key: item.get(key) for key in _MLST_DETAIL_KEYS}
                         for item in (result.get("candidates") or [])[:_MAX_DETAIL]],
            "identification_status": result.get("identification_status"),
            "best_scheme_path": result.get("best_scheme_path")}


def _species_quarantine(result: dict) -> str:
    """Name the reason from the engine's own numbers, never by reading its prose."""
    thresholds = result.get("thresholds") or {}
    gap, minimum = result.get("gap_ani"), thresholds.get("min_species_gap_ani_pct")
    if result.get("discordant_low_fraction_hits"):
        return "conflicting_evidence"
    if result.get("runner_up") and gap is not None and minimum is not None and gap < minimum:
        return "conflicting_evidence"
    return "not_in_reference_panel"


def _blank(path: Path, digest: str = "", kind: str = "") -> dict:
    return {"format_version": EVIDENCE_FORMAT_VERSION, "input_path": str(path),
            "input_sha256": digest, "kind": kind, "basis": "none", "confidence": "unresolved",
            "status": "quarantined", "proposed": {"genus": "", "species": ""},
            "accepted": {"genus": "", "species": ""},
            "quarantine_reason": "awaiting_identification", "confirmed_by": None,
            "confirmed_utc": None, "panel": {}, "detail": {}, "margin_ani": None,
            "runner_up": None, "reason": "", "notes": [], "limitations": [], "errors": []}


def _adopt_species(verdict: dict, result: dict, basis: str, root: Path, kind: str) -> None:
    """Copy an ANI outcome into the verdict without strengthening a single word."""
    species = result.get("species", "")
    if result.get("confidence") == "genus_only":
        species = ""  # A genus-level reference never licenses a species name.
    verdict.update(basis=basis, confidence=result["confidence"], status="proposed",
                   proposed={"genus": result.get("genus", ""), "species": species},
                   quarantine_reason=None, reason=result.get("reason", ""),
                   notes=list(result.get("notes") or []),
                   limitations=list(result.get("limitations") or []),
                   detail=_ani_detail(result), margin_ani=result.get("gap_ani"),
                   runner_up=_taxon(result.get("runner_up")),
                   panel={"kind": kind, "path": str(root),
                          "reference_digest": result.get("reference_digest"),
                          "reference_count": result.get("reference_count")})


def _adopt_panel(verdict: dict, result: dict, scheme_paths) -> None:
    """Copy an MLST outcome in, capped at panel compatibility however good the match."""
    organism = result.get("organism") or {}
    verdict.update(basis="mlst_panel", confidence="panel_compatibility", status="proposed",
                   proposed={"genus": organism.get("genus", ""), "species": organism.get("species", "")},
                   quarantine_reason=None, notes=list(result.get("notes") or []),
                   limitations=[], detail=_mlst_detail(result), margin_ani=None, runner_up=None,
                   panel={"kind": "mlst_schemes", "path": result.get("best_scheme_path"),
                          "reference_digest": None, "reference_count": len(list(scheme_paths))})
    verdict["reason"] = verdict["notes"][0] if verdict["notes"] else ""


def _query_panel(path, root, cancelled, progress):
    try:
        return species_evidence.identify_species(path, root, cancelled, progress), ""
    except AnalysisCancelled:
        raise
    except (OSError, ValueError) as error:
        return None, str(error)


def identify_input(path, *, species_panel_root=None, kpsc_panel_root=None, scheme_paths=(),
                   cancelled=None, progress=None) -> dict:
    """Layered organism evidence for one FASTA or FASTQ, run on the user's own file.

    Tried in order, first label wins: whole-genome ANI against the installed
    species panel; a second ANI query against a focused complex panel when the
    first says Klebsiella; then MLST panel compatibility, which can never exceed
    'panel_compatibility' however complete the allele match is. An input this
    cannot identify returns an 'unresolved' verdict with a quarantine reason; it
    never raises for that, and it never invents a genus. Cancellation propagates.
    """
    check_cancelled(cancelled)
    path = Path(path).expanduser().resolve()
    scheme_paths = list(scheme_paths)
    verdict = _blank(path)
    try:
        digest = file_sha256(path, cancelled)
        with SequenceReader(path, cancelled) as reader:
            kind = reader.kind
    except AnalysisCancelled:
        raise
    except (OSError, ValueError) as error:
        verdict["errors"].append({"tier": "input", "message": str(error)})
        verdict["reason"] = str(error)
        return verdict
    verdict.update(input_sha256=digest, kind=kind)
    if kind != "fasta":
        verdict.update(quarantine_reason="reads_not_assembled",
                       reason="Reads were not identified because both installed engines "
                              "require an assembled FASTA. Assemble this sample first.")
        return verdict

    panels = []
    if species_panel_root:
        panels.append((Path(species_panel_root).expanduser().resolve(), "species_panel", "genomic_ani"))
    if kpsc_panel_root:
        focused = Path(kpsc_panel_root).expanduser().resolve()
        if all(focused != existing for existing, _, _ in panels):
            panels.append((focused, "characterization_starter", "genomic_ani_kpsc"))
    species = None
    for root, kind_token, basis in panels:
        check_cancelled(cancelled)
        if verdict["basis"] != "none" and verdict["proposed"]["genus"] != "Klebsiella":
            break  # A second panel refines a complex or fills a gap, never retries a call.
        result, error = _query_panel(path, root, cancelled, progress)
        if error:
            verdict["errors"].append({"tier": basis, "message": error})
            continue
        if result.get("confidence") in _LABELLED and (
                species is None
                or confidence_rank(result["confidence"]) >= confidence_rank(verdict["confidence"])):
            _adopt_species(verdict, result, basis, root, kind_token)
        species = species or result

    if verdict["basis"] == "none" and scheme_paths:
        check_cancelled(cancelled)
        try:
            panel = identification.identify_assembly(path, scheme_paths, cancelled, progress)
        except AnalysisCancelled:
            raise
        except (OSError, ValueError) as error:
            verdict["errors"].append({"tier": "mlst_panel", "message": str(error)})
        else:
            organism = panel.get("organism") or {}
            if panel.get("identification_status") == "assigned" and organism.get("genus"):
                _adopt_panel(verdict, panel, scheme_paths)
            else:
                verdict["detail"] = {**verdict["detail"], **_mlst_detail(panel)}
                verdict["notes"] = list(panel.get("notes") or [])

    if verdict["basis"] == "none":
        if species is not None:
            verdict.update(quarantine_reason=_species_quarantine(species),
                           reason=species.get("reason", ""),
                           limitations=list(species.get("limitations") or []),
                           detail={**verdict["detail"], **_ani_detail(species)},
                           margin_ani=species.get("gap_ani"), runner_up=_taxon(species.get("runner_up")))
        elif not panels:
            verdict["reason"] = ("No species reference panel is installed, so no genomic "
                                 "comparison was attempted.")
    return verdict


def identify_batch(paths, *, allocation=None, cancelled=None, progress=None, on_result=None,
                   **options) -> list[dict]:
    """Identify a batch of original files, returning one verdict per input in order.

    This runs before anything is copied, so a rejected or cancelled identification
    leaves the project exactly as it was and no file is filed twice.
    """
    items = [Path(value).expanduser().resolve() for value in paths]
    if not items:
        return []

    def operation(item, plan, is_cancelled, report):
        return identify_input(item, cancelled=is_cancelled, progress=report, **options)

    completed = run_bounded(items, operation,
                            allocation or resources_for_run({"threads": 2, "memory_gb": 2}),
                            cancelled=cancelled, progress=progress, on_result=on_result)
    verdicts = {item: result for item, result in completed}
    return [verdicts[item] for item in items if item in verdicts]


def resolve_policy(source=None) -> dict:
    """Read the filing policy from a project or a mapping; default to review-everything."""
    if source is None:
        values = {}
    elif hasattr(source, "get_setting"):
        values = {"auto_confirm": source.get_setting(SETTING_AUTO_CONFIRM, DEFAULT_POLICY["auto_confirm"]),
                  "min_confidence": source.get_setting(SETTING_MIN_CONFIDENCE, DEFAULT_POLICY["min_confidence"])}
    else:
        values = dict(source)
    minimum = str(values.get("min_confidence") or DEFAULT_POLICY["min_confidence"])
    if minimum not in CONFIDENCE_ORDER:
        raise ValueError(f"Unknown filing confidence floor: {minimum}")
    return {"auto_confirm": bool(values.get("auto_confirm", DEFAULT_POLICY["auto_confirm"])),
            "min_confidence": minimum}


def meets_policy(verdict: dict, policy=None) -> bool:
    """Whether this verdict may be filed without a person looking at it."""
    policy = resolve_policy(policy)
    if not verdict.get("proposed", {}).get("genus"):
        return False
    rank = confidence_rank(verdict.get("confidence"))
    return bool(policy["auto_confirm"] and rank >= confidence_rank(policy["min_confidence"])
                and rank >= confidence_rank(AUTO_CONFIRM_FLOOR))


def apply_policy(verdict: dict, policy=None) -> dict:
    """Confirm a verdict only where the policy allows; otherwise leave it proposed."""
    if verdict.get("status") in {"confirmed", "quarantined"}:
        return dict(verdict)
    updated = dict(verdict)
    if meets_policy(verdict, policy):
        updated.update(status="confirmed", confirmed_by="auto_policy",
                       accepted=dict(verdict["proposed"]))
    return updated


def proposed_destination(verdict: dict, *, policy=None) -> tuple[str | None, str, str]:
    """Where this verdict files right now: (quarantine bucket or None, genus, species).

    A proposal nobody has accepted goes to the needs-review tree, because the
    folder a file sits in must reflect a decision somebody actually made.
    """
    verdict = apply_policy(verdict, policy)
    if verdict["status"] == "quarantined":
        return verdict.get("quarantine_reason") or "awaiting_identification", "", ""
    if verdict["status"] == "confirmed":
        accepted = verdict.get("accepted") or {}
        if accepted.get("genus"):
            return None, accepted["genus"], accepted.get("species", "")
        return "user_deferred", "", ""
    return "awaiting_identification", "", ""


def assignment_for(verdict: dict, *, policy=None, name=None, scheme_path=None) -> dict:
    """Turn a verdict into the assignment dict storage.import_samples consumes."""
    quarantine, genus, species = proposed_destination(verdict, policy=policy)
    assignment = {"path": verdict["input_path"], "genus": genus, "species": species,
                  "typing_mode": "manual" if genus else "auto", "quarantine": quarantine,
                  "organism_evidence": apply_policy(verdict, policy)}
    if name:
        assignment["name"] = name
    if scheme_path:
        assignment["scheme_path"] = scheme_path
    return assignment


def evidence_sentences(verdict: dict) -> list[str]:
    """Every sentence the engines wrote about this verdict, verbatim and in order."""
    sentences, seen = [], set()
    for value in [verdict.get("reason") or "", *(verdict.get("notes") or []),
                  *(verdict.get("limitations") or [])]:
        if value and value not in seen:
            seen.add(value)
            sentences.append(value)
    return sentences


def confidence_label(verdict_or_token) -> tuple[str, str]:
    """Plain-language wording first, with the technical explanation for a tooltip."""
    token = (verdict_or_token.get("confidence") if isinstance(verdict_or_token, dict)
             else verdict_or_token)
    return CONFIDENCE_LABELS.get(str(token), CONFIDENCE_LABELS["unresolved"])
