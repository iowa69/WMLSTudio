"""Read-only outbreak workflow summaries: what is known, missing, and actionable.

The map does not run analyses or infer a clinical verdict. Counts are derived
from stable sample records; linked read mates are not counted as extra isolates.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from typing import Any, Mapping, Sequence


@dataclass(frozen=True)
class JourneyStep:
    key: str
    title: str
    question: str
    action: str
    action_label: str
    summary: str
    state: str = "ready"


def isolate_records(samples: Sequence[Mapping[str, Any]]) -> list[Mapping[str, Any]]:
    return [sample for sample in samples
            if (sample.get("metadata") or {}).get("workflow", {}).get("source_kind") != "read_mate"]


def characterization_state(sample: Mapping[str, Any]) -> tuple[str, Mapping[str, Any]]:
    """Do not display a previous assembly's characterization as current evidence."""
    metadata = sample.get("metadata") or {}
    evidence = metadata.get("characterization") or {}
    if not isinstance(evidence, dict) or not evidence:
        return "not_run", {}
    from wmlstudio.sample_workflow import current_input_sha256
    expected = current_input_sha256(sample)
    observed = evidence.get("input_sha256")
    if expected and observed and expected != observed:
        return "stale", evidence
    if not observed or not expected:
        return "unverified", evidence
    return str(evidence.get("status", "completed")), evidence


def journey_summary(
    samples: Sequence[Mapping[str, Any]], *, selected_ids: set[str] | None = None,
    comparison_ids: set[str] | None = None,
) -> dict[str, Any]:
    """Summarize saved evidence without opening genomes or mutating project state."""
    from wmlstudio.sample_workflow import hydra_evidence_status

    isolates = isolate_records(samples)
    selected_ids = selected_ids or set()
    counts: Counter[str] = Counter(total=len(isolates), linked_mates=len(samples) - len(isolates))
    review_ids = []
    pending_ids = []
    for sample in isolates:
        result = sample.get("result") or {}
        counts["selected"] += sample["id"] in selected_ids
        counts["comparison"] += comparison_ids is not None and sample["id"] in comparison_ids
        if sample.get("status") in {"queued", "interrupted"}:
            pending_ids.append(sample["id"])
        if sample.get("status") in {"failed", "interrupted"} or result.get("status") in {
            "mixed", "ambiguous", "incomplete", "failed",
        }:
            review_ids.append(sample["id"])
        counts["unavailable"] += bool(sample.get("missing_input"))
        counts["profile_only"] += bool(sample.get("profile_only"))
        counts["typed"] += bool(result.get("st") is not None and result.get("status") == "complete"
                                 and sample.get("status") == "completed")
        counts["profiled"] += bool(result.get("alleles") and sample.get("status") == "completed")
        counts["qc"] += bool(result.get("qc"))
        counts["amr"] += hydra_evidence_status(sample)["status"] == "current"
        state, evidence = characterization_state(sample)
        species = evidence.get("species_evidence") or {}
        counts["characterized"] += state == "completed"
        counts["species_reviewed"] += state == "completed" and species.get("status") in {
            "completed", "resolved", "supported", "identified",
        }
        counts["characterization_stale"] += state == "stale"
    counts.update(pending=len(pending_ids), review=len(review_ids))
    total = counts["total"]
    steps = [
        JourneyStep("inputs", "1  Define your isolates", "Which files belong to each isolate?",
                    "samples", "Review samples && reads",
                    f"{total} isolates · {counts['linked_mates']} linked read mates · "
                    f"{counts['pending']} pending", "ready" if total else "empty"),
        JourneyStep("identity", "2  Check identity & quality", "Are these the organisms I expect?",
                    "characterize", "Review characterization…",
                    f"{counts['characterized']} characterization results · {counts['review']} need review. "
                    "MLST lineage alone is not independent species confirmation.",
                    "review" if counts["review"] or counts["characterization_stale"] else "ready"),
        JourneyStep("typing", "3  Establish the baseline", "What are the STs and allele profiles?",
                    "analyse", "Review typing run…",
                    f"{counts['typed']} registered STs · {counts['profiled']} saved primary profiles. "
                    "Additional cgMLST/wgMLST calls preserve classical MLST."),
        JourneyStep("features", "4  Inspect resistance & virulence", "Which determinants were detected, and how?",
                    "features", "Open linked evidence",
                    f"{counts['amr']} current linked AMR results. Drug associations are not measured "
                    "susceptibility; plasmid markers are hypotheses, not transmission proof."),
        JourneyStep("relatedness", "5  Investigate relatedness", "Who is close, and is the comparison reliable?",
                    "compare", "Open investigation && graph",
                    f"{counts['comparison']} isolates in the explicit comparison cohort. "
                    "Choose one scheme snapshot; inspect shared loci and cluster chaining."),
        JourneyStep("report", "6  Freeze & communicate", "Which findings belong in this report?",
                    "reports", "Choose report cohort",
                    "Select isolates or reviewed clusters, preserve the snapshot, and add next week's "
                    "samples without discarding earlier evidence."),
    ]
    return {"counts": dict(counts), "review_ids": review_ids, "pending_ids": pending_ids, "steps": steps}
