"""Dated, reviewable publication evidence, never an organism-wide clinical rule.

The catalog is deliberately separate from allele databases and contains no
licensed scheme contents. A matching taxon or locus count does not validate a
caller, missing-data policy, or clustering threshold.
"""

from __future__ import annotations

import hashlib
import json
from copy import deepcopy
from datetime import date, datetime, timezone

CATALOG_VERSION = "2026-09-12.1"
REVIEWED_ON = "2026-09-12"
INTERPRETATION = (
    "A genomic cluster is a research signal for epidemiological review, not proof of direct transmission, "
    "direction of transmission, or clinical causation. Single linkage can join distant endpoints through "
    "intermediate isolates. Missing calls are unknown, not identical. A published cutoff is not validated "
    "for WMLSTudio merely because the organism or locus count matches."
)

SOURCES = {
    "glasgow2025": {"citation": "Glasgow et al. (2025). Comparison of core genome multi-locus sequencing typing pipelines for hospital outbreak detection of common bacterial pathogens.",
        "doi": "10.1128/jcm.00646-25", "url": "https://journals.asm.org/doi/10.1128/jcm.00646-25", "published": "2025-08-27",
        "locator": "Materials and methods; Supplemental Tables S1–S2", "type": "peer_reviewed_pipeline_comparison"},
    "siddall2025": {"citation": "Siddall, Starkey and Patel (2025). Automated whole genome sequencing platform for bacterial strain typing in clinical microbiology laboratories.",
        "doi": "10.1128/jcm.00178-25", "url": "https://journals.asm.org/doi/10.1128/jcm.00178-25", "published": "2025-04-22",
        "locator": "Table 1 and whole genome sequencing analysis methods", "type": "peer_reviewed_local_validation"},
    "kampmeier2022": {"citation": "Kampmeier et al. (2022). Development and Evaluation of a Core Genome Multilocus Sequencing Typing (cgMLST) Scheme for Serratia marcescens Molecular Surveillance and Outbreak Investigations.",
        "doi": "10.1128/jcm.01196-22", "url": "https://journals.asm.org/doi/10.1128/jcm.01196-22", "published": "2022-10-10",
        "locator": "Results and discussion; Figures 1–4", "type": "peer_reviewed_scheme_evaluation"},
    "debeen2015": {"citation": "de Been et al. (2015). Core Genome Multilocus Sequence Typing Scheme for High-Resolution Typing of Enterococcus faecium.",
        "doi": "10.1128/jcm.01946-15", "url": "https://journals.asm.org/doi/10.1128/jcm.01946-15", "published": "2015-11-18",
        "locator": "Results: definition and evaluation of clonally related isolates", "type": "peer_reviewed_scheme_evaluation"},
    "higgs2022": {"citation": "Higgs et al. (2022). Optimising genomic approaches for identifying vancomycin-resistant Enterococcus faecium transmission in healthcare settings.",
        "doi": "10.1038/s41467-022-28156-4", "url": "https://www.nature.com/articles/s41467-022-28156-4", "published": "2022-01-26",
        "locator": "Results; Figures 1, 2 and 4; Methods", "type": "peer_reviewed_genomic_epidemiology"},
    "tonnies2021": {"citation": "Tönnies et al. (2021). Establishment and Evaluation of a Core Genome Multilocus Sequence Typing Scheme for Whole-Genome Sequence-Based Typing of Pseudomonas aeruginosa.",
        "doi": "10.1128/jcm.01987-20", "url": "https://journals.asm.org/doi/10.1128/jcm.01987-20", "published": "2021-02-18",
        "locator": "Scheme definition and outbreak evaluation", "type": "peer_reviewed_scheme_evaluation"},
    "citrobacter2025": {"citation": "Kieninger et al. (2025). Development and validation of a core genome multilocus sequence typing scheme for Citrobacter freundii: application in outbreak investigations and comparative analysis across the Citrobacter genus.",
        "doi": "10.1128/jcm.00860-25", "url": "https://journals.asm.org/doi/10.1128/jcm.00860-25", "published": "2025-09-19",
        "locator": "Results: species-specific scheme and outbreak evaluation; Figure 4", "type": "peer_reviewed_scheme_evaluation"},
    "moura2016": {"citation": "Moura, Criscuolo, Pouseele, Maury, Leclercq, Tarr, Björkman, Dallman, Reimer, Enouf, Larsonneur, Carleton, Bracq-Dieye, Katz, Jones, Touchon, Tourdjman, Walker, Stroika, Cantinelli, Chenal-Francisque, Kucerova, Rocha, Nadon, Grant, Nielsen, Pot, Gerner-Smidt, Lecuit and Brisse (2016). Whole genome-based population biology and epidemiological surveillance of Listeria monocytogenes.",
        "doi": "10.1038/nmicrobiol.2016.185", "url": "https://www.nature.com/articles/nmicrobiol2016185", "published": "2016-10-10",
        "locator": "Results: nomenclature of Lm cgMLST profiles; Figure 1a; Supplementary Sections 2.1 and 2.7", "type": "peer_reviewed_scheme_definition"},
    "vanwalle2018": {"citation": "Van Walle, Björkman, Cormican, Dallman, Mossong, Moura, Pietzka, Ruppitsch, Takkinen and the European Listeria WGS typing group (2018). Retrospective validation of whole genome sequencing-enhanced surveillance of listeriosis in Europe, 2010 to 2015.",
        "doi": "10.2807/1560-7917.ES.2018.23.33.1700798", "url": "https://www.eurosurveillance.org/content/10.2807/1560-7917.ES.2018.23.33.1700798", "published": "2018-08-16",
        "locator": "Results: epidemiological validation; Figure 2b; Table 3", "type": "peer_reviewed_multi_country_evaluation"},
    "bletz2018": {"citation": "Bletz, Janezic, Harmsen, Rupnik and Mellmann (2018). Defining and Evaluating a Core Genome Multilocus Sequence Typing Scheme for Genome-Wide Typing of Clostridium difficile.",
        "doi": "10.1128/JCM.01987-17", "url": "https://journals.asm.org/doi/10.1128/JCM.01987-17", "published": "2018-04-04",
        "locator": "Results: outbreak evaluation preceding Figure 2; Figure 2 legend; Discussion", "type": "peer_reviewed_scheme_evaluation"},
    "kohl2018": {"citation": "Kohl, Harmsen, Rothgänger, Walker, Diel and Niemann (2018). Harmonized Genome Wide Typing of Tubercle Bacilli Using a Web-Based Gene-By-Gene Nomenclature System.",
        "doi": "10.1016/j.ebiom.2018.07.030", "url": "https://www.sciencedirect.com/science/article/pii/S2352396418302779", "published": "2018-08-11",
        "locator": "Research in Context; Methods 2.5; Results and Figure 2a; Table 3", "type": "peer_reviewed_scheme_evaluation"},
    "leeper2023": {"citation": "Leeper, Tolar, Griswold, Vidyaprakash, Hise, Williams, Im, Chen, Pouseele and Carleton (2023). Evaluation of whole and core genome multilocus sequence typing allele schemes for Salmonella enterica outbreak detection in a national surveillance network, PulseNet USA.",
        "doi": "10.3389/fmicb.2023.1254777", "url": "https://www.frontiersin.org/journals/microbiology/articles/10.3389/fmicb.2023.1254777/full", "published": "2023-10-19",
        "locator": "Introduction: PulseNet cluster definition and scheme locus counts; Table 1; Discussion limitations", "type": "peer_reviewed_surveillance_evaluation"},
    "morangilad2015": {"citation": "Moran-Gilad, Prior, Yakunin, Harrison, Underwood, Lazarovitch, Valinsky, Lück, Krux, Agmon, Grotto and Harmsen (2015). Design and application of a core genome multilocus sequence typing scheme for investigation of Legionnaires' disease incidents.",
        "doi": "10.2807/1560-7917.ES2015.20.28.21186", "url": "https://www.eurosurveillance.org/content/10.2807/1560-7917.ES2015.20.28.21186", "published": "2015-07-16",
        "locator": "Results: humidifier-associated incident and cluster-type calibration; Discussion, opening paragraph", "type": "peer_reviewed_scheme_definition"},
}

# This is a surveillance work-list, not a ranking of danger or prevalence.
ORGANISMS = (
    "Acinetobacter baumannii", "Escherichia coli", "Enterococcus faecium", "Enterococcus faecalis",
    "Klebsiella pneumoniae", "Klebsiella variicola", "Klebsiella quasipneumoniae", "Klebsiella oxytoca complex",
    "Klebsiella aerogenes", "Enterobacter cloacae complex", "Citrobacter freundii", "Serratia marcescens",
    "Pseudomonas aeruginosa", "Stenotrophomonas maltophilia", "Burkholderia cepacia complex",
    "Proteus mirabilis", "Morganella morganii", "Providencia stuartii", "Staphylococcus aureus",
    "Staphylococcus epidermidis", "Staphylococcus lugdunensis", "Staphylococcus capitis",
    "Streptococcus pneumoniae", "Streptococcus pyogenes", "Streptococcus agalactiae",
    "Clostridioides difficile", "Legionella pneumophila", "Cutibacterium acnes", "Salmonella enterica",
    "Listeria monocytogenes", "Mycobacterium tuberculosis complex",
)

# Fields are facts extracted from the cited study, not default application values.
# Explicit None means the exact target set is not established by this entry.
_STANDARD = (
    ("Acinetobacter baumannii", 9, 2390, "cgmlst.org:abaumannii-2390"),
    ("Escherichia coli", 10, 2513, "cgmlst.org:ecoli-2513"),
    ("Enterococcus faecalis", 7, 1972, "cgmlst.org:efaecalis-1972"),
    ("Enterococcus faecium", 20, 1423, "cgmlst.org:efaecium-1423"),
    ("Klebsiella pneumoniae", 15, 2358, "cgmlst.org:kpneumoniae-2358"),
    ("Staphylococcus aureus", 24, 1861, "cgmlst.org:saureus-1861"),
)
_MAYO = (
    ("Staphylococcus aureus", 8), ("Acinetobacter baumannii", 9), ("Klebsiella pneumoniae", 15),
    ("Legionella pneumophila", 4), ("Clostridioides difficile", 6), ("Escherichia coli", 10),
    ("Enterobacter cloacae complex", 15), ("Enterococcus faecium", 7), ("Enterococcus faecalis", 7),
    ("Streptococcus pyogenes", 20), ("Serratia marcescens", 12), ("Pseudomonas aeruginosa", 6),
    ("Streptococcus agalactiae", 20), ("Staphylococcus lugdunensis", 8),
    ("Staphylococcus epidermidis", 8), ("Cutibacterium acnes", 5),
)


def _entry(organism, value, source, *, suffix="", scheme_key=None, loci=None,
           method="cgmlst", scope="", missing_policy="Protocol-specific; review original methods.", note=""):
    return {"id": organism.lower().replace(" ", "-") + "-" + source + suffix,
            "organism": organism, "method": method, "published_threshold": value,
            "operator": "<=", "unit": "allele differences" if method == "cgmlst" else "SNPs",
            "scheme_key": scheme_key, "locus_count": loci, "source_id": source,
            "scope": scope, "missing_policy": missing_policy, "limitations": note,
            "auto_apply": False, "reviewed_on": REVIEWED_ON}


def catalog_entries():
    entries = [_entry(organism, value, "glasgow2025", scheme_key=key, loci=loci,
        scope="SeqSphere+ public-scheme comparison in a selected hospital isolate cohort; not WMLSTudio calibration.",
        note="Selection was based on prior clustering; pipeline agreement is not an independent proof of transmission. Exact snapshot/caller settings still require review.")
        for organism, value, loci, key in _STANDARD]
    entries.extend(_entry(organism, value, "siddall2025", suffix="-local",
        scope="Mayo Clinic local related-isolate criterion; SeqSphere+ 10.0.5 / SKESA 2.3.0.",
        note="Exact scheme identity is not curated here. A local validation criterion is not transferable merely by species name. See study for its intermediate category.")
        for organism, value in _MAYO)
    entries.append(_entry("Serratia marcescens", 12, "kampmeier2022", loci=2692,
        scheme_key="cgmlst.org:smarcescens-2692", scope="Db11-derived 2,692-locus scheme; outbreak evaluation.",
        missing_policy="Pairwise ignore missing targets.", note="Prompts epidemiological investigation; incomplete epidemiological information in source. Not a direct-transmission diagnosis."))
    entries.append(_entry("Enterococcus faecium", 20, "debeen2015", loci=1423,
        scheme_key="cgmlst.org:efaecium-1423", scope="de Been 1,423-locus scheme; relatedness within the evaluated outbreaks.",
        note="Later work uses different cutoffs for different questions. Keep study context and isolate QC with the threshold."))
    entries.append(_entry("Enterococcus faecium", 25, "higgs2022", suffix="-population", loci=1423,
        scheme_key="cgmlst.org:efaecium-1423", scope="COREugate / de Been scheme; single-linkage population grouping before fine-scale analysis.",
        note="This is a population grouping cutoff, not the study's fine-scale transmission cutoff."))
    entries.append(_entry("Enterococcus faecium", 7, "higgs2022", suffix="-ska", method="snp",
        scheme_key="higgs2022:ska-short-reads-k15", scope="Original SKA short-read protocol, k=15, within-patient calibrated single linkage.",
        note="Not validated for SKA2 assembly distance, reference-mapped SNPs, coding-only SNPs or changed masks. No automatic transfer."))
    entries.append(_entry("Pseudomonas aeruginosa", None, "tonnies2021", loci=3867,
        scheme_key="cgmlst.org:paeruginosa-3867", scope="3,867-locus public scheme; do not substitute the 2025 comparison's ad hoc scheme.",
        note="No universal numeric rule curated from this source. Review lineage, hypermutation, environmental persistence and the exact protocol."))
    entries.append(_entry("Citrobacter freundii", 10, "citrobacter2025", loci=3250,
        scheme_key="kieninger2025:cfreundii-3250", scope="Species-specific 3,250-target scheme; maximum intracluster distance in evaluated outbreaks.",
        note="Not the paper's separate combined-species scheme or its 8-allele criterion. Single-linkage endpoints may exceed the published within-cluster bound; review cluster diameters."))
    entries.append(_entry("Listeria monocytogenes", 7, "moura2016", loci=1748,
        scheme_key="pasteur:lmonocytogenes-1748", scope="Institut Pasteur BIGSdb-Lm 1,748-locus scheme; cgMLST type (CT) definition by single linkage.",
        missing_policy="Mismatches counted among loci called in both profiles.",
        note="A CT is a surveillance grouping, not a transmission finding. The authors state CTs diversify slowly (about 0.2 alleles per year), so short-term transmission may need finer methods. The separate 150-mismatch sublineage cutoff is not an outbreak cutoff."))
    entries.append(_entry("Listeria monocytogenes", 7, "vanwalle2018", suffix="-ruppitsch", loci=1701,
        scheme_key="cgmlst.org:lmonocytogenes-1701", scope="ECDC-led retrospective evaluation across 27 EU/EEA countries and 19 confirmed outbreaks; Ruppitsch 1,701-locus scheme in SeqSphere+.",
        note="Confirmed useful for cluster detection at 7 allele differences, with positive predictive value near 58-68%: most detected clusters were not confirmed outbreaks. The authors offer 4 allele differences as a stricter option for more compelling microbiological evidence. Case definitions must state the scheme, its full locus count and the cutoff."))
    entries.append(_entry("Clostridioides difficile", 6, "bletz2018", loci=2270,
        scheme_key="cgmlst.org:cdifficile-2270", scope="Bletz 2,270-target scheme seeded on strain 630 (NC_009089.1); cluster type grouping on a minimum spanning tree.",
        missing_policy="Pairwise ignore missing targets.",
        note="The authors observed at most 3 allele differences among linked isolates and doubled it to 6 as a precaution, so this is precautionary rather than statistically calibrated. SCHEME VERSION MATTERS: the provider now lists this 2,270-target scheme as deprecated v1, with a current v2 of 2,147 targets. Bind the version you actually ran. A separate local 6-allele criterion exists under siddall2025 with no bound scheme; the two must not be read as independent confirmation."))
    entries.append(_entry("Mycobacterium tuberculosis complex", 5, "kohl2018", loci=2891,
        scheme_key="cgmlst.org:mtbc-2891", scope="MTBC 2,891-target scheme seeded on H37Rv (NC_000962.3), PE/PPE and repetitive genes excluded; recent transmission judged likely at or below this distance.",
        note="Banded, not binary: above 12 allele differences recent transmission is judged unlikely, and 6 to 12 is explicitly INDETERMINATE and must never render as unrelated. Derived by transferring the SNP thresholds of Walker et al. 2013 onto the same 390-genome UK cohort, so it is a method transfer, not an independent derivation. A reduced, extended or re-derived target set voids the cutoff. No paediatric derivation or validation exists."))
    entries.append(_entry("Salmonella enterica", 10, "leeper2023", loci=3002,
        scheme_key="enterobase:senterica-cgmlst-3002", scope="PulseNet USA national surveillance cluster-detection criterion on the EnteroBase-derived 3,002-locus core scheme, evaluated but not derived by this study.",
        note="An operational surveillance criterion, not a validated outbreak or transmission cutoff. It is only half the rule: PulseNet also requires at least seven clinical cases, or three for rarer serotypes, within 60 days. The evaluation clustered by UPGMA, so do not assume single linkage. EnteroBase HC5 is a hierarchical clustering level, not this criterion and not a validated threshold. ECDC and EFSA set their Salmonella cutoff per outbreak, so no single European number is curated here."))
    entries.append(_entry("Legionella pneumophila", 4, "morangilad2015", loci=1521,
        scheme_key="cgmlst.org:lpneumophila-1521", scope="ESGLI-associated 1,521-target scheme seeded on Philadelphia-1 (NC_002942.5); preliminary cluster-type distance.",
        missing_policy="Calibration compared the 1,446 of 1,521 targets shared by all analysed genomes; missing targets excluded rather than counted as differences.",
        note="PRELIMINARY by the authors' own wording, from only 17 genomes, and never independently revalidated; they state it should be further evaluated and fine-tuned. Linked isolates actually differed by 0, 1 and 3 alleles, so 4 is a rounded bound, not a measured maximum. David et al. 2016 (10.1128/JCM.00432-16) question this scheme's resolution for outbreak use. It is a linkage distance, so report cluster diameter. A separate local 4-allele criterion exists under siddall2025 with no bound scheme; the numeric coincidence is not independent confirmation."))
    return deepcopy(entries)


def guidance_for(organism, method=None):
    """Exact taxon matching only: no genus fallback, taxonomic or method inference."""
    text = " ".join(str(organism).split()).casefold()
    entries = [entry for entry in catalog_entries() if entry["organism"].casefold() == text
               and (method is None or entry["method"] == method)]
    for entry in entries:
        entry["source"] = deepcopy(SOURCES[entry["source_id"]])
    return {"organism": str(organism), "catalog_version": CATALOG_VERSION, "reviewed_on": REVIEWED_ON,
            "status": "published_contexts_require_review" if entries else "no_curated_transferable_cutoff",
            "message": "No reviewed transferable cutoff in this catalog. This is an evidence gap, not proof that no publications exist." if not entries else INTERPRETATION,
            "entries": entries, "interpretation": INTERPRETATION,
            "review_scope": "Targeted primary-source review, including 2025 literature; not a systematic or continuously updated review."}


def review_age_days(today=None):
    return ((today or date.today()) - date.fromisoformat(REVIEWED_ON)).days


def record_decision(entry_id, context, *, selected_threshold=None, justification="",
                    protocol_reviewed=False, schema_reviewed=False, epi_reviewed=False):
    """Freeze citation + actual local protocol. Never call an adaptation validated.

    An entry without an exact curated target set can be cited for context, but
    cannot supply an approved threshold. Application requires matching method,
    explicit scheme binding, matching target count, and a local justification.
    """
    entry = next((entry for entry in catalog_entries() if entry["id"] == entry_id), None)
    if entry is None:
        raise ValueError("Unknown publication-guidance entry.")
    context = deepcopy(context)
    source = deepcopy(SOURCES[entry["source_id"]])
    approved = None
    if selected_threshold is not None:
        if isinstance(selected_threshold, bool) or not isinstance(selected_threshold, int) or selected_threshold < 0:
            raise ValueError("Threshold must be a nonnegative integer.")
        if entry["published_threshold"] is None:
            raise ValueError("This entry does not supply a reviewed numeric cutoff.")
        if context.get("method") != entry["method"] or context.get("organism") != entry["organism"]:
            raise ValueError("Organism and comparison method must exactly match the cited context.")
        if not entry["scheme_key"] or context.get("scheme_key") != entry["scheme_key"]:
            raise ValueError("Bind the exact published scheme/protocol first; an organism or locus count alone is insufficient.")
        if entry["locus_count"] is not None and context.get("locus_count") != entry["locus_count"]:
            raise ValueError("The full published target set is required; no automatic scaling for missing or accessory loci.")
        if not context.get("scheme_digest") or not context.get("missing_policy") or not context.get("caller"):
            raise ValueError("Record reference fingerprint, caller and missing-data policy.")
        if not all((protocol_reviewed, schema_reviewed, epi_reviewed)) or len(justification.strip()) < 20:
            raise ValueError("Review scheme, protocol and epidemiology, and record a meaningful local justification.")
        approved = selected_threshold
    payload = {"format_version": 1, "catalog_version": CATALOG_VERSION, "reviewed_on": REVIEWED_ON,
               "created_at": datetime.now(timezone.utc).isoformat(), "entry_id": entry_id,
               "citation": source["citation"], "doi": source["doi"], "url": source["url"], "source": source,
               "published_threshold": entry["published_threshold"], "approved_threshold": approved,
               "scheme_scope": entry["scheme_key"] or "Study-specific scheme not curated",
               "protocol_scope": entry["scope"], "context": context,
               "justification": justification.strip(), "status": "local_adaptation_requires_validation" if approved is not None else "citation_only_no_threshold_change",
               "departs_from_published_number": approved is not None and approved != entry["published_threshold"],
               "review_attestations": {"scheme": schema_reviewed, "protocol": protocol_reviewed, "epidemiology": epi_reviewed},
               "limitations": entry["limitations"], "interpretation": INTERPRETATION}
    payload["evidence_digest"] = hashlib.sha256(json.dumps(payload, sort_keys=True, ensure_ascii=False).encode()).hexdigest()
    return payload
