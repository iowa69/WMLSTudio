"""Conservative automatic selection of installed classical MLST panels.

A shared seed index screens every allele, then the normal complete-allele engine
confirms candidate schemes. A seed match alone never assigns an organism or ST.
The result is compatibility with an installed panel, not a species diagnostic.
"""

from __future__ import annotations

import math
import re
import threading
import urllib.parse
from collections import OrderedDict
from pathlib import Path

import ahocorasick

from .reference_catalog import organism_parts
from .sequence import SequenceError, SequenceReader, check_cancelled, file_signature
from .typing import ALLELE_SUFFIXES, call_assembly, load_scheme, reverse_complement

# Public database labels verified against https://rest.pubmlst.org/db on
# 2026-09-12. These describe schema scope, not a validated species classifier.
_SOURCE_LABELS = {
    "aactinomycetemcomitans": "Aggregatibacter actinomycetemcomitans",
    "abaumannii": "Acinetobacter baumannii", "achromobacter": "Achromobacter spp.",
    "actinobacillus": "Actinobacillus spp.", "aeromonas": "Aeromonas spp.",
    "afumigatus": "Aspergillus fumigatus", "aparagallinarum": "Avibacterium paragallinarum",
    "aphagocytophilum": "Anaplasma phagocytophilum", "arcobacter": "Arcobacter spp.",
    "bbacilliformis": "Bartonella bacilliformis", "bcc": "Burkholderia cepacia complex",
    "bcereus": "Bacillus cereus", "bfragilis": "Bacteroides fragilis",
    "bhenselae": "Bartonella henselae", "blastocystis": "Blastocystis spp.",
    "blicheniformis": "Bacillus licheniformis", "bmallei": "Burkholderia mallei",
    "borrelia": "Borrelia spp.", "bpseudomallei": "Burkholderia pseudomallei",
    "brachyspira": "Brachyspira spp.", "brucella": "Brucella spp.",
    "bsubtilis": "Bacillus subtilis", "bwashoensis": "Bartonella washoensis",
    "calbicans": "Candida albicans", "campylobacter_nonjejuni": "Campylobacter non-jejuni/coli",
    "campylobacter": "Campylobacter jejuni/coli", "cauris": "Candida auris",
    "cavidum": "Cutibacterium avidum", "cbotulinum": "Clostridium botulinum",
    "cchauvoei": "Clostridium chauvoei", "cdifficile": "Clostridioides difficile",
    "cfreundii": "Citrobacter spp.", "cglabrata": "Candida glabrata",
    "chlamydiales": "Chlamydiales spp.", "ckrusei": "Candida krusei",
    "cmaltaromaticum": "Carnobacterium maltaromaticum", "cperfringens": "Clostridium perfringens",
    "cronobacter": "Cronobacter spp.", "csepticum": "Clostridium septicum",
    "csinensis": "Clonorchis sinensis", "ctropicalis": "Candida tropicalis",
    "dnodosus": "Dichelobacter nodosus", "ecloacae": "Enterobacter spp.",
    "edwardsiella": "Edwardsiella spp.", "efaecalis": "Enterococcus faecalis",
    "efaecium": "Enterococcus faecium", "escherichia": "Escherichia spp.",
    "fpsychrophilum": "Flavobacterium psychrophilum", "gallibacterium": "Gallibacterium anatis",
    "geotrichum": "Geotrichum spp.", "hcinaedi": "Helicobacter cinaedi",
    "helicobacter": "Helicobacter pylori", "hinfluenzae": "Haemophilus influenzae",
    "hparasuis": "Glaesserella parasuis", "hsuis": "Helicobacter suis",
    "kaerogenes": "Klebsiella aerogenes", "koxytoca": "Klebsiella oxytoca",
    "kseptempunctata": "Kudoa septempunctata", "leptospira": "Leptospira spp.",
    "lgarvieae": "Lactococcus garvieae", "liberibacter": "Candidatus Liberibacter solanacearum",
    "lsalivarius": "Lactobacillus salivarius", "mabscessus": "Mycobacteroides abscessus complex",
    "magalactiae": "Mycoplasma agalactiae", "manserisalpingitidis": "Mycoplasma anserisalpingitidis",
    "mbovis": "Mycoplasma bovis", "mcanis": "Macrococcus canis",
    "mcaseolyticus": "Macrococcus caseolyticus", "mflocculare": "Mycoplasma flocculare",
    "mgallisepticum": "Mycoplasma gallisepticum", "mgenitalium": "Mycoplasma genitalium",
    "mhaemolytica": "Mannheimia haemolytica", "mhominis": "Mycoplasma hominis",
    "mhyopneumoniae": "Mycoplasma hyopneumoniae", "mhyorhinis": "Mycoplasma hyorhinis",
    "mhyosynoviae": "Mycoplasma hyosynoviae", "miowae": "Mycoplasma iowae",
    "mplutonius": "Melissococcus plutonius", "mpneumoniae": "Mycoplasma pneumoniae",
    "msciuri": "Mammaliicoccus sciuri", "msynoviae": "Mycoplasma synoviae",
    "mycobacteria": "Mycobacteria spp.", "neisseria": "Neisseria spp.",
    "oralstrep": "Oral Streptococcus spp.", "orhinotracheale": "Ornithobacterium rhinotracheale",
    "otsutsugamushi": "Orientia tsutsugamushi", "pacnes": "Cutibacterium acnes",
    "paeruginosa": "Pseudomonas aeruginosa", "pdamselae": "Photobacterium damselae",
    "pfluorescens": "Pseudomonas fluorescens", "pgingivalis": "Porphyromonas gingivalis",
    "plarvae": "Paenibacillus larvae", "pmultocida": "Pasteurella multocida",
    "ppentosaceus": "Pediococcus pentosaceus", "pputida": "Pseudomonas putida",
    "proteus": "Proteus spp.", "providencia": "Providencia spp.",
    "psalmonis": "Piscirickettsia salmonis", "ranatipestifer": "Riemerella anatipestifer",
    "rhodococcus": "Rhodococcus spp.", "sagalactiae": "Streptococcus agalactiae",
    "salmonella": "Salmonella spp.", "saureus": "Staphylococcus aureus",
    "sbsec": "Streptococcus bovis/equinus complex (SBSEC)", "scanis": "Streptococcus canis",
    "scapitis": "Staphylococcus capitis", "schromogenes": "Staphylococcus chromogenes",
    "sdysgalactiae": "Streptococcus dysgalactiae", "sepidermidis": "Staphylococcus epidermidis",
    "serratia": "Serratia spp.", "sgallolyticus": "Streptococcus gallolyticus",
    "shaemolyticus": "Staphylococcus haemolyticus", "shewanella": "Shewanella spp.",
    "shominis": "Staphylococcus hominis", "siniae": "Streptococcus iniae",
    "sinorhizobium": "Sinorhizobium spp.", "smaltophilia": "Stenotrophomonas maltophilia",
    "smitis": "Streptococcus mitis", "sparasitica": "Saprolegnia parasitica",
    "spneumoniae": "Streptococcus pneumoniae", "spseudintermedius": "Staphylococcus pseudintermedius",
    "spyogenes": "Streptococcus pyogenes", "ssuis": "Streptococcus suis",
    "sthermophilus": "Streptococcus thermophilus", "streptomyces": "Streptomyces spp",
    "suberis": "Streptococcus uberis", "szooepidemicus": "Streptococcus zooepidemicus",
    "taylorella": "Taylorella spp.", "tenacibaculum": "Tenacibaculum spp.",
    "tpallidum": "Treponema pallidum", "tpyogenes": "Trueperella pyogenes",
    "ureaplasma": "Ureaplasma spp.", "vcholerae": "Vibrio cholerae", "vibrio": "Vibrio spp.",
    "vparahaemolyticus": "Vibrio parahaemolyticus", "vtapetis": "Vibrio tapetis",
    "vvulnificus": "Vibrio vulnificus", "wolbachia": "Wolbachia spp.",
    "xcitri": "Xanthomonas citri", "xfastidiosa": "Xylella fastidiosa",
    "ypseudotuberculosis_achtman": "Yersinia pseudotuberculosis", "yruckeri": "Yersinia ruckeri",
}
_INDEX_CACHE: OrderedDict = OrderedDict()
_SCHEME_CACHE: OrderedDict = OrderedDict()
_CACHE_LOCK = threading.Lock()


def _signature(path: Path) -> tuple:
    if not path.is_dir():
        return (str(path), "missing")
    return (str(path), tuple((p.name, file_signature(p)) for p in sorted(path.iterdir())
                            if p.is_file() and p.suffix.casefold() in
                            ALLELE_SUFFIXES | {".json", ".txt", ".tsv", ".gz", ".bz2"}))


def cached_scheme(path, cancelled=None):
    path = Path(path).resolve()
    signature = _signature(path)
    check_cancelled(cancelled)
    with _CACHE_LOCK:
        if signature in _SCHEME_CACHE:
            _SCHEME_CACHE.move_to_end(signature)
            return _SCHEME_CACHE[signature]
    scheme = load_scheme(path, cancelled=cancelled)
    if _signature(path) != signature:
        raise SequenceError("Reference changed while preparing the typing cache.")
    with _CACHE_LOCK:
        if len(scheme.loci) > 30:
            # Keep only one large reference: cgMLST allele snapshots may occupy
            # hundreds of MB before their exact-matching index is constructed.
            for old_key, old_scheme in list(_SCHEME_CACHE.items()):
                if len(old_scheme.loci) > 30:
                    del _SCHEME_CACHE[old_key]
        _SCHEME_CACHE[signature] = scheme
        while len(_SCHEME_CACHE) > 4:
            _SCHEME_CACHE.popitem(last=False)
    return scheme


def scheme_organism(scheme) -> tuple[str, dict]:
    metadata = scheme.metadata
    label = metadata.get("organism")
    if isinstance(label, dict):
        parts = {key: str(label.get(key) or "") for key in ("genus", "species")}
        return " ".join(value for value in parts.values() if value), parts
    if isinstance(label, str) and label:
        return label, organism_parts(label)
    if metadata.get("genus"):
        label = " ".join(str(metadata.get(key) or "") for key in ("genus", "species")).strip()
        return label, organism_parts(label)
    api = str(metadata.get("API") or "")
    match = re.search(r"/db/pubmlst_(.+?)_seqdef/", urllib.parse.urlsplit(api).path)
    label = _SOURCE_LABELS.get(match.group(1), "") if match else ""
    return label or scheme.name, organism_parts(label)


def clear_identification_cache():
    with _CACHE_LOCK:
        _INDEX_CACHE.clear()
        _SCHEME_CACHE.clear()


def _build_index(paths, cancelled, progress, min_loci, max_loci):
    key = (min_loci, max_loci, tuple(_signature(path) for path in paths))
    with _CACHE_LOCK:
        if key in _INDEX_CACHE:
            return _INDEX_CACHE[key]
    automaton = ahocorasick.Automaton()
    panels, excluded = [], []
    digests = set()
    for index, path in enumerate(paths, 1):
        check_cancelled(cancelled)
        try:
            # Count before parsing large cgMLST references: discovery uses small
            # typing panels; large schemes remain available for selected typing.
            count = sum((p.with_suffix("") if p.suffix.casefold() in {".gz", ".bz2"} else p)
                        .suffix.casefold() in ALLELE_SUFFIXES for p in path.iterdir() if p.is_file())
            if not min_loci <= count <= max_loci:
                excluded.append({"scheme_path": str(path), "reason": f"Discovery uses {min_loci}–{max_loci} locus panels; this scheme has {count}."})
                continue
            scheme = cached_scheme(path, cancelled)
            if scheme.digest in digests:
                continue
            label, organism = scheme_organism(scheme)
            biological_scope = ' '.join(str(scheme.metadata.get(key, ''))
                                        for key in ('API', 'name', 'organism', 'type', 'scope'))
            if "plasmid" in biological_scope.casefold():
                excluded.append({"scheme_path": str(path), "reason": "Plasmid panels do not establish host organism identity."})
                continue
            panel = {"scheme": scheme.name, "scheme_path": str(path), "scheme_digest": scheme.digest,
                     "locus_count": len(scheme.loci), "organism_label": label, "organism": organism}
            number = len(panels)
            for locus, alleles in scheme.alleles.items():
                seeds = set()
                for sequence in alleles.values():
                    check_cancelled(cancelled)
                    if set(sequence) - set("ACGT"):
                        continue
                    size = min(31, len(sequence))
                    start = (len(sequence) - size) // 2
                    seed = sequence[start:start + size]
                    seeds.add(seed)
                    seeds.add(reverse_complement(seed))
                for seed in seeds:
                    members = automaton.get(seed, set())
                    members.add((number, locus))
                    automaton.add_word(seed, members)
            panels.append(panel)
            digests.add(scheme.digest)
        except (OSError, ValueError) as error:
            excluded.append({"scheme_path": str(path), "reason": str(error)})
        if progress:
            progress(index, len(paths), f"Preparing automatic typing: {path.name}")
    if len(automaton):
        automaton.make_automaton()
    if key != (min_loci, max_loci, tuple(_signature(path) for path in paths)):
        raise SequenceError("Installed schemes changed while building the discovery index; retry.")
    result = (automaton, panels, excluded)
    with _CACHE_LOCK:
        _INDEX_CACHE.clear()
        _INDEX_CACHE[key] = result
    return result


def identify_assembly(path, scheme_paths, cancelled=None, progress=None, *,
                      min_loci=4, min_coverage=0.8, min_margin=0.15, max_panel_loci=30) -> dict:
    """Return provisional panel compatibility, selecting only a strong unique scheme.

    Coverage and margin are fractions, not posterior probabilities. No taxonomic
    species claim follows from a match; group/complex references stay genus-only.
    """
    if (not isinstance(min_loci, int) or isinstance(min_loci, bool) or min_loci < 1
            or not 0 < min_coverage <= 1 or not 0 <= min_margin <= 1
            or not math.isfinite(min_coverage) or not math.isfinite(min_margin)):
        raise ValueError("Discovery requires positive locus counts and finite coverage/margin fractions.")
    path = Path(path).resolve()
    input_signature = file_signature(path)
    paths = sorted({Path(value).resolve() for value in scheme_paths}, key=str)
    check_cancelled(cancelled)
    automaton, panels, excluded = _build_index(paths, cancelled, progress, min_loci, max_panel_loci)
    supported = [set() for _ in panels]
    with SequenceReader(path, cancelled) as reader:
        if reader.kind != "fasta":
            raise SequenceError("Automatic MLST discovery needs a FASTA assembly; reads must be assembled first.")
        for record in reader:
            check_cancelled(cancelled)
            if not len(automaton):
                continue
            iterator = automaton.iter("")
            for start in range(0, len(record.sequence), 65536):
                check_cancelled(cancelled)
                iterator.set(record.sequence[start:start + 65536], False)
                for index, (_end, members) in enumerate(iterator):
                    if index % 4096 == 0:
                        check_cancelled(cancelled)
                    for number, locus in members:
                        supported[number].add(locus)
    candidates, confirmed = [], {}
    shortlist = [(index, panel) for index, panel in enumerate(panels)
                 if len(supported[index]) >= min_loci
                 and len(supported[index]) / panel["locus_count"] >= max(0, min_coverage - min_margin)]
    for done, (index, panel) in enumerate(shortlist, 1):
        check_cancelled(cancelled)
        if progress:
            progress(done, len(shortlist), f"Confirming full alleles: {panel['scheme']}")
        result = call_assembly(path, cached_scheme(panel["scheme_path"], cancelled), cancelled=cancelled)
        if result["scheme_digest"] != panel["scheme_digest"]:
            raise SequenceError("A candidate scheme changed during automatic identification; retry.")
        count = sum(call["status"] == "exact" for call in result["calls"])
        candidate = {**panel, "matched_loci": count, "coverage": count / panel["locus_count"],
                     "seed_loci": len(supported[index]), "status": result["status"], "st": result["st"]}
        candidates.append(candidate)
        confirmed[panel["scheme_path"]] = result
    candidates.sort(key=lambda item: (-item["coverage"], -item["matched_loci"], item["scheme_path"]))
    eligible = [item for item in candidates if item["matched_loci"] >= min_loci
                and item["coverage"] >= min_coverage and item["status"] not in {"mixed", "ambiguous"}]
    chosen = eligible[0] if eligible else None
    status = "insufficient" if panels else "no_schemes"
    notes = ["Organism labels are provisional compatibility with installed MLST references, not species confirmation.",
             "Absence of mixed allele evidence does not exclude contamination; assembly-based panel matching is not a taxonomic purity assessment."]
    if chosen:
        competing = [item for item in candidates if item is not chosen
                     and item["coverage"] >= chosen["coverage"] - min_margin]
        if competing:
            chosen = None
            status = "ambiguous"
            notes.append("Several schemes have similar exact-locus support; choose a scheme explicitly.")
        else:
            status = "assigned"
    elif any(item["status"] in {"mixed", "ambiguous"} for item in candidates):
        status = "ambiguous"
        notes.append("Mixed or ambiguous allele evidence prevents automatic organism assignment.")
    else:
        notes.append("No installed panel has enough complete, unambiguous allele evidence for assignment.")
    if file_signature(path) != input_signature:
        raise SequenceError("Input changed during automatic scheme discovery; retry.")
    organism = chosen["organism"] if chosen else {"genus": "", "species": ""}
    response = {
        "identification_status": status, "provisional_label": chosen["organism_label"] if chosen else "Unknown",
        **organism, "organism": organism, "best_scheme_path": chosen["scheme_path"] if chosen else None,
        "candidates": candidates, "excluded_schemes": excluded, "notes": notes,
        "parameters": {"min_loci": min_loci, "min_coverage": min_coverage, "min_margin": min_margin,
                       "max_panel_loci": max_panel_loci, "seed_length": 31, "method": "seed-screen-full-exact-confirmation"},
    }
    if chosen:
        response["typing_result"] = confirmed[chosen["scheme_path"]]
    return response
