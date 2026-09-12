"""Native cgMLST: exact reference calls and guarded full-CDS novel inference.

Novel alleles are local sequence identities, never registered allele numbers or
nearest-reference STs. BLAST searches reference alleles against the small predicted
CDS database, so a highly variable locus cannot exhaust another locus's hit cap.
"""

from __future__ import annotations

import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
import time
from collections import defaultdict
from pathlib import Path

import pyrodigal

from . import __version__
from .sequence import SequenceReader, check_cancelled, file_sha256, sample_name
from .typing import Scheme, call_assembly, load_scheme


class CGTypingError(ValueError):
    """The cgMLST analysis could not produce a defensible result."""


def resolve_blast(tool, explicit=None):
    if tool not in {"blastn", "makeblastdb"}:
        raise ValueError("Unsupported nucleotide typing executable.")
    if explicit is not None:
        path = Path(explicit)
        if not path.is_file():
            raise CGTypingError(f"The selected {tool} executable does not exist: {path}")
        return str(path.resolve())
    suffix = ".exe" if os.name == "nt" else ""
    roots = [
        Path(getattr(sys, "_MEIPASS", Path(sys.executable).parent)) / "Tools" / "blast" / "bin",
        Path(sys.executable).parent / "Tools" / "blast" / "bin",
        Path(__file__).resolve().parent / "resources" / "tools" / "blast" / "bin",
    ]
    for root in roots:
        candidate = root / (tool + suffix)
        if candidate.is_file():
            return str(candidate)
    found = shutil.which(tool)
    if found:
        return found
    raise CGTypingError(f"Native NCBI BLAST+ ({tool}) is required for full-CDS cgMLST inference.")


def _run(command, *, cancelled=None, output=None):
    check_cancelled(cancelled)
    with tempfile.TemporaryFile() as errors, tempfile.TemporaryFile() as captured:
        flags = subprocess.CREATE_NO_WINDOW if os.name == "nt" else 0
        process = subprocess.Popen(command, stdout=captured if output is None else output,
                                   stderr=errors, creationflags=flags)
        try:
            while process.poll() is None:
                check_cancelled(cancelled)
                time.sleep(0.05)
        except BaseException:
            process.terminate()
            try:
                process.wait(timeout=2)
            except subprocess.TimeoutExpired:
                process.kill()
                process.wait()
            raise
        errors.seek(0)
        message = errors.read(12000).decode("utf-8", errors="replace")
        if process.returncode:
            raise CGTypingError(f"{Path(command[0]).name} failed ({process.returncode}): {message.strip()}")
        captured.seek(0)
        return captured.read().decode("utf-8", errors="replace")


def full_cds_qc(sequence, *, partial_begin=False, partial_end=False, genetic_code=11):
    sequence = sequence.upper()
    reasons = []
    if partial_begin or partial_end:
        reasons.append("contig-edge partial CDS")
    if not sequence or set(sequence) - set("ACGT"):
        reasons.append("ambiguous or empty nucleotide sequence")
    if len(sequence) % 3:
        reasons.append("length is not a multiple of three")
    if sequence[:3] not in {"ATG", "GTG", "TTG"}:
        reasons.append("no supported complete start codon")
    if genetic_code not in {4, 11}:
        reasons.append("unsupported genetic code for full-CDS validation")
    stops = {"TAA", "TAG"} if genetic_code == 4 else {"TAA", "TAG", "TGA"}
    if sequence[-3:] not in stops:
        reasons.append("no complete terminal stop codon")
    if any(sequence[index:index + 3] in stops for index in range(3, max(3, len(sequence) - 3), 3)):
        reasons.append("internal in-frame stop codon")
    return {"valid": not reasons, "reasons": reasons, "genetic_code": genetic_code,
            "length": len(sequence), "partial_begin": bool(partial_begin),
            "partial_end": bool(partial_end)}


def predict_cds(path, *, genetic_code=11, cancelled=None, progress=None):
    if isinstance(genetic_code, bool) or genetic_code not in {4, 11}:
        raise ValueError('Full-CDS typing currently supports genetic codes 11 and 4 only.')
    records = []
    with SequenceReader(path, cancelled) as reader:
        if reader.kind != "fasta":
            raise CGTypingError("cgMLST requires a FASTA assembly.")
        records = list(reader)
    total_bases = sum(len(record.sequence) for record in records)
    meta = total_bases < 20000
    bins = pyrodigal.MetagenomicBins([item for item in pyrodigal.METAGENOMIC_BINS
                                    if item.training_info.translation_table == genetic_code])
    finder = pyrodigal.GeneFinder(meta=meta, closed=False, metagenomic_bins=bins)
    if not meta:
        check_cancelled(cancelled)
        finder.train(*(record.sequence.encode("ascii") for record in records),
                     translation_table=genetic_code)
    genes = []
    for index, record in enumerate(records, 1):
        check_cancelled(cancelled)
        predicted = finder.find_genes(record.sequence.encode("ascii"))
        for gene in predicted:
            sequence = gene.sequence().upper()
            code = gene.translation_table if meta else genetic_code
            genes.append({
                "id": f"cds{len(genes)}", "contig": record.identifier,
                "start": gene.begin, "end": gene.end, "strand": "+" if gene.strand == 1 else "-",
                "sequence": sequence, "sequence_sha256": hashlib.sha256(sequence.encode()).hexdigest(),
                "cds_qc": full_cds_qc(sequence, partial_begin=gene.partial_begin,
                                      partial_end=gene.partial_end, genetic_code=code),
            })
        if progress:
            progress(index, len(records), f"Predicted coding sequences in {record.identifier}")
    return genes, {"tool": "pyrodigal", "version": pyrodigal.__version__,
                   "mode": "meta-small-input" if meta else "single", "genetic_code": genetic_code}


def _apply_inference(result, scheme, genes, alignments, min_identity, min_coverage):
    """Interpret all qualifying CDS/locus relationships, preserving ambiguity."""
    by_gene = {gene["id"]: gene for gene in genes}
    accepted = defaultdict(dict)
    rejected = {}
    for hit in alignments:
        gene = by_gene[hit["gene_id"]]
        coverage = min(hit["query_coverage"], hit["subject_coverage"])
        if hit["identity"] < min_identity or coverage < min_coverage:
            if hit["bitscore"] > rejected.get(hit["locus"], {}).get("bitscore", -1):
                rejected[hit["locus"]] = hit
            continue
        if not gene["cds_qc"]["valid"]:
            if hit["bitscore"] > rejected.get(hit["locus"], {}).get("bitscore", -1):
                rejected[hit["locus"]] = hit
            continue
        previous = accepted[hit["locus"]].get(gene["id"])
        if previous is None or hit["bitscore"] > previous["bitscore"]:
            accepted[hit["locus"]][gene["id"]] = hit
    gene_loci = defaultdict(set)
    for locus, genes_for_locus in accepted.items():
        for gene_id in genes_for_locus:
            gene_loci[gene_id].add(locus)
    novel_sequences = []
    for call in result["calls"]:
        locus = call["locus"]
        candidates = [by_gene[identifier] for identifier in accepted[locus]]
        if call["status"] in {"mixed", "ambiguous"}:
            continue
        if call["status"] == "exact":
            known = scheme.alleles[locus][call["allele"]]
            # A second full homologous CDS is a copy-number ambiguity, even if
            # its sequence is novel and therefore invisible to exact matching.
            extra = [gene for gene in candidates if gene["sequence"] != known and not any(
                gene["contig"] == hit["contig"] and gene["start"] <= hit["end"]
                and hit["start"] <= gene["end"] for hit in call["hits"])]
            if extra:
                call.update(allele=None, status="mixed",
                            reason="A distinct full-length homologous CDS accompanies the exact allele.")
                result["alleles"][locus] = None
            continue
        if not candidates:
            if locus in rejected:
                best = rejected[locus]
                gene = by_gene[best["gene_id"]]
                call.update(status="partial" if not gene["cds_qc"]["valid"] else "low_similarity",
                            cds_qc=gene["cds_qc"], best_alignment=best,
                            reason="A homolog was found, but complete-CDS or alignment criteria failed.")
            continue
        if len(candidates) != 1 or len(gene_loci[candidates[0]["id"]]) != 1:
            call.update(allele=None, status="ambiguous",
                        reason="Multiple CDS copies or cross-locus homology prevent unique assignment.",
                        cds_candidates=[gene["id"] for gene in candidates])
            continue
        gene = candidates[0]
        hit = accepted[locus][gene["id"]]
        # Full-CDS identity is independent of the nearest known allele number.
        digest = gene["sequence_sha256"]
        token = "NOVEL_" + digest
        call.update(allele=token, status="novel_validated", sequence_sha256=digest,
                    hit_count=1, candidates=[token], evidence_truncated=False,
                    cds_qc=gene["cds_qc"], identity_percent=hit["identity"] * 100,
                    query_coverage=hit["query_coverage"], subject_coverage=hit["subject_coverage"],
                    nearest_reference_allele=hit["reference_allele"],
                    hits=[{key: gene[key] for key in ("contig", "start", "end", "strand")}],
                    reason="Unambiguous complete CDS homolog; local sequence identity, not a registered allele.")
        result["alleles"][locus] = token
        novel_sequences.append({"locus": locus, "allele": token,
                                "sequence_sha256": digest, "sequence": gene["sequence"]})
    statuses = {call["status"] for call in result["calls"]}
    result["novel_sequences"] = novel_sequences
    if "mixed" in statuses:
        result["status"], result["st"] = "mixed", None
    elif "ambiguous" in statuses:
        result["status"], result["st"] = "ambiguous", None
    elif any(value is None for value in result["alleles"].values()):
        result["status"], result["st"] = "incomplete", None
    elif novel_sequences:
        result["status"], result["st"] = "novel_alleles", None
    vector = [(locus, result["alleles"][locus]) for locus in sorted(result["alleles"])]
    result["cg_profile_digest"] = hashlib.sha256(
        json.dumps({"scheme": scheme.digest, "alleles": vector}, sort_keys=True).encode()).hexdigest()
    result["cg_profile_id"] = "LOCAL_" + result["cg_profile_digest"][:16]
    result["cg_profile_complete"] = all(value is not None for value in result["alleles"].values())
    result.setdefault("notes", []).append(
        "Novel CDS alleles use full SHA-256 sequence identities. LOCAL profile identifiers are not registered cgSTs.")
    return result


def call_cgassembly(path, scheme, cancelled=None, progress=None, *, cache_dir=None,
                    blastn_path=None, makeblastdb_path=None, threads=2,
                    min_identity=0.90, min_coverage=0.98, genetic_code=11):
    """Exact allele calls followed by full-CDS inference, with local provenance."""
    if not 0 < min_identity <= 1 or not 0 < min_coverage <= 1:
        raise ValueError("cgMLST identity and coverage must be fractions in (0, 1].")
    if not isinstance(threads, int) or threads < 1:
        raise ValueError("The BLAST thread count must be a positive integer.")
    path = Path(path).resolve()
    if not isinstance(scheme, Scheme):
        scheme = load_scheme(scheme, cancelled=cancelled)
    blastn, makeblastdb = resolve_blast("blastn", blastn_path), resolve_blast("makeblastdb", makeblastdb_path)
    blast_version = _run([blastn, "-version"], cancelled=cancelled).splitlines()[0]
    database_version = _run([makeblastdb, "-version"], cancelled=cancelled).splitlines()[0]
    input_digest = file_sha256(path, cancelled)
    parameters = {"method": "full-cds-cgmlst-v2", "min_identity": min_identity,
                  "min_coverage": min_coverage, "genetic_code": genetic_code,
                  "pyrodigal": pyrodigal.__version__, "blastn": blast_version,
                  "makeblastdb": database_version, "threads": threads,
                  "blastn_sha256": file_sha256(blastn, cancelled),
                  "makeblastdb_sha256": file_sha256(makeblastdb, cancelled),
                  "input_sha256": input_digest, "scheme_digest": scheme.digest,
                  "software_version": __version__}
    key = hashlib.sha256(json.dumps(parameters, sort_keys=True).encode()).hexdigest()
    cache_root = Path(cache_dir or Path(tempfile.gettempdir()) / "wmlstudio-cgcache").resolve()
    cache_root.mkdir(parents=True, exist_ok=True)
    cache_path = cache_root / f"{key}.json"
    if cache_path.is_file():
        try:
            cached = json.loads(cache_path.read_text(encoding="utf-8"))
            encoded = json.dumps(cached["result"], sort_keys=True).encode()
            if cached["parameters"] == parameters and hashlib.sha256(encoded).hexdigest() == cached["sha256"]:
                cached["result"]["input_path"] = str(path)
                cached["result"]["sample_name"] = sample_name(path)
                check_cancelled(cancelled)
                return cached["result"]
        except (OSError, ValueError, KeyError, TypeError):
            pass
    result = call_assembly(path, scheme, cancelled=cancelled, progress=progress)
    result['notes'] = [note for note in result['notes']
                       if note != 'Exact nucleotide matching only; missing matches do not establish novel alleles.']
    result['notes'].append('Known alleles require exact assembly matches; novel inference additionally requires an unambiguous complete predicted CDS.')
    genes, prediction = predict_cds(path, genetic_code=genetic_code, cancelled=cancelled, progress=progress)
    if not genes:
        result["notes"].append("No CDS predictions were available; only exact assembly matches are reported.")
        result = _apply_inference(result, scheme, [], [], min_identity, min_coverage)
    else:
        with tempfile.TemporaryDirectory(prefix="wmlstudio-cg-", dir=cache_root) as work:
            work = Path(work)
            cds_path, query_path = work / "cds.fasta", work / "references.fasta"
            with cds_path.open("w", encoding="ascii") as handle:
                for gene in genes:
                    handle.write(f">{gene['id']}\n{gene['sequence']}\n")
            queries = {}
            with query_path.open("w", encoding="ascii") as handle:
                for locus in scheme.loci:
                    check_cancelled(cancelled)
                    for allele, sequence in scheme.alleles[locus].items():
                        if set(sequence) - set("ACGT"):
                            continue
                        identifier = f"ref{len(queries)}"
                        queries[identifier] = (locus, allele)
                        handle.write(f">{identifier}\n{sequence}\n")
            database = work / "cdsdb"
            if progress:
                progress(0, 1, "Indexing predicted CDS for cgMLST")
            _run([makeblastdb, "-in", str(cds_path), "-dbtype", "nucl",
                  "-out", str(database), "-blastdb_version", "4"], cancelled=cancelled)
            output = work / "alignments.tsv"
            if progress:
                progress(0, 1, "Comparing complete reference alleles with predicted CDS")
            _run([blastn, "-task", "megablast", "-query", str(query_path), "-db", str(database),
                  "-out", str(output), "-outfmt",
                  "6 qseqid sseqid pident length qstart qend sstart send bitscore qlen slen",
                  "-num_threads", str(threads), "-perc_identity", str(min_identity * 100),
                  "-qcov_hsp_perc", str(min_coverage * 100), "-max_target_seqs", str(len(genes)),
                  "-max_hsps", "1", "-dust", "no", "-evalue", "1e-10"], cancelled=cancelled)

            def alignments():
                with output.open(encoding="ascii") as handle:
                    for index, line in enumerate(handle):
                        if index % 4096 == 0:
                            check_cancelled(cancelled)
                        fields = line.rstrip("\n").split("\t")
                        if len(fields) != 11 or fields[0] not in queries:
                            raise CGTypingError("BLAST returned malformed or unexpected alignment records.")
                        locus, allele = queries[fields[0]]
                        yield {"locus": locus, "reference_allele": allele, "gene_id": fields[1],
                               "identity": float(fields[2]) / 100,
                               "query_coverage": (abs(int(fields[5]) - int(fields[4])) + 1) / int(fields[9]),
                               "subject_coverage": (abs(int(fields[7]) - int(fields[6])) + 1) / int(fields[10]),
                               "bitscore": float(fields[8])}
            result = _apply_inference(result, scheme, genes, alignments(), min_identity, min_coverage)
    if file_sha256(path, cancelled) != input_digest:
        raise CGTypingError("Input changed during cgMLST analysis; retry.")
    result["parameters"] = {**result["parameters"], **parameters, "threads": threads}
    result["gene_prediction"] = prediction
    result["tools"] = {"blastn": {"path": blastn, "version": blast_version},
                       "makeblastdb": {"path": makeblastdb, "version": database_version},
                       "pyrodigal": pyrodigal.__version__}
    payload = {"parameters": parameters, "result": result,
               "sha256": hashlib.sha256(json.dumps(result, sort_keys=True).encode()).hexdigest()}
    temporary = None
    try:
        with tempfile.NamedTemporaryFile("w", encoding="utf-8", dir=cache_root,
                                         suffix=".part", delete=False) as handle:
            temporary = Path(handle.name)
            json.dump(payload, handle)
            handle.flush()
            os.fsync(handle.fileno())
        check_cancelled(cancelled)
        os.replace(temporary, cache_path)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)
    return result
