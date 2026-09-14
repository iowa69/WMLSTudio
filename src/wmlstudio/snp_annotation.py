"""Reference-bound SNP coding context and explicit interval masking.

Consequences describe a single nucleotide change in an annotated/predicted CDS,
not its phenotype, pathogenicity or transmission relevance. Coding selection,
repeat/IS masking and recombination handling are distinct recorded decisions.
"""

from __future__ import annotations

import bisect
import copy
import hashlib
import io
import itertools
import json
from collections import defaultdict
from pathlib import Path
from urllib.parse import unquote

from . import __version__
from .cgtyping import predict_cds
from .sequence import SequenceReader, check_cancelled, file_sha256, file_signature
from .typing import reverse_complement

_CODONS = dict(zip((''.join(bases) for bases in itertools.product('TCAG', repeat=3)),
                   'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG', strict=True))
_STARTS = {11: {'ATG', 'GTG', 'TTG', 'CTG', 'ATT', 'ATC', 'ATA'},
           4: {'ATG', 'GTG', 'TTG', 'CTG', 'ATT', 'ATC', 'ATA', 'TTA'}}


def _sequences(path, cancelled):
    sequences, total = {}, 0
    with SequenceReader(path, cancelled) as reader:
        if reader.kind != 'fasta':
            raise ValueError('SNP annotation requires an exact reference FASTA.')
        for record in reader:
            total += len(record.sequence)
            if total > 100_000_000:
                raise ValueError('Reference exceeds the explicit 100 million base bacterial annotation bound.')
            sequences[record.identifier] = record.sequence
    return sequences


def _attributes(value):
    fields = {}
    for field in value.split(';'):
        if not field:
            continue
        if '=' not in field:
            raise ValueError('GFF3 attributes require key=value syntax; GTF is not supported.')
        key, text = field.split('=', 1)
        if key in fields:
            raise ValueError(f'GFF3 contains a repeated {key!r} attribute.')
        fields[key] = unquote(text)
    return fields


def _embedded_fasta(text):
    sequences, name, parts = {}, None, []
    for line in io.StringIO(text):
        line = line.strip()
        if not line:
            continue
        if line.startswith('>'):
            if name is not None:
                sequences[name] = ''.join(parts).upper()
            if not line[1:].strip():
                raise ValueError('GFF3 embedded FASTA contains an empty identifier.')
            name, parts = line[1:].split()[0], []
            if name in sequences:
                raise ValueError('GFF3 embedded FASTA has duplicate reference identifiers.')
        elif name is None:
            raise ValueError('Malformed GFF3 embedded FASTA.')
        else:
            parts.append(line)
    if name is not None:
        if name in sequences:
            raise ValueError('GFF3 embedded FASTA has duplicate reference identifiers.')
        sequences[name] = ''.join(parts).upper()
    return sequences


def _gff_features(path, sequences, *, reference_sha256, annotation_reference_sha256, cancelled):
    if Path(path).suffix.lower() not in {'.gff', '.gff3'}:
        raise ValueError('Supplied CDS annotation currently supports GFF3 only, not GenBank/GTF. Export verified GFF3 with embedded FASTA or use explicit Prodigal prediction.')
    if Path(path).stat().st_size > 64 * 1024 * 1024:
        raise ValueError('GFF3 exceeds the 64 MiB safety bound.')
    text = Path(path).read_text(encoding='utf-8')
    annotation_text, separator, fasta_text = text.partition('##FASTA')
    if separator:
        embedded = _embedded_fasta(fasta_text)
        if embedded != sequences:
            raise ValueError('GFF3 embedded reference sequences do not exactly match the selected reference FASTA.')
        binding = 'embedded_reference_sequence_verified'
    elif annotation_reference_sha256 is not None:
        if annotation_reference_sha256 != reference_sha256:
            raise ValueError('Declared annotation reference SHA-256 does not match the selected reference FASTA.')
        binding = 'explicit_user_reference_hash_binding'
    else:
        binding = 'coordinate_only_unverified'
    features, groups = [], defaultdict(list)
    for number, line in enumerate(annotation_text.splitlines(), 1):
        check_cancelled(cancelled)
        if not line or line.startswith('#'):
            continue
        fields = line.split('\t')
        if len(fields) != 9:
            raise ValueError(f'GFF3 line {number} does not contain nine tab-separated columns.')
        contig, source, kind, start, end, score, strand, phase, attributes = fields
        if kind != 'CDS':
            continue
        contig = unquote(contig)
        try:
            start, end, phase = int(start), int(end), int(phase)
        except ValueError as error:
            raise ValueError(f'GFF3 CDS line {number} requires integer start/end and phase 0, 1 or 2.') from error
        if contig not in sequences or not 1 <= start <= end <= len(sequences[contig]) or strand not in {'+', '-'} or phase not in {0, 1, 2}:
            raise ValueError(f'GFF3 CDS line {number} is not compatible with the exact reference coordinates.')
        attrs = _attributes(attributes)
        code = int(attrs.get('transl_table') or 11)
        if code not in {4, 11}:
            raise ValueError('SNP consequence annotation supports bacterial genetic codes 11 and 4 only.')
        identifier = attrs.get('ID') or attrs.get('Parent') or f'CDS-line-{number}'
        sequence = sequences[contig][start - 1:end]
        if strand == '-':
            sequence = reverse_complement(sequence)
        partial = attrs.get('partial', '').lower() in {'true', '1'} or any('.' in attrs.get(key, '') for key in ('start_range', 'end_range'))
        feature = {'id': identifier, 'contig': contig, 'start': start, 'end': end, 'strand': strand,
                   'phase': phase, 'sequence': sequence, 'genetic_code': code, 'partial': partial,
                   'gene': attrs.get('gene') or attrs.get('Name') or '', 'locus_tag': attrs.get('locus_tag') or '',
                   'product': attrs.get('product') or '', 'source': source, 'origin': 'supplied_GFF3',
                   'annotation_binding': binding, 'compound': False}
        groups[identifier].append(feature)
        features.append(feature)
        if len(features) > 100_000:
            raise ValueError('GFF3 exceeds the 100,000 CDS feature safety bound.')
    for group in groups.values():
        if len(group) > 1:
            for feature in group:
                feature['compound'] = True
    return features, {'kind': 'supplied_GFF3', 'binding': binding, 'sha256': file_sha256(path, cancelled),
                      'path': str(Path(path).resolve()), 'cds_count': len(features),
                      'multipart_cds_policy': 'Preserved as coding context; consequences unresolved, never silently concatenated.'}


def _bed_intervals(path, sequences, cancelled):
    if Path(path).stat().st_size > 32 * 1024 * 1024:
        raise ValueError('BED mask exceeds the 32 MiB safety bound.')
    intervals = []
    with Path(path).open(encoding='utf-8') as handle:
        for number, line in enumerate(handle, 1):
            check_cancelled(cancelled)
            if not line.strip() or line.startswith(('#', 'track ', 'browser ')):
                continue
            values = line.rstrip().split('\t')
            if len(values) < 3:
                raise ValueError(f'BED line {number} must contain at least three tab-separated fields.')
            contig, start, end = values[:3]
            try:
                start, end = int(start), int(end)
            except ValueError as error:
                raise ValueError(f'BED line {number} has invalid integer coordinates.') from error
            if contig not in sequences or not 0 <= start < end <= len(sequences[contig]):
                raise ValueError(f'BED line {number} is outside the exact reference; expected zero-based, half-open coordinates.')
            intervals.append({'contig': contig, 'start': start + 1, 'end': end,
                              'label': values[3] if len(values) > 3 else f'BED-line-{number}'})
            if len(intervals) > 100_000:
                raise ValueError('BED mask exceeds the 100,000 interval safety bound.')
    return intervals


def _index_intervals(entries):
    groups = defaultdict(list)
    for entry in entries:
        groups[entry['contig']].append(entry)
    indexed = {}
    for contig, rows in groups.items():
        rows.sort(key=lambda row: (row['start'], row['end']))
        maximum, prefix = 0, []
        for row in rows:
            maximum = max(maximum, row['end'])
            prefix.append(maximum)
        indexed[contig] = rows, [row['start'] for row in rows], prefix
    return indexed


def _overlaps(index, contig, position):
    rows, starts, prefix = index.get(contig, ([], [], []))
    i = bisect.bisect_right(starts, position) - 1
    matches = []
    while i >= 0 and prefix[i] >= position:
        if rows[i]['end'] >= position:
            matches.append(rows[i])
        i -= 1
    return matches


def eligible_reference_intervals(sequences, features, masks, *, coding_only=False):
    """Merged CDS/all-reference minus BED masks, zero-based half-open intervals.

    These are protocol-eligible coordinates, not necessarily callable sites.
    The distance engine must still intersect them with actual per-sample calls.
    """
    def merge(rows):
        merged = []
        for start, end in sorted(rows):
            if merged and start <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], end)
            else:
                merged.append([start, end])
        return merged

    output, feature_groups, mask_groups = [], defaultdict(list), defaultdict(list)
    for entry in features:
        feature_groups[entry['contig']].append((entry['start'] - 1, entry['end']))
    for entry in masks:
        mask_groups[entry['contig']].append((entry['start'] - 1, entry['end']))
    for contig, sequence in sorted(sequences.items()):
        included = merge(feature_groups[contig] if coding_only else [(0, len(sequence))])
        excluded = merge(mask_groups[contig])
        cursor = 0
        for start, end in included:
            while cursor < len(excluded) and excluded[cursor][1] <= start:
                cursor += 1
            index, current = cursor, start
            while index < len(excluded) and excluded[index][0] < end:
                left, right = excluded[index]
                if left > current:
                    output.append({'contig': contig, 'start': current, 'end': min(end, left)})
                current = max(current, right)
                if current >= end:
                    break
                index += 1
            if current < end:
                output.append({'contig': contig, 'start': current, 'end': end})
    return output


def _consequence(feature, position, alternate):
    result = {key: feature[key] for key in ('id', 'gene', 'locus_tag', 'product', 'start', 'end', 'strand', 'origin', 'annotation_binding')}
    result.update(effect='unresolved', confidence='coding_context_only', reason='',
                  genetic_code=feature['genetic_code'], protein_effect='not_inferred')
    if feature['compound']:
        return dict(result, reason='Multipart/compound CDS is not interpreted by this bounded consequence engine.')
    if feature['partial'] or feature['phase'] != 0:
        return dict(result, reason='Partial CDS or nonzero initial phase: coding context retained without claiming a complete codon frame.')
    offset = position - feature['start'] if feature['strand'] == '+' else feature['end'] - position
    sequence = feature['sequence']
    start = (offset // 3) * 3
    codon = sequence[start:start + 3]
    if len(codon) != 3 or set(codon) - set('ACGT'):
        return dict(result, reason='Reference codon is incomplete or contains an ambiguous base.')
    changed = alternate if feature['strand'] == '+' else reverse_complement(alternate)
    alternative = codon[:offset % 3] + changed + codon[offset % 3 + 1:]
    code = feature['genetic_code']

    def translate(triplet):
        if start == 0 and triplet in _STARTS[code]:
            return 'M'
        return 'W' if code == 4 and triplet == 'TGA' else _CODONS[triplet]

    aa, alt_aa = translate(codon), translate(alternative)
    effect = ('start_lost' if start == 0 and codon in _STARTS[code] and alternative not in _STARTS[code]
              else 'stop_gained' if aa != '*' and alt_aa == '*' else 'stop_lost' if aa == '*' and alt_aa != '*'
              else 'synonymous' if aa == alt_aa else 'missense')
    return dict(result, effect=effect, confidence='predicted_CDS' if feature['origin'] == 'predicted_Prodigal' else feature['annotation_binding'],
                cds_position=offset + 1, codon_number=offset // 3 + 1, ref_codon=codon, alt_codon=alternative,
                ref_amino_acid=aa, alt_amino_acid=alt_aa,
                reason='Single-variant coding consequence only; effects of nearby variants are not haplotype-phased.')


def annotate_snps(reference_path, variants, *, annotation_path=None, annotation_reference_sha256=None,
                   mask_bed=None, coding_only=False, predict_coding=True, genetic_code=11,
                   max_sites=100_000, cancelled=None, progress=None):
    """Annotate retained and excluded SNPs against one exact, hashed reference.

    Variants have ``contig``, one-based ``pos``, ``ref``, and ``alt`` (one base,
    comma-delimited bases, or a list). Indels/symbolic alleles are retained as
    unsupported rows, never converted to SNPs. Plain GFF3 binding is explicitly
    unverified unless the caller supplies a matching reference SHA attestation.
    """
    check_cancelled(cancelled)
    if isinstance(max_sites, bool) or not isinstance(max_sites, int) or not 1 <= max_sites <= 1_000_000:
        raise ValueError('Maximum sites must be an integer from 1 to 1,000,000.')
    reference = Path(reference_path).resolve()
    tracked = [reference] + ([Path(annotation_path).resolve()] if annotation_path else []) + ([Path(mask_bed).resolve()] if mask_bed else [])
    signatures = [file_signature(path) for path in tracked]
    reference_sha = file_sha256(reference, cancelled)
    sequences = _sequences(reference, cancelled)
    if annotation_path is not None:
        features, annotation = _gff_features(annotation_path, sequences, reference_sha256=reference_sha,
                                            annotation_reference_sha256=annotation_reference_sha256, cancelled=cancelled)
    elif predict_coding:
        genes, provenance = predict_cds(reference, genetic_code=genetic_code, cancelled=cancelled, progress=progress)
        features = [dict(gene, phase=0, genetic_code=gene['cds_qc']['genetic_code'], partial=gene['cds_qc']['partial_begin'] or gene['cds_qc']['partial_end'],
                         gene='', locus_tag='', product='', origin='predicted_Prodigal', annotation_binding='exact_reference_prediction', compound=False)
                    for gene in genes]
        annotation = {'kind': 'predicted_Prodigal', 'binding': 'exact_reference_prediction', 'cds_count': len(features),
                      'provenance': provenance, 'reference_sha256': reference_sha}
    else:
        features, annotation = [], {'kind': 'not_run', 'binding': 'none', 'cds_count': 0}
    if coding_only and (annotation['kind'] == 'not_run' or annotation['binding'] == 'coordinate_only_unverified'):
        raise ValueError('Coding-only filtering requires a reference-bound CDS annotation or explicit Prodigal prediction. Coordinate-only GFF3 is unverified.')
    intervals = _bed_intervals(mask_bed, sequences, cancelled) if mask_bed else []
    eligible = eligible_reference_intervals(sequences, features, intervals, coding_only=coding_only)
    feature_index, mask_index = _index_intervals(features), _index_intervals(intervals)
    rows = []
    for index, variant in enumerate(variants):
        check_cancelled(cancelled)
        if index >= max_sites:
            raise ValueError('SNP input exceeds the explicit site bound; no truncated report was accepted.')
        contig, position, ref = variant.get('contig'), variant.get('pos'), str(variant.get('ref') or '').upper()
        if contig not in sequences or isinstance(position, bool) or not isinstance(position, int) or not 1 <= position <= len(sequences[contig]):
            raise ValueError(f'Variant {index + 1} has an unknown contig or invalid one-based reference position.')
        if not ref or sequences[contig][position - 1:position - 1 + len(ref)] != ref:
            raise ValueError(f'Variant {index + 1} REF does not exactly match the selected reference FASTA.')
        alternatives = variant.get('alt')
        alternatives = alternatives.split(',') if isinstance(alternatives, str) else alternatives
        if not isinstance(alternatives, (list, tuple)) or not alternatives or any(not isinstance(value, str) or not value for value in alternatives):
            raise ValueError(f'Variant {index + 1} requires at least one explicit ALT allele.')
        if len(alternatives) > 20:
            raise ValueError('Variant exceeds the explicit 20 alternate allele bound.')
        overlaps, masks = _overlaps(feature_index, contig, position), _overlaps(mask_index, contig, position)
        annotations, excluded = [], []
        if masks:
            excluded.append('explicit_BED_mask')
        if coding_only and not overlaps:
            excluded.append('outside_selected_CDS_annotation')
        for alternate in alternatives:
            alternate = alternate.upper()
            if len(ref) != 1 or ref not in 'ACGT' or len(alternate) != 1 or alternate not in 'ACGT':
                annotations.append({'alt': alternate, 'status': 'unsupported_non_SNP', 'consequences': []})
                continue
            if alternate == ref:
                annotations.append({'alt': alternate, 'status': 'reference_allele', 'consequences': []})
                continue
            effects = [_consequence(feature, position, alternate) for feature in overlaps]
            annotations.append({'alt': alternate, 'status': 'annotated' if overlaps else 'noncoding' if annotation['kind'] != 'not_run' else 'annotation_not_run',
                                'consequences': effects})
        if all(item['status'] in {'unsupported_non_SNP', 'reference_allele'} for item in annotations):
            excluded.append('not_a_supported_SNP')
        rows.append({'variant': copy.deepcopy(variant), 'contig': contig, 'pos': position, 'ref': ref,
                     'coding_context': bool(overlaps), 'included': not excluded, 'exclusion_reasons': excluded,
                     'mask_intervals': masks, 'alternatives': annotations})
        if progress and index % 1000 == 0:
            progress(index + 1, 0, f'Annotating reference SNP site {index + 1:,}')
    if signatures != [file_signature(path) for path in tracked]:
        raise ValueError('Reference, annotation or BED mask changed during annotation; no result was accepted.')
    protocol = {'reference_sha256': reference_sha, 'annotation': annotation,
                'mask': {'status': 'applied' if mask_bed else 'not_run', 'sha256': file_sha256(mask_bed, cancelled) if mask_bed else None,
                         'path': str(Path(mask_bed).resolve()) if mask_bed else None, 'interval_count': len(intervals),
                         'coordinate_system': 'BED zero-based half-open; reported intervals one-based inclusive',
                         'interpretation': 'User-supplied intervals; repeats/IS identity is not inferred or independently validated.'},
                'coding_only': bool(coding_only), 'recombination_masking': 'not_run', 'max_sites': max_sites,
                'eligible_reference_intervals': eligible, 'eligible_intervals_coordinate_system': 'zero-based half-open',
                'eligible_reference_bases': sum(row['end'] - row['start'] for row in eligible),
                'callability_policy': 'Protocol-eligible intervals must additionally intersect actual per-sample callable bases; eligibility is not callability.',
                'single_variant_codon_effects': True, 'noncoding_exclusion_is_not_functional_classification': True,
                'method_sources': ['https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi',
                                   'https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md',
                                   'https://genome.ucsc.edu/FAQ/FAQformat.html#format1'],
                'software': 'WMLSTudio', 'version': __version__}
    protocol_sha = hashlib.sha256(json.dumps(protocol, sort_keys=True, separators=(',', ':')).encode()).hexdigest()
    if signatures != [file_signature(path) for path in tracked]:
        raise ValueError('Reference, annotation or BED mask changed while hashing the protocol; no result was accepted.')
    return {'format_version': 1, 'status': 'completed', 'reference_path': str(reference), 'reference_sha256': reference_sha,
            'protocol': protocol, 'protocol_sha256': protocol_sha, 'sites': rows,
            'summary': {'sites': len(rows), 'included': sum(row['included'] for row in rows),
                        'excluded': sum(not row['included'] for row in rows)},
            'limitations': ['Coding consequences are sequence-level annotations, not proof of functional effect, pathogenicity or antimicrobial phenotype.',
                            'Noncoding variants may alter regulatory function; coding-only filtering can remove biologically important signals.',
                            'Repeat/IS BED masking is independent of coding selection. Coding genes may lie in repeats or mobile elements.',
                            'Recombination has not been inferred or masked by this service. No Gubbins-equivalent analysis is claimed.',
                            'SNP distances and annotations do not establish direct transmission, direction or source attribution.',
                            'Single-variant codon effects are not phased with nearby variants; partial and compound CDS consequences remain unresolved.']}
