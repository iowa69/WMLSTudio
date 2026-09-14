"""Registry mapping an organism match rule to the organism-specific assays offered.

The match rule never gates execution. It records whether a module's reference
panel was curated for this isolate's organism, so a result run outside those
taxa is reported as such instead of being quietly presented as validated. No
module assigns metadata['organism']; organism assignment stays an explicit user
action through storage.assign_organism.
"""

from __future__ import annotations

import json
import threading
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable

# Reserved by characterize_assembly itself; a module may never shadow one.
CORE_SECTIONS = frozenset({'species_evidence', 'virulence', 'drug_associations', 'plasmid_hypotheses'})
APPLICABILITY_ORDER = ('recommended', 'possible', 'unknown_organism', 'off_panel')
_ASSAY_MODULES = ('sccmec_evidence', 'klebsiella_evidence')
_LOCK = threading.Lock()
_LOADED = False


@dataclass(frozen=True)
class OrganismMatch:
    """Three-state applicability: curated for, plausible for, or outside the panel."""

    genera: frozenset
    species: frozenset = frozenset()
    exclude_species: frozenset = frozenset()

    def evaluate(self, genus, species):
        """Pure; returns (state, reason). Never gates execution, only recommendation."""
        name, epithet = (genus or '').strip(), (species or '').strip()
        folded, folded_species = name.casefold(), epithet.casefold()
        if not folded and not folded_species:
            return 'unknown_organism', 'No genus or species is assigned or detected for this isolate.'
        if folded not in {value.casefold() for value in self.genera}:
            return 'off_panel', f'{name or "Unknown"} is outside the taxa this reference panel was curated for.'
        if folded_species and folded_species in {value.casefold() for value in self.exclude_species}:
            return 'off_panel', f'{name} {epithet} is explicitly excluded from this panel.'
        if self.species and folded_species not in {value.casefold() for value in self.species}:
            return 'possible', ('The genus matches, but the species is outside the taxa this panel was curated '
                                'against. Results are reported as outside the primary validation taxa.')
        return 'recommended', 'Genus and species match the taxa this reference panel was curated for.'


@dataclass(frozen=True)
class OrganismModule:
    """One offerable assay: what it needs, how it runs and how it is presented."""

    key: str
    title: str
    column_title: str
    match: OrganismMatch
    runner: Callable
    option_keys: tuple
    manifest_sections: tuple
    summary: Callable
    detail_html: Callable
    report_default: bool = True
    locus_st_provider: Callable | None = None
    limitations: tuple = field(default_factory=tuple)


REGISTRY: dict[str, OrganismModule] = {}


def register(module):
    """Refuse a key that would shadow a core characterization section or another module."""
    if not isinstance(module, OrganismModule) or not module.key:
        raise ValueError('An organism module must be an OrganismModule with a key.')
    if module.key in CORE_SECTIONS:
        raise ValueError(f'Organism module key {module.key!r} collides with a core characterization section.')
    if module.key in REGISTRY and REGISTRY[module.key] is not module:
        raise ValueError(f'Organism module key {module.key!r} is already registered.')
    unsupported = set(module.option_keys) - {'threads', 'blastn_path', 'makeblastdb_path'}
    if unsupported:
        raise ValueError(f'Organism module {module.key!r} requests unsupported options: {sorted(unsupported)}.')
    REGISTRY[module.key] = module
    return module


def _load():
    """Import the assay modules lazily so importing the registry stays cheap."""
    global _LOADED
    with _LOCK:
        if _LOADED:
            return REGISTRY
        _LOADED = True
    for name in _ASSAY_MODULES:
        __import__(f'{__package__}.{name}')
    return REGISTRY


def registered_modules():
    return dict(_load())


def organism_of(record):
    """Assigned organism first, local detection second; mirrors ui_common.organism_for."""
    metadata = record.get('metadata') or {}
    assigned = metadata.get('organism') or {}
    if isinstance(assigned, str):
        parts = assigned.split(maxsplit=1)
        assigned = {'genus': parts[0] if parts else '', 'species': parts[1] if len(parts) > 1 else ''}
    if assigned.get('genus') or assigned.get('species'):
        return assigned.get('genus', ''), assigned.get('species', ''), 'Assigned'
    identification = (record.get('result') or {}).get('identification') or {}
    detected = identification.get('organism') or {}
    genus = detected.get('genus') or identification.get('genus') or ''
    species = detected.get('species') or identification.get('species') or ''
    return genus, species, 'Provisional' if genus else 'Unknown'


def modules_for(genus, species):
    """Every registered module with its applicability, most applicable first."""
    rows = [(module.match.evaluate(genus, species), module) for module in _load().values()]
    return sorted(((state, reason, module) for (state, reason), module in rows),
                  key=lambda row: (APPLICABILITY_ORDER.index(row[0]), row[2].key))


def cohort_applicability(module, records):
    """Counter of applicability states over a selected cohort; never a filter."""
    states = Counter()
    for record in records:
        genus, species, _ = organism_of(record)
        states[module.match.evaluate(genus, species)[0]] += 1
    return states


def _manifest(reference_root):
    if reference_root is None:
        return None
    path = Path(reference_root) / 'manifest.json'
    try:
        if path.stat().st_size > 2 * 1024 * 1024:
            return None
        manifest = json.loads(path.read_text(encoding='utf-8'))
    except (OSError, ValueError):
        return None
    return manifest if isinstance(manifest, dict) else None


def _section_present(manifest, section):
    node = manifest
    for part in section.split('.'):
        if not isinstance(node, dict) or not node.get(part):
            return False
        node = node[part]
    return True


def module_tasks(selected, reference_root, *, threads=2, blastn_path=None, makeblastdb_path=None):
    """Assay tasks for the selected modules in characterize_assembly's task shape.

    A selected module whose reference sections are absent is reported as not run
    with the missing section named, so an out-of-date reference snapshot is
    diagnosable from the stored evidence instead of looking like a negative.
    """
    available = {'threads': threads, 'blastn_path': blastn_path, 'makeblastdb_path': makeblastdb_path}
    selected = selected or {}
    manifest = _manifest(reference_root)
    tasks = []
    for key, module in _load().items():
        options = {name: available[name] for name in module.option_keys}
        if not selected.get(key):
            tasks.append((key, False, module.runner, options, 'Assay not selected.'))
            continue
        missing = [section for section in module.manifest_sections if not _section_present(manifest or {}, section)]
        if manifest is None:
            tasks.append((key, False, module.runner, options,
                          'No readable characterization reference manifest is installed; install or select a verified snapshot first.'))
        elif missing:
            tasks.append((key, False, module.runner, options,
                          f'This reference snapshot (format {manifest.get("format_version")}) carries no '
                          f'{", ".join(missing)} section for {module.title}. Install an updated snapshot from '
                          'Characterization, Install / update.'))
        else:
            tasks.append((key, True, module.runner, options, None))
    return tasks


def stamp_applicability(result, organism):
    """Record, never gate: every module block says which taxa its panel covers."""
    genus, species = (organism or ('', ''))[:2]
    for key, module in _load().items():
        block = result.get(key)
        if isinstance(block, dict):
            block['applicability'], block['applicability_reason'] = module.match.evaluate(genus, species)
    return result


def summarize_record(evidence):
    """One flat line of organism-specific calls for exports and at-a-glance rows."""
    parts = []
    for key, module in _load().items():
        block = (evidence or {}).get(key)
        if not isinstance(block, dict) or block.get('status') in {None, 'not_run'}:
            continue
        text = module.summary(block)
        if text:
            parts.append(text)
    return '; '.join(parts)
