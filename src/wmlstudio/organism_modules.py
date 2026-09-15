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
from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import Callable

# Reserved by characterize_assembly itself; a module may never shadow one.
CORE_SECTIONS = frozenset({'species_evidence', 'virulence', 'drug_associations', 'plasmid_hypotheses'})
APPLICABILITY_ORDER = ('recommended', 'possible', 'unknown_organism', 'off_panel')
# Said wherever a user chooses or reads these assays, not only in a docstring.
BOUNDARY = ('Each organism-specific tool is a WMLSTudio BLAST or exact-allele screen over pinned public '
            'reference data, with the source repository, revision and licence recorded in the snapshot. '
            'None of them is Kleborate, Kaptive, AMRFinderPlus or SCCmecFinder, and none is equivalent to '
            'those tools: no susceptibility, phenotype, serotype or transmission conclusion is produced.')
APPLICABILITY_TITLES = {'recommended': 'within the curated taxa',
                        'possible': 'same genus, other species',
                        'unknown_organism': 'no organism assigned yet',
                        'off_panel': 'outside these taxa'}
# Why an assay produced nothing. These are four different answers and a column
# that showed the same 'not_run' token for all of them is what leaves a user
# hunting for output that was never asked for. The strings below are the ones
# actually written into a result, kept here so a table cell can say which kind
# of blank it is without re-deriving it from prose.
NOT_SELECTED = 'Assay not selected.'
NO_MANIFEST = ('No readable characterization reference manifest is installed; install or select a '
               'verified snapshot first.')
# characterize_assembly substitutes this one when no snapshot is chosen at all;
# test_organism_modules pins it, so a change upstream fails there rather than
# silently turning a missing reference into an unexplained blank.
NO_SNAPSHOT = 'Install or select a verified characterization reference snapshot first.'
SKIPPED_PREFIX = 'Not run for this isolate.'
# How every "the installed snapshot does not carry this" sentence opens, whether
# it is written here or by an assay that read the manifest itself.
SNAPSHOT_PREFIX = 'This reference snapshot'
BLANK_SUMMARIES = {'not_recorded': 'not recorded',
                   'not_selected': 'not run: not selected',
                   'reference_missing': 'not run: reference not installed',
                   'off_panel': 'not run: outside these taxa',
                   'unrecorded': 'not run'}
BLANK_SENTENCES = {
    'not_recorded': ('No result is recorded for this assay on this isolate: it was never run, or its result '
                     'belongs to an earlier assembly. Nothing was screened, so this is not a negative result.'),
    'not_selected': ('This assay was not selected for the run, so nothing was screened for this isolate. '
                     'That is not a negative result.'),
    'reference_missing': ('This assay could not run because the reference data it needs is not installed. '
                          'Install or update the characterization reference snapshot and run it again. '
                          'Nothing was screened, so this is not a negative result.'),
    'off_panel': ('This assay was not run for this isolate: it is outside the taxa its reference panel covers. '
                  'Nothing was screened, so this is not a negative result.'),
    'unrecorded': ('This assay produced no result and did not record why. Nothing was screened, so this is not '
                   'a negative result.'),
}
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
    # Describes this assay's own results only. Registration wraps it so that a
    # blank — no result recorded, or a not-run one — is answered identically for
    # every assay, saying which kind of blank it is rather than 'not_run'.
    summary: Callable
    detail_html: Callable
    report_default: bool = True
    # One plain sentence, for someone who is not a bioinformatician: what a result
    # from this assay establishes, and what it explicitly does not.
    purpose: str = ''
    locus_st_provider: Callable | None = None
    limitations: tuple = field(default_factory=tuple)


REGISTRY: dict[str, OrganismModule] = {}
_ORDER: dict[str, int] = {}


def blank_kind(block):
    """Which kind of blank a module result is, or None when it is a result.

    Four blanks that used to read alike: nothing recorded at all, an assay the
    run did not select, an assay whose reference data is not installed, and an
    isolate outside the taxa the panel covers. A user who cannot tell them apart
    reads every one of them as 'we looked and found nothing'.
    """
    if not isinstance(block, dict) or not block:
        return 'not_recorded'
    if block.get('status') != 'not_run':
        return None
    if block.get('not_run_kind') in BLANK_SUMMARIES:
        return block['not_run_kind']
    reason = str(block.get('reason') or '')
    if reason.startswith(SKIPPED_PREFIX):
        return 'off_panel'
    if reason == NOT_SELECTED:
        return 'not_selected'
    if reason in {NO_MANIFEST, NO_SNAPSHOT} or reason.startswith(SNAPSHOT_PREFIX):
        return 'reference_missing'
    return 'unrecorded'


def blank_summary(block):
    """The short phrase a table cell shows for a blank, or '' when there is a result."""
    kind = blank_kind(block)
    return BLANK_SUMMARIES[kind] if kind else ''


def blank_sentence(block):
    """Why there is no result, as a sentence a reader can act on, or '' for a result."""
    kind = blank_kind(block)
    if not kind:
        return ''
    reason = str(block.get('reason') or '').strip() if isinstance(block, dict) else ''
    if kind in {'off_panel', 'unrecorded'} and reason:
        return reason
    # The recorded reason is repeated only where it names what the headline
    # cannot: which section of the installed snapshot is absent. 'Install a
    # snapshot first' twice in one paragraph helps nobody.
    if kind == 'reference_missing' and reason.startswith(SNAPSHOT_PREFIX):
        return f'{BLANK_SENTENCES[kind]} {reason}'
    return BLANK_SENTENCES[kind]


def _blank_aware(summary):
    """Answer blanks here so every assay says which kind of blank it is.

    Each module describes its own results; none of them should have to reason
    about the four ways a result can be absent, and none of them should answer
    that question differently from the others.
    """
    def summarize(block):
        return blank_summary(block) or summary(block)

    return summarize


def register(module):
    """Refuse a key that would shadow a core characterization section or another module.

    The stored module is a copy whose ``summary`` answers blanks uniformly; it
    is also what is returned, so the registry and the caller hold one object.
    """
    if not isinstance(module, OrganismModule) or not module.key:
        raise ValueError('An organism module must be an OrganismModule with a key.')
    if module.key in CORE_SECTIONS:
        raise ValueError(f'Organism module key {module.key!r} collides with a core characterization section.')
    if module.key in REGISTRY and REGISTRY[module.key] is not module:
        raise ValueError(f'Organism module key {module.key!r} is already registered.')
    unsupported = set(module.option_keys) - {'threads', 'blastn_path', 'makeblastdb_path'}
    if unsupported:
        raise ValueError(f'Organism module {module.key!r} requests unsupported options: {sorted(unsupported)}.')
    module = replace(module, summary=_blank_aware(module.summary))
    _ORDER.setdefault(module.key, len(_ORDER))
    REGISTRY[module.key] = module
    return module


def _declared_order(item):
    """Order modules by _ASSAY_MODULES, then by their order within their own file.

    Whichever assay module Python happens to import first must not decide the
    column order a user sees, or the order of the columns in an exported table.
    """
    key, module = item
    origin = (getattr(module.runner, '__module__', '') or '').rsplit('.', 1)[-1]
    position = _ASSAY_MODULES.index(origin) if origin in _ASSAY_MODULES else len(_ASSAY_MODULES)
    return (position, _ORDER.get(key, len(_ORDER)))


def _load():
    """Import the assay modules lazily so importing the registry stays cheap."""
    global _LOADED
    with _LOCK:
        if _LOADED:
            return REGISTRY
        _LOADED = True
    for name in _ASSAY_MODULES:
        __import__(f'{__package__}.{name}')
    # An assay already imported by another path registered ahead of its turn, so
    # settle the order once here rather than leaving it to import timing.
    ordered = sorted(REGISTRY.items(), key=_declared_order)
    REGISTRY.clear()
    REGISTRY.update(ordered)
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


def _isolates(count):
    return f'{count} isolate' if count == 1 else f'{count} isolates'


def cohort_sentence(module, counts, *, off_panel_included=False):
    """Name, in plain words, exactly which selected isolates this assay will and will not run on.

    A mixed-genus cohort is the normal case in an outbreak workspace, so the
    dialog states the split rather than leaving the user to infer it from a
    badge. Skipping an off-panel isolate records nothing for it, and saying so
    is what keeps the blank from reading as a negative result.
    """
    total = sum(counts.values())
    off_panel = counts.get('off_panel', 0)
    running = total if off_panel_included else total - off_panel
    detail = [f'{counts[state]} {APPLICABILITY_TITLES[state]}' for state in APPLICABILITY_ORDER
              if counts.get(state) and (state != 'off_panel' or off_panel_included)]
    sentence = f'Will run on {_isolates(running)} of {total} selected'
    if detail:
        sentence += ' (' + '; '.join(detail) + ')'
    if not off_panel:
        return sentence + '.'
    if off_panel_included:
        return (sentence + f'. The {_isolates(off_panel)} outside these taxa are included at your request and '
                'their results are labelled off-panel, which is weak evidence in either direction.')
    return (sentence + f'. Will not run on the {_isolates(off_panel)} outside these taxa; nothing is recorded '
            'for them, which is not a negative result.')


def selection_for_record(record, selected, *, include_off_panel=False):
    """Split one isolate's chosen assays into those that run and those skipped, with the reason.

    The registry never gates on the match rule; this is the caller's choice, made
    once here so the plan dialog and the command line make it the same way.
    """
    genus, species, _ = organism_of(record)
    running, skipped = {}, {}
    for key, module in _load().items():
        if not (selected or {}).get(key):
            continue
        state, reason = module.match.evaluate(genus, species)
        if state == 'off_panel' and not include_off_panel:
            skipped[key] = reason
        else:
            running[key] = True
    return running, skipped


def record_skipped(result, skipped):
    """Replace 'Assay not selected.' with why this isolate in particular was skipped."""
    for key, reason in (skipped or {}).items():
        block = result.get(key)
        if isinstance(block, dict) and block.get('status') == 'not_run':
            block['reason'] = (SKIPPED_PREFIX + ' ' + reason + ' Nothing was screened, so this is not '
                               'a negative result; include isolates outside a tool\'s reference taxa to run it '
                               'anyway.')
            block['not_run_kind'] = 'off_panel'
    return result


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
            tasks.append((key, False, module.runner, options, NOT_SELECTED))
            continue
        missing = [section for section in module.manifest_sections if not _section_present(manifest or {}, section)]
        if manifest is None:
            tasks.append((key, False, module.runner, options, NO_MANIFEST))
        elif missing:
            tasks.append((key, False, module.runner, options,
                          f'This reference snapshot (format {manifest.get("format_version")}) carries no '
                          f'{", ".join(missing)} section for {module.title}. Install an updated snapshot from '
                          'Characterization, Install / update.'))
        else:
            tasks.append((key, True, module.runner, options, None))
    return tasks


def stamp_applicability(result, organism):
    """Record, never gate: every module block says which taxa its panel covers.

    A block with no result also records which kind of blank it is, so a table
    cell and an exported row can say whether the assay was not selected, could
    not run for want of a reference, or was skipped for this isolate's organism
    without any of them having to read the sentence back.
    """
    genus, species = (organism or ('', ''))[:2]
    for key, module in _load().items():
        block = result.get(key)
        if isinstance(block, dict):
            block['applicability'], block['applicability_reason'] = module.match.evaluate(genus, species)
            kind = blank_kind(block)
            if kind and kind != 'not_recorded':
                block['not_run_kind'] = kind
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
