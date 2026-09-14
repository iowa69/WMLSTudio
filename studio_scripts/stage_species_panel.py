"""Stage the explicit pinned broad species panel; never called at startup."""

from __future__ import annotations

import argparse
import gzip
import json
import os
import shutil
import tempfile
from pathlib import Path

from wmlstudio.characterization_refs import validate_characterization_references
from wmlstudio.organism_panel import (
    SPECIES_PANEL,
    provision_species_panel,
    verify_pins,
)

ROOT = Path(__file__).resolve().parents[1]


def stage_snapshot(source, destination):
    source, destination = Path(source).resolve(), Path(destination).resolve()
    manifest = validate_characterization_references(source)
    if source == destination:
        return manifest
    if destination.exists():
        existing = validate_characterization_references(destination)
        if existing['reference_digest'] != manifest['reference_digest']:
            raise ValueError('Refusing to replace a different installed species panel; stage to a new destination.')
        return existing
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='.species-panel-stage-', dir=destination.parent) as temporary:
        stage = Path(temporary) / 'snapshot'
        stage.mkdir()
        for relative in ['manifest.json', *(entry['path'] for entry in manifest['files'])]:
            target = stage / relative
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(source / relative, target)
        validate_characterization_references(stage)
        os.replace(stage, destination)
    return manifest


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--destination', type=Path,
                        help='Where to place the verified snapshot; defaults to a temporary install under --root')
    parser.add_argument('--root', type=Path, default=ROOT / 'artifacts' / 'species-panel',
                        help='Directory that receives species-panel-<digest> when no --destination is given')
    parser.add_argument('--source-root', type=Path,
                        help='Optional local mirror of the pinned NCBI genomes/all tree')
    parser.add_argument('--snapshot', type=Path,
                        help='Existing verified immutable species panel, staged without downloading')
    parser.add_argument('--verify-pins', action='store_true',
                        help='Only re-read the published md5checksums.txt for every pinned assembly')
    args = parser.parse_args(argv)
    if args.verify_pins:
        report = verify_pins(source_root=args.source_root)
        print(json.dumps({'pinned': len(SPECIES_PANEL), **report}, indent=2))
        return 0 if not report['missing'] else 1
    source = args.snapshot
    if source is None:
        source = Path(provision_species_panel(args.root, source_root=args.source_root)['path'])
    manifest = stage_snapshot(source, args.destination) if args.destination else \
        validate_characterization_references(source)
    destination = Path(args.destination or source).resolve()
    uncompressed = 0
    for entry in manifest['files']:
        path = destination / entry['path']
        if path.suffix == '.gz':
            with gzip.open(path, 'rb') as handle:
                while chunk := handle.read(1024 * 1024):
                    uncompressed += len(chunk)
        else:
            uncompressed += entry['bytes']
    print(json.dumps({'path': str(destination), 'reference_digest': manifest['reference_digest'],
                      'species_references': len(manifest['species']),
                      'taxa': [f"{entry['genus']} {entry['species']}".strip() for entry in manifest['species']],
                      'stored_bytes': sum(entry['bytes'] for entry in manifest['files']),
                      'uncompressed_reference_bytes': uncompressed}, indent=2))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
