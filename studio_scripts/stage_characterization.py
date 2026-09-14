"""Stage an explicit pinned species/virulence starter; never called at startup."""

from __future__ import annotations

import argparse
import gzip
import json
import os
import shutil
import tempfile
from pathlib import Path

from wmlstudio.characterization_refs import (
    provision_characterization_references,
    validate_characterization_references,
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
            raise ValueError('Refusing to replace a different installed characterization snapshot; stage to a new destination.')
        return existing
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='.characterization-stage-', dir=destination.parent) as temporary:
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
    parser.add_argument('--destination', type=Path, default=ROOT / 'src/wmlstudio/resources/characterization/starter')
    parser.add_argument('--source-root', type=Path, help='Optional local checkout of the pinned official Kleborate source')
    parser.add_argument('--snapshot', type=Path, help='Existing verified immutable characterization snapshot, without downloading')
    args = parser.parse_args(argv)
    with tempfile.TemporaryDirectory(prefix='wmlstudio-characterization-stage-') as temporary:
        source = args.snapshot
        if source is None:
            source = Path(provision_characterization_references(temporary, source_root=args.source_root)['path'])
        manifest = stage_snapshot(source, args.destination)
    uncompressed = 0
    for entry in manifest['files']:
        path = args.destination / entry['path']
        if path.suffix == '.gz':
            with gzip.open(path, 'rb') as handle:
                while chunk := handle.read(1024 * 1024):
                    uncompressed += len(chunk)
        else:
            uncompressed += entry['bytes']
    print(json.dumps({'path': str(args.destination.resolve()), 'reference_digest': manifest['reference_digest'],
                      'species_references': len(manifest['species']), 'virulence_loci': list(manifest['virulence']),
                      'stored_bytes': sum(entry['bytes'] for entry in manifest['files']),
                      'uncompressed_reference_bytes': uncompressed}, indent=2))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
