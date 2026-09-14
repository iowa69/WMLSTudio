"""Fetch a pinned public practice cohort on request; never invoked at app startup.

Sequence data is never committed to this repository or shipped in the portable ZIP.
The cohorts are teaching material: no expected ST, cluster or threshold is supplied,
and a folder the app creates from an organism proposal is a filing decision, not a
laboratory identification.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from wmlstudio.practice_cohorts import (
    cohort_names,
    describe_cohorts,
    download_cohort,
    verify_cohort,
    verify_pins,
)


def _report_progress(index, total, message):
    print(f'[{index}/{total}] {message}', flush=True)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--list', action='store_true', help='Describe the available cohorts')
    parser.add_argument('--cohort', choices=cohort_names(), help='Cohort to fetch or check')
    parser.add_argument('--destination', type=Path,
                        help='Directory to create for the cohort files')
    parser.add_argument('--cache-root', type=Path,
                        help='Shared content-addressed cache; defaults to _cache beside the cohort')
    parser.add_argument('--source-root', type=Path,
                        help='Optional local mirror of the NCBI genomes/all tree')
    parser.add_argument('--verify-pins', action='store_true',
                        help='Compare pinned MD5s with NCBI published checksums, no genomes')
    parser.add_argument('--verify', type=Path,
                        help='Re-hash an installed cohort directory against the pinned literals')
    parser.add_argument('--quiet', action='store_true', help='Do not print per-genome progress')
    args = parser.parse_args(argv)
    progress = None if args.quiet else _report_progress
    if args.list:
        print(json.dumps(describe_cohorts(), indent=2))
        return 0
    if args.verify is not None:
        manifest = verify_cohort(args.verify)
        print(json.dumps({'path': str(Path(args.verify).resolve()), 'cohort': manifest['cohort'],
                          'cohort_digest': manifest['cohort_digest'],
                          'genomes': len(manifest['genomes']),
                          'total_bytes': manifest['total_bytes']}, indent=2))
        return 0
    if args.cohort is None:
        parser.error('Choose --cohort, --list or --verify.')
    if args.verify_pins:
        report = verify_pins(args.cohort, source_root=args.source_root, progress=progress)
        print(json.dumps(report, indent=2))
        return 1 if report['stale'] else 0
    if args.destination is None:
        parser.error('--destination is required when fetching a cohort.')
    result = download_cohort(args.cohort, args.destination, cache_root=args.cache_root,
                             source_root=args.source_root, progress=progress)
    print(json.dumps(result, indent=2))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
