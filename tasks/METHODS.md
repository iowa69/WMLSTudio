# WMLSTudio methods and build record

## 2026-09-12

Host: Linux x86_64. Workspace: `/home/iowa/Desktop/WMLSTudio`.
Original HEAD: `6cad46f` from iowa69/WMLST. Legacy sources are already deleted.
References: https://github.com/iowa69/mlstudio and https://github.com/iowa69/WMLST.
Qt deployment: https://doc.qt.io/qtforpython-6/deployment/index.html.
Competitor baseline: https://www.ridom.de/seqsphere/ug/v105/User_Guide.html.

The new analysis engine will use exact nucleotide matches on both strands. Unknown,
ambiguous, mixed and missing states must stay distinct from assigned alleles. Read
QC must say whether it sampled records or scanned the complete file. MST edges must
include their shared-locus denominators; insufficient overlap is not zero distance.

Existing sequence data and databases remain read only. Local validation outputs go
to `results/2026-09-12_validation`. Synthetic fixtures carry explicit labels.
No study-level biological inference or patient interpretation is intended.

## Environment setup and commands

Installed isolated uv 0.10.11 from its official GitHub release archive into `/tmp`,
then `uv sync --python /usr/bin/python3`. No base conda/system environment changes.
`uv.lock` pins PySide6-Essentials 6.11.2, Shiboken 6.11.2, pyahocorasick 2.3.1,
PyInstaller 6.22.2, pytest 9.1.1 and the resolved supporting dependencies.
Linux interpreter: CPython 3.12.3. Windows build interpreter: CPython 3.12.10.

Reference staging:

```sh
uv run python studio_scripts/stage_schemes.py --source ../wmlst/db/pubmlst
```

162 schemes, 1,432 data files, 116,388,783 bytes. Staged manifest SHA-256:
`65293f60de5dc9d55e3a8e20177a05f7220337888e517ab1a0fca6c919e6150d`.
Immediate source commit: `6cad46ffd9dfddfaa55f7993cf391f80a70556d3`.
All 162 schemas loaded; schema warnings are retained in results. Reference data
were not modified. Both-strand exact matching, no nucleotide ambiguity wildcards.

Verification commands:

```sh
QT_QPA_PLATFORM=offscreen uv run pytest studio_tests -q --junitxml=artifacts/pytest-linux.xml
uv run ruff check src/wmlstudio studio_tests studio_scripts studio_packaging
QT_QPA_PLATFORM=offscreen uv run python studio_scripts/capture_desktop.py
uv run python studio_scripts/validate_local_cohort.py
uv run python studio_scripts/write_validation_manifest.py
```

`validate_local_cohort.py` reruns the four local assemblies and 10,000-record read
QC, checks all P. aeruginosa alleles against the preserved older WMLST/BLAST
output, and captures final source hashes. The older comparator used WMLST 1.2.1
and BLAST 2.12.0 against the same PubMLST database. This is regression agreement,
not independent clinical validation or an estimate of population accuracy.

S. epidermidis input SHA-256:
`a18f883d81d4a09350526706698a9e2a0637e3240f2acd638a60fa854c8ded39`.
Expected ST184; observed all seven alleles exact. P. aeruginosa accessions
GCF059737285v1, GCF059544455v1 and GCF049679795v1 gave ST155, ST2952 and ST1858.
FASTQ QC explicitly samples a prefix; the hash covers all compressed file bytes.

The practice exercise uses seed 42 and seven synthetic loci. No stochastic steps
are used in exact matching or graph construction; graph tie ordering is deterministic.
Generated screenshots and PDF were inspected visually; report layout was corrected
for Qt's supported rich-text subset. Local Linux X11 launch also passed after
extracting libxcb-cursor0 0.1.4-1build1 into a task-specific temporary directory.
No system packages were installed.

HYDRA importer contract pinned to upstream commit
`6d36c109491c16544e8919fe6962b4b62e97d3d7` (1.4.0). Import preserves evidence,
parameters, database information and source report SHA-256. Upstream performance
claims were not reproduced and are not presented as WMLSTudio results.

Windows packaging uses separate Windows wheels and interpreter under an isolated
Wine prefix. The exact final build outcome is recorded in SUMMARY.md and package
artifacts. Wine execution is distinct from clean-machine Windows 11 acceptance.

Additional real read-pair inspection:

```sh
uv run wmlstudio-cli check-pair \
  /media/iowa/u/tesseract_fastq/saureus/SRR12343864_1.fastq.gz \
  /media/iowa/u/tesseract_fastq/saureus/SRR12343864_2.fastq.gz --max-reads 10000
```

10,000 corresponding identifiers were checked. The result is `unmarked` because
the headers lack explicit mate indicators; matching file names alone are not used
to claim validated pair orientation. The remaining records were not inspected.

Final Windows artifact: `dist/WMLSTudio-Windows-x64.zip`, 46,294,204 bytes,
SHA-256 `9470e4e1e27666eec42dfea83d5c09b71ada25b38933dc83bd4648cca34422f3`.
Archive integrity and relocated execution from a directory containing spaces and
`é` passed under Wine: the extracted CLI reported ST184 and the expected input
hash. The final frozen GUI passed both offscreen and Windows Qt platform smoke
checks after a screenshot-only animation-settling adjustment. The 198-test source
suites preceded that adjustment. Detailed context and evidence are retained in
`artifacts/windows-build-context.json`, `artifacts/frozen-windows-wine-check.json`
and `artifacts/frozen-windows-relocated.json`. No clean Windows 11 host was available.
