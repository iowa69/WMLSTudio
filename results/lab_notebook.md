# WMLSTudio validation notebook

## 2026-09-12 — Native application and regression checks

Goal: replace the browser workflow with a portable desktop and retain transparent
scientific evidence. Used a new package; preserved existing legacy deletions.

The synthetic practice scheme uses seed 42 and seven 420-base loci. Known profiles,
an exact duplicate and a partial assembly test both positive and negative states.
`studio_scripts/capture_desktop.py` reruns this exercise and renders native pages,
an interactive forest and a PDF. Images are generated evidence, not UI mockups.

The real S. epidermidis positive control returned ST184. Three distinct P. aeruginosa
assemblies returned ST155, ST2952 and ST1858 and agreed with the older WMLST/BLAST
caller using the same allele database. A real S. aureus FASTQ file exercised a
10,000-record prefix, with full compressed-file hashing. Exact local paths and
results remain in ignored validation artifacts; source sequences are unmodified.

Corrections found through tests: project switching could retain old HYDRA text;
practice projects could duplicate samples; interrupted results needed job-state
filtering; file extensions were insufficient for choosing read versus assembly
processing. Regression tests now cover those cases. Imported schemes are checked
again after copying and when a cached snapshot is reused.

Native Linux launch initially failed because libxcb-cursor0 was missing. An Ubuntu
runtime package was extracted into an isolated temporary directory, and an actual
X11 smoke run then succeeded without changing system packages. Windows build and
Wine validation are recorded separately from Windows 11 acceptance.
