# WMLSTudio 0.1.0 — development preview

Native Qt application implemented, including exact assembly typing, read QC,
saved projects, interactive comparisons, reports and HYDRA result integration.
No web UI, web server, WSL, or end-user Python installation is used by the frozen app.

## Observed validation

- Linux: **198 tests passed**, 18.14 seconds in the final full run.
- Windows Python 3.12.10 / Qt 6.11.2 under Wine 11.17: **198 tests passed**, 85.54 seconds.
- All **162** staged PubMLST schemas load. Reference inconsistencies are disclosed.
- S. epidermidis assembly: expected and observed **ST184**, all seven loci exact.
- Three P. aeruginosa assemblies: **ST155**, **ST2952**, **ST1858**; all allele calls
  agree with WMLST 1.2.1 / BLAST 2.12.0 using the same reference snapshot.
- Real compressed FASTQ: **10,000** records QC-checked, explicitly reported as a
  prefix; the SHA-256 covers the entire compressed input.
- Native desktop, compact layout, comparison graph and PDF inspected visually.
- Source distribution built at `dist/source/wmlstudio-0.1.0.tar.gz`.

Exact commands and caveats: [methods](tasks/METHODS.md).
Generated evidence: `results/2026-09-12_validation/`, Linux JUnit
`artifacts/pytest-linux.xml`, Windows/Wine JUnit `artifacts/windows-wine-pytest.xml`.
The source/environment manifest is `results/2026-09-12_validation/source-manifest.json`.

## Remaining product work

The Windows portable build is at
[`dist/WMLSTudio-Windows-x64.zip`](dist/WMLSTudio-Windows-x64.zip), approximately
46.3 MB compressed / 210 MB extracted. Frozen GUI startup and practice analysis
passed on both offscreen and Windows Qt platforms under Wine; the frozen CLI
typed the bundled positive control and the real ST184 assembly correctly.
The final ZIP passed integrity checks and real ST184 typing after extraction into
a directory containing spaces and `é`. SHA-256:
`9470e4e1e27666eec42dfea83d5c09b71ada25b38933dc83bd4648cca34422f3`.
Build context: `artifacts/windows-build-context.json`; native Windows Qt preview:
`artifacts/frozen-windows-desktop.png`.
**Clean Windows 11 acceptance has not been performed.** Wine success is a separate test.

HYDRA integration currently imports and displays its JSON evidence with provenance;
it does not execute HYDRA. Read assembly, approximate/novel allele inference,
full cgMLST benchmark validation and broader MLSTudio/SeqSphere+ parity remain
future work. This preview is not a production or diagnostic validation claim.
See [the capability matrix](docs/STUDIO_CAPABILITIES.md) and
[HYDRA integration](docs/HYDRA_INTEGRATION.md).
