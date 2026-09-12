# SKESA native Windows port

This recipe targets SKESA 2.4.0, exact upstream commit
c1413581e4f37211892d3c4310d01f3d9a9b3490, with the upstream NO_NGS
configuration. It does not use WSL, Docker, or Linux emulation at runtime.
Source: https://github.com/ncbi/SKESA/tree/c1413581e4f37211892d3c4310d01f3d9a9b3490

Use the Native Windows SKESA GitHub Actions workflow or run
bash studio_packaging/build_skesa_windows.sh in an MSYS2 UCRT64 shell.
Required packages are declared in .github/workflows/skesa-windows.yml.
The output directory must not already exist. Compiler/dependency versions,
source commit, patch hash, binary hashes and exact modified source are packaged.
Build dependencies are resolved from MSYS2 and their exact installed versions
are recorded; this is not a bit-for-bit reproducible dependency lock.

## Reviewed portability changes

* Add BSD unsigned integer aliases on Windows, preserving exact bit widths.
* Include stdint.h explicitly in glb_align.cpp; its fixed-width integer types
  were previously obtained indirectly from Linux standard-library headers.
* Use GCC's __builtin_ffsll instead of the unavailable POSIX ffsll;
  both return the one-based index of the first set bit, or zero for zero.
* Compile the two upstream NO_NGS translation units with UCRT64 libraries,
  omitting Linux-only rt and dl linker dependencies. No assembly algorithm
  parameters or source logic are otherwise changed.
* Use C++14 for contemporary Boost compatibility and its native -mt DLL import
  libraries. Current Boost.System and Boost.Regex are header-only.

The build retains upstream's SSE4.2 CPU requirement and disables online NGS/SRA
access. Its directory includes the recursive non-system DLL closure; the build
fails on unresolved imports. A real paired-read assembly of a seeded 12 kbp
synthetic genome must recover a >=9 kbp contig and all returned contigs must
match ground truth before the artifact is published. This is a platform smoke
test, not biological or clinical validation. Windows execution remains
unverified until this workflow has passed on a Windows runner.

## Licensing and distribution

Read SKESA-LICENSE.txt: most upstream code carries NCBI public-domain notices,
but the GATB integer implementation carries AGPL version 3 terms.
Do not describe the entire executable as public-domain or omit these files.
The package includes exact modified source, the portability patch, build
recipe and available installed dependency-license texts. Distributors must
review and meet the applicable source/license obligations, including those of
the dynamically linked compiler/runtime libraries; this note is not a legal
certification. Keep the source archive beside the executable in redistribution.

WMLSTudio invokes this separately distributed engine as a subprocess. Original
FASTQ files are never modified. Basic contiguity and read-QC metrics do not
establish genome completeness, species purity, or absence of contamination.
QUAST and QuickClade are not bundled/executed by this adapter; these biological
QC stages must be reported as not performed, not as passed.
