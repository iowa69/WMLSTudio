# Decisions and prevention rules

- Preserve pre-existing deletions; build new code under a separate package.
- Never claim an absent locus is a novel allele, or FASTQ QC is assembly/typing.
- Never interpret unshared loci as zero differences or assign outbreak causality.
- All visible metrics must derive from actual project data, including the demo.
- A Linux test run does not validate a Windows executable.
- Measure what an addition costs the portable package, and gate the total. A size
  promise nobody measures stops being true one commit at a time.
- A tool with no official binary for a platform is absent on that platform and
  says so. Never substitute a different program, and never offer to run anyway.
- Bundle a core reference whole or not at all: shipping the genes without the
  point mutations makes an unperformed search read like a negative result.
- An automatic "install everything" must not accept a licence for the user.
