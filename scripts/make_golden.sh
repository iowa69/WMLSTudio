#!/usr/bin/env bash
# Regenerate the golden reference corpus using the ORIGINAL tseemann/mlst (Perl).
# These goldens are the byte-for-byte compatibility contract that WMLST must satisfy.
# Requires: perl + Moo/List::MoreUtils/JSON/Path::Tiny, blastn, any2fasta, and a
# checkout of https://github.com/tseemann/mlst in $REF_MLST.
set -euo pipefail
REF_MLST="${REF_MLST:?set REF_MLST to a tseemann/mlst checkout}"
export MLST_DBDIR="$(cd "$(dirname "$0")/.." && pwd)/db"
cd "$(dirname "$0")/../tests/data"
G=../golden; mkdir -p $G
REF="$REF_MLST/bin/mlst --quiet --skipcheck"
run(){ n="$1"; shift; $REF "$@" > "$G/$n.out" 2> "$G/$n.err" || true; }
run default example.fna
run full --full example.fna
run csv --csv example.fna
run nopath --nopath --full example.fna
run label --label MYSAMPLE --full example.fna
run gz --full example.fna.gz
run gbk --full example.gbk.gz
run messy --full messy.fa
run mixedzip --full mixed.fa.zip
run nullfa --full null.fa
run novelfa --full novel.fa
run nonefa --full none.fa
run emptyfa --full empty.fa
run issue146 --full issue146.fa
run equality --full equality.fa.gz
run novelbz2 --full novel.fasta.bz2
run scheme_forced --scheme saureus --full example.fna
run legacy --scheme sepidermidis --legacy example.fna
run multi --full example.fna novel.fa null.fa none.fa messy.fa
$REF --json $G/example.json example.fna >/dev/null 2>&1 || true
$REF --novel $G/messy_novel.fa --full messy.fa >/dev/null 2>&1 || true
$REF --novel $G/lepto_novel.fa --full novel.fasta.bz2 >/dev/null 2>&1 || true
echo "golden corpus regenerated in $G"
