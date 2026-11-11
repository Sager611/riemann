#!/bin/bash
# run_spec.sh: small-lambda Galerkin spectra (N = 120, 240) for the LA / gap / count tests; one process at a time.
cd "$(dirname "$0")"
OUT=../logs/spec.jsonl
for spec in "2 192" "3 192" "4 224" "5 256" "6 256" "8 288" "10 320" "12 352" "16 416" "20 512"; do
  set -- $spec
  for N in 120 240; do
    /usr/bin/time -l perl -e 'alarm 1200; exec @ARGV' uv run --no-project --with python-flint==0.9.0 python3 spec_small.py $1 $N $2 30 >> $OUT 2> ../logs/spec_c$1_N$N.time
    echo "c=$1 N=$N done $(date -u +%H:%M:%S) rss=$(rg -o '[0-9]+  maximum resident' ../logs/spec_c$1_N$N.time)" >> ../logs/spec_progress.log
  done
done
echo ALLDONE >> ../logs/spec_progress.log
