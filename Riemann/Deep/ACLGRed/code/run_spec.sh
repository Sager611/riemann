#!/bin/bash
# full low spectra (both sectors, all eigenvalues) of QW_lambda^N for eigenvalue counts; one process at a time
cd /Users/adriansager/playground/riemann/Riemann/Deep/ACLGFals/code
OUT=/Users/adriansager/playground/riemann/.scratch/aclgred/spec_all.jsonl
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
for spec in "5 256 240" "10 320 240" "16 416 240" "20 512 240" "25 608 300" "30 704 300"; do
  set -- $spec
  /usr/bin/time -l perl -e 'alarm 900; exec @ARGV' uv run --no-project --with python-flint==0.9.0 python3 spec_small.py $1 $3 $2 400 >> $OUT 2> /Users/adriansager/playground/riemann/.scratch/aclgred/spec_c$1.time
  echo "c=$1 N=$3 done $(date +%H:%M:%S) rss=$(rg -o '[0-9]+  maximum resident' /Users/adriansager/playground/riemann/.scratch/aclgred/spec_c$1.time)" >> /Users/adriansager/playground/riemann/.scratch/aclgred/spec_progress.log
done
echo ALLDONE >> /Users/adriansager/playground/riemann/.scratch/aclgred/spec_progress.log
