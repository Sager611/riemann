#!/bin/bash
# npaframe queue 1: fixed-theta inside gains G(N, n0) = T<n0> - P (one process at a time, memcap 290 MB)
# args: N SET iters B seed steps (npz path appended)
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
cd "$(dirname "$0")"
run() { N=$1; S=$2; tag=N${N}_${S}; WALL=${WALL:-850} /usr/bin/time -l bash ../../../.scratch/memcap.sh 290 1800 uv run --no-project --with numpy --with scipy python3 cg_theta.py "$@" logs/$tag.npz > logs/$tag.log 2> logs/$tag.err; echo "$tag :: $(tail -1 logs/$tag.log) :: $(rg -o 'memcap.*' logs/$tag.err) $(rg -o '[0-9]+  maximum resident' logs/$tag.err)" >> logs/summary.log; }
run 4 T2 150 384 1 300
run 6 T3 150 384 1 300
run 8 T4 150 192 1 150
run 10 P 150 192 1 150
run 10 T5 150 192 1 150
WALL=1500 run 12 P 150 192 1 150
WALL=1500 run 12 T6 150 192 1 150
WALL=1500 run 12 T9 150 192 1 150
