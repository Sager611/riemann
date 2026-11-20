#!/bin/bash
# npaframe queue 16 (02:14 UTC): cheap theta = 1/4 open-band check (To2 at N = 8, To3 at N = 12) against the closed 0.000020 / 0.000108
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export PRUNE=9000 NCUE=1000
cd "$(dirname "$0")"
run() { prog=$1; shift; N=$1; S=$2; tag=N${N}_${S}${TAG}; WALL=${WALL:-850} /usr/bin/time -l bash ../../../.scratch/memcap.sh 290 1800 uv run --no-project --with numpy --with scipy python3 $prog "$@" logs/$tag.npz > logs/$tag.log 2> logs/$tag.err; echo "$tag :: $(tail -1 logs/$tag.log) :: $(rg -o 'memcap.*' logs/$tag.err) $(rg -o '[0-9]+  maximum resident' logs/$tag.err)" >> logs/summary.log; }
WALL=500 run cg_law.py 8 To2 150 192 1 150
WALL=600 run cg_law.py 12 To3 150 192 1 150
