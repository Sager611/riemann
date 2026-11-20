#!/bin/bash
# npaframe queue 2 (replaces the tail of run_q1.sh): fixed-theta series; one process at a time, memcap 290 MB
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
cd "$(dirname "$0")"
run() { prog=$1; shift; N=$1; S=$2; tag=N${N}_${S}${TAG}; WALL=${WALL:-850} /usr/bin/time -l bash ../../../.scratch/memcap.sh 290 1800 uv run --no-project --with numpy --with scipy python3 $prog "$@" logs/$tag.npz > logs/$tag.log 2> logs/$tag.err; echo "$tag :: $(tail -1 logs/$tag.log) :: $(rg -o 'memcap.*' logs/$tag.err) $(rg -o '[0-9]+  maximum resident' logs/$tag.err)" >> logs/summary.log; }
run cg_law.py 9 T6 150 192 1 150
run cg_law.py 8 T2 150 192 1 150
WALL=1500 run cg_law.py 12 P 150 192 1 150
WALL=1500 run cg_law.py 12 T3 150 192 1 150
WALL=1500 run cg_law.py 12 T6 150 192 1 150
TAG=_law run cg_law.py 8 P 150 192 1 150
WALL=1500 run cg_law.py 12 T9 150 192 1 150
