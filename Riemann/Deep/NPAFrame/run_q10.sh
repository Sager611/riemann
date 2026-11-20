#!/bin/bash
# npaframe queue 10 (01:39 UTC, replaces run_q9.sh after N4/N12 Po): OPEN-BAND fixed-theta series, cheapest first; one process at a time.
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export PRUNE=15000
cd "$(dirname "$0")"
run() { prog=$1; shift; N=$1; S=$2; tag=N${N}_${S}${TAG}; WALL=${WALL:-850} /usr/bin/time -l bash ../../../.scratch/memcap.sh 290 1800 uv run --no-project --with numpy --with scipy python3 $prog "$@" logs/$tag.npz > logs/$tag.log 2> logs/$tag.err; echo "$tag :: $(tail -1 logs/$tag.log) :: $(rg -o 'memcap.*' logs/$tag.err) $(rg -o '[0-9]+  maximum resident' logs/$tag.err)" >> logs/summary.log; }
run cg_law.py 6 Po 150 384 1 300
run cg_law.py 6 To3 150 384 1 300
run cg_law.py 6 To4 150 384 1 300
run cg_law.py 10 Po 150 192 1 150
run cg_law.py 10 To5 150 192 1 150
run cg_law.py 3 Po 150 384 1 300
run cg_law.py 3 To2 150 384 1 300
run cg_law.py 9 Po 150 192 1 150
run cg_law.py 9 To6 150 192 1 150
WALL=900 run cg_law.py 12 To6 150 192 1 150
run cg_law.py 3 P 150 384 1 300
run cg_law.py 3 T2 150 384 1 300
run cg_law.py 5 P 150 384 1 300
run cg_law.py 5 T4 150 384 1 300
TAG=_cheap WALL=1100 run cg_law.py 12 To9 150 96 1 100
