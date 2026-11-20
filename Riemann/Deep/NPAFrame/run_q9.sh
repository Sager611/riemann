#!/bin/bash
# npaframe queue 9 (01:36 UTC, replaces run_q8.sh): OPEN-BAND model (pair words k <= N-1, AF's finite programme), fixed-theta gains
# To<n0> - Po, theta = 3/4 first; then theta = 1/2, 2/3; then npafals reproduction checks. One process at a time, memcap 290 MB.
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export PRUNE=15000
cd "$(dirname "$0")"
run() { prog=$1; shift; N=$1; S=$2; tag=N${N}_${S}${TAG}; WALL=${WALL:-850} /usr/bin/time -l bash ../../../.scratch/memcap.sh 290 1800 uv run --no-project --with numpy --with scipy python3 $prog "$@" logs/$tag.npz > logs/$tag.log 2> logs/$tag.err; echo "$tag :: $(tail -1 logs/$tag.log) :: $(rg -o 'memcap.*' logs/$tag.err) $(rg -o '[0-9]+  maximum resident' logs/$tag.err)" >> logs/summary.log; }
run cg_law.py 4 Po 150 384 1 300
run cg_law.py 4 To3 150 384 1 300
run cg_law.py 4 To2 150 384 1 300
WALL=1000 run cg_law.py 12 Po 150 192 1 150
WALL=1200 run cg_law.py 12 To9 150 192 1 150
run cg_law.py 6 Po 150 384 1 300
run cg_law.py 6 To3 150 384 1 300
run cg_law.py 6 To4 150 384 1 300
run cg_law.py 10 Po 150 192 1 150
run cg_law.py 10 To5 150 192 1 150
run cg_law.py 3 P 150 384 1 300
run cg_law.py 3 T2 150 384 1 300
run cg_law.py 5 P 150 384 1 300
run cg_law.py 5 T4 150 384 1 300
WALL=1000 run cg_law.py 12 To6 150 192 1 150
run cg_law.py 10 To7 150 192 1 150
