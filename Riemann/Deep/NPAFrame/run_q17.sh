#!/bin/bash
# npaframe queue 17 (02:16 UTC): theta = 1/2 at N = 14 (open To7 and closed T7; P/Po first), meeting test; cg_lite, PRUNE 9000
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export PRUNE=9000 NCUE=1000
cd "$(dirname "$0")"
run() { prog=$1; shift; N=$1; S=$2; tag=N${N}_${S}${TAG}; WALL=${WALL:-850} /usr/bin/time -l bash ../../../.scratch/memcap.sh 290 1800 uv run --no-project --with numpy --with scipy python3 $prog "$@" logs/$tag.npz > logs/$tag.log 2> logs/$tag.err; echo "$tag :: $(tail -1 logs/$tag.log) :: $(rg -o 'memcap.*' logs/$tag.err) $(rg -o '[0-9]+  maximum resident' logs/$tag.err)" >> logs/summary.log; }
WALL=300 run cg_lite.py 14 Po 150 192 1 150
WALL=300 run cg_lite.py 14 P 150 192 1 150
WALL=700 run cg_lite.py 14 To7 150 192 1 150
WALL=700 run cg_lite.py 14 T7 150 192 1 150
