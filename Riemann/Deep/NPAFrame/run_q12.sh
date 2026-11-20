#!/bin/bash
# npaframe queue 12: closed-band theta = 3/4 at N = 12 with the same cheap pricing as N12_To9_cheap (meeting test), after q10/q11
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export PRUNE=15000
cd "$(dirname "$0")"
while pgrep -f 'run_q10.sh|run_q11.sh' > /dev/null; do sleep 5; done
run() { prog=$1; shift; N=$1; S=$2; tag=N${N}_${S}${TAG}; WALL=${WALL:-850} /usr/bin/time -l bash ../../../.scratch/memcap.sh 290 1800 uv run --no-project --with numpy --with scipy python3 $prog "$@" logs/$tag.npz > logs/$tag.log 2> logs/$tag.err; echo "$tag :: $(tail -1 logs/$tag.log) :: $(rg -o 'memcap.*' logs/$tag.err) $(rg -o '[0-9]+  maximum resident' logs/$tag.err)" >> logs/summary.log; }
TAG=_cheap WALL=900 run cg_law.py 12 T9 150 96 1 100
