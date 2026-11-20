#!/bin/bash
# npaframe queue 3: which words carry the fixed-theta gain? Dg = diagonal quartics |p_a p_b|^2 only; Ph = phase words (lam != mu) only
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
cd "$(dirname "$0")"
while pgrep -f run_q2.sh > /dev/null; do sleep 10; done
run() { prog=$1; shift; N=$1; S=$2; tag=N${N}_${S}${TAG}; WALL=${WALL:-850} /usr/bin/time -l bash ../../../.scratch/memcap.sh 290 1800 uv run --no-project --with numpy --with scipy python3 $prog "$@" logs/$tag.npz > logs/$tag.log 2> logs/$tag.err; echo "$tag :: $(tail -1 logs/$tag.log) :: $(rg -o 'memcap.*' logs/$tag.err) $(rg -o '[0-9]+  maximum resident' logs/$tag.err)" >> logs/summary.log; }
run cg_sets.py 8 Dg4 150 192 1 150
run cg_sets.py 8 Ph4 150 192 1 150
run cg_sets.py 10 Dg5 150 192 1 150
run cg_sets.py 10 Ph5 150 192 1 150
run cg_sets.py 8 Dg6 150 192 1 150
run cg_sets.py 8 Ph6 150 192 1 150
