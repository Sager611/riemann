#!/bin/bash
# npaframe queue 5 (replaces run_q3.sh and run_q4.sh): after queue 2, one process at a time, memcap 290 MB.
# (a) pricing robustness: rerun N = 10 P and T5 with another seed and a 2x larger Adam batch / 2x steps;
# (b) which words carry the fixed-theta gain (Dg = diagonal quartics only, Ph = phase words only);
# (c) off-line immunity of fixed-theta sets on lattices Z_M with depths and rays (npa_cg_lean.py, exhaustive pricing).
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
cd "$(dirname "$0")"
while pgrep -f 'run_q2.sh' > /dev/null; do sleep 10; done
run() { prog=$1; shift; N=$1; S=$2; tag=N${N}_${S}${TAG}; WALL=${WALL:-850} /usr/bin/time -l bash ../../../.scratch/memcap.sh 290 1800 uv run --no-project --with numpy --with scipy python3 $prog "$@" logs/$tag.npz > logs/$tag.log 2> logs/$tag.err; echo "$tag :: $(tail -1 logs/$tag.log) :: $(rg -o 'memcap.*' logs/$tag.err) $(rg -o '[0-9]+  maximum resident' logs/$tag.err)" >> logs/summary.log; }
D=0.05,0.2,0.6,1.5,4.0
lrun() { tag=L_N$1_M$2_$4_$5${6:+$6}; /usr/bin/time -l bash ../../../.scratch/memcap.sh 290 1500 uv run --no-project --with numpy --with scipy python3 npa_cg_lean.py "$@" > logs/$tag.log 2>&1; echo "$tag :: $(rg 'RESULT' logs/$tag.log | cut -c1-170) :: $(rg -o 'memcap.*' logs/$tag.log)" >> logs/summary.log; }
TAG=_s2 run cg_law.py 10 P 150 384 2 300
TAG=_s2 run cg_law.py 10 T5 150 384 2 300
run cg_sets.py 8 Dg4 150 192 1 150
run cg_sets.py 8 Ph4 150 192 1 150
run cg_sets.py 10 Dg5 150 192 1 150
run cg_sets.py 10 Ph5 150 192 1 150
lrun 6 18 $D P on
lrun 6 18 $D T3 on
lrun 6 18 $D T3 off -rays
lrun 6 18 $D T4 on
lrun 6 18 $D T4 off -rays
lrun 7 21 $D P on
lrun 7 21 $D T5 on
lrun 7 21 $D T5 off -rays
WALL=1500 run cg_law.py 14 P 150 192 1 150
WALL=1500 run cg_law.py 14 T7 150 192 1 150
