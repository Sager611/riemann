#!/bin/bash
# npaframe master queue (replaces the tails of run_q2/q5/q6.sh, killed 01:03 UTC): one process at a time, memcap 290 MB.
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export PRUNE=40000
cd "$(dirname "$0")"
run() { prog=$1; shift; N=$1; S=$2; tag=N${N}_${S}${TAG}; WALL=${WALL:-850} /usr/bin/time -l bash ../../../.scratch/memcap.sh 290 1800 uv run --no-project --with numpy --with scipy python3 $prog "$@" logs/$tag.npz > logs/$tag.log 2> logs/$tag.err; echo "$tag :: $(tail -1 logs/$tag.log) :: $(rg -o 'memcap.*' logs/$tag.err) $(rg -o '[0-9]+  maximum resident' logs/$tag.err)" >> logs/summary.log; }
D=0.05,0.2,0.6,1.5,4.0
lrun() { tag=L_N$1_M$2_$4_$5${6:+$6}; /usr/bin/time -l bash ../../../.scratch/memcap.sh 290 1500 uv run --no-project --with numpy --with scipy python3 npa_cg_lean.py "$@" > logs/$tag.log 2>&1; echo "$tag :: $(rg 'RESULT' logs/$tag.log | cut -c1-170) :: $(rg -o 'memcap.*' logs/$tag.log)" >> logs/summary.log; }
WALL=1200 run cg_law.py 12 T6 150 192 1 150
TAG=_s2 run cg_law.py 10 P 150 384 2 300
PRUNE=15000 TAG=_s2 run cg_law.py 10 T5 150 384 2 300
lrun 6 18 $D P on
lrun 6 18 $D T3 on
lrun 6 18 $D T3 off -rays
lrun 6 18 $D T4 on
lrun 6 18 $D T4 off -rays
lrun 7 21 $D P on
lrun 7 21 $D T5 on
lrun 7 21 $D T5 off -rays
PRUNE=15000 run cg_sets.py 8 Dg4 150 192 1 150
lrun 6 18 $D T3 on
lrun 6 18 $D P on
PRUNE=15000 run cg_sets.py 8 Ph4 150 192 1 150
PRUNE=15000 run cg_sets.py 10 Dg5 150 192 1 150
PRUNE=15000 run cg_sets.py 10 Ph5 150 192 1 150
PRUNE=15000 run cg_law.py 8 Po 150 192 1 150
PRUNE=15000 run cg_law.py 8 To4 150 192 1 150
PRUNE=15000 run cg_law.py 8 To6 150 192 1 150
PRUNE=15000 run cg_law.py 3 P 150 384 1 300
PRUNE=15000 run cg_law.py 3 T2 150 384 1 300
PRUNE=15000 run cg_law.py 5 P 150 384 1 300
PRUNE=15000 run cg_law.py 5 T4 150 384 1 300
PRUNE=15000 TAG=_law run cg_law.py 8 P 150 192 1 150
uv run --no-project --with numpy python3 analyze_law.py logs/N12_P.npz logs/N12_T3.npz logs/N12_T6.npz logs/N8_P_law.npz > logs/analyze_laws.log 2>&1
PRUNE=15000 WALL=1200 run cg_law.py 12 T9 150 192 1 150
PRUNE=15000 WALL=1200 run cg_law.py 14 P 150 192 1 150
PRUNE=15000 WALL=1200 run cg_law.py 14 T7 150 192 1 150
PRUNE=15000 WALL=1200 run cg_law.py 12 T8 150 192 1 150
