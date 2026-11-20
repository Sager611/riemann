#!/bin/bash
# npaframe queue 13 (01:58 UTC; replaces the last job of run_q10.sh and run_q11/q12.sh): theta = 3/4 at N = 12, open and closed, same
# cheap pricing (B 96, 100 Adam steps), PRUNE 9000 keeping the newest 4500 columns (memory), then the law analyses. One process at a time.
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export PRUNE=9000
cd "$(dirname "$0")"
while ! rg -q 'FINAL|Traceback' logs/N5_T4.log 2>/dev/null; do sleep 2; done
pkill -f run_q10.sh; sleep 2; pkill -f 'cg_law.py 12 To9'; sleep 2
run() { prog=$1; shift; N=$1; S=$2; tag=N${N}_${S}${TAG}; WALL=${WALL:-850} /usr/bin/time -l bash ../../../.scratch/memcap.sh 290 1800 uv run --no-project --with numpy --with scipy python3 $prog "$@" logs/$tag.npz > logs/$tag.log 2> logs/$tag.err; echo "$tag :: $(tail -1 logs/$tag.log) :: $(rg -o 'memcap.*' logs/$tag.err) $(rg -o '[0-9]+  maximum resident' logs/$tag.err)" >> logs/summary.log; }
TAG=_cheap WALL=900 run cg_law.py 12 To9 150 96 1 100
uv run --no-project --with numpy python3 formfactor.py logs/N8_Po.npz logs/N8_To6.npz logs/N12_Po.npz logs/N12_P.npz logs/N12_To6.npz logs/N12_T6.npz > logs/formfactor.log 2>&1
TAG=_cheap WALL=900 run cg_law.py 12 T9 150 96 1 100
