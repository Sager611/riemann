#!/bin/bash
# npaframe queue 6: reproduction checks of npafals' L4o - P at N = 3, 5 (T_{N-1} = L4o), after queue 5
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
cd "$(dirname "$0")"
while pgrep -f 'run_q2.sh|run_q5.sh' > /dev/null; do sleep 10; done
run() { prog=$1; shift; N=$1; S=$2; tag=N${N}_${S}${TAG}; WALL=${WALL:-850} /usr/bin/time -l bash ../../../.scratch/memcap.sh 290 1800 uv run --no-project --with numpy --with scipy python3 $prog "$@" logs/$tag.npz > logs/$tag.log 2> logs/$tag.err; echo "$tag :: $(tail -1 logs/$tag.log) :: $(rg -o 'memcap.*' logs/$tag.err) $(rg -o '[0-9]+  maximum resident' logs/$tag.err)" >> logs/summary.log; }
run cg_law.py 3 P 150 384 1 300
run cg_law.py 3 T2 150 384 1 300
run cg_law.py 5 P 150 384 1 300
run cg_law.py 5 T4 150 384 1 300
# law structure (types, word violations by level, between-type covariances)
uv run --no-project --with numpy python3 analyze_law.py logs/N12_P.npz logs/N12_T3.npz logs/N12_T6.npz logs/N8_P_law.npz > logs/analyze_laws.log 2>&1
# open band (AF's finite programme: pair words k <= N-1 only) and fixed-theta gains relative to it
run cg_law.py 8 Po 150 192 1 150
run cg_law.py 8 To4 150 192 1 150
run cg_law.py 8 To6 150 192 1 150
