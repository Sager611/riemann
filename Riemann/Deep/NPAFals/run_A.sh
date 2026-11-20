#!/bin/bash
# npafals (A): continuous-position colgen, RH world, min E[#simple]/N; P vs L4o (inside, immune) vs L4, N = 4..; one process at a time.
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
for job in "$@"; do
  set -- $job; N=$1; S=$2; IT=${3:-50}; B=${4:-256}; ST=${5:-200}
  /usr/bin/time -l bash ../../../.scratch/memcap.sh 290 900 uv run --no-project --with numpy --with scipy python3 colgen_simple.py $N $S $IT $B 1 $ST > logs/A_N${N}_${S}.log 2>&1
  rg 'FINAL|memcap|maximum resident' logs/A_N${N}_${S}.log | tr '\n' ' ' >> logs/A_summary.log; echo >> logs/A_summary.log
done
