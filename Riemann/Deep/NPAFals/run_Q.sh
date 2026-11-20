#!/bin/bash
# npafals queue: one process at a time. Jobs: "C N S [IT B ST]" = colgen_simple.py (continuous, RH world);
# "L N M depths S world [-rays]" = npa_cg_lean.py (lattice Z_M, off-line pairs at the listed depths, rays = infinitely deep pairs).
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
UV="uv run --no-project --with numpy --with scipy python3"
for job in "$@"; do
  set -- $job; kind=$1; shift
  if [ "$kind" = C ]; then
    N=$1; S=$2; IT=${3:-60}; B=${4:-256}; ST=${5:-200}; log=logs/A_N${N}_${S}.log
    /usr/bin/time -l bash ../../../.scratch/memcap.sh 290 1000 $UV colgen_simple.py $N $S $IT $B 1 $ST > $log 2>&1
    echo "C $job :: $(rg 'FINAL' $log) $(rg -o 'peak [0-9]+ MB.*' $log) $(rg 'maximum resident' $log | awk '{print $1}')" >> logs/summary.log
  else
    log=logs/L_N$1_M$2_$4_$5$6.log
    /usr/bin/time -l bash ../../../.scratch/memcap.sh 290 1500 $UV npa_cg_lean.py "$@" > $log 2>&1
    echo "L $job :: $(rg 'RESULT' $log) $(rg -o 'peak [0-9]+ MB.*' $log) $(rg 'maximum resident' $log | awk '{print $1}')" >> logs/summary.log
  fi
done
echo QUEUE-DONE >> logs/summary.log
