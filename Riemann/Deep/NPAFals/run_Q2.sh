#!/bin/bash
# npafals queue 2 (replaces run_Q.sh's remaining list): waits for the running job, then runs jobs one at a time.
# Job kinds as run_Q.sh plus "B N J s0list words..." = b_cluster_search.py.
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
UV="uv run --no-project --with numpy --with scipy python3"
while pgrep -f 'python3 (colgen_simple|npa_cg_lean|b_cluster)' > /dev/null; do sleep 3; done
for job in "$@"; do
  set -- $job; kind=$1; shift
  if [ "$kind" = C ]; then
    N=$1; S=$2; IT=${3:-60}; B=${4:-256}; ST=${5:-200}; log=logs/A_N${N}_${S}.log
    /usr/bin/time -l bash ../../../.scratch/memcap.sh 290 1000 $UV colgen_simple.py $N $S $IT $B 1 $ST > $log 2>&1
    echo "C $job :: $(rg 'FINAL' $log) $(rg -o 'peak [0-9]+ MB.*' $log) $(rg 'maximum resident' $log | awk '{print $1}')" >> logs/summary.log
  elif [ "$kind" = B ]; then
    /usr/bin/time -l bash ../../../.scratch/memcap.sh 290 900 $UV b_cluster_search.py "$@" >> logs/B_cluster.log 2>&1
    echo "B $job :: done $(tail -2 logs/B_cluster.log | rg -o 'peak [0-9]+ MB.*')" >> logs/summary.log
  else
    log=logs/L_N$1_M$2_$4_$5$6.log
    /usr/bin/time -l bash ../../../.scratch/memcap.sh 290 1500 $UV npa_cg_lean.py "$@" > $log 2>&1
    echo "L $job :: $(rg 'RESULT' $log) $(rg -o 'peak [0-9]+ MB.*' $log) $(rg 'maximum resident' $log | awk '{print $1}')" >> logs/summary.log
  fi
done
echo QUEUE2-DONE >> logs/summary.log
