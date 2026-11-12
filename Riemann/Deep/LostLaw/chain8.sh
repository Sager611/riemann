#!/bin/bash
cd /Users/adriansager/playground/riemann/.scratch/lostlaw
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
watch() { while kill -0 $1 2>/dev/null; do for pid in $(pgrep -f "$2"); do r=$(ps -o rss= -p $pid 2>/dev/null | tr -d ' '); if [ -n "$r" ] && [ "$r" -gt 290000 ]; then kill -9 $pid; echo "$3 KILLED rss=$r" >> logs/chain.status; fi; done; sleep 0.3; done; }
job() { tag=$1; al=$2; shift 2
  /usr/bin/time -l perl -e "alarm $al; exec @ARGV" uv run --no-project --with mpmath --with python-flint python3 "$@" > logs/$tag.log 2> logs/$tag.time &
  tp=$!; watch $tp "$1 $2 $3" $tag; wait $tp; echo "$tag exit=$?" >> logs/chain.status; }
job zhp3 900 zhp.py 310 1,2,3,10,20,30,50,74,100
KS=1,2,3,10,20,27,30,50,74,100,129,150
ZFILE=logs/zhp_all.log KLIST=$KS job gal_c40_N400 900 gal.py 40 400 896 1.9 -
ZFILE=logs/zhp_all.log KLIST=$KS job gal_c45_N420 900 gal.py 45 420 1024 1.9 -
ZFILE=logs/zhp3.log KLIST=1,2,3,10,20,30,50,74,100 job gal_c50_N400 900 gal.py 50 400 1152 1.9 -
echo chain8 done >> logs/chain.status
