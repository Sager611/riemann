#!/bin/bash
cd /Users/adriansager/playground/riemann/.scratch/lostlaw
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
watch() { while kill -0 $1 2>/dev/null; do for pid in $(pgrep -f "$2"); do r=$(ps -o rss= -p $pid 2>/dev/null | tr -d ' '); if [ -n "$r" ] && [ "$r" -gt 290000 ]; then kill -9 $pid; echo "$3 KILLED rss=$r" >> logs/chain.status; fi; done; sleep 0.3; done; }
job() { tag=$1; al=$2; shift 2
  /usr/bin/time -l perl -e "alarm $al; exec @ARGV" uv run --no-project --with mpmath --with python-flint python3 "$@" > logs/$tag.log 2> logs/$tag.time &
  tp=$!; watch $tp "$1 $2 $3" $tag; wait $tp; echo "$tag exit=$?" >> logs/chain.status; }
job tailritz_c10 900 tailritz.py 10 520 1000 416
job gal_c5_N500 300 gal.py 5 500 320 0.01 - rig
job tailritz_c16 900 tailritz.py 16 520 900 512
echo chain6 done >> logs/chain.status
