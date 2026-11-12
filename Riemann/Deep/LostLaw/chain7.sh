#!/bin/bash
cd /Users/adriansager/playground/riemann/.scratch/lostlaw
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
watch() { while kill -0 $1 2>/dev/null; do for pid in $(pgrep -f "$2"); do r=$(ps -o rss= -p $pid 2>/dev/null | tr -d ' '); if [ -n "$r" ] && [ "$r" -gt 290000 ]; then kill -9 $pid; echo "$3 KILLED rss=$r" >> logs/chain.status; fi; done; sleep 0.3; done; }
job() { tag=$1; al=$2; shift 2
  /usr/bin/time -l perl -e "alarm $al; exec @ARGV" uv run --no-project --with mpmath --with python-flint python3 "$@" > logs/$tag.log 2> logs/$tag.time &
  tp=$!; watch $tp "$1 $2 $3" $tag; wait $tp; echo "$tag exit=$?" >> logs/chain.status; }
job gal_c12_N520 600 gal.py 12 520 416 2.0 sets_c12.json rig
job tailritz_c12 600 tailritz.py 12 520 1000 416
K2MUL=2.2 job cert_c12_K22 1200 certmulti.py 12 416 sets_c12.json logs/gal_c12_N520.log
echo chain7 done >> logs/chain.status
