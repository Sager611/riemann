#!/bin/bash
# sequential chain, one compute process at a time; alarm 1500 s; external RSS watchdog (kill > 290 MB)
cd /Users/adriansager/playground/riemann/.scratch/lostlaw
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
watch() { # tp pattern tag
  while kill -0 $1 2>/dev/null; do
    for pid in $(pgrep -f "$2"); do
      r=$(ps -o rss= -p $pid 2>/dev/null | tr -d ' ')
      if [ -n "$r" ] && [ "$r" -gt 290000 ]; then kill -9 $pid; echo "$3 KILLED rss=$r" >> logs/chain.status; fi
    done
    sleep 0.3
  done
}
run() { # c N bits rig ymax
  tag=gal_c$1_N$2
  /usr/bin/time -l perl -e 'alarm 1500; exec @ARGV' uv run --no-project --with mpmath --with python-flint python3 gal.py $1 $2 $3 $5 sets_c$1.json $4 > logs/$tag.log 2> logs/$tag.time &
  tp=$!
  watch $tp "gal.py $1 $2 $3" $tag
  wait $tp; echo "$tag exit=$?" >> logs/chain.status
}
run 16 200 512 norig 2.3
run 20 300 576 norig 2.3
run 20 480 576 rig 2.3
run 25 360 640 norig 2.3
run 30 450 704 norig 1.75
run 16 520 512 rig 2.3
echo chain2 done >> logs/chain.status
