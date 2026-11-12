#!/bin/bash
cd /Users/adriansager/playground/riemann/.scratch/lostlaw
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
watch() { while kill -0 $1 2>/dev/null; do for pid in $(pgrep -f "$2"); do r=$(ps -o rss= -p $pid 2>/dev/null | tr -d ' '); if [ -n "$r" ] && [ "$r" -gt 290000 ]; then kill -9 $pid; echo "$3 KILLED rss=$r" >> logs/chain.status; fi; done; sleep 0.3; done; }
cert() { # c prec gallog
  tag=cert_c$1
  /usr/bin/time -l perl -e 'alarm 1700; exec @ARGV' uv run --no-project --with mpmath --with python-flint python3 certmulti.py $1 $2 sets_c$1.json $3 > logs/$tag.log 2> logs/$tag.time &
  tp=$!; watch $tp "certmulti.py $1 $2" $tag; wait $tp; echo "$tag exit=$?" >> logs/chain.status
}
gal() { # c N bits ymax setsfile
  tag=gal_c$1_N$2
  /usr/bin/time -l perl -e 'alarm 1500; exec @ARGV' uv run --no-project --with mpmath --with python-flint python3 gal.py $1 $2 $3 $4 $5 norig > logs/$tag.log 2> logs/$tag.time &
  tp=$!; watch $tp "gal.py $1 $2 $3" $tag; wait $tp; echo "$tag exit=$?" >> logs/chain.status
}
cert 10 416 logs/gal_c10_N400.log
cert 16 512 logs/gal_c16_N520.log
gal 8 200 320 2.05 sets_min.json
gal 12 200 384 2.05 sets_min.json
gal 14 240 416 2.05 sets_min.json
gal 18 300 480 2.05 sets_min.json
gal 22 320 576 2.05 sets_min.json
gal 24 340 608 2.05 sets_min.json
gal 26 380 640 2.05 sets_min.json
gal 28 420 672 2.05 sets_min.json
cert 20 576 logs/gal_c20_N480.log
echo chain3 done >> logs/chain.status
