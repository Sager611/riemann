#!/bin/bash
cd /Users/adriansager/playground/riemann/.scratch/lostlaw
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
# wait for the running cert_c16 (chain3 child) to finish: one compute process at a time
while pgrep -f "certmulti.py 16 512" > /dev/null; do sleep 2; done
echo "cert_c16 finished (seen by chain4)" >> logs/chain.status
watch() { while kill -0 $1 2>/dev/null; do for pid in $(pgrep -f "$2"); do r=$(ps -o rss= -p $pid 2>/dev/null | tr -d ' '); if [ -n "$r" ] && [ "$r" -gt 290000 ]; then kill -9 $pid; echo "$3 KILLED rss=$r" >> logs/chain.status; fi; done; sleep 0.3; done; }
cert() { # c prec gallog tag
  tag=$4
  /usr/bin/time -l perl -e 'alarm 1750; exec @ARGV' uv run --no-project --with mpmath --with python-flint python3 certmulti.py $1 $2 sets_c$1.json $3 > logs/$tag.log 2> logs/$tag.time &
  tp=$!; watch $tp "certmulti.py $1 $2" $tag; wait $tp; echo "$tag exit=$?" >> logs/chain.status
}
gal() { # c N bits ymax setsfile rig
  tag=gal_c$1_N$2
  /usr/bin/time -l perl -e 'alarm 1500; exec @ARGV' uv run --no-project --with mpmath --with python-flint python3 gal.py $1 $2 $3 $4 $5 $6 > logs/$tag.log 2> logs/$tag.time &
  tp=$!; watch $tp "gal.py $1 $2 $3" $tag; wait $tp; echo "$tag exit=$?" >> logs/chain.status
}
trig() { # c N bits ks
  tag=theta_c$1_N$2
  /usr/bin/time -l perl -e 'alarm 900; exec @ARGV' uv run --no-project --with mpmath --with python-flint python3 theta_rig.py $1 $2 $3 $4 > logs/$tag.log 2> logs/$tag.time &
  tp=$!; watch $tp "theta_rig.py $1 $2" $tag; wait $tp; echo "$tag exit=$?" >> logs/chain.status
}
for c in 20.03 20.06 20.09 20.12 20.15 20.18; do gal $c 320 576 1.9 - norig; done
gal 8 200 320 2.05 sets_min.json norig
gal 12 200 384 2.05 sets_min.json norig
gal 14 240 416 2.05 sets_min.json norig
gal 18 300 480 2.05 sets_min.json norig
gal 22 320 576 2.05 sets_min.json norig
gal 24 340 608 2.05 sets_min.json norig
gal 26 380 640 2.05 sets_min.json norig
gal 28 420 672 2.05 sets_min.json norig
gal 20 320 576 1.9 - norig
trig 16 300 512 1,10,30,69
trig 30 300 704 1,27,74,129
gal 10 520 416 2.3 sets_c10.json rig
echo chain4 done >> logs/chain.status
