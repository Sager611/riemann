#!/bin/bash
cd /Users/adriansager/playground/riemann/.scratch/lostlaw
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
while pgrep -f "chain4.sh" > /dev/null; do sleep 3; done
watch() { while kill -0 $1 2>/dev/null; do for pid in $(pgrep -f "$2"); do r=$(ps -o rss= -p $pid 2>/dev/null | tr -d ' '); if [ -n "$r" ] && [ "$r" -gt 290000 ]; then kill -9 $pid; echo "$3 KILLED rss=$r" >> logs/chain.status; fi; done; sleep 0.3; done; }
job() { # tag alarm cmd...
  tag=$1; al=$2; shift 2
  /usr/bin/time -l perl -e "alarm $al; exec @ARGV" uv run --no-project --with mpmath --with python-flint python3 "$@" > logs/$tag.log 2> logs/$tag.time &
  tp=$!; watch $tp "$1 $2 $3" $tag; wait $tp; echo "$tag exit=$?" >> logs/chain.status
}
for c in 10.00 10.04 10.08 10.12 10.16 10.20 10.24; do job gal_c${c}_N200 600 gal.py $c 200 416 1.9 sets_scan10.json; done
for c in 30.00 30.02 30.04 30.06 30.08 30.10; do job gal_c${c}_N400 900 gal.py $c 400 704 1.75 sets_scan30.json; done
for c in 20.00 20.06 20.12 20.18 20.21 20.24 20.27 20.30; do job gal_c${c}_N320s 600 gal.py $c 320 576 1.9 sets_scan20.json; done
job zhp 900 zhp.py 250 1,2,3,10,20,27,30,50,74,100,129,150
KS=1,2,3,10,20,27,30,50,74,100,129,150
ZFILE=logs/zhp.log KLIST=$KS job gal_c35_N460 1500 gal.py 35 460 800 1.9 -
ZFILE=logs/zhp.log KLIST=$KS job gal_c40_N460 1500 gal.py 40 460 896 1.9 -
job cert_c10_new 600 certmulti.py 10 416 sets_c10.json logs/gal_c10_N400.log full,d1
K2MUL=3 job cert_c10_K3 900 certmulti.py 10 416 sets_c10.json logs/gal_c10_N400.log full,d1,d14,d12
job gal_c16_N520x 900 gal.py 16 520 512 0.01 sets_c16x.json rig
job cert_c16x 1200 certmulti.py 16 512 sets_c16x.json logs/gal_c16_N520x.log
K2MUL=3 job cert_c10_K3full 900 certmulti.py 10 416 sets_c10.json logs/gal_c10_N520.log
SOLVEMODE=resid job cert_c10_resid 600 certmulti.py 10 416 sets_c10.json logs/gal_c10_N400.log full,d1,d14
job zhp2 600 zhp.py 250 34,41,61,91,109,123,147
cat logs/zhp.log logs/zhp2.log > logs/zhp_all.log
job gal_c10_add 600 gal.py 10 400 416 0.01 sets_add_c10.json rig
job gal_c16_add 900 gal.py 16 520 512 0.01 sets_add_c16.json rig
job gal_c20_add 900 gal.py 20 480 576 0.01 sets_add_c20.json rig
job gal_c25_add 900 gal.py 25 360 640 0.01 sets_add_c25.json norig
job gal_c30_add 900 gal.py 30 450 704 0.01 sets_add_c30.json norig
ZFILE=logs/zhp_all.log job gal_c35_add 900 gal.py 35 460 800 0.01 sets_add_c35.json norig
ZFILE=logs/zhp_all.log job gal_c40_add 900 gal.py 40 460 896 0.01 sets_add_c40.json norig
echo waiting-go >> logs/chain.status
while [ ! -f logs/go ]; do sleep 2; done
K2MUL=3 job cert_c10_K3add 900 certmulti.py 10 416 sets_add_c10.json logs/gal_c10_add.log
cat logs/gal_c20_N480.log logs/gal_c20_add.log > logs/gal_c20_comb.log
SOLVEMODE=resid job cert_c20 1780 certmulti.py 20 512 sets_c20comb.json logs/gal_c20_comb.log full,d1,d14,d41,d94,d12,d123,d1_41,d1_41_94,A3,A5,pair27_1e-40,pair27_1e-20,d105
job tailritz_c10 900 tailritz.py 10 520 1000 416
job tailritz_c5 600 tailritz.py 5 200 700 320
echo chain5 done >> logs/chain.status
