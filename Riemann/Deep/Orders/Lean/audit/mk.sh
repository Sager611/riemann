#!/bin/bash
# usage: bash .scratch/orderslean/mk.sh <Module> [timeout_s]   (run from repo root)
cd /Users/adriansager/playground/riemann
m=$1; to=${2:-1800}
LEANB=~/.elan/toolchains/leanprover--lean4---v4.35.0-rc2/bin/lean
OL=.scratch/orderslean/olean; d=Riemann/Deep/Orders/Lean
mkdir -p $OL/$d
export OMP_NUM_THREADS=1
t0=$(date +%s)
LEAN_PATH=$OL bash .scratch/memcap.sh 290 $to /usr/bin/time -l $LEANB -j1 -o $OL/$d/$m.olean $d/$m.lean > .scratch/orderslean/log_$m.txt 2> .scratch/orderslean/err_$m.txt
rc=$?
rss=$(rg -o '[0-9]+ +peak memory footprint' .scratch/orderslean/err_$m.txt | awk '{printf "%d", $1/1048576}')
pk=$(rg -o 'memcap: peak [0-9]+ MB.*' .scratch/orderslean/err_$m.txt)
echo "$(date -u '+%H:%M') $m rc=$rc $(( $(date +%s) - t0 ))s footprint=${rss}MB $pk" | tee -a .scratch/orderslean/peaks.log
head -c 6000 .scratch/orderslean/log_$m.txt

