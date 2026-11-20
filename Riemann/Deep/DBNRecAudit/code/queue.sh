#!/bin/bash
# dbnrecaudit queue.sh: the falsification runs, one process at a time, each under memcap 290 MB and /usr/bin/time -l, single thread.
cd /Users/adriansager/playground/riemann
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
L=Riemann/Deep/DBNRecAudit/logs; C=Riemann/Deep/DBNRecAudit/code
UV="uv run --no-project --with mpmath --with numpy --with python-flint python3"
log() { echo "$(date -u +%H:%M:%S) $*" >> $L/queue.log; }
step() { name=$1; to=$2; shift 2; log "start $name"
  /usr/bin/time -l bash .scratch/memcap.sh 290 $to "$@" > $L/$name.out 2> $L/$name.err; rc=$?
  log "end $name rc=$rc $(rg -o 'memcap: peak [0-9]+ MB.*' $L/$name.err | tail -1) $(rg -o '[0-9]+  maximum resident' $L/$name.err | tail -1)"; }
step dcoef 600 $UV $C/falsify_dcoef.py
step cmpdat_r8_00 600 $UV $C/cmp_dat.py /tmp/dbnra_keep/r8_00.dat Riemann/Deep/DBNRec2/logs/r8_00.prep.json 8
step dense_r1 1500 $UV $C/falsify_dense.py 1
step dense_r8 1500 $UV $C/falsify_dense.py 8
step scan_r8 1500 $UV $C/scan_wrap.py 8 300 8
for r in 2 3 4 5 6 7; do step dense_r$r 1200 $UV $C/falsify_dense.py $r; done
log "queue done"
