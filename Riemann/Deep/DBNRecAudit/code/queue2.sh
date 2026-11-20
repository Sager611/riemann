#!/bin/bash
# dbnrecaudit queue2.sh: after queue.sh (waits for it), rerun the dense tests of rows 1 and 8 (their first runs died on a KeyError fixed at 16:54).
cd /Users/adriansager/playground/riemann
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
L=Riemann/Deep/DBNRecAudit/logs; C=Riemann/Deep/DBNRecAudit/code
UV="uv run --no-project --with mpmath --with numpy --with python-flint python3"
while pgrep -f 'DBNRecAudit/code/queue.sh' > /dev/null; do sleep 5; done
log() { echo "$(date -u +%H:%M:%S) $*" >> $L/queue.log; }
step() { name=$1; to=$2; shift 2; log "start $name"
  /usr/bin/time -l bash .scratch/memcap.sh 290 $to "$@" > $L/$name.out 2> $L/$name.err; rc=$?
  log "end $name rc=$rc $(rg -o 'memcap: peak [0-9]+ MB.*' $L/$name.err | tail -1) $(rg -o '[0-9]+  maximum resident' $L/$name.err | tail -1)"; }
step dense_r1 1500 $UV $C/falsify_dense.py 1
step dense_r8 1500 $UV $C/falsify_dense.py 8
log "queue2 done"
