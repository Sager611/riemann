#!/bin/bash
# dbnrecaudit queue4.sh: after queue3.sh, the anchor scan of row 8 shell 0 again with 32 in-anchor offsets (the 8-offset run was too coarse:
# row 8's 2-group varies by < 0.08 over an anchor, so the 8 offsets were 2.27 apart). One process at a time, memcap 290 MB.
cd /Users/adriansager/playground/riemann
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
L=Riemann/Deep/DBNRecAudit/logs; C=Riemann/Deep/DBNRecAudit/code
UV="uv run --no-project --with mpmath --with numpy --with python-flint python3"
while pgrep -f 'DBNRecAudit/code/queue[23]?.sh' > /dev/null; do sleep 5; done
log() { echo "$(date -u +%H:%M:%S) $*" >> $L/queue.log; }
log "start scan_r8_32off"
/usr/bin/time -l bash .scratch/memcap.sh 290 1500 $UV $C/scan_wrap.py 8 300 32 > $L/scan_r8_32off.out 2> $L/scan_r8_32off.err
log "end scan_r8_32off rc=$? $(rg -o 'memcap: peak [0-9]+ MB.*' $L/scan_r8_32off.err | tail -1) $(rg -o '[0-9]+  maximum resident' $L/scan_r8_32off.err | tail -1)"
log "queue4 done"
