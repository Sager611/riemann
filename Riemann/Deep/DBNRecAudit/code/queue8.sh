#!/bin/bash
# dbnrecaudit queue8.sh: after queue7.sh (edges), the anchor scan of row 8 shell 0 with 16 offsets, 2-group held constant (the 32-offset run of
# queue6 was stopped at 17:33: ~30 min under load). memcap 290 MB, one process.
cd /Users/adriansager/playground/riemann
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
L=Riemann/Deep/DBNRecAudit/logs; C=Riemann/Deep/DBNRecAudit/code
while pgrep -f 'DBNRecAudit/code/queue7.sh' > /dev/null; do sleep 3; done
echo "$(date -u +%H:%M:%S) start scan_r8_16off" >> $L/queue.log
/usr/bin/time -l bash .scratch/memcap.sh 290 1500 uv run --no-project --with mpmath --with numpy --with python-flint python3 $C/scan_wrap.py 8 300 16 > $L/scan_r8_16off.out 2> $L/scan_r8_16off.err
echo "$(date -u +%H:%M:%S) end scan_r8_16off rc=$? $(rg -o 'memcap: peak [0-9]+ MB.*' $L/scan_r8_16off.err | tail -1) $(rg -o '[0-9]+  maximum resident' $L/scan_r8_16off.err | tail -1)" >> $L/queue.log
