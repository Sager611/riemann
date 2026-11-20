#!/bin/bash
# dbnrecaudit queue7.sh: after queue6.sh, the edge test (falsify_edges.py rows 1, 8). memcap 290 MB, one process.
cd /Users/adriansager/playground/riemann
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
L=Riemann/Deep/DBNRecAudit/logs; C=Riemann/Deep/DBNRecAudit/code
while pgrep -f 'DBNRecAudit/code/queue6.sh' > /dev/null; do sleep 3; done
echo "$(date -u +%H:%M:%S) start edges" >> $L/queue.log
/usr/bin/time -l bash .scratch/memcap.sh 290 900 uv run --no-project --with mpmath --with numpy --with python-flint python3 $C/falsify_edges.py 1 8 > $L/edges.out 2> $L/edges.err
echo "$(date -u +%H:%M:%S) end edges rc=$? $(rg -o 'memcap: peak [0-9]+ MB.*' $L/edges.err | tail -1) $(rg -o '[0-9]+  maximum resident' $L/edges.err | tail -1)" >> $L/queue.log
