#!/bin/bash
# padhi2: queue after 'all' (one process at a time): repulsion off, N3 removed, KLN sensitivities.
cd "$(dirname "$0")/.."
while pgrep -f "qlab.py all 797" > /dev/null; do sleep 5; done
CELLS_OVERRIDE=logs/norep_lo.json bash code/run.sh 900 qv_norep_797 code/qlab.py published 797 1000
bash code/run.sh 900 qv_none_797 code/qlab.py none 797 1000
bash code/run.sh 600 qv_kln2_867 code/qlab.py kln2 867 887
bash code/run.sh 600 qv_nokln_867 code/qlab.py nokln 867 887
