#!/bin/bash
# padhi2: N3 variants and sensitivities on the binding stretch, one process at a time.
cd "$(dirname "$0")/.."
lo=${1:-797}; hi=${2:-1000}
for v in both all none kln2 nokln; do bash code/run.sh 900 qv_${v}_${lo} code/qlab.py $v $lo $hi; done
CELLS_OVERRIDE=logs/norep_lo.json bash code/run.sh 900 qv_norep_${lo} code/qlab.py published $lo $hi
