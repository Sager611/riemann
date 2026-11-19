#!/bin/bash
# padhi2: d* halved (cells lk <= 32, k <= 140), then q on the binding stretch.
cd "$(dirname "$0")/.."
DSTAR_SCALE=0.5 bash code/run.sh 1200 dstar_half code/dstar0.py 29,29.5,30,30.5,31,32 dhalf_lo 140
CELLS_OVERRIDE=logs/dhalf_lo.json bash code/run.sh 600 qv_dhalf_847 code/qlab.py published 847 947
