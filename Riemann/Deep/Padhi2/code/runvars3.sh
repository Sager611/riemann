#!/bin/bash
# padhi2: after runvars2: per-row N3 report, then a sample (every 6th interval) of [1100, 3159) with N3 deleted.
cd "$(dirname "$0")/.."
while pgrep -f "runvars2.sh" > /dev/null; do sleep 5; done
bash code/run.sh 300 n3bind code/n3bind.py
bash code/run.sh 1200 qv_none_hi code/qlab.py none 1100 3159 6
