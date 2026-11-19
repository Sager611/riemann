#!/bin/bash
# padhi2 Thm 5.2: re-certify binding cells one at a time, 420 s each.
cd "$(dirname "$0")/.."
for c in 32,68 30,72 31,71 31,69; do bash code/run.sh 420 repcert_$(echo $c | tr , _) code/repcert_cells.py $c; done
