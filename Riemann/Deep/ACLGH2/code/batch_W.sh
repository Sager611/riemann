#!/bin/bash
# batch_W.sh: shifted-band model at single coupling pockets of S = {inf,2,3,5,7,11} (pockets.py), window +-75, depth 16.
cd "$(dirname "$0")"
P=2,3,5,7,11
for tc in 3099.82 2447.58 16751.62 19851.48 2350.48 2877.02 3471.68 4352.44 2001.88 1500000; do
  ./runw.sh W_p235711_c$tc $P $tc 100 16
done
echo BATCH_DONE $(date -u +%H:%M:%S)
