#!/bin/bash
# batch_L.sh: Laurent (band-limited, unitary) model, t = 1, R = 5, 20, depth 30, tmax = 100, 200; one run at a time.
cd "$(dirname "$0")"
for P in 2,3 2,3,5 2,3,5,7 2,3,5,7,11; do
  tag=L_p$(echo $P | tr -d ,)
  for tm in 100 200; do
    [ -s ../logs/${tag}_t${tm}_d30.json ] && continue
    echo "== $tag $tm $(date -u +%H:%M:%S)"; ./runl.sh ${tag}_t${tm}_d30 $P $tm 30 5,20
  done
done
echo BATCH_DONE $(date -u +%H:%M:%S)
