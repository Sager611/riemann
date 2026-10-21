#!/bin/sh
# e3 (half-width 3/2) pipeline, one heavy process at a time, niced. Log: e3run.log
cd /Users/adriansager/playground/riemann/.scratch/weilrec3
PY=../weil/venv/bin/python; T="/usr/bin/time -l"
for k in 0 1 2; do ln -sf ../weilrec/gramspart_e3_3000_384_${k}of3.bin .; done
while [ "$(ls ../weilrec/gramspart_e3_3000_384_*of3.bin 2>/dev/null | wc -l | tr -d ' ')" -lt 3 ] || pgrep -f 'bc.py grams e3' >/dev/null; do sleep 60; done
echo "parts done $(date)"
nice -n 10 $T $PY bc.py merge e3 3000 384 3 > merge_e3.log 2>&1 || { echo "merge failed"; exit 1; }
echo "merge done $(date)"
nice -n 10 $T $PY bc.py band e3 3000 1216 103/20 384 > band_e3.log 2>&1 || { echo "band failed"; exit 1; }
echo "band done $(date)"
TAU="1/1$(printf '0%.0s' $(seq 1 100))"
nice -n 10 $T $PY cert3.py e3 grams_e3_3000_384.bin band_e3_3000_1216_103o20_384.bin w_e3_32000.txt 42959/10000 1/2 $TAU 660 219 384 - cert_e3_summary.json > cert_e3.log 2>&1
echo "cert exit $? $(date): $(rg RESULT cert_e3.log)"
