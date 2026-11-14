#!/bin/bash
# scan2.sh (bothaaudit): point.py over logs/grid.txt (c N bits L_botha), one point at a time, each under memcap 290 MB / 1700 s
# and /usr/bin/time -l; resumable (skips c already in logs/points.jsonl).
cd "$(dirname "$0")"; export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
while read c N bits L; do
  rg -q "\"c\": \"$c\"" logs/points.jsonl 2>/dev/null && continue
  /usr/bin/time -l bash ../../../.scratch/memcap.sh 290 1700 uv run --no-project --with python-flint --with mpmath python3 point.py logs/points.jsonl $c $N $bits $L > /dev/null 2> logs/pt_$c.time
  echo "$c $(rg -o 'memcap: peak [0-9]+ MB[^r]*' logs/pt_$c.time) maxrss_MB $(rg -o '[0-9]+  maximum resident' logs/pt_$c.time | awk '{printf "%.0f", $1/1048576}')" >> logs/points_mem.txt
done < logs/grid.txt
echo DONE >> logs/points_mem.txt
