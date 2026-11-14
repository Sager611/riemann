#!/bin/bash
# scan_grid.sh (bothaaudit): float scan over logs/grid.txt (c N bits L_botha), one point at a time (one.sh), resumable.
cd "$(dirname "$0")"
while read c N bits L; do
  rg -q "\"c\": \"$c\"" logs/scan.jsonl 2>/dev/null && continue
  ./one.sh logs/scan.jsonl $c $N $bits 3 2 4 $L > /dev/null 2>> logs/scan.err
done < logs/grid.txt
echo DONE >> logs/scan.err
