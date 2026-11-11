#!/bin/bash
# batch_L3.sh: Laurent model for S = {inf} u {p <= P}, P = 13 .. 97 (|S| = 7 .. 26), t = 1, tmax = 100, depth 40 + #primes;
# tmax 200 (float32) checks at P = 13, 23, 47.
cd "$(dirname "$0")"
PR=(2 3 5 7 11 13 17 19 23 29 31 37 41 43 47 53 59 61 67 71 73 79 83 89 97)
for i in $(seq 5 24); do
  P=$(IFS=,; echo "${PR[*]:0:$((i+1))}"); top=${PR[$i]}; d=$((40 + i + 1))
  ./runl.sh L_P${top}_t100 $P 100 $d 5,20
  case $top in 13|23|47) LA_F32=1 ./runl.sh L_P${top}_t200 $P 200 $d 5,20;; esac
done
echo BATCH_DONE $(date -u +%H:%M:%S)
