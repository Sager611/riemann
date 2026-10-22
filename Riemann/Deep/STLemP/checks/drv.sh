#!/bin/bash
# runs famqE (compile: gcc -O2 -o famqE famqE.c -lm) for X = 20000 and the 25 (N, mode) pairs of NOTE §8; 3 jobs in parallel, 400 s alarm each
cd "$(dirname "$0")"
run() { perl -e 'alarm 400; exec @ARGV' ./famqE 20000 $1 $2 $3 > out_$2_$1.txt 2>&1; }
export -f run
for N in 2500 10000 20000 40000 80000; do
  for ms in "3 0.5" "4 0.5" "5 0.8333333" "6 0.8333333" "2 0.8333333"; do echo "$N $ms"; done
done | xargs -P 3 -L 1 bash -c 'run $0 $1 $2'
