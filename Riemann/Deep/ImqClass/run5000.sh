#!/bin/bash
# run5000.sh: the N = 5000 sieve in four sequential chunks (one process at a time, single-threaded).
cd "$(dirname "$0")"
run() { perl -e 'alarm 3000; exec @ARGV' ./isieve 5000 $1 $2 runs/n5000_$3 || echo "FAIL chunk $3 exit $?" >> runs/n5000.status; echo "chunk $3 end $(date +%H:%M:%S)" >> runs/n5000.status; }
echo "start $(date +%H:%M:%S)" > runs/n5000.status
run 15032385536 49999999999 a
run 50000000000 99999999999 b
run 100000000000 152499999999 c
run 152500000000 190000000000 d
echo "ALL DONE $(date +%H:%M:%S)" >> runs/n5000.status
