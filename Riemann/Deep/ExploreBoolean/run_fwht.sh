#!/bin/bash
# run_fwht.sh — exact spectra sweep (sequential, single-threaded)
B=/Users/adriansager/playground/riemann/.scratch/explore33/fwht
echo "## system 0 (primes), all x <= 96"
for x in $(seq 2 96); do $B 0 $x 1 0 | head -2; done
for x in 23 47 89; do
 echo "## x=$x adversarial"
 $B 1 $x 1 0 | head -3
 for s in 1 2 3 4 5; do $B 2 $x $s 0.5 | head -2; done
 for s in 1 2 3 4 5; do $B 3 $x $s 0 | head -2; done
done
