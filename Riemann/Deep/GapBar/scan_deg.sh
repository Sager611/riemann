#!/bin/bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
for d in 4 8 12 16 20 24 32; do
  perl -e 'alarm 600; exec @ARGV' uv run --no-project --with numpy --with scipy python3 copt.py 2 0.2499 $d 3.30 scan
done
for d in 4 8 12 16; do
  perl -e 'alarm 600; exec @ARGV' uv run --no-project --with numpy --with scipy python3 copt.py 1 0.2499 $d 3.28 scan
done
