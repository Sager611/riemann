#!/bin/sh
# sequential, single-threaded rerun of all tests (about 12 min)
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
R="uv run --no-project --with mpmath --with numpy --with scipy python3"
for s in "g1.py" "g2.py C" "g3.py" "g4.py" "g5.py"; do
  o=$(echo $s | cut -d. -f1)_out.txt
  perl -e 'alarm 1500; exec @ARGV' $R $s > $o 2>&1
done
