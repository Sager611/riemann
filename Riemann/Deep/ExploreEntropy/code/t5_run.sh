#!/bin/bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
R="uv run --no-project --with mpmath --with python-flint python3 t5_movedpair.py"
$R 1 0.003,0.03,0.3 1.5 3.2
$R 2 0.03,0.3 1.5 4.0
$R 3 0.03,0.3 1.5 4.5
$R 5 0.003,0.03 2.6 4.6
$R 10 0.03,0.3 3.5 5.4
$R 20 0.03,0.3 5.0 7.6
