#!/bin/bash
cd /Users/adriansager/playground/riemann/Riemann/Deep/GapCeil
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
mkdir -p logs/evid
for p in a b c; do perl -e 'alarm 900; exec @ARGV' /usr/bin/time -l uv run --no-project --with numpy --with scipy python3 evid.py $p > logs/evid/evid_$p.txt 2>&1; done
for k in 10 14; do perl -e 'alarm 900; exec @ARGV' /usr/bin/time -l uv run --no-project --with python-flint --with numpy --with scipy python3 evid_arb.py $k 16 > logs/evid/evid_arb_$k.txt 2>&1; done
echo EVID DONE > logs/evid/done.txt
