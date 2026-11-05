#!/bin/sh
# queue2.sh: first scaling runs c = 20, 25 (eps ~ Chekmachev's eps_c / 100 or below)
cd /Users/adriansager/playground/riemann/Riemann/Deep/WeilZeta
export WZ_RSCAN=/Users/adriansager/playground/riemann/.scratch/weilzeta/rscan WZ_SCAN=/Users/adriansager/playground/riemann/.scratch/weilzeta/scan WZ_CHOL_BS=64
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
R="nice -n 10 uv run --no-project --with mpmath --with numpy --with python-flint python3"
for s in even odd; do perl -e 'alarm 1800; exec @ARGV' /usr/bin/time -l $R cert_zeta.py 20 $s 2^-327 448 >> logs/z20.log 2>> logs/z20.time; done
for s in even odd; do perl -e 'alarm 2400; exec @ARGV' /usr/bin/time -l $R cert_zeta.py 25 $s 2^-416 576 >> logs/z25.log 2>> logs/z25.time; done
echo QUEUE2 DONE >> logs/z25.log
