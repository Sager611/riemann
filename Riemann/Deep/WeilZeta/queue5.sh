#!/bin/sh
# queue5.sh: independent certificates at Chekmachev's own eps: c = 25 odd, c = 29 even and odd (his margins: 8.078775e-120, 1.426111e-146, 2.257925e-141)
cd /Users/adriansager/playground/riemann/Riemann/Deep/WeilZeta
export WZ_RSCAN=/Users/adriansager/playground/riemann/.scratch/weilzeta/rscan WZ_SCAN=/Users/adriansager/playground/riemann/.scratch/weilzeta/scan WZ_CHOL_BS=64
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
R="nice -n 10 uv run --no-project --with mpmath --with numpy --with python-flint python3"
perl -e 'alarm 1800; exec @ARGV' /usr/bin/time -l $R cert_zeta.py 25 odd 4.0546561810437569560e-119 576 >> logs/z25.log 2>> logs/z25.time
perl -e 'alarm 2400; exec @ARGV' /usr/bin/time -l $R cert_zeta.py 29 even 3.3887641037680338350e-145 640 >> logs/z29.log 2>> logs/z29.time
perl -e 'alarm 2400; exec @ARGV' /usr/bin/time -l $R cert_zeta.py 29 odd 1.5329466444386175069e-140 640 >> logs/z29.log 2>> logs/z29.time
echo QUEUE5 DONE >> logs/z29.log
