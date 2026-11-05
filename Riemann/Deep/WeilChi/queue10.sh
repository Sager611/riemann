#!/bin/sh
# queue10.sh: chi_-3 at c = 40, both parities, eps = 2^-226 (fallback window); then chi_-3 at c = 47, odd f, K2 = 601 (320 bits)
cd /Users/adriansager/playground/riemann/Riemann/Deep/WeilChi
export WEILCHI_SCAN=/Users/adriansager/playground/riemann/.scratch/weilchi/scan WEILCHI_CSCAN=/Users/adriansager/playground/riemann/.scratch/weilchi/cscan
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
R="uv run --no-project --with mpmath --with numpy --with python-flint python3"
for s in even odd; do perl -e 'alarm 1500; exec @ARGV' /usr/bin/time -l $R cert_chi.py m3 40 $s 2^-226 384 601 >> logs/m3_40_cert.log 2>> logs/m3_40.time; done
perl -e 'alarm 2400; exec @ARGV' /usr/bin/time -l $R cert_chi.py m3 47 odd 2^-275 320 601 >> logs/m3_47_cert.log 2>> logs/m3_47_odd4.time
echo QUEUE10 DONE >> logs/m3_47_cert.log
