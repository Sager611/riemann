#!/bin/sh
# queue12.sh: fallback chi_-3 window c = 30 (both parities, eps = 2^-170, K2 = 201), then c = 40 odd with a single Woodbury solve
cd /Users/adriansager/playground/riemann/Riemann/Deep/WeilChi
export WEILCHI_SCAN=/Users/adriansager/playground/riemann/.scratch/weilchi/scan WEILCHI_CSCAN=/Users/adriansager/playground/riemann/.scratch/weilchi/cscan
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
while pgrep -f 'cert_chi.py m3 40 even' > /dev/null; do sleep 3; done
R="uv run --no-project --with mpmath --with numpy --with python-flint python3"
for s in even odd; do perl -e 'alarm 900; exec @ARGV' /usr/bin/time -l $R cert_chi.py m3 30 $s 2^-170 384 201 >> logs/m3_30_cert.log 2>> logs/m3_30.time; done
perl -e 'alarm 1500; exec @ARGV' /usr/bin/time -l $R cert_chi.py m3 40 odd 2^-226 384 601 >> logs/m3_40_cert.log 2>> logs/m3_40.time
echo QUEUE12 DONE >> logs/m3_40_cert.log
