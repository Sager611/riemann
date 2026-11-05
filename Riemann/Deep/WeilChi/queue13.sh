#!/bin/sh
# queue13.sh: chi_-3 at c = 35, both parities, eps = 2^-200 (K2 = 201, then 401 for any failing parity)
cd /Users/adriansager/playground/riemann/Riemann/Deep/WeilChi
export WEILCHI_SCAN=/Users/adriansager/playground/riemann/.scratch/weilchi/scan WEILCHI_CSCAN=/Users/adriansager/playground/riemann/.scratch/weilchi/cscan
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
R="uv run --no-project --with mpmath --with numpy --with python-flint python3"
for s in even odd; do
  perl -e 'alarm 900; exec @ARGV' /usr/bin/time -l $R cert_chi.py m3 35 $s 2^-200 384 201 >> logs/m3_35_cert.log 2>> logs/m3_35.time
  if ! tail -1 logs/m3_35_cert.log | rg -q '"rigorous": true'; then
    perl -e 'alarm 900; exec @ARGV' /usr/bin/time -l $R cert_chi.py m3 35 $s 2^-200 384 401 >> logs/m3_35_cert.log 2>> logs/m3_35.time
  fi
done
echo QUEUE13 DONE >> logs/m3_35_cert.log
