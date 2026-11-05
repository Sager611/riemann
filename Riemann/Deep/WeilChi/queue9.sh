#!/bin/sh
# queue9.sh: test the robust Lemma 5.5 code (m4 c=30 even: smallest Ritz value must be 2.2313345e-36), then p8 c=40 odd, then chi_-3 c=47 odd
cd /Users/adriansager/playground/riemann/Riemann/Deep/WeilChi
export WEILCHI_SCAN=/Users/adriansager/playground/riemann/.scratch/weilchi/scan WEILCHI_CSCAN=/Users/adriansager/playground/riemann/.scratch/weilchi/cscan
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
true
R="uv run --no-project --with mpmath --with numpy --with python-flint python3"
perl -e 'alarm 600; exec @ARGV' $R cert_chi.py m4 30 even 0 448 201 >> logs/leancheck.log 2>&1
perl -e 'alarm 600; exec @ARGV' $R cert_chi.py m4 16 even 2.4e-18 256 601 >> logs/leancheck.log 2>&1
perl -e 'alarm 2400; exec @ARGV' /usr/bin/time -l $R cert_chi.py p8 40 odd 2^-91 384 601 >> logs/p8_40_cert.log 2>> logs/p8_40.time
perl -e 'alarm 3500; exec @ARGV' /usr/bin/time -l $R cert_chi.py m3 47 odd 2^-275 320 201 >> logs/m3_47_cert.log 2>> logs/m3_47_odd3.time
echo QUEUE9 DONE >> logs/m3_47_cert.log
