#!/bin/sh
# queue7.sh: m3 c=47 odd (t1 = 1/2, K2 = 201, NBO = 8, 320 bits), Galerkin bracket tests, even characters
cd /Users/adriansager/playground/riemann/Riemann/Deep/WeilChi
export WEILCHI_SCAN=/Users/adriansager/playground/riemann/.scratch/weilchi/scan WEILCHI_CSCAN=/Users/adriansager/playground/riemann/.scratch/weilchi/cscan
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
R="uv run --no-project --with mpmath --with numpy --with python-flint python3"
WEILCHI_T1=0.5 WEILCHI_NBO=8 perl -e 'alarm 3500; exec @ARGV' /usr/bin/time -l $R cert_chi.py m3 47 odd 2^-275 320 201 >> logs/m3_47_cert.log 2>> logs/m3_47_odd.time
for e in 4.8e-18 2.4e-18; do perl -e 'alarm 600; exec @ARGV' $R cert_chi.py m4 16 even $e 256 601 >> logs/bracket.log 2>&1; done
for e in 5.4e-25 2.6e-25; do perl -e 'alarm 600; exec @ARGV' $R cert_chi.py m3 16 even $e 256 601 >> logs/bracket.log 2>&1; done
for s in even odd; do
  perl -e 'alarm 600; exec @ARGV' $R cert_chi.py p5 2 $s 2^-4 128 >> logs/trivial.log 2>&1
  perl -e 'alarm 600; exec @ARGV' $R cert_chi.py p8 3 $s 2^-4 128 >> logs/trivial.log 2>&1
done
for s in even odd; do
  perl -e 'alarm 2400; exec @ARGV' /usr/bin/time -l $R cert_chi.py p5 40 $s 2^-140 384 601 >> logs/p5_40_cert.log 2>> logs/p5_40.time
  perl -e 'alarm 2400; exec @ARGV' /usr/bin/time -l $R cert_chi.py p8 40 $s 2^-91 384 601 >> logs/p8_40_cert.log 2>> logs/p8_40.time
done
echo QUEUE7 DONE >> logs/p8_40_cert.log
