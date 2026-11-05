#!/bin/sh
# queue4.sh: chi_-3 at c = 47 (K2 = 1001), trivial windows, Galerkin cross-checks, even-character sweeps
cd /Users/adriansager/playground/riemann/Riemann/Deep/WeilChi
export WEILCHI_SCAN=/Users/adriansager/playground/riemann/.scratch/weilchi/scan WEILCHI_CSCAN=/Users/adriansager/playground/riemann/.scratch/weilchi/cscan
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
R="uv run --no-project --with mpmath --with numpy --with python-flint python3"
for s in even odd; do perl -e 'alarm 3500; exec @ARGV' $R cert_chi.py m3 47 $s 2^-275 448 1001 >> logs/m3_47_cert.log 2>&1; done
for s in even odd; do
  perl -e 'alarm 600; exec @ARGV' $R cert_chi.py m4 3 $s 0 128 >> logs/trivial.log 2>&1
  perl -e 'alarm 600; exec @ARGV' $R cert_chi.py m4 3 $s 2^-4 128 >> logs/trivial.log 2>&1
  perl -e 'alarm 600; exec @ARGV' $R cert_chi.py m3 2 $s 0 128 >> logs/trivial.log 2>&1
  perl -e 'alarm 600; exec @ARGV' $R cert_chi.py m3 2 $s 2^-4 128 >> logs/trivial.log 2>&1
  perl -e 'alarm 600; exec @ARGV' $R cert_chi.py p5 2 $s 0 128 >> logs/trivial.log 2>&1
  perl -e 'alarm 600; exec @ARGV' $R cert_chi.py p8 3 $s 0 128 >> logs/trivial.log 2>&1
done
perl -e 'alarm 1200; exec @ARGV' $R galerkin_chi.py m4 16 60,120,200 256 >> logs/galerkin.log 2>&1
perl -e 'alarm 1200; exec @ARGV' $R galerkin_chi.py m3 16 60,120,200 256 >> logs/galerkin.log 2>&1
for c in 16; do for ch in m4 m3; do perl -e 'alarm 600; exec @ARGV' $R cert_chi.py $ch $c even 0 256 >> logs/galerkin.log 2>&1; done; done
for c in 10 20 30; do for s in even odd; do
  perl -e 'alarm 1200; exec @ARGV' $R cert_chi.py p5 $c $s 0 384 601 >> logs/p5_sweep.log 2>&1
  perl -e 'alarm 1200; exec @ARGV' $R cert_chi.py p8 $c $s 0 384 601 >> logs/p8_sweep.log 2>&1
done; done
echo QUEUE4 DONE >> logs/p8_sweep.log
