#!/bin/sh
# queue3.sh: why does m3 c=40 even fail at eps=0? knobs: K2, neighbour windows, t1
cd /Users/adriansager/playground/riemann/Riemann/Deep/WeilChi
export WEILCHI_SCAN=/Users/adriansager/playground/riemann/.scratch/weilchi/scan WEILCHI_CSCAN=/Users/adriansager/playground/riemann/.scratch/weilchi/cscan
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
while pgrep -f queue1.sh > /dev/null; do sleep 3; done
R="uv run --no-project --with mpmath --with numpy --with python-flint python3"
perl -e 'alarm 1200; exec @ARGV' $R cert_chi.py m3 40 even 0 384 601 >> logs/m3_40_tests.log 2>&1
WEILCHI_NBEK=3 WEILCHI_NBEJ=8 perl -e 'alarm 1200; exec @ARGV' $R cert_chi.py m3 40 even 0 384 201 >> logs/m3_40_tests.log 2>&1
WEILCHI_T1=2 perl -e 'alarm 1200; exec @ARGV' $R cert_chi.py m3 40 even 0 384 201 >> logs/m3_40_tests.log 2>&1
WEILCHI_NBO=64 perl -e 'alarm 1200; exec @ARGV' $R cert_chi.py m3 40 even 0 384 201 >> logs/m3_40_tests.log 2>&1
echo QUEUE3 DONE >> logs/m3_40_tests.log
