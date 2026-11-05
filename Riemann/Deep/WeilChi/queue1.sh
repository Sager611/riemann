#!/bin/sh
# queue1.sh: sequential jobs (one compute process at a time), started after the running m4 c=47 eps=0 job
cd /Users/adriansager/playground/riemann/Riemann/Deep/WeilChi
export WEILCHI_SCAN=/Users/adriansager/playground/riemann/.scratch/weilchi/scan WEILCHI_CSCAN=/Users/adriansager/playground/riemann/.scratch/weilchi/cscan
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
while kill -0 98712 2>/dev/null; do sleep 5; done
R="uv run --no-project --with mpmath --with numpy --with python-flint python3"
perl -e 'alarm 1200; exec @ARGV' $R scan_check.py m4 47 201 >> logs/scan_check.log 2>&1
for s in even odd; do perl -e 'alarm 3000; exec @ARGV' $R cert_chi.py m4 47 $s 2^-200 448 201 >> logs/m4_47_cert.log 2>&1; done
for c in 20 25 30 40; do for s in even odd; do perl -e 'alarm 1800; exec @ARGV' $R cert_chi.py m3 $c $s 0 384 201 >> logs/m3_sweep.log 2>&1; done; done
echo QUEUE1 DONE >> logs/m3_sweep.log
