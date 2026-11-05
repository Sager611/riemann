#!/bin/sh
# queue3.sh: (1) re-validate the column-blocked Woodbury code (c = 16 both sectors, c = 20 even); (2) early-exit scan rscan2 vs rscan (c = 20, 25)
cd /Users/adriansager/playground/riemann/Riemann/Deep/WeilZeta
export WZ_RSCAN=/Users/adriansager/playground/riemann/.scratch/weilzeta/rscan WZ_SCAN=/Users/adriansager/playground/riemann/.scratch/weilzeta/scan WZ_CHOL_BS=64
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
R="nice -n 10 uv run --no-project --with mpmath --with numpy --with python-flint python3"
perl -e 'alarm 900; exec @ARGV' /usr/bin/time -l $R cert_zeta.py 16 even 1.512800554e-75 512 >> logs/val2.log 2>> logs/val2.time
perl -e 'alarm 900; exec @ARGV' /usr/bin/time -l $R cert_zeta.py 16 odd 2.105612047e-71 512 >> logs/val2.log 2>> logs/val2.time
perl -e 'alarm 900; exec @ARGV' /usr/bin/time -l $R cert_zeta.py 20 even 2^-327 448 >> logs/val2.log 2>> logs/val2.time
mkdir -p /Users/adriansager/playground/riemann/.scratch/weilzeta/scan2
export WZ_RSCAN2=/Users/adriansager/playground/riemann/.scratch/weilzeta/rscan2 WZ_SCAN=/Users/adriansager/playground/riemann/.scratch/weilzeta/scan2
for c in 20 25; do perl -e 'alarm 900; exec @ARGV' /usr/bin/time -l $R scan_zeta.py $c - 4 >> logs/scancmp.log 2>> logs/scancmp.time; done
echo QUEUE3 DONE >> logs/scancmp.log
