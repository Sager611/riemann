#!/bin/sh
# queue14.sh: certificates re-run with the RIGOROUS scan (rscan.c, RIGOR.md): chi_-3 c = 30, 35; chi_-4 c = 47
cd /Users/adriansager/playground/riemann/Riemann/Deep/WeilChi
export WEILCHI_SCAN=/Users/adriansager/playground/riemann/.scratch/weilchi/scan WEILCHI_RSCAN=/Users/adriansager/playground/riemann/.scratch/weilchi/rscan
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
R="nice -n 10 uv run --no-project --with mpmath --with numpy --with python-flint python3"
perl -e 'alarm 900; exec @ARGV' /usr/bin/time -l $R cert_chi.py m3 30 odd 2^-170 384 201 >> logs/rig_m3_30.log 2>> logs/rig.time
for s in even odd; do perl -e 'alarm 1500; exec @ARGV' /usr/bin/time -l $R cert_chi.py m3 35 $s 2^-200 384 201 >> logs/rig_m3_35.log 2>> logs/rig.time; done
for s in even odd; do perl -e 'alarm 1800; exec @ARGV' /usr/bin/time -l $R cert_chi.py m4 47 $s 2^-200 448 201 >> logs/rig_m4_47.log 2>> logs/rig.time; done
echo QUEUE14 DONE >> logs/rig_m4_47.log
