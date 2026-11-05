#!/bin/sh
# queue1.sh: validation runs c = 16 (ChekAudit: even margin 7.2936e-77 at eps_16 = 1.512800554e-75; odd 2.993e-72 at 2.105612047e-71)
cd /Users/adriansager/playground/riemann/Riemann/Deep/WeilZeta
export WZ_RSCAN=/Users/adriansager/playground/riemann/.scratch/weilzeta/rscan WZ_SCAN=/Users/adriansager/playground/riemann/.scratch/weilzeta/scan
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
R="nice -n 10 uv run --no-project --with mpmath --with numpy --with python-flint python3"
perl -e 'alarm 900; exec @ARGV' /usr/bin/time -l $R cert_zeta.py 5 even 7.447489666e-18 256 >> logs/val.log 2>> logs/val.time
perl -e 'alarm 900; exec @ARGV' /usr/bin/time -l $R cert_zeta.py 5 odd 5.693186289e-15 256 >> logs/val.log 2>> logs/val.time
perl -e 'alarm 900; exec @ARGV' /usr/bin/time -l $R cert_zeta.py 10 even 9.133835606e-44 384 >> logs/val.log 2>> logs/val.time
perl -e 'alarm 900; exec @ARGV' /usr/bin/time -l $R cert_zeta.py 16 even 1.512800554e-75 512 >> logs/val.log 2>> logs/val.time
perl -e 'alarm 900; exec @ARGV' /usr/bin/time -l $R cert_zeta.py 16 odd 2.105612047e-71 512 >> logs/val.log 2>> logs/val.time
echo QUEUE1 DONE >> logs/val.log
