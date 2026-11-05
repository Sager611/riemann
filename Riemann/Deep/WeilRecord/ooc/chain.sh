#!/bin/sh
# chain runner: one oocert process at a time, each under a hard 2400 s alarm; exit 3 = resume.
cd /Users/adriansager/playground/riemann/.scratch/weilooc
export VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 OOC_BIN=$PWD
RD=$1
while true; do
  /usr/bin/time -l perl -e 'alarm 2400; exec @ARGV' ./venv/bin/python ../../Riemann/Deep/WeilRecord/ooc/oocert.py $RD >> $RD/chain.out 2>&1
  rc=$?
  echo "$(date '+%H:%M:%S') oocert exit $rc" >> $RD/chain.out
  [ $rc -eq 3 ] || break
done
echo "$(date '+%H:%M:%S') chain done rc=$rc" >> $RD/chain.out
