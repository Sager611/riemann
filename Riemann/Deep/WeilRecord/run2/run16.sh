#!/bin/sh
# End-to-end re-run of the mu=16 (half-width log 4) certificate with the built-in checks (2026-10-20).
# One process at a time, niced; peak memory via /usr/bin/time -l.  Summary log: run16.log
cd /Users/adriansager/playground/riemann/.scratch/weilrec2
PY=../weil/venv/bin/python
T="/usr/bin/time -l"
echo "start $(date)"
nice -n 10 $T $PY bc.py schur 16 w_16_32000.txt > schur_16.log 2>&1 || { echo "schur failed"; exit 1; }
echo "schur done $(date)"
nice -n 10 $T $PY bc.py grams 16 1400 384 > grams_16.log 2>&1 || { echo "grams failed"; exit 1; }
echo "grams done $(date)"
nice -n 10 $T $PY bc.py band 16 1400 640 9/2 384 > band_16.log 2>&1 || { echo "band failed"; exit 1; }
echo "band done $(date)"
BC_WFILE=w_16_32000.txt BC_SUMMARY=cert_16_summary.json nice -n 10 $T $PY bc.py cert 16 grams_16_1400_384.bin band_16_1400_640_9o2_384.bin 36421/10000 1/2 1/10000000000000000000000000000000000000000000000000000000000000000000000000000 343 114 384 - 1400 > cert_16.log 2>&1
echo "cert exit $? $(date)"
