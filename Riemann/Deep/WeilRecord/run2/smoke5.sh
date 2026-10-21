#!/bin/sh
# smoke test of the patched bc.py at mu=5 (small N); log: smoke5.log
cd /Users/adriansager/playground/riemann/.scratch/weilrec2
PY=../weil/venv/bin/python
nice -n 10 $PY bc.py grams 5 140 320 > grams_5.log 2>&1 || { echo "grams failed"; exit 1; }
nice -n 10 $PY bc.py band 5 140 48 7/4 320 > band_5.log 2>&1 || { echo "band failed"; exit 1; }
BC_WFILE=w_5_8000.txt BC_SUMMARY=cert_5_summary.json nice -n 10 $PY bc.py cert 5 grams_5_140_320.bin band_5_140_48_7o4_320.bin 12192/10000 1/4 1/100000000000000000000 33 10 320 - 140 > cert_5.log 2>&1
echo "cert exit $?"
