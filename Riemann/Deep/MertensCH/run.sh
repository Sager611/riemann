#!/bin/sh
# run.sh -- full certificate (about 40 min, one single-threaded process at a time, < 120 MB)
set -e
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cd "$(dirname "$0")"; mkdir -p logs; rm -f logs/JL0.json
UV="uv run --no-project --with python-flint --with numpy python3"
$UV bprime.py > logs/bprime.log 2>&1                              # routine rerun with B' (Mertens78 pipeline)
$UV cert_ch.py > logs/cert_ch.log 2>&1                             # new bound, Mertens78 (7/8) certificates
$UV cert_ch.py --bp > logs/cert_ch_bp.log 2>&1                     # new bound, + B' certificates (headline)
MCH_TP=1000 $UV cert_ch.py --bp > logs/cert_ch_bp_T1000.log 2>&1   # routine baseline: same pipeline without the RH strip
python3 table.py > logs/table.md                                   # comparison table (floating point)
