#!/bin/bash
# queue (ccmfmm): the window c = 35 (early-exit scan via ../WeilZeta/scan_zeta.py + rscan2.c, then run_ccm.py); one process at a time
cd "$(dirname "$0")"
S=../../../.scratch/ccmfmm
UVP="uv run --no-project --with mpmath --with numpy --with python-flint python3"
export WZ_RSCAN=$PWD/$S/rscan WZ_SCAN=$PWD/$S/scan OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
echo "scan 35 start $(date -u +%H:%M:%S)"
( cd ../WeilZeta && PYTHONDONTWRITEBYTECODE=1 WZ_RSCAN2=$OLDPWD/$S/rscan2 WZ_SCAN=$OLDPWD/$S/scan /usr/bin/time -l perl -e 'alarm 1790; exec @ARGV' $UVP scan_zeta.py 35 - 8 ) > logs/scan35.out 2> logs/scan35.time
tail -1 logs/scan35.out
echo "int35 start $(date -u +%H:%M:%S)"
$UVP run_ccm.py logs/int.jsonl 35
echo "q5 done $(date -u +%H:%M:%S)"
