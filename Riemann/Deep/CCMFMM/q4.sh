#!/bin/bash
# queue (ccmfmm), one compute process at a time: waits for run_ccm.py 26..30, then chain [25, 30]; early-exit scan (../WeilZeta/scan_zeta.py
# + rscan2.c) and the window c = 35. (c = 40 dropped: time.)
cd "$(dirname "$0")"
S=../../../.scratch/ccmfmm; R=$S/runs
UVP="uv run --no-project --with mpmath --with numpy --with python-flint python3"
export WZ_RSCAN=$PWD/$S/rscan WZ_SCAN=$PWD/$S/scan OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
while pgrep -f 'run_ccm.py logs/int.jsonl 26' > /dev/null; do sleep 5; done
echo "chain25 start $(date -u +%H:%M:%S)"
$UVP chain_ccm.py logs/int.jsonl logs/chain.jsonl 25 30
for c in 35; do
  echo "scan $c start $(date -u +%H:%M:%S)"
  ( cd ../WeilZeta && PYTHONDONTWRITEBYTECODE=1 WZ_RSCAN2=$OLDPWD/$S/rscan2 WZ_SCAN=$OLDPWD/$S/scan /usr/bin/time -l perl -e 'alarm 1790; exec @ARGV' $UVP scan_zeta.py $c - 8 ) > logs/scan$c.out 2> logs/scan$c.time
  tail -1 logs/scan$c.out
done
echo "int35 start $(date -u +%H:%M:%S)"
$UVP run_ccm.py logs/int.jsonl 35
echo "q4 done $(date -u +%H:%M:%S)"
