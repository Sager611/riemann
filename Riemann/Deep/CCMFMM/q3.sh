#!/bin/bash
# queue (ccmfmm), one compute process at a time: waits for run_ccm.py 21..25, then chain [20, 25]; adversarial test (c = 10, 21);
# integer windows 26..30; chain [25, 30]; early-exit scans (../WeilZeta/scan_zeta.py + rscan2.c) and windows c = 35, 40.
cd "$(dirname "$0")"
S=../../../.scratch/ccmfmm; R=$S/runs
UVP="uv run --no-project --with mpmath --with numpy --with python-flint python3"
export WZ_RSCAN=$PWD/$S/rscan WZ_SCAN=$PWD/$S/scan OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
while pgrep -f 'run_ccm.py logs/int.jsonl 21' > /dev/null; do sleep 5; done
echo "chain20 start $(date -u +%H:%M:%S)"
CCM_U0=1.65743315856e-96 CCM_O0=3.76226271626722e-92 $UVP chain_ccm.py logs/int.jsonl logs/chain.jsonl 20 25
echo "adv start $(date -u +%H:%M:%S)"
mkdir -p logs/adv
for spec in "10 7.348097e-40 1.349542e-41 320 201" "21 1.195093e-97 1.618853e-99 640 735"; do
  set -- $spec; c=$1; e2=$2; eo=$3; bits=$4; k2=$5
  CCM_G=$R/g$c.json /usr/bin/time -l perl -e 'alarm 1790; exec @ARGV' $UVP prep_ccm.py $c defl $e2 $bits $R/adv${c}d $k2 1 > logs/adv/adv${c}d.prep.json 2> logs/adv/adv${c}d.prep.time
  SLT_BRD=$R/adv${c}d.brd SLT_POINT=1 /usr/bin/time -l perl -e 'alarm 1790; exec @ARGV' $S/slt $R/adv${c}d $bits -1 -0.1 0 $e2 > logs/adv/adv${c}d.slt.json 2> logs/adv/adv${c}d.slt.time
  /usr/bin/time -l perl -e 'alarm 1790; exec @ARGV' $UVP prep_ccm.py $c odd $eo $bits $R/adv${c}o $k2 1 > logs/adv/adv${c}o.prep.json 2> logs/adv/adv${c}o.prep.time
  SLT_POINT=1 /usr/bin/time -l perl -e 'alarm 1790; exec @ARGV' $S/slt $R/adv${c}o $bits -1 -0.1 0 $eo > logs/adv/adv${c}o.slt.json 2> logs/adv/adv${c}o.slt.time
done
echo "int26 start $(date -u +%H:%M:%S)"
$UVP run_ccm.py logs/int.jsonl 26 27 28 29 30
echo "chain25 start $(date -u +%H:%M:%S)"
$UVP chain_ccm.py logs/int.jsonl logs/chain.jsonl 25 30
for c in 35 40; do
  echo "scan $c start $(date -u +%H:%M:%S)"
  ( cd ../WeilZeta && PYTHONDONTWRITEBYTECODE=1 WZ_RSCAN2=$OLDPWD/$S/rscan2 WZ_SCAN=$OLDPWD/$S/scan /usr/bin/time -l perl -e 'alarm 1790; exec @ARGV' $UVP scan_zeta.py $c - 8 ) > logs/scan$c.out 2> logs/scan$c.time
  tail -1 logs/scan$c.out
done
echo "int35 start $(date -u +%H:%M:%S)"
$UVP run_ccm.py logs/int.jsonl 35 40
echo "q3 done $(date -u +%H:%M:%S)"
