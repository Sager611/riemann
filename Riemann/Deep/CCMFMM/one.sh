#!/bin/bash
# one certificate run (ccmfmm): prep (Python) + structured LDL (C, point mode), one process at a time, each under alarm + /usr/bin/time -l
# usage: one.sh <tag> <c> <sector even|odd|defl> <eps> <bits> <K2|-> <nu1> [nu2 ...]     (defl needs CCM_G=<gstate json>)
cd "$(dirname "$0")"
R=${CCM_RUNS:-../../../.scratch/ccmfmm/runs}; X=${CCM_BIN:-../../../.scratch/ccmfmm}
export WZ_RSCAN=$X/rscan WZ_SCAN=$X/scan OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
tag=$1; c=$2; sec=$3; eps=$4; bits=$5; k2=$6; shift 6
b=$R/$tag
[ -s $b.idx ] || /usr/bin/time -l perl -e 'alarm 1790; exec @ARGV' uv run --no-project --with mpmath --with numpy --with python-flint python3 prep_ccm.py $c $sec $eps $bits $b $k2 > $b.prep.json 2> $b.prep.time
BRD=""; [ "$sec" = defl ] && BRD=$b.brd
SLT_BRD=$BRD SLT_POINT=1 /usr/bin/time -l perl -e 'alarm 1790; exec @ARGV' $X/slt $b $bits -g $b.gen "$@" > $b.slt.json 2> $b.slt.time
cat $b.prep.json $b.slt.json
echo "RSS_MB prep $(rg -o '[0-9]+  maximum resident' $b.prep.time | awk '{printf "%.0f", $1/1048576}') slt $(rg -o '[0-9]+  maximum resident' $b.slt.time | awk '{printf "%.0f", $1/1048576}') ; time prep $(rg -o '[0-9.]+ real' $b.prep.time | head -1) slt $(rg -o '[0-9.]+ real' $b.slt.time | head -1)"
