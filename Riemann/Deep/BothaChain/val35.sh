#!/bin/bash
# val35.sh (bothachain): scan premise at c = 35 by rscan2 (compared with CCMFMM's cache), then val35.py <bits> (U(35) and the odd certificate at
# CCMFMM's level with the reduced precision). Each process under memcap 290 MB.   usage: val35.sh <bits>
cd "$(dirname "$0")"; export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
R=$(cd ../../..; pwd)/.scratch/bothachain/ccm; V=$(cd ../../..; pwd)/.scratch/bothachain/val35; mkdir -p $V/runs
export CCM_RUNS=$V/runs CCM_BIN=$R
U="uv run --no-project --with mpmath --with numpy --with python-flint python3"
[ -s $R/scan/scan_zeta_35_2043_rig.json ] || (cd ccm && WZ_RSCAN2=$R/rscan2 WZ_SCAN=$R/scan PYTHONPATH=../../WeilZeta bash ../../../../.scratch/memcap.sh 290 600 $U scan_zeta.py 35 2043 > $V/scan35.out 2>&1)
python3 -c "
import json
x=json.load(open('$R/scan/scan_zeta_35_2043_rig.json')); y=json.load(open('../CCMFMM/logs/scan/scan_zeta_35_2043_rig.json'))
print(json.dumps(dict(scan35_equal_ccmfmm=all(x[k]==y[k] for k in ('ks','kse','selo','sele')), selo=len(x['selo']), sele=len(x['sele']))))" > logs/val35.jsonl
bash ../../../.scratch/memcap.sh 290 1500 $U val35.py $1 >> logs/val35.jsonl 2> $V/val35.err
tail -1 $V/val35.err >> logs/val35.jsonl
