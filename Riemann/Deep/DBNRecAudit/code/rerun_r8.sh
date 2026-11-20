#!/bin/bash
# dbnrecaudit rerun_r8.sh: from-scratch rerun of DBNRec2 row 8 (partial, shells u 0 -> 0.09) in /tmp/dbnra, then certify; compare with the logged cert.
# Copies: prep_end.py, bb_rec2.c, certify.py (paths sed'ed to /tmp/dbnra); inputs = audit's own Zenodo fetch (sha256-pinned). One process at a time, memcap 290 MB.
set -u
R=/Users/adriansager/playground/riemann; V=/tmp/dbnra
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
mkdir -p $V/code $V/logs $V/dl $V/ws/lambda_round15/flow_main_v2
cp $R/Riemann/Deep/DBNRec2/code/prep_end.py $R/Riemann/Deep/DBNRec2/code/bb_rec2.c $R/Riemann/Deep/DBNRec2/code/certify.py $V/code/
cp $R/.scratch/memcap.sh $V/memcap.sh
cp $R/.scratch/dbnrecaudit/dl/*ACTUAL_67_ROW_ARC_INPUTS.json $R/.scratch/dbnrecaudit/dl/*ALL_67_ORIGINAL_AFFINE_NATIVE.json $V/dl/
cp $R/.scratch/dbnrecaudit/dl/*MAIN_CLOCK_SOURCE.json $V/ws/lambda_round15/flow_main_v2/MAIN_CLOCK_SOURCE.json
sed -i '' -e "s#L = 'Riemann/Deep/DBNRec2/logs/'; S_ = '.scratch/dbnrec2/'#L = '/tmp/dbnra/logs/'; S_ = '/tmp/dbnra/'#" $V/code/certify.py
left=$(rg -n 'Riemann/|\.scratch|/Users/' $V/code); [ -z "$left" ] || { echo "repo paths left: $left"; exit 1; }
cd $V
cc --version | head -1 > logs/cc.txt; cc -O2 -o bb_rec code/bb_rec2.c -lm; shasum -a 256 bb_rec code/* >> logs/cc.txt
ARC=$V/dl/Proof_and_Replay__Lambda_Peer_Review__verification__arcs__ACTUAL_67_ROW_ARC_INPUTS.json
AFF=$V/dl/Proof_and_Replay__Lambda_Peer_Review__verification__referee__ALL_67_ORIGINAL_AFFINE_NATIVE.json
row=8; us=(0 0.01 0.02 0.035 0.05 0.07 0.09)
for ((j=0; j<6; j++)); do
  a=${us[$j]}; b=${us[$((j+1))]}; tag=r8_$(printf %02d $j)
  if [ $j -eq 0 ]; then T=0.913118; else T=0.913695; fi
  /usr/bin/time -l bash $V/memcap.sh 290 600 uv run --no-project --with mpmath --with python-flint python3 code/prep_end.py $ARC $AFF $row $a $b $T $V/$tag.dat > logs/$tag.prep.json 2> logs/$tag.prep.err
  /usr/bin/time -l bash $V/memcap.sh 290 1800 ./bb_rec $V/$tag.dat 48 400 0 0 0.001 logs/$tag.unres.txt > logs/$tag.bb.out 2> logs/$tag.bb.err
  echo "$(date -u +%H:%M:%S) $tag done $(rg -o 'CERT.*' logs/$tag.bb.out)" >> logs/progress.txt
done
/usr/bin/time -l bash $V/memcap.sh 290 600 uv run --no-project --with mpmath --with python-flint python3 code/certify.py 8 r8_00 r8_01 r8_02 r8_03 r8_04 r8_05 > logs/cert_r8.json 2> logs/cert_r8.err
echo "$(date -u +%H:%M:%S) certify rc=$?" >> logs/progress.txt
rm -f $V/r8_0*.dat
echo "$(date -u +%H:%M:%S) ALL DONE" >> logs/progress.txt

