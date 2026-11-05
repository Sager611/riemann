#!/bin/sh
# queue6.sh: full rigorous early-exit scan at c = 45 (K2 = 3375, Chekmachev's), then comparison with his retained nodes
cd /Users/adriansager/playground/riemann/Riemann/Deep/WeilZeta
export WZ_MARGIN=1e-5 WZ_RSCAN2=/Users/adriansager/playground/riemann/.scratch/weilzeta/rscan2 WZ_SCAN=/Users/adriansager/playground/riemann/.scratch/weilzeta/scan2
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
R="nice -n 10 uv run --no-project --with mpmath --with numpy --with python-flint python3"
perl -e 'alarm 3000; exec @ARGV' /usr/bin/time -l $R scan_zeta.py 45 3375 24 >> logs/scan45.log 2>> logs/scan45.time
C=/Users/adriansager/playground/riemann/.scratch/weilzeta/chek
$R cmp_chek.py /Users/adriansager/playground/riemann/.scratch/weilzeta/scan2/scan_zeta_45_3375_rig.json $C/certificates__c45_even_eps__nodes_odd.csv $C/certificates__c45_even_eps__nodes_even.csv 3375 >> logs/scan45_cmp.log 2>&1
echo QUEUE6 DONE >> logs/scan45.log
