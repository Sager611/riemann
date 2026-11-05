#!/bin/sh
# queue4.sh: full rigorous early-exit scan at c = 47 (K2 = 3683, Chekmachev's), resumable chunks; then comparison with his retained nodes
cd /Users/adriansager/playground/riemann/Riemann/Deep/WeilZeta
export WZ_MARGIN=1e-5 WZ_RSCAN2=/Users/adriansager/playground/riemann/.scratch/weilzeta/rscan2 WZ_SCAN=/Users/adriansager/playground/riemann/.scratch/weilzeta/scan2
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
R="nice -n 10 uv run --no-project --with mpmath --with numpy --with python-flint python3"
perl -e 'alarm 3500; exec @ARGV' /usr/bin/time -l $R scan_zeta.py 47 3683 24 >> logs/scan47.log 2>> logs/scan47.time
C=/Users/adriansager/playground/riemann/.scratch/weilzeta/chek
$R cmp_chek.py /Users/adriansager/playground/riemann/.scratch/weilzeta/scan2/scan_zeta_47_3683_rig.json $C/certificates__c47_even_eps__nodes_odd.csv $C/certificates__c47_even_eps__nodes_even.csv 3683 >> logs/scan47_cmp.log 2>&1
echo QUEUE4 DONE >> logs/scan47.log
