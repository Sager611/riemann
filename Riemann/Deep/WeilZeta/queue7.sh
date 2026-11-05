#!/bin/sh
# queue7.sh: (0) checks of Chekmachev's c = 45, 47 retained sets (initial segments, symbol spot check);
# after the memory fixes (no midpoint copy of Id - W_up, no copy of L21 L21^T): re-validate (c = 16 even, c = 20 odd: same margins
# expected), then c = 29 at Chekmachev's eps with 576 bits (the 640-bit even run passed 1.5 GB at 17 min and was stopped)
cd /Users/adriansager/playground/riemann/Riemann/Deep/WeilZeta
export WZ_RSCAN=/Users/adriansager/playground/riemann/.scratch/weilzeta/rscan WZ_SCAN=/Users/adriansager/playground/riemann/.scratch/weilzeta/scan WZ_CHOL_BS=64
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
R="nice -n 10 uv run --no-project --with mpmath --with numpy --with python-flint python3"
C=/Users/adriansager/playground/riemann/.scratch/weilzeta/chek
perl -e 'alarm 600; exec @ARGV' $R chek_spot.py 47 3683 $C/certificates__c47_even_eps__nodes_odd.csv $C/certificates__c47_even_eps__nodes_even.csv 1000 >> logs/chek_spot.log 2>&1
perl -e 'alarm 600; exec @ARGV' $R chek_spot.py 45 3375 $C/certificates__c45_even_eps__nodes_odd.csv $C/certificates__c45_even_eps__nodes_even.csv 300 >> logs/chek_spot.log 2>&1
./memwatch.sh &
perl -e 'alarm 900; exec @ARGV' /usr/bin/time -l $R cert_zeta.py 16 even 1.512800554e-75 512 >> logs/val3.log 2>> logs/val3.time
perl -e 'alarm 900; exec @ARGV' /usr/bin/time -l $R cert_zeta.py 20 odd 2^-327 448 >> logs/val3.log 2>> logs/val3.time
perl -e 'alarm 2700; exec @ARGV' /usr/bin/time -l $R cert_zeta.py 29 even 3.3887641037680338350e-145 576 >> logs/z29.log 2>> logs/z29.time
perl -e 'alarm 2700; exec @ARGV' /usr/bin/time -l $R cert_zeta.py 29 odd 1.5329466444386175069e-140 576 >> logs/z29.log 2>> logs/z29.time
echo QUEUE7 DONE >> logs/z29.log
