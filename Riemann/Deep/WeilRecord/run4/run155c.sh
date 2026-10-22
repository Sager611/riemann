#!/bin/sh
# 1.55 grams, part A rerun (the first g155_A.bin was truncated by a full disk): one parity per process.
cd /Users/adriansager/playground/riemann/.scratch/weilrec4
T="/usr/bin/time -l"
GRAMS_PARITY=0 GRAMS_PIECES=0,2,7,8,9,10,11 nice -n 10 $T perl -e 'alarm shift; exec @ARGV' 18000 ./grams_c L31o20 4200 448 0/2 g155_A0.bin bound > g155_A0.log 2>&1 &
GRAMS_PARITY=1 GRAMS_PIECES=0,2,7,8,9,10,11 nice -n 10 $T perl -e 'alarm shift; exec @ARGV' 18000 ./grams_c L31o20 4200 448 0/2 g155_A1.bin bound > g155_A1.log 2>&1 &
wait
rg -q 'grams_c done' g155_A0.log && rg -q 'grams_c done' g155_A1.log || { echo "A rerun failed $(date)"; exit 1; }
echo "A rerun done $(date)"
TAU=1/1$(printf '0%.0s' $(seq 1 112)); PY=../weil/venv/bin/python; G=g155_B.bin,g155_A0.bin,g155_A1.bin
DIAG=1 nice -n 10 $T perl -e 'alarm shift; exec @ARGV' 18000 $PY cert4.py L31o20 $G 1640 109/20 w_L31o20_32000.txt 45913/10000 1/2 $TAU 913 303 448 3400 c155_R1.json > c155_R1.log 2>&1 &
DIAG=1 nice -n 10 $T perl -e 'alarm shift; exec @ARGV' 18000 $PY cert4.py L31o20 $G 1340 131/25 w_L31o20_32000.txt 45913/10000 1/2 $TAU 798 266 448 2825 c155_R2.json > c155_R2.log 2>&1 &
wait
echo "R1 $(date): $(rg -o 'RESULT.*' c155_R1.log)"
echo "R2: $(rg -o 'RESULT.*' c155_R2.log)"
