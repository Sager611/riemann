#!/bin/sh
# half-width 1.55 grams (weilrec4): two full-parity parts with balanced explicit piece lists, niced, 4 h timeout each.
cd /Users/adriansager/playground/riemann/.scratch/weilrec4
TO="perl -e 'alarm shift; exec @ARGV' 14400"
T="/usr/bin/time -l"
GRAMS_PIECES=0,2,7,8,9,10,11 nice -n 10 perl -e 'alarm shift; exec @ARGV' 14400 $T ./grams_c L31o20 4200 448 0/2 g155_A.bin bound > g155_A.log 2>&1 &
GRAMS_PIECES=1,3,4,5,6,12 nice -n 10 perl -e 'alarm shift; exec @ARGV' 14400 $T ./grams_c L31o20 4200 448 1/2 g155_B.bin bound > g155_B.log 2>&1 &
wait
echo "grams 1.55 done $(date): $(tail -n 30 g155_A.log | rg -o 'grams_c done.*') | $(tail -n 30 g155_B.log | rg -o 'grams_c done.*')"
