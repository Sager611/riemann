#!/bin/sh
# log 5 (x = 25) phase (weilrec4): waits for run155c.sh; picks parameters from the 1.55 results:
#  B (gap 0.65, Omega 1920, N 4210) if R2 (1.55 at gap 0.65) certified, else A (gap 0.86, Omega 2380, N 5210) if R1 certified.
cd /Users/adriansager/playground/riemann/.scratch/weilrec4
PY=./venv/bin/python; T="/usr/bin/time -l"
while pgrep -f run155c.sh >/dev/null; do sleep 120; done
if rg -q 'RESULT.*CERTIFIED' c155_R2.log; then SEL=B; N=4210; OM=1920; GA=561/100; U0=1087; U1=362;
elif rg -q 'RESULT.*CERTIFIED' c155_R1.log; then SEL=A; N=5210; OM=2380; GA=291/50; U0=1405; U1=468;
else echo "no 1.55 certificate: stop $(date)"; exit 1; fi
echo "log5 selection $SEL N=$N Omega=$OM gamma=$GA u=($U0,$U1) $(date)"
rm -f g155_B.bin g155_A0.bin g155_A1.bin
GRAMS_PARITY=0 nice -n 10 $T perl -e 'alarm shift; exec @ARGV' 64800 ./grams_c 25 $N 480 0/1 g25_p0.bin bound > g25_p0.log 2>&1 &
GRAMS_PARITY=1 nice -n 10 $T perl -e 'alarm shift; exec @ARGV' 64800 ./grams_c 25 $N 480 0/1 g25_p1.bin bound > g25_p1.log 2>&1 &
wait
rg -q 'grams_c done' g25_p0.log && rg -q 'grams_c done' g25_p1.log || { echo "log5 grams failed $(date)"; exit 1; }
echo "log5 grams done $(date)"
TAU=1/1$(printf '0%.0s' $(seq 1 127))
PREC_C=640 nice -n 10 $T perl -e 'alarm shift; exec @ARGV' 64800 $PY cert4.py 25 g25_p0.bin,g25_p1.bin $OM $GA w_25_32000.txt 49604/10000 1/2 $TAU $U0 $U1 480 $N c25_$SEL.json > c25_$SEL.log 2>&1
echo "log5 cert exit $? $(date): $(rg -o 'RESULT.*' c25_$SEL.log)"
