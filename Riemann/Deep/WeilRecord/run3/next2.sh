#!/bin/sh
# Disk-lean follow-up to e3run.sh (scratch kept < 2 GB).  One heavy process at a time, niced.  Log: next2.log
cd /Users/adriansager/playground/riemann/.scratch/weilrec3
PY=../weil/venv/bin/python; T="/usr/bin/time -l"; W2=../weilrec2
while pgrep -f e3run.sh >/dev/null; do sleep 60; done
echo "e3run finished $(date): $(rg -o 'RESULT.*' cert_e3.log | head -1)"
rm -f band_e3_3000_1216_103o20_384.bin
nice -n 10 $T $PY cert3.py 16 $W2/grams_16_1400_384.bin $W2/band_16_1400_640_9o2_384.bin $W2/w_16_32000.txt 36421/10000 1/2 1/1$(printf '0%.0s' $(seq 1 76)) 343 114 384 - cert3_16_summary.json > cert3_16.log 2>&1
echo "cert3 mu=16 exit $? $(date): $(rg -o 'RESULT.*' cert3_16.log | head -1)"
rg -q 'RESULT.*CERTIFIED' cert_e3.log || { echo "e3 not certified: stop"; exit 1; }
nice -n 10 $T ./grams_c e3 3000 384 0/1 ctest/gpe3_bound.bin bound > grams_c_e3_bound.log 2>&1 || { echo "grams_c e3 failed"; exit 1; }
nice -n 10 $PY cmp_grams.py ctest/gpe3_bound.bin grams_e3_3000_384.bin 384 > cmp_e3_bound.log 2>&1
NOK=$(rg -c 'non-overlapping=0 ' cmp_e3_bound.log)
echo "bound-vs-exact e3 $(date): $NOK/4 keys overlap"
rm -f ctest/gpe3_bound.bin grams_e3_3000_384.bin
[ "$NOK" = "4" ] || { echo "bound mode disagrees: stop"; exit 1; }
nice -n 10 $T ./grams_c L31o20 4200 448 0/1 gramspart_L31o20_4200_448_0of1.bin bound > grams_L155.log 2>&1 || { echo "grams L155 failed"; exit 1; }
echo "grams L155 done $(date)"
nice -n 10 $T $PY cert3.py L31o20 gramspart_L31o20_4200_448_0of1.bin compute:1640:109/20 w_L31o20_32000.txt 45913/10000 1/2 1/1$(printf '0%.0s' $(seq 1 112)) 913 303 448 - cert_L155_summary.json > cert_L155.log 2>&1
echo "cert L155 exit $? $(date): $(rg -o 'RESULT.*' cert_L155.log | head -1)"
rm -f gramspart_L31o20_4200_448_0of1.bin
