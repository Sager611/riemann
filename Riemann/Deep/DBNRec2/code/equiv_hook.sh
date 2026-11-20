#!/bin/bash
# dbnrec2 equivalence test of bb_rec2 (prefilter) against tranche 1's output and the unchanged bb_rec (run by master2.sh between rows: one compute process at a time).
cd /Users/adriansager/playground/riemann
L=Riemann/Deep/DBNRec2/logs; S=.scratch/dbnrec2; E=$L/equiv.txt
strip() { sed -E 's/\| [0-9]+ s$//; s/tables [0-9.]+ s//' "$1" | rg -v 'memcap|maximum resident|peak memory'; }
echo "$(date -u +%H:%M:%S) equivalence test start" >> $E
/usr/bin/time -l bash .scratch/memcap.sh 290 600 $S/bb_new $S/r4_00.dat 48 400 0 0 0.001 $S/eq_r4_00.unres.txt > $S/eq_r4_00.bb.out 2> $S/eq_r4_00.err
echo "new on r4_00 (full shell, thr 0.832387): rc=$? $(rg -o '[0-9]+ s$' $S/eq_r4_00.bb.out | tail -1) $(rg 'memcap' $S/eq_r4_00.err) $(rg 'maximum resident' $S/eq_r4_00.err)" >> $E
if diff <(strip $S/eq_r4_00.bb.out) <(strip $L/r4_00.bb.out) > /dev/null && cmp -s $S/eq_r4_00.unres.txt $L/r4_00.unres.txt; then echo "  r4_00: bb.out (times stripped) and unres dump IDENTICAL to tranche 1 ($(wc -l < $S/eq_r4_00.unres.txt) cells)" >> $E
else echo "  r4_00: DIFFERS from tranche 1" >> $E; diff <(strip $S/eq_r4_00.bb.out) <(strip $L/r4_00.bb.out) | head -5 >> $E; fi
for b in bb_rec bb_new; do
  /usr/bin/time -l bash .scratch/memcap.sh 290 300 $S/$b $S/r4_06.dat 48 400 300000000 0 0.001 $S/eq_$b.unres.txt > $S/eq_$b.bb.out 2> $S/eq_$b.err
  echo "$b on r4_06 (K = 3e8): $(rg -o '[0-9]+ s$' $S/eq_$b.bb.out | tail -1) $(rg 'real' $S/eq_$b.err | sed -E 's/ +/ /g')" >> $E
done
if diff <(strip $S/eq_bb_rec.bb.out) <(strip $S/eq_bb_new.bb.out) > /dev/null && cmp -s $S/eq_bb_rec.unres.txt $S/eq_bb_new.unres.txt; then echo "  r4_06 K=3e8: old and new IDENTICAL" >> $E; else echo "  r4_06 K=3e8: DIFFER" >> $E; fi
echo "$(date -u +%H:%M:%S) equivalence test end" >> $E
