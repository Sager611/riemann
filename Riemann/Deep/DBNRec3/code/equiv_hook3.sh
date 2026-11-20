#!/bin/bash
# dbnrec3 equiv_hook3.sh: equivalence test of bb_rec3 (two-group prefilter, .scratch/dbnrec3/bb_new) against the bb_rec2 build, run by master3.sh as a hook
# between rows (one compute process at a time): (1) old vs new on shell r8_05 (outer, thr 0.939035) with K = 6e8 anchors (4 progress lines + DONE), speed;
# (2) new on the full shell r8_00 (shell 0, thr 0.913118, 350 dumped cells) against tranche 2's bb.out and unres dump (bb_rec2 build), same arguments as
# run_row.sh; times stripped; (3) mkrec.py path check; (4) adoption only if (1) and (2) are IDENTICAL and new is faster in (1) (bb's own CPU seconds, DONE line).
cd /Users/adriansager/playground/riemann
L=Riemann/Deep/DBNRec3/logs; S=.scratch/dbnrec3; E=$L/equiv.txt
strip() { sed -E 's/\| [0-9]+ s$//; s/tables [0-9.]+ s//' "$1" | rg -v 'memcap|maximum resident|peak memory'; }
cpu() { awk '/ user /{for(i=1;i<NF;i++){if($(i+1)=="user")u=$i; if($(i+1)=="sys")s=$i}} END{printf "%.2f", u+s}' "$1"; }
echo "$(date -u +%H:%M:%S) equivalence test start (bb_new sha256 $(shasum -a 256 $S/bb_new | cut -c1-16), bb_rec sha256 $(shasum -a 256 $S/bb_rec | cut -c1-16))" >> $E
for b in bb_rec bb_new; do
  /usr/bin/time -l bash .scratch/memcap.sh 290 300 $S/$b $S/r8_05.dat 48 400 600000000 0 0.001 $S/eq_k_$b.unres.txt > $S/eq_k_$b.bb.out 2> $S/eq_k_$b.err
  echo "$b on r8_05 (K = 6e8): $(rg -o '[0-9]+ s$' $S/eq_k_$b.bb.out | tail -1) | user+sys $(cpu $S/eq_k_$b.err) s | $(rg 'memcap' $S/eq_k_$b.err)" >> $E
done
if diff <(strip $S/eq_k_bb_rec.bb.out) <(strip $S/eq_k_bb_new.bb.out) > /dev/null && cmp -s $S/eq_k_bb_rec.unres.txt $S/eq_k_bb_new.unres.txt; then echo "  r8_05 K=6e8: old and new IDENTICAL ($(wc -l < $S/eq_k_bb_new.bb.out | tr -d ' ') lines)" >> $E; else echo "  r8_05 K=6e8: DIFFER" >> $E; fi
bbs() { rg -o '[0-9]+ s$' "$1" | tail -1 | tr -dc 0-9; }   # bb's own CPU seconds (DONE line, clock())
FASTER=$(awk -v o=$(bbs $S/eq_k_bb_rec.bb.out) -v n=$(bbs $S/eq_k_bb_new.bb.out) 'BEGIN{print (n + 0 < o + 0) ? 1 : 0}')
tag=r8_00
/usr/bin/time -l bash .scratch/memcap.sh 290 900 $S/bb_new $S/$tag.dat 48 400 0 0 0.001 $S/eq_$tag.unres.txt > $S/eq_$tag.bb.out 2> $S/eq_$tag.err
echo "new on $tag (full shell): rc=$? $(rg -o '[0-9]+ s$' $S/eq_$tag.bb.out | tail -1) (tranche 2, bb_rec2 build: $(rg -o '[0-9]+ s$' $L/$tag.bb.out | tail -1)) | user+sys $(cpu $S/eq_$tag.err) s | $(rg 'memcap' $S/eq_$tag.err) $(rg 'maximum resident' $S/eq_$tag.err | sed -E 's/ +/ /g')" >> $E
if diff <(strip $S/eq_$tag.bb.out) <(strip $L/$tag.bb.out) > /dev/null && cmp -s $S/eq_$tag.unres.txt $L/$tag.unres.txt; then echo "  $tag: bb.out (times stripped) and unres dump IDENTICAL to tranche 2 ($(wc -l < $S/eq_$tag.unres.txt | tr -d ' ') cells)" >> $E
else echo "  $tag: DIFFERS from tranche 2" >> $E; diff <(strip $S/eq_$tag.bb.out) <(strip $L/$tag.bb.out) | head -5 >> $E; fi
# mkrec.py path check: tranche 1-2 certs (rows 1-8) must rebuild chain F's source (sha256 6321ae69...9aab)
UVF="uv run --no-project --with mpmath --with python-flint"
bash .scratch/memcap.sh 290 300 $UVF python3 Riemann/Deep/DBNRec3/code/mkrec.py source Riemann/Deep/DBNRec/logs/cert_r1.json Riemann/Deep/DBNRec/logs/cert_r2.json Riemann/Deep/DBNRec/logs/cert_r3.json \
  Riemann/Deep/DBNRec2/logs/cert_r4.json Riemann/Deep/DBNRec2/logs/cert_r5.json Riemann/Deep/DBNRec2/logs/cert_r6.json Riemann/Deep/DBNRec2/logs/cert_r7.json Riemann/Deep/DBNRec2/logs/cert_r8.json > $S/mkrec_F.json 2>&1
echo "mkrec rc=$? sha256 $(shasum -a 256 $S/ws/dbnrec/MAIN_CLOCK_SOURCE_REC.json | cut -c1-64) (tranche 2 chain F: 6321ae690ef0760caa98808bb5fc05af5e08a0113b96046046a538e5e4109aab)" >> $E
# adoption (reviewer pass 2 CORRECT, 15:03 UTC): atomic rename, old build kept as bb_rec2_build
if [ $(rg -c 'IDENTICAL' $E) -ge 2 ] && ! rg -q 'DIFFER' $E && [ "$FASTER" = 1 ]; then cp $S/bb_rec $S/bb_rec2_build && cp $S/bb_new $S/bb_rec.tmp && mv $S/bb_rec.tmp $S/bb_rec && echo "$(date -u +%H:%M:%S) ADOPTED: .scratch/dbnrec3/bb_rec is now the bb_rec3 build (sha256 $(shasum -a 256 $S/bb_rec | cut -c1-16))" >> $E
else echo "$(date -u +%H:%M:%S) NOT adopted (identical: $(rg -c 'IDENTICAL' $E), faster: $FASTER); bb_rec unchanged" >> $E; fi
echo "$(date -u +%H:%M:%S) equivalence test end" >> $E
