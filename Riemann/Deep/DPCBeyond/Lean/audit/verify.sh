#!/bin/bash
# dpcblean final verification: every module in dependency order, one at a time, lean -j1 directly,
# memcap 290 MB / 1800 s, /usr/bin/time -l; pplean's FinSum, Model, Liminf and dpclean's modules first (imported, unmodified).
cd /Users/adriansager/playground/riemann
out=.scratch/dpcblean/verify.log; : > $out
rm -rf .scratch/dpcblean/olean; mkdir -p .scratch/dpcblean/olean
for m in FinSum Model Liminf; do
  LINES_MAX=5 bash .scratch/dpcblean/mk.sh $m Riemann/Deep/ExplorePointProc/Lean >> $out 2>&1
done
for m in Params Gaps Tail Phi Step7 TheoremU DLam Hairpin Fine; do
  LINES_MAX=5 bash .scratch/dpcblean/mk.sh $m Riemann/Deep/DPCUncond/Lean >> $out 2>&1
done
for m in ParamsB Rpow Step7B HairpinB SumB DLamB TheoremUB FinalB LemmaZ LemmaM LemmaM4 LemmaM3 FineB LemmaZ2 LemmaZW CheckB WiredB; do
  LINES_MAX=200 bash .scratch/dpcblean/mk.sh $m >> $out 2>&1
done
echo VERIFY-DONE >> $out

