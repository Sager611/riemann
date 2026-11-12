#!/bin/bash
# build the imported project modules (pplean, dpclean) one at a time
cd /Users/adriansager/playground/riemann
for m in FinSum Model Liminf; do LINES_MAX=5 bash .scratch/dpcblean/mk.sh $m Riemann/Deep/ExplorePointProc/Lean; done
for m in Params Gaps Tail Phi Step7 TheoremU DLam Hairpin Fine; do LINES_MAX=5 bash .scratch/dpcblean/mk.sh $m Riemann/Deep/DPCUncond/Lean; done
echo DEPS-DONE
