#!/bin/bash
# b3lean verify.sh: clean rebuild of E1Lean's modules (read-only sources) and E1PLean into a fresh olean dir, one lean at a time.
cd /Users/adriansager/playground/riemann
O=/Users/adriansager/playground/riemann/.scratch/b3lean/olean_v
rm -rf $O; mkdir -p $O
LP="$O:$(cat Riemann/Deep/B3Cert/Lean/audit/leanpath.txt)"
LEANBIN=/Users/adriansager/.elan/toolchains/leanprover--lean4---v4.35.0-rc2/bin/lean
: > Riemann/Deep/B3Cert/Lean/audit/verify.log
for dm in SantiAudit:Defs SantiAudit:Detector SantiAudit:AlgebraH SantiAudit:AlgebraC SantiAudit:Numeric SantiAudit:Kappa SantiAudit:Arch SantiAudit:Parity SantiAudit:Bounds2 B3Cert:Trig B3Cert:Exp B3Cert:Num B3Cert:Arch B3Cert:Window B3Cert:Parity B3Cert:DetTrig B3Cert:Det B3Cert:Check; do
  d=Riemann/Deep/${dm%%:*}/Lean; m=${dm##*:}
  mkdir -p $O/$d; log=Riemann/Deep/B3Cert/Lean/audit/logs/vlog_${dm%%:*}_$m.txt; t0=$(date +%s)
  LEAN_PATH="$LP" OMP_NUM_THREADS=1 bash .scratch/memcap.sh 290 1800 /usr/bin/time -l $LEANBIN -j1 -o $O/$d/$m.olean $d/$m.lean > $log 2>&1
  rc=$?
  pk=$(rg -o '[0-9]+  peak memory footprint' $log | awk '{printf "%d", $1/1048576}')
  mc=$(rg -o 'memcap: peak [0-9]+ MB' $log)
  echo "$d/$m rc=$rc $(( $(date +%s) - t0 ))s footprint=${pk}MiB $mc" >> Riemann/Deep/B3Cert/Lean/audit/verify.log
  [ $rc -ne 0 ] && { echo FAIL >> Riemann/Deep/B3Cert/Lean/audit/verify.log; exit 1; }
done
echo DONE >> Riemann/Deep/B3Cert/Lean/audit/verify.log
