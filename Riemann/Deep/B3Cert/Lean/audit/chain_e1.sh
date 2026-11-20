#!/bin/bash
cd /Users/adriansager/playground/riemann
for m in Defs Detector AlgebraH AlgebraC Numeric Kappa Arch Parity Bounds2; do
  bash Riemann/Deep/B3Cert/Lean/audit/mk.sh Riemann/Deep/SantiAudit/Lean $m > Riemann/Deep/B3Cert/Lean/audit/chain_$m.out 2>&1 || { echo "FAIL $m" >> Riemann/Deep/B3Cert/Lean/audit/chain.log; exit 1; }
  echo "ok $m" >> Riemann/Deep/B3Cert/Lean/audit/chain.log
done
echo DONE >> Riemann/Deep/B3Cert/Lean/audit/chain.log
