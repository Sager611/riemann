#!/bin/bash
# T1 import-memory probes (2026-11-22). Run from the repo root: bash Riemann/Deep/LeanNBRate/audit/probe.sh
# Every lean/lake run under .scratch/memcap.sh 290 (kills on ps RSS, 0.5 s samples); for direct lean runs
# fp.py samples lean's RSS and physical footprint (top's MEM) every 20 ms alongside.
A=Riemann/Deep/LeanNBRate/audit; T=$(mktemp -d /tmp/nbprobe.XXXX); Z=Mathlib.NumberTheory.LSeries.RiemannZeta
L=$HOME/.elan/toolchains/leanprover--lean4---v4.35.0-rc2/bin/lean
LP=$PWD/.lake/build/lib/lean; for d in .lake/packages/*/.lake/build/lib/lean; do LP=$LP:$PWD/$d; done
cap() { bash .scratch/memcap.sh 290 600 "$@" 2>&1 | tail -1; }
mk() { printf "$2" > $T/$1.lean; echo "== $1: $(tr '\n' ';' < $T/$1.lean)"; }
direct() { mk "$@"
  (bash .scratch/memcap.sh 60 120 python3 $A/fp.py "^$L --threads=1 $T/$1.lean" > $T/$1.fp 2>&1 &); sleep 0.3
  LEAN_PATH=$LP cap $L --threads=1 $T/$1.lean; sleep 1; head -1 $T/$1.fp; }
viaLake() { mk "$@"; cap lake env lean --threads=1 $T/$1.lean; }
date -u; echo "== lake env true"; cap lake env true
viaLake lake_zeta "import $Z\n#check RiemannHypothesis\n"
viaLake lake_zeta_mod "module\nimport $Z\n#check RiemannHypothesis\n"
direct init "#check Nat\n"
direct init_mod "module\n#check Nat\n"
direct lean_mod "module\nimport Lean\n"
direct order_mod "module\nimport Mathlib.Order.Basic\n"
direct zeta "import $Z\n#check RiemannHypothesis\n"
direct zeta_mod "module\nimport $Z\n#check RiemannHypothesis\n"
for m in Mathlib.MeasureTheory.Function.LpSpace.Basic Mathlib.NumberTheory.Harmonic.EulerMascheroni Mathlib.Topology.MetricSpace.HausdorffDistance; do
  direct "$(echo $m | tr . _)_mod" "module\nimport $m\n"; done
rm -rf $T
