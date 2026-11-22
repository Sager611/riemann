#!/bin/bash
# Reproduces the BLOCKED-MEMCAP measurements. Every run is under .scratch/memcap.sh 290 (RSS cap);
# tools/fp.c samples RSS and physical footprint of the run's process every 10 ms until memcap stops it.
cd "$(git rev-parse --show-toplevel)" || exit 1
. Riemann/Deep/LeanLi/tools/env.sh
d=Riemann/Deep/LeanLi
cc -O2 -o /tmp/leanli_fp $d/tools/fp.c || exit 1
run() { # label cmd...
  local label=$1; shift
  bash .scratch/memcap.sh 290 600 "$@" > /tmp/leanli_run.out 2>&1 &
  local mc=$! pid=""
  for i in 1 2 3 4 5 6 7 8 9 10; do pid=$(pgrep -P $mc | head -1); [ -n "$pid" ] && break; sleep 0.01; done
  echo "== $label"
  [ -n "$pid" ] && /tmp/leanli_fp "$pid" 10 1000000 | tail -1
  wait $mc; echo "exit $?; $(tail -1 /tmp/leanli_run.out)"
}
mkdir -p $d/probe
for m in Init Mathlib.Order.Basic Mathlib.Data.Real.Basic Mathlib.NumberTheory.LSeries.RiemannZeta; do
  printf 'module\npublic import %s\n#check Nat\n' "$m" > $d/probe/Import.lean
  run "lean --threads=1, module file: public import $m" "$LEAN" --threads=1 $d/probe/Import.lean
done
printf 'import Mathlib.NumberTheory.LSeries.RiemannZeta\n#check RiemannHypothesis\n' > $d/probe/Import.lean
run "lean --threads=1, non-module file: import Mathlib.NumberTheory.LSeries.RiemannZeta" "$LEAN" --threads=1 $d/probe/Import.lean
run "lake env true" lake env true
for m in Init Mathlib.Order.Basic Mathlib.Data.Real.Basic Mathlib.NumberTheory.LSeries.RiemannZeta; do
  python3 $d/tools/cone.py $m
done
