#!/bin/bash
# usage (repo root): bash Riemann/Deep/LeanPrefix/tools/mk.sh Mod [Mod ...]
# Compiles Riemann/Deep/LeanPrefix/<Mod>.lean one at a time (lake env lean --threads=1) under
# .scratch/lean4/leanrun.sh; writes the .olean to .lake/build/lib/lean; log in /tmp/lp/logs/<Mod>.log.
cd /Users/adriansager/playground/riemann
mkdir -p .lake/build/lib/lean/Riemann/Deep/LeanPrefix /tmp/lp/logs
for M in "$@"; do
  src=Riemann/Deep/LeanPrefix/$M.lean
  out=.lake/build/lib/lean/Riemann/Deep/LeanPrefix/$M.olean
  s=$SECONDS
  bash .scratch/lean4/leanrun.sh -- lake env lean --threads=1 -o "$out" "$src" > /tmp/lp/logs/$M.log 2>&1
  rc=$?
  echo "$M rc=$rc $((SECONDS - s))s $(rg -o 'memcap: .*' /tmp/lp/logs/$M.log | tail -1)"
  if [ $rc -ne 0 ]; then rg -v '^memcap' /tmp/lp/logs/$M.log | head -80; exit $rc; fi
done
