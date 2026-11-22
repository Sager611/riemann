#!/bin/bash
# T4 memory probes (rerun from the repo root): bash Riemann/Deep/LeanPrefix/tools/probe.sh
# Each probe: a 2-line file importing one module, compiled by lean --threads=1 under .scratch/memcap.sh 290,
# sampled by tools/fpsamp.py (RSS and physical footprint = top's MEM).
cd /Users/adriansager/playground/riemann
L=$HOME/.elan/toolchains/leanprover--lean4---v4.35.0-rc2/bin/lean
export LEAN_PATH=$(ls -d $PWD/.lake/packages/*/.lake/build/lib/lean | paste -sd: -):$HOME/.elan/toolchains/leanprover--lean4---v4.35.0-rc2/lib/lean
d=/tmp/lp_probe; mkdir -p $d; date -u
for spec in "M:Init" "M:Mathlib.GroupTheory.Perm.Sign" "N:Mathlib.GroupTheory.Perm.Sign" "N:Mathlib.GroupTheory.SpecificGroups.Alternating" "N:Mathlib.GroupTheory.GroupAction.Primitive" "N:Mathlib.GroupTheory.GroupAction.Jordan" "N:Mathlib.GroupTheory.SpecificGroups.Dihedral" "N:Mathlib.GroupTheory.SemidirectProduct"; do
  kind=${spec%%:*}; mod=${spec#*:}; f=$d/probe_${kind}_${mod//./_}.lean
  if [ $kind = M ]; then printf 'module\n\npublic import %s\n\n#check 1\n' $mod > $f; else printf 'import %s\n\n#check 1\n' $mod > $f; fi
  echo "== $([ $kind = M ] && echo module || echo non-module) file: import $mod"
  python3 Riemann/Deep/LeanPrefix/tools/fpsamp.py "^$L --threads=1 $f" > $d/samp.out 2>&1 &
  sp=$!; sleep 0.3
  bash .scratch/memcap.sh 290 600 $L --threads=1 $f 2>&1 | rg -v '^1 : Nat|^Nat : Type'
  wait $sp; cat $d/samp.out
done

