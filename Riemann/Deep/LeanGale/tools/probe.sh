#!/bin/bash
# Import probes: physical footprint (top's MEM metric, lifetime max) of `lean -j1` on a file that only imports M,
# plain and as a `module`, each under fpcap.py's 290 MiB cap.  Run from anywhere: bash Riemann/Deep/LeanGale/tools/probe.sh
cd "$(dirname "$0")/../../../.." || exit 1
TC=$(ls -d ~/.elan/toolchains/leanprover--lean4---v4.35.0-rc2)
LP=$PWD/.lake/build/lib/lean; for p in .lake/packages/*; do LP="$LP:$PWD/$p/.lake/build/lib/lean"; done
d=$(mktemp -d)
for m in ${@:-Mathlib.Order.Basic Mathlib.Data.Real.Basic Mathlib.Topology.MetricSpace.Bounded Mathlib.Analysis.Real.Sqrt Mathlib.Analysis.InnerProductSpace.Basic Mathlib.Analysis.Normed.Lp.PiLp Mathlib.Analysis.InnerProductSpace.PiL2 Mathlib.NumberTheory.LSeries.RiemannZeta}; do
  printf 'import %s\n' "$m" > "$d/a.lean"; printf 'module\n\nimport %s\n' "$m" > "$d/b.lean"
  for f in a b; do
    r=$(LEAN_PATH=$LP uv run --no-project python3 Riemann/Deep/LeanGale/tools/fpcap.py 290 600 "$TC/bin/lean" -j1 "$d/$f.lean" 2>&1 | rg -o 'fpcap: .*')
    echo "$m $([ $f = a ] && echo plain || echo module): $r"
  done
done
rm -rf "$d"
