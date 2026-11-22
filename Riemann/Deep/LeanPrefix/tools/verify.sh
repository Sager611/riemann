#!/bin/bash
# Clean rebuild of Riemann/Deep/LeanPrefix (one lake env lean at a time, under leanrun.sh), then the
# axiom audit and the static scan. Usage (repo root): bash Riemann/Deep/LeanPrefix/tools/verify.sh
cd /Users/adriansager/playground/riemann
P=Riemann/Deep/LeanPrefix
rm -f .lake/build/lib/lean/Riemann/Deep/LeanPrefix/*.olean
MODS="Defs CycleA CycleB CycleC IsCycleP Parity BlocksA BlocksB Orbits Prim0 Bridge OrbitsB PrimA PrimB1 PrimB2 Prop5A Prop5B Prop5C Small Small23 Lemma7A Lemma7B Lemma7C Inputs Alt Wreath WreathA WreathB WreathC B2 Main Check"
date -u
bash $P/tools/mk.sh $MODS || exit 1
echo "== axioms"
rg "depends on axioms" /tmp/lp/logs/Check.log
echo "== static scan"
rg -n 'sorry|admit|native_decide|axiom|implemented_by|extern|unsafe' $P --glob '*.lean' || echo "(no matches)"
date -u

