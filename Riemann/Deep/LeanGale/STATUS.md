REVIEW
UTC: 2026-11-22T13:16:42Z
sorry count: 0 (rg 'sorry|admit|native_decide|axiom|implemented_by|extern|unsafe' Riemann/Deep/LeanGale/*.lean: only the #print axioms lines of Check.lean)
files: Gale.lean (definitions: E, beta, GaleConjecture, GaleConjectureBall; bridge to EuclideanSpace/Metric.diam), Main.lean (claims and main theorems), Family.lean, Row.lean, Splx.lean, Pts.lean, Cert.lean, P3/P4/P7/P8.lean, Data/*.lean (generated, tools/gen.py, tools/rows.py; tools/checkdata.py: all equal to the source files), Check.lean
main theorem: LeanGale.gale_fails (hHin : Hinrichs2002) : ∀ n ≥ 3, ¬ GaleConjecture n; also gale_fails_le_1000 (n = 3…1000, no hypotheses), claim1_n7, claim1_n8, claim2, gale_ball_fails
literature inputs: Hinrichs2002 (A. Hinrichs, Spherical codes and Borsuk's conjecture, Discrete Math. 243 (2002) 253–256, Thm 1: Borsuk fails in all dimensions > 322), only for n > 1000; SimplexPartitionBound (Voronov arXiv 2610.28828 l. 96–107, unrefereed), only for the ball-form corollary gale_ball_fails
#print axioms (all six main theorems and splx_simplex): [propext, Classical.choice, Quot.sound]
NOT-FORMALIZED: none
last reviewer verdict: none yet
next step: reviewer rounds (statement faithfulness, literature inputs, vacuity)

## Build (one module at a time, .scratch/lean4/leanrun.sh, lake build Riemann.Deep.LeanGale.<M>)
Peaks: Rows00–19 2.4 GB, 13–18 s each; P3 4.7 GB 282 s; P4 5.1 GB 152 s; P7 3.2 GB 93 s; P8 2.5 GB 41 s; Gale (PiL2) 2.4 GB; Family 1.7 GB; Main 2.4 GB.
