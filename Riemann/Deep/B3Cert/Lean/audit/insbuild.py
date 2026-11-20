
p = "Riemann/Deep/B3Cert/LEAN.md"; s = open(p).read()
tab = open("Riemann/Deep/B3Cert/Lean/audit/build_table.md").read().replace("Scripts (Python, exact Fractions): ≤ 1 s each, below memcap's first sample.", "Scripts (Python, exact Fractions, single-threaded) each ended before memcap's first sample (peak 0 MB).")
block = """## Build
Mathlib: \`rm .lake && bash .scratch/memcap.sh 290 1500 lake exe cache get\` (20 Nov 05:42 UTC; 52 s, peak 246 MB, rc = 0). LEAN_PATH = \`lake env printenv LEAN_PATH\` (Riemann/Deep/B3Cert/Lean/audit/leanpath.txt). Imports: E1Lean's modules only (Numeric, Bounds2, Parity) plus our own; no new Mathlib import. Each module separately, one at a time, lean directly (no lake build), oleans to .scratch/b3lean/olean (deleted at the end):

    LEAN_PATH=.scratch/b3lean/olean:$(cat Riemann/Deep/B3Cert/Lean/audit/leanpath.txt) OMP_NUM_THREADS=1 bash .scratch/memcap.sh 290 1800 \\\\
      /usr/bin/time -l ~/.elan/toolchains/leanprover--lean4---v4.35.0-rc2/bin/lean -j1 \\\\
      -o .scratch/b3lean/olean/$d/$m.olean $d/$m.lean

Scripts: Riemann/Deep/B3Cert/Lean/audit/mk.sh <dir> <Module> (one module), verify.sh (clean rebuild into a fresh olean dir, E1Lean's 9 modules then ours, order Trig, Exp, Num, Arch, Window, Parity, DetTrig, Det, Check). Final clean rebuild (20 Nov 06:05–06:15 UTC, system load ≈ 10 from other agents, background QoS), every module rc = 0:

""" + tab + """
## Static scan
\`rg -n -w 'sorry|admit|native_decide|skipKernelTC|decide' Riemann/Deep/B3Cert/Lean/*.lean\` and \`rg -n '^\\\\s*axiom\\\\s' Riemann/Deep/B3Cert/Lean/*.lean\`: no match (rc = 1). No errors or warnings in any compile log.

## #print axioms (Check.lean, 37 lines, all identical)

    'E1PLean.theoremE1P' depends on axioms: [propext, Classical.choice, Quot.sound]
    (same line for theoremE1P_gauss, theoremE1P_parity, theoremE1P_fine, kappa2_gt_kstar, e1p_enclosures, detector_threshold, inputs_consistentP, gauss_consistent, H1P_of_gauss, H3P_of_inputs, eta_edge, eta_step, nu_edge_iff, nu_edge_eq, exponent_edge, cW_edge_ge, cW_edge_le, kappa2_bounds, kappa2_upper, h_edge_gt, h_edge_le, kstar_edge_bounds, exp_Lhi2, exp_Llo2, sq_lo, sq_up, base_sin, cos_two, trig_xlo2, trig_xhi2, trig_det_lo, trig_det_hi, cW_detp_ge, cW_detm_le, kstar_anti, kstar_det_bounds)

"""
i = s.index("## Review")
s = s[:i] + block + s[i:]
open(p, "w").write(s)
print("inserted", block.count("\n"))

