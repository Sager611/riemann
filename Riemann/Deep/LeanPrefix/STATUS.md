REVIEW
UTC: 2026-11-22 13:58
Sorry count: 0 (rg 'sorry|admit|native_decide|axiom|implemented_by|extern|unsafe': only the #print axioms lines of Check.lean)
Files: Riemann/Deep/LeanPrefix/*.lean (32 modules; Main.lean = statement), tools/verify.sh (clean rebuild + axioms + scan), audit/ (review prompts and verdicts), PLAN.md
Main theorem: Riemann.Deep.LeanPrefix.prefix_reversals (Main.lean); #print axioms: [propext, Classical.choice, Quot.sound]
Literature inputs (hypotheses, Inputs.lean): LZZW2026_Thm3_1 = Lei-Zhao-Zhu-Wu, arXiv 2610.30892v1 (unrefereed), Thm 3.1, first sentence (TeX l. 221-223); Jones2014_Thm1_2_k0 = G. A. Jones, Bull. Aust. Math. Soc. 89 (2014) 159-165, arXiv 1209.5169, Thm 1.2 case k = 0 (TeX l. 54-72), conclusions weakened (lower bounds C_p ≤ G, PGL_d(q) ≤ G and the identification in (c) dropped)
NOT-FORMALIZED: none
Proved (ours, no inputs): Prop 5 (Steps 1-4, Prop5A-C), Lemma 7 for involutions (Lemma7A-C), Prop 8 from Jones (Alt), B.1 orbits for every g (OrbitsB), B.2 (B2), B.3 incl. G = imprimitive D_{2α} ≀ D_{2q} and G ≃* Wreath, |G| (Wreath, WreathA-C), Theorem A, B.4 (Main); kernel G = S_n for (7,4,2), (11,6,4), (13,10,2), (23,18,4), (23,16,6), (23,12,10)
Issue (Claim wording): the Claim's second sentence is literally false without NOTE Thm B's hypotheses: (10,6,3) gcd 1, α ∤ n, imprimitive; (11,8,2) n odd, α ∤ n, 2 orbits; (6,3,2) α ∣ n, |G| = 72 ≠ 144 = |D_6 ≀ D_4|. Lean proves it with them (B.3: n odd; B.4: n odd, gcd = 1), as the T4 brief says ("as NOTE.md states it").
Last reviewer verdict (round 1, audit/r1_*_out.md): literature inputs "VERDICT: FULLY-PROVED-AND-FAITHFUL"; vacuity "VERDICT: FULLY-PROVED-AND-FAITHFUL"; statement "VERDICT: GAPS — 1. statement mismatch — the Claim's second sentence omits necessary hypotheses … Fix: Amend the Claim … Keep the Lean hypotheses. With this correction, the supplied theorem is faithful."
Next step: round-2 statement reviewer against the proposed amended Claim (audit/r2_stmt_prompt.txt); then BLOCKED-MATH pending Adrian's approval of the Claim wording.
Build: bash Riemann/Deep/LeanPrefix/tools/verify.sh (one lake env lean at a time via .scratch/lean4/leanrun.sh; 32 modules rc 0, peak 1.9 GB, ~8 min)
