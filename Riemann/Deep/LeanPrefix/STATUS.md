REVIEW
UTC: 2026-11-22 13:41
Sorry count: 0 (rg 'sorry|admit|native_decide|axiom|implemented_by|extern|unsafe': only the #print axioms lines of Check.lean)
Files: Riemann/Deep/LeanPrefix/*.lean (32 modules; Main.lean = statement), tools/verify.sh (clean rebuild + axioms + scan), PLAN.md
Main theorem: Riemann.Deep.LeanPrefix.prefix_reversals (Main.lean); #print axioms: [propext, Classical.choice, Quot.sound]
Literature inputs (hypotheses, Inputs.lean): LZZW2026_Thm3_1 = Lei-Zhao-Zhu-Wu, arXiv 2610.30892v1 (unrefereed), Thm 3.1, first sentence (TeX l. 221-223); Jones2014_Thm1_2_k0 = G. A. Jones, Bull. Aust. Math. Soc. 89 (2014) 159-165, arXiv 1209.5169, Thm 1.2 case k = 0 (TeX l. 54-72), conclusions weakened (lower bounds C_p ≤ G, PGL_d(q) ≤ G and the identification in (c) dropped)
NOT-FORMALIZED: none
Proved (ours, no inputs): Prop 5 (Steps 1-4, Prop5A-C), Lemma 7 for involutions (Lemma7A-C), Prop 8 from Jones (Alt), B.1 orbits for every g (OrbitsB), B.2 (B2), B.3 incl. G = imprimitive D_{2α} ≀ D_{2q} and G ≃* Wreath, |G| (Wreath, WreathA-C), Theorem A, B.4 (Main); kernel G = S_n for (7,4,2), (11,6,4), (13,10,2), (23,18,4), (23,16,6), (23,12,10)
Reading of the Claim (deviation to review): "when α ∣ n … otherwise a primitive group containing A_n" carries NOTE Thm B's hypotheses (B.3: n odd; B.4: n odd, gcd = 1), as the T4 brief says ("as NOTE.md states it"); literally false without them: (12,8,3) α ∣ n, |G| = 768 ≠ 3072 = |D_8 ≀ D_6|; (10,6,3) gcd = 1, α ∤ n, imprimitive; (11,8,2) gcd = 3, intransitive (sympy). B.2 (n even) added so every line triple is covered.
Last reviewer verdict: none yet
Next step: reviewer rounds (statement faithfulness; literature inputs; vacuity).
Build: bash Riemann/Deep/LeanPrefix/tools/verify.sh (one lake env lean at a time via .scratch/lean4/leanrun.sh; 32 modules rc 0, peak 1.9 GB, ~8 min)
