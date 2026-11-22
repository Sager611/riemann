REVIEW
UTC: 2026-11-22 19:05
Sorry count: 0 (rg 'sorry|admit|native_decide|axiom|implemented_by|extern|unsafe' on *.lean: only the word "axioms" in a Lit.lean docstring)
Files: Riemann/Deep/LeanLi/*.lean (54 modules incl. LR/*), entry Main.lean; PLAN.md
Main theorem: LeanLi.li_remainder_laws : LiteratureInputs → LiRemainderLaws (Main.lean; Statement.lean Part0–Part5); #print axioms: propext, Classical.choice, Quot.sound. Also LeanLi.part1_lower_rh_only : LitCore → RH → (Claim 1 lower bound), no MontgomeryPC.
Build: lake build Riemann.Deep.LeanLi.Main — Build completed successfully (3714 jobs), peak 1636 MB
Literature inputs (Lit.lean; LitCore = all but Montgomery): ζ ≠ 0 on (0,1) (Nakamura 1405.1504; Apostol 12.21); θ_RS Stirling O(1/t) (2004.00926; Edwards §6.5); DLMF 5.11.2 digamma O(|z|^-2); S = O(log T) (Titchmarsh 9.3+9.4; BPT 2009.13791); S₁ = O(log T) (Titchmarsh 9.9); Selberg ∫S² ≪ T log log T (Selberg 1946 Thms 6–7); Voros (TS) S̄_n = t(n)+3/4+O(1/n) (math/0506326); Montgomery 1973 pair correlation (RH; α ∈ [0,1/2]; Goldston math/0412313).
Proved (no input): Li's formula under RH, Lemma 0 (¬RH ⇒ |λ_n| ≥ e^{cn} i.o.), Lemma A with 7/4, Theorems 1, 2, 2+Montgomery, 4, Corollaries.
NOT-FORMALIZED: Claim 5's "Numerically a = 2 first fails at n = 5083" (certified λ_n for all n ≤ 5083 to ~0.09 absolute accuracy: no certified numerics for λ_n in Mathlib, native_decide forbidden).
Last reviewer verdict: round 1 (statement only): GAPS — (1) 5083 not encoded → now declared NOT-FORMALIZED in Statement.lean; (2) Part 4 integer scales → now real N ≥ 2; (3) Montgomery shared by all clauses → LitCore split, part1_lower_rh_only.
Next step: reviewer round 2 (faithfulness + vacuity; literature inputs), fix gaps, then page Lean line, DONE.

## History
- 2026-11-22 11:40 BLOCKED-MEMCAP under memcap.sh 290 (RSS): no Mathlib import fits. Lifted by Adrian 12:08.
- 2026-11-22 15:15 Theorem 1 proved (thm1_li), axioms standard.
- 2026-11-22 17:27 li_remainder_laws assembled, 0 sorry, axioms standard.
- 2026-11-22 19:00 Lemma 0 proved (LiExpGrowth no longer an input); LitCore split; Claim 4 at real scales.

