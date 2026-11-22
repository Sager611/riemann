WORKING
UTC: 2026-11-22 17:10
Sorry count: 0 in compiled modules (Defs, Lit, Statement, ZerosBasic, Abel, ZerosCount, CountBounds, ZeroSums, LiCoeff, LiRH, Smooth, LR/*)
Files: Riemann/Deep/LeanLi/ (probe, tools); PLAN.md next
Main theorem: (planned) Riemann.Deep.LeanLi.li_remainder_laws (the five parts of discoveries/li-remainder.md)
Literature inputs: (planned) Riemann–von Mangoldt (Titchmarsh Thm 9.4), S₁ = O(log T) (Titchmarsh 9.9), Selberg 1946 (∫S²), Montgomery 1973 (pair correlation, RH), Stirling/digamma asymptotics (DLMF 5.11), Voros math/0506326 (TS), ζ ≠ 0 on (0, 1)
NOT-FORMALIZED: part 5's "a = 2 first fails at n = 5083" (certified λ_n for n ≤ 5083 out of reach in Lean)
Last reviewer verdict: none
Next step: Part0 proved (λ_n converges, real); Li's formula under RH; smooth exponential-sum decay. Now: I(n) closed form (Fourier–log integrals), weights, trend lemma, Theorem 1.

## History
- 2026-11-22 11:40 BLOCKED-MEMCAP under memcap.sh 290 (RSS): no Mathlib import fits; lake alone 620–661 MB RSS (memprobe.txt). Lifted by Adrian 12:08.
