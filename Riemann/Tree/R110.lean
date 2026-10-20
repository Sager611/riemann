/-!
# Round 110 — records for certified Weil positivity on windows (scouts; 4 of 8 answered)

| source | window and rigor |
|--------|------------------|
| Chuk 2026, arXiv 2608.24827 | certified numerics; convention: test function in [−L, L] (autocorrelation [−2L, 2L]) |
| Mori | autocorrelation in [−log 15, log 15]: heuristic numerics, not a proof |
| 2026 analytic work | no rigorous extension of Yoshida's analytic window beyond log 2 found |
| Li coefficients | λ_n > 0 certified for n ≤ 10⁶ (Johansson 2023), not a window result |
| formal verification | none found except ours (weilPos_lt_log2); a GitHub note (Kuberwastaken 2026) is heuristic |

So a Lean proof beyond log 2 (yoshida5) would be the first formal one; certified-numerics windows exist (Chuk; ours 1.75 in Arb).
-/

/-! Late scouts: Connes–Consani 2020 (arXiv 2006.13771) window [−1,1] heuristic; Garrett heuristic; Suzuki arXiv 2607.24830 numerical; a claimed certified autocorrelation window [−3,3] by "Philip Gordon 2026, id CC+" is unverified (no usable id). -/
