/-!
# Round 99 — tooling to scale certified zeros to height 100–1000 in Lean

| item | verdict |
|------|---------|
| Riemann–Siegel with Gabcke's explicit error (1979) | LATER (saddle-point remainder is a large formalization) |
| Turing's method for exact N(T) (Turing 1953, Trudgian's constants) | USE (avoids covering the strip) |
| explicit S(T) bounds (Backlund 1918, Trudgian 2014) | LATER |
| Lean interval / Taylor-model libraries | LATER (porting cost high) |
| kernel-checked certificates without native_decide (GMP-backed Nat, staged norm_num) | USE |
| Euler–Maclaurin vs Riemann–Siegel crossover | LATER (EM is fine to height ~100) |
| Gram points, Gram's law, Rosser's rule | USE (cheap separation certificate) |
| Hardy Z(t) sign changes | USE (Λ sign route already in Lean) |

Plan: keep Euler–Maclaurin to t ≈ 100, add Gram-point sign certificates and Turing's method for exact counts.
-/
