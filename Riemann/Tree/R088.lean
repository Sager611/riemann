/-!
# Round 88 — GRH-side targets (Dirichlet L-functions)

| item | known | Lean cost |
|------|-------|-----------|
| first zeros of L(s, χ), q ≤ 10 | Rumely 1993; Platt 2013 (GRH verified to height 10⁸/q, q ≤ 400,000) | high (beta zero: ultraworker running) |
| Hardy's theorem for L(s, χ) | Hardy–Littlewood 1921; Davenport–Heilbronn 1936 | high (reuse our Hardy files) |
| zero-free regions / Siegel zeros, small q | de la Vallée Poussin; Siegel 1935 (ineffective) | high |
| GRH-verified heights | Platt | high |
| Siegel's theorem (statement) | ineffective constant | high |
| L(1, χ₋₄) = π/4, L(1, χ₋₃) = π/(3√3) | classical | check Mathlib (Leibniz series) |
| simplicity of first zeros | open in general | high |
| ζ_{ℚ(i)} = ζ · L(·, χ₋₄) zeros + Weil positivity | — | high |
-/
