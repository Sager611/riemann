import Mathlib

/-!
# Round 38 — next audit targets: 2025–26 claimed RH proofs

| area | claim | key step | testable |
|------|-------|----------|----------|
| Weil positivity | S. Velez, Zenodo 21173721 (Jul 2026) | positivity of a "completed arithmetic Weil kernel" | yes (kernel eigenvalues) |
| Li coefficients | Zenodo 16900110 (Aug 2025) | phase alignment forces λ_n ≥ 0 | yes (λ_n from zeros) |
| de Bruijn–Newman | A. Gershon, Preprints 202604.1513 (Apr 2026) | strict log-concavity of the kernel + finite Jensen checks | yes → Orca audit |
| Hilbert–Pólya | E. Yakaboylu, arXiv 2408.15135v11 | self-adjoint H with spectrum i(1/2 - ρ) | yes (truncations) |
| Nyman–Beurling | Bhattacharjee–Nandi–Frederick, Preprints 202506.0772 | partition/denseness argument | yes (d_N) |
| Robin | Zenodo 21869761 (Aug 2026) | Robin above a threshold | yes |
| quadratic form | Zenodo 18923112 (Mar 2026) | Q(A) ≪ A^ε | yes |
| argument principle | Y. Abe, Cambridge Open Engage v25 (Jul 2026) | symmetric counting forces Re = 1/2 | yes |

Lean below: algebraic core of a counterexample to "strictly log-concave kernel ⇒ real zeros".
K(t) = e^{-t²}(1 + t⁴/4) is even, positive and strictly log-concave (with x = t²,
(log K)'' < 0 ⇔ x⁴/8 + x³/4 + x² - 3x + 2 > 0). Its cosine transform is
√π e^{-z²/4} (19/16 - 3z²/16 + z⁴/64) (Gaussian moments), whose factor has no real root w = z²,
so all four zeros are non-real.
-/

namespace Fayed.R38

/-- Strict log-concavity certificate for K(t) = e^{-t²}(1 + t⁴/4), with x = t². -/
theorem logconcave_poly (x : ℝ) (hx : 0 ≤ x) : 0 < x ^ 4 / 8 + x ^ 3 / 4 + x ^ 2 - 3 * x + 2 := by
  have hq : 0 ≤ x ^ 2 / 8 + x / 2 + 15 / 8 := by nlinarith [sq_nonneg x]
  nlinarith [mul_nonneg (sq_nonneg (x - 1)) hq]

/-- The transform's polynomial factor in w = z² has no real root. -/
theorem quartic_factor_pos (w : ℝ) : 0 < 19 / 16 - 3 / 16 * w + w ^ 2 / 64 := by
  nlinarith [sq_nonneg (w / 8 - 3 / 4)]

end Fayed.R38

#print axioms Fayed.R38.logconcave_poly
