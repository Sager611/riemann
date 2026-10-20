import Mathlib

/-!
# Round 34 — real-rootedness toolbox applied to Ξ: all dead or circular

| tool | proves | needed for Ξ | verdict |
|------|--------|--------------|---------|
| Lee–Yang / Newman 1976 | ferromagnetic zeros on a circle / line | Lee–Yang measure for Ξ | CIRCULAR |
| Pólya 1927, de Bruijn 1950 | ∫ e^{-t^{2n}} cos zt real-rooted | Pólya kernel hypotheses for Riemann's kernel | DEAD |
| Jensen polynomials (GORZ 2019), Turán (CNV 1986) | each degree hyperbolic for large shifts | all degrees, all shifts | DEAD (asymptotic only) |
| Borcea–Brändén preservers | classification of real-rootedness preservers | Ξ as image of a real-rooted object | DEAD |
| interlacing families (MSS 2015) | one good leaf | Ξ itself controlled | DEAD |
| Hermite–Biehler | A = (E+E#)/2 real-rooted | admissible E with A = Ξ | CIRCULAR |
| Ki–Kim–Lee 2009, Cardon | Λ < 1/2; Gaussian weights preserve real zeros | Λ ≤ 0 | CIRCULAR |
| Jacobi matrices | char. polys real-rooted | Jacobi approximants of Ξ | CIRCULAR |

Lean below: the degree-2 case — a real-rooted quadratic Jensen polynomial forces the Turán inequality.
-/

namespace Fayed.R34

/-- If a z² + 2 b z + c factors over ℝ as a (z - x)(z - y), then a c ≤ b² (Turán inequality). -/
theorem turan_of_real_roots (a b c x y : ℝ)
    (h : ∀ z : ℝ, a * z ^ 2 + 2 * b * z + c = a * (z - x) * (z - y)) : a * c ≤ b ^ 2 := by
  have h0 := h 0
  have h1 := h 1
  have h2 := h (-1)
  have hc : c = a * x * y := by linarith
  have hb : b = -(a * x + a * y) / 2 := by linarith
  subst hc
  subst hb
  nlinarith [sq_nonneg (a * x - a * y)]

end Fayed.R34

#print axioms Fayed.R34.turan_of_real_roots
