import Riemann.Fayed.StripReflect

/-!
# Round 7 — transfer from proven RH analogs

| id | analog | status |
|----|--------|--------|
| R7.1 | Ihara zeta, Ramanujan graphs | DEAD: no graph/building with zeta = ζ (Ihara zetas rational) |
| R7.2 | interlacing families (MSS) | DEAD: no finite random-matrix model with E[char poly] → Ξ |
| R7.3 | Selberg zeta, compact surfaces | DEAD: no compact quotient; ζ zeros only resonances |
| R7.4 | Goss zeta via Newton polygons | DEAD: no archimedean valuation forcing |
| R7.5 | period polynomials (Cohn criterion) | DEAD: no certified self-inversive approximants of ξ |
| R7.6 | Pólya / de Bruijn cosine kernels | DEAD: Φ not PF∞; limits need uniform control (real_of_approx_real) |
| R7.7 | Stepanov–Bombieri auxiliary polynomials | DEAD: no degree counting / Frobenius over ℤ (aux_poly_vanishing) |
| R7.8 | fractal uncertainty principle | OPEN, weaker goal (QuasiRH δ strip); mechanism absent |
| R7.9 | Lee–Yang / primon gas | DEAD: no ferromagnetic model with partition function ξ |
| R7.10 | Deligne: monodromy, families | DEAD: equidistribution is statistical, not per member |
-/

open scoped ComplexConjugate

namespace Fayed.R7

/-- Quasi-RH: no zeros with Re s > 1 - δ (s ≠ 1). Open for every δ > 0; δ = 1/2 is RH. -/
def QuasiRH (δ : ℝ) : Prop := ∀ s : ℂ, 1 - δ < s.re → s ≠ 1 → riemannZeta s ≠ 0

/-- RH ⟺ QuasiRH (1/2), using Mathlib's functional equation (strip_holds, reflect_holds). -/
theorem RH_iff_QuasiRH_half : RiemannHypothesis ↔ QuasiRH (1 / 2) := by
  constructor
  · intro h s hs hs1 hz
    rcases le_or_gt 1 s.re with h1 | h1
    · exact riemannZeta_ne_zero_of_one_le_re h1 hz
    · have hnt : ¬∃ n : ℕ, s = -2 * (n + 1) := by
        rintro ⟨n, rfl⟩
        have : (-2 * ((n : ℂ) + 1)).re = -2 * ((n : ℝ) + 1) := by simp
        rw [this] at hs
        have : (0 : ℝ) ≤ n := Nat.cast_nonneg n
        linarith
      have := h s hz hnt hs1
      linarith
  · intro h s hz hnt hs1
    obtain ⟨h0, h1⟩ := Fayed.strip_holds s hz hnt hs1
    rcases lt_trichotomy s.re (1 / 2) with hlt | heq | hgt
    · exfalso
      have hz' := Fayed.reflect_holds s hz h0 h1
      have hre : (1 - conj s).re = 1 - s.re := by simp
      refine h (1 - conj s) (by rw [hre]; linarith) ?_ hz'
      intro he
      have := congrArg Complex.re he
      rw [hre] at this
      simp at this
      linarith
    · exact heq
    · exact absurd hz (h s (by linarith) hs1)

/-- R7.7: Stepanov's principle at its core — a polynomial divisible by (X - a)^n with degree
below n vanishes. Over ℤ there is no degree count to make this bite. -/
theorem aux_poly_vanishing {K : Type*} [Field K] {P : Polynomial K} {a : K} {n : ℕ}
    (hdiv : (Polynomial.X - Polynomial.C a) ^ n ∣ P) (hdeg : P.natDegree < n) : P = 0 := by
  by_contra hP
  have := Polynomial.natDegree_le_of_dvd hdiv hP
  rw [Polynomial.natDegree_pow, Polynomial.natDegree_X_sub_C, mul_one] at this
  omega

/-- R7.6: Hurwitz/Pólya limit principle — a limit of real roots is real. The hard part
(uniform approximation of Ξ by real-rooted functions) is what is missing. -/
theorem real_of_approx_real {w : ℂ}
    (h : ∀ n : ℕ, ∃ z : ℂ, z.im = 0 ∧ ‖z - w‖ < 1 / ((n : ℝ) + 1)) : w.im = 0 := by
  by_contra hw
  obtain ⟨n, hn⟩ := exists_nat_one_div_lt (abs_pos.2 hw)
  obtain ⟨z, hz, hzw⟩ := h n
  have h1 : |(z - w).im| ≤ ‖z - w‖ := Complex.abs_im_le_norm _
  have h2 : (z - w).im = -w.im := by simp [hz]
  rw [h2, abs_neg] at h1
  linarith

end Fayed.R7
