import Riemann.Deep.PNTError.LogDeriv

/-!
# Power-saving PNT, part 3: ψ(x) = x + O(x^{σ₀+ε}) from a zero-free half-plane

`psi_error_of_zeroFree`: if ζ has no zero in Re s > σ₀ (1/2 ≤ σ₀ < 1), then for every ε > 0
there is C with |ψ(x) - x| ≤ C x^{σ₀+ε} for all x ≥ 2.

Proof. For x ≥ 2 at distance ≥ 1/2 from the integers, take T = x⁵ and σ₁ = σ₀ + ε/2:
ψ(x) = (1/2π)∫_{-T}^{T} F(4+iy) dy + O(1) (`perron_psi`); move the line to Re s = σ₁, picking up
the residue x at s = 1 (`shift_psi`); the horizontal sides are O(log x / x) and the new vertical
line is O(x^{σ₁} log² x) by `logDeriv_bound` (|t| ≥ 3) and compactness (|t| ≤ 3). General x: ψ is
constant on [N, N+1), so compare with N + 1/2.
-/

open Complex Set MeasureTheory intervalIntegral Filter Topology
open scoped Real Interval

namespace Fayed.PNTError

noncomputable section

lemma norm_Fψ {x : ℝ} (hx : 0 < x) (s : ℂ) :
    ‖Fψ x s‖ = ‖deriv riemannZeta s / riemannZeta s‖ * x ^ s.re / ‖s‖ := by
  rw [Fψ, norm_div, norm_mul, neg_div, norm_neg, Complex.norm_cpow_eq_rpow_re_of_pos hx]

/-- The horizontal sides of the contour. -/
lemma horiz_psi {σ₁ x T C : ℝ} (hC : ∀ s : ℂ, σ₁ ≤ s.re → 3 ≤ |s.im| →
      ‖deriv riemannZeta s / riemannZeta s‖ ≤ C * Real.log |s.im|) (hC0 : 0 ≤ C)
    (hσ₁0 : 0 ≤ σ₁) (hσ₁ : σ₁ < 4) (hx : 1 ≤ x) (hT : 3 ≤ T) {y : ℝ} (hy : |y| = T) :
    ‖∫ t in σ₁..4, Fψ x (t + y * I)‖ ≤ 4 * (C * Real.log T * x ^ 4 / T) := by
  have hx0 : 0 < x := by linarith
  have hT0 : 0 < T := by linarith
  have hlog : 0 ≤ Real.log T := Real.log_nonneg (by linarith)
  have hK0 : 0 ≤ C * Real.log T * x ^ 4 / T :=
    div_nonneg (mul_nonneg (mul_nonneg hC0 hlog) (by positivity)) hT0.le
  have hpt : ∀ t ∈ Ι σ₁ 4, ‖Fψ x (t + y * I)‖ ≤ C * Real.log T * x ^ 4 / T := by
    intro t ht
    rw [uIoc_of_le hσ₁.le] at ht
    obtain ⟨ht1, ht2⟩ := ht
    have hre : ((t : ℂ) + y * I).re = t := by simp
    have him : ((t : ℂ) + y * I).im = y := by simp
    rw [norm_Fψ hx0, hre]
    have h1 := hC ((t : ℂ) + y * I) (by rw [hre]; linarith) (by rw [him, hy]; exact hT)
    rw [him, hy] at h1
    have h2 : x ^ t ≤ x ^ (4 : ℕ) := by
      have := Real.rpow_le_rpow_of_exponent_le hx ht2
      rwa [show (4 : ℝ) = ((4 : ℕ) : ℝ) by norm_num, Real.rpow_natCast] at this
    have h3 : T ≤ ‖(t : ℂ) + y * I‖ := by
      rw [← hy]
      have := Complex.abs_im_le_norm ((t : ℂ) + y * I)
      rwa [him] at this
    apply div_le_div₀ (mul_nonneg (mul_nonneg hC0 hlog) (by positivity)) _ hT0 h3
    exact mul_le_mul h1 h2 (Real.rpow_nonneg hx0.le t) (mul_nonneg hC0 hlog)
  refine (intervalIntegral.norm_integral_le_of_norm_le_const hpt).trans ?_
  have hσ₁' : |(4 : ℝ) - σ₁| ≤ 4 := by
    rw [abs_le]
    constructor <;> linarith
  calc C * Real.log T * x ^ 4 / T * |4 - σ₁| ≤ C * Real.log T * x ^ 4 / T * 4 :=
        mul_le_mul_of_nonneg_left hσ₁' hK0
    _ = 4 * (C * Real.log T * x ^ 4 / T) := by ring

/-- The vertical line Re s = σ₁. -/
lemma vert_psi {σ₁ x T C B : ℝ} (hC : ∀ s : ℂ, σ₁ ≤ s.re → 3 ≤ |s.im| →
      ‖deriv riemannZeta s / riemannZeta s‖ ≤ C * Real.log |s.im|) (hC0 : 0 ≤ C)
    (hB : ∀ y : ℝ, |y| ≤ 3 →
      ‖deriv riemannZeta (σ₁ + y * I) / riemannZeta (σ₁ + y * I)‖ ≤ B)
    (hB0 : 0 ≤ B) (hσ₁0 : 0 < σ₁) (hx : 1 ≤ x) (hT : 3 ≤ T) :
    ‖∫ y in (-T)..T, Fψ x ((σ₁ : ℂ) + y * I)‖ ≤
      x ^ σ₁ * (4 * B / σ₁ + 2 * C * Real.log T) * (2 * Real.log (1 + T)) := by
  have hx0 : 0 < x := by linarith
  have hT0 : 0 < T := by linarith
  have hlogT : 0 ≤ Real.log T := Real.log_nonneg (by linarith)
  have hxs : 0 < x ^ σ₁ := Real.rpow_pos_of_pos hx0 σ₁
  have hK1 : 4 * B / σ₁ ≤ 4 * B / σ₁ + 2 * C * Real.log T := by
    have : 0 ≤ 2 * C * Real.log T := mul_nonneg (mul_nonneg (by norm_num) hC0) hlogT
    linarith
  have hK2 : 2 * C * Real.log T ≤ 4 * B / σ₁ + 2 * C * Real.log T := by
    have : 0 ≤ 4 * B / σ₁ := div_nonneg (mul_nonneg (by norm_num) hB0) hσ₁0.le
    linarith
  have hpt : ∀ y ∈ Ioc (-T) T, ‖Fψ x ((σ₁ : ℂ) + y * I)‖ ≤
      x ^ σ₁ * (4 * B / σ₁ + 2 * C * Real.log T) * (1 + |y|)⁻¹ := by
    intro y hy
    have hre : ((σ₁ : ℂ) + y * I).re = σ₁ := by simp
    have him : ((σ₁ : ℂ) + y * I).im = y := by simp
    rw [norm_Fψ hx0, hre]
    have hn1 : σ₁ ≤ ‖(σ₁ : ℂ) + y * I‖ := by
      have := Complex.abs_re_le_norm ((σ₁ : ℂ) + y * I)
      rwa [hre, abs_of_pos hσ₁0] at this
    have hn2 : |y| ≤ ‖(σ₁ : ℂ) + y * I‖ := by
      have := Complex.abs_im_le_norm ((σ₁ : ℂ) + y * I)
      rwa [him] at this
    have hy1 : 0 < 1 + |y| := by positivity
    by_cases h3 : 3 ≤ |y|
    · have h1 := hC ((σ₁ : ℂ) + y * I) (by rw [hre]) (by rw [him]; exact h3)
      rw [him] at h1
      have hyT : |y| ≤ T := abs_le.mpr ⟨by linarith [hy.1], hy.2⟩
      have hlog : Real.log |y| ≤ Real.log T := Real.log_le_log (by linarith) hyT
      have h1' : ‖deriv riemannZeta ((σ₁ : ℂ) + y * I) / riemannZeta ((σ₁ : ℂ) + y * I)‖ ≤
          C * Real.log T := h1.trans (mul_le_mul_of_nonneg_left hlog hC0)
      have hs2 : (1 + |y|) / 2 ≤ ‖(σ₁ : ℂ) + y * I‖ := by linarith
      calc ‖deriv riemannZeta ((σ₁ : ℂ) + y * I) / riemannZeta ((σ₁ : ℂ) + y * I)‖ * x ^ σ₁ /
            ‖(σ₁ : ℂ) + y * I‖
          ≤ C * Real.log T * x ^ σ₁ / ((1 + |y|) / 2) :=
            div_le_div₀ (mul_nonneg (mul_nonneg hC0 hlogT) hxs.le)
              (mul_le_mul_of_nonneg_right h1' hxs.le) (by positivity) hs2
        _ = x ^ σ₁ * (2 * C * Real.log T) * (1 + |y|)⁻¹ := by
            field_simp
        _ ≤ x ^ σ₁ * (4 * B / σ₁ + 2 * C * Real.log T) * (1 + |y|)⁻¹ :=
            mul_le_mul_of_nonneg_right (mul_le_mul_of_nonneg_left hK2 hxs.le)
              (inv_nonneg.mpr hy1.le)
    · rw [not_le] at h3
      have h1 := hB y h3.le
      calc ‖deriv riemannZeta ((σ₁ : ℂ) + y * I) / riemannZeta ((σ₁ : ℂ) + y * I)‖ * x ^ σ₁ /
            ‖(σ₁ : ℂ) + y * I‖
          ≤ B * x ^ σ₁ / σ₁ :=
            div_le_div₀ (mul_nonneg hB0 hxs.le) (mul_le_mul_of_nonneg_right h1 hxs.le) hσ₁0 hn1
        _ = x ^ σ₁ * (4 * B / σ₁) * (1 / 4) := by ring
        _ ≤ x ^ σ₁ * (4 * B / σ₁ + 2 * C * Real.log T) * (1 + |y|)⁻¹ := by
            apply mul_le_mul (mul_le_mul_of_nonneg_left hK1 hxs.le) _ (by norm_num)
              (mul_nonneg hxs.le (le_trans (div_nonneg (mul_nonneg (by norm_num) hB0)
                hσ₁0.le) hK1))
            rw [one_div]
            exact inv_anti₀ hy1 (by linarith)
  have hint : IntervalIntegrable
      (fun y : ℝ => x ^ σ₁ * (4 * B / σ₁ + 2 * C * Real.log T) * (1 + |y|)⁻¹) volume (-T) T :=
    (continuous_const.mul (Continuous.inv₀ (by fun_prop) (fun y => by positivity))).intervalIntegrable
      _ _
  refine (intervalIntegral.norm_integral_le_of_norm_le (by linarith) (ae_of_all _ hpt) hint).trans
    (le_of_eq ?_)
  rw [intervalIntegral.integral_const_mul, integral_inv_one_add_abs hT0.le]

/-- ζ'/ζ is bounded on the segment Re s = σ₁, |Im s| ≤ 3. -/
lemma vert_compact {σ₀ σ₁ : ℝ} (hZ : ZeroFree σ₀) (h01 : σ₀ < σ₁) (hσ₁ : σ₁ < 1) :
    ∃ B : ℝ, 0 ≤ B ∧ ∀ y : ℝ, |y| ≤ 3 →
      ‖deriv riemannZeta (σ₁ + y * I) / riemannZeta (σ₁ + y * I)‖ ≤ B := by
  have hK : IsCompact ((fun y : ℝ => (σ₁ : ℂ) + y * I) '' Icc (-3) 3) :=
    isCompact_Icc.image (by fun_prop)
  obtain ⟨B, hB⟩ := logDeriv_bdd_compact hZ hK (by
    rintro s ⟨y, -, rfl⟩
    refine ⟨by simpa using h01, fun h => ?_⟩
    have := congrArg Complex.re h
    simp at this
    linarith)
  refine ⟨max B 0, le_max_right _ _, fun y hy => ?_⟩
  exact (hB _ ⟨y, abs_le.mp hy, rfl⟩).trans (le_max_left _ _)


/-- Combination of Perron, the contour shift and the three side bounds, for general T ≥ 3. -/
lemma psi_bound_aux {σ₀ σ₁ x T C B : ℝ} (hZ : ZeroFree σ₀) (h01 : σ₀ < σ₁) (hσ₁0 : 0 < σ₁)
    (hσ₁ : σ₁ < 1) (hC : ∀ s : ℂ, σ₁ ≤ s.re → 3 ≤ |s.im| →
      ‖deriv riemannZeta s / riemannZeta s‖ ≤ C * Real.log |s.im|) (hC0 : 0 ≤ C)
    (hB : ∀ y : ℝ, |y| ≤ 3 →
      ‖deriv riemannZeta (σ₁ + y * I) / riemannZeta (σ₁ + y * I)‖ ≤ B)
    (hB0 : 0 ≤ B) (hx : 1 ≤ x) (hT : 3 ≤ T) (hfar : ∀ n : ℕ, 1 / 2 ≤ |x - n|) :
    |Chebyshev.psi x - x| ≤ 4 * x ^ 4 * (x + 1) / (π * T) * (π ^ 2 / 6) +
      1 / (2 * π) * (4 * (C * Real.log T * x ^ 4 / T) + 4 * (C * Real.log T * x ^ 4 / T)) +
      1 / (2 * π) * (x ^ σ₁ * (4 * B / σ₁ + 2 * C * Real.log T) * (2 * Real.log (1 + T))) := by
  have hx0 : 0 < x := by linarith
  have hT0 : 0 < T := by linarith
  have hP := perron_psi hx hT0 hfar
  have hS := shift_psi hZ h01 hσ₁0 hσ₁ hx0 hT0
  have hHb := horiz_psi hC hC0 hσ₁0.le (by linarith) hx hT (y := -T)
    (by rw [abs_neg, abs_of_pos hT0])
  have hHt := horiz_psi hC hC0 hσ₁0.le (by linarith) hx hT (y := T) (abs_of_pos hT0)
  have hV := vert_psi hC hC0 hB hB0 hσ₁0 hx hT
  have hc : ‖((1 / (2 * π) : ℝ) : ℂ)‖ = 1 / (2 * π) := by
    rw [Complex.norm_real, Real.norm_eq_abs, abs_of_pos (by positivity)]
  have e : ((Chebyshev.psi x - x : ℝ) : ℂ) =
      ((Chebyshev.psi x : ℂ) - ((1 / (2 * π) : ℝ) : ℂ) * ∫ y in (-T)..T, Fψ x ((4 : ℝ) + y * I)) +
      ((1 / (2 * π) : ℝ) : ℂ) * (I * (∫ t in σ₁..4, Fψ x (t + (-T : ℝ) * I)) -
        I * (∫ t in σ₁..4, Fψ x (t + (T : ℝ) * I))) +
      ((1 / (2 * π) : ℝ) : ℂ) * ∫ y in (-T)..T, Fψ x ((σ₁ : ℂ) + y * I) := by
    rw [← hS]
    push_cast
    ring
  have hnorm : |Chebyshev.psi x - x| = ‖((Chebyshev.psi x - x : ℝ) : ℂ)‖ := by
    rw [Complex.norm_real, Real.norm_eq_abs]
  rw [hnorm, e]
  refine (norm_add_le _ _).trans (add_le_add ((norm_add_le _ _).trans (add_le_add hP ?_)) ?_)
  · rw [norm_mul, hc]
    refine mul_le_mul_of_nonneg_left ((norm_sub_le _ _).trans ?_) (by positivity)
    simp only [norm_mul, Complex.norm_I, one_mul]
    exact add_le_add hHb hHt
  · rw [norm_mul, hc]
    exact mul_le_mul_of_nonneg_left hV (by positivity)

lemma num_b1 {x : ℝ} (hx : 2 ≤ x) : 4 * x ^ 4 * (x + 1) / (π * x ^ 5) * (π ^ 2 / 6) ≤ 2 * π := by
  have hx0 : 0 < x := by linarith
  have hπ : 0 < π := Real.pi_pos
  have e1 : 4 * x ^ 4 * (x + 1) / (π * x ^ 5) * (π ^ 2 / 6) = 2 * π / 3 * ((x + 1) / x) := by
    field_simp
    ring
  have hq : (x + 1) / x ≤ 3 / 2 := by
    rw [div_le_iff₀ hx0]
    linarith
  rw [e1]
  have : 2 * π / 3 * ((x + 1) / x) ≤ 2 * π / 3 * (3 / 2) :=
    mul_le_mul_of_nonneg_left hq (by positivity)
  linarith

lemma num_b2 {C x : ℝ} (hC : 0 ≤ C) (hx : 1 ≤ x) :
    1 / (2 * π) * (4 * (C * Real.log (x ^ 5) * x ^ 4 / x ^ 5) +
      4 * (C * Real.log (x ^ 5) * x ^ 4 / x ^ 5)) ≤ 20 * C / π := by
  have hx0 : 0 < x := by linarith
  have hπ : 0 < π := Real.pi_pos
  have hLx : Real.log x ≤ x := Real.log_le_self hx0.le
  have e2 : 1 / (2 * π) * (4 * (C * Real.log (x ^ 5) * x ^ 4 / x ^ 5) +
      4 * (C * Real.log (x ^ 5) * x ^ 4 / x ^ 5)) = 20 * C / π * (Real.log x / x) := by
    rw [Real.log_pow]
    field_simp
    ring
  have hq : Real.log x / x ≤ 1 := by
    rw [div_le_iff₀ hx0]
    linarith
  rw [e2]
  have h0 : 0 ≤ 20 * C / π := div_nonneg (by linarith) hπ.le
  calc 20 * C / π * (Real.log x / x) ≤ 20 * C / π * 1 := mul_le_mul_of_nonneg_left hq h0
    _ = 20 * C / π := mul_one _

lemma num_step2 {a b L : ℝ} (ha : 0 ≤ a) (hb : 0 ≤ b) (hL : 0 ≤ L) :
    (a + 10 * b * L) * L ≤ (a + 10 * b) * (1 + L) ^ 2 := by
  nlinarith [mul_nonneg ha hL, mul_nonneg (mul_nonneg ha hL) hL, mul_nonneg hb hL]

lemma num_b3 {σ₁ B C x : ℝ} (hσ₁0 : 0 < σ₁) (hB : 0 ≤ B) (hC : 0 ≤ C) (hx : 2 ≤ x) :
    1 / (2 * π) * (x ^ σ₁ * (4 * B / σ₁ + 2 * C * Real.log (x ^ 5)) *
      (2 * Real.log (1 + x ^ 5))) ≤
      6 / π * (4 * B / σ₁ + 10 * C) * x ^ σ₁ * (1 + Real.log x) ^ 2 := by
  have hx0 : 0 < x := by linarith
  have hπ : 0 < π := Real.pi_pos
  have hxs : 0 < x ^ σ₁ := Real.rpow_pos_of_pos hx0 σ₁
  have hL2 : Real.log 2 ≤ Real.log x := Real.log_le_log (by norm_num) hx
  have hlog2 : 0 < Real.log 2 := Real.log_pos (by norm_num)
  have hL0 : 0 ≤ Real.log x := by linarith
  have hT0 : 0 < x ^ 5 := by positivity
  have hlog1T : Real.log (1 + x ^ 5) ≤ 6 * Real.log x := by
    have h1 : (1 : ℝ) ≤ x ^ 5 := one_le_pow₀ (by linarith)
    calc Real.log (1 + x ^ 5) ≤ Real.log (2 * x ^ 5) := Real.log_le_log (by linarith) (by linarith)
      _ = Real.log 2 + 5 * Real.log x := by
          rw [Real.log_mul (by norm_num) hT0.ne', Real.log_pow]
          norm_num
      _ ≤ 6 * Real.log x := by linarith
  have hB' : 0 ≤ 4 * B / σ₁ := div_nonneg (by linarith) hσ₁0.le
  have hK : 4 * B / σ₁ + 2 * C * Real.log (x ^ 5) = 4 * B / σ₁ + 10 * C * Real.log x := by
    rw [Real.log_pow]
    ring
  have hK0 : 0 ≤ 4 * B / σ₁ + 10 * C * Real.log x := add_nonneg hB' (by positivity)
  have step1 : x ^ σ₁ * (4 * B / σ₁ + 2 * C * Real.log (x ^ 5)) * (2 * Real.log (1 + x ^ 5)) ≤
      x ^ σ₁ * (4 * B / σ₁ + 10 * C * Real.log x) * (12 * Real.log x) := by
    rw [hK]
    apply mul_le_mul_of_nonneg_left _ (mul_nonneg hxs.le hK0)
    linarith
  have step2 := num_step2 hB' hC hL0
  calc 1 / (2 * π) * (x ^ σ₁ * (4 * B / σ₁ + 2 * C * Real.log (x ^ 5)) *
        (2 * Real.log (1 + x ^ 5)))
      ≤ 1 / (2 * π) * (x ^ σ₁ * (4 * B / σ₁ + 10 * C * Real.log x) * (12 * Real.log x)) :=
        mul_le_mul_of_nonneg_left step1 (by positivity)
    _ = 6 / π * x ^ σ₁ * ((4 * B / σ₁ + 10 * C * Real.log x) * Real.log x) := by
        field_simp
        ring
    _ ≤ 6 / π * x ^ σ₁ * ((4 * B / σ₁ + 10 * C) * (1 + Real.log x) ^ 2) :=
        mul_le_mul_of_nonneg_left step2 (mul_nonneg (by positivity) hxs.le)
    _ = 6 / π * (4 * B / σ₁ + 10 * C) * x ^ σ₁ * (1 + Real.log x) ^ 2 := by ring

/-- The estimate at points at distance ≥ 1/2 from the integers:
|ψ(x) - x| ≤ A + B x^{σ₁} (1 + log x)². -/
theorem psi_halfint {σ₀ : ℝ} (hZ : ZeroFree σ₀) (h12 : 1 / 2 ≤ σ₀) {σ₁ : ℝ} (h01 : σ₀ < σ₁)
    (hσ₁ : σ₁ < 1) :
    ∃ A B : ℝ, 0 ≤ A ∧ 0 ≤ B ∧ ∀ x : ℝ, 2 ≤ x → (∀ n : ℕ, 1 / 2 ≤ |x - n|) →
      |Chebyshev.psi x - x| ≤ A + B * x ^ σ₁ * (1 + Real.log x) ^ 2 := by
  have h1 : σ₀ < 1 := lt_trans h01 hσ₁
  have hσ₁0 : 0 < σ₁ := by linarith
  have hπ : 0 < π := Real.pi_pos
  obtain ⟨C, hC0, hC⟩ := logDeriv_bound hZ h12 h1 (sub_pos.mpr h01)
  have hC' : ∀ s : ℂ, σ₁ ≤ s.re → 3 ≤ |s.im| →
      ‖deriv riemannZeta s / riemannZeta s‖ ≤ C * Real.log |s.im| :=
    fun s hs ht => hC s (by linarith) ht
  obtain ⟨B₀, hB₀0, hB₀⟩ := vert_compact hZ h01 hσ₁
  have hA0 : 0 ≤ 2 * π + 20 * C / π := add_nonneg (by positivity) (div_nonneg (by linarith) hπ.le)
  have hD0 : 0 ≤ 4 * B₀ / σ₁ + 10 * C :=
    add_nonneg (div_nonneg (by linarith) hσ₁0.le) (by linarith)
  refine ⟨2 * π + 20 * C / π, 6 / π * (4 * B₀ / σ₁ + 10 * C), hA0,
    mul_nonneg (div_nonneg (by norm_num) hπ.le) hD0, ?_⟩
  intro x hx hfar
  have hx1 : 1 ≤ x := by linarith
  have hT3 : 3 ≤ x ^ 5 := by
    have := pow_le_pow_left₀ (by norm_num : (0 : ℝ) ≤ 2) hx 5
    norm_num at this
    linarith
  have hm := psi_bound_aux hZ h01 hσ₁0 hσ₁ hC' hC0 hB₀ hB₀0 hx1 hT3 hfar
  have b1 := num_b1 hx
  have b2 := num_b2 hC0 hx1
  have b3 := num_b3 hσ₁0 hB₀0 hC0 hx
  linarith only [hm, b1, b2, b3]

/-- **Rung 2 (general form).** If ζ has no zero in Re s > σ₀ with 1/2 ≤ σ₀ < 1, then for every
ε > 0 there is C with |ψ(x) - x| ≤ C x^{σ₀+ε} for all x ≥ 2. -/
theorem psi_error_of_zeroFree {σ₀ : ℝ} (hZ : ZeroFree σ₀) (h12 : 1 / 2 ≤ σ₀) (h1 : σ₀ < 1)
    {ε : ℝ} (hε : 0 < ε) :
    ∃ C : ℝ, ∀ x : ℝ, 2 ≤ x → |Chebyshev.psi x - x| ≤ C * x ^ (σ₀ + ε) := by
  by_cases hbig : 1 ≤ σ₀ + ε
  · refine ⟨Real.log 4 + 5, fun x hx => ?_⟩
    have hx0 : 0 ≤ x := by linarith
    have hψ := Chebyshev.psi_le_const_mul_self hx0
    have hψ0 := Chebyshev.psi_nonneg x
    have hxe : x ≤ x ^ (σ₀ + ε) := by
      calc x = x ^ (1 : ℝ) := (Real.rpow_one x).symm
        _ ≤ x ^ (σ₀ + ε) := Real.rpow_le_rpow_of_exponent_le (by linarith) hbig
    have hl4 : 0 ≤ Real.log 4 := Real.log_nonneg (by norm_num)
    rw [abs_le]
    constructor <;> nlinarith
  · rw [not_le] at hbig
    obtain ⟨A, B, hA, hB, hmain⟩ :=
      psi_halfint hZ h12 (σ₁ := σ₀ + ε / 2) (by linarith) (by linarith)
    set D := A + B * (1 + 4 / ε) ^ 2 with hD
    have hD0 : 0 ≤ D := add_nonneg hA (mul_nonneg hB (by positivity))
    -- the bound at half-integers
    have hhalf : ∀ y : ℝ, 2 ≤ y → (∀ n : ℕ, 1 / 2 ≤ |y - n|) →
        |Chebyshev.psi y - y| ≤ D * y ^ (σ₀ + ε) := by
      intro y hy hfar
      have hy0 : 0 < y := by linarith
      have hy1 : 1 ≤ y := by linarith
      have hm := hmain y hy hfar
      have hq : 0 < ε / 4 := by linarith
      have hlog : Real.log y ≤ y ^ (ε / 4) / (ε / 4) := Real.log_le_rpow_div hy0.le hq
      have hy4 : 1 ≤ y ^ (ε / 4) := Real.one_le_rpow hy1 hq.le
      have hlog0 : 0 ≤ Real.log y := Real.log_nonneg hy1
      have hl1 : 1 + Real.log y ≤ (1 + 4 / ε) * y ^ (ε / 4) := by
        have : y ^ (ε / 4) / (ε / 4) = 4 / ε * y ^ (ε / 4) := by field_simp
        rw [this] at hlog
        nlinarith
      have hsq : (1 + Real.log y) ^ 2 ≤ (1 + 4 / ε) ^ 2 * y ^ (ε / 2) := by
        have e2 : y ^ (ε / 2) = y ^ (ε / 4) * y ^ (ε / 4) := by
          rw [← Real.rpow_add hy0]
          ring_nf
        rw [e2]
        have h0 : 0 ≤ 1 + Real.log y := by linarith
        nlinarith [mul_le_mul hl1 hl1 h0 (by positivity)]
      have hys : 0 < y ^ (σ₀ + ε / 2) := Real.rpow_pos_of_pos hy0 _
      have hyy : y ^ (σ₀ + ε / 2) * y ^ (ε / 2) = y ^ (σ₀ + ε) := by
        rw [← Real.rpow_add hy0]
        ring_nf
      have hy1' : 1 ≤ y ^ (σ₀ + ε) := Real.one_le_rpow hy1 (by linarith)
      calc |Chebyshev.psi y - y| ≤ A + B * y ^ (σ₀ + ε / 2) * (1 + Real.log y) ^ 2 := hm
        _ ≤ A * y ^ (σ₀ + ε) + B * y ^ (σ₀ + ε / 2) * ((1 + 4 / ε) ^ 2 * y ^ (ε / 2)) := by
            apply add_le_add
            · nlinarith
            · exact mul_le_mul_of_nonneg_left hsq (mul_nonneg hB hys.le)
        _ = D * y ^ (σ₀ + ε) := by
            rw [hD, ← hyy]
            ring
    refine ⟨2 * D + 1, fun x hx => ?_⟩
    have hx0 : 0 < x := by linarith
    set N := ⌊x⌋₊ with hN
    have hN2 : 2 ≤ N := by
      rw [hN]
      exact Nat.le_floor (by exact_mod_cast hx)
    have hNx : (N : ℝ) ≤ x := Nat.floor_le hx0.le
    have hxN : x < N + 1 := Nat.lt_floor_add_one x
    set y : ℝ := N + 1 / 2 with hy
    have hy2 : 2 ≤ y := by
      have : (2 : ℝ) ≤ N := by exact_mod_cast hN2
      rw [hy]; linarith
    have hfar : ∀ n : ℕ, 1 / 2 ≤ |y - n| := by
      intro n
      rcases le_or_gt n N with h | h
      · have : (n : ℝ) ≤ N := by exact_mod_cast h
        rw [abs_of_pos (by rw [hy]; linarith)]
        rw [hy]; linarith
      · have : (N : ℝ) + 1 ≤ n := by exact_mod_cast h
        rw [abs_of_neg (by rw [hy]; linarith)]
        rw [hy]; linarith
    have hfloor : ⌊y⌋₊ = N := by
      rw [Nat.floor_eq_iff (by linarith)]
      constructor <;> rw [hy] <;> linarith
    have hψ : Chebyshev.psi x = Chebyshev.psi y := by
      rw [Chebyshev.psi_eq_psi_coe_floor x, Chebyshev.psi_eq_psi_coe_floor y, hfloor]
    have hyx : |y - x| ≤ 1 / 2 := by
      rw [abs_le]; constructor <;> rw [hy] <;> linarith
    have hθ0 : 0 < σ₀ + ε := by linarith
    have hθ1 : σ₀ + ε ≤ 1 := hbig.le
    have hy2x : y ≤ 2 * x := by rw [hy]; linarith
    have hpow : y ^ (σ₀ + ε) ≤ 2 * x ^ (σ₀ + ε) := by
      calc y ^ (σ₀ + ε) ≤ (2 * x) ^ (σ₀ + ε) :=
            Real.rpow_le_rpow (by linarith) hy2x hθ0.le
        _ = 2 ^ (σ₀ + ε) * x ^ (σ₀ + ε) := Real.mul_rpow (by norm_num) hx0.le
        _ ≤ 2 * x ^ (σ₀ + ε) := by
            apply mul_le_mul_of_nonneg_right _ (Real.rpow_nonneg hx0.le _)
            calc (2 : ℝ) ^ (σ₀ + ε) ≤ 2 ^ (1 : ℝ) :=
                  Real.rpow_le_rpow_of_exponent_le (by norm_num) hθ1
              _ = 2 := Real.rpow_one 2
    have hx1 : 1 ≤ x ^ (σ₀ + ε) := Real.one_le_rpow (by linarith) hθ0.le
    have hmain' := hhalf y hy2 hfar
    have htri : |Chebyshev.psi x - x| ≤ |Chebyshev.psi y - y| + |y - x| := by
      rw [hψ]
      calc |Chebyshev.psi y - x| = |(Chebyshev.psi y - y) + (y - x)| := by ring_nf
        _ ≤ |Chebyshev.psi y - y| + |y - x| := abs_add_le _ _
    calc |Chebyshev.psi x - x| ≤ D * y ^ (σ₀ + ε) + 1 / 2 := by linarith
      _ ≤ D * (2 * x ^ (σ₀ + ε)) + 1 * x ^ (σ₀ + ε) := by
          apply add_le_add (mul_le_mul_of_nonneg_left hpow hD0)
          linarith
      _ = (2 * D + 1) * x ^ (σ₀ + ε) := by ring

end

end Fayed.PNTError
