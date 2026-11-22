import Riemann.Deep.LeanNBRate.Multiplier
import Riemann.Deep.LeanNBRate.Windows
import Mathlib.Analysis.SpecialFunctions.ImproperIntegrals

/-!
# The multiplier bound in terms of sup-norms of the jets of `G`
-/

open MeasureTheory Filter Topology Set Complex

noncomputable section

namespace Riemann.Deep.LeanNBRate

/-- `∫ (R² + v²)⁻¹ dv = π / R`. -/
theorem integral_inv_sq_add_sq {R : ℝ} (hR : 0 < R) :
    ∫ v : ℝ, (R ^ 2 + v ^ 2)⁻¹ = Real.pi / R := by
  have h := Measure.integral_comp_mul_left (fun u : ℝ => (R ^ 2 * (1 + u ^ 2))⁻¹) R⁻¹
  have e : (fun v : ℝ => (R ^ 2 * (1 + (R⁻¹ * v) ^ 2))⁻¹) = fun v => (R ^ 2 + v ^ 2)⁻¹ := by
    ext v; congr 1; field_simp
  rw [e] at h
  rw [h, inv_inv, abs_of_pos hR, smul_eq_mul]
  have : (fun u : ℝ => (R ^ 2 * (1 + u ^ 2))⁻¹) = fun u => (R ^ 2)⁻¹ * (1 + u ^ 2)⁻¹ := by
    ext u; rw [mul_inv]
  rw [this, integral_const_mul, integral_univ_inv_one_add_sq]
  field_simp

theorem integrable_inv_sq_add_sq {R : ℝ} (hR : 0 < R) :
    Integrable fun v : ℝ => (R ^ 2 + v ^ 2)⁻¹ := by
  have h := (integrable_inv_one_add_sq.comp_mul_left' (inv_ne_zero hR.ne')).const_mul (R ^ 2)⁻¹
  refine h.congr (Eventually.of_forall fun v => ?_)
  simp only
  rw [← mul_inv]; congr 1; field_simp

/-- `∫ |v|^j |kern G v| ≤ R A_j + A_{j+2}/R`, `A_n = ∫ |G^{(n)}|`. -/
theorem IsTest.integral_pow_mul_norm_kern_le {G : ℝ → ℂ} (hG : IsTest G) (j : ℕ) {R : ℝ}
    (hR : 0 < R) :
    ∫ v : ℝ, |v| ^ j * ‖kern G v‖ ≤
      R * (∫ θ : ℝ, ‖iteratedDeriv j G θ‖) + (∫ θ : ℝ, ‖iteratedDeriv (j + 2) G θ‖) / R := by
  set Aj := ∫ θ : ℝ, ‖iteratedDeriv j G θ‖
  set Aj2 := ∫ θ : ℝ, ‖iteratedDeriv (j + 2) G θ‖
  have hπ : 0 < Real.pi := Real.pi_pos
  have hAj : 0 ≤ Aj := integral_nonneg fun _ => norm_nonneg _
  have hAj2 : 0 ≤ Aj2 := integral_nonneg fun _ => norm_nonneg _
  set a := (2 * Real.pi)⁻¹ * Aj
  set b := (2 * Real.pi)⁻¹ * Aj2
  have ha : 0 ≤ a := by positivity
  have hb : 0 ≤ b := by positivity
  have hpt : ∀ v, |v| ^ j * ‖kern G v‖ ≤ 2 * (a + b / R ^ 2) * R ^ 2 * (R ^ 2 + v ^ 2)⁻¹ := by
    intro v
    have h1 : |v| ^ j * ‖kern G v‖ ≤ a := hG.norm_pow_mul_kern_le j v
    have h2 : |v| ^ (j + 2) * ‖kern G v‖ ≤ b := hG.norm_pow_mul_kern_le (j + 2) v
    have hR2 : 0 < R ^ 2 := by positivity
    have hden : 0 < R ^ 2 + v ^ 2 := by positivity
    rw [← div_eq_mul_inv, le_div_iff₀ hden]
    set X := |v| ^ j * ‖kern G v‖
    have hX : 0 ≤ X := by positivity
    have h2' : X * v ^ 2 ≤ b := by
      calc X * v ^ 2 = |v| ^ (j + 2) * ‖kern G v‖ := by simp only [X]; rw [pow_add, ← sq_abs v]; ring
        _ ≤ b := h2
    have : X * (R ^ 2 + v ^ 2) ≤ a * R ^ 2 + b := by nlinarith
    calc X * (R ^ 2 + v ^ 2) ≤ a * R ^ 2 + b := this
      _ ≤ 2 * (a + b / R ^ 2) * R ^ 2 := by
          rw [show 2 * (a + b / R ^ 2) * R ^ 2 = 2 * a * R ^ 2 + 2 * b by field_simp]; nlinarith
  have hint := hG.integrable_norm_pow_mul_kern j
  have hint2 := (integrable_inv_sq_add_sq hR).const_mul (2 * (a + b / R ^ 2) * R ^ 2)
  calc ∫ v : ℝ, |v| ^ j * ‖kern G v‖ ≤ ∫ v : ℝ, 2 * (a + b / R ^ 2) * R ^ 2 * (R ^ 2 + v ^ 2)⁻¹ :=
        integral_mono hint hint2 hpt
    _ = 2 * (a + b / R ^ 2) * R ^ 2 * (Real.pi / R) := by
        rw [integral_const_mul, integral_inv_sq_add_sq hR]
    _ = R * Aj + Aj2 / R := by
        simp only [a, b]; field_simp

/-- `(1 + t)^n ≤ 2^n (1 + t^n)` for `t ≥ 0`. -/
theorem one_add_pow_le (n : ℕ) {t : ℝ} (ht : 0 ≤ t) : (1 + t) ^ n ≤ 2 ^ n * (1 + t ^ n) := by
  rcases le_total t 1 with h | h
  · calc (1 + t) ^ n ≤ 2 ^ n := by gcongr; linarith
      _ ≤ 2 ^ n * (1 + t ^ n) := by
          have : 0 ≤ t ^ n := by positivity
          nlinarith [pow_nonneg (by norm_num : (0:ℝ) ≤ 2) n]
  · calc (1 + t) ^ n ≤ (2 * t) ^ n := by gcongr; linarith
      _ = 2 ^ n * t ^ n := by rw [mul_pow]
      _ ≤ 2 ^ n * (1 + t ^ n) := by gcongr; linarith

/-- **Multiplier bound.** `‖Λ_r(e^{iθx}G)‖ ≤ c^{M−1}·K·2^{M−1}·(R₀A₀ + A₂/R₀ + ℓ^{1−M}(R₁A_{M−1} + A_{M+1}/R₁))`,
`A_n = ∫|G^{(n)}|`. -/
theorem ml_bound (hT : Nazarov1993_Turan) : ∃ c : ℝ, 0 < c ∧
    ∀ (r : ExpSum) {G : ℝ → ℂ} (_ : IsTest G) (ℓ K x R₀ R₁ : ℝ), 0 < ℓ → 1 ≤ r.M →
      x ∈ Icc 0 ℓ → 0 < R₀ → 0 < R₁ → (∀ u ∈ Icc 0 ℓ, ‖r.eval u‖ ≤ K) →
      ‖r.Lam x G‖ ≤ c ^ (r.M - 1) * K * 2 ^ (r.M - 1) *
        (R₀ * (∫ θ : ℝ, ‖G θ‖) + (∫ θ : ℝ, ‖iteratedDeriv 2 G θ‖) / R₀ +
          (ℓ ^ (r.M - 1))⁻¹ * (R₁ * (∫ θ : ℝ, ‖iteratedDeriv (r.M - 1) G θ‖) +
            (∫ θ : ℝ, ‖iteratedDeriv (r.M - 1 + 2) G θ‖) / R₁)) := by
  obtain ⟨c, hc, hTs⟩ := turan_shift hT
  refine ⟨c, hc, fun r G hG ℓ K x R₀ R₁ hℓ hM hx hR₀ hR₁ hK => ?_⟩
  set m := r.M - 1
  have hK0 : 0 ≤ K := (norm_nonneg _).trans (hK 0 ⟨le_rfl, hℓ.le⟩)
  rw [← r.integral_kern_mul_eval hG x]
  have hpt : ∀ v : ℝ, ‖kern G v * r.eval (x + v)‖ ≤
      c ^ m * K * 2 ^ m * (‖kern G v‖ + (ℓ ^ m)⁻¹ * (|v| ^ m * ‖kern G v‖)) := by
    intro v
    rw [norm_mul]
    have h1 := hTs r ℓ K x v hℓ hM hx hK
    have h2 : (c * (1 + |v| / ℓ)) ^ m ≤ c ^ m * (2 ^ m * (1 + (|v| / ℓ) ^ m)) := by
      rw [mul_pow]; gcongr; exact one_add_pow_le m (by positivity)
    have h3 : ‖r.eval (x + v)‖ ≤ c ^ m * (2 ^ m * (1 + (|v| / ℓ) ^ m)) * K :=
      h1.trans (mul_le_mul_of_nonneg_right h2 hK0)
    calc ‖kern G v‖ * ‖r.eval (x + v)‖ ≤ ‖kern G v‖ * (c ^ m * (2 ^ m * (1 + (|v| / ℓ) ^ m)) * K) := by
          gcongr
      _ = c ^ m * K * 2 ^ m * (‖kern G v‖ + (ℓ ^ m)⁻¹ * (|v| ^ m * ‖kern G v‖)) := by
          rw [div_pow]; field_simp
  have hint0 : Integrable fun v : ℝ => ‖kern G v‖ := hG.integrable_kern.norm
  have hintm := hG.integrable_norm_pow_mul_kern m
  have hmaj : Integrable fun v : ℝ =>
      c ^ m * K * 2 ^ m * (‖kern G v‖ + (ℓ ^ m)⁻¹ * (|v| ^ m * ‖kern G v‖)) :=
    (hint0.add (hintm.const_mul _)).const_mul _
  calc ‖∫ v : ℝ, kern G v * r.eval (x + v)‖ ≤ ∫ v : ℝ, ‖kern G v * r.eval (x + v)‖ :=
        norm_integral_le_integral_norm _
    _ ≤ ∫ v : ℝ, c ^ m * K * 2 ^ m * (‖kern G v‖ + (ℓ ^ m)⁻¹ * (|v| ^ m * ‖kern G v‖)) :=
        integral_mono_of_nonneg (Eventually.of_forall fun v => norm_nonneg _) hmaj
          (Eventually.of_forall hpt)
    _ = c ^ m * K * 2 ^ m * ((∫ v : ℝ, ‖kern G v‖) +
          (ℓ ^ m)⁻¹ * ∫ v : ℝ, |v| ^ m * ‖kern G v‖) := by
        rw [integral_const_mul, integral_add hint0 (hintm.const_mul _), integral_const_mul]
    _ ≤ c ^ m * K * 2 ^ m * ((R₀ * (∫ θ : ℝ, ‖G θ‖) + (∫ θ : ℝ, ‖iteratedDeriv 2 G θ‖) / R₀) +
          (ℓ ^ m)⁻¹ * (R₁ * (∫ θ : ℝ, ‖iteratedDeriv m G θ‖) +
            (∫ θ : ℝ, ‖iteratedDeriv (m + 2) G θ‖) / R₁)) := by
        have h0 := hG.integral_pow_mul_norm_kern_le 0 hR₀
        simp only [pow_zero, one_mul, iteratedDeriv_zero, zero_add] at h0
        have hm := hG.integral_pow_mul_norm_kern_le m hR₁
        gcongr
    _ = _ := by ring

end Riemann.Deep.LeanNBRate
