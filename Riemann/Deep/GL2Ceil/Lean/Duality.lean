import Riemann.Deep.GL2Ceil.Lean.Defs
import Mathlib.MeasureTheory.Integral.DivergenceTheorem

/-!
# GL2Ceil (gl2lean): Lemma D, weak duality (DERIVATION §2)

For R, α, γ > 0, σ continuous on [0,1] with derivative σ' off a countable set (continuous piecewise
polynomials included), σ', μ integrable with (σ' + μ e^{-Ry})² integrable, μ(y) = μ(1-y), M = ∫₀¹ μ:
  J(Q) = (1 + w(1)²)/2 + α∫w² + γ∫w'² ≥ 1/2 + M - 2σ(0) - (2σ(1) + M e^{-R})²/2
          - (1/α)∫(σ' + μe^{-Ry})² - (1/γ)∫σ²   for every Q in the real class (w = e^{Ry}Q).
Proof: pointwise Fenchel, integration by parts (FTC off a countable set), folding y ↦ 1-y, minimum
over u = β - 1 (a square).
-/

open Real MeasureTheory Set

noncomputable section

namespace GL2Ceil

/-- p = σ' + μ e^{-Ry}. -/
def pfun (R : ℝ) (σ' μ : ℝ → ℝ) (y : ℝ) : ℝ := σ' y + μ y * Real.exp (-(R * y))

/-- The dual value D of Lemma D (DERIVATION §2). -/
def Dual (R α γ : ℝ) (σ σ' μ : ℝ → ℝ) : ℝ :=
  1 / 2 + (∫ y in (0:ℝ)..1, μ y) - 2 * σ 0
    - (2 * σ 1 + (∫ y in (0:ℝ)..1, μ y) * Real.exp (-R)) ^ 2 / 2
    - (1 / α) * (∫ y in (0:ℝ)..1, pfun R σ' μ y ^ 2)
    - (1 / γ) * ∫ y in (0:ℝ)..1, σ y ^ 2

/-- Admissible dual pair (σ, μ), with σ' the derivative of σ off a countable set. -/
structure DualAdm (σ σ' μ : ℝ → ℝ) : Prop where
  cont : ContinuousOn σ (Icc 0 1)
  hasDeriv : ∃ s : Set ℝ, s.Countable ∧ ∀ x ∈ Ioo (0:ℝ) 1 \ s, HasDerivAt σ (σ' x) x
  int_σ' : IntervalIntegrable σ' volume 0 1
  int_μ : IntervalIntegrable μ volume 0 1
  int_sq : ∀ R : ℝ, IntervalIntegrable (fun y => pfun R σ' μ y ^ 2) volume 0 1
  symm : ∀ y ∈ Icc (0:ℝ) 1, μ (1 - y) = μ y

/-- Lemma D (weak duality): M(R, λ) ≥ D. -/
theorem weak_duality {R α γ : ℝ} (hα : 0 < α) (hγ : 0 < γ) {σ σ' μ : ℝ → ℝ}
    (hd : DualAdm σ σ' μ) {Q : ℝ → ℝ} (hQ : QAdm Q) :
    Dual R α γ σ σ' μ ≤ Jag R α γ Q := by
  obtain ⟨hQd, hQc⟩ := contDiff_one_iff_deriv.mp hQ.smooth
  obtain ⟨β, hβ⟩ := hQ.symm
  have hQ0c : Continuous Q := hQd.continuous
  have hwc := w_cont (R := R) hQ0c
  have hwdc := wd_cont (R := R) hQ0c hQc
  have hQ1 : Q 1 = β - 1 := by
    have := hβ 0 ⟨le_refl 0, zero_le_one⟩
    rw [sub_zero, hQ.zero] at this
    linarith
  have hw0 : w R Q 0 = 1 := by simp [w, hQ.zero]
  have hw1 : w R Q 1 * Real.exp (-R) = β - 1 := by
    simp only [w, mul_one]
    rw [hQ1]
    have : Real.exp R * Real.exp (-R) = 1 := by rw [← Real.exp_add]; simp
    linear_combination (β - 1) * this
  -- integrability
  have hexpc : Continuous (fun y : ℝ => Real.exp (-(R * y))) := by fun_prop
  have hμe : IntervalIntegrable (fun y => μ y * Real.exp (-(R * y))) volume 0 1 :=
    hd.int_μ.mul_continuousOn hexpc.continuousOn
  have hp : IntervalIntegrable (pfun R σ' μ) volume 0 1 := hd.int_σ'.add hμe
  have hpw : IntervalIntegrable (fun y => pfun R σ' μ y * w R Q y) volume 0 1 :=
    hp.mul_continuousOn hwc.continuousOn
  have hp2 := hd.int_sq R
  have hw2 : IntervalIntegrable (fun y => w R Q y ^ 2) volume 0 1 :=
    (hwc.fun_pow 2).intervalIntegrable 0 1
  have hwd2 : IntervalIntegrable (fun y => wd R Q y ^ 2) volume 0 1 :=
    (hwdc.fun_pow 2).intervalIntegrable 0 1
  have hσwd : IntervalIntegrable (fun y => σ y * wd R Q y) volume 0 1 :=
    (hd.cont.fun_mul hwdc.continuousOn).intervalIntegrable_of_Icc zero_le_one
  have hσ2 : IntervalIntegrable (fun y => σ y ^ 2) volume 0 1 :=
    (hd.cont.fun_pow 2).intervalIntegrable_of_Icc zero_le_one
  have hσ'w : IntervalIntegrable (fun y => σ' y * w R Q y) volume 0 1 :=
    hd.int_σ'.mul_continuousOn hwc.continuousOn
  have hμQ : IntervalIntegrable (fun y => μ y * Q y) volume 0 1 :=
    hd.int_μ.mul_continuousOn hQ0c.continuousOn
  have hQr : Continuous (fun y => Q (1 - y)) := by fun_prop
  have hμQ' : IntervalIntegrable (fun y => μ y * Q (1 - y)) volume 0 1 :=
    hd.int_μ.mul_continuousOn hQr.continuousOn
  -- Fenchel for α ∫ w²
  have F1 : 2 * (∫ y in (0:ℝ)..1, pfun R σ' μ y * w R Q y)
      - (1 / α) * (∫ y in (0:ℝ)..1, pfun R σ' μ y ^ 2) ≤ α * W0 R Q := by
    have key : ∫ y in (0:ℝ)..1, (2 * (pfun R σ' μ y * w R Q y) - (1 / α) * pfun R σ' μ y ^ 2)
        ≤ ∫ y in (0:ℝ)..1, α * w R Q y ^ 2 := by
      apply intervalIntegral.integral_mono_on zero_le_one
        ((hpw.const_mul 2).sub (hp2.const_mul _)) (hw2.const_mul _)
      intro y _
      have h1 : 0 ≤ (α * w R Q y - pfun R σ' μ y) ^ 2 / α := div_nonneg (sq_nonneg _) hα.le
      have e : (α * w R Q y - pfun R σ' μ y) ^ 2 / α
          = α * w R Q y ^ 2 - (2 * (pfun R σ' μ y * w R Q y) - (1 / α) * pfun R σ' μ y ^ 2) := by
        field_simp
        ring
      linarith
    rw [intervalIntegral.integral_sub (hpw.const_mul 2) (hp2.const_mul _),
      intervalIntegral.integral_const_mul, intervalIntegral.integral_const_mul,
      intervalIntegral.integral_const_mul] at key
    unfold W0
    linarith
  -- Fenchel for γ ∫ w'²
  have F2 : 2 * (∫ y in (0:ℝ)..1, σ y * wd R Q y)
      - (1 / γ) * (∫ y in (0:ℝ)..1, σ y ^ 2) ≤ γ * W1 R Q := by
    have key : ∫ y in (0:ℝ)..1, (2 * (σ y * wd R Q y) - (1 / γ) * σ y ^ 2)
        ≤ ∫ y in (0:ℝ)..1, γ * wd R Q y ^ 2 := by
      apply intervalIntegral.integral_mono_on zero_le_one
        ((hσwd.const_mul 2).sub (hσ2.const_mul _)) (hwd2.const_mul _)
      intro y _
      have h1 : 0 ≤ (γ * wd R Q y - σ y) ^ 2 / γ := div_nonneg (sq_nonneg _) hγ.le
      have e : (γ * wd R Q y - σ y) ^ 2 / γ
          = γ * wd R Q y ^ 2 - (2 * (σ y * wd R Q y) - (1 / γ) * σ y ^ 2) := by
        field_simp
        ring
      linarith
    rw [intervalIntegral.integral_sub (hσwd.const_mul 2) (hσ2.const_mul _),
      intervalIntegral.integral_const_mul, intervalIntegral.integral_const_mul,
      intervalIntegral.integral_const_mul] at key
    unfold W1
    linarith
  -- integration by parts: ∫ σ' w + ∫ σ w' = σ(1) w(1) - σ(0) w(0)
  have IBP : (∫ y in (0:ℝ)..1, σ' y * w R Q y) + (∫ y in (0:ℝ)..1, σ y * wd R Q y)
      = σ 1 * w R Q 1 - σ 0 * w R Q 0 := by
    obtain ⟨s, hs, hsd⟩ := hd.hasDeriv
    rw [← intervalIntegral.integral_add hσ'w hσwd]
    exact MeasureTheory.integral_eq_of_hasDerivAt_off_countable_of_le (fun y => σ y * w R Q y) _
      zero_le_one hs (hd.cont.fun_mul hwc.continuousOn)
      (fun x hx => (hsd x hx).fun_mul (hasDerivAt_w (R := R) hQd x)) (hσ'w.add hσwd)
  rw [hw0, mul_one] at IBP
  -- ∫ p w = ∫ σ' w + ∫ μ Q
  have Split : (∫ y in (0:ℝ)..1, pfun R σ' μ y * w R Q y)
      = (∫ y in (0:ℝ)..1, σ' y * w R Q y) + ∫ y in (0:ℝ)..1, μ y * Q y := by
    rw [← intervalIntegral.integral_add hσ'w hμQ]
    apply intervalIntegral.integral_congr
    intro y _
    simp only [pfun, w]
    have : Real.exp (-(R * y)) * Real.exp (R * y) = 1 := by rw [← Real.exp_add]; simp
    linear_combination (μ y * Q y) * this
  -- folding y ↦ 1 - y: 2 ∫ μ Q = β M
  have Fold : 2 * (∫ y in (0:ℝ)..1, μ y * Q y) = β * ∫ y in (0:ℝ)..1, μ y := by
    have hrefl : (∫ y in (0:ℝ)..1, μ y * Q y) = ∫ y in (0:ℝ)..1, μ y * Q (1 - y) := by
      have h := intervalIntegral.integral_comp_sub_left (fun y => μ y * Q y) (a := 0) (b := 1) 1
      simp only [sub_self, sub_zero] at h
      rw [← h]
      apply intervalIntegral.integral_congr
      intro y hy
      rw [uIcc_of_le zero_le_one] at hy
      simp only [hd.symm y hy]
    calc 2 * (∫ y in (0:ℝ)..1, μ y * Q y)
        = (∫ y in (0:ℝ)..1, μ y * Q y) + ∫ y in (0:ℝ)..1, μ y * Q (1 - y) := by
          rw [← hrefl]; ring
      _ = ∫ y in (0:ℝ)..1, (μ y * Q y + μ y * Q (1 - y)) :=
          (intervalIntegral.integral_add hμQ hμQ').symm
      _ = ∫ y in (0:ℝ)..1, β * μ y := by
          apply intervalIntegral.integral_congr
          intro y hy
          rw [uIcc_of_le zero_le_one] at hy
          simp only
          rw [← mul_add, hβ y hy, mul_comm]
      _ = β * ∫ y in (0:ℝ)..1, μ y := intervalIntegral.integral_const_mul _ _
  have hw1M : (β - 1) * (∫ y in (0:ℝ)..1, μ y)
      = w R Q 1 * Real.exp (-R) * (∫ y in (0:ℝ)..1, μ y) := by rw [hw1]
  have E : Jag R α γ Q - Dual R α γ σ σ' μ
      = (α * W0 R Q - (2 * (∫ y in (0:ℝ)..1, pfun R σ' μ y * w R Q y)
          - (1 / α) * (∫ y in (0:ℝ)..1, pfun R σ' μ y ^ 2)))
        + (γ * W1 R Q - (2 * (∫ y in (0:ℝ)..1, σ y * wd R Q y)
          - (1 / γ) * (∫ y in (0:ℝ)..1, σ y ^ 2)))
        + (w R Q 1 + 2 * σ 1 + (∫ y in (0:ℝ)..1, μ y) * Real.exp (-R)) ^ 2 / 2 := by
    unfold Jag Dual
    linear_combination 2 * Split + 2 * IBP + Fold + hw1M
  nlinarith [F1, F2, E, sq_nonneg (w R Q 1 + 2 * σ 1 + (∫ y in (0:ℝ)..1, μ y) * Real.exp (-R))]

end GL2Ceil

