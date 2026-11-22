import Riemann.Deep.LeanNBRate.KJ2

/-!
# K = zeta J to the right of the contour (Re s > 3)
-/

open Complex Filter Topology Set MeasureTheory intervalIntegral Metric
open scoped Real Interval

noncomputable section

namespace Riemann.Deep.LeanNBRate

variable {H : ℂ → ℂ} {M β D : ℝ} {N : ℕ}

/-- Integrability on the three pieces of the upper half-strip boundary. -/
structure UInt (f : ℂ → ℂ) (a b T : ℝ) : Prop where
  bot : IntervalIntegrable (fun x : ℝ => f (x + T * I)) volume a b
  rgt : IntegrableOn (fun y : ℝ => f (b + y * I)) (Ioi T)
  lft : IntegrableOn (fun y : ℝ => f (a + y * I)) (Ioi T)

/-- Integrability on the three pieces of the lower half-strip boundary. -/
structure LInt (f : ℂ → ℂ) (a b T : ℝ) : Prop where
  bot : IntervalIntegrable (fun x : ℝ => f (x + ((-T : ℝ) : ℂ) * I)) volume a b
  rgt : IntegrableOn (fun y : ℝ => f (b + y * I)) (Iic (-T))
  lft : IntegrableOn (fun y : ℝ => f (a + y * I)) (Iic (-T))

theorem UpperU_sub {f g : ℂ → ℂ} {a b T : ℝ} (hf : UInt f a b T) (hg : UInt g a b T) :
    UpperU (fun w => f w - g w) a b T = UpperU f a b T - UpperU g a b T := by
  simp only [UpperU]
  rw [intervalIntegral.integral_sub hf.bot hg.bot, integral_sub hf.rgt hg.rgt, integral_sub hf.lft hg.lft]
  ring

theorem LowerU_sub {f g : ℂ → ℂ} {a b T : ℝ} (hf : LInt f a b T) (hg : LInt g a b T) :
    LowerU (fun w => f w - g w) a b T = LowerU f a b T - LowerU g a b T := by
  simp only [LowerU]
  rw [intervalIntegral.integral_sub hf.bot hg.bot, integral_sub hf.rgt hg.rgt, integral_sub hf.lft hg.lft]
  ring

theorem UpperU_const_mul (c : ℂ) (f : ℂ → ℂ) (a b T : ℝ) :
    UpperU (fun w => c * f w) a b T = c * UpperU f a b T := by
  simp only [UpperU, intervalIntegral.integral_const_mul, MeasureTheory.integral_const_mul]
  ring

theorem LowerU_const_mul (c : ℂ) (f : ℂ → ℂ) (a b T : ℝ) :
    LowerU (fun w => c * f w) a b T = c * LowerU f a b T := by
  simp only [LowerU, intervalIntegral.integral_const_mul, MeasureTheory.integral_const_mul]
  ring

theorem UpperU_congr {f g : ℂ → ℂ} {a b T : ℝ} (hab : a ≤ b)
    (hbot : ∀ x ∈ Icc a b, f (x + T * I) = g (x + T * I))
    (hr : ∀ y ∈ Ioi T, f (b + y * I) = g (b + y * I)) (hl : ∀ y ∈ Ioi T, f (a + y * I) = g (a + y * I)) :
    UpperU f a b T = UpperU g a b T := by
  simp only [UpperU]
  rw [intervalIntegral.integral_congr (fun x hx => hbot x (by rwa [Set.uIcc_of_le hab] at hx)),
    setIntegral_congr_fun measurableSet_Ioi hr, setIntegral_congr_fun measurableSet_Ioi hl]

theorem LowerU_congr {f g : ℂ → ℂ} {a b T : ℝ} (hab : a ≤ b)
    (hbot : ∀ x ∈ Icc a b, f (x + ((-T : ℝ) : ℂ) * I) = g (x + ((-T : ℝ) : ℂ) * I))
    (hr : ∀ y ∈ Iic (-T), f (b + y * I) = g (b + y * I)) (hl : ∀ y ∈ Iic (-T), f (a + y * I) = g (a + y * I)) :
    LowerU f a b T = LowerU g a b T := by
  simp only [LowerU]
  rw [intervalIntegral.integral_congr (fun x hx => hbot x (by rwa [Set.uIcc_of_le hab] at hx)),
    setIntegral_congr_fun measurableSet_Iic hr, setIntegral_congr_fun measurableSet_Iic hl]

/-- The difference zeta(s) Psi_s - kint_s on the contour. -/
def hdiff (H : ℂ → ℂ) (L : ℝ) (s w : ℂ) : ℂ := H w * kap L (s - w) / (s - w)

theorem zeta_Psi_eq {L : ℝ} {s w : ℂ} (hw : riemannZeta w ≠ 0) (hsw : s ≠ w) :
    riemannZeta s * Psi H L s w = kint H L s w + hdiff H L s w := by
  have hsw' : s - w ≠ 0 := sub_ne_zero.mpr hsw
  rw [Psi, Fq, kint, hdiff, dslope_of_ne _ hsw, slope_def_field]
  field_simp
  ring

/-- Bound for functions of the type H(w) * (bounded) on the contour pieces. -/
theorem norm_comb_le (hH : HBound H M β D N) {g : ℂ → ℂ} {C : ℝ} (hC : 0 ≤ C)
    (hg : ∀ w : ℂ, -1 ≤ w.re → w.re ≤ 2 → ‖g w‖ ≤ C * ‖H w‖) {σ : ℝ} (h1 : -1 ≤ σ) (h2 : σ ≤ 2) (y : ℝ) :
    ‖g (σ + y * I)‖ ≤ C * (M * Real.exp (β * max 0 (1 / 2 - σ)) * (D / min D 1) ^ 2 * ((1 + |y|) ^ 2)⁻¹) := by
  have h := hg (σ + y * I) (by simp; exact h1) (by simp; exact h2)
  exact h.trans (mul_le_mul_of_nonneg_left (norm_H_le_line hH σ y) hC)

end Riemann.Deep.LeanNBRate

