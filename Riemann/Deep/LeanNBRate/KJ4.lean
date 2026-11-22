import Riemann.Deep.LeanNBRate.KJ3

/-!
# K = zeta J for Re s > 3
-/

open Complex Filter Topology Set MeasureTheory intervalIntegral Metric
open scoped Real Interval

noncomputable section

namespace Riemann.Deep.LeanNBRate

variable {H : ℂ → ℂ} {M β D : ℝ} {N : ℕ}

theorem norm_Psi_le_on {L : ℝ} (hL : 0 ≤ L) {s w : ℂ} (hs : 3 < s.re) (hw : w.re ≤ 2) {B : ℝ}
    (hB : ‖(riemannZeta w)⁻¹‖ ≤ B) : ‖Psi H L s w‖ ≤ B * ‖H w‖ := by
  rw [norm_Psi_eq]
  have hsw : 1 ≤ ‖s - w‖ := by
    have := Complex.abs_re_le_norm (s - w); simp at this
    have h' : s.re - w.re ≤ |s.re - w.re| := le_abs_self _
    linarith
  have hk : Real.exp (-(L * (s - w).re)) ≤ 1 := by
    rw [Real.exp_le_one_iff]; simp; nlinarith
  have hB0 : 0 ≤ B := (norm_nonneg _).trans hB
  calc ‖H w‖ * ‖(riemannZeta w)⁻¹‖ * Real.exp (-(L * (s - w).re)) / ‖s - w‖
      ≤ ‖H w‖ * ‖(riemannZeta w)⁻¹‖ * Real.exp (-(L * (s - w).re)) := div_le_self (by positivity) hsw
    _ ≤ ‖H w‖ * B * 1 := by gcongr
    _ = B * ‖H w‖ := by ring

theorem norm_hdiff_le {L : ℝ} (hL : 0 ≤ L) {s w : ℂ} (hs : 3 < s.re) (hw : w.re ≤ 2) :
    ‖hdiff H L s w‖ ≤ 1 * ‖H w‖ := by
  rw [hdiff, norm_div, norm_mul, one_mul]
  have hsw : 1 ≤ ‖s - w‖ := by
    have := Complex.abs_re_le_norm (s - w); simp at this
    have h' : s.re - w.re ≤ |s.re - w.re| := le_abs_self _
    linarith
  have hk : ‖kap L (s - w)‖ ≤ 1 := norm_kap_le_one hL (by simp; linarith)
  calc ‖H w‖ * ‖kap L (s - w)‖ / ‖s - w‖ ≤ ‖H w‖ * ‖kap L (s - w)‖ := div_le_self (by positivity) hsw
    _ ≤ ‖H w‖ * 1 := mul_le_mul_of_nonneg_left hk (norm_nonneg _)
    _ = ‖H w‖ := mul_one _

theorem integrableOn_line_of_le (hH : HBound H M β D N) {f : ℂ → ℂ} {σ : ℝ} {S : Set ℝ}
    (hS : MeasurableSet S) (hc : ContinuousOn (fun y : ℝ => f (σ + y * I)) S) {B : ℝ}
    (hb : ∀ y ∈ S, ‖f (σ + y * I)‖ ≤ B * ‖H (σ + y * I)‖) : IntegrableOn (fun y : ℝ => f (σ + y * I)) S := by
  refine ((integrable_norm_H_line hH σ).integrableOn.const_mul B).mono' (hc.aestronglyMeasurable hS) ?_
  filter_upwards [ae_restrict_mem hS] with y hy
  exact hb y hy

theorem continuousAt_Psi_contour (hH : HBound H M β D N) (L : ℝ) {s w : ℂ} (hs : 3 < s.re) (hw : w.re ≤ 2)
    (hw1 : w ≠ 1) (hw0 : riemannZeta w ≠ 0) : ContinuousAt (Psi H L s) w :=
  continuousAt_Psi hH L hw1 hw0 (fun h => by rw [h] at hw; linarith)

theorem continuousAt_hdiff (hH : HBound H M β D N) (L : ℝ) {s w : ℂ} (hsw : w ≠ s) :
    ContinuousAt (hdiff H L s) w := by
  unfold hdiff kap
  exact ((hH.diff w).continuousAt.mul (by fun_prop)).div (continuousAt_const.sub continuousAt_id)
    (sub_ne_zero.mpr (Ne.symm hsw))

theorem differentiableAt_hdiff (hH : HBound H M β D N) (L : ℝ) {s w : ℂ} (hsw : w ≠ s) :
    DifferentiableAt ℂ (hdiff H L s) w := by
  unfold hdiff kap
  exact ((hH.diff w).mul (by fun_prop)).div ((differentiableAt_const _).sub differentiableAt_id)
    (sub_ne_zero.mpr (Ne.symm hsw))

theorem ne_of_re_lt {s w : ℂ} (h : w.re < s.re) : w ≠ s := fun e => by rw [e] at h; exact lt_irrefl _ h

/-- **K = zeta J** for Re s > 3. -/
theorem Kf_eq_zeta_Jf (hH : HBound H M β D N) {T₂ : ℝ} (hZ : HZeros H T₂) (hT1 : 1 ≤ T₂) {L : ℝ}
    (hL : 0 ≤ L) {s : ℂ} (hs : 3 < s.re) : Kf H L T₂ s = riemannZeta s * Jf H L T₂ s := by
  have hT0 : T₂ ≠ 0 := by linarith
  set ζs := riemannZeta s
  obtain ⟨B₁, hB₁0, hB₁⟩ := exists_bound_bottom hT0 hZ.good
  obtain ⟨B₂, hB₂0, hB₂⟩ := exists_bound_bottom (T := -T₂) (by linarith) (by
    intro σ hσ; have := hZ.good_neg σ hσ; simpa using this)
  set Bl := 60 * π * Real.sqrt (2 * π) with hBl
  have hBl0 : 0 ≤ Bl := by positivity
  -- pointwise identity on the contour
  have hid : ∀ w : ℂ, w.re ≤ 2 → riemannZeta w ≠ 0 →
      kint H L s w = ζs * Psi H L s w - hdiff H L s w := by
    intro w hw hw0
    rw [zeta_Psi_eq hw0 (Ne.symm (ne_of_re_lt (by linarith)))]; ring
  -- facts on the pieces
  have hbotU : ∀ x ∈ Icc (-1 : ℝ) 2, riemannZeta ((x : ℂ) + T₂ * I) ≠ 0 := hZ.good
  have hbotL : ∀ x ∈ Icc (-1 : ℝ) 2, riemannZeta ((x : ℂ) + ((-T₂ : ℝ) : ℂ) * I) ≠ 0 := hZ.good_neg
  have hz2 : ∀ y : ℝ, riemannZeta (((2 : ℝ) : ℂ) + y * I) ≠ 0 := fun y =>
    riemannZeta_ne_zero_of_one_lt_re (by simp)
  have hzl : ∀ y : ℝ, y ≠ 0 → riemannZeta (((-1 : ℝ) : ℂ) + y * I) ≠ 0 := fun y hy =>
    zeta_ne_zero_neg_one_line hy
  -- integrability: Psi
  have hcontPsi : ∀ w : ℂ, w.re ≤ 2 → w ≠ 1 → riemannZeta w ≠ 0 → ContinuousAt (fun w => ζs * Psi H L s w) w :=
    fun w hw hw1 hw0 => continuousAt_const.mul (continuousAt_Psi_contour hH L hs hw hw1 hw0)
  have hUP : UInt (fun w => ζs * Psi H L s w) (-1) 2 T₂ := by
    refine ⟨?_, ?_, ?_⟩
    · refine ContinuousOn.intervalIntegrable fun x hx => ?_
      rw [Set.uIcc_of_le (by norm_num)] at hx
      exact ContinuousAt.continuousWithinAt (ContinuousAt.comp (f := fun x : ℝ => (x : ℂ) + T₂ * I) (x := x)
        (hcontPsi _ (by simp; exact hx.2) (line_ne_one hT0) (hbotU x hx)) (by fun_prop))
    · refine integrableOn_line_of_le hH (f := fun w => ζs * Psi H L s w) (σ := 2) measurableSet_Ioi
        (fun y _ => ?_) (B := ‖ζs‖ * 3) fun y _ => ?_
      · exact ContinuousAt.continuousWithinAt (ContinuousAt.comp (f := fun y : ℝ => ((2 : ℝ) : ℂ) + y * I)
          (x := y) (hcontPsi _ (by simp) (by intro h; have := congrArg Complex.re h; simp at this) (hz2 y))
          (by fun_prop))
      · rw [norm_mul, mul_assoc]
        exact mul_le_mul_of_nonneg_left (norm_Psi_le_on hL hs (by simp) (norm_inv_zeta_le_three (by simp)))
          (norm_nonneg _)
    · refine integrableOn_line_of_le hH (f := fun w => ζs * Psi H L s w) (σ := -1) measurableSet_Ioi
        (fun y hy => ?_) (B := ‖ζs‖ * Bl) fun y hy => ?_
      · have hy' : (0 : ℝ) < y := by have : T₂ < y := hy; linarith
        exact ContinuousAt.continuousWithinAt (ContinuousAt.comp (f := fun y : ℝ => ((-1 : ℝ) : ℂ) + y * I)
          (x := y) (hcontPsi _ (by simp; norm_num) (line_ne_one hy'.ne') (hzl y hy'.ne')) (by fun_prop))
      · have hy' : T₂ < y := hy
        rw [norm_mul, mul_assoc]
        exact mul_le_mul_of_nonneg_left (norm_Psi_le_on hL hs (by simp; norm_num)
          (norm_inv_zeta_neg_one_line (by rw [abs_of_pos (by linarith)]; linarith))) (norm_nonneg _)
  have hLP : LInt (fun w => ζs * Psi H L s w) (-1) 2 T₂ := by
    refine ⟨?_, ?_, ?_⟩
    · refine ContinuousOn.intervalIntegrable fun x hx => ?_
      rw [Set.uIcc_of_le (by norm_num)] at hx
      exact ContinuousAt.continuousWithinAt (ContinuousAt.comp
        (f := fun x : ℝ => (x : ℂ) + ((-T₂ : ℝ) : ℂ) * I) (x := x)
        (hcontPsi _ (by simp; exact hx.2) (line_ne_one (by linarith)) (hbotL x hx)) (by fun_prop))
    · refine integrableOn_line_of_le hH (f := fun w => ζs * Psi H L s w) (σ := 2) measurableSet_Iic
        (fun y _ => ?_) (B := ‖ζs‖ * 3) fun y _ => ?_
      · exact ContinuousAt.continuousWithinAt (ContinuousAt.comp (f := fun y : ℝ => ((2 : ℝ) : ℂ) + y * I)
          (x := y) (hcontPsi _ (by simp) (by intro h; have := congrArg Complex.re h; simp at this) (hz2 y))
          (by fun_prop))
      · rw [norm_mul, mul_assoc]
        exact mul_le_mul_of_nonneg_left (norm_Psi_le_on hL hs (by simp) (norm_inv_zeta_le_three (by simp)))
          (norm_nonneg _)
    · refine integrableOn_line_of_le hH (f := fun w => ζs * Psi H L s w) (σ := -1) measurableSet_Iic
        (fun y hy => ?_) (B := ‖ζs‖ * Bl) fun y hy => ?_
      · have hy' : y < 0 := by have : y ≤ -T₂ := hy; linarith
        exact ContinuousAt.continuousWithinAt (ContinuousAt.comp (f := fun y : ℝ => ((-1 : ℝ) : ℂ) + y * I)
          (x := y) (hcontPsi _ (by simp; norm_num) (line_ne_one hy'.ne) (hzl y hy'.ne)) (by fun_prop))
      · have hy' : y ≤ -T₂ := hy
        rw [norm_mul, mul_assoc]
        exact mul_le_mul_of_nonneg_left (norm_Psi_le_on hL hs (by simp; norm_num)
          (norm_inv_zeta_neg_one_line (by rw [abs_of_neg (by linarith)]; linarith))) (norm_nonneg _)
  -- integrability: hdiff
  have hcd : ∀ w : ℂ, w.re ≤ 2 → ContinuousAt (hdiff H L s) w := fun w hw =>
    continuousAt_hdiff hH L (ne_of_re_lt (by linarith))
  have hUH : UInt (hdiff H L s) (-1) 2 T₂ := by
    refine ⟨?_, ?_, ?_⟩
    · refine ContinuousOn.intervalIntegrable fun x hx => ?_
      rw [Set.uIcc_of_le (by norm_num)] at hx
      exact ContinuousAt.continuousWithinAt (ContinuousAt.comp (f := fun x : ℝ => (x : ℂ) + T₂ * I) (x := x)
        (hcd _ (by simp; exact hx.2)) (by fun_prop))
    · refine integrableOn_line_of_le hH (f := hdiff H L s) (σ := 2) measurableSet_Ioi (fun y _ => ?_)
        (B := 1) fun y _ => norm_hdiff_le hL hs (by simp)
      exact ContinuousAt.continuousWithinAt (ContinuousAt.comp (f := fun y : ℝ => ((2 : ℝ) : ℂ) + y * I) (x := y)
        (hcd _ (by simp)) (by fun_prop))
    · refine integrableOn_line_of_le hH (f := hdiff H L s) (σ := -1) measurableSet_Ioi (fun y _ => ?_)
        (B := 1) fun y _ => norm_hdiff_le hL hs (by simp; norm_num)
      exact ContinuousAt.continuousWithinAt (ContinuousAt.comp (f := fun y : ℝ => ((-1 : ℝ) : ℂ) + y * I) (x := y)
        (hcd _ (by simp; norm_num)) (by fun_prop))
  have hLH : LInt (hdiff H L s) (-1) 2 T₂ := by
    refine ⟨?_, ?_, ?_⟩
    · refine ContinuousOn.intervalIntegrable fun x hx => ?_
      rw [Set.uIcc_of_le (by norm_num)] at hx
      exact ContinuousAt.continuousWithinAt (ContinuousAt.comp
        (f := fun x : ℝ => (x : ℂ) + ((-T₂ : ℝ) : ℂ) * I) (x := x) (hcd _ (by simp; exact hx.2)) (by fun_prop))
    · refine integrableOn_line_of_le hH (f := hdiff H L s) (σ := 2) measurableSet_Iic (fun y _ => ?_)
        (B := 1) fun y _ => norm_hdiff_le hL hs (by simp)
      exact ContinuousAt.continuousWithinAt (ContinuousAt.comp (f := fun y : ℝ => ((2 : ℝ) : ℂ) + y * I) (x := y)
        (hcd _ (by simp)) (by fun_prop))
    · refine integrableOn_line_of_le hH (f := hdiff H L s) (σ := -1) measurableSet_Iic (fun y _ => ?_)
        (B := 1) fun y _ => norm_hdiff_le hL hs (by simp; norm_num)
      exact ContinuousAt.continuousWithinAt (ContinuousAt.comp (f := fun y : ℝ => ((-1 : ℝ) : ℂ) + y * I) (x := y)
        (hcd _ (by simp; norm_num)) (by fun_prop))
  -- the hdiff half-strip integrals vanish
  have hdecH : ∀ ε > 0, ∀ᶠ R : ℝ in atTop, ∀ x ∈ Icc (-1 : ℝ) 2, ‖hdiff H L s (x + R * I)‖ ≤ ε ∧
      ‖hdiff H L s (x + ((-R : ℝ) : ℂ) * I)‖ ≤ ε := by
    refine decay_of_le_inv_sq (C := 1 * (M * Real.exp (β * (3 / 2)) * (D / min D 1) ^ 2)) fun z h1 h2 _ => ?_
    have hz : z = (z.re : ℂ) + z.im * I := (Complex.re_add_im z).symm
    have := norm_comb_le hH (g := hdiff H L s) zero_le_one (fun w hw1 hw2 => norm_hdiff_le hL hs hw2) h1 h2 z.im
    rw [← hz] at this
    have hM := hH.M_nonneg
    refine this.trans ?_
    have hexp : Real.exp (β * max 0 (1 / 2 - z.re)) ≤ Real.exp (β * (3 / 2)) := by
      apply Real.exp_le_exp.mpr
      exact mul_le_mul_of_nonneg_left (max_le (by norm_num) (by linarith)) hH.β_nonneg
    have : (0 : ℝ) ≤ ((1 + |z.im|) ^ 2)⁻¹ := by positivity
    nlinarith [mul_le_mul_of_nonneg_right hexp (by positivity : (0 : ℝ) ≤ M * (D / min D 1) ^ 2 * ((1 + |z.im|) ^ 2)⁻¹)]
  have hU0 : UpperU (hdiff H L s) (-1) 2 T₂ = 0 :=
    upperU_eq_zero (by norm_num) (fun x _ hx2 _ => hcd x hx2)
      (fun x _ hx2 _ => differentiableAt_hdiff hH L (ne_of_re_lt (by linarith))) hUH.lft hUH.rgt
      (fun ε hε => (hdecH ε hε).mono fun R h x hx => (h x hx).1)
  have hL0 : LowerU (hdiff H L s) (-1) 2 T₂ = 0 :=
    lowerU_eq_zero (by norm_num) (fun x _ hx2 _ => hcd x hx2)
      (fun x _ hx2 _ => differentiableAt_hdiff hH L (ne_of_re_lt (by linarith))) hLH.lft hLH.rgt
      (fun ε hε => (hdecH ε hε).mono fun R h x hx => (h x hx).2)
  -- assemble
  have hKU : UpperU (kint H L s) (-1) 2 T₂ = ζs * UpperU (Psi H L s) (-1) 2 T₂ := by
    rw [UpperU_congr (g := fun w => ζs * Psi H L s w - hdiff H L s w) (by norm_num)
      (fun x hx => hid _ (by simp; exact hx.2) (hbotU x hx))
      (fun y _ => hid _ (by simp) (hz2 y))
      (fun y hy => hid _ (by simp; norm_num) (hzl y (by have : T₂ < y := hy; linarith)))]
    rw [UpperU_sub hUP hUH, hU0, sub_zero, UpperU_const_mul]
  have hKL : LowerU (kint H L s) (-1) 2 T₂ = ζs * LowerU (Psi H L s) (-1) 2 T₂ := by
    rw [LowerU_congr (g := fun w => ζs * Psi H L s w - hdiff H L s w) (by norm_num)
      (fun x hx => hid _ (by simp; exact hx.2) (hbotL x hx))
      (fun y _ => hid _ (by simp) (hz2 y))
      (fun y hy => hid _ (by simp; norm_num) (hzl y (by have : y ≤ -T₂ := hy; linarith)))]
    rw [LowerU_sub hLP hLH, hL0, sub_zero, LowerU_const_mul]
  rw [Kf, Jf, hKU, hKL]
  ring

end Riemann.Deep.LeanNBRate

