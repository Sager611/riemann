import Riemann.Deep.LeanNBRate.KJ1

/-!
# K is analytic off s = 1
-/

open Complex Filter Topology Set MeasureTheory intervalIntegral Metric
open scoped Real Interval

noncomputable section

namespace Riemann.Deep.LeanNBRate

variable {H : ℂ → ℂ} {M β D : ℝ} {N : ℕ}

theorem norm_inv_zeta_neg_one_line {y : ℝ} (hy : 1 ≤ |y|) :
    ‖(riemannZeta (((-1 : ℝ) : ℂ) + y * I))⁻¹‖ ≤ 60 * π * Real.sqrt (2 * π) := by
  set z : ℂ := ((-1 : ℝ) : ℂ) + y * I with hz
  have hre : z.re = -1 := by simp [hz]
  have him : z.im = y := by simp [hz]
  have hzi : 1 ≤ |z.im| := by rw [him]; exact hy
  have h := norm_inv_zeta_left (s := z) (by rw [hre]) (not_nat_of_im hzi) 2 (by rw [hre]; norm_num)
    (by rw [hre]; norm_num)
  rw [hre] at h
  have hP : |y| * |y| ≤ ∏ j ∈ Finset.range 2, ‖z + j‖ := by
    simp only [Finset.prod_range_succ, Finset.prod_range_zero, one_mul, Nat.cast_zero, add_zero, Nat.cast_one]
    have h1 : |y| ≤ ‖z‖ := by rw [← him]; exact Complex.abs_im_le_norm z
    have h2 : |y| ≤ ‖z + 1‖ := by
      have := Complex.abs_im_le_norm (z + 1); simp [him] at this; exact this
    exact mul_le_mul h1 h2 (abs_nonneg _) (norm_nonneg _)
  have hnz : ‖z‖ ≤ 1 + |y| := by
    calc ‖z‖ ≤ |z.re| + |z.im| := Complex.norm_le_abs_re_add_abs_im z
      _ = 1 + |y| := by rw [hre, him]; norm_num
  have hy0 : 0 < |y| := by linarith
  have hyy : 0 < |y| * |y| := by positivity
  have hsq : Real.sqrt (2 * π) ≥ 0 := Real.sqrt_nonneg _
  have h2pi : (2 * π) ^ (-(-1 : ℝ)) = 2 * π := by simp
  rw [h2pi] at h
  have h3 : ‖(riemannZeta z)⁻¹‖ * (|y| * |y|) ≤ 6 * Real.sqrt (2 * π) * (2 * π) * (‖z‖ + 2 + 1) := by
    calc ‖(riemannZeta z)⁻¹‖ * (|y| * |y|) ≤ ‖(riemannZeta z)⁻¹‖ * ∏ j ∈ Finset.range 2, ‖z + j‖ :=
          mul_le_mul_of_nonneg_left hP (norm_nonneg _)
      _ ≤ _ := by simpa using h
  have h4 : 6 * Real.sqrt (2 * π) * (2 * π) * (‖z‖ + 2 + 1) ≤ 60 * π * Real.sqrt (2 * π) * (|y| * |y|) := by
    have : ‖z‖ + 2 + 1 ≤ 5 * (|y| * |y|) := by nlinarith
    have hpi := Real.pi_pos
    calc 6 * Real.sqrt (2 * π) * (2 * π) * (‖z‖ + 2 + 1) ≤ 6 * Real.sqrt (2 * π) * (2 * π) * (5 * (|y| * |y|)) := by
          gcongr
      _ = _ := by ring
  exact le_of_mul_le_mul_right (h3.trans h4) hyy

theorem norm_H_le_line (hH : HBound H M β D N) (σ y : ℝ) :
    ‖H (σ + y * I)‖ ≤ M * Real.exp (β * max 0 (1 / 2 - σ)) * (D / min D 1) ^ 2 * ((1 + |y|) ^ 2)⁻¹ := by
  have h := hH.bound (σ + y * I)
  have hw := weight_mono hH.D_pos N (z := (σ : ℂ) + y * I)
  have hwl := weight_le_inv_sq hH.D_pos (by linarith [hH.N_ge] : 2 ≤ N) y
  have him : ((σ : ℂ) + y * I).im = y := by simp
  have hre : ((σ : ℂ) + y * I).re = σ := by simp
  rw [him] at hw; rw [hre] at h
  have hM := hH.M_nonneg
  calc ‖H (σ + y * I)‖ ≤ M * Real.exp (β * max 0 (1 / 2 - σ)) * (D / (D + ‖(σ : ℂ) + y * I - 1 / 2‖)) ^ N := h
    _ ≤ M * Real.exp (β * max 0 (1 / 2 - σ)) * ((D / min D 1) ^ 2 * ((1 + |y|) ^ 2)⁻¹) := by
        gcongr; exact hw.trans hwl
    _ = _ := by ring

theorem integrable_norm_H_line (hH : HBound H M β D N) (σ : ℝ) :
    Integrable fun y : ℝ => ‖H (σ + y * I)‖ := by
  have hc : Continuous fun y : ℝ => ‖H (σ + y * I)‖ := (hH.diff.continuous.comp (by fun_prop)).norm
  refine (integrable_inv_one_add_abs_sq.const_mul (M * Real.exp (β * max 0 (1 / 2 - σ)) * (D / min D 1) ^ 2)).mono'
    hc.aestronglyMeasurable (Filter.Eventually.of_forall fun y => ?_)
  rw [norm_norm]
  exact norm_H_le_line hH σ y

theorem exists_bound_bottom {T : ℝ} (hT : T ≠ 0) (hgood : ∀ σ ∈ Icc (-1 : ℝ) 2, riemannZeta (σ + T * I) ≠ 0) :
    ∃ B : ℝ, 0 ≤ B ∧ ∀ x ∈ Icc (-1 : ℝ) 2, ‖(riemannZeta (x + T * I))⁻¹‖ ≤ B := by
  have hc : ContinuousOn (fun x : ℝ => (riemannZeta (x + T * I))⁻¹) (Icc (-1) 2) := by
    intro x hx
    have h1 : (x : ℂ) + T * I ≠ 1 := fun h => by
      have := congrArg Complex.im h; simp at this; exact hT this
    refine ContinuousAt.continuousWithinAt ?_
    exact ContinuousAt.comp (f := fun x : ℝ => (x : ℂ) + T * I) (x := x)
      ((differentiableAt_riemannZeta h1).continuousAt.inv₀ (hgood x hx))
      (by fun_prop : Continuous fun x : ℝ => (x : ℂ) + T * I).continuousAt
  obtain ⟨B, hB⟩ := isCompact_Icc.exists_bound_of_continuousOn hc
  exact ⟨max B 0, le_max_right _ _, fun x hx => (hB x hx).trans (le_max_left _ _)⟩

theorem zeta_ne_zero_neg_one_line {y : ℝ} (hy : y ≠ 0) : riemannZeta (((-1 : ℝ) : ℂ) + y * I) ≠ 0 := by
  refine zeta_ne_zero_left (by simp) fun k h => ?_
  have := congrArg Complex.im h; simp [tz] at this; exact hy this

theorem line_ne_one {σ y : ℝ} (hy : y ≠ 0) : (σ : ℂ) + y * I ≠ 1 := fun h => by
  have := congrArg Complex.im h; simp at this; exact hy this

/-- **K is analytic on C minus {1}.** -/
theorem differentiableAt_Kf (hH : HBound H M β D N) {T₂ : ℝ} (hZ : HZeros H T₂) (hT1 : 1 ≤ T₂)
    {L : ℝ} (hL : 0 ≤ L) {s₀ : ℂ} (hs₀ : s₀ ≠ 1) : DifferentiableAt ℂ (Kf H L T₂) s₀ := by
  have hT0 : T₂ ≠ 0 := by linarith
  obtain ⟨B₁, hB₁0, hB₁⟩ := exists_bound_bottom hT0 hZ.good
  obtain ⟨B₂, hB₂0, hB₂⟩ := exists_bound_bottom (T := -T₂) (by linarith) (by
    intro σ hσ; have := hZ.good_neg σ hσ; simpa using this)
  have hHc : Continuous H := hH.diff.continuous
  have hBl : 0 ≤ 60 * π * Real.sqrt (2 * π) := by positivity
  -- the six pieces
  have hP1 := differentiableAt_piece (H := H) hL (S := Ioc (-1) 2) measurableSet_Ioc
    (γ := fun x : ℝ => (x : ℂ) + T₂ * I) (by fun_prop) hB₁0 (fun x _ => line_ne_one hT0)
    (fun x hx => hZ.good x ⟨hx.1.le, hx.2⟩) (fun x hx => hB₁ x ⟨hx.1.le, hx.2⟩)
    (fun x hx => by simp; exact hx.2) hHc
    ((hHc.comp (by fun_prop : Continuous fun x : ℝ => (x : ℂ) + T₂ * I)).norm.continuousOn.integrableOn_Icc.mono_set
      Ioc_subset_Icc_self) hs₀
  have hP4 := differentiableAt_piece (H := H) hL (S := Ioc (-1) 2) measurableSet_Ioc
    (γ := fun x : ℝ => (x : ℂ) + ((-T₂ : ℝ) : ℂ) * I) (by fun_prop) hB₂0 (fun x _ => line_ne_one (by linarith))
    (fun x hx => hZ.good_neg x ⟨hx.1.le, hx.2⟩) (fun x hx => by
      have := hB₂ x ⟨hx.1.le, hx.2⟩; simpa using this)
    (fun x hx => by simp; exact hx.2) hHc
    ((hHc.comp (by fun_prop : Continuous fun x : ℝ => (x : ℂ) + ((-T₂ : ℝ) : ℂ) * I)).norm.continuousOn.integrableOn_Icc.mono_set
      Ioc_subset_Icc_self) hs₀
  have hP2 := differentiableAt_piece (H := H) hL (S := Ioi T₂) measurableSet_Ioi
    (γ := fun y : ℝ => ((2 : ℝ) : ℂ) + y * I) (by fun_prop) (by norm_num : (0:ℝ) ≤ 3)
    (fun y hy => line_ne_one (by have : T₂ < y := hy; linarith))
    (fun y _ => riemannZeta_ne_zero_of_one_lt_re (by simp))
    (fun y _ => norm_inv_zeta_le_three (by simp)) (fun y _ => by simp) hHc
    (integrable_norm_H_line hH 2).integrableOn hs₀
  have hP3 := differentiableAt_piece (H := H) hL (S := Ioi T₂) measurableSet_Ioi
    (γ := fun y : ℝ => ((-1 : ℝ) : ℂ) + y * I) (by fun_prop) hBl
    (fun y hy => line_ne_one (by have : T₂ < y := hy; linarith))
    (fun y hy => zeta_ne_zero_neg_one_line (by have : T₂ < y := hy; linarith))
    (fun y hy => norm_inv_zeta_neg_one_line (by have : T₂ < y := hy; rw [abs_of_pos (by linarith)]; linarith))
    (fun y _ => by simp; norm_num) hHc (integrable_norm_H_line hH (-1)).integrableOn hs₀
  have hP5 := differentiableAt_piece (H := H) hL (S := Iic (-T₂)) measurableSet_Iic
    (γ := fun y : ℝ => ((2 : ℝ) : ℂ) + y * I) (by fun_prop) (by norm_num : (0:ℝ) ≤ 3)
    (fun y hy => line_ne_one (by have : y ≤ -T₂ := hy; linarith))
    (fun y _ => riemannZeta_ne_zero_of_one_lt_re (by simp))
    (fun y _ => norm_inv_zeta_le_three (by simp)) (fun y _ => by simp) hHc
    (integrable_norm_H_line hH 2).integrableOn hs₀
  have hP6 := differentiableAt_piece (H := H) hL (S := Iic (-T₂)) measurableSet_Iic
    (γ := fun y : ℝ => ((-1 : ℝ) : ℂ) + y * I) (by fun_prop) hBl
    (fun y hy => line_ne_one (by have : y ≤ -T₂ := hy; linarith))
    (fun y hy => zeta_ne_zero_neg_one_line (by have : y ≤ -T₂ := hy; linarith))
    (fun y hy => norm_inv_zeta_neg_one_line (by have : y ≤ -T₂ := hy; rw [abs_of_neg (by linarith)]; linarith))
    (fun y _ => by simp; norm_num) hHc (integrable_norm_H_line hH (-1)).integrableOn hs₀
  have e : Kf H L T₂ = fun s => (2 * π * I)⁻¹ *
      (((∫ x in Ioc (-1 : ℝ) 2, kint H L s ((x : ℂ) + T₂ * I)) +
        I * (∫ y in Ioi T₂, kint H L s (((2 : ℝ) : ℂ) + y * I)) -
        I * ∫ y in Ioi T₂, kint H L s (((-1 : ℝ) : ℂ) + y * I)) -
      ((∫ x in Ioc (-1 : ℝ) 2, kint H L s ((x : ℂ) + ((-T₂ : ℝ) : ℂ) * I)) -
        I * (∫ y in Iic (-T₂), kint H L s (((2 : ℝ) : ℂ) + y * I)) +
        I * ∫ y in Iic (-T₂), kint H L s (((-1 : ℝ) : ℂ) + y * I))) := by
    funext s
    simp only [Kf, UpperU, LowerU, intervalIntegral.integral_of_le (by norm_num : (-1 : ℝ) ≤ 2)]
  rw [e]
  exact (((hP1.add (hP2.const_mul I)).sub (hP3.const_mul I)).sub
    ((hP4.sub (hP5.const_mul I)).add (hP6.const_mul I))).const_mul _

end Riemann.Deep.LeanNBRate

