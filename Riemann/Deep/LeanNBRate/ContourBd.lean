import Riemann.Deep.LeanNBRate.KJ5

/-!
# Norm bounds for the half-strip boundary integrals (contour mass of H)
-/

open Complex Filter Topology Set MeasureTheory intervalIntegral
open scoped Real

noncomputable section

namespace Riemann.Deep.LeanNBRate

variable {H : ℂ → ℂ} {M β D : ℝ} {N : ℕ}

theorem hasDerivAt_neg_sq_div {D x : ℝ} (hx : D + x ≠ 0) :
    HasDerivAt (fun y : ℝ => -(D ^ 2) / (D + y)) (D ^ 2 / (D + x) ^ 2) x := by
  have h := (((hasDerivAt_id x).const_add D).inv hx).const_mul (-(D ^ 2))
  convert h using 1
  · funext y; simp [div_eq_mul_inv]
  · simp only [id]; ring

theorem integral_Ioi_sq_div {D T : ℝ} (hD : 0 < D) (hT : 0 ≤ T) :
    IntegrableOn (fun y : ℝ => D ^ 2 / (D + y) ^ 2) (Ioi T) ∧
    ∫ y in Ioi T, D ^ 2 / (D + y) ^ 2 = D ^ 2 / (D + T) := by
  have hlim : Tendsto (fun y : ℝ => -(D ^ 2) / (D + y)) atTop (𝓝 0) :=
    tendsto_const_nhds.div_atTop (tendsto_atTop_add_const_left _ D tendsto_id)
  have hT' : D + T ≠ 0 := by linarith
  have hcont : ContinuousWithinAt (fun y : ℝ => -(D ^ 2) / (D + y)) (Ici T) T :=
    (hasDerivAt_neg_sq_div hT').continuousAt.continuousWithinAt
  have hder : ∀ x ∈ Ioi T, HasDerivAt (fun y : ℝ => -(D ^ 2) / (D + y)) (D ^ 2 / (D + x) ^ 2) x := by
    intro x hx
    have : T < x := hx
    exact hasDerivAt_neg_sq_div (by linarith)
  have hpos : ∀ x ∈ Ioi T, 0 ≤ D ^ 2 / (D + x) ^ 2 := fun x _ => by positivity
  refine ⟨integrableOn_Ioi_deriv_of_nonneg hcont hder hpos hlim, ?_⟩
  rw [integral_Ioi_of_hasDerivAt_of_nonneg hcont hder hpos hlim]
  ring

/-- Half-line bound: if φ ≤ K (D/(D+y))^N on y > T, then the integral is ≤ K (D/(D+T))^(N-2) D. -/
theorem setIntegral_Ioi_le_weight {φ : ℝ → ℝ} {K D T : ℝ} {N : ℕ} (hD : 0 < D) (hT : 0 ≤ T) (hK : 0 ≤ K)
    (hN : 2 ≤ N) (hφ : IntegrableOn φ (Ioi T)) (hb : ∀ y ∈ Ioi T, φ y ≤ K * (D / (D + y)) ^ N) :
    ∫ y in Ioi T, φ y ≤ K * (D / (D + T)) ^ (N - 2) * D := by
  obtain ⟨hint, hval⟩ := integral_Ioi_sq_div hD hT
  calc ∫ y in Ioi T, φ y ≤ ∫ y in Ioi T, K * (D / (D + T)) ^ (N - 2) * (D ^ 2 / (D + y) ^ 2) := by
        refine setIntegral_mono_on hφ (hint.const_mul _) measurableSet_Ioi fun y hy => ?_
        refine (hb y hy).trans ?_
        have hy' : T < y := hy
        have h1 : D / (D + y) ≤ D / (D + T) := div_le_div_of_nonneg_left hD.le (by linarith) (by linarith)
        have h0 : 0 ≤ D / (D + y) := div_nonneg hD.le (by linarith)
        have e : (D / (D + y)) ^ N = (D / (D + y)) ^ (N - 2) * (D ^ 2 / (D + y) ^ 2) := by
          rw [← div_pow, ← pow_add]; congr 1; omega
        rw [e]
        have h2 : (D / (D + y)) ^ (N - 2) ≤ (D / (D + T)) ^ (N - 2) := pow_le_pow_left₀ h0 h1 _
        have hb0 : 0 ≤ D ^ 2 / (D + y) ^ 2 := by positivity
        calc K * ((D / (D + y)) ^ (N - 2) * (D ^ 2 / (D + y) ^ 2))
            ≤ K * ((D / (D + T)) ^ (N - 2) * (D ^ 2 / (D + y) ^ 2)) :=
              mul_le_mul_of_nonneg_left (mul_le_mul_of_nonneg_right h2 hb0) hK
          _ = _ := by ring
    _ = K * (D / (D + T)) ^ (N - 2) * (D ^ 2 / (D + T)) := by
        rw [MeasureTheory.integral_const_mul, hval]
    _ ≤ K * (D / (D + T)) ^ (N - 2) * D := by
        have h3 : D ^ 2 / (D + T) ≤ D := by
          rw [div_le_iff₀ (by linarith)]; nlinarith
        have h4 : 0 ≤ K * (D / (D + T)) ^ (N - 2) := by
          have : 0 ≤ D / (D + T) := div_nonneg hD.le (by linarith)
          positivity
        exact mul_le_mul_of_nonneg_left h3 h4

theorem norm_H_contour (hH : HBound H M β D N) {w : ℂ} (h1 : -1 ≤ w.re) :
    ‖H w‖ ≤ M * Real.exp (β * (3 / 2)) * (D / (D + |w.im|)) ^ N := by
  have hb := hH.bound w
  have hw := weight_mono hH.D_pos N (z := w)
  have hexp : Real.exp (β * max 0 (1 / 2 - w.re)) ≤ Real.exp (β * (3 / 2)) :=
    Real.exp_le_exp.mpr (mul_le_mul_of_nonneg_left (max_le (by norm_num) (by linarith)) hH.β_nonneg)
  have hM := hH.M_nonneg
  have hD := hH.D_pos
  have h0 : 0 ≤ (D / (D + ‖w - 1 / 2‖)) ^ N := by positivity
  exact hb.trans (mul_le_mul (mul_le_mul_of_nonneg_left hexp hM) hw h0 (by positivity))

/-- Contour mass of H on the upper half-strip boundary. -/
def ZU (H : ℂ → ℂ) (T : ℝ) : ℝ :=
  (∫ x in (-1 : ℝ)..2, ‖H (x + T * I)‖) + (∫ y in Ioi T, ‖H (((2 : ℝ) : ℂ) + y * I)‖) +
    ∫ y in Ioi T, ‖H (((-1 : ℝ) : ℂ) + y * I)‖

/-- Contour mass of H on the lower half-strip boundary. -/
def ZL (H : ℂ → ℂ) (T : ℝ) : ℝ :=
  (∫ x in (-1 : ℝ)..2, ‖H (x + ((-T : ℝ) : ℂ) * I)‖) + (∫ y in Iic (-T), ‖H (((2 : ℝ) : ℂ) + y * I)‖) +
    ∫ y in Iic (-T), ‖H (((-1 : ℝ) : ℂ) + y * I)‖

theorem norm_UpperU_le (hH : HBound H M β D N) {f : ℂ → ℂ} {C T : ℝ} (hC : 0 ≤ C)
    (hbot : ∀ x ∈ Icc (-1 : ℝ) 2, ‖f (x + T * I)‖ ≤ C * ‖H (x + T * I)‖)
    (hr : ∀ y ∈ Ioi T, ‖f (((2 : ℝ) : ℂ) + y * I)‖ ≤ C * ‖H (((2 : ℝ) : ℂ) + y * I)‖)
    (hl : ∀ y ∈ Ioi T, ‖f (((-1 : ℝ) : ℂ) + y * I)‖ ≤ C * ‖H (((-1 : ℝ) : ℂ) + y * I)‖) :
    ‖UpperU f (-1) 2 T‖ ≤ C * ZU H T := by
  have hcb : Continuous fun x : ℝ => ‖H (x + T * I)‖ := (hH.diff.continuous.comp (by fun_prop)).norm
  have hA : ‖∫ x in (-1 : ℝ)..2, f (x + T * I)‖ ≤ ∫ x in (-1 : ℝ)..2, C * ‖H (x + T * I)‖ :=
    intervalIntegral.norm_integral_le_of_norm_le (by norm_num)
      (Filter.Eventually.of_forall fun x hx => hbot x (Ioc_subset_Icc_self hx))
      ((continuous_const.mul hcb).intervalIntegrable _ _)
  have hB : ‖∫ y in Ioi T, f (((2 : ℝ) : ℂ) + y * I)‖ ≤ ∫ y in Ioi T, C * ‖H (((2 : ℝ) : ℂ) + y * I)‖ :=
    norm_integral_le_of_norm_le ((integrable_norm_H_line hH 2).integrableOn.const_mul C)
      ((ae_restrict_iff' measurableSet_Ioi).mpr (Filter.Eventually.of_forall hr))
  have hC' : ‖∫ y in Ioi T, f (((-1 : ℝ) : ℂ) + y * I)‖ ≤ ∫ y in Ioi T, C * ‖H (((-1 : ℝ) : ℂ) + y * I)‖ :=
    norm_integral_le_of_norm_le ((integrable_norm_H_line hH (-1)).integrableOn.const_mul C)
      ((ae_restrict_iff' measurableSet_Ioi).mpr (Filter.Eventually.of_forall hl))
  rw [intervalIntegral.integral_const_mul] at hA
  rw [MeasureTheory.integral_const_mul] at hB hC'
  unfold UpperU ZU
  have n1 := norm_sub_le ((∫ x in (-1 : ℝ)..2, f (x + T * I)) + I * (∫ y in Ioi T, f (((2 : ℝ) : ℂ) + y * I)))
    (I * ∫ y in Ioi T, f (((-1 : ℝ) : ℂ) + y * I))
  have n2 := norm_add_le (∫ x in (-1 : ℝ)..2, f (x + T * I)) (I * (∫ y in Ioi T, f (((2 : ℝ) : ℂ) + y * I)))
  simp only [norm_mul, Complex.norm_I, one_mul] at n1 n2
  linarith

theorem norm_LowerU_le (hH : HBound H M β D N) {f : ℂ → ℂ} {C T : ℝ} (hC : 0 ≤ C)
    (hbot : ∀ x ∈ Icc (-1 : ℝ) 2, ‖f (x + ((-T : ℝ) : ℂ) * I)‖ ≤ C * ‖H (x + ((-T : ℝ) : ℂ) * I)‖)
    (hr : ∀ y ∈ Iic (-T), ‖f (((2 : ℝ) : ℂ) + y * I)‖ ≤ C * ‖H (((2 : ℝ) : ℂ) + y * I)‖)
    (hl : ∀ y ∈ Iic (-T), ‖f (((-1 : ℝ) : ℂ) + y * I)‖ ≤ C * ‖H (((-1 : ℝ) : ℂ) + y * I)‖) :
    ‖LowerU f (-1) 2 T‖ ≤ C * ZL H T := by
  have hcb : Continuous fun x : ℝ => ‖H (x + ((-T : ℝ) : ℂ) * I)‖ := (hH.diff.continuous.comp (by fun_prop)).norm
  have hA : ‖∫ x in (-1 : ℝ)..2, f (x + ((-T : ℝ) : ℂ) * I)‖ ≤ ∫ x in (-1 : ℝ)..2, C * ‖H (x + ((-T : ℝ) : ℂ) * I)‖ :=
    intervalIntegral.norm_integral_le_of_norm_le (by norm_num)
      (Filter.Eventually.of_forall fun x hx => hbot x (Ioc_subset_Icc_self hx))
      ((continuous_const.mul hcb).intervalIntegrable _ _)
  have hB : ‖∫ y in Iic (-T), f (((2 : ℝ) : ℂ) + y * I)‖ ≤ ∫ y in Iic (-T), C * ‖H (((2 : ℝ) : ℂ) + y * I)‖ :=
    norm_integral_le_of_norm_le ((integrable_norm_H_line hH 2).integrableOn.const_mul C)
      ((ae_restrict_iff' measurableSet_Iic).mpr (Filter.Eventually.of_forall hr))
  have hC' : ‖∫ y in Iic (-T), f (((-1 : ℝ) : ℂ) + y * I)‖ ≤ ∫ y in Iic (-T), C * ‖H (((-1 : ℝ) : ℂ) + y * I)‖ :=
    norm_integral_le_of_norm_le ((integrable_norm_H_line hH (-1)).integrableOn.const_mul C)
      ((ae_restrict_iff' measurableSet_Iic).mpr (Filter.Eventually.of_forall hl))
  rw [intervalIntegral.integral_const_mul] at hA
  rw [MeasureTheory.integral_const_mul] at hB hC'
  unfold LowerU ZL
  have n1 := norm_add_le ((∫ x in (-1 : ℝ)..2, f (x + ((-T : ℝ) : ℂ) * I)) -
    I * (∫ y in Iic (-T), f (((2 : ℝ) : ℂ) + y * I))) (I * ∫ y in Iic (-T), f (((-1 : ℝ) : ℂ) + y * I))
  have n2 := norm_sub_le (∫ x in (-1 : ℝ)..2, f (x + ((-T : ℝ) : ℂ) * I))
    (I * (∫ y in Iic (-T), f (((2 : ℝ) : ℂ) + y * I)))
  simp only [norm_mul, Complex.norm_I, one_mul] at n1 n2
  linarith

theorem ZU_le (hH : HBound H M β D N) {T : ℝ} (hT : 0 ≤ T) :
    ZU H T ≤ M * Real.exp (β * (3 / 2)) * (D / (D + T)) ^ (N - 2) * (3 + 2 * D) := by
  have hD := hH.D_pos
  have hM := hH.M_nonneg
  have hN : 2 ≤ N := by linarith [hH.N_ge]
  obtain ⟨K, hK⟩ : ∃ K : ℝ, K = M * Real.exp (β * (3 / 2)) := ⟨_, rfl⟩
  have hK0 : 0 ≤ K := by rw [hK]; positivity
  have hq0 : 0 ≤ D / (D + T) := div_nonneg hD.le (by linarith)
  have hq1 : D / (D + T) ≤ 1 := (div_le_one (by linarith)).mpr (by linarith)
  have hcb : Continuous fun x : ℝ => ‖H (x + T * I)‖ := (hH.diff.continuous.comp (by fun_prop)).norm
  have hbot : ∫ x in (-1 : ℝ)..2, ‖H (x + T * I)‖ ≤ 3 * (K * (D / (D + T)) ^ (N - 2)) := by
    have h := intervalIntegral.integral_mono_on (by norm_num : (-1 : ℝ) ≤ 2) (hcb.intervalIntegrable _ _)
      (intervalIntegrable_const (μ := volume) (c := K * (D / (D + T)) ^ (N - 2)))
      (fun x hx => by
        have h1 := norm_H_contour hH (w := (x : ℂ) + T * I) (by simp; exact hx.1)
        have him : |((x : ℂ) + T * I).im| = T := by simp [abs_of_nonneg hT]
        rw [him, ← hK] at h1
        exact h1.trans (mul_le_mul_of_nonneg_left (pow_le_pow_of_le_one hq0 hq1 (Nat.sub_le N 2)) hK0))
    rw [intervalIntegral.integral_const, smul_eq_mul, show (2 : ℝ) - -1 = 3 by norm_num] at h
    exact h
  have hvert : ∀ σ : ℝ, -1 ≤ σ → ∫ y in Ioi T, ‖H ((σ : ℂ) + y * I)‖ ≤ K * (D / (D + T)) ^ (N - 2) * D := by
    intro σ hσ
    refine setIntegral_Ioi_le_weight hD hT hK0 hN (integrable_norm_H_line hH σ).integrableOn fun y hy => ?_
    have hy' : T < y := hy
    have h1 := norm_H_contour hH (w := (σ : ℂ) + y * I) (by simp; exact hσ)
    have him : |((σ : ℂ) + y * I).im| = y := by simp [abs_of_nonneg (by linarith : (0 : ℝ) ≤ y)]
    rwa [him, ← hK] at h1
  have h2 := hvert 2 (by norm_num)
  have h3 := hvert (-1) le_rfl
  unfold ZU
  rw [← hK]
  nlinarith

theorem ZL_le (hH : HBound H M β D N) {T : ℝ} (hT : 0 ≤ T) :
    ZL H T ≤ M * Real.exp (β * (3 / 2)) * (D / (D + T)) ^ (N - 2) * (3 + 2 * D) := by
  have hD := hH.D_pos
  have hM := hH.M_nonneg
  have hN : 2 ≤ N := by linarith [hH.N_ge]
  obtain ⟨K, hK⟩ : ∃ K : ℝ, K = M * Real.exp (β * (3 / 2)) := ⟨_, rfl⟩
  have hK0 : 0 ≤ K := by rw [hK]; positivity
  have hq0 : 0 ≤ D / (D + T) := div_nonneg hD.le (by linarith)
  have hq1 : D / (D + T) ≤ 1 := (div_le_one (by linarith)).mpr (by linarith)
  have hcb : Continuous fun x : ℝ => ‖H (x + ((-T : ℝ) : ℂ) * I)‖ := (hH.diff.continuous.comp (by fun_prop)).norm
  have hbot : ∫ x in (-1 : ℝ)..2, ‖H (x + ((-T : ℝ) : ℂ) * I)‖ ≤ 3 * (K * (D / (D + T)) ^ (N - 2)) := by
    have h := intervalIntegral.integral_mono_on (by norm_num : (-1 : ℝ) ≤ 2) (hcb.intervalIntegrable _ _)
      (intervalIntegrable_const (μ := volume) (c := K * (D / (D + T)) ^ (N - 2)))
      (fun x hx => by
        have h1 := norm_H_contour hH (w := (x : ℂ) + ((-T : ℝ) : ℂ) * I) (by simp; exact hx.1)
        have him : |((x : ℂ) + ((-T : ℝ) : ℂ) * I).im| = T := by simp [abs_of_nonneg hT]
        rw [him, ← hK] at h1
        exact h1.trans (mul_le_mul_of_nonneg_left (pow_le_pow_of_le_one hq0 hq1 (Nat.sub_le N 2)) hK0))
    rw [intervalIntegral.integral_const, smul_eq_mul, show (2 : ℝ) - -1 = 3 by norm_num] at h
    exact h
  have hvert : ∀ σ : ℝ, -1 ≤ σ → ∫ y in Iic (-T), ‖H ((σ : ℂ) + y * I)‖ ≤ K * (D / (D + T)) ^ (N - 2) * D := by
    intro σ hσ
    rw [← integral_comp_neg_Ioi]
    have hi : IntegrableOn (fun y : ℝ => ‖H ((σ : ℂ) + ((-y : ℝ) : ℂ) * I)‖) (Ioi T) :=
      ((integrable_norm_H_line hH σ).comp_neg).integrableOn
    refine setIntegral_Ioi_le_weight hD hT hK0 hN hi fun y hy => ?_
    have hy' : T < y := hy
    have h1 := norm_H_contour hH (w := (σ : ℂ) + ((-y : ℝ) : ℂ) * I) (by simp; exact hσ)
    have him : |((σ : ℂ) + ((-y : ℝ) : ℂ) * I).im| = y := by
      simp [abs_of_nonneg (by linarith : (0 : ℝ) ≤ y)]
    rwa [him, ← hK] at h1
  have h2 := hvert 2 (by norm_num)
  have h3 := hvert (-1) le_rfl
  unfold ZL
  rw [← hK]
  nlinarith

end Riemann.Deep.LeanNBRate

