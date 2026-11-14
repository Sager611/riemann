/-
E1Lean.Kappa — the detector side of E1.md §4: e^{8.10242} ≥ c_hi, so κ₁ = 1 − R⁻¹ log c ≥ 37379/442500 for
c ∈ [c_lo, c_hi] (the H4 interval), and then h(1/2; κ₁) > 0.0055, k*(1/2) < κ₁.
-/
import Riemann.Deep.SantiAudit.Lean.Numeric

namespace E1Lean

variable (F : Fns)

/-- c_hi = 3302.286923204427565640 ≤ e^{8 + 2·0.05121} (Real.exp_one_gt_d9, quadratic_le_exp_of_nonneg). -/
theorem exp_810242 : (3302286923204427565640 / 1000000000000000000 : ℝ) ≤ F.exp (810242 / 100000) := by
  have e1 := F.exp_one_gt_d9
  have h8 : F.exp 8 = F.exp 1 ^ 8 := by exact_mod_cast exp_nat F 8
  have hq := F.quadratic_le_exp_of_nonneg (5121 / 100000) (by norm_num)
  have hsplit : F.exp (810242 / 100000) = F.exp 8 * (F.exp (5121 / 100000) * F.exp (5121 / 100000)) := by
    rw [← F.exp_add, ← F.exp_add]; norm_num
  rw [hsplit, h8]
  have hpos : (0 : ℝ) ≤ 1 + 5121 / 100000 + (5121 / 100000) ^ 2 / 2 := by norm_num
  have hq2 := mul_le_mul hq hq hpos (le_trans hpos hq)
  have he8 := pow_le_pow_left₀ (by norm_num) e1.le 8
  calc (3302286923204427565640 / 1000000000000000000 : ℝ)
        ≤ (2.7182818283 : ℝ) ^ 8 * ((1 + 5121 / 100000 + (5121 / 100000) ^ 2 / 2) *
          (1 + 5121 / 100000 + (5121 / 100000) ^ 2 / 2)) := by norm_num
    _ ≤ F.exp 1 ^ 8 * (F.exp (5121 / 100000) * F.exp (5121 / 100000)) :=
        mul_le_mul he8 hq2 (by norm_num) (pow_nonneg (F.exp_pos 1).le 8)

/-- For c in the H4 interval: 37379/442500 ≤ κ = 1 − (177/20)⁻¹ log c < 1. -/
theorem kappa_bounds {c : ℝ} (hlo : (3302286923204427565638 / 1000000000000000000 : ℝ) ≤ c)
    (hhi : c ≤ 3302286923204427565640 / 1000000000000000000) :
    (37379 / 442500 : ℝ) ≤ kappa F c (177 / 20) ∧ kappa F c (177 / 20) < 1 := by
  have hc0 : 0 < c := lt_of_lt_of_le (by norm_num) hlo
  have hlog : F.log c ≤ 810242 / 100000 := (F.log_le_iff_le_exp c _ hc0).2 (hhi.trans (exp_810242 F))
  have hlog0 : 0 < F.log c := by
    by_contra hneg
    rw [not_lt] at hneg
    have h1 := (F.log_le_iff_le_exp c 0 hc0).1 hneg
    rw [F.exp_zero] at h1
    linarith
  unfold kappa
  constructor
  · have h1 : (177 / 20 : ℝ)⁻¹ * F.log c ≤ (177 / 20 : ℝ)⁻¹ * (810242 / 100000) :=
      mul_le_mul_of_nonneg_left hlog (by norm_num)
    have h2 : (37379 / 442500 : ℝ) = 1 - (177 / 20 : ℝ)⁻¹ * (810242 / 100000) := by norm_num
    linarith
  · have h1 : 0 < (177 / 20 : ℝ)⁻¹ * F.log c := mul_pos (by norm_num) hlog0
    linarith

/-- E1.md §4: h(1/2; κ) > 0.0055 for every κ ∈ [37379/442500, 1). -/
theorem h_half_gt {κ : ℝ} (hk : (37379 / 442500 : ℝ) ≤ κ) (hk1 : κ < 1) : (0.0055 : ℝ) < h F (1 / 2) κ := by
  have hcw := cW_half_ge F
  have hc1 := cW_lt_one F (by norm_num : (0 : ℝ) < 1 / 2) (by norm_num)
  have m1 := hc_mono_k F hc1 hk hk1
  have m2 := hc_mono_c F (by norm_num : (37379 / 442500 : ℝ) < 1) hcw hc1
  have base : (0.0055 : ℝ) < hc F (-16596708521464472031 / 100000000000000000000 : ℝ) (37379 / 442500) := by
    unfold hc
    have hA : (0 : ℝ) ≤ (1 - 37379 / 442500) * (9 - 37379 / 442500 - 8 * (-16596708521464472031 / 100000000000000000000 : ℝ)) := by norm_num
    have := sqrt_lt_of_lt_sq F hA (B := 3 + 37379 / 442500 - 0.022) (by norm_num) (by norm_num)
    linarith
  unfold h
  linarith

/-- E1.md §4: k*(1/2) < κ for every κ ∈ [37379/442500, 1) (so κ₁ − k*(1/2) > 0). -/
theorem kstar_half_lt {κ : ℝ} (hk : (37379 / 442500 : ℝ) ≤ κ) (hk1 : κ < 1) : kstar F (1 / 2) < κ :=
  (h_pos_iff F (by norm_num) (by norm_num) hk1).1 (lt_trans (by norm_num) (h_half_gt F hk hk1))

end E1Lean
