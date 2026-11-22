import Riemann.Deep.LeanNBRate.KJ5

/-!
# R_tr on Re s ≥ 1/2 from a bound on H at the trivial zeros
-/

open Complex Filter Topology Set Function
open scoped Real Nat

noncomputable section

namespace Riemann.Deep.LeanNBRate

theorem norm_rk_le' {H : ℂ → ℂ} {M' β' : ℝ}
    (hb : ∀ j : ℕ, ‖H (tz j)‖ ≤ M' * Real.exp (β' * (2 * j + 5 / 2))) (j : ℕ) :
    ‖rk H j‖ ≤ M' * Real.exp (β' * (2 * j + 5 / 2)) * (6 * (2 * π) ^ (2 * j + 2) / (2 * j + 2).factorial) := by
  have hd := norm_deriv_zeta_neg_two_ge j
  have hfpos : (0 : ℝ) < ((2 * j + 2).factorial : ℝ) := by exact_mod_cast Nat.factorial_pos _
  have hdpos : (0 : ℝ) < ((2 * j + 2).factorial : ℝ) / (6 * (2 * π) ^ (2 * j + 2)) := by positivity
  have h1 := hb j
  unfold tz at h1
  rw [rk, norm_div, div_le_iff₀ (lt_of_lt_of_le hdpos hd)]
  calc ‖H (-2 * ((j : ℂ) + 1))‖ ≤ M' * Real.exp (β' * (2 * j + 5 / 2)) := h1
    _ = M' * Real.exp (β' * (2 * j + 5 / 2)) * (6 * (2 * π) ^ (2 * j + 2) / (2 * j + 2).factorial) *
        (((2 * j + 2).factorial : ℝ) / (6 * (2 * π) ^ (2 * j + 2))) := by field_simp
    _ ≤ _ := by
        have hM : 0 ≤ M' * Real.exp (β' * (2 * j + 5 / 2)) * (6 * (2 * π) ^ (2 * j + 2) / (2 * j + 2).factorial) := by
          have h0 : 0 ≤ M' * Real.exp (β' * (2 * j + 5 / 2)) := (norm_nonneg _).trans h1
          positivity
        exact mul_le_mul_of_nonneg_left hd hM

theorem M'_nonneg {H : ℂ → ℂ} {M' β' : ℝ}
    (hb : ∀ j : ℕ, ‖H (tz j)‖ ≤ M' * Real.exp (β' * (2 * j + 5 / 2))) : 0 ≤ M' := by
  by_contra h
  push Not at h
  have h1 := hb 0
  have h2 : M' * Real.exp (β' * (2 * ((0 : ℕ) : ℝ) + 5 / 2)) < 0 := mul_neg_of_neg_of_pos h (Real.exp_pos _)
  linarith [norm_nonneg (H (tz 0))]

theorem norm_Ares_le' {H : ℂ → ℂ} {M' β' L : ℝ} (hL : 0 ≤ L) (hβ : β' ≤ L)
    (hb : ∀ j : ℕ, ‖H (tz j)‖ ≤ M' * Real.exp (β' * (2 * j + 5 / 2))) {s : ℂ} (hs : 1 / 2 ≤ s.re) (j : ℕ) :
    ‖Ares H L s j‖ ≤ M' * Real.exp (5 / 2 * (β' - L)) * (12 / 5 * ((2 * π) ^ (2 * j + 2) / (2 * j + 2).factorial)) := by
  have hj0 : (0 : ℝ) ≤ j := Nat.cast_nonneg j
  have htz : (tz j).re = -2 * j - 2 := tz_re j
  have hsz : 5 / 2 ≤ ‖s - tz j‖ := by
    have h1 := Complex.abs_re_le_norm (s - tz j)
    rw [sub_re, htz] at h1
    have h2 : (0 : ℝ) < s.re - (-2 * j - 2) := by linarith
    rw [abs_of_pos h2] at h1
    linarith
  have hk : ‖kap L (s - tz j)‖ ≤ Real.exp (-(L * (2 * j + 5 / 2))) := by
    rw [norm_kap, Real.exp_le_exp, sub_re, htz]; nlinarith
  have hrk := norm_rk_le' hb j
  have hM' := M'_nonneg hb
  have hfpos : (0 : ℝ) < ((2 * j + 2).factorial : ℝ) := by exact_mod_cast Nat.factorial_pos _
  rw [Ares, norm_div, norm_mul]
  calc ‖rk H j‖ * ‖kap L (s - tz j)‖ / ‖s - tz j‖ ≤ ‖rk H j‖ * ‖kap L (s - tz j)‖ / (5 / 2) :=
        div_le_div_of_nonneg_left (by positivity) (by norm_num) hsz
    _ ≤ (M' * Real.exp (β' * (2 * j + 5 / 2)) * (6 * (2 * π) ^ (2 * j + 2) / (2 * j + 2).factorial)) *
          Real.exp (-(L * (2 * j + 5 / 2))) / (5 / 2) := by gcongr
    _ = M' * Real.exp ((β' - L) * (2 * j + 5 / 2)) *
          (12 / 5 * ((2 * π) ^ (2 * j + 2) / (2 * j + 2).factorial)) := by
        rw [show (β' - L) * (2 * j + 5 / 2) = β' * (2 * j + 5 / 2) + -(L * (2 * j + 5 / 2)) by ring,
          Real.exp_add]
        ring
    _ ≤ _ := by
        have h3 : Real.exp ((β' - L) * (2 * j + 5 / 2)) ≤ Real.exp (5 / 2 * (β' - L)) :=
          Real.exp_le_exp.mpr (by nlinarith)
        have h4 : 0 ≤ 12 / 5 * ((2 * π) ^ (2 * j + 2) / ((2 * j + 2).factorial : ℝ)) := by positivity
        exact mul_le_mul_of_nonneg_right (mul_le_mul_of_nonneg_left h3 hM') h4

theorem summable_even_exp (x : ℝ) :
    Summable fun j : ℕ => x ^ (2 * j + 2) / ((2 * j + 2).factorial : ℝ) := by
  have hinj : Injective fun j : ℕ => 2 * j + 2 := fun a b h => by simp at h; omega
  exact (Real.summable_pow_div_factorial x).comp_injective hinj

theorem tsum_even_exp_le {x : ℝ} (hx : 0 ≤ x) :
    ∑' j : ℕ, x ^ (2 * j + 2) / ((2 * j + 2).factorial : ℝ) ≤ Real.exp x := by
  have hinj : Injective fun j : ℕ => 2 * j + 2 := fun a b h => by simp at h; omega
  calc ∑' j : ℕ, x ^ (2 * j + 2) / ((2 * j + 2).factorial : ℝ) ≤ ∑' n : ℕ, x ^ n / (n ! : ℝ) :=
        Summable.tsum_le_tsum_of_inj (fun j : ℕ => 2 * j + 2) hinj (fun c _ => by positivity)
          (fun i => le_rfl) (summable_even_exp x) (Real.summable_pow_div_factorial x)
    _ ≤ Real.exp x := Real.tsum_le_of_sum_range_le (fun n => by positivity)
        (fun n => Real.sum_le_exp_of_nonneg hx n)

/-- **R_tr on Re s ≥ 1/2.** -/
theorem norm_Rtr_le {H : ℂ → ℂ} {M' β' L : ℝ} (hL : 0 ≤ L) (hβ : β' ≤ L)
    (hb : ∀ j : ℕ, ‖H (tz j)‖ ≤ M' * Real.exp (β' * (2 * j + 5 / 2))) {s : ℂ} (hs : 1 / 2 ≤ s.re) :
    ‖Rtr H L s‖ ≤ M' * Real.exp (5 / 2 * (β' - L)) * (12 / 5 * Real.exp (2 * π)) := by
  rw [Rtr_eq]
  have hM' := M'_nonneg hb
  have hbd := fun j => norm_Ares_le' hL hβ hb hs j
  have hS : Summable fun j : ℕ => M' * Real.exp (5 / 2 * (β' - L)) *
      (12 / 5 * ((2 * π) ^ (2 * j + 2) / ((2 * j + 2).factorial : ℝ))) :=
    ((summable_even_exp (2 * π)).mul_left (12 / 5)).mul_left _
  have hS' : Summable fun j : ℕ => ‖Ares H L s j‖ := hS.of_nonneg_of_le (fun j => norm_nonneg _) hbd
  calc ‖∑' j, Ares H L s j‖ ≤ ∑' j, ‖Ares H L s j‖ := norm_tsum_le_tsum_norm hS'
    _ ≤ ∑' j : ℕ, M' * Real.exp (5 / 2 * (β' - L)) *
          (12 / 5 * ((2 * π) ^ (2 * j + 2) / ((2 * j + 2).factorial : ℝ))) :=
        Summable.tsum_le_tsum hbd hS' hS
    _ = M' * Real.exp (5 / 2 * (β' - L)) *
          (12 / 5 * ∑' j : ℕ, (2 * π) ^ (2 * j + 2) / ((2 * j + 2).factorial : ℝ)) := by
        rw [tsum_mul_left, tsum_mul_left]
    _ ≤ _ := by
        have h := tsum_even_exp_le (x := 2 * π) (by positivity)
        have h0 : 0 ≤ M' * Real.exp (5 / 2 * (β' - L)) := by positivity
        exact mul_le_mul_of_nonneg_left (mul_le_mul_of_nonneg_left h (by norm_num)) h0

end Riemann.Deep.LeanNBRate

