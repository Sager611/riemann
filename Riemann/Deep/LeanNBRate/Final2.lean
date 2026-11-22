import Riemann.Deep.LeanNBRate.Final1

/-!
# Asymptotic lemmas and the explicit error constants for the final estimate
-/

open MeasureTheory Filter Topology Set Complex
open scoped Real

noncomputable section

namespace Riemann.Deep.LeanNBRate

theorem tendsto_L_div_ell : Tendsto (fun L : ℝ => L / (L - Real.sqrt L)) atTop (𝓝 1) := by
  have h1 : Tendsto (fun L : ℝ => (Real.sqrt L)⁻¹) atTop (𝓝 0) :=
    tendsto_inv_atTop_zero.comp Real.tendsto_sqrt_atTop
  have h2 : Tendsto (fun L : ℝ => (1 - (Real.sqrt L)⁻¹)⁻¹) atTop (𝓝 1) := by
    have := (tendsto_const_nhds (x := (1 : ℝ))).sub h1
    have := this.inv₀ (by norm_num : (1 : ℝ) - 0 ≠ 0)
    simpa using this
  refine h2.congr' ?_
  filter_upwards [eventually_gt_atTop (1 : ℝ)] with L hL
  have hs : 1 < Real.sqrt L := by
    rw [show (1 : ℝ) = Real.sqrt 1 by simp]; exact Real.sqrt_lt_sqrt (by norm_num) hL
  have hsq : Real.sqrt L ^ 2 = L := Real.sq_sqrt (by linarith)
  obtain ⟨x, hx⟩ : ∃ x, x = Real.sqrt L := ⟨_, rfl⟩
  rw [← hx] at hs hsq ⊢
  rw [← hsq]
  have hne : x - 1 ≠ 0 := by linarith
  have hs0 : x ≠ 0 := by linarith
  have hne' : x ^ 2 - x ≠ 0 := by
    have : x ^ 2 - x = x * (x - 1) := by ring
    rw [this]; exact mul_ne_zero hs0 hne
  field_simp

theorem tendsto_L_exp_neg4 : Tendsto (fun L : ℝ => L * Real.exp (-(4 * L))) atTop (𝓝 0) := by
  have h := Real.tendsto_pow_mul_exp_neg_atTop_nhds_zero 1
  simp only [pow_one] at h
  refine squeeze_zero' ?_ ?_ h
  · filter_upwards [eventually_ge_atTop (0 : ℝ)] with L hL
    positivity
  · filter_upwards [eventually_ge_atTop (0 : ℝ)] with L hL
    exact mul_le_mul_of_nonneg_left (Real.exp_le_exp.mpr (by linarith)) hL

theorem tendsto_L3_exp_sqrt : Tendsto (fun L : ℝ => L ^ 3 * Real.exp (-(5 * Real.sqrt L))) atTop (𝓝 0) := by
  have h := (Real.tendsto_pow_mul_exp_neg_atTop_nhds_zero 6).comp Real.tendsto_sqrt_atTop
  refine squeeze_zero' ?_ ?_ h
  · filter_upwards [eventually_ge_atTop (0 : ℝ)] with L hL
    positivity
  · filter_upwards [eventually_ge_atTop (0 : ℝ)] with L hL
    have hs := Real.sqrt_nonneg L
    have e : L ^ 3 = Real.sqrt L ^ 6 := by
      rw [show (6 : ℕ) = 2 * 3 by norm_num, pow_mul, Real.sq_sqrt hL]
    rw [e]
    simp only [Function.comp]
    exact mul_le_mul_of_nonneg_left (Real.exp_le_exp.mpr (by linarith)) (by positivity)

theorem pow_four_div_le {L : ℝ} {m p : ℕ} (hL : 2 ≤ L) (h5 : 4 * Real.exp 5 ≤ L ^ 3) (hm : L ≤ m)
    (hp : (p : ℝ) ≤ L) : (4 / L ^ 3) ^ m ≤ Real.exp (-(5 * (L - p))) / L ^ p := by
  have hL0 : 0 < L := by linarith
  have hq0 : 0 ≤ 4 / L ^ 3 := by positivity
  have hq5 : 4 / L ^ 3 ≤ Real.exp (-5) := by
    rw [div_le_iff₀ (by positivity), Real.exp_neg]
    have := Real.exp_pos 5
    rw [← div_eq_inv_mul, le_div_iff₀ this]
    linarith
  have hqL : 4 / L ^ 3 ≤ 1 / L := by
    rw [div_le_div_iff₀ (by positivity) hL0]
    have h4 : 4 ≤ L ^ 2 := by nlinarith
    have h6 : 4 * L ≤ L ^ 2 * L := by nlinarith
    linarith [show L ^ 3 = L ^ 2 * L by ring]
  have hpm : p ≤ m := by exact_mod_cast hp.trans hm
  calc (4 / L ^ 3) ^ m = (4 / L ^ 3) ^ (m - p) * (4 / L ^ 3) ^ p := by
        rw [← pow_add, Nat.sub_add_cancel hpm]
    _ ≤ Real.exp (-5) ^ (m - p) * (1 / L) ^ p :=
        mul_le_mul (pow_le_pow_left₀ hq0 hq5 _) (pow_le_pow_left₀ hq0 hqL _) (by positivity) (by positivity)
    _ = Real.exp (-(5 * ((m : ℝ) - p))) / L ^ p := by
        rw [← Real.exp_nat_mul, Nat.cast_sub hpm, one_div_pow, mul_one_div]
        congr 2; ring
    _ ≤ Real.exp (-(5 * (L - p))) / L ^ p := by
        gcongr

/-- The exponent of the 1/ζ bound, rounded up. -/
def kA (A : ℝ) : ℕ := ⌈max A 0⌉₊

/-- The polynomial degree collected in Q. -/
def pE (A : ℝ) : ℕ := 6 * kA A + 13

def KB (C A : ℝ) : ℝ := |C| * 2 ^ kA A + 723

def KQ (C A : ℝ) : ℝ := (2 * π)⁻¹ * (8 + 200 * KB C A) * 80 * (1 + 2 * Real.sqrt (zeroSum + 1))

/-- The explicit bound for Q. -/
def Qb (C A L : ℝ) : ℝ := KQ C A * Real.exp (5 * pE A) * Real.exp (-(2 * L))

def KR : ℝ := (1 + Real.sqrt (zeroSum + 1)) * Real.exp 5 * (12 / 5 * Real.exp (2 * π))

/-- The explicit bound for R_tr. -/
def Rb (L : ℝ) : ℝ := KR * L * Real.exp (-(5 / 2 * Real.sqrt L))

theorem tendsto_L_Qb_sq (C A : ℝ) : Tendsto (fun L : ℝ => L * Qb C A L ^ 2) atTop (𝓝 0) := by
  have h := tendsto_L_exp_neg4.const_mul ((KQ C A * Real.exp (5 * pE A)) ^ 2)
  rw [mul_zero] at h
  refine h.congr' (Filter.Eventually.of_forall fun L => ?_)
  unfold Qb
  rw [show -(4 * L) = -(2 * L) + -(2 * L) by ring, Real.exp_add]
  ring

theorem tendsto_L_Rb_sq : Tendsto (fun L : ℝ => L * Rb L ^ 2) atTop (𝓝 0) := by
  have h := tendsto_L3_exp_sqrt.const_mul (KR ^ 2)
  rw [mul_zero] at h
  refine h.congr' (Filter.Eventually.of_forall fun L => ?_)
  unfold Rb
  rw [show -(5 * Real.sqrt L) = -(5 / 2 * Real.sqrt L) + -(5 / 2 * Real.sqrt L) by ring, Real.exp_add]
  ring

/-- The error term tends to zero. -/
theorem tendsto_err (C A θ : ℝ) :
    Tendsto (fun L : ℝ => 3 * (1 + θ⁻¹) * (2 * π)⁻¹ * (8 * π / L + (4 * π + 4 * Zz) * (L * Qb C A L ^ 2) +
      8 * Zz * (L * Rb L ^ 2))) atTop (𝓝 0) := by
  have h1 : Tendsto (fun L : ℝ => 8 * π / L) atTop (𝓝 0) :=
    tendsto_const_nhds.div_atTop tendsto_id
  have h := ((h1.add ((tendsto_L_Qb_sq C A).const_mul (4 * π + 4 * Zz))).add
    (tendsto_L_Rb_sq.const_mul (8 * Zz))).const_mul (3 * (1 + θ⁻¹) * (2 * π)⁻¹)
  simpa using h

end Riemann.Deep.LeanNBRate

