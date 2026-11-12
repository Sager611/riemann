module

public import Riemann.Deep.Combos.Lean.Expo1pC3

@[expose] public section

/-! # Theorem 1′(a): class I exponent (RCApps2 l.50–55) -/

open Real

namespace RingClass

/-- class I: (2Y log²Y)^{2−2σ+2ε₀} ≤ (8/ε₀²) W^{e(σ)+3ε₀}, Y = W^{1/(3−2σ)}, W ≥ 1 -/
lemma classI_exp (W σ ε₀ : ℝ) (hW : 1 ≤ W) (h0 : 0.6 ≤ σ) (h1 : σ ≤ 1) (he : 0 < ε₀)
    (he1 : ε₀ ≤ 1 / 10) :
    (2 * W ^ (1 / (3 - 2 * σ)) * Real.log (W ^ (1 / (3 - 2 * σ))) ^ 2) ^ (2 - 2 * σ + 2 * ε₀)
      ≤ 8 / ε₀ ^ 2 * W ^ (ee σ + 3 * ε₀) := by
  have hW0 : 0 < W := lt_of_lt_of_le one_pos hW
  have h06 : (0.6 : ℝ) = 3 / 5 := by norm_num
  have hd : 0 < 3 - 2 * σ := by linarith only [h1]
  have hy0 : 0 ≤ 1 / (3 - 2 * σ) := (one_div_pos.mpr hd).le
  have hy1 : 1 / (3 - 2 * σ) ≤ 1 := by rw [div_le_one hd]; linarith only [h1]
  have hp0 : 0 ≤ 2 - 2 * σ + 2 * ε₀ := by linarith only [h1, he]
  have hp1 : 2 - 2 * σ + 2 * ε₀ ≤ 1 := by linarith only [h0, h06, he1]
  have step1 := classI_step1 W (1 / (3 - 2 * σ)) ε₀ (2 - 2 * σ + 2 * ε₀) hW hy0 hy1 he hp0
  have he2 : 0 < ε₀ ^ 2 := pow_pos he 2
  have hc1 : 1 ≤ 8 / ε₀ ^ 2 := by
    rw [le_div_iff₀ he2, one_mul]
    have h := pow_le_pow_left₀ he.le he1 2
    have h' : ((1 : ℝ) / 10) ^ 2 = 1 / 100 := by norm_num
    linarith only [h, h']
  have a1 : (8 / ε₀ ^ 2) ^ (2 - 2 * σ + 2 * ε₀) ≤ 8 / ε₀ ^ 2 := by
    have h := Real.rpow_le_rpow_of_exponent_le hc1 hp1
    rwa [Real.rpow_one] at h
  have a2 : (W ^ (1 / (3 - 2 * σ))) ^ (2 - 2 * σ + 2 * ε₀) ≤ W ^ (ee σ + 2 * ε₀) := by
    rw [← Real.rpow_mul hW0.le]
    exact Real.rpow_le_rpow_of_exponent_le hW (exp_ineq σ ε₀ he hy1)
  have a3 : (W ^ ε₀) ^ (2 - 2 * σ + 2 * ε₀) ≤ W ^ ε₀ := by
    rw [← Real.rpow_mul hW0.le]
    apply Real.rpow_le_rpow_of_exponent_le hW
    have h := mul_le_mul_of_nonneg_left hp1 he.le
    linarith only [h]
  have b1 : 0 ≤ (W ^ (1 / (3 - 2 * σ))) ^ (2 - 2 * σ + 2 * ε₀) :=
    Real.rpow_nonneg (Real.rpow_pos_of_pos hW0 _).le _
  have b2 : 0 ≤ (W ^ ε₀) ^ (2 - 2 * σ + 2 * ε₀) := Real.rpow_nonneg (Real.rpow_pos_of_pos hW0 _).le _
  have b4 : 0 ≤ W ^ (ee σ + 2 * ε₀) := (Real.rpow_pos_of_pos hW0 _).le
  have hc0 : 0 ≤ 8 / ε₀ ^ 2 := le_trans zero_le_one hc1
  have m1 := mul_le_mul a1 a2 b1 hc0
  have m2 := mul_le_mul m1 a3 b2 (mul_nonneg hc0 b4)
  have e : 8 / ε₀ ^ 2 * W ^ (ee σ + 2 * ε₀) * W ^ ε₀ = 8 / ε₀ ^ 2 * W ^ (ee σ + 3 * ε₀) := by
    rw [mul_assoc, ← Real.rpow_add hW0]; congr 2; ring
  rw [← e]
  exact step1.trans m2

end RingClass

