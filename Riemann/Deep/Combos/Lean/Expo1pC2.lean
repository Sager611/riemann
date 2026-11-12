import Riemann.Deep.Combos.Lean.Expo1pA

/-! # Theorem 1′(a): class I exponent, the exponent inequality (RCApps2 l.55, 59) -/

open Real

namespace RingClass

/-- y(2 − 2σ + 2ε₀) ≤ e(σ) + 2ε₀ for y = 1/(3 − 2σ) ≤ 1 -/
lemma exp_ineq (σ ε₀ : ℝ) (he : 0 < ε₀) (hy1 : 1 / (3 - 2 * σ) ≤ 1) :
    1 / (3 - 2 * σ) * (2 - 2 * σ + 2 * ε₀) ≤ ee σ + 2 * ε₀ := by
  have hee : ee σ = (2 * (1 - σ)) * (1 / (3 - 2 * σ)) := by
    unfold ee; exact div_eq_mul_one_div _ _
  rw [hee]
  have h := mul_le_mul_of_nonneg_left hy1 (le_of_lt (mul_pos two_pos he))
  linarith

/-- (abc)^p = a^p b^p c^p -/
lemma rpow_mul3 (a b c p : ℝ) (ha : 0 ≤ a) (hb : 0 ≤ b) (hc : 0 ≤ c) :
    (a * b * c) ^ p = a ^ p * b ^ p * c ^ p := by
  rw [Real.mul_rpow (mul_nonneg ha hb) hc, Real.mul_rpow ha hb]

end RingClass

