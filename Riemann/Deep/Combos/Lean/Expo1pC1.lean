module

public import Riemann.Deep.Combos.Lean.Expo1pA

@[expose] public section

/-! # Theorem 1′(a): class I exponent, logarithms (RCApps2 l.50–55) -/

open Real

namespace RingClass

/-- log²Y ≤ (4/ε₀²) W^{ε₀} for Y = W^y, W ≥ 1, 0 ≤ y ≤ 1 -/
lemma logY_sq_le (W y ε₀ : ℝ) (hW : 1 ≤ W) (hy0 : 0 ≤ y) (hy1 : y ≤ 1) (he : 0 < ε₀) :
    Real.log (W ^ y) ^ 2 ≤ 4 / ε₀ ^ 2 * W ^ ε₀ := by
  have hW0 : 0 < W := by linarith
  have hlW : 0 ≤ Real.log W := Real.log_nonneg hW
  have hlY : Real.log (W ^ y) = y * Real.log W := Real.log_rpow hW0 y
  have hlY0 : 0 ≤ Real.log (W ^ y) := by rw [hlY]; exact mul_nonneg hy0 hlW
  have hlYW : Real.log (W ^ y) ≤ Real.log W := by
    rw [hlY]; exact mul_le_of_le_one_left hlW hy1
  have hlog : Real.log W ≤ W ^ (ε₀ / 2) / (ε₀ / 2) := Real.log_le_rpow_div hW0.le (by linarith)
  have hsq : (W ^ (ε₀ / 2)) ^ 2 = W ^ ε₀ := by
    rw [← Real.rpow_natCast, ← Real.rpow_mul hW0.le]; congr 1; push_cast; ring
  have h1 : Real.log (W ^ y) ^ 2 ≤ (W ^ (ε₀ / 2) / (ε₀ / 2)) ^ 2 :=
    pow_le_pow_left₀ hlY0 (hlYW.trans hlog) 2
  have h2 : (W ^ (ε₀ / 2) / (ε₀ / 2)) ^ 2 = 4 / ε₀ ^ 2 * (W ^ (ε₀ / 2)) ^ 2 := by
    rw [div_pow]; ring
  rw [h2, hsq] at h1
  exact h1

/-- Z = 2Y log²Y ≤ (8/ε₀²) Y W^{ε₀} -/
lemma Z_le (W y ε₀ : ℝ) (hW : 1 ≤ W) (hy0 : 0 ≤ y) (hy1 : y ≤ 1) (he : 0 < ε₀) :
    2 * W ^ y * Real.log (W ^ y) ^ 2 ≤ 8 / ε₀ ^ 2 * W ^ y * W ^ ε₀ := by
  have hW0 : 0 < W := by linarith
  have hl2 := logY_sq_le W y ε₀ hW hy0 hy1 he
  have hY : 0 ≤ 2 * W ^ y := by have := Real.rpow_pos_of_pos hW0 y; linarith
  have h := mul_le_mul_of_nonneg_left hl2 hY
  have e : 2 * W ^ y * (4 / ε₀ ^ 2 * W ^ ε₀) = 8 / ε₀ ^ 2 * W ^ y * W ^ ε₀ := by ring
  linarith

end RingClass

