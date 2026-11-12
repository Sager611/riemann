import Riemann.Deep.Combos.Lean.Expo1pC1
import Riemann.Deep.Combos.Lean.Expo1pC2

/-! # Theorem 1′(a): class I exponent, first step (RCApps2 l.50–55) -/

open Real

namespace RingClass

/-- Z^p ≤ (8/ε₀²)^p (W^y)^p (W^{ε₀})^p for Z = 2W^y log²(W^y) -/
lemma classI_step1 (W y ε₀ p : ℝ) (hW : 1 ≤ W) (hy0 : 0 ≤ y) (hy1 : y ≤ 1) (he : 0 < ε₀)
    (hp0 : 0 ≤ p) :
    (2 * W ^ y * Real.log (W ^ y) ^ 2) ^ p
      ≤ (8 / ε₀ ^ 2) ^ p * (W ^ y) ^ p * (W ^ ε₀) ^ p := by
  have hW0 : 0 < W := by linarith
  have hZ := Z_le W y ε₀ hW hy0 hy1 he
  have hY : 0 < W ^ y := Real.rpow_pos_of_pos hW0 y
  have hZ0 : 0 ≤ 2 * W ^ y * Real.log (W ^ y) ^ 2 :=
    mul_nonneg (mul_nonneg zero_le_two hY.le) (sq_nonneg _)
  have h := Real.rpow_le_rpow hZ0 hZ hp0
  rw [rpow_mul3 _ _ _ _ (div_nonneg (by norm_num) (sq_nonneg ε₀)) hY.le
    (Real.rpow_pos_of_pos hW0 ε₀).le] at h
  exact h

end RingClass

