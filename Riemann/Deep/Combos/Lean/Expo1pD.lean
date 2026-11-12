import Riemann.Deep.Combos.Lean.Expo1pA

/-! # Theorem 1′(a): window count and last step (RCApps2 l.57–61) -/

open Real

namespace RingClass

/-- the window count (RCApps2 l.57): log(|D|T²) ≤ 2 log W ≤ (2/ε₀) W^{ε₀}, W = 9|D|T -/
lemma window_log (absD T ε₀ : ℝ) (hD : 1 ≤ absD) (hT : 1 ≤ T) (he : 0 < ε₀) :
    Real.log (absD * T ^ 2) ≤ 2 / ε₀ * (9 * absD * T) ^ ε₀ := by
  have hW : 0 < 9 * absD * T := by positivity
  have h1 : absD * T ^ 2 ≤ (9 * absD * T) ^ 2 := by
    have h0 : absD ≤ absD ^ 2 := by rw [sq]; exact le_mul_of_one_le_right (by linarith) hD
    have : absD ≤ 81 * absD ^ 2 := by have := sq_nonneg absD; linarith
    have hT2 : 0 ≤ T ^ 2 := sq_nonneg T
    have e : (9 * absD * T) ^ 2 = 81 * absD ^ 2 * T ^ 2 := by ring
    rw [e]; exact mul_le_mul_of_nonneg_right this hT2
  have h2 : Real.log (absD * T ^ 2) ≤ Real.log ((9 * absD * T) ^ 2) :=
    Real.log_le_log (by positivity) h1
  rw [Real.log_pow] at h2
  have h3 := Real.log_le_rpow_div hW.le he
  have e : 2 / ε₀ * (9 * absD * T) ^ ε₀ = 2 * ((9 * absD * T) ^ ε₀ / ε₀) := by ring
  rw [e]; push_cast at h2; linarith

/-- the last step: W^{e+ε/2} ≤ 81 (|D|T)^{e+ε}, W = 9|D|T, for 0 < ε ≤ 1, 0 ≤ e ≤ 1 -/
lemma final_exp (absD T e ε : ℝ) (hD : 1 ≤ absD) (hT : 1 ≤ T) (he0 : 0 ≤ e) (he1 : e ≤ 1)
    (hε : 0 < ε) (hε1 : ε ≤ 1) :
    (9 * absD * T) ^ (e + ε / 2) ≤ 81 * (absD * T) ^ (e + ε) := by
  have hDT : 1 ≤ absD * T := one_le_mul_of_one_le_of_one_le hD hT
  rw [mul_assoc, Real.mul_rpow (by norm_num) (by linarith)]
  have h9 : (9 : ℝ) ^ (e + ε / 2) ≤ 81 := by
    calc (9 : ℝ) ^ (e + ε / 2) ≤ (9 : ℝ) ^ (2 : ℝ) :=
          Real.rpow_le_rpow_of_exponent_le (by norm_num) (by linarith)
      _ = 81 := by rw [Real.rpow_two]; norm_num
  have hb : (absD * T) ^ (e + ε / 2) ≤ (absD * T) ^ (e + ε) :=
    Real.rpow_le_rpow_of_exponent_le hDT (by linarith)
  exact mul_le_mul h9 hb (by positivity) (by norm_num)

end RingClass
