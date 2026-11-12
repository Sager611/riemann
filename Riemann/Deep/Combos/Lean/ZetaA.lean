import Mathlib.Analysis.SpecialFunctions.Pow.Real
import Mathlib.Algebra.Order.BigOperators.Ring.Finset

/-!
# ζ(1 + r) ≤ 1 + 1/r, in the finite form used by Lemma 2 (NOTE l.89)

For r > 0 and every finite set K of positive integers: Σ_{k∈K} k^{−(1+r)} ≤ 1 + 1/r.
Elementary: ((m+1)/m)^r ≥ 1 + r log((m+1)/m) ≥ 1 + r/(m+1), which telescopes.
-/

open Finset

namespace RingClass

lemma ratio_bound (r m : ℝ) (hr : 0 < r) (hm : 1 ≤ m) :
    1 + r / (m + 1) ≤ (m + 1) ^ r / m ^ r := by
  have hm0 : 0 < m := by linarith
  have hm1 : 0 < m + 1 := by linarith
  rw [← Real.div_rpow hm1.le hm0.le, Real.rpow_def_of_pos (div_pos hm1 hm0)]
  have h1 : 1 - ((m + 1) / m)⁻¹ ≤ Real.log ((m + 1) / m) :=
    Real.one_sub_inv_le_log_of_pos (div_pos hm1 hm0)
  have h2 : 1 - ((m + 1) / m)⁻¹ = 1 / (m + 1) := by
    rw [inv_div, eq_div_iff hm1.ne', sub_mul, div_mul_cancel₀ _ hm1.ne']; ring
  have h3 := Real.add_one_le_exp (Real.log ((m + 1) / m) * r)
  have h4 : 1 / (m + 1) * r ≤ Real.log ((m + 1) / m) * r :=
    mul_le_mul_of_nonneg_right (h2 ▸ h1) hr.le
  have h5 : r / (m + 1) = 1 / (m + 1) * r := by ring
  linarith

/-- (m+1)^{−(1+r)} ≤ (m^{−r} − (m+1)^{−r})/r for m ≥ 1 -/
lemma tele_step (r m : ℝ) (hr : 0 < r) (hm : 1 ≤ m) :
    ((m + 1) ^ (1 + r))⁻¹ ≤ ((m ^ r)⁻¹ - ((m + 1) ^ r)⁻¹) / r := by
  have hm0 : 0 < m := by linarith
  have hm1 : 0 < m + 1 := by linarith
  have ha : 0 < m ^ r := Real.rpow_pos_of_pos hm0 r
  have hb : 0 < (m + 1) ^ r := Real.rpow_pos_of_pos hm1 r
  have hsplit : (m + 1) ^ (1 + r) = (m + 1) * (m + 1) ^ r := by
    rw [Real.rpow_add hm1, Real.rpow_one]
  have hq := ratio_bound r m hr hm
  rw [le_div_iff₀ ha] at hq
  -- hq : (1 + r/(m+1)) * m^r ≤ (m+1)^r
  have key : r * m ^ r ≤ (m + 1) * ((m + 1) ^ r - m ^ r) := by
    have e : (1 + r / (m + 1)) * m ^ r = m ^ r + r * m ^ r / (m + 1) := by ring
    rw [e] at hq
    have h6 : r * m ^ r / (m + 1) ≤ (m + 1) ^ r - m ^ r := by linarith
    rw [div_le_iff₀ hm1] at h6
    linarith
  rw [hsplit, inv_eq_one_div, inv_eq_one_div, inv_eq_one_div, div_sub_div _ _ ha.ne' hb.ne',
    div_div, div_le_div_iff₀ (by positivity) (by positivity)]
  have key2 := mul_le_mul_of_nonneg_right key hb.le
  linarith [key2]

end RingClass
