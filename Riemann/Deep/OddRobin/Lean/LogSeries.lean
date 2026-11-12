import Riemann.Deep.OddRobin.Lean.LogLike
import Mathlib.Analysis.SpecialFunctions.Log.Deriv
import Mathlib.Analysis.SpecificLimits.Basic

/-!
# The series bounds for Real.log (fields series_lo, series_hi of LogLike Real.log)

From Mathlib's Real.hasSum_log_one_add: log (1 + a) = sum 2 y^(2k+1)/(2k+1), y = a/(a+2).
Heavy imports (Log.Deriv); compiled at 09:03 UTC (watchdog 275 MB) with lsum/ltail defined locally;
this version takes the identical definitions from LogLike.lean and is not recompiled (memory cap).
-/

namespace OddRobin

open Real Finset

theorem lsum_le_log {a : ℝ} (ha : 0 ≤ a) (m : ℕ) : lsum a m ≤ log (1 + a) := by
  unfold lsum
  refine sum_le_hasSum (range m) (fun k _ => ?_) (hasSum_log_one_add ha)
  have : 0 ≤ a / (a + 2) := div_nonneg ha (by linarith)
  positivity

theorem log_le_lsum {a : ℝ} (ha : 0 ≤ a) (m : ℕ) : log (1 + a) ≤ lsum a m + ltail a m := by
  have hs := hasSum_log_one_add ha
  unfold lsum ltail
  set y := a / (a + 2) with hy
  have hy0 : 0 ≤ y := div_nonneg ha (by linarith)
  have hy1 : y < 1 := by rw [hy, div_lt_one (by linarith)]; linarith
  have hy2 : y ^ 2 < 1 := by nlinarith
  have ht := (hasSum_nat_add_iff' m).mpr hs
  have hg : HasSum (fun k : ℕ => (2 / (2 * (m : ℝ) + 1) * y ^ (2 * m + 1)) * (y ^ 2) ^ k)
      (2 / (2 * (m : ℝ) + 1) * y ^ (2 * m + 1) * (1 - y ^ 2)⁻¹) :=
    (hasSum_geometric_of_lt_one (by positivity) hy2).mul_left _
  have hle := hasSum_le (fun k => ?_) ht hg
  · have h1 : 0 < 1 - y ^ 2 := by linarith
    have e : 2 / (2 * (m : ℝ) + 1) * y ^ (2 * m + 1) * (1 - y ^ 2)⁻¹
        = 2 * y ^ (2 * m + 1) / ((2 * m + 1) * (1 - y ^ 2)) := by
      field_simp
    linarith
  · have e : y ^ (2 * (k + m) + 1) = y ^ (2 * m + 1) * (y ^ 2) ^ k := by
      rw [← pow_mul, ← pow_add]; ring_nf
    rw [e]
    have hk : 2 * (m : ℝ) + 1 ≤ 2 * ((k + m : ℕ) : ℝ) + 1 := by
      push_cast; linarith [(k.cast_nonneg : (0:ℝ) ≤ k)]
    have hk' : 1 / (2 * ((k + m : ℕ) : ℝ) + 1) ≤ 1 / (2 * (m : ℝ) + 1) :=
      one_div_le_one_div_of_le (by positivity) hk
    have hp : 0 ≤ y ^ (2 * m + 1) * (y ^ 2) ^ k := by positivity
    calc 2 * (1 / (2 * ((k + m : ℕ) : ℝ) + 1)) * (y ^ (2 * m + 1) * (y ^ 2) ^ k)
        ≤ 2 * (1 / (2 * (m : ℝ) + 1)) * (y ^ (2 * m + 1) * (y ^ 2) ^ k) := by gcongr
      _ = 2 / (2 * (m : ℝ) + 1) * y ^ (2 * m + 1) * (y ^ 2) ^ k := by ring

end OddRobin

