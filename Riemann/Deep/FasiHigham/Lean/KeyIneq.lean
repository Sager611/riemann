import Mathlib.Algebra.Order.Archimedean.Real.Basic
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.Positivity
import Mathlib.Tactic.Ring

/-! # Fasi–Higham, m = 3: the key inequality (no use of the cubic relation) -/

namespace FasiLean

/-- Key inequality a²(a² − 2ua + v) < v(v − u²) with a = N + d, u = N − d/2, v = N² − dN + d² + 3N,
for all N ≥ 1, d ∈ (1/2, 2/3): the difference is (9/4)((4 − 4d − d²)N² + (d² − 3d³)N − d⁴). -/
theorem key_ineq (N d : ℝ) (hN : 1 ≤ N) (hd1 : 1 / 2 < d) (hd2 : d < 2 / 3) :
    (N + d) ^ 2 * ((N + d) ^ 2 - 2 * (N - d / 2) * (N + d) + (N ^ 2 - d * N + d ^ 2 + 3 * N))
      < (N ^ 2 - d * N + d ^ 2 + 3 * N) * ((N ^ 2 - d * N + d ^ 2 + 3 * N) - (N - d / 2) ^ 2) := by
  have hd0 : (0 : ℝ) ≤ d := by linarith
  have hd2' : d ^ 2 < 4 / 9 := by
    calc d ^ 2 = d * d := sq d
      _ < 2 / 3 * (2 / 3) := mul_self_lt_mul_self hd0 hd2
      _ = 4 / 9 := by norm_num
  have hA : 8 / 9 < 4 - 4 * d - d ^ 2 := by linarith
  have hB : -(4 / 9) < d ^ 2 - 3 * d ^ 3 := by
    have h1 : d ^ 2 * (3 * d - 1) < 4 / 9 * (3 * d - 1) :=
      mul_lt_mul_of_pos_right hd2' (by linarith)
    linarith
  have hD : d ^ 4 < 16 / 81 := by
    calc d ^ 4 = d ^ 2 * d ^ 2 := by ring
      _ < 4 / 9 * (4 / 9) := mul_self_lt_mul_self (sq_nonneg d) hd2'
      _ = 16 / 81 := by norm_num
  have hN0 : 0 < N := by linarith
  have e : (N ^ 2 - d * N + d ^ 2 + 3 * N) * ((N ^ 2 - d * N + d ^ 2 + 3 * N) - (N - d / 2) ^ 2)
      - (N + d) ^ 2 * ((N + d) ^ 2 - 2 * (N - d / 2) * (N + d) + (N ^ 2 - d * N + d ^ 2 + 3 * N))
      = 9 / 4 * ((4 - 4 * d - d ^ 2) * N ^ 2 + (d ^ 2 - 3 * d ^ 3) * N - d ^ 4) := by ring
  have h1 : 8 / 9 * N ^ 2 < (4 - 4 * d - d ^ 2) * N ^ 2 := mul_lt_mul_of_pos_right hA (by positivity)
  have h2 : -(4 / 9) * N < (d ^ 2 - 3 * d ^ 3) * N := mul_lt_mul_of_pos_right hB hN0
  have h3 : 0 ≤ (N - 1) * (2 * N + 1) := mul_nonneg (by linarith) (by linarith)
  linarith

end FasiLean
