import Riemann.Deep.FasiHigham.Lean.Main

/-!
# Fasi–Higham, m = 3: the quantitative form c_n ≥ (R₀ − 2|R₊|) z₀^{−n} (real form)

qR k z = Σ_{j ≤ 3} qcoef k 3 j z^j is q_{k,3} on ℝ and qR' its derivative. For N = k + 1 and the root d:
z₀ = N + d is a zero of q, N(N+1)(N+2) q(z) = (z₀ − z)(z² − 2uz + v) with u² < v (the complex pair z_±,
|z_±|² = v) and z₀² < v (D1), R₀ = −1/(z₀ q'(z₀)), K = R₀² · (NOTE.md's closed form of (2|R₊|/R₀)²),
0 ≤ K < R₀² and (c_n − R₀ z₀^{−n})² ≤ K z₀^{−2n}; hence c_n ≥ (R₀ − r) z₀^{−n} for every r ≥ 0 with r² ≥ K,
in particular for r = √K = 2|R₊| < R₀.
-/

namespace FasiLean

open Finset

/-- q_{k,3}(z) for real z, from the general coefficient formula. -/
noncomputable def qR (k : ℕ) (z : ℝ) : ℝ := ∑ j ∈ range 4, (qcoef k 3 j : ℝ) * z ^ j

/-- The derivative q_{k,3}'(z) (termwise). -/
noncomputable def qR' (k : ℕ) (z : ℝ) : ℝ := ∑ j ∈ range 4, (j : ℝ) * (qcoef k 3 j : ℝ) * z ^ (j - 1)

theorem qR_eq (k : ℕ) (N z : ℝ) (hN : N = k + 1) : N * (N + 1) * (N + 2) * qR k z
    = N * (N + 1) * (N + 2) - 3 * N * (N + 1) * z + 3 * N * z ^ 2 - z ^ 3 := by
  subst hN
  have h1 : (k : ℝ) + 1 ≠ 0 := by positivity
  have h2 : (k : ℝ) + 2 ≠ 0 := by positivity
  have h3 : (k : ℝ) + 3 ≠ 0 := by positivity
  simp only [qR, Finset.sum_range_succ, Finset.sum_range_zero, qcoef_zero, qcoef_one, qcoef_two,
    qcoef_three]
  push_cast
  field_simp
  ring

theorem qR'_eq (k : ℕ) (N z : ℝ) (hN : N = k + 1) : N * (N + 1) * (N + 2) * qR' k z
    = -3 * N * (N + 1) + 6 * N * z - 3 * z ^ 2 := by
  subst hN
  have h1 : (k : ℝ) + 1 ≠ 0 := by positivity
  have h2 : (k : ℝ) + 2 ≠ 0 := by positivity
  have h3 : (k : ℝ) + 3 ≠ 0 := by positivity
  simp only [qR', Finset.sum_range_succ, Finset.sum_range_zero, qcoef_zero, qcoef_one, qcoef_two,
    qcoef_three]
  norm_num
  field_simp
  ring

/-- K/R₀² equals NOTE.md §2's closed form of (2|R₊|/R₀)² on the curve d³ + 3Nd = 2N. -/
theorem K_note (N d : ℝ) (hN : 1 ≤ N) (hd1 : 1 / 2 < d) (hd2 : d < 2 / 3)
    (hr : d ^ 3 + 3 * N * d = 2 * N) :
    3 * (N + d) ^ 2 * (d ^ 2 + N) / (vR N d * (vR N d - uR N d ^ 2))
      = 8 * N * (N + d) ^ 2 * (1 - d) ^ 2
        / (d ^ 3 * ((N + d) ^ 2 + 3 * N * (1 - d)) * (3 * N + d ^ 2 + d + 1)) := by
  have hN0 : 0 < N := by linarith
  have hd0 : 0 < d := by linarith
  have hv : vR N d = (N - d / 2) ^ 2 + 3 / 4 * d ^ 2 + 3 * N := by simp only [vR]; ring
  have hvu : vR N d - uR N d ^ 2 = 3 / 4 * d ^ 2 + 3 * N := by simp only [vR, uR]; ring
  have hv0 : 0 < vR N d := by rw [hv]; positivity
  have hvu0 : 0 < vR N d - uR N d ^ 2 := by rw [hvu]; positivity
  have hw : 0 < (N + d) ^ 2 + 3 * N * (1 - d) := by
    have : 0 < 3 * N * (1 - d) := mul_pos (by linarith) (by linarith)
    positivity
  rw [div_eq_div_iff (mul_pos hv0 hvu0).ne' (by positivity), hvu]
  simp only [vR]
  linear_combination (3 * (N + d) ^ 2 * (N ^ 2 - d * N + d ^ 2 + 3 * N)
    * (N * d ^ 2 - 2 * N * d + 4 * N + d ^ 4 + d ^ 3 + d ^ 2)) * hr

end FasiLean
