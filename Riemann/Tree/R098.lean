import Mathlib

/-!
# Round 98 — analytic Weil positivity just beyond log 2

g = f * f~ supported in (−L, L). For L < log 2 there is no prime term (Yoshida); for log 2 < L < log 3
only n = 2 enters, through g(log 2) = ∫ over an overlap of length δ = L − log 2.

| idea | verdict |
|------|---------|
| bound g(log 2) by δ·‖f‖∞² or overlap L² mass | dead: mass can sit at both support ends, g(log 2) = ‖f‖²/2 with no δ-decay |
| uniform margin of the archimedean + pole form at L = log 2 | plausible (numerics: positivity far beyond log 2) |
| log-energy coercivity: endpoint mass ≲ E(f)/log(1/δ) | gives an analytic window log 2 + δ, δ tiny |
| sign of Re g(log 2) | either sign (Bochner): the prime term can help or hurt |
| Burnol S = {∞, 2} framework | natural setting for log 2 < L < log 3; no published proof found |
| novelty | an analytic δ-window would be weaker than certified numerical windows (ours 1.75) |

Lean below: the support lemmas behind "no prime term below log 2" and "only n = 2 below log 3".
-/

namespace Fayed.R98

open MeasureTheory

/-- If f vanishes outside [−a, a] and x > 2a, the autocorrelation at lag x vanishes. -/
theorem autocorr_eq_zero {f : ℝ → ℝ} {a x : ℝ} (hf : ∀ y, a < |y| → f y = 0)
    (hx : 2 * a < x) : ∫ y, f y * f (y - x) = 0 := by
  have h0 : ∀ y, f y * f (y - x) = 0 := by
    intro y
    by_cases h : |y| ≤ a
    · have hy := (abs_le.mp h).2
      have h2 : a < |y - x| := by
        rw [abs_sub_comm]
        exact lt_of_lt_of_le (by linarith) (le_abs_self (x - y))
      rw [hf (y - x) h2, mul_zero]
    · rw [hf y (not_le.mp h), zero_mul]
  simp [h0]

/-- Windows shorter than log 3: every lag log n with n ≥ 3 is inactive, only n = 2 can enter. -/
theorem only_two_active {f : ℝ → ℝ} {a : ℝ} (hf : ∀ y, a < |y| → f y = 0)
    (ha : 2 * a < Real.log 3) {n : ℕ} (hn : 3 ≤ n) : ∫ y, f y * f (y - Real.log n) = 0 := by
  apply autocorr_eq_zero hf
  have h3 : Real.log 3 ≤ Real.log n := Real.log_le_log (by norm_num) (by exact_mod_cast hn)
  linarith

end Fayed.R98
