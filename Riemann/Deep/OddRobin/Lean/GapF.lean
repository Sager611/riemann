import Riemann.Deep.OddRobin.Lean.SigVals

namespace OddRobin

open ArithmeticFunction
open scoped ArithmeticFunction.sigma

variable {L : ℝ → ℝ}

/-- 7^2 * 13 * 17 * ... * 61. -/
def KK : ℕ := 7 ^ 2 * 13 * 17 * 19 * 23 * 29 * 31 * 37 * 41 * 43 * 47 * 53 * 59 * 61

lemma final_cases (x3 x5 x11 x67 x71 : ℕ) (h3 : x3 = 4 ∨ x3 = 5) (h5 : x5 = 3 ∨ x5 = 2)
    (h11 : x11 = 1 ∨ x11 = 2) (h67 : x67 = 0 ∨ x67 = 1) (h71 : x71 = 0 ∨ x71 = 1)
    (x35 : ¬ (x3 = 5 ∧ x5 = 2)) (x311 : ¬ (x3 = 5 ∧ x11 = 2)) (x371 : ¬ (x3 = 5 ∧ x71 = 1))
    (x511 : ¬ (x5 = 2 ∧ x11 = 2)) (x571 : ¬ (x5 = 2 ∧ x71 = 1))
    (x1171 : ¬ (x11 = 2 ∧ x71 = 1))
    (hlo : NN < 3 ^ x3 * (5 ^ x5 * (11 ^ x11 * (67 ^ x67 * (71 ^ x71 * KK)))))
    (hhi : 3 ^ x3 * (5 ^ x5 * (11 ^ x11 * (67 ^ x67 * (71 ^ x71 * KK)))) < c0) :
    3 ^ x3 * (5 ^ x5 * (11 ^ x11 * (67 ^ x67 * (71 ^ x71 * KK)))) = m3 ∨
      3 ^ x3 * (5 ^ x5 * (11 ^ x11 * (67 ^ x67 * (71 ^ x71 * KK)))) = m11 ∨
      3 ^ x3 * (5 ^ x5 * (11 ^ x11 * (67 ^ x67 * (71 ^ x71 * KK)))) = m67 := by
  rcases h3 with rfl | rfl <;> rcases h5 with rfl | rfl <;> rcases h11 with rfl | rfl <;>
    rcases h67 with rfl | rfl <;> rcases h71 with rfl | rfl <;>
    norm_num [KK, NN, c0, m3, m11, m67] at *


end OddRobin
