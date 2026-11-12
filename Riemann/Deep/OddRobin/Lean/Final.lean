import Riemann.Deep.OddRobin.Lean.BigN
import Riemann.Deep.OddRobin.Lean.Big3
import Riemann.Deep.OddRobin.Lean.Big11
import Riemann.Deep.OddRobin.Lean.Big67

namespace OddRobin

open ArithmeticFunction
open scoped ArithmeticFunction.sigma

variable {L : ℝ → ℝ}

theorem delta0_lt (hL : LogLike L) {g : ℝ} (hg : (5772145 / 10000000 : ℝ) ≤ g) :
    L ((σ 1 NN : ℝ) / (NN : ℝ)) - g + L 2 - L (L (L (NN : ℝ))) < 13 / 10000 := by
  linarith [LQNN_le hL, L2_le hL, LLLNN_ge hL]
theorem OR_m3 (hL : LogLike L) {E g : ℝ} (hE : 0 < E) (hLE : L E = g)
    (hg : (5772145 / 10000000 : ℝ) ≤ g) : ORf L E m3 := by
  rw [ORf_iff_log hL hE (by norm_num [m3]) (by linarith [LLm3_ge hL]), hLE]
  linarith [LQm3_le hL, L2_le hL, LLLm3_ge hL]
theorem OR_m11 (hL : LogLike L) {E g : ℝ} (hE : 0 < E) (hLE : L E = g)
    (hg : (5772145 / 10000000 : ℝ) ≤ g) : ORf L E m11 := by
  rw [ORf_iff_log hL hE (by norm_num [m11]) (by linarith [LLm11_ge hL]), hLE]
  linarith [LQm11_le hL, L2_le hL, LLLm11_ge hL]
theorem OR_m67 (hL : LogLike L) {E g : ℝ} (hE : 0 < E) (hLE : L E = g)
    (hg : (5772145 / 10000000 : ℝ) ≤ g) : ORf L E m67 := by
  rw [ORf_iff_log hL hE (by norm_num [m67]) (by linarith [LLm67_ge hL]), hLE]
  linarith [LQm67_le hL, L2_le hL, LLLm67_ge hL]
theorem not_OR_NN (hL : LogLike L) {E g : ℝ} (hE : 0 < E) (hLE : L E = g)
    (hg : g ≤ (5776 / 10000 : ℝ)) : ¬ ORf L E NN := by
  rw [ORf_iff_log hL hE (by norm_num [NN]) (by linarith [LLNN_ge hL]), hLE]
  exact not_lt.2 (by linarith [LQNN_ge hL, L2_ge hL, LLLNN_le hL])
theorem Lc0_le (hL : LogLike L) : L (c0 : ℝ) ≤ (16289 / 250 : ℝ) := by
  have h1 : (c0 : ℝ) ≤ (2 : ℝ) ^ 94 * 1 := by norm_num [c0]
  have h2 := hL.mono _ _ (by norm_num [c0]) h1
  rw [hL.two_pow_mul 94 (by norm_num), hL.map_one] at h2
  push_cast at h2; linarith [L2_le hL]
theorem q_T_eq (hL : LogLike L) : L (16289 / 250 : ℝ) = 6 * L 2 + L (16289 / 16000 : ℝ) := by
  have h := hL.two_pow_mul 6 (r := (16289 / 16000 : ℝ)) (by norm_num)
  rw [show (2 : ℝ) ^ 6 * (16289 / 16000 : ℝ) = (16289 / 250 : ℝ) by norm_num] at h
  exact_mod_cast h
theorem q_T_ge (hL : LogLike L) : (835356878899 / 200000000000 : ℝ) ≤ L (16289 / 250 : ℝ) := by
  have h : (55941597299 / 3125000000000 : ℝ) ≤ L (16289 / 16000 : ℝ) := hL.ge_of (q := (16289 / 16000 : ℝ)) (a := (289 / 16000 : ℝ)) (by norm_num) (by norm_num) 3 (by norm_num [lsum, Finset.sum_range_succ])
  rw [q_T_eq hL]; linarith [L2_ge hL]
theorem q_T_le (hL : LogLike L) : L (16289 / 250 : ℝ) ≤ (41767843944957 / 10000000000000 : ℝ) := by
  have h : L (16289 / 16000 : ℝ) ≤ (1790131113569 / 100000000000000 : ℝ) := hL.le_of (q := (16289 / 16000 : ℝ)) (a := (289 / 16000 : ℝ)) (by norm_num) (by norm_num) 3 (by norm_num [lsum, ltail, Finset.sum_range_succ])
  rw [q_T_eq hL]; linarith [L2_le hL]
theorem LT_le (hL : LogLike L) : L (16289 / 250 : ℝ) ≤ (835357 / 200000 : ℝ) := by linarith [q_T_le hL]
theorem eps_TT2 (hL : LogLike L) : eps L * ((16289 / 250 : ℝ) * (835357 / 200000 : ℝ)) ≤ 1 := by
  have := eps_le hL
  nlinarith
theorem LNN_gt_one (hL : LogLike L) : 1 < L (NN : ℝ) := by linarith [LNN_ge hL]

end OddRobin
