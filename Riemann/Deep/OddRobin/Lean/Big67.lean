import Riemann.Deep.OddRobin.Lean.NumA
import Riemann.Deep.OddRobin.Lean.SigVals

namespace OddRobin

open ArithmeticFunction
open scoped ArithmeticFunction.sigma

variable {L : ℝ → ℝ}

theorem Lm67_ge (hL : LogLike L) : (634816530193 / 10000000000 : ℝ) ≤ L (m67 : ℝ) := by
  have h1 : (2 : ℝ) ^ 91 * (14996917520451 / 10000000000000 : ℝ) ≤ (m67 : ℝ) := by norm_num [m67]
  have h2 := hL.mono _ _ (by norm_num) h1
  rw [hL.two_pow_mul 91 (by norm_num)] at h2
  have h3 : (20262979417691 / 50000000000000 : ℝ) ≤ L (14996917520451 / 10000000000000 : ℝ) := hL.ge_of (q := (14996917520451 / 10000000000000 : ℝ)) (a := (4996917520451 / 10000000000000 : ℝ)) (by norm_num) (by norm_num) 9 (by norm_num [lsum, Finset.sum_range_succ])
  push_cast at h2; linarith [L2_ge hL]
theorem Lm67_le (hL : LogLike L) : L (m67 : ℝ) ≤ (317408265097 / 5000000000 : ℝ) := by
  have h1 : (m67 : ℝ) ≤ (2 : ℝ) ^ 91 * (3749229380113 / 2500000000000 : ℝ) := by norm_num [m67]
  have h2 := hL.mono _ _ (by norm_num [m67]) h1
  rw [hL.two_pow_mul 91 (by norm_num)] at h2
  have h3 : L (3749229380113 / 2500000000000 : ℝ) ≤ (4052595883539 / 10000000000000 : ℝ) := hL.le_of (q := (3749229380113 / 2500000000000 : ℝ)) (a := (1249229380113 / 2500000000000 : ℝ)) (by norm_num) (by norm_num) 9 (by norm_num [lsum, ltail, Finset.sum_range_succ])
  push_cast at h2; linarith [L2_le hL]
theorem q_m67_A_lo_eq (hL : LogLike L) : L (634816530193 / 10000000000 : ℝ) = 5 * L 2 + L (634816530193 / 320000000000 : ℝ) := by
  have h := hL.two_pow_mul 5 (r := (634816530193 / 320000000000 : ℝ)) (by norm_num)
  rw [show (2 : ℝ) ^ 5 * (634816530193 / 320000000000 : ℝ) = (634816530193 / 10000000000 : ℝ) by norm_num] at h
  exact_mod_cast h
theorem q_m67_A_lo_ge (hL : LogLike L) : (10376877338301 / 2500000000000 : ℝ) ≤ L (634816530193 / 10000000000 : ℝ) := by
  have h : (34250751626047 / 50000000000000 : ℝ) ≤ L (634816530193 / 320000000000 : ℝ) := hL.ge_of (q := (634816530193 / 320000000000 : ℝ)) (a := (314816530193 / 320000000000 : ℝ)) (by norm_num) (by norm_num) 12 (by norm_num [lsum, Finset.sum_range_succ])
  rw [q_m67_A_lo_eq hL]; linarith [L2_ge hL]
theorem q_m67_A_lo_le (hL : LogLike L) : L (634816530193 / 10000000000 : ℝ) ≤ (41507509353211 / 10000000000000 : ℝ) := by
  have h : L (634816530193 / 320000000000 : ℝ) ≤ (68501503252103 / 100000000000000 : ℝ) := hL.le_of (q := (634816530193 / 320000000000 : ℝ)) (a := (314816530193 / 320000000000 : ℝ)) (by norm_num) (by norm_num) 12 (by norm_num [lsum, ltail, Finset.sum_range_succ])
  rw [q_m67_A_lo_eq hL]; linarith [L2_le hL]
theorem q_m67_A_hi_eq (hL : LogLike L) : L (317408265097 / 5000000000 : ℝ) = 5 * L 2 + L (317408265097 / 160000000000 : ℝ) := by
  have h := hL.two_pow_mul 5 (r := (317408265097 / 160000000000 : ℝ)) (by norm_num)
  rw [show (2 : ℝ) ^ 5 * (317408265097 / 160000000000 : ℝ) = (317408265097 / 5000000000 : ℝ) by norm_num] at h
  exact_mod_cast h
theorem q_m67_A_hi_ge (hL : LogLike L) : (2075375467661 / 500000000000 : ℝ) ≤ L (317408265097 / 5000000000 : ℝ) := by
  have h : (68501503252251 / 100000000000000 : ℝ) ≤ L (317408265097 / 160000000000 : ℝ) := hL.ge_of (q := (317408265097 / 160000000000 : ℝ)) (a := (157408265097 / 160000000000 : ℝ)) (by norm_num) (by norm_num) 12 (by norm_num [lsum, Finset.sum_range_succ])
  rw [q_m67_A_hi_eq hL]; linarith [L2_ge hL]
theorem q_m67_A_hi_le (hL : LogLike L) : L (317408265097 / 5000000000 : ℝ) ≤ (20753754676613 / 5000000000000 : ℝ) := by
  have h : L (317408265097 / 160000000000 : ℝ) ≤ (3425075162613 / 5000000000000 : ℝ) := hL.le_of (q := (317408265097 / 160000000000 : ℝ)) (a := (157408265097 / 160000000000 : ℝ)) (by norm_num) (by norm_num) 12 (by norm_num [lsum, ltail, Finset.sum_range_succ])
  rw [q_m67_A_hi_eq hL]; linarith [L2_le hL]
theorem LLm67_ge (hL : LogLike L) : (41507509353 / 10000000000 : ℝ) ≤ L (L (m67 : ℝ)) := by
  have h := hL.mono _ _ (by norm_num) (Lm67_ge hL)
  linarith [q_m67_A_lo_ge hL]
theorem LLm67_le (hL : LogLike L) : L (L (m67 : ℝ)) ≤ (20753754677 / 5000000000 : ℝ) := by
  have h := hL.mono _ _ (by linarith [Lm67_ge hL]) (Lm67_le hL)
  linarith [q_m67_A_hi_le hL]
theorem q_m67_B_lo_eq (hL : LogLike L) : L (41507509353 / 10000000000 : ℝ) = 2 * L 2 + L (41507509353 / 40000000000 : ℝ) := by
  have h := hL.two_pow_mul 2 (r := (41507509353 / 40000000000 : ℝ)) (by norm_num)
  rw [show (2 : ℝ) ^ 2 * (41507509353 / 40000000000 : ℝ) = (41507509353 / 10000000000 : ℝ) by norm_num] at h
  exact_mod_cast h
theorem q_m67_B_lo_ge (hL : LogLike L) : (1779111582673 / 1250000000000 : ℝ) ≤ L (41507509353 / 10000000000 : ℝ) := by
  have h : (3699490501861 / 100000000000000 : ℝ) ≤ L (41507509353 / 40000000000 : ℝ) := hL.ge_of (q := (41507509353 / 40000000000 : ℝ)) (a := (1507509353 / 40000000000 : ℝ)) (by norm_num) (by norm_num) 4 (by norm_num [lsum, Finset.sum_range_succ])
  rw [q_m67_B_lo_eq hL]; linarith [L2_ge hL]
theorem q_m67_B_lo_le (hL : LogLike L) : L (41507509353 / 10000000000 : ℝ) ≤ (14232892661387 / 10000000000000 : ℝ) := by
  have h : L (41507509353 / 40000000000 : ℝ) ≤ (1849745250931 / 50000000000000 : ℝ) := hL.le_of (q := (41507509353 / 40000000000 : ℝ)) (a := (1507509353 / 40000000000 : ℝ)) (by norm_num) (by norm_num) 4 (by norm_num [lsum, ltail, Finset.sum_range_succ])
  rw [q_m67_B_lo_eq hL]; linarith [L2_le hL]
theorem q_m67_B_hi_eq (hL : LogLike L) : L (20753754677 / 5000000000 : ℝ) = 2 * L 2 + L (20753754677 / 20000000000 : ℝ) := by
  have h := hL.two_pow_mul 2 (r := (20753754677 / 20000000000 : ℝ)) (by norm_num)
  rw [show (2 : ℝ) ^ 2 * (20753754677 / 20000000000 : ℝ) = (20753754677 / 5000000000 : ℝ) by norm_num] at h
  exact_mod_cast h
theorem q_m67_B_hi_ge (hL : LogLike L) : (113863141293 / 80000000000 : ℝ) ≤ L (20753754677 / 5000000000 : ℝ) := by
  have h : (369949050427 / 10000000000000 : ℝ) ≤ L (20753754677 / 20000000000 : ℝ) := hL.ge_of (q := (20753754677 / 20000000000 : ℝ)) (a := (753754677 / 20000000000 : ℝ)) (by norm_num) (by norm_num) 4 (by norm_num [lsum, Finset.sum_range_succ])
  rw [q_m67_B_hi_eq hL]; linarith [L2_ge hL]
theorem q_m67_B_hi_le (hL : LogLike L) : L (20753754677 / 5000000000 : ℝ) ≤ (3558223165407 / 2500000000000 : ℝ) := by
  have h : L (20753754677 / 20000000000 : ℝ) ≤ (3699490504271 / 100000000000000 : ℝ) := hL.le_of (q := (20753754677 / 20000000000 : ℝ)) (a := (753754677 / 20000000000 : ℝ)) (by norm_num) (by norm_num) 4 (by norm_num [lsum, ltail, Finset.sum_range_succ])
  rw [q_m67_B_hi_eq hL]; linarith [L2_le hL]
theorem LLLm67_ge (hL : LogLike L) : (14232892661 / 10000000000 : ℝ) ≤ L (L (L (m67 : ℝ))) := by
  have h := hL.mono _ _ (by norm_num) (LLm67_ge hL)
  linarith [q_m67_B_lo_ge hL]
theorem LLLm67_le (hL : LogLike L) : L (L (L (m67 : ℝ))) ≤ (7116446331 / 5000000000 : ℝ) := by
  have h := hL.mono _ _ (by linarith [LLm67_ge hL]) (LLm67_le hL)
  linarith [q_m67_B_hi_le hL]
theorem q_m67_Q_lo_eq (hL : LogLike L) : L (3694511153223 / 1000000000000 : ℝ) = 1 * L 2 + L (3694511153223 / 2000000000000 : ℝ) := by
  have h := hL.two_pow_mul 1 (r := (3694511153223 / 2000000000000 : ℝ)) (by norm_num)
  rw [show (2 : ℝ) ^ 1 * (3694511153223 / 2000000000000 : ℝ) = (3694511153223 / 1000000000000 : ℝ) by norm_num] at h
  exact_mod_cast h
theorem q_m67_Q_lo_ge (hL : LogLike L) : (13068482461139 / 10000000000000 : ℝ) ≤ L (3694511153223 / 1000000000000 : ℝ) := by
  have h : (61370106555409 / 100000000000000 : ℝ) ≤ L (3694511153223 / 2000000000000 : ℝ) := hL.ge_of (q := (3694511153223 / 2000000000000 : ℝ)) (a := (1694511153223 / 2000000000000 : ℝ)) (by norm_num) (by norm_num) 11 (by norm_num [lsum, Finset.sum_range_succ])
  rw [q_m67_Q_lo_eq hL]; linarith [L2_ge hL]
theorem q_m67_Q_lo_le (hL : LogLike L) : L (3694511153223 / 1000000000000 : ℝ) ≤ (6534241230571 / 5000000000000 : ℝ) := by
  have h : L (3694511153223 / 2000000000000 : ℝ) ≤ (30685053277709 / 50000000000000 : ℝ) := hL.le_of (q := (3694511153223 / 2000000000000 : ℝ)) (a := (1694511153223 / 2000000000000 : ℝ)) (by norm_num) (by norm_num) 11 (by norm_num [lsum, ltail, Finset.sum_range_succ])
  rw [q_m67_Q_lo_eq hL]; linarith [L2_le hL]
theorem q_m67_Q_hi_eq (hL : LogLike L) : L (461813894153 / 125000000000 : ℝ) = 1 * L 2 + L (461813894153 / 250000000000 : ℝ) := by
  have h := hL.two_pow_mul 1 (r := (461813894153 / 250000000000 : ℝ)) (by norm_num)
  rw [show (2 : ℝ) ^ 1 * (461813894153 / 250000000000 : ℝ) = (461813894153 / 125000000000 : ℝ) by norm_num] at h
  exact_mod_cast h
theorem q_m67_Q_hi_ge (hL : LogLike L) : (6534241230571 / 5000000000000 : ℝ) ≤ L (461813894153 / 125000000000 : ℝ) := by
  have h : (15342526638859 / 25000000000000 : ℝ) ≤ L (461813894153 / 250000000000 : ℝ) := hL.ge_of (q := (461813894153 / 250000000000 : ℝ)) (a := (211813894153 / 250000000000 : ℝ)) (by norm_num) (by norm_num) 11 (by norm_num [lsum, Finset.sum_range_succ])
  rw [q_m67_Q_hi_eq hL]; linarith [L2_ge hL]
theorem q_m67_Q_hi_le (hL : LogLike L) : L (461813894153 / 125000000000 : ℝ) ≤ (2613696492229 / 2000000000000 : ℝ) := by
  have h : L (461813894153 / 250000000000 : ℝ) ≤ (12274021311089 / 20000000000000 : ℝ) := hL.le_of (q := (461813894153 / 250000000000 : ℝ)) (a := (211813894153 / 250000000000 : ℝ)) (by norm_num) (by norm_num) 11 (by norm_num [lsum, ltail, Finset.sum_range_succ])
  rw [q_m67_Q_hi_eq hL]; linarith [L2_le hL]
theorem LQm67_le (hL : LogLike L) : L ((σ 1 m67 : ℝ) / (m67 : ℝ)) ≤ (2613696492229 / 2000000000000 : ℝ) := by
  rw [sig_m67]
  have h1 : ((13717930251099546889224192000 : ℕ) : ℝ) / (m67 : ℝ) ≤ (461813894153 / 125000000000 : ℝ) := by norm_num [m67]
  have h2 := hL.mono _ _ (by norm_num [m67]) h1
  linarith [q_m67_Q_hi_le hL]
theorem LQm67_ge (hL : LogLike L) : (13068482461139 / 10000000000000 : ℝ) ≤ L ((σ 1 m67 : ℝ) / (m67 : ℝ)) := by
  rw [sig_m67]
  have h1 : (3694511153223 / 1000000000000 : ℝ) ≤ ((13717930251099546889224192000 : ℕ) : ℝ) / (m67 : ℝ) := by norm_num [m67]
  have h2 := hL.mono _ _ (by norm_num) h1
  linarith [q_m67_Q_lo_ge hL]

end OddRobin
