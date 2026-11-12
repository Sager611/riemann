import Riemann.Deep.OddRobin.Lean.NumA
import Riemann.Deep.OddRobin.Lean.SigVals

namespace OddRobin

open ArithmeticFunction
open scoped ArithmeticFunction.sigma

variable {L : ℝ → ℝ}

theorem LNN_ge (hL : LogLike L) : (608863983123 / 10000000000 : ℝ) ≤ L (NN : ℝ) := by
  have h1 : (2 : ℝ) ^ 87 * (8953383594299 / 5000000000000 : ℝ) ≤ (NN : ℝ) := by norm_num [NN]
  have h2 := hL.mono _ _ (by norm_num) h1
  rw [hL.two_pow_mul 87 (by norm_num)] at h2
  have h3 : (5825936036367 / 10000000000000 : ℝ) ≤ L (8953383594299 / 5000000000000 : ℝ) := hL.ge_of (q := (8953383594299 / 5000000000000 : ℝ)) (a := (3953383594299 / 5000000000000 : ℝ)) (by norm_num) (by norm_num) 11 (by norm_num [lsum, Finset.sum_range_succ])
  push_cast at h2; linarith [L2_ge hL]
theorem LNN_le (hL : LogLike L) : L (NN : ℝ) ≤ (152215995781 / 2500000000 : ℝ) := by
  have h1 : (NN : ℝ) ≤ (2 : ℝ) ^ 87 * (17906767188599 / 10000000000000 : ℝ) := by norm_num [NN]
  have h2 := hL.mono _ _ (by norm_num [NN]) h1
  rw [hL.two_pow_mul 87 (by norm_num)] at h2
  have h3 : L (17906767188599 / 10000000000000 : ℝ) ≤ (58259360363679 / 100000000000000 : ℝ) := hL.le_of (q := (17906767188599 / 10000000000000 : ℝ)) (a := (7906767188599 / 10000000000000 : ℝ)) (by norm_num) (by norm_num) 11 (by norm_num [lsum, ltail, Finset.sum_range_succ])
  push_cast at h2; linarith [L2_le hL]
theorem q_NN_A_lo_eq (hL : LogLike L) : L (608863983123 / 10000000000 : ℝ) = 5 * L 2 + L (608863983123 / 320000000000 : ℝ) := by
  have h := hL.two_pow_mul 5 (r := (608863983123 / 320000000000 : ℝ)) (by norm_num)
  rw [show (2 : ℝ) ^ 5 * (608863983123 / 320000000000 : ℝ) = (608863983123 / 10000000000 : ℝ) by norm_num] at h
  exact_mod_cast h
theorem q_NN_A_lo_ge (hL : LogLike L) : (5136262256431 / 1250000000000 : ℝ) ≤ L (608863983123 / 10000000000 : ℝ) := by
  have h : (8040923779317 / 12500000000000 : ℝ) ≤ L (608863983123 / 320000000000 : ℝ) := hL.ge_of (q := (608863983123 / 320000000000 : ℝ)) (a := (288863983123 / 320000000000 : ℝ)) (by norm_num) (by norm_num) 12 (by norm_num [lsum, Finset.sum_range_succ])
  rw [q_NN_A_lo_eq hL]; linarith [L2_ge hL]
theorem q_NN_A_lo_le (hL : LogLike L) : L (608863983123 / 10000000000 : ℝ) ≤ (20545049025727 / 5000000000000 : ℝ) := by
  have h : L (608863983123 / 320000000000 : ℝ) ≤ (64327390234539 / 100000000000000 : ℝ) := hL.le_of (q := (608863983123 / 320000000000 : ℝ)) (a := (288863983123 / 320000000000 : ℝ)) (by norm_num) (by norm_num) 12 (by norm_num [lsum, ltail, Finset.sum_range_succ])
  rw [q_NN_A_lo_eq hL]; linarith [L2_le hL]
theorem q_NN_A_hi_eq (hL : LogLike L) : L (152215995781 / 2500000000 : ℝ) = 5 * L 2 + L (152215995781 / 80000000000 : ℝ) := by
  have h := hL.two_pow_mul 5 (r := (152215995781 / 80000000000 : ℝ)) (by norm_num)
  rw [show (2 : ℝ) ^ 5 * (152215995781 / 80000000000 : ℝ) = (152215995781 / 2500000000 : ℝ) by norm_num] at h
  exact_mod_cast h
theorem q_NN_A_hi_ge (hL : LogLike L) : (8218019610293 / 2000000000000 : ℝ) ≤ L (152215995781 / 2500000000 : ℝ) := by
  have h : (643273902347 / 1000000000000 : ℝ) ≤ L (152215995781 / 80000000000 : ℝ) := hL.ge_of (q := (152215995781 / 80000000000 : ℝ)) (a := (72215995781 / 80000000000 : ℝ)) (by norm_num) (by norm_num) 12 (by norm_num [lsum, Finset.sum_range_succ])
  rw [q_NN_A_hi_eq hL]; linarith [L2_ge hL]
theorem q_NN_A_hi_le (hL : LogLike L) : L (152215995781 / 2500000000 : ℝ) ≤ (41090098051471 / 10000000000000 : ℝ) := by
  have h : L (152215995781 / 80000000000 : ℝ) ≤ (64327390234703 / 100000000000000 : ℝ) := hL.le_of (q := (152215995781 / 80000000000 : ℝ)) (a := (72215995781 / 80000000000 : ℝ)) (by norm_num) (by norm_num) 12 (by norm_num [lsum, ltail, Finset.sum_range_succ])
  rw [q_NN_A_hi_eq hL]; linarith [L2_le hL]
theorem LLNN_ge (hL : LogLike L) : (41090098051 / 10000000000 : ℝ) ≤ L (L (NN : ℝ)) := by
  have h := hL.mono _ _ (by norm_num) (LNN_ge hL)
  linarith [q_NN_A_lo_ge hL]
theorem LLNN_le (hL : LogLike L) : L (L (NN : ℝ)) ≤ (10272524513 / 2500000000 : ℝ) := by
  have h := hL.mono _ _ (by linarith [LNN_ge hL]) (LNN_le hL)
  linarith [q_NN_A_hi_le hL]
theorem q_NN_B_lo_eq (hL : LogLike L) : L (41090098051 / 10000000000 : ℝ) = 2 * L 2 + L (41090098051 / 40000000000 : ℝ) := by
  have h := hL.two_pow_mul 2 (r := (41090098051 / 40000000000 : ℝ)) (by norm_num)
  rw [show (2 : ℝ) ^ 2 * (41090098051 / 40000000000 : ℝ) = (41090098051 / 10000000000 : ℝ) by norm_num] at h
  exact_mod_cast h
theorem q_NN_B_lo_ge (hL : LogLike L) : (353295519037 / 250000000000 : ℝ) ≤ L (41090098051 / 10000000000 : ℝ) := by
  have h : (134438575141 / 5000000000000 : ℝ) ≤ L (41090098051 / 40000000000 : ℝ) := hL.ge_of (q := (41090098051 / 40000000000 : ℝ)) (a := (1090098051 / 40000000000 : ℝ)) (by norm_num) (by norm_num) 3 (by norm_num [lsum, Finset.sum_range_succ])
  rw [q_NN_B_lo_eq hL]; linarith [L2_ge hL]
theorem q_NN_B_lo_le (hL : LogLike L) : L (41090098051 / 10000000000 : ℝ) ≤ (14131820761483 / 10000000000000 : ℝ) := by
  have h : L (41090098051 / 40000000000 : ℝ) ≤ (336096437853 / 12500000000000 : ℝ) := hL.le_of (q := (41090098051 / 40000000000 : ℝ)) (a := (1090098051 / 40000000000 : ℝ)) (by norm_num) (by norm_num) 3 (by norm_num [lsum, ltail, Finset.sum_range_succ])
  rw [q_NN_B_lo_eq hL]; linarith [L2_le hL]
theorem q_NN_B_hi_eq (hL : LogLike L) : L (10272524513 / 2500000000 : ℝ) = 2 * L 2 + L (10272524513 / 10000000000 : ℝ) := by
  have h := hL.two_pow_mul 2 (r := (10272524513 / 10000000000 : ℝ)) (by norm_num)
  rw [show (2 : ℝ) ^ 2 * (10272524513 / 10000000000 : ℝ) = (10272524513 / 2500000000 : ℝ) by norm_num] at h
  exact_mod_cast h
theorem q_NN_B_hi_ge (hL : LogLike L) : (14131820761723 / 10000000000000 : ℝ) ≤ L (10272524513 / 2500000000 : ℝ) := by
  have h : (1344385752627 / 50000000000000 : ℝ) ≤ L (10272524513 / 10000000000 : ℝ) := hL.ge_of (q := (10272524513 / 10000000000 : ℝ)) (a := (272524513 / 10000000000 : ℝ)) (by norm_num) (by norm_num) 3 (by norm_num [lsum, Finset.sum_range_succ])
  rw [q_NN_B_hi_eq hL]; linarith [L2_ge hL]
theorem q_NN_B_hi_le (hL : LogLike L) : L (10272524513 / 2500000000 : ℝ) ≤ (7065910380863 / 5000000000000 : ℝ) := by
  have h : L (10272524513 / 10000000000 : ℝ) ≤ (2688771505257 / 100000000000000 : ℝ) := hL.le_of (q := (10272524513 / 10000000000 : ℝ)) (a := (272524513 / 10000000000 : ℝ)) (by norm_num) (by norm_num) 3 (by norm_num [lsum, ltail, Finset.sum_range_succ])
  rw [q_NN_B_hi_eq hL]; linarith [L2_le hL]
theorem LLLNN_ge (hL : LogLike L) : (14131820761 / 10000000000 : ℝ) ≤ L (L (L (NN : ℝ))) := by
  have h := hL.mono _ _ (by norm_num) (LLNN_ge hL)
  linarith [q_NN_B_lo_ge hL]
theorem LLLNN_le (hL : LogLike L) : L (L (L (NN : ℝ))) ≤ (7065910381 / 5000000000 : ℝ) := by
  have h := hL.mono _ _ (by linarith [LLNN_ge hL]) (LLNN_le hL)
  linarith [q_NN_B_hi_le hL]
theorem q_NN_Q_lo_eq (hL : LogLike L) : L (18318325699 / 5000000000 : ℝ) = 1 * L 2 + L (18318325699 / 10000000000 : ℝ) := by
  have h := hL.two_pow_mul 1 (r := (18318325699 / 10000000000 : ℝ)) (by norm_num)
  rw [show (2 : ℝ) ^ 1 * (18318325699 / 10000000000 : ℝ) = (18318325699 / 5000000000 : ℝ) by norm_num] at h
  exact_mod_cast h
theorem q_NN_Q_lo_ge (hL : LogLike L) : (12984640506593 / 10000000000000 : ℝ) ≤ L (18318325699 / 5000000000 : ℝ) := by
  have h : (3026584350497 / 5000000000000 : ℝ) ≤ L (18318325699 / 10000000000 : ℝ) := hL.ge_of (q := (18318325699 / 10000000000 : ℝ)) (a := (8318325699 / 10000000000 : ℝ)) (by norm_num) (by norm_num) 11 (by norm_num [lsum, Finset.sum_range_succ])
  rw [q_NN_Q_lo_eq hL]; linarith [L2_ge hL]
theorem q_NN_Q_lo_le (hL : LogLike L) : L (18318325699 / 5000000000 : ℝ) ≤ (2596928101319 / 2000000000000 : ℝ) := by
  have h : L (18318325699 / 10000000000 : ℝ) ≤ (30265843504973 / 50000000000000 : ℝ) := hL.le_of (q := (18318325699 / 10000000000 : ℝ)) (a := (8318325699 / 10000000000 : ℝ)) (by norm_num) (by norm_num) 11 (by norm_num [lsum, ltail, Finset.sum_range_succ])
  rw [q_NN_Q_lo_eq hL]; linarith [L2_le hL]
theorem q_NN_Q_hi_eq (hL : LogLike L) : L (3663665139801 / 1000000000000 : ℝ) = 1 * L 2 + L (3663665139801 / 2000000000000 : ℝ) := by
  have h := hL.two_pow_mul 1 (r := (3663665139801 / 2000000000000 : ℝ)) (by norm_num)
  rw [show (2 : ℝ) ^ 1 * (3663665139801 / 2000000000000 : ℝ) = (3663665139801 / 1000000000000 : ℝ) by norm_num] at h
  exact_mod_cast h
theorem q_NN_Q_hi_ge (hL : LogLike L) : (2596928101319 / 2000000000000 : ℝ) ≤ L (3663665139801 / 1000000000000 : ℝ) := by
  have h : (60531687009967 / 100000000000000 : ℝ) ≤ L (3663665139801 / 2000000000000 : ℝ) := hL.ge_of (q := (3663665139801 / 2000000000000 : ℝ)) (a := (1663665139801 / 2000000000000 : ℝ)) (by norm_num) (by norm_num) 11 (by norm_num [lsum, Finset.sum_range_succ])
  rw [q_NN_Q_hi_eq hL]; linarith [L2_ge hL]
theorem q_NN_Q_hi_le (hL : LogLike L) : L (3663665139801 / 1000000000000 : ℝ) ≤ (6492320253299 / 5000000000000 : ℝ) := by
  have h : L (3663665139801 / 2000000000000 : ℝ) ≤ (60531687009973 / 100000000000000 : ℝ) := hL.le_of (q := (3663665139801 / 2000000000000 : ℝ)) (a := (1663665139801 / 2000000000000 : ℝ)) (by norm_num) (by norm_num) 11 (by norm_num [lsum, ltail, Finset.sum_range_succ])
  rw [q_NN_Q_hi_eq hL]; linarith [L2_le hL]
theorem LQNN_le (hL : LogLike L) : L ((σ 1 NN : ℝ) / (NN : ℝ)) ≤ (6492320253299 / 5000000000000 : ℝ) := by
  rw [sig_NN]
  have h1 : ((1015178899037727378898944000 : ℕ) : ℝ) / (NN : ℝ) ≤ (3663665139801 / 1000000000000 : ℝ) := by norm_num [NN]
  have h2 := hL.mono _ _ (by norm_num [NN]) h1
  linarith [q_NN_Q_hi_le hL]
theorem LQNN_ge (hL : LogLike L) : (12984640506593 / 10000000000000 : ℝ) ≤ L ((σ 1 NN : ℝ) / (NN : ℝ)) := by
  rw [sig_NN]
  have h1 : (18318325699 / 5000000000 : ℝ) ≤ ((1015178899037727378898944000 : ℕ) : ℝ) / (NN : ℝ) := by norm_num [NN]
  have h2 := hL.mono _ _ (by norm_num) h1
  linarith [q_NN_Q_lo_ge hL]

end OddRobin
