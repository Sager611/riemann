import Riemann.Deep.OddRobin.Lean.NumA
import Riemann.Deep.OddRobin.Lean.SigVals

namespace OddRobin

open ArithmeticFunction
open scoped ArithmeticFunction.sigma

variable {L : ℝ → ℝ}

theorem Lm11_ge (hL : LogLike L) : (632842935851 / 10000000000 : ℝ) ≤ L (m11 : ℝ) := by
  have h1 : (2 : ℝ) ^ 91 * (12310902442161 / 10000000000000 : ℝ) ≤ (m11 : ℝ) := by norm_num [m11]
  have h2 := hL.mono _ _ (by norm_num) h1
  rw [hL.two_pow_mul 91 (by norm_num)] at h2
  have h3 : (2598751927441 / 12500000000000 : ℝ) ≤ L (12310902442161 / 10000000000000 : ℝ) := hL.ge_of (q := (12310902442161 / 10000000000000 : ℝ)) (a := (2310902442161 / 10000000000000 : ℝ)) (by norm_num) (by norm_num) 6 (by norm_num [lsum, Finset.sum_range_succ])
  push_cast at h2; linarith [L2_ge hL]
theorem Lm11_le (hL : LogLike L) : L (m11 : ℝ) ≤ (158210733963 / 2500000000 : ℝ) := by
  have h1 : (m11 : ℝ) ≤ (2 : ℝ) ^ 91 * (6155451221081 / 5000000000000 : ℝ) := by norm_num [m11]
  have h2 := hL.mono _ _ (by norm_num [m11]) h1
  rw [hL.two_pow_mul 91 (by norm_num)] at h2
  have h3 : L (6155451221081 / 5000000000000 : ℝ) ≤ (1039500770977 / 5000000000000 : ℝ) := hL.le_of (q := (6155451221081 / 5000000000000 : ℝ)) (a := (1155451221081 / 5000000000000 : ℝ)) (by norm_num) (by norm_num) 6 (by norm_num [lsum, ltail, Finset.sum_range_succ])
  push_cast at h2; linarith [L2_le hL]
theorem q_m11_A_lo_eq (hL : LogLike L) : L (632842935851 / 10000000000 : ℝ) = 5 * L 2 + L (632842935851 / 320000000000 : ℝ) := by
  have h := hL.two_pow_mul 5 (r := (632842935851 / 320000000000 : ℝ)) (by norm_num)
  rw [show (2 : ℝ) ^ 5 * (632842935851 / 320000000000 : ℝ) = (632842935851 / 10000000000 : ℝ) by norm_num] at h
  exact_mod_cast h
theorem q_m11_A_lo_ge (hL : LogLike L) : (10369092929359 / 2500000000000 : ℝ) ≤ L (632842935851 / 10000000000 : ℝ) := by
  have h : (68190126894413 / 100000000000000 : ℝ) ≤ L (632842935851 / 320000000000 : ℝ) := hL.ge_of (q := (632842935851 / 320000000000 : ℝ)) (a := (312842935851 / 320000000000 : ℝ)) (by norm_num) (by norm_num) 12 (by norm_num [lsum, Finset.sum_range_succ])
  rw [q_m11_A_lo_eq hL]; linarith [L2_ge hL]
theorem q_m11_A_lo_le (hL : LogLike L) : L (632842935851 / 10000000000 : ℝ) ≤ (41476371717443 / 10000000000000 : ℝ) := by
  have h : L (632842935851 / 320000000000 : ℝ) ≤ (68190126894421 / 100000000000000 : ℝ) := hL.le_of (q := (632842935851 / 320000000000 : ℝ)) (a := (312842935851 / 320000000000 : ℝ)) (by norm_num) (by norm_num) 12 (by norm_num [lsum, ltail, Finset.sum_range_succ])
  rw [q_m11_A_lo_eq hL]; linarith [L2_le hL]
theorem q_m11_A_hi_eq (hL : LogLike L) : L (158210733963 / 2500000000 : ℝ) = 5 * L 2 + L (158210733963 / 80000000000 : ℝ) := by
  have h := hL.two_pow_mul 5 (r := (158210733963 / 80000000000 : ℝ)) (by norm_num)
  rw [show (2 : ℝ) ^ 5 * (158210733963 / 80000000000 : ℝ) = (158210733963 / 2500000000 : ℝ) by norm_num] at h
  exact_mod_cast h
theorem q_m11_A_hi_ge (hL : LogLike L) : (10369092929363 / 2500000000000 : ℝ) ≤ L (158210733963 / 2500000000 : ℝ) := by
  have h : (68190126894571 / 100000000000000 : ℝ) ≤ L (158210733963 / 80000000000 : ℝ) := hL.ge_of (q := (158210733963 / 80000000000 : ℝ)) (a := (78210733963 / 80000000000 : ℝ)) (by norm_num) (by norm_num) 12 (by norm_num [lsum, Finset.sum_range_succ])
  rw [q_m11_A_hi_eq hL]; linarith [L2_ge hL]
theorem q_m11_A_hi_le (hL : LogLike L) : L (158210733963 / 2500000000 : ℝ) ≤ (20738185858729 / 5000000000000 : ℝ) := by
  have h : L (158210733963 / 80000000000 : ℝ) ≤ (68190126894579 / 100000000000000 : ℝ) := hL.le_of (q := (158210733963 / 80000000000 : ℝ)) (a := (78210733963 / 80000000000 : ℝ)) (by norm_num) (by norm_num) 12 (by norm_num [lsum, ltail, Finset.sum_range_succ])
  rw [q_m11_A_hi_eq hL]; linarith [L2_le hL]
theorem LLm11_ge (hL : LogLike L) : (41476371717 / 10000000000 : ℝ) ≤ L (L (m11 : ℝ)) := by
  have h := hL.mono _ _ (by norm_num) (Lm11_ge hL)
  linarith [q_m11_A_lo_ge hL]
theorem LLm11_le (hL : LogLike L) : L (L (m11 : ℝ)) ≤ (20738185859 / 5000000000 : ℝ) := by
  have h := hL.mono _ _ (by linarith [Lm11_ge hL]) (Lm11_le hL)
  linarith [q_m11_A_hi_le hL]
theorem q_m11_B_lo_eq (hL : LogLike L) : L (41476371717 / 10000000000 : ℝ) = 2 * L 2 + L (41476371717 / 40000000000 : ℝ) := by
  have h := hL.two_pow_mul 2 (r := (41476371717 / 40000000000 : ℝ)) (by norm_num)
  rw [show (2 : ℝ) ^ 2 * (41476371717 / 40000000000 : ℝ) = (41476371717 / 10000000000 : ℝ) by norm_num] at h
  exact_mod_cast h
theorem q_m11_B_lo_ge (hL : LogLike L) : (1778173519851 / 1250000000000 : ℝ) ≤ L (41476371717 / 10000000000 : ℝ) := by
  have h : (453055684513 / 12500000000000 : ℝ) ≤ L (41476371717 / 40000000000 : ℝ) := hL.ge_of (q := (41476371717 / 40000000000 : ℝ)) (a := (1476371717 / 40000000000 : ℝ)) (by norm_num) (by norm_num) 4 (by norm_num [lsum, Finset.sum_range_succ])
  rw [q_m11_B_lo_eq hL]; linarith [L2_ge hL]
theorem q_m11_B_lo_le (hL : LogLike L) : L (41476371717 / 10000000000 : ℝ) ≤ (14225388158811 / 10000000000000 : ℝ) := by
  have h : L (41476371717 / 40000000000 : ℝ) ≤ (724889095221 / 20000000000000 : ℝ) := hL.le_of (q := (41476371717 / 40000000000 : ℝ)) (a := (1476371717 / 40000000000 : ℝ)) (by norm_num) (by norm_num) 4 (by norm_num [lsum, ltail, Finset.sum_range_succ])
  rw [q_m11_B_lo_eq hL]; linarith [L2_le hL]
theorem q_m11_B_hi_eq (hL : LogLike L) : L (20738185859 / 5000000000 : ℝ) = 2 * L 2 + L (20738185859 / 20000000000 : ℝ) := by
  have h := hL.two_pow_mul 2 (r := (20738185859 / 20000000000 : ℝ)) (by norm_num)
  rw [show (2 : ℝ) ^ 2 * (20738185859 / 20000000000 : ℝ) = (20738185859 / 5000000000 : ℝ) by norm_num] at h
  exact_mod_cast h
theorem q_m11_B_hi_ge (hL : LogLike L) : (14225388159049 / 10000000000000 : ℝ) ≤ L (20738185859 / 5000000000 : ℝ) := by
  have h : (724889095703 / 20000000000000 : ℝ) ≤ L (20738185859 / 20000000000 : ℝ) := hL.ge_of (q := (20738185859 / 20000000000 : ℝ)) (a := (738185859 / 20000000000 : ℝ)) (by norm_num) (by norm_num) 4 (by norm_num [lsum, Finset.sum_range_succ])
  rw [q_m11_B_hi_eq hL]; linarith [L2_ge hL]
theorem q_m11_B_hi_le (hL : LogLike L) : L (20738185859 / 5000000000 : ℝ) ≤ (3556347039763 / 2500000000000 : ℝ) := by
  have h : L (20738185859 / 20000000000 : ℝ) ≤ (906111369629 / 25000000000000 : ℝ) := hL.le_of (q := (20738185859 / 20000000000 : ℝ)) (a := (738185859 / 20000000000 : ℝ)) (by norm_num) (by norm_num) 4 (by norm_num [lsum, ltail, Finset.sum_range_succ])
  rw [q_m11_B_hi_eq hL]; linarith [L2_le hL]
theorem LLLm11_ge (hL : LogLike L) : (7112694079 / 5000000000 : ℝ) ≤ L (L (L (m11 : ℝ))) := by
  have h := hL.mono _ _ (by norm_num) (LLm11_ge hL)
  linarith [q_m11_B_lo_ge hL]
theorem LLLm11_le (hL : LogLike L) : L (L (L (m11 : ℝ))) ≤ (22227169 / 15625000 : ℝ) := by
  have h := hL.mono _ _ (by linarith [LLm11_ge hL]) (LLm11_le hL)
  linarith [q_m11_B_hi_le hL]
theorem q_m11_Q_lo_eq (hL : LogLike L) : L (1845710089369 / 500000000000 : ℝ) = 1 * L 2 + L (1845710089369 / 1000000000000 : ℝ) := by
  have h := hL.two_pow_mul 1 (r := (1845710089369 / 1000000000000 : ℝ)) (by norm_num)
  rw [show (2 : ℝ) ^ 1 * (1845710089369 / 1000000000000 : ℝ) = (1845710089369 / 500000000000 : ℝ) by norm_num] at h
  exact_mod_cast h
theorem q_m11_Q_lo_ge (hL : LogLike L) : (6530056281473 / 5000000000000 : ℝ) ≤ L (1845710089369 / 500000000000 : ℝ) := by
  have h : (61286407573479 / 100000000000000 : ℝ) ≤ L (1845710089369 / 1000000000000 : ℝ) := hL.ge_of (q := (1845710089369 / 1000000000000 : ℝ)) (a := (845710089369 / 1000000000000 : ℝ)) (by norm_num) (by norm_num) 11 (by norm_num [lsum, Finset.sum_range_succ])
  rw [q_m11_Q_lo_eq hL]; linarith [L2_ge hL]
theorem q_m11_Q_lo_le (hL : LogLike L) : L (1845710089369 / 500000000000 : ℝ) ≤ (13060112562949 / 10000000000000 : ℝ) := by
  have h : L (1845710089369 / 1000000000000 : ℝ) ≤ (3830400473343 / 6250000000000 : ℝ) := hL.le_of (q := (1845710089369 / 1000000000000 : ℝ)) (a := (845710089369 / 1000000000000 : ℝ)) (by norm_num) (by norm_num) 11 (by norm_num [lsum, ltail, Finset.sum_range_succ])
  rw [q_m11_Q_lo_eq hL]; linarith [L2_le hL]
theorem q_m11_Q_hi_eq (hL : LogLike L) : L (3691420178739 / 1000000000000 : ℝ) = 1 * L 2 + L (3691420178739 / 2000000000000 : ℝ) := by
  have h := hL.two_pow_mul 1 (r := (3691420178739 / 2000000000000 : ℝ)) (by norm_num)
  rw [show (2 : ℝ) ^ 1 * (3691420178739 / 2000000000000 : ℝ) = (3691420178739 / 1000000000000 : ℝ) by norm_num] at h
  exact_mod_cast h
theorem q_m11_Q_hi_ge (hL : LogLike L) : (13060112562949 / 10000000000000 : ℝ) ≤ L (3691420178739 / 1000000000000 : ℝ) := by
  have h : (30643203786753 / 50000000000000 : ℝ) ≤ L (3691420178739 / 2000000000000 : ℝ) := hL.ge_of (q := (3691420178739 / 2000000000000 : ℝ)) (a := (1691420178739 / 2000000000000 : ℝ)) (by norm_num) (by norm_num) 11 (by norm_num [lsum, Finset.sum_range_succ])
  rw [q_m11_Q_hi_eq hL]; linarith [L2_ge hL]
theorem q_m11_Q_hi_le (hL : LogLike L) : L (3691420178739 / 1000000000000 : ℝ) ≤ (1632514070369 / 1250000000000 : ℝ) := by
  have h : L (3691420178739 / 2000000000000 : ℝ) ≤ (12257281514703 / 20000000000000 : ℝ) := hL.le_of (q := (3691420178739 / 2000000000000 : ℝ)) (a := (1691420178739 / 2000000000000 : ℝ)) (by norm_num) (by norm_num) 11 (by norm_num [lsum, ltail, Finset.sum_range_succ])
  rw [q_m11_Q_hi_eq hL]; linarith [L2_le hL]
theorem LQm11_le (hL : LogLike L) : L ((σ 1 m11 : ℝ) / (m11 : ℝ)) ≤ (1632514070369 / 1250000000000 : ℝ) := by
  rw [sig_m11]
  have h1 : ((11251566131001478449463296000 : ℕ) : ℝ) / (m11 : ℝ) ≤ (3691420178739 / 1000000000000 : ℝ) := by norm_num [m11]
  have h2 := hL.mono _ _ (by norm_num [m11]) h1
  linarith [q_m11_Q_hi_le hL]
theorem LQm11_ge (hL : LogLike L) : (6530056281473 / 5000000000000 : ℝ) ≤ L ((σ 1 m11 : ℝ) / (m11 : ℝ)) := by
  rw [sig_m11]
  have h1 : (1845710089369 / 500000000000 : ℝ) ≤ ((11251566131001478449463296000 : ℕ) : ℝ) / (m11 : ℝ) := by norm_num [m11]
  have h2 := hL.mono _ _ (by norm_num) h1
  linarith [q_m11_Q_lo_ge hL]

end OddRobin
