import Riemann.Deep.OddRobin.Lean.NumA
import Riemann.Deep.OddRobin.Lean.SigVals

namespace OddRobin

open ArithmeticFunction
open scoped ArithmeticFunction.sigma

variable {L : ℝ → ℝ}

theorem Lm3_ge (hL : LogLike L) : (61985010601 / 1000000000 : ℝ) ≤ L (m3 : ℝ) := by
  have h1 : (2 : ℝ) ^ 89 * (1678759423931 / 1250000000000 : ℝ) ≤ (m3 : ℝ) := by norm_num [m3]
  have h2 := hL.mono _ _ (by norm_num) h1
  rw [hL.two_pow_mul 89 (by norm_num)] at h2
  have h3 : (29491153118487 / 100000000000000 : ℝ) ≤ L (1678759423931 / 1250000000000 : ℝ) := hL.ge_of (q := (1678759423931 / 1250000000000 : ℝ)) (a := (428759423931 / 1250000000000 : ℝ)) (by norm_num) (by norm_num) 7 (by norm_num [lsum, Finset.sum_range_succ])
  push_cast at h2; linarith [L2_ge hL]
theorem Lm3_le (hL : LogLike L) : L (m3 : ℝ) ≤ (619850106011 / 10000000000 : ℝ) := by
  have h1 : (m3 : ℝ) ≤ (2 : ℝ) ^ 89 * (13430075391449 / 10000000000000 : ℝ) := by norm_num [m3]
  have h2 := hL.mono _ _ (by norm_num [m3]) h1
  rw [hL.two_pow_mul 89 (by norm_num)] at h2
  have h3 : L (13430075391449 / 10000000000000 : ℝ) ≤ (29491153118499 / 100000000000000 : ℝ) := hL.le_of (q := (13430075391449 / 10000000000000 : ℝ)) (a := (3430075391449 / 10000000000000 : ℝ)) (by norm_num) (by norm_num) 7 (by norm_num [lsum, ltail, Finset.sum_range_succ])
  push_cast at h2; linarith [L2_le hL]
theorem q_m3_A_lo_eq (hL : LogLike L) : L (61985010601 / 1000000000 : ℝ) = 5 * L 2 + L (61985010601 / 32000000000 : ℝ) := by
  have h := hL.two_pow_mul 5 (r := (61985010601 / 32000000000 : ℝ)) (by norm_num)
  rw [show (2 : ℝ) ^ 5 * (61985010601 / 32000000000 : ℝ) = (61985010601 / 1000000000 : ℝ) by norm_num] at h
  exact_mod_cast h
theorem q_m3_A_lo_ge (hL : LogLike L) : (825378518263 / 200000000000 : ℝ) ≤ L (61985010601 / 1000000000 : ℝ) := by
  have h : (33057834425779 / 50000000000000 : ℝ) ≤ L (61985010601 / 32000000000 : ℝ) := hL.ge_of (q := (61985010601 / 32000000000 : ℝ)) (a := (29985010601 / 32000000000 : ℝ)) (by norm_num) (by norm_num) 12 (by norm_num [lsum, Finset.sum_range_succ])
  rw [q_m3_A_lo_eq hL]; linarith [L2_ge hL]
theorem q_m3_A_lo_le (hL : LogLike L) : L (61985010601 / 1000000000 : ℝ) ≤ (41268925913157 / 10000000000000 : ℝ) := by
  have h : L (61985010601 / 32000000000 : ℝ) ≤ (33057834425781 / 50000000000000 : ℝ) := hL.le_of (q := (61985010601 / 32000000000 : ℝ)) (a := (29985010601 / 32000000000 : ℝ)) (by norm_num) (by norm_num) 12 (by norm_num [lsum, ltail, Finset.sum_range_succ])
  rw [q_m3_A_lo_eq hL]; linarith [L2_le hL]
theorem q_m3_A_hi_eq (hL : LogLike L) : L (619850106011 / 10000000000 : ℝ) = 5 * L 2 + L (619850106011 / 320000000000 : ℝ) := by
  have h := hL.two_pow_mul 5 (r := (619850106011 / 320000000000 : ℝ)) (by norm_num)
  rw [show (2 : ℝ) ^ 5 * (619850106011 / 320000000000 : ℝ) = (619850106011 / 10000000000 : ℝ) by norm_num] at h
  exact_mod_cast h
theorem q_m3_A_hi_ge (hL : LogLike L) : (20634462956583 / 5000000000000 : ℝ) ≤ L (619850106011 / 10000000000 : ℝ) := by
  have h : (66115668851719 / 100000000000000 : ℝ) ≤ L (619850106011 / 320000000000 : ℝ) := hL.ge_of (q := (619850106011 / 320000000000 : ℝ)) (a := (299850106011 / 320000000000 : ℝ)) (by norm_num) (by norm_num) 12 (by norm_num [lsum, Finset.sum_range_succ])
  rw [q_m3_A_hi_eq hL]; linarith [L2_ge hL]
theorem q_m3_A_hi_le (hL : LogLike L) : L (619850106011 / 10000000000 : ℝ) ≤ (41268925913173 / 10000000000000 : ℝ) := by
  have h : L (619850106011 / 320000000000 : ℝ) ≤ (16528917212931 / 25000000000000 : ℝ) := hL.le_of (q := (619850106011 / 320000000000 : ℝ)) (a := (299850106011 / 320000000000 : ℝ)) (by norm_num) (by norm_num) 12 (by norm_num [lsum, ltail, Finset.sum_range_succ])
  rw [q_m3_A_hi_eq hL]; linarith [L2_le hL]
theorem LLm3_ge (hL : LogLike L) : (41268925913 / 10000000000 : ℝ) ≤ L (L (m3 : ℝ)) := by
  have h := hL.mono _ _ (by norm_num) (Lm3_ge hL)
  linarith [q_m3_A_lo_ge hL]
theorem LLm3_le (hL : LogLike L) : L (L (m3 : ℝ)) ≤ (20634462957 / 5000000000 : ℝ) := by
  have h := hL.mono _ _ (by linarith [Lm3_ge hL]) (Lm3_le hL)
  linarith [q_m3_A_hi_le hL]
theorem q_m3_B_lo_eq (hL : LogLike L) : L (41268925913 / 10000000000 : ℝ) = 2 * L 2 + L (41268925913 / 40000000000 : ℝ) := by
  have h := hL.two_pow_mul 2 (r := (41268925913 / 40000000000 : ℝ)) (by norm_num)
  rw [show (2 : ℝ) ^ 2 * (41268925913 / 40000000000 : ℝ) = (41268925913 / 10000000000 : ℝ) by norm_num] at h
  exact_mod_cast h
theorem q_m3_B_lo_ge (hL : LogLike L) : (7087623622887 / 5000000000000 : ℝ) ≤ L (41268925913 / 10000000000 : ℝ) := by
  have h : (3123036345767 / 100000000000000 : ℝ) ≤ L (41268925913 / 40000000000 : ℝ) := hL.ge_of (q := (41268925913 / 40000000000 : ℝ)) (a := (1268925913 / 40000000000 : ℝ)) (by norm_num) (by norm_num) 3 (by norm_num [lsum, Finset.sum_range_succ])
  rw [q_m3_B_lo_eq hL]; linarith [L2_ge hL]
theorem q_m3_B_lo_le (hL : LogLike L) : L (41268925913 / 10000000000 : ℝ) ≤ (7087623622889 / 5000000000000 : ℝ) := by
  have h : L (41268925913 / 40000000000 : ℝ) ≤ (124921453831 / 4000000000000 : ℝ) := hL.le_of (q := (41268925913 / 40000000000 : ℝ)) (a := (1268925913 / 40000000000 : ℝ)) (by norm_num) (by norm_num) 3 (by norm_num [lsum, ltail, Finset.sum_range_succ])
  rw [q_m3_B_lo_eq hL]; linarith [L2_le hL]
theorem q_m3_B_hi_eq (hL : LogLike L) : L (20634462957 / 5000000000 : ℝ) = 2 * L 2 + L (20634462957 / 20000000000 : ℝ) := by
  have h := hL.two_pow_mul 2 (r := (20634462957 / 20000000000 : ℝ)) (by norm_num)
  rw [show (2 : ℝ) ^ 2 * (20634462957 / 20000000000 : ℝ) = (20634462957 / 5000000000 : ℝ) by norm_num] at h
  exact_mod_cast h
theorem q_m3_B_hi_ge (hL : LogLike L) : (14175247246017 / 10000000000000 : ℝ) ≤ L (20634462957 / 5000000000 : ℝ) := by
  have h : (3123036348191 / 100000000000000 : ℝ) ≤ L (20634462957 / 20000000000 : ℝ) := hL.ge_of (q := (20634462957 / 20000000000 : ℝ)) (a := (634462957 / 20000000000 : ℝ)) (by norm_num) (by norm_num) 3 (by norm_num [lsum, Finset.sum_range_succ])
  rw [q_m3_B_hi_eq hL]; linarith [L2_ge hL]
theorem q_m3_B_hi_le (hL : LogLike L) : L (20634462957 / 5000000000 : ℝ) ≤ (708762362301 / 500000000000 : ℝ) := by
  have h : L (20634462957 / 20000000000 : ℝ) ≤ (1561518174099 / 50000000000000 : ℝ) := hL.le_of (q := (20634462957 / 20000000000 : ℝ)) (a := (634462957 / 20000000000 : ℝ)) (by norm_num) (by norm_num) 3 (by norm_num [lsum, ltail, Finset.sum_range_succ])
  rw [q_m3_B_hi_eq hL]; linarith [L2_le hL]
theorem LLLm3_ge (hL : LogLike L) : (2835049449 / 2000000000 : ℝ) ≤ L (L (L (m3 : ℝ))) := by
  have h := hL.mono _ _ (by norm_num) (LLm3_ge hL)
  linarith [q_m3_B_lo_ge hL]
theorem LLLm3_le (hL : LogLike L) : L (L (L (m3 : ℝ))) ≤ (14175247247 / 10000000000 : ℝ) := by
  have h := hL.mono _ _ (by linarith [LLm3_ge hL]) (LLm3_le hL)
  linarith [q_m3_B_hi_le hL]
theorem q_m3_Q_lo_eq (hL : LogLike L) : L (229609867577 / 62500000000 : ℝ) = 1 * L 2 + L (229609867577 / 125000000000 : ℝ) := by
  have h := hL.two_pow_mul 1 (r := (229609867577 / 125000000000 : ℝ)) (by norm_num)
  rw [show (2 : ℝ) ^ 1 * (229609867577 / 125000000000 : ℝ) = (229609867577 / 62500000000 : ℝ) by norm_num] at h
  exact_mod_cast h
theorem q_m3_Q_lo_ge (hL : LogLike L) : (13012150840311 / 10000000000000 : ℝ) ≤ L (229609867577 / 62500000000 : ℝ) := by
  have h : (60806790347129 / 100000000000000 : ℝ) ≤ L (229609867577 / 125000000000 : ℝ) := hL.ge_of (q := (229609867577 / 125000000000 : ℝ)) (a := (104609867577 / 125000000000 : ℝ)) (by norm_num) (by norm_num) 11 (by norm_num [lsum, Finset.sum_range_succ])
  rw [q_m3_Q_lo_eq hL]; linarith [L2_ge hL]
theorem q_m3_Q_lo_le (hL : LogLike L) : L (229609867577 / 62500000000 : ℝ) ≤ (6506075420157 / 5000000000000 : ℝ) := by
  have h : L (229609867577 / 125000000000 : ℝ) ≤ (475053049587 / 781250000000 : ℝ) := hL.le_of (q := (229609867577 / 125000000000 : ℝ)) (a := (104609867577 / 125000000000 : ℝ)) (by norm_num) (by norm_num) 11 (by norm_num [lsum, ltail, Finset.sum_range_succ])
  rw [q_m3_Q_lo_eq hL]; linarith [L2_le hL]
theorem q_m3_Q_hi_eq (hL : LogLike L) : L (3673757881233 / 1000000000000 : ℝ) = 1 * L 2 + L (3673757881233 / 2000000000000 : ℝ) := by
  have h := hL.two_pow_mul 1 (r := (3673757881233 / 2000000000000 : ℝ)) (by norm_num)
  rw [show (2 : ℝ) ^ 1 * (3673757881233 / 2000000000000 : ℝ) = (3673757881233 / 1000000000000 : ℝ) by norm_num] at h
  exact_mod_cast h
theorem q_m3_Q_hi_ge (hL : LogLike L) : (6506075420157 / 5000000000000 : ℝ) ≤ L (3673757881233 / 1000000000000 : ℝ) := by
  have h : (15201697586789 / 25000000000000 : ℝ) ≤ L (3673757881233 / 2000000000000 : ℝ) := hL.ge_of (q := (3673757881233 / 2000000000000 : ℝ)) (a := (1673757881233 / 2000000000000 : ℝ)) (by norm_num) (by norm_num) 11 (by norm_num [lsum, Finset.sum_range_succ])
  rw [q_m3_Q_hi_eq hL]; linarith [L2_ge hL]
theorem q_m3_Q_hi_le (hL : LogLike L) : L (3673757881233 / 1000000000000 : ℝ) ≤ (13012150840317 / 10000000000000 : ℝ) := by
  have h : L (3673757881233 / 2000000000000 : ℝ) ≤ (60806790347163 / 100000000000000 : ℝ) := hL.le_of (q := (3673757881233 / 2000000000000 : ℝ)) (a := (1673757881233 / 2000000000000 : ℝ)) (by norm_num) (by norm_num) 11 (by norm_num [lsum, ltail, Finset.sum_range_succ])
  rw [q_m3_Q_hi_eq hL]; linarith [L2_le hL]
theorem LQm3_le (hL : LogLike L) : L ((σ 1 m3 : ℝ) / (m3 : ℝ)) ≤ (13012150840317 / 10000000000000 : ℝ) := by
  rw [sig_m3]
  have h1 : ((3053926605369692280324096000 : ℕ) : ℝ) / (m3 : ℝ) ≤ (3673757881233 / 1000000000000 : ℝ) := by norm_num [m3]
  have h2 := hL.mono _ _ (by norm_num [m3]) h1
  linarith [q_m3_Q_hi_le hL]
theorem LQm3_ge (hL : LogLike L) : (13012150840311 / 10000000000000 : ℝ) ≤ L ((σ 1 m3 : ℝ) / (m3 : ℝ)) := by
  rw [sig_m3]
  have h1 : (229609867577 / 62500000000 : ℝ) ≤ ((3053926605369692280324096000 : ℕ) : ℝ) / (m3 : ℝ) := by norm_num [m3]
  have h2 := hL.mono _ _ (by norm_num) h1
  linarith [q_m3_Q_lo_ge hL]

end OddRobin
