import Riemann.Deep.OddRobin.Lean.Cost

namespace OddRobin

open ArithmeticFunction
open scoped ArithmeticFunction.sigma

variable {L : ℝ → ℝ}

theorem L2_ge (hL : LogLike L) : (6931471805599 / 10000000000000 : ℝ) ≤ L 2 := hL.ge_of (q := (2 : ℝ)) (a := (1 : ℝ)) (by norm_num) (by norm_num) 14 (by norm_num [lsum, Finset.sum_range_succ])
theorem L2_le (hL : LogLike L) : L 2 ≤ (8664339757 / 12500000000 : ℝ) := hL.le_of (q := (2 : ℝ)) (a := (1 : ℝ)) (by norm_num) (by norm_num) 14 (by norm_num [lsum, ltail, Finset.sum_range_succ])
theorem L3_eq (hL : LogLike L) : L (3 : ℝ) = 1 * L 2 + L (3 / 2 : ℝ) := by
  have h := hL.two_pow_mul 1 (r := (3 / 2 : ℝ)) (by norm_num)
  rw [show (2 : ℝ) ^ 1 * (3 / 2 : ℝ) = 3 by norm_num] at h
  exact_mod_cast h
theorem L3_ge (hL : LogLike L) : (274653072167 / 250000000000 : ℝ) ≤ L (3 : ℝ) := by
  have h : (8109302162163 / 20000000000000 : ℝ) ≤ L (3 / 2 : ℝ) := hL.ge_of (q := (3 / 2 : ℝ)) (a := (1 / 2 : ℝ)) (by norm_num) (by norm_num) 9 (by norm_num [lsum, Finset.sum_range_succ])
  rw [L3_eq hL]; linarith [L2_ge hL]
theorem L3_le (hL : LogLike L) : L (3 : ℝ) ≤ (1098612288669 / 1000000000000 : ℝ) := by
  have h : L (3 / 2 : ℝ) ≤ (40546510810817 / 100000000000000 : ℝ) := hL.le_of (q := (3 / 2 : ℝ)) (a := (1 / 2 : ℝ)) (by norm_num) (by norm_num) 9 (by norm_num [lsum, ltail, Finset.sum_range_succ])
  rw [L3_eq hL]; linarith [L2_le hL]
theorem L5_eq (hL : LogLike L) : L (5 : ℝ) = 2 * L 2 + L (5 / 4 : ℝ) := by
  have h := hL.two_pow_mul 2 (r := (5 / 4 : ℝ)) (by norm_num)
  rw [show (2 : ℝ) ^ 2 * (5 / 4 : ℝ) = 5 by norm_num] at h
  exact_mod_cast h
theorem L5_ge (hL : LogLike L) : (1609437912433 / 1000000000000 : ℝ) ≤ L (5 : ℝ) := by
  have h : (11157177565707 / 50000000000000 : ℝ) ≤ L (5 / 4 : ℝ) := hL.ge_of (q := (5 / 4 : ℝ)) (a := (1 / 4 : ℝ)) (by norm_num) (by norm_num) 6 (by norm_num [lsum, Finset.sum_range_succ])
  rw [L5_eq hL]; linarith [L2_ge hL]
theorem L5_le (hL : LogLike L) : L (5 : ℝ) ≤ (321887582487 / 200000000000 : ℝ) := by
  have h : L (5 / 4 : ℝ) ≤ (22314355131421 / 100000000000000 : ℝ) := hL.le_of (q := (5 / 4 : ℝ)) (a := (1 / 4 : ℝ)) (by norm_num) (by norm_num) 6 (by norm_num [lsum, ltail, Finset.sum_range_succ])
  rw [L5_eq hL]; linarith [L2_le hL]
theorem L7_eq (hL : LogLike L) : L (7 : ℝ) = 2 * L 2 + L (7 / 4 : ℝ) := by
  have h := hL.two_pow_mul 2 (r := (7 / 4 : ℝ)) (by norm_num)
  rw [show (2 : ℝ) ^ 2 * (7 / 4 : ℝ) = 7 by norm_num] at h
  exact_mod_cast h
theorem L7_ge (hL : LogLike L) : (389182029811 / 200000000000 : ℝ) ≤ L (7 : ℝ) := by
  have h : (55961578793541 / 100000000000000 : ℝ) ≤ L (7 / 4 : ℝ) := hL.ge_of (q := (7 / 4 : ℝ)) (a := (3 / 4 : ℝ)) (by norm_num) (by norm_num) 11 (by norm_num [lsum, Finset.sum_range_succ])
  rw [L7_eq hL]; linarith [L2_ge hL]
theorem L7_le (hL : LogLike L) : L (7 : ℝ) ≤ (30404846079 / 15625000000 : ℝ) := by
  have h : L (7 / 4 : ℝ) ≤ (55961578793543 / 100000000000000 : ℝ) := hL.le_of (q := (7 / 4 : ℝ)) (a := (3 / 4 : ℝ)) (by norm_num) (by norm_num) 11 (by norm_num [lsum, ltail, Finset.sum_range_succ])
  rw [L7_eq hL]; linarith [L2_le hL]
theorem L11_eq (hL : LogLike L) : L (11 : ℝ) = 3 * L 2 + L (11 / 8 : ℝ) := by
  have h := hL.two_pow_mul 3 (r := (11 / 8 : ℝ)) (by norm_num)
  rw [show (2 : ℝ) ^ 3 * (11 / 8 : ℝ) = 11 by norm_num] at h
  exact_mod_cast h
theorem L11_ge (hL : LogLike L) : (1198947636399 / 500000000000 : ℝ) ≤ L (11 : ℝ) := by
  have h : (31845373111853 / 100000000000000 : ℝ) ≤ L (11 / 8 : ℝ) := hL.ge_of (q := (11 / 8 : ℝ)) (a := (3 / 8 : ℝ)) (by norm_num) (by norm_num) 8 (by norm_num [lsum, Finset.sum_range_succ])
  rw [L11_eq hL]; linarith [L2_ge hL]
theorem L11_le (hL : LogLike L) : L (11 : ℝ) ≤ (2397895272799 / 1000000000000 : ℝ) := by
  have h : L (11 / 8 : ℝ) ≤ (15922686555927 / 50000000000000 : ℝ) := hL.le_of (q := (11 / 8 : ℝ)) (a := (3 / 8 : ℝ)) (by norm_num) (by norm_num) 8 (by norm_num [lsum, ltail, Finset.sum_range_succ])
  rw [L11_eq hL]; linarith [L2_le hL]
theorem L13_eq (hL : LogLike L) : L (13 : ℝ) = 3 * L 2 + L (13 / 8 : ℝ) := by
  have h := hL.two_pow_mul 3 (r := (13 / 8 : ℝ)) (by norm_num)
  rw [show (2 : ℝ) ^ 3 * (13 / 8 : ℝ) = 13 by norm_num] at h
  exact_mod_cast h
theorem L13_ge (hL : LogLike L) : (2564949357461 / 1000000000000 : ℝ) ≤ L (13 : ℝ) := by
  have h : (48550781578169 / 100000000000000 : ℝ) ≤ L (13 / 8 : ℝ) := hL.ge_of (q := (13 / 8 : ℝ)) (a := (5 / 8 : ℝ)) (by norm_num) (by norm_num) 10 (by norm_num [lsum, Finset.sum_range_succ])
  rw [L13_eq hL]; linarith [L2_ge hL]
theorem L13_le (hL : LogLike L) : L (13 : ℝ) ≤ (1282474678731 / 500000000000 : ℝ) := by
  have h : L (13 / 8 : ℝ) ≤ (48550781578171 / 100000000000000 : ℝ) := hL.le_of (q := (13 / 8 : ℝ)) (a := (5 / 8 : ℝ)) (by norm_num) (by norm_num) 10 (by norm_num [lsum, ltail, Finset.sum_range_succ])
  rw [L13_eq hL]; linarith [L2_le hL]
theorem L61_eq (hL : LogLike L) : L (61 : ℝ) = 5 * L 2 + L (61 / 32 : ℝ) := by
  have h := hL.two_pow_mul 5 (r := (61 / 32 : ℝ)) (by norm_num)
  rw [show (2 : ℝ) ^ 5 * (61 / 32 : ℝ) = 61 by norm_num] at h
  exact_mod_cast h
theorem L61_ge (hL : LogLike L) : (4110873864173 / 1000000000000 : ℝ) ≤ L (61 : ℝ) := by
  have h : (16128449034339 / 25000000000000 : ℝ) ≤ L (61 / 32 : ℝ) := hL.ge_of (q := (61 / 32 : ℝ)) (a := (29 / 32 : ℝ)) (by norm_num) (by norm_num) 12 (by norm_num [lsum, Finset.sum_range_succ])
  rw [L61_eq hL]; linarith [L2_ge hL]
theorem L61_le (hL : LogLike L) : L (61 : ℝ) ≤ (2055436932087 / 500000000000 : ℝ) := by
  have h : L (61 / 32 : ℝ) ≤ (64513796137359 / 100000000000000 : ℝ) := hL.le_of (q := (61 / 32 : ℝ)) (a := (29 / 32 : ℝ)) (by norm_num) (by norm_num) 12 (by norm_num [lsum, ltail, Finset.sum_range_succ])
  rw [L61_eq hL]; linarith [L2_le hL]
theorem L67_eq (hL : LogLike L) : L (67 : ℝ) = 6 * L 2 + L (67 / 64 : ℝ) := by
  have h := hL.two_pow_mul 6 (r := (67 / 64 : ℝ)) (by norm_num)
  rw [show (2 : ℝ) ^ 6 * (67 / 64 : ℝ) = 67 by norm_num] at h
  exact_mod_cast h
theorem L67_ge (hL : LogLike L) : (420469261939 / 100000000000 : ℝ) ≤ L (67 : ℝ) := by
  have h : (4580953603129 / 100000000000000 : ℝ) ≤ L (67 / 64 : ℝ) := hL.ge_of (q := (67 / 64 : ℝ)) (a := (3 / 64 : ℝ)) (by norm_num) (by norm_num) 4 (by norm_num [lsum, Finset.sum_range_succ])
  rw [L67_eq hL]; linarith [L2_ge hL]
theorem L67_le (hL : LogLike L) : L (67 : ℝ) ≤ (32849161089 / 7812500000 : ℝ) := by
  have h : L (67 / 64 : ℝ) ≤ (458095360313 / 10000000000000 : ℝ) := hL.le_of (q := (67 / 64 : ℝ)) (a := (3 / 64 : ℝ)) (by norm_num) (by norm_num) 4 (by norm_num [lsum, ltail, Finset.sum_range_succ])
  rw [L67_eq hL]; linarith [L2_le hL]
theorem L71_eq (hL : LogLike L) : L (71 : ℝ) = 6 * L 2 + L (71 / 64 : ℝ) := by
  have h := hL.two_pow_mul 6 (r := (71 / 64 : ℝ)) (by norm_num)
  rw [show (2 : ℝ) ^ 6 * (71 / 64 : ℝ) = 71 by norm_num] at h
  exact_mod_cast h
theorem L71_ge (hL : LogLike L) : (4262679877041 / 1000000000000 : ℝ) ≤ L (71 : ℝ) := by
  have h : (2594919842041 / 25000000000000 : ℝ) ≤ L (71 / 64 : ℝ) := hL.ge_of (q := (71 / 64 : ℝ)) (a := (7 / 64 : ℝ)) (by norm_num) (by norm_num) 5 (by norm_num [lsum, Finset.sum_range_succ])
  rw [L71_eq hL]; linarith [L2_ge hL]
theorem L71_le (hL : LogLike L) : L (71 : ℝ) ≤ (2131339938521 / 500000000000 : ℝ) := by
  have h : L (71 / 64 : ℝ) ≤ (2075935873633 / 20000000000000 : ℝ) := hL.le_of (q := (71 / 64 : ℝ)) (a := (7 / 64 : ℝ)) (by norm_num) (by norm_num) 5 (by norm_num [lsum, ltail, Finset.sum_range_succ])
  rw [L71_eq hL]; linarith [L2_le hL]
theorem L73_eq (hL : LogLike L) : L (73 : ℝ) = 6 * L 2 + L (73 / 64 : ℝ) := by
  have h := hL.two_pow_mul 6 (r := (73 / 64 : ℝ)) (by norm_num)
  rw [show (2 : ℝ) ^ 6 * (73 / 64 : ℝ) = 73 by norm_num] at h
  exact_mod_cast h
theorem L73_ge (hL : LogLike L) : (1072614860287 / 250000000000 : ℝ) ≤ L (73 : ℝ) := by
  have h : (1315763577887 / 10000000000000 : ℝ) ≤ L (73 / 64 : ℝ) := hL.ge_of (q := (73 / 64 : ℝ)) (a := (9 / 64 : ℝ)) (by norm_num) (by norm_num) 5 (by norm_num [lsum, Finset.sum_range_succ])
  rw [L73_eq hL]; linarith [L2_ge hL]
theorem L73_le (hL : LogLike L) : L (73 : ℝ) ≤ (4290459441149 / 1000000000000 : ℝ) := by
  have h : L (73 / 64 : ℝ) ≤ (1644704472359 / 12500000000000 : ℝ) := hL.le_of (q := (73 / 64 : ℝ)) (a := (9 / 64 : ℝ)) (by norm_num) (by norm_num) 5 (by norm_num [lsum, ltail, Finset.sum_range_succ])
  rw [L73_eq hL]; linarith [L2_le hL]
theorem L6867_ge (hL : LogLike L) : (74075428925703 / 5000000000000000 : ℝ) ≤ L (68 / 67 : ℝ) := hL.ge_of (q := (68 / 67 : ℝ)) (a := (1 / 67 : ℝ)) (by norm_num) (by norm_num) 4 (by norm_num [lsum, Finset.sum_range_succ])
theorem L6867_le (hL : LogLike L) : L (68 / 67 : ℝ) ≤ (148150857851407 / 10000000000000000 : ℝ) := hL.le_of (q := (68 / 67 : ℝ)) (a := (1 / 67 : ℝ)) (by norm_num) (by norm_num) 4 (by norm_num [lsum, ltail, Finset.sum_range_succ])
theorem L67_pos (hL : LogLike L) : 0 < L (67 : ℝ) := by linarith [L67_ge hL]
theorem eps_ge (hL : LogLike L) : (352346464443 / 100000000000000 : ℝ) ≤ eps L := by
  unfold eps; rw [le_div_iff₀ (L67_pos hL)]; linarith [L6867_ge hL, L67_le hL]
theorem eps_le (hL : LogLike L) : eps L ≤ (88086616111 / 25000000000000 : ℝ) := by
  unfold eps; rw [div_le_iff₀ (L67_pos hL)]; linarith [L6867_le hL, L67_ge hL]
theorem eps_pos (hL : LogLike L) : 0 < eps L := by linarith [eps_ge hL]

end OddRobin
