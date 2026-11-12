import Riemann.Deep.OddRobin.Lean.NumA

namespace OddRobin

open ArithmeticFunction
open scoped ArithmeticFunction.sigma

variable {L : ℝ → ℝ}

theorem sig_3_3 : σ 1 (3 ^ 3) = 40 := by
  rw [sigma_one_apply_prime_pow (by norm_num)]; norm_num [Finset.sum_range_succ]
theorem inc_3_3_eq : inc L 3 3 = L (121 / 120 : ℝ) - eps L * L (3 : ℝ) := by
  rw [inc_val 3 3 40 sig_3_3]; norm_num
theorem inc_3_3_ge (hL : LogLike L) : (221394062881 / 50000000000000 : ℝ) ≤ inc L 3 3 := by
  have h1 : (1659760562939 / 200000000000000 : ℝ) ≤ L (121 / 120 : ℝ) := hL.ge_of (q := (121 / 120 : ℝ)) (a := (1 / 120 : ℝ)) (by norm_num) (by norm_num) 3 (by norm_num [lsum, Finset.sum_range_succ])
  have h2 : eps L * L (3 : ℝ) ≤ (88086616111 / 25000000000000 : ℝ) * (1098612288669 / 1000000000000 : ℝ) :=
    mul_le_mul (eps_le hL) (L3_le hL) (by linarith [L3_ge hL]) (by norm_num)
  rw [inc_3_3_eq]; linarith
theorem inc_3_3_le (hL : LogLike L) : inc L 3 3 ≤ (110697031441 / 25000000000000 : ℝ) := by
  have h1 : L (121 / 120 : ℝ) ≤ (82988028146951 / 10000000000000000 : ℝ) := hL.le_of (q := (121 / 120 : ℝ)) (a := (1 / 120 : ℝ)) (by norm_num) (by norm_num) 3 (by norm_num [lsum, ltail, Finset.sum_range_succ])
  have h2 : (352346464443 / 100000000000000 : ℝ) * (274653072167 / 250000000000 : ℝ) ≤ eps L * L (3 : ℝ) :=
    mul_le_mul (eps_ge hL) (L3_ge hL) (by norm_num) (eps_pos hL).le
  rw [inc_3_3_eq]; linarith
theorem sig_3_4 : σ 1 (3 ^ 4) = 121 := by
  rw [sigma_one_apply_prime_pow (by norm_num)]; norm_num [Finset.sum_range_succ]
theorem inc_3_4_eq : inc L 3 4 = L (364 / 363 : ℝ) - eps L * L (3 : ℝ) := by
  rw [inc_val 3 4 121 sig_3_4]; norm_num
theorem inc_3_4_ge (hL : LogLike L) : (-111988818519 / 100000000000000 : ℝ) ≤ inc L 3 4 := by
  have h1 : (13755166859449 / 5000000000000000 : ℝ) ≤ L (364 / 363 : ℝ) := hL.ge_of (q := (364 / 363 : ℝ)) (a := (1 / 363 : ℝ)) (by norm_num) (by norm_num) 3 (by norm_num [lsum, Finset.sum_range_succ])
  have h2 : eps L * L (3 : ℝ) ≤ (88086616111 / 25000000000000 : ℝ) * (1098612288669 / 1000000000000 : ℝ) :=
    mul_le_mul (eps_le hL) (L3_le hL) (by linarith [L3_ge hL]) (by norm_num)
  rw [inc_3_4_eq]; linarith
theorem inc_3_4_le (hL : LogLike L) : inc L 3 4 ≤ (-27997204629 / 25000000000000 : ℝ) := by
  have h1 : L (364 / 363 : ℝ) ≤ (27510333718899 / 10000000000000000 : ℝ) := hL.le_of (q := (364 / 363 : ℝ)) (a := (1 / 363 : ℝ)) (by norm_num) (by norm_num) 3 (by norm_num [lsum, ltail, Finset.sum_range_succ])
  have h2 : (352346464443 / 100000000000000 : ℝ) * (274653072167 / 250000000000 : ℝ) ≤ eps L * L (3 : ℝ) :=
    mul_le_mul (eps_ge hL) (L3_ge hL) (by norm_num) (eps_pos hL).le
  rw [inc_3_4_eq]; linarith
theorem sig_3_5 : σ 1 (3 ^ 5) = 364 := by
  rw [sigma_one_apply_prime_pow (by norm_num)]; norm_num [Finset.sum_range_succ]
theorem inc_3_5_eq : inc L 3 5 = L (1093 / 1092 : ℝ) - eps L * L (3 : ℝ) := by
  rw [inc_val 3 5 364 sig_3_5]; norm_num
theorem inc_3_5_ge (hL : LogLike L) : (-295558968539 / 100000000000000 : ℝ) ≤ inc L 3 5 := by
  have h1 : (4576659358441 / 5000000000000000 : ℝ) ≤ L (1093 / 1092 : ℝ) := hL.ge_of (q := (1093 / 1092 : ℝ)) (a := (1 / 1092 : ℝ)) (by norm_num) (by norm_num) 2 (by norm_num [lsum, Finset.sum_range_succ])
  have h2 : eps L * L (3 : ℝ) ≤ (88086616111 / 25000000000000 : ℝ) * (1098612288669 / 1000000000000 : ℝ) :=
    mul_le_mul (eps_le hL) (L3_le hL) (by linarith [L3_ge hL]) (by norm_num)
  rw [inc_3_5_eq]; linarith
theorem inc_3_5_le (hL : LogLike L) : inc L 3 5 ≤ (-36944871067 / 12500000000000 : ℝ) := by
  have h1 : L (1093 / 1092 : ℝ) ≤ (9153318716883 / 10000000000000000 : ℝ) := hL.le_of (q := (1093 / 1092 : ℝ)) (a := (1 / 1092 : ℝ)) (by norm_num) (by norm_num) 2 (by norm_num [lsum, ltail, Finset.sum_range_succ])
  have h2 : (352346464443 / 100000000000000 : ℝ) * (274653072167 / 250000000000 : ℝ) ≤ eps L * L (3 : ℝ) :=
    mul_le_mul (eps_ge hL) (L3_ge hL) (by norm_num) (eps_pos hL).le
  rw [inc_3_5_eq]; linarith
theorem sig_5_1 : σ 1 (5 ^ 1) = 6 := by
  rw [sigma_one_apply_prime_pow (by norm_num)]; norm_num [Finset.sum_range_succ]
theorem inc_5_1_eq : inc L 5 1 = L (31 / 30 : ℝ) - eps L * L (5 : ℝ) := by
  rw [inc_val 5 1 6 sig_5_1]; norm_num
theorem inc_5_1_ge (hL : LogLike L) : (271190252411 / 10000000000000 : ℝ) ≤ inc L 5 1 := by
  have h1 : (81974557057477 / 2500000000000000 : ℝ) ≤ L (31 / 30 : ℝ) := hL.ge_of (q := (31 / 30 : ℝ)) (a := (1 / 30 : ℝ)) (by norm_num) (by norm_num) 4 (by norm_num [lsum, Finset.sum_range_succ])
  have h2 : eps L * L (5 : ℝ) ≤ (88086616111 / 25000000000000 : ℝ) * (321887582487 / 200000000000 : ℝ) :=
    mul_le_mul (eps_le hL) (L5_le hL) (by linarith [L5_ge hL]) (by norm_num)
  rw [inc_5_1_eq]; linarith
theorem inc_5_1_le (hL : LogLike L) : inc L 5 1 ≤ (2711902524113 / 100000000000000 : ℝ) := by
  have h1 : L (31 / 30 : ℝ) ≤ (327898228229909 / 10000000000000000 : ℝ) := hL.le_of (q := (31 / 30 : ℝ)) (a := (1 / 30 : ℝ)) (by norm_num) (by norm_num) 4 (by norm_num [lsum, ltail, Finset.sum_range_succ])
  have h2 : (352346464443 / 100000000000000 : ℝ) * (1609437912433 / 1000000000000 : ℝ) ≤ eps L * L (5 : ℝ) :=
    mul_le_mul (eps_ge hL) (L5_ge hL) (by norm_num) (eps_pos hL).le
  rw [inc_5_1_eq]; linarith
theorem sig_5_2 : σ 1 (5 ^ 2) = 31 := by
  rw [sigma_one_apply_prime_pow (by norm_num)]; norm_num [Finset.sum_range_succ]
theorem inc_5_2_eq : inc L 5 2 = L (156 / 155 : ℝ) - eps L * L (5 : ℝ) := by
  rw [inc_val 5 2 31 sig_5_2]; norm_num
theorem inc_5_2_ge (hL : LogLike L) : (1900231871 / 2500000000000 : ℝ) ≤ inc L 5 2 := by
  have h1 : (8038612912863 / 1250000000000000 : ℝ) ≤ L (156 / 155 : ℝ) := hL.ge_of (q := (156 / 155 : ℝ)) (a := (1 / 155 : ℝ)) (by norm_num) (by norm_num) 3 (by norm_num [lsum, Finset.sum_range_succ])
  have h2 : eps L * L (5 : ℝ) ≤ (88086616111 / 25000000000000 : ℝ) * (321887582487 / 200000000000 : ℝ) :=
    mul_le_mul (eps_le hL) (L5_le hL) (by linarith [L5_ge hL]) (by norm_num)
  rw [inc_5_2_eq]; linarith
theorem inc_5_2_le (hL : LogLike L) : inc L 5 2 ≤ (76009274843 / 100000000000000 : ℝ) := by
  have h1 : L (156 / 155 : ℝ) ≤ (12861780660581 / 2000000000000000 : ℝ) := hL.le_of (q := (156 / 155 : ℝ)) (a := (1 / 155 : ℝ)) (by norm_num) (by norm_num) 3 (by norm_num [lsum, ltail, Finset.sum_range_succ])
  have h2 : (352346464443 / 100000000000000 : ℝ) * (1609437912433 / 1000000000000 : ℝ) ≤ eps L * L (5 : ℝ) :=
    mul_le_mul (eps_ge hL) (L5_ge hL) (by norm_num) (eps_pos hL).le
  rw [inc_5_2_eq]; linarith
theorem sig_5_3 : σ 1 (5 ^ 3) = 156 := by
  rw [sigma_one_apply_prime_pow (by norm_num)]; norm_num [Finset.sum_range_succ]
theorem inc_5_3_eq : inc L 5 3 = L (781 / 780 : ℝ) - eps L * L (5 : ℝ) := by
  rw [inc_val 5 3 156 sig_5_3]; norm_num
theorem inc_5_3_ge (hL : LogLike L) : (-54869592823 / 12500000000000 : ℝ) ≤ inc L 5 3 := by
  have h1 : (2562460312097 / 2000000000000000 : ℝ) ≤ L (781 / 780 : ℝ) := hL.ge_of (q := (781 / 780 : ℝ)) (a := (1 / 780 : ℝ)) (by norm_num) (by norm_num) 2 (by norm_num [lsum, Finset.sum_range_succ])
  have h2 : eps L * L (5 : ℝ) ≤ (88086616111 / 25000000000000 : ℝ) * (321887582487 / 200000000000 : ℝ) :=
    mul_le_mul (eps_le hL) (L5_le hL) (by linarith [L5_ge hL]) (by norm_num)
  rw [inc_5_3_eq]; linarith
theorem inc_5_3_le (hL : LogLike L) : inc L 5 3 ≤ (-438956742581 / 100000000000000 : ℝ) := by
  have h1 : L (781 / 780 : ℝ) ≤ (6406150780243 / 5000000000000000 : ℝ) := hL.le_of (q := (781 / 780 : ℝ)) (a := (1 / 780 : ℝ)) (by norm_num) (by norm_num) 2 (by norm_num [lsum, ltail, Finset.sum_range_succ])
  have h2 : (352346464443 / 100000000000000 : ℝ) * (1609437912433 / 1000000000000 : ℝ) ≤ eps L * L (5 : ℝ) :=
    mul_le_mul (eps_ge hL) (L5_ge hL) (by norm_num) (eps_pos hL).le
  rw [inc_5_3_eq]; linarith
theorem sig_7_1 : σ 1 (7 ^ 1) = 8 := by
  rw [sigma_one_apply_prime_pow (by norm_num)]; norm_num [Finset.sum_range_succ]
theorem inc_7_1_eq : inc L 7 1 = L (57 / 56 : ℝ) - eps L * L (7 : ℝ) := by
  rw [inc_val 7 1 8 sig_7_1]; norm_num
theorem inc_7_1_ge (hL : LogLike L) : (542161574397 / 50000000000000 : ℝ) ≤ inc L 7 1 := by
  have h1 : (176995770994009 / 10000000000000000 : ℝ) ≤ L (57 / 56 : ℝ) := hL.ge_of (q := (57 / 56 : ℝ)) (a := (1 / 56 : ℝ)) (by norm_num) (by norm_num) 4 (by norm_num [lsum, Finset.sum_range_succ])
  have h2 : eps L * L (7 : ℝ) ≤ (88086616111 / 25000000000000 : ℝ) * (30404846079 / 15625000000 : ℝ) :=
    mul_le_mul (eps_le hL) (L7_le hL) (by linarith [L7_ge hL]) (by norm_num)
  rw [inc_7_1_eq]; linarith
theorem inc_7_1_le (hL : LogLike L) : inc L 7 1 ≤ (1084323148797 / 100000000000000 : ℝ) := by
  have h1 : L (57 / 56 : ℝ) ≤ (17699577099401 / 1000000000000000 : ℝ) := hL.le_of (q := (57 / 56 : ℝ)) (a := (1 / 56 : ℝ)) (by norm_num) (by norm_num) 4 (by norm_num [lsum, ltail, Finset.sum_range_succ])
  have h2 : (352346464443 / 100000000000000 : ℝ) * (389182029811 / 200000000000 : ℝ) ≤ eps L * L (7 : ℝ) :=
    mul_le_mul (eps_ge hL) (L7_ge hL) (by norm_num) (eps_pos hL).le
  rw [inc_7_1_eq]; linarith
theorem sig_7_2 : σ 1 (7 ^ 2) = 57 := by
  rw [sigma_one_apply_prime_pow (by norm_num)]; norm_num [Finset.sum_range_succ]
theorem inc_7_2_eq : inc L 7 2 = L (400 / 399 : ℝ) - eps L * L (7 : ℝ) := by
  rw [inc_val 7 2 57 sig_7_2]; norm_num
theorem inc_7_2_ge (hL : LogLike L) : (-217660769667 / 50000000000000 : ℝ) ≤ inc L 7 2 := by
  have h1 : (5006260436237 / 2000000000000000 : ℝ) ≤ L (400 / 399 : ℝ) := hL.ge_of (q := (400 / 399 : ℝ)) (a := (1 / 399 : ℝ)) (by norm_num) (by norm_num) 3 (by norm_num [lsum, Finset.sum_range_succ])
  have h2 : eps L * L (7 : ℝ) ≤ (88086616111 / 25000000000000 : ℝ) * (30404846079 / 15625000000 : ℝ) :=
    mul_le_mul (eps_le hL) (L7_le hL) (by linarith [L7_ge hL]) (by norm_num)
  rw [inc_7_2_eq]; linarith
theorem inc_7_2_le (hL : LogLike L) : inc L 7 2 ≤ (-435321539331 / 100000000000000 : ℝ) := by
  have h1 : L (400 / 399 : ℝ) ≤ (12515651090593 / 5000000000000000 : ℝ) := hL.le_of (q := (400 / 399 : ℝ)) (a := (1 / 399 : ℝ)) (by norm_num) (by norm_num) 3 (by norm_num [lsum, ltail, Finset.sum_range_succ])
  have h2 : (352346464443 / 100000000000000 : ℝ) * (389182029811 / 200000000000 : ℝ) ≤ eps L * L (7 : ℝ) :=
    mul_le_mul (eps_ge hL) (L7_ge hL) (by norm_num) (eps_pos hL).le
  rw [inc_7_2_eq]; linarith

end OddRobin
