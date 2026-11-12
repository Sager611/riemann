import Riemann.Deep.OddRobin.Lean.NumA

namespace OddRobin

open ArithmeticFunction
open scoped ArithmeticFunction.sigma

variable {L : ℝ → ℝ}

theorem sig_11_0 : σ 1 (11 ^ 0) = 1 := by
  rw [sigma_one_apply_prime_pow (by norm_num)]; norm_num [Finset.sum_range_succ]
theorem inc_11_0_eq : inc L 11 0 = L (12 / 11 : ℝ) - eps L * L (11 : ℝ) := by
  rw [inc_val 11 0 1 sig_11_0]; norm_num
theorem inc_11_0_ge (hL : LogLike L) : (1571249555497 / 20000000000000 : ℝ) ≤ inc L 11 0 := by
  have h1 : (174022753979259 / 2000000000000000 : ℝ) ≤ L (12 / 11 : ℝ) := hL.ge_of (q := (12 / 11 : ℝ)) (a := (1 / 11 : ℝ)) (by norm_num) (by norm_num) 5 (by norm_num [lsum, Finset.sum_range_succ])
  have h2 : eps L * L (11 : ℝ) ≤ (88086616111 / 25000000000000 : ℝ) * (2397895272799 / 1000000000000 : ℝ) :=
    mul_le_mul (eps_le hL) (L11_le hL) (by linarith [L11_ge hL]) (by norm_num)
  rw [inc_11_0_eq]; linarith
theorem inc_11_0_le (hL : LogLike L) : inc L 11 0 ≤ (7856247777489 / 100000000000000 : ℝ) := by
  have h1 : L (12 / 11 : ℝ) ≤ (435056884948149 / 5000000000000000 : ℝ) := hL.le_of (q := (12 / 11 : ℝ)) (a := (1 / 11 : ℝ)) (by norm_num) (by norm_num) 5 (by norm_num [lsum, ltail, Finset.sum_range_succ])
  have h2 : (352346464443 / 100000000000000 : ℝ) * (1198947636399 / 500000000000 : ℝ) ≤ eps L * L (11 : ℝ) :=
    mul_le_mul (eps_ge hL) (L11_ge hL) (by norm_num) (eps_pos hL).le
  rw [inc_11_0_eq]; linarith
theorem sig_11_1 : σ 1 (11 ^ 1) = 12 := by
  rw [sigma_one_apply_prime_pow (by norm_num)]; norm_num [Finset.sum_range_succ]
theorem inc_11_1_eq : inc L 11 1 = L (133 / 132 : ℝ) - eps L * L (11 : ℝ) := by
  rw [inc_val 11 1 12 sig_11_1]; norm_num
theorem inc_11_1_ge (hL : LogLike L) : (-4508467897 / 5000000000000 : ℝ) ≤ inc L 11 1 := by
  have h1 : (75472056353829 / 10000000000000000 : ℝ) ≤ L (133 / 132 : ℝ) := hL.ge_of (q := (133 / 132 : ℝ)) (a := (1 / 132 : ℝ)) (by norm_num) (by norm_num) 3 (by norm_num [lsum, Finset.sum_range_succ])
  have h2 : eps L * L (11 : ℝ) ≤ (88086616111 / 25000000000000 : ℝ) * (2397895272799 / 1000000000000 : ℝ) :=
    mul_le_mul (eps_le hL) (L11_le hL) (by linarith [L11_ge hL]) (by norm_num)
  rw [inc_11_1_eq]; linarith
theorem inc_11_1_le (hL : LogLike L) : inc L 11 1 ≤ (-5635584871 / 6250000000000 : ℝ) := by
  have h1 : L (133 / 132 : ℝ) ≤ (7547205635383 / 1000000000000000 : ℝ) := hL.le_of (q := (133 / 132 : ℝ)) (a := (1 / 132 : ℝ)) (by norm_num) (by norm_num) 3 (by norm_num [lsum, ltail, Finset.sum_range_succ])
  have h2 : (352346464443 / 100000000000000 : ℝ) * (1198947636399 / 500000000000 : ℝ) ≤ eps L * L (11 : ℝ) :=
    mul_le_mul (eps_ge hL) (L11_ge hL) (by norm_num) (eps_pos hL).le
  rw [inc_11_1_eq]; linarith
theorem sig_11_2 : σ 1 (11 ^ 2) = 133 := by
  rw [sigma_one_apply_prime_pow (by norm_num)]; norm_num [Finset.sum_range_succ]
theorem inc_11_2_eq : inc L 11 2 = L (1464 / 1463 : ℝ) - eps L * L (11 : ℝ) := by
  rw [inc_val 11 2 133 sig_11_2]; norm_num
theorem inc_11_2_ge (hL : LogLike L) : (-155312114273 / 20000000000000 : ℝ) ≤ inc L 11 2 := by
  have h1 : (273317400453 / 400000000000000 : ℝ) ≤ L (1464 / 1463 : ℝ) := hL.ge_of (q := (1464 / 1463 : ℝ)) (a := (1 / 1463 : ℝ)) (by norm_num) (by norm_num) 2 (by norm_num [lsum, Finset.sum_range_succ])
  have h2 : eps L * L (11 : ℝ) ≤ (88086616111 / 25000000000000 : ℝ) * (2397895272799 / 1000000000000 : ℝ) :=
    mul_le_mul (eps_le hL) (L11_le hL) (by linarith [L11_ge hL]) (by norm_num)
  rw [inc_11_2_eq]; linarith
theorem inc_11_2_le (hL : LogLike L) : inc L 11 2 ≤ (-776560571361 / 100000000000000 : ℝ) := by
  have h1 : L (1464 / 1463 : ℝ) ≤ (3416467505663 / 5000000000000000 : ℝ) := hL.le_of (q := (1464 / 1463 : ℝ)) (a := (1 / 1463 : ℝ)) (by norm_num) (by norm_num) 2 (by norm_num [lsum, ltail, Finset.sum_range_succ])
  have h2 : (352346464443 / 100000000000000 : ℝ) * (1198947636399 / 500000000000 : ℝ) ≤ eps L * L (11 : ℝ) :=
    mul_le_mul (eps_ge hL) (L11_ge hL) (by norm_num) (eps_pos hL).le
  rw [inc_11_2_eq]; linarith
theorem sig_67_1 : σ 1 (67 ^ 1) = 68 := by
  rw [sigma_one_apply_prime_pow (by norm_num)]; norm_num [Finset.sum_range_succ]
theorem inc_67_1_eq : inc L 67 1 = L (4557 / 4556 : ℝ) - eps L * L (67 : ℝ) := by
  rw [inc_val 67 1 68 sig_67_1]; norm_num
theorem inc_67_1_ge (hL : LogLike L) : (-364890477209 / 25000000000000 : ℝ) ≤ inc L 67 1 := by
  have h1 : (2194666968097 / 10000000000000000 : ℝ) ≤ L (4557 / 4556 : ℝ) := hL.ge_of (q := (4557 / 4556 : ℝ)) (a := (1 / 4556 : ℝ)) (by norm_num) (by norm_num) 2 (by norm_num [lsum, Finset.sum_range_succ])
  have h2 : eps L * L (67 : ℝ) ≤ (88086616111 / 25000000000000 : ℝ) * (32849161089 / 7812500000 : ℝ) :=
    mul_le_mul (eps_le hL) (L67_le hL) (by linarith [L67_ge hL]) (by norm_num)
  rw [inc_67_1_eq]; linarith
theorem inc_67_1_le (hL : LogLike L) : inc L 67 1 ≤ (-145956190883 / 10000000000000 : ℝ) := by
  have h1 : L (4557 / 4556 : ℝ) ≤ (1097333484049 / 5000000000000000 : ℝ) := hL.le_of (q := (4557 / 4556 : ℝ)) (a := (1 / 4556 : ℝ)) (by norm_num) (by norm_num) 2 (by norm_num [lsum, ltail, Finset.sum_range_succ])
  have h2 : (352346464443 / 100000000000000 : ℝ) * (420469261939 / 100000000000 : ℝ) ≤ eps L * L (67 : ℝ) :=
    mul_le_mul (eps_ge hL) (L67_ge hL) (by norm_num) (eps_pos hL).le
  rw [inc_67_1_eq]; linarith
theorem sig_71_0 : σ 1 (71 ^ 0) = 1 := by
  rw [sigma_one_apply_prime_pow (by norm_num)]; norm_num [Finset.sum_range_succ]
theorem inc_71_0_eq : inc L 71 0 = L (72 / 71 : ℝ) - eps L * L (71 : ℝ) := by
  rw [inc_val 71 0 1 sig_71_0]; norm_num
theorem inc_71_0_ge (hL : LogLike L) : (-103315986259 / 100000000000000 : ℝ) ≤ inc L 71 0 := by
  have h1 : (34965604936849 / 2500000000000000 : ℝ) ≤ L (72 / 71 : ℝ) := hL.ge_of (q := (72 / 71 : ℝ)) (a := (1 / 71 : ℝ)) (by norm_num) (by norm_num) 3 (by norm_num [lsum, Finset.sum_range_succ])
  have h2 : eps L * L (71 : ℝ) ≤ (88086616111 / 25000000000000 : ℝ) * (2131339938521 / 500000000000 : ℝ) :=
    mul_le_mul (eps_le hL) (L71_le hL) (by linarith [L71_ge hL]) (by norm_num)
  rw [inc_71_0_eq]; linarith
theorem inc_71_0_le (hL : LogLike L) : inc L 71 0 ≤ (-103315986253 / 100000000000000 : ℝ) := by
  have h1 : L (72 / 71 : ℝ) ≤ (139862419747399 / 10000000000000000 : ℝ) := hL.le_of (q := (72 / 71 : ℝ)) (a := (1 / 71 : ℝ)) (by norm_num) (by norm_num) 3 (by norm_num [lsum, ltail, Finset.sum_range_succ])
  have h2 : (352346464443 / 100000000000000 : ℝ) * (4262679877041 / 1000000000000 : ℝ) ≤ eps L * L (71 : ℝ) :=
    mul_le_mul (eps_ge hL) (L71_ge hL) (by norm_num) (eps_pos hL).le
  rw [inc_71_0_eq]; linarith
theorem sig_71_1 : σ 1 (71 ^ 1) = 72 := by
  rw [sigma_one_apply_prime_pow (by norm_num)]; norm_num [Finset.sum_range_succ]
theorem inc_71_1_eq : inc L 71 1 = L (5113 / 5112 : ℝ) - eps L * L (71 : ℝ) := by
  rw [inc_val 71 1 72 sig_71_1]; norm_num
theorem inc_71_1_ge (hL : LogLike L) : (-148238028147 / 10000000000000 : ℝ) ≤ inc L 71 1 := by
  have h1 : (391198045257 / 2000000000000000 : ℝ) ≤ L (5113 / 5112 : ℝ) := hL.ge_of (q := (5113 / 5112 : ℝ)) (a := (1 / 5112 : ℝ)) (by norm_num) (by norm_num) 2 (by norm_num [lsum, Finset.sum_range_succ])
  have h2 : eps L * L (71 : ℝ) ≤ (88086616111 / 25000000000000 : ℝ) * (2131339938521 / 500000000000 : ℝ) :=
    mul_le_mul (eps_le hL) (L71_le hL) (by linarith [L71_ge hL]) (by norm_num)
  rw [inc_71_1_eq]; linarith
theorem inc_71_1_le (hL : LogLike L) : inc L 71 1 ≤ (-185297535183 / 12500000000000 : ℝ) := by
  have h1 : L (5113 / 5112 : ℝ) ≤ (977995113143 / 5000000000000000 : ℝ) := hL.le_of (q := (5113 / 5112 : ℝ)) (a := (1 / 5112 : ℝ)) (by norm_num) (by norm_num) 2 (by norm_num [lsum, ltail, Finset.sum_range_succ])
  have h2 : (352346464443 / 100000000000000 : ℝ) * (4262679877041 / 1000000000000 : ℝ) ≤ eps L * L (71 : ℝ) :=
    mul_le_mul (eps_ge hL) (L71_ge hL) (by norm_num) (eps_pos hL).le
  rw [inc_71_1_eq]; linarith
theorem sig_67_0 : σ 1 (67 ^ 0) = 1 := by
  rw [sigma_one_apply_prime_pow (by norm_num)]; norm_num
theorem inc_67_0 (hL : LogLike L) : inc L 67 0 = 0 := by
  rw [inc_val 67 0 1 sig_67_0]
  have : (1 : ℝ) + 1 / (((67 : ℕ) : ℝ) * ((1 : ℕ) : ℝ)) = 68 / 67 := by norm_num
  rw [this]; unfold eps; push_cast; rw [div_mul_cancel₀ _ (L67_pos hL).ne']; ring
theorem L_u61_ge (hL : LogLike L) : (40651302179449 / 2500000000000000 : ℝ) ≤ L (62 / 61 : ℝ) := hL.ge_of (q := (62 / 61 : ℝ)) (a := (1 / 61 : ℝ)) (by norm_num) (by norm_num) 3 (by norm_num [lsum, Finset.sum_range_succ])
theorem L_u61_le (hL : LogLike L) : L (62 / 61 : ℝ) ≤ (40651302179451 / 2500000000000000 : ℝ) := hL.le_of (q := (62 / 61 : ℝ)) (a := (1 / 61 : ℝ)) (by norm_num) (by norm_num) 3 (by norm_num [lsum, ltail, Finset.sum_range_succ])
theorem L_u182_ge (hL : LogLike L) : (10958931529251 / 2000000000000000 : ℝ) ≤ L (183 / 182 : ℝ) := hL.ge_of (q := (183 / 182 : ℝ)) (a := (1 / 182 : ℝ)) (by norm_num) (by norm_num) 3 (by norm_num [lsum, Finset.sum_range_succ])
theorem L_u182_le (hL : LogLike L) : L (183 / 182 : ℝ) ≤ (3424666102891 / 625000000000000 : ℝ) := hL.le_of (q := (183 / 182 : ℝ)) (a := (1 / 182 : ℝ)) (by norm_num) (by norm_num) 3 (by norm_num [lsum, ltail, Finset.sum_range_succ])
theorem L_u73_ge (hL : LogLike L) : (17007065069723 / 1250000000000000 : ℝ) ≤ L (74 / 73 : ℝ) := hL.ge_of (q := (74 / 73 : ℝ)) (a := (1 / 73 : ℝ)) (by norm_num) (by norm_num) 3 (by norm_num [lsum, Finset.sum_range_succ])
theorem L_u73_le (hL : LogLike L) : L (74 / 73 : ℝ) ≤ (136056520557787 / 10000000000000000 : ℝ) := hL.le_of (q := (74 / 73 : ℝ)) (a := (1 / 73 : ℝ)) (by norm_num) (by norm_num) 3 (by norm_num [lsum, ltail, Finset.sum_range_succ])

end OddRobin
