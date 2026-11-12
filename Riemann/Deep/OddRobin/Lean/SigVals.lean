import Riemann.Deep.OddRobin.Lean.Defs
import Mathlib.Tactic.NormNum.Prime
import Mathlib.Tactic.NormNum.GCD

namespace OddRobin

open ArithmeticFunction
open scoped ArithmeticFunction.sigma

variable {L : ℝ → ℝ}

lemma sig_mul {a b sa sb : ℕ} (h : Nat.Coprime a b) (ha : σ 1 a = sa) (hb : σ 1 b = sb) :
    σ 1 (a * b) = sa * sb := by
  rw [isMultiplicative_sigma.map_mul_of_coprime h, ha, hb]
lemma sv_3_4 : σ 1 81 = 121 := by
  have h := sigma_one_apply_prime_pow (p := 3) (i := 4) (by norm_num)
  norm_num [Finset.sum_range_succ] at h
  exact h
lemma sv_3_5 : σ 1 243 = 364 := by
  have h := sigma_one_apply_prime_pow (p := 3) (i := 5) (by norm_num)
  norm_num [Finset.sum_range_succ] at h
  exact h
lemma sv_5_2 : σ 1 25 = 31 := by
  have h := sigma_one_apply_prime_pow (p := 5) (i := 2) (by norm_num)
  norm_num [Finset.sum_range_succ] at h
  exact h
lemma sv_5_3 : σ 1 125 = 156 := by
  have h := sigma_one_apply_prime_pow (p := 5) (i := 3) (by norm_num)
  norm_num [Finset.sum_range_succ] at h
  exact h
lemma sv_7_2 : σ 1 49 = 57 := by
  have h := sigma_one_apply_prime_pow (p := 7) (i := 2) (by norm_num)
  norm_num [Finset.sum_range_succ] at h
  exact h
lemma sv_11_1 : σ 1 11 = 12 := by
  have h := sigma_one_apply_prime_pow (p := 11) (i := 1) (by norm_num)
  norm_num [Finset.sum_range_succ] at h
  exact h
lemma sv_11_2 : σ 1 121 = 133 := by
  have h := sigma_one_apply_prime_pow (p := 11) (i := 2) (by norm_num)
  norm_num [Finset.sum_range_succ] at h
  exact h
lemma sv_13_1 : σ 1 13 = 14 := by
  have h := sigma_one_apply_prime_pow (p := 13) (i := 1) (by norm_num)
  norm_num [Finset.sum_range_succ] at h
  exact h
lemma sv_17_1 : σ 1 17 = 18 := by
  have h := sigma_one_apply_prime_pow (p := 17) (i := 1) (by norm_num)
  norm_num [Finset.sum_range_succ] at h
  exact h
lemma sv_19_1 : σ 1 19 = 20 := by
  have h := sigma_one_apply_prime_pow (p := 19) (i := 1) (by norm_num)
  norm_num [Finset.sum_range_succ] at h
  exact h
lemma sv_23_1 : σ 1 23 = 24 := by
  have h := sigma_one_apply_prime_pow (p := 23) (i := 1) (by norm_num)
  norm_num [Finset.sum_range_succ] at h
  exact h
lemma sv_29_1 : σ 1 29 = 30 := by
  have h := sigma_one_apply_prime_pow (p := 29) (i := 1) (by norm_num)
  norm_num [Finset.sum_range_succ] at h
  exact h
lemma sv_31_1 : σ 1 31 = 32 := by
  have h := sigma_one_apply_prime_pow (p := 31) (i := 1) (by norm_num)
  norm_num [Finset.sum_range_succ] at h
  exact h
lemma sv_37_1 : σ 1 37 = 38 := by
  have h := sigma_one_apply_prime_pow (p := 37) (i := 1) (by norm_num)
  norm_num [Finset.sum_range_succ] at h
  exact h
lemma sv_41_1 : σ 1 41 = 42 := by
  have h := sigma_one_apply_prime_pow (p := 41) (i := 1) (by norm_num)
  norm_num [Finset.sum_range_succ] at h
  exact h
lemma sv_43_1 : σ 1 43 = 44 := by
  have h := sigma_one_apply_prime_pow (p := 43) (i := 1) (by norm_num)
  norm_num [Finset.sum_range_succ] at h
  exact h
lemma sv_47_1 : σ 1 47 = 48 := by
  have h := sigma_one_apply_prime_pow (p := 47) (i := 1) (by norm_num)
  norm_num [Finset.sum_range_succ] at h
  exact h
lemma sv_53_1 : σ 1 53 = 54 := by
  have h := sigma_one_apply_prime_pow (p := 53) (i := 1) (by norm_num)
  norm_num [Finset.sum_range_succ] at h
  exact h
lemma sv_59_1 : σ 1 59 = 60 := by
  have h := sigma_one_apply_prime_pow (p := 59) (i := 1) (by norm_num)
  norm_num [Finset.sum_range_succ] at h
  exact h
lemma sv_61_1 : σ 1 61 = 62 := by
  have h := sigma_one_apply_prime_pow (p := 61) (i := 1) (by norm_num)
  norm_num [Finset.sum_range_succ] at h
  exact h
lemma sv_67_1 : σ 1 67 = 68 := by
  have h := sigma_one_apply_prime_pow (p := 67) (i := 1) (by norm_num)
  norm_num [Finset.sum_range_succ] at h
  exact h
lemma chNN_c15 : σ 1 (59 * 61) = 60 * 62 :=
  sig_mul (by norm_num) sv_59_1 sv_61_1
lemma chNN_c14 : σ 1 (53 * (59 * 61)) = 54 * (60 * 62) :=
  sig_mul (by norm_num) sv_53_1 chNN_c15
lemma chNN_c13 : σ 1 (47 * (53 * (59 * 61))) = 48 * (54 * (60 * 62)) :=
  sig_mul (by norm_num) sv_47_1 chNN_c14
lemma chNN_c12 : σ 1 (43 * (47 * (53 * (59 * 61)))) = 44 * (48 * (54 * (60 * 62))) :=
  sig_mul (by norm_num) sv_43_1 chNN_c13
lemma chNN_c11 : σ 1 (41 * (43 * (47 * (53 * (59 * 61))))) = 42 * (44 * (48 * (54 * (60 * 62)))) :=
  sig_mul (by norm_num) sv_41_1 chNN_c12
lemma chNN_c10 : σ 1 (37 * (41 * (43 * (47 * (53 * (59 * 61)))))) = 38 * (42 * (44 * (48 * (54 * (60 * 62))))) :=
  sig_mul (by norm_num) sv_37_1 chNN_c11
lemma chNN_c9 : σ 1 (31 * (37 * (41 * (43 * (47 * (53 * (59 * 61))))))) = 32 * (38 * (42 * (44 * (48 * (54 * (60 * 62)))))) :=
  sig_mul (by norm_num) sv_31_1 chNN_c10
lemma chNN_c8 : σ 1 (29 * (31 * (37 * (41 * (43 * (47 * (53 * (59 * 61)))))))) = 30 * (32 * (38 * (42 * (44 * (48 * (54 * (60 * 62))))))) :=
  sig_mul (by norm_num) sv_29_1 chNN_c9
lemma chNN_c7 : σ 1 (23 * (29 * (31 * (37 * (41 * (43 * (47 * (53 * (59 * 61))))))))) = 24 * (30 * (32 * (38 * (42 * (44 * (48 * (54 * (60 * 62)))))))) :=
  sig_mul (by norm_num) sv_23_1 chNN_c8
lemma chNN_c6 : σ 1 (19 * (23 * (29 * (31 * (37 * (41 * (43 * (47 * (53 * (59 * 61)))))))))) = 20 * (24 * (30 * (32 * (38 * (42 * (44 * (48 * (54 * (60 * 62))))))))) :=
  sig_mul (by norm_num) sv_19_1 chNN_c7
lemma chNN_c5 : σ 1 (17 * (19 * (23 * (29 * (31 * (37 * (41 * (43 * (47 * (53 * (59 * 61))))))))))) = 18 * (20 * (24 * (30 * (32 * (38 * (42 * (44 * (48 * (54 * (60 * 62)))))))))) :=
  sig_mul (by norm_num) sv_17_1 chNN_c6
lemma chNN_c4 : σ 1 (13 * (17 * (19 * (23 * (29 * (31 * (37 * (41 * (43 * (47 * (53 * (59 * 61)))))))))))) = 14 * (18 * (20 * (24 * (30 * (32 * (38 * (42 * (44 * (48 * (54 * (60 * 62))))))))))) :=
  sig_mul (by norm_num) sv_13_1 chNN_c5
lemma chNN_c3 : σ 1 (11 * (13 * (17 * (19 * (23 * (29 * (31 * (37 * (41 * (43 * (47 * (53 * (59 * 61))))))))))))) = 12 * (14 * (18 * (20 * (24 * (30 * (32 * (38 * (42 * (44 * (48 * (54 * (60 * 62)))))))))))) :=
  sig_mul (by norm_num) sv_11_1 chNN_c4
lemma chNN_c2 : σ 1 (49 * (11 * (13 * (17 * (19 * (23 * (29 * (31 * (37 * (41 * (43 * (47 * (53 * (59 * 61)))))))))))))) = 57 * (12 * (14 * (18 * (20 * (24 * (30 * (32 * (38 * (42 * (44 * (48 * (54 * (60 * 62))))))))))))) :=
  sig_mul (by norm_num) sv_7_2 chNN_c3
lemma chNN_c1 : σ 1 (125 * (49 * (11 * (13 * (17 * (19 * (23 * (29 * (31 * (37 * (41 * (43 * (47 * (53 * (59 * 61))))))))))))))) = 156 * (57 * (12 * (14 * (18 * (20 * (24 * (30 * (32 * (38 * (42 * (44 * (48 * (54 * (60 * 62)))))))))))))) :=
  sig_mul (by norm_num) sv_5_3 chNN_c2
lemma chNN_c0 : σ 1 (81 * (125 * (49 * (11 * (13 * (17 * (19 * (23 * (29 * (31 * (37 * (41 * (43 * (47 * (53 * (59 * 61)))))))))))))))) = 121 * (156 * (57 * (12 * (14 * (18 * (20 * (24 * (30 * (32 * (38 * (42 * (44 * (48 * (54 * (60 * 62))))))))))))))) :=
  sig_mul (by norm_num) sv_3_4 chNN_c1
lemma sig_NN : σ 1 NN = 1015178899037727378898944000 := by
  rw [show NN = 81 * (125 * (49 * (11 * (13 * (17 * (19 * (23 * (29 * (31 * (37 * (41 * (43 * (47 * (53 * (59 * 61))))))))))))))) by norm_num [NN], chNN_c0] <;> norm_num
lemma chm3_c15 : σ 1 (59 * 61) = 60 * 62 :=
  sig_mul (by norm_num) sv_59_1 sv_61_1
lemma chm3_c14 : σ 1 (53 * (59 * 61)) = 54 * (60 * 62) :=
  sig_mul (by norm_num) sv_53_1 chm3_c15
lemma chm3_c13 : σ 1 (47 * (53 * (59 * 61))) = 48 * (54 * (60 * 62)) :=
  sig_mul (by norm_num) sv_47_1 chm3_c14
lemma chm3_c12 : σ 1 (43 * (47 * (53 * (59 * 61)))) = 44 * (48 * (54 * (60 * 62))) :=
  sig_mul (by norm_num) sv_43_1 chm3_c13
lemma chm3_c11 : σ 1 (41 * (43 * (47 * (53 * (59 * 61))))) = 42 * (44 * (48 * (54 * (60 * 62)))) :=
  sig_mul (by norm_num) sv_41_1 chm3_c12
lemma chm3_c10 : σ 1 (37 * (41 * (43 * (47 * (53 * (59 * 61)))))) = 38 * (42 * (44 * (48 * (54 * (60 * 62))))) :=
  sig_mul (by norm_num) sv_37_1 chm3_c11
lemma chm3_c9 : σ 1 (31 * (37 * (41 * (43 * (47 * (53 * (59 * 61))))))) = 32 * (38 * (42 * (44 * (48 * (54 * (60 * 62)))))) :=
  sig_mul (by norm_num) sv_31_1 chm3_c10
lemma chm3_c8 : σ 1 (29 * (31 * (37 * (41 * (43 * (47 * (53 * (59 * 61)))))))) = 30 * (32 * (38 * (42 * (44 * (48 * (54 * (60 * 62))))))) :=
  sig_mul (by norm_num) sv_29_1 chm3_c9
lemma chm3_c7 : σ 1 (23 * (29 * (31 * (37 * (41 * (43 * (47 * (53 * (59 * 61))))))))) = 24 * (30 * (32 * (38 * (42 * (44 * (48 * (54 * (60 * 62)))))))) :=
  sig_mul (by norm_num) sv_23_1 chm3_c8
lemma chm3_c6 : σ 1 (19 * (23 * (29 * (31 * (37 * (41 * (43 * (47 * (53 * (59 * 61)))))))))) = 20 * (24 * (30 * (32 * (38 * (42 * (44 * (48 * (54 * (60 * 62))))))))) :=
  sig_mul (by norm_num) sv_19_1 chm3_c7
lemma chm3_c5 : σ 1 (17 * (19 * (23 * (29 * (31 * (37 * (41 * (43 * (47 * (53 * (59 * 61))))))))))) = 18 * (20 * (24 * (30 * (32 * (38 * (42 * (44 * (48 * (54 * (60 * 62)))))))))) :=
  sig_mul (by norm_num) sv_17_1 chm3_c6
lemma chm3_c4 : σ 1 (13 * (17 * (19 * (23 * (29 * (31 * (37 * (41 * (43 * (47 * (53 * (59 * 61)))))))))))) = 14 * (18 * (20 * (24 * (30 * (32 * (38 * (42 * (44 * (48 * (54 * (60 * 62))))))))))) :=
  sig_mul (by norm_num) sv_13_1 chm3_c5
lemma chm3_c3 : σ 1 (11 * (13 * (17 * (19 * (23 * (29 * (31 * (37 * (41 * (43 * (47 * (53 * (59 * 61))))))))))))) = 12 * (14 * (18 * (20 * (24 * (30 * (32 * (38 * (42 * (44 * (48 * (54 * (60 * 62)))))))))))) :=
  sig_mul (by norm_num) sv_11_1 chm3_c4
lemma chm3_c2 : σ 1 (49 * (11 * (13 * (17 * (19 * (23 * (29 * (31 * (37 * (41 * (43 * (47 * (53 * (59 * 61)))))))))))))) = 57 * (12 * (14 * (18 * (20 * (24 * (30 * (32 * (38 * (42 * (44 * (48 * (54 * (60 * 62))))))))))))) :=
  sig_mul (by norm_num) sv_7_2 chm3_c3
lemma chm3_c1 : σ 1 (125 * (49 * (11 * (13 * (17 * (19 * (23 * (29 * (31 * (37 * (41 * (43 * (47 * (53 * (59 * 61))))))))))))))) = 156 * (57 * (12 * (14 * (18 * (20 * (24 * (30 * (32 * (38 * (42 * (44 * (48 * (54 * (60 * 62)))))))))))))) :=
  sig_mul (by norm_num) sv_5_3 chm3_c2
lemma chm3_c0 : σ 1 (243 * (125 * (49 * (11 * (13 * (17 * (19 * (23 * (29 * (31 * (37 * (41 * (43 * (47 * (53 * (59 * 61)))))))))))))))) = 364 * (156 * (57 * (12 * (14 * (18 * (20 * (24 * (30 * (32 * (38 * (42 * (44 * (48 * (54 * (60 * 62))))))))))))))) :=
  sig_mul (by norm_num) sv_3_5 chm3_c1
/-- 3N. -/
def m3 : ℕ := 831281402884796906843926125
lemma sig_m3 : σ 1 m3 = 3053926605369692280324096000 := by
  rw [show m3 = 243 * (125 * (49 * (11 * (13 * (17 * (19 * (23 * (29 * (31 * (37 * (41 * (43 * (47 * (53 * (59 * 61))))))))))))))) by norm_num [m3], chm3_c0] <;> norm_num
lemma chm11_c15 : σ 1 (59 * 61) = 60 * 62 :=
  sig_mul (by norm_num) sv_59_1 sv_61_1
lemma chm11_c14 : σ 1 (53 * (59 * 61)) = 54 * (60 * 62) :=
  sig_mul (by norm_num) sv_53_1 chm11_c15
lemma chm11_c13 : σ 1 (47 * (53 * (59 * 61))) = 48 * (54 * (60 * 62)) :=
  sig_mul (by norm_num) sv_47_1 chm11_c14
lemma chm11_c12 : σ 1 (43 * (47 * (53 * (59 * 61)))) = 44 * (48 * (54 * (60 * 62))) :=
  sig_mul (by norm_num) sv_43_1 chm11_c13
lemma chm11_c11 : σ 1 (41 * (43 * (47 * (53 * (59 * 61))))) = 42 * (44 * (48 * (54 * (60 * 62)))) :=
  sig_mul (by norm_num) sv_41_1 chm11_c12
lemma chm11_c10 : σ 1 (37 * (41 * (43 * (47 * (53 * (59 * 61)))))) = 38 * (42 * (44 * (48 * (54 * (60 * 62))))) :=
  sig_mul (by norm_num) sv_37_1 chm11_c11
lemma chm11_c9 : σ 1 (31 * (37 * (41 * (43 * (47 * (53 * (59 * 61))))))) = 32 * (38 * (42 * (44 * (48 * (54 * (60 * 62)))))) :=
  sig_mul (by norm_num) sv_31_1 chm11_c10
lemma chm11_c8 : σ 1 (29 * (31 * (37 * (41 * (43 * (47 * (53 * (59 * 61)))))))) = 30 * (32 * (38 * (42 * (44 * (48 * (54 * (60 * 62))))))) :=
  sig_mul (by norm_num) sv_29_1 chm11_c9
lemma chm11_c7 : σ 1 (23 * (29 * (31 * (37 * (41 * (43 * (47 * (53 * (59 * 61))))))))) = 24 * (30 * (32 * (38 * (42 * (44 * (48 * (54 * (60 * 62)))))))) :=
  sig_mul (by norm_num) sv_23_1 chm11_c8
lemma chm11_c6 : σ 1 (19 * (23 * (29 * (31 * (37 * (41 * (43 * (47 * (53 * (59 * 61)))))))))) = 20 * (24 * (30 * (32 * (38 * (42 * (44 * (48 * (54 * (60 * 62))))))))) :=
  sig_mul (by norm_num) sv_19_1 chm11_c7
lemma chm11_c5 : σ 1 (17 * (19 * (23 * (29 * (31 * (37 * (41 * (43 * (47 * (53 * (59 * 61))))))))))) = 18 * (20 * (24 * (30 * (32 * (38 * (42 * (44 * (48 * (54 * (60 * 62)))))))))) :=
  sig_mul (by norm_num) sv_17_1 chm11_c6
lemma chm11_c4 : σ 1 (13 * (17 * (19 * (23 * (29 * (31 * (37 * (41 * (43 * (47 * (53 * (59 * 61)))))))))))) = 14 * (18 * (20 * (24 * (30 * (32 * (38 * (42 * (44 * (48 * (54 * (60 * 62))))))))))) :=
  sig_mul (by norm_num) sv_13_1 chm11_c5
lemma chm11_c3 : σ 1 (121 * (13 * (17 * (19 * (23 * (29 * (31 * (37 * (41 * (43 * (47 * (53 * (59 * 61))))))))))))) = 133 * (14 * (18 * (20 * (24 * (30 * (32 * (38 * (42 * (44 * (48 * (54 * (60 * 62)))))))))))) :=
  sig_mul (by norm_num) sv_11_2 chm11_c4
lemma chm11_c2 : σ 1 (49 * (121 * (13 * (17 * (19 * (23 * (29 * (31 * (37 * (41 * (43 * (47 * (53 * (59 * 61)))))))))))))) = 57 * (133 * (14 * (18 * (20 * (24 * (30 * (32 * (38 * (42 * (44 * (48 * (54 * (60 * 62))))))))))))) :=
  sig_mul (by norm_num) sv_7_2 chm11_c3
lemma chm11_c1 : σ 1 (125 * (49 * (121 * (13 * (17 * (19 * (23 * (29 * (31 * (37 * (41 * (43 * (47 * (53 * (59 * 61))))))))))))))) = 156 * (57 * (133 * (14 * (18 * (20 * (24 * (30 * (32 * (38 * (42 * (44 * (48 * (54 * (60 * 62)))))))))))))) :=
  sig_mul (by norm_num) sv_5_3 chm11_c2
lemma chm11_c0 : σ 1 (81 * (125 * (49 * (121 * (13 * (17 * (19 * (23 * (29 * (31 * (37 * (41 * (43 * (47 * (53 * (59 * 61)))))))))))))))) = 121 * (156 * (57 * (133 * (14 * (18 * (20 * (24 * (30 * (32 * (38 * (42 * (44 * (48 * (54 * (60 * 62))))))))))))))) :=
  sig_mul (by norm_num) sv_3_4 chm11_c1
/-- 11N. -/
def m11 : ℕ := 3048031810577588658427729125
lemma sig_m11 : σ 1 m11 = 11251566131001478449463296000 := by
  rw [show m11 = 81 * (125 * (49 * (121 * (13 * (17 * (19 * (23 * (29 * (31 * (37 * (41 * (43 * (47 * (53 * (59 * 61))))))))))))))) by norm_num [m11], chm11_c0] <;> norm_num
lemma chm67_c16 : σ 1 (61 * 67) = 62 * 68 :=
  sig_mul (by norm_num) sv_61_1 sv_67_1
lemma chm67_c15 : σ 1 (59 * (61 * 67)) = 60 * (62 * 68) :=
  sig_mul (by norm_num) sv_59_1 chm67_c16
lemma chm67_c14 : σ 1 (53 * (59 * (61 * 67))) = 54 * (60 * (62 * 68)) :=
  sig_mul (by norm_num) sv_53_1 chm67_c15
lemma chm67_c13 : σ 1 (47 * (53 * (59 * (61 * 67)))) = 48 * (54 * (60 * (62 * 68))) :=
  sig_mul (by norm_num) sv_47_1 chm67_c14
lemma chm67_c12 : σ 1 (43 * (47 * (53 * (59 * (61 * 67))))) = 44 * (48 * (54 * (60 * (62 * 68)))) :=
  sig_mul (by norm_num) sv_43_1 chm67_c13
lemma chm67_c11 : σ 1 (41 * (43 * (47 * (53 * (59 * (61 * 67)))))) = 42 * (44 * (48 * (54 * (60 * (62 * 68))))) :=
  sig_mul (by norm_num) sv_41_1 chm67_c12
lemma chm67_c10 : σ 1 (37 * (41 * (43 * (47 * (53 * (59 * (61 * 67))))))) = 38 * (42 * (44 * (48 * (54 * (60 * (62 * 68)))))) :=
  sig_mul (by norm_num) sv_37_1 chm67_c11
lemma chm67_c9 : σ 1 (31 * (37 * (41 * (43 * (47 * (53 * (59 * (61 * 67)))))))) = 32 * (38 * (42 * (44 * (48 * (54 * (60 * (62 * 68))))))) :=
  sig_mul (by norm_num) sv_31_1 chm67_c10
lemma chm67_c8 : σ 1 (29 * (31 * (37 * (41 * (43 * (47 * (53 * (59 * (61 * 67))))))))) = 30 * (32 * (38 * (42 * (44 * (48 * (54 * (60 * (62 * 68)))))))) :=
  sig_mul (by norm_num) sv_29_1 chm67_c9
lemma chm67_c7 : σ 1 (23 * (29 * (31 * (37 * (41 * (43 * (47 * (53 * (59 * (61 * 67)))))))))) = 24 * (30 * (32 * (38 * (42 * (44 * (48 * (54 * (60 * (62 * 68))))))))) :=
  sig_mul (by norm_num) sv_23_1 chm67_c8
lemma chm67_c6 : σ 1 (19 * (23 * (29 * (31 * (37 * (41 * (43 * (47 * (53 * (59 * (61 * 67))))))))))) = 20 * (24 * (30 * (32 * (38 * (42 * (44 * (48 * (54 * (60 * (62 * 68)))))))))) :=
  sig_mul (by norm_num) sv_19_1 chm67_c7
lemma chm67_c5 : σ 1 (17 * (19 * (23 * (29 * (31 * (37 * (41 * (43 * (47 * (53 * (59 * (61 * 67)))))))))))) = 18 * (20 * (24 * (30 * (32 * (38 * (42 * (44 * (48 * (54 * (60 * (62 * 68))))))))))) :=
  sig_mul (by norm_num) sv_17_1 chm67_c6
lemma chm67_c4 : σ 1 (13 * (17 * (19 * (23 * (29 * (31 * (37 * (41 * (43 * (47 * (53 * (59 * (61 * 67))))))))))))) = 14 * (18 * (20 * (24 * (30 * (32 * (38 * (42 * (44 * (48 * (54 * (60 * (62 * 68)))))))))))) :=
  sig_mul (by norm_num) sv_13_1 chm67_c5
lemma chm67_c3 : σ 1 (11 * (13 * (17 * (19 * (23 * (29 * (31 * (37 * (41 * (43 * (47 * (53 * (59 * (61 * 67)))))))))))))) = 12 * (14 * (18 * (20 * (24 * (30 * (32 * (38 * (42 * (44 * (48 * (54 * (60 * (62 * 68))))))))))))) :=
  sig_mul (by norm_num) sv_11_1 chm67_c4
lemma chm67_c2 : σ 1 (49 * (11 * (13 * (17 * (19 * (23 * (29 * (31 * (37 * (41 * (43 * (47 * (53 * (59 * (61 * 67))))))))))))))) = 57 * (12 * (14 * (18 * (20 * (24 * (30 * (32 * (38 * (42 * (44 * (48 * (54 * (60 * (62 * 68)))))))))))))) :=
  sig_mul (by norm_num) sv_7_2 chm67_c3
lemma chm67_c1 : σ 1 (25 * (49 * (11 * (13 * (17 * (19 * (23 * (29 * (31 * (37 * (41 * (43 * (47 * (53 * (59 * (61 * 67)))))))))))))))) = 31 * (57 * (12 * (14 * (18 * (20 * (24 * (30 * (32 * (38 * (42 * (44 * (48 * (54 * (60 * (62 * 68))))))))))))))) :=
  sig_mul (by norm_num) sv_5_2 chm67_c2
lemma chm67_c0 : σ 1 (81 * (25 * (49 * (11 * (13 * (17 * (19 * (23 * (29 * (31 * (37 * (41 * (43 * (47 * (53 * (59 * (61 * 67))))))))))))))))) = 121 * (31 * (57 * (12 * (14 * (18 * (20 * (24 * (30 * (32 * (38 * (42 * (44 * (48 * (54 * (60 * (62 * 68)))))))))))))))) :=
  sig_mul (by norm_num) sv_3_4 chm67_c1
/-- 67N/5. -/
def m67 : ℕ := 3713056932885426183902870025
lemma sig_m67 : σ 1 m67 = 13717930251099546889224192000 := by
  rw [show m67 = 81 * (25 * (49 * (11 * (13 * (17 * (19 * (23 * (29 * (31 * (37 * (41 * (43 * (47 * (53 * (59 * (61 * 67)))))))))))))))) by norm_num [m67], chm67_c0] <;> norm_num

end OddRobin
