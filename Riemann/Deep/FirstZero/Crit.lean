import Riemann.Deep.FirstZero.EM
import Riemann.Deep.FirstZero.Ball
import Mathlib.Analysis.Real.Pi.Bounds
import Mathlib.Analysis.Complex.ExponentialBounds

/-!
# ζ(1/2 + it) through the order-8 Euler–Maclaurin formula with N = 6

ζ(s) ≈ 1 + 2^{-s} + 3^{-s} + 4^{-s} + 5^{-s} + 6^{-s}·K(s), with 4^{-s} = (2^{-s})², 6^{-s} = 2^{-s}3^{-s},
so only p^{-s} for p = 2, 3, 5 need enclosures: p^{-s} = p^{-1/2} e^{-i t log p}.
-/

open Complex

namespace Fayed.FirstZero

/-- The coefficient of 6^{-s} in the N = 6 Euler–Maclaurin main term. -/
noncomputable def EMK (s : ℂ) : ℂ := 1 - 6 / (1 - s) - 1 / 2 +
  s * (1 / 12 * gc s 0 / 6 - 1 / 720 * gc s 2 / 216 + 1 / 30240 * gc s 4 / 7776
    - 1 / 1209600 * gc s 6 / 279936)

lemma cpow_shift6 (s : ℂ) (j : ℕ) :
    ((((6 : ℕ) : ℝ) : ℂ)) ^ (-s - 1 - (j : ℂ)) = (6 : ℂ) ^ (-s) * ((6 : ℂ) ^ (j + 1))⁻¹ := by
  have h6 : ((((6 : ℕ) : ℝ) : ℂ)) = (6 : ℂ) := by norm_num
  rw [h6, show -s - 1 - (j : ℂ) = -s + (-((j + 1 : ℕ) : ℂ)) by push_cast; ring,
    cpow_add _ _ (by norm_num), cpow_neg (6 : ℂ) ((j + 1 : ℕ) : ℂ), cpow_natCast]

lemma em_main_eq (s : ℂ) (hs : s ≠ 0) (hs1 : s ≠ 1) :
    (∑ n ∈ Finset.range (6 + 1), 1 / (n : ℂ) ^ s) + (-((6 : ℕ) : ℂ) ^ (1 - s)) / (1 - s)
      + (-((6 : ℕ) : ℂ) ^ (-s)) / 2 + s * EMPhi s ((6 : ℕ) : ℝ) =
    1 + (2 : ℂ) ^ (-s) + (3 : ℂ) ^ (-s) + (2 : ℂ) ^ (-s) * (2 : ℂ) ^ (-s) + (5 : ℂ) ^ (-s)
      + ((2 : ℂ) ^ (-s) * (3 : ℂ) ^ (-s)) * EMK s := by
  have h4 : (4 : ℂ) ^ (-s) = (2 : ℂ) ^ (-s) * (2 : ℂ) ^ (-s) := by
    have := natCast_mul_natCast_cpow 2 2 (-s); norm_num at this; exact this
  have h6 : (6 : ℂ) ^ (-s) = (2 : ℂ) ^ (-s) * (3 : ℂ) ^ (-s) := by
    have := natCast_mul_natCast_cpow 2 3 (-s); norm_num at this; exact this
  have h6' : (6 : ℂ) ^ (1 - s) = 6 * (6 : ℂ) ^ (-s) := by
    rw [sub_eq_add_neg, cpow_add _ _ (by norm_num), cpow_one]
  have hs1' : (1 : ℂ) - s ≠ 0 := sub_ne_zero.mpr (Ne.symm hs1)
  simp only [EMPhi, gf, cpow_shift6]
  simp only [Finset.sum_range_succ, Finset.sum_range_zero, Nat.cast_zero, Nat.cast_one,
    Nat.cast_ofNat, zero_cpow (neg_ne_zero.mpr hs), one_cpow, zero_add, one_div, ← cpow_neg]
  rw [h6', h4, h6]
  simp only [EMK]
  field_simp
  ring

/-- p^{-(1/2 + i t)} = p^{-1/2} e^{-i t log p}. -/
lemma cpow_crit {x : ℝ} (hx : 0 < x) (t : ℝ) :
    (x : ℂ) ^ (-((1 / 2 : ℂ) + (t : ℂ) * I)) =
      (((Real.sqrt x)⁻¹ : ℝ) : ℂ) * cexp (((-(t * Real.log x)) : ℝ) * I) := by
  have hx0 : (x : ℂ) ≠ 0 := by exact_mod_cast hx.ne'
  rw [cpow_def_of_ne_zero hx0, ← ofReal_log hx.le]
  have hs : (Real.sqrt x)⁻¹ = Real.exp (-(Real.log x / 2)) := by
    rw [Real.sqrt_eq_rpow, Real.rpow_def_of_pos hx, ← Real.exp_neg]
    ring_nf
  rw [hs, ofReal_exp, ← Complex.exp_add]
  congr 1
  push_cast
  ring

lemma log2_near : |Real.log 2 - 69314718055 / 10 ^ 11| ≤ 3 / 10 ^ 10 := by
  have := Real.log_two_gt_d9; have := Real.log_two_lt_d9
  rw [abs_le]; constructor <;> norm_num at * <;> linarith

lemma log3_near : |Real.log 3 - 10986122887 / 10 ^ 10| ≤ 1 / 10 ^ 9 := by
  have h := Real.abs_log_sub_add_sum_range_le (x := (1 / 4 : ℝ)) (by norm_num [abs_of_pos]) 20
  have e : Real.log 3 = 2 * Real.log 2 + Real.log (3 / 4) := by
    rw [show (3 : ℝ) = 2 ^ 2 * (3 / 4) by norm_num, Real.log_mul (by norm_num) (by norm_num),
      Real.log_pow]
    push_cast; ring
  have h2 := log2_near
  rw [abs_le] at h h2 ⊢
  norm_num [Finset.sum_range_succ] at h h2 ⊢
  rw [e]
  constructor <;> linarith [h.1, h.2, h2.1, h2.2]

lemma log5_near : |Real.log 5 - 16094379124 / 10 ^ 10| ≤ 1 / 10 ^ 9 := by
  have h := Real.abs_log_sub_add_sum_range_le (x := (1 / 5 : ℝ)) (by norm_num [abs_of_pos]) 20
  have e : Real.log 5 = 2 * Real.log 2 - Real.log (4 / 5) := by
    rw [show (5 : ℝ) = 2 ^ 2 / (4 / 5) by norm_num, Real.log_div (by norm_num) (by norm_num),
      Real.log_pow]
    push_cast; ring
  have h2 := log2_near
  rw [abs_le] at h h2 ⊢
  norm_num [Finset.sum_range_succ] at h h2 ⊢
  rw [e]
  constructor <;> linarith [h.1, h.2, h2.1, h2.2]

lemma norm_gc_le (s : ℂ) (c : ℕ → ℝ) (n : ℕ) (hc : ∀ j < n, ‖s + 1 + j‖ ≤ c j) :
    ‖gc s n‖ ≤ ∏ j ∈ Finset.range n, c j := by
  induction n with
  | zero => simp [gc]
  | succ n ih =>
    rw [Finset.prod_range_succ]
    have h1 : ‖gc s (n + 1)‖ = ‖gc s n‖ * ‖s + 1 + n‖ := by
      rw [gc, norm_mul, show -s - 1 - (n : ℂ) = -(s + 1 + n) by ring, norm_neg]
    rw [h1]
    exact mul_le_mul (ih fun j hj => hc j (by omega)) (hc n (by omega)) (norm_nonneg _)
      (le_trans (norm_nonneg _) (ih fun j hj => hc j (by omega)))

lemma six_rpow_le : (6 : ℝ) ^ (-((1 / 2 : ℝ) + 7)) ≤ 1 / (279936 * (2449 / 1000)) := by
  have h6 : (0 : ℝ) < 6 := by norm_num
  rw [Real.rpow_neg h6.le, Real.rpow_add h6, show (7 : ℝ) = ((7 : ℕ) : ℝ) by norm_num,
    Real.rpow_natCast, ← Real.sqrt_eq_rpow, one_div]
  have hs : (2449 / 1000 : ℝ) ≤ Real.sqrt 6 := Real.le_sqrt_of_sq_le (by norm_num)
  apply inv_anti₀ (by norm_num)
  norm_num
  nlinarith

/-- Turn a ball for the Euler–Maclaurin main term into a ball for ζ(1/2 + it). -/
lemma zeta_ball (t : ℝ) (ht : t ≠ 0) {A B R1 R2 : ℝ}
    (hmain : ‖1 + (2 : ℂ) ^ (-((1 / 2 : ℂ) + (t : ℂ) * I)) + (3 : ℂ) ^ (-((1 / 2 : ℂ) + (t : ℂ) * I))
      + (2 : ℂ) ^ (-((1 / 2 : ℂ) + (t : ℂ) * I)) * (2 : ℂ) ^ (-((1 / 2 : ℂ) + (t : ℂ) * I))
      + (5 : ℂ) ^ (-((1 / 2 : ℂ) + (t : ℂ) * I))
      + ((2 : ℂ) ^ (-((1 / 2 : ℂ) + (t : ℂ) * I)) * (3 : ℂ) ^ (-((1 / 2 : ℂ) + (t : ℂ) * I)))
        * EMK ((1 / 2 : ℂ) + (t : ℂ) * I) - (A + B * I)‖ ≤ R1)
    (hem : ‖(1 / 2 : ℂ) + (t : ℂ) * I‖ * (‖gc ((1 / 2 : ℂ) + (t : ℂ) * I) 7‖ / 1000000 *
      ((1 / (279936 * (2449 / 1000))) / (1 / 2 + 7))) ≤ R2) :
    ‖riemannZeta ((1 / 2 : ℂ) + (t : ℂ) * I) - (A + B * I)‖ ≤ R1 + R2 := by
  set s : ℂ := (1 / 2 : ℂ) + (t : ℂ) * I with hsdef
  have hre : s.re = 1 / 2 := by simp [hsdef]
  have hs0 : 0 < s.re := by rw [hre]; norm_num
  have hs1 : s ≠ 1 := by
    intro h; have := congrArg Complex.im h; simp [hsdef] at this; exact ht this
  have hsz : s ≠ 0 := by
    intro h; have := congrArg Complex.re h; rw [hre] at this; norm_num at this
  have hEM := zeta_EM8 hs0 hs1 (N := 6) (by norm_num)
  rw [em_main_eq s hsz hs1] at hEM
  have hb : ‖s‖ * (‖gc s 7‖ / 1000000 * (((6 : ℕ) : ℝ) ^ (-(s.re + 7)) / (s.re + 7))) ≤ R2 := by
    refine le_trans ?_ hem
    rw [hre]
    gcongr
    have := six_rpow_le
    norm_num at this ⊢
    exact this
  calc ‖riemannZeta s - (A + B * I)‖ ≤ ‖riemannZeta s - (1 + (2 : ℂ) ^ (-s) + (3 : ℂ) ^ (-s)
        + (2 : ℂ) ^ (-s) * (2 : ℂ) ^ (-s) + (5 : ℂ) ^ (-s) + ((2 : ℂ) ^ (-s) * (3 : ℂ) ^ (-s)) * EMK s)‖
        + ‖(1 + (2 : ℂ) ^ (-s) + (3 : ℂ) ^ (-s) + (2 : ℂ) ^ (-s) * (2 : ℂ) ^ (-s) + (5 : ℂ) ^ (-s)
        + ((2 : ℂ) ^ (-s) * (3 : ℂ) ^ (-s)) * EMK s) - (A + B * I)‖ := norm_sub_le_norm_sub_add_norm_sub _ _ _
    _ ≤ R2 + R1 := add_le_add (hEM.trans hb) hmain
    _ = R1 + R2 := add_comm _ _

end Fayed.FirstZero
