import Riemann.Deep.ConjFals.Lean.Identity
import Mathlib.Analysis.SpecialFunctions.Pow.Real
import Mathlib.Analysis.SpecialFunctions.Log.Base
import Mathlib.Analysis.PSeries
import Mathlib.Algebra.Field.GeomSum

/-!
# Explicit bounds for Σ (−a)^{Ω(n)}

For `a > 1` with `α = log₂ a ≥ 6` (e.g. `a ≥ 64`):
* `abs_WN_le`: `|WN a N| ≤ (a/(a−1)) (1 + 9a^{−1/3}) N^α` for every `N`;
* `abs_WN_two_pow_ge`: `|WN a 2^k| ≥ (a^{k+1} − 1)/(a+1) − (a/(a−1)) 9a^{−1/3} a^k`.
Key inputs: the 2-adic identity `WN_eq`, `Σ_{j∈F} a^j ≤ X a/(a−1)` when every `a^j ≤ X`, and
`a^{Ω(m)} m^{−α} ≤ 9 a^{−1/3} m^{−2}` for odd `m ≥ 3` (from `2^{3Ω(m)} ≤ m²`).
-/

open Finset Real
open scoped ArithmeticFunction.Omega

namespace OmegaLean

/-- A sum of distinct powers `a^j ≤ X` is at most `X · a/(a−1)`. -/
lemma sum_pow_le {a X : ℝ} (ha : 1 < a) (hX : 0 ≤ X) (F : Finset ℕ)
    (hF : ∀ j ∈ F, a ^ j ≤ X) : ∑ j ∈ F, a ^ j ≤ X * (a / (a - 1)) := by
  have ha0 : 0 < a := by linarith
  have ha1 : 0 < a - 1 := by linarith
  have hq : 0 ≤ a / (a - 1) := div_nonneg ha0.le ha1.le
  rcases F.eq_empty_or_nonempty with rfl | hne
  · rw [sum_empty]; exact mul_nonneg hX hq
  · have hsub : F ⊆ range (F.max' hne + 1) := fun j hj =>
      mem_range.mpr (Nat.lt_succ_of_le (F.le_max' j hj))
    have hK : a ^ (F.max' hne) ≤ X := hF _ (F.max'_mem hne)
    calc ∑ j ∈ F, a ^ j ≤ ∑ j ∈ range (F.max' hne + 1), a ^ j :=
          sum_le_sum_of_subset_of_nonneg hsub (fun j _ _ => by positivity)
      _ = (a ^ (F.max' hne + 1) - 1) / (a - 1) := geom_sum_eq (ne_of_gt ha) _
      _ ≤ a ^ (F.max' hne + 1) / (a - 1) := div_le_div_of_nonneg_right (by linarith) ha1.le
      _ = a ^ (F.max' hne) * (a / (a - 1)) := by rw [pow_succ]; ring
      _ ≤ X * (a / (a - 1)) := mul_le_mul_of_nonneg_right hK hq

/-- `(2^j)^{log₂ a} = a^j`. -/
lemma two_pow_rpow_logb {a : ℝ} (ha : 0 < a) (j : ℕ) :
    ((2:ℝ) ^ j) ^ (logb 2 a) = a ^ j := by
  have h2a : (2:ℝ) ^ (logb 2 a) = a := Real.rpow_logb (by norm_num) (by norm_num) ha
  rw [← Real.rpow_natCast (2:ℝ) j, ← Real.rpow_mul (by norm_num : (0:ℝ) ≤ 2), mul_comm,
    Real.rpow_mul (by norm_num : (0:ℝ) ≤ 2), h2a, Real.rpow_natCast]

/-- If `m 2^j ≤ N` then `a^j ≤ N^α m^{−α}`. -/
lemma pow_le_of_mul_two_pow_le {a : ℝ} (ha : 1 ≤ a) {N m j : ℕ} (hm : 1 ≤ m)
    (h : m * 2 ^ j ≤ N) :
    a ^ j ≤ (N : ℝ) ^ (logb 2 a) * (m : ℝ) ^ (-logb 2 a) := by
  have ha0 : 0 < a := by linarith
  have hα : 0 ≤ logb 2 a := Real.logb_nonneg (by norm_num) ha
  have hm0 : (0:ℝ) < m := by exact_mod_cast hm
  have hcast : (m : ℝ) * 2 ^ j ≤ N := by exact_mod_cast h
  have h1 : ((m : ℝ) * 2 ^ j) ^ (logb 2 a) ≤ (N : ℝ) ^ (logb 2 a) :=
    Real.rpow_le_rpow (by positivity) hcast hα
  rw [Real.mul_rpow hm0.le (by positivity), two_pow_rpow_logb ha0] at h1
  have hmα : 0 < (m : ℝ) ^ (logb 2 a) := Real.rpow_pos_of_pos hm0 _
  rw [Real.rpow_neg hm0.le, ← div_eq_mul_inv, le_div_iff₀ hmα]
  linarith [h1]

/-- `|E a N m| ≤ N^α m^{−α} · a/(a−1)`. -/
lemma abs_E_le {a : ℝ} (ha : 1 < a) {N m : ℕ} (hm : 1 ≤ m) :
    |E a N m| ≤ ((N : ℝ) ^ (logb 2 a) * (m : ℝ) ^ (-logb 2 a)) * (a / (a - 1)) := by
  unfold E
  refine (abs_sum_le_sum_abs _ _).trans ?_
  have habs : ∀ j ∈ (range N).filter (fun j => m * 2 ^ j ≤ N), |(-a) ^ j| = a ^ j := by
    intro j _
    rw [abs_pow, abs_neg, abs_of_pos (by linarith)]
  rw [sum_congr rfl habs]
  apply sum_pow_le ha (by positivity)
  intro j hj
  exact pow_le_of_mul_two_pow_le ha.le hm (mem_filter.mp hj).2

/-- One term of the 2-adic identity. -/
lemma abs_term_le {a : ℝ} (ha : 1 < a) {N m : ℕ} (hm : 1 ≤ m) :
    |(-a) ^ (Ω m) * E a N m| ≤
      (a ^ (Ω m) * (m : ℝ) ^ (-logb 2 a)) * ((N : ℝ) ^ (logb 2 a) * (a / (a - 1))) := by
  have ha0 : 0 < a := by linarith
  rw [abs_mul, abs_pow, abs_neg, abs_of_pos ha0]
  calc a ^ (Ω m) * |E a N m|
      ≤ a ^ (Ω m) * (((N : ℝ) ^ (logb 2 a) * (m : ℝ) ^ (-logb 2 a)) * (a / (a - 1))) :=
        mul_le_mul_of_nonneg_left (abs_E_le ha hm) (by positivity)
    _ = (a ^ (Ω m) * (m : ℝ) ^ (-logb 2 a)) * ((N : ℝ) ^ (logb 2 a) * (a / (a - 1))) := by ring

/-- For odd `m ≥ 3` and `log₂ a ≥ 6`: `a^{Ω(m)} m^{−α} ≤ 9 a^{−1/3} m^{−2}`. -/
lemma term_le {a : ℝ} (ha : 1 < a) (hα : 6 ≤ logb 2 a) {m : ℕ} (hm : Odd m) (hm3 : 3 ≤ m) :
    a ^ (Ω m) * (m : ℝ) ^ (-logb 2 a) ≤ 9 * a ^ (-(1/3 : ℝ)) * ((m : ℝ) ^ 2)⁻¹ := by
  have ha0 : 0 < a := by linarith
  have hm0 : (0 : ℝ) < m := by exact_mod_cast (show 0 < m by omega)
  have hl2 : 0 < Real.log 2 := Real.log_pos (by norm_num)
  have hL : Real.log a = logb 2 a * Real.log 2 := by
    rw [Real.logb, div_mul_cancel₀ _ hl2.ne']
  have h8 : 3 * (Ω m : ℝ) * Real.log 2 ≤ 2 * Real.log m := by
    have h' : ((2 : ℝ) ^ (3 * Ω m)) ≤ (m : ℝ) ^ 2 := by exact_mod_cast two_pow_three_mul_le hm
    have := Real.log_le_log (by positivity) h'
    rw [Real.log_pow, Real.log_pow] at this
    push_cast at this
    linarith
  have hu : Real.log 3 ≤ Real.log m := Real.log_le_log (by norm_num) (by exact_mod_cast hm3)
  have h23 : Real.log 2 < Real.log 3 := Real.log_lt_log (by norm_num) (by norm_num)
  have hlog9 : Real.log 9 = 2 * Real.log 3 := by
    rw [show (9:ℝ) = 3 ^ 2 by norm_num, Real.log_pow]; norm_num
  have p1 : 0 < a ^ (Ω m) := pow_pos ha0 _
  have p2 : 0 < (m : ℝ) ^ (-logb 2 a) := Real.rpow_pos_of_pos hm0 _
  have p3 : 0 < a ^ (-(1/3 : ℝ)) := Real.rpow_pos_of_pos ha0 _
  have p4 : 0 < ((m : ℝ) ^ 2)⁻¹ := inv_pos.mpr (pow_pos hm0 2)
  have p5 : (0 : ℝ) < 9 := by norm_num
  rw [← Real.log_le_log_iff (mul_pos p1 p2) (mul_pos (mul_pos p5 p3) p4),
    Real.log_mul p1.ne' p2.ne', Real.log_mul (mul_pos p5 p3).ne' p4.ne',
    Real.log_mul p5.ne' p3.ne', Real.log_pow, Real.log_rpow ha0, Real.log_rpow hm0,
    Real.log_inv, Real.log_pow, hlog9, hL]
  have key : logb 2 a * ((Ω m : ℝ) * Real.log 2) ≤ logb 2 a * (2 * Real.log m / 3) :=
    mul_le_mul_of_nonneg_left (by linarith) (by linarith)
  have hprod : 0 ≤ (logb 2 a / 3 - 2) * (Real.log m - Real.log 2) :=
    mul_nonneg (by linarith) (by linarith)
  push_cast
  linarith [key, hprod, h23]

/-- `Σ_{m=2}^{N} m^{−2} ≤ 1`. -/
lemma sum_Icc_two_inv_sq_le (N : ℕ) : ∑ m ∈ Icc 2 N, ((m : ℝ) ^ 2)⁻¹ ≤ 1 := by
  have h := sum_Ioo_inv_sq_le (α := ℝ) 1 (N + 1)
  have he : Ioo 1 (N + 1) = Icc 2 N := by
    ext x; simp only [mem_Ioo, mem_Icc]; omega
  rw [he] at h
  calc ∑ m ∈ Icc 2 N, ((m : ℝ) ^ 2)⁻¹ ≤ 2 / ((1 : ℕ) + 1 : ℝ) := h
    _ = 1 := by norm_num

/-- The odd terms `m ≥ 3` contribute at most `9 a^{−1/3}`. -/
lemma tail_le {a : ℝ} (ha : 1 < a) (hα : 6 ≤ logb 2 a) (N : ℕ) :
    ∑ m ∈ Icc 2 N with Odd m, a ^ (Ω m) * (m : ℝ) ^ (-logb 2 a) ≤ 9 * a ^ (-(1/3 : ℝ)) := by
  have ha0 : 0 < a := by linarith
  have hc : 0 ≤ 9 * a ^ (-(1/3 : ℝ)) := by positivity
  calc ∑ m ∈ Icc 2 N with Odd m, a ^ (Ω m) * (m : ℝ) ^ (-logb 2 a)
      ≤ ∑ m ∈ Icc 2 N with Odd m, 9 * a ^ (-(1/3 : ℝ)) * ((m : ℝ) ^ 2)⁻¹ := by
        apply sum_le_sum
        intro m hm
        rw [mem_filter, mem_Icc] at hm
        have h3 : 3 ≤ m := by
          obtain ⟨⟨h1, _⟩, ⟨r, hr⟩⟩ := hm
          omega
        exact term_le ha hα hm.2 h3
    _ ≤ ∑ m ∈ Icc 2 N, 9 * a ^ (-(1/3 : ℝ)) * ((m : ℝ) ^ 2)⁻¹ :=
        sum_le_sum_of_subset_of_nonneg (filter_subset _ _) (fun m _ _ => by positivity)
    _ = 9 * a ^ (-(1/3 : ℝ)) * ∑ m ∈ Icc 2 N, ((m : ℝ) ^ 2)⁻¹ := by rw [mul_sum]
    _ ≤ 9 * a ^ (-(1/3 : ℝ)) * 1 := mul_le_mul_of_nonneg_left (sum_Icc_two_inv_sq_le N) hc
    _ = 9 * a ^ (-(1/3 : ℝ)) := mul_one _

/-- `Σ_{m ≤ N odd} a^{Ω(m)} m^{−α} ≤ 1 + 9 a^{−1/3}`. -/
lemma sum_odd_le {a : ℝ} (ha : 1 < a) (hα : 6 ≤ logb 2 a) (N : ℕ) :
    ∑ m ∈ Icc 1 N with Odd m, a ^ (Ω m) * (m : ℝ) ^ (-logb 2 a) ≤ 1 + 9 * a ^ (-(1/3 : ℝ)) := by
  have ha0 : 0 < a := by linarith
  have hsub : (Icc 1 N).filter (fun m => Odd m) ⊆ insert 1 ((Icc 2 N).filter (fun m => Odd m)) := by
    intro m hm
    rw [mem_filter, mem_Icc] at hm
    rw [mem_insert, mem_filter, mem_Icc]
    by_cases h : m = 1
    · left; exact h
    · right; exact ⟨⟨by omega, hm.1.2⟩, hm.2⟩
  have h1 : (1:ℕ) ∉ (Icc 2 N).filter (fun m => Odd m) := by simp
  calc ∑ m ∈ Icc 1 N with Odd m, a ^ (Ω m) * (m : ℝ) ^ (-logb 2 a)
      ≤ ∑ m ∈ insert 1 ((Icc 2 N).filter (fun m => Odd m)), a ^ (Ω m) * (m : ℝ) ^ (-logb 2 a) :=
        sum_le_sum_of_subset_of_nonneg hsub (fun m _ _ => by positivity)
    _ = 1 + ∑ m ∈ Icc 2 N with Odd m, a ^ (Ω m) * (m : ℝ) ^ (-logb 2 a) := by
        rw [sum_insert h1]; simp
    _ ≤ 1 + 9 * a ^ (-(1/3 : ℝ)) := by linarith [tail_le ha hα N]

/-- Upper bound: `|WN a N| ≤ (a/(a−1)) (1 + 9a^{−1/3}) N^α`. -/
theorem abs_WN_le {a : ℝ} (ha : 1 < a) (hα : 6 ≤ logb 2 a) (N : ℕ) :
    |WN a N| ≤ (a / (a - 1)) * (1 + 9 * a ^ (-(1/3 : ℝ))) * (N : ℝ) ^ (logb 2 a) := by
  rw [WN_eq]
  have hC : 0 ≤ (N : ℝ) ^ (logb 2 a) * (a / (a - 1)) :=
    mul_nonneg (by positivity) (div_nonneg (by linarith) (by linarith))
  calc |∑ m ∈ Icc 1 N with Odd m, (-a) ^ (Ω m) * E a N m|
      ≤ ∑ m ∈ Icc 1 N with Odd m, |(-a) ^ (Ω m) * E a N m| := abs_sum_le_sum_abs _ _
    _ ≤ ∑ m ∈ Icc 1 N with Odd m,
          (a ^ (Ω m) * (m : ℝ) ^ (-logb 2 a)) * ((N : ℝ) ^ (logb 2 a) * (a / (a - 1))) := by
        apply sum_le_sum
        intro m hm
        exact abs_term_le ha (mem_Icc.mp (mem_filter.mp hm).1).1
    _ = (∑ m ∈ Icc 1 N with Odd m, a ^ (Ω m) * (m : ℝ) ^ (-logb 2 a)) *
          ((N : ℝ) ^ (logb 2 a) * (a / (a - 1))) := by rw [sum_mul]
    _ ≤ (1 + 9 * a ^ (-(1/3 : ℝ))) * ((N : ℝ) ^ (logb 2 a) * (a / (a - 1))) :=
        mul_le_mul_of_nonneg_right (sum_odd_le ha hα N) hC
    _ = (a / (a - 1)) * (1 + 9 * a ^ (-(1/3 : ℝ))) * (N : ℝ) ^ (logb 2 a) := by ring

/-- `E a 2^k 1 = Σ_{j=0}^{k} (−a)^j`. -/
lemma E_one_two_pow (a : ℝ) (k : ℕ) : E a (2 ^ k) 1 = ∑ j ∈ range (k + 1), (-a) ^ j := by
  unfold E
  congr 1
  ext j
  simp only [mem_filter, mem_range, one_mul]
  constructor
  · rintro ⟨_, h⟩
    have := (Nat.pow_le_pow_iff_right (by norm_num : 1 < 2)).mp h
    omega
  · intro h
    have hk : k < 2 ^ k := Nat.lt_two_pow_self
    exact ⟨by omega, Nat.pow_le_pow_right (by norm_num) (by omega)⟩

/-- `|E a 2^k 1| ≥ (a^{k+1} − 1)/(a+1)`. -/
lemma abs_E_one_ge {a : ℝ} (ha : 1 < a) (k : ℕ) :
    (a ^ (k + 1) - 1) / (a + 1) ≤ |E a (2 ^ k) 1| := by
  have ha1 : 0 < a + 1 := by linarith
  rw [E_one_two_pow, geom_sum_eq (ne_of_lt (by linarith : -a < 1)), abs_div,
    show (-a - 1) = -(a + 1) by ring, abs_neg, abs_of_pos ha1]
  apply div_le_div_of_nonneg_right _ ha1.le
  have := abs_sub_abs_le_abs_sub ((-a) ^ (k + 1)) 1
  rw [abs_pow, abs_neg, abs_of_pos (by linarith : (0:ℝ) < a), abs_one] at this
  exact this

/-- Lower bound at `N = 2^k`. -/
theorem abs_WN_two_pow_ge {a : ℝ} (ha : 1 < a) (hα : 6 ≤ logb 2 a) (k : ℕ) :
    (a ^ (k + 1) - 1) / (a + 1) - (a / (a - 1)) * (9 * a ^ (-(1/3 : ℝ))) * a ^ k ≤
      |WN a (2 ^ k)| := by
  have hN : 1 ≤ 2 ^ k := Nat.one_le_two_pow
  have hsplit : (Icc 1 (2 ^ k)).filter (fun m => Odd m) =
      insert 1 ((Icc 2 (2 ^ k)).filter (fun m => Odd m)) := by
    ext m
    simp only [mem_filter, mem_Icc, mem_insert]
    constructor
    · rintro ⟨⟨h1, h2⟩, ho⟩
      by_cases h : m = 1
      · left; exact h
      · right; exact ⟨⟨by omega, h2⟩, ho⟩
    · rintro (rfl | ⟨⟨h1, h2⟩, ho⟩)
      · exact ⟨⟨le_refl _, hN⟩, odd_one⟩
      · exact ⟨⟨by omega, h2⟩, ho⟩
  have h1 : (1:ℕ) ∉ (Icc 2 (2 ^ k)).filter (fun m => Odd m) := by simp
  rw [WN_eq, hsplit, sum_insert h1, ArithmeticFunction.cardFactors_one, pow_zero, one_mul]
  set B := ∑ m ∈ Icc 2 (2 ^ k) with Odd m, (-a) ^ (Ω m) * E a (2 ^ k) m with hBdef
  have hC : 0 ≤ ((2 ^ k : ℕ) : ℝ) ^ (logb 2 a) * (a / (a - 1)) :=
    mul_nonneg (by positivity) (div_nonneg (by linarith) (by linarith))
  have hpow : ((2 ^ k : ℕ) : ℝ) ^ (logb 2 a) = a ^ k := by
    push_cast; exact two_pow_rpow_logb (by linarith) k
  have hB : |B| ≤ (a / (a - 1)) * (9 * a ^ (-(1/3 : ℝ))) * a ^ k := by
    calc |B| ≤ ∑ m ∈ Icc 2 (2 ^ k) with Odd m, |(-a) ^ (Ω m) * E a (2 ^ k) m| :=
          abs_sum_le_sum_abs _ _
      _ ≤ ∑ m ∈ Icc 2 (2 ^ k) with Odd m, (a ^ (Ω m) * (m : ℝ) ^ (-logb 2 a)) *
            (((2 ^ k : ℕ) : ℝ) ^ (logb 2 a) * (a / (a - 1))) := by
          apply sum_le_sum
          intro m hm
          have hm1 : 1 ≤ m := by
            have := (mem_Icc.mp (mem_filter.mp hm).1).1
            omega
          exact abs_term_le ha hm1
      _ = (∑ m ∈ Icc 2 (2 ^ k) with Odd m, a ^ (Ω m) * (m : ℝ) ^ (-logb 2 a)) *
            (((2 ^ k : ℕ) : ℝ) ^ (logb 2 a) * (a / (a - 1))) := by rw [sum_mul]
      _ ≤ (9 * a ^ (-(1/3 : ℝ))) * (((2 ^ k : ℕ) : ℝ) ^ (logb 2 a) * (a / (a - 1))) :=
          mul_le_mul_of_nonneg_right (tail_le ha hα _) hC
      _ = (a / (a - 1)) * (9 * a ^ (-(1/3 : ℝ))) * a ^ k := by rw [hpow]; ring
  have hA := abs_E_one_ge ha k
  have htri := abs_sub_abs_le_abs_sub (E a (2 ^ k) 1) (-B)
  rw [abs_neg, sub_neg_eq_add] at htri
  linarith [hA, hB, htri]

end OmegaLean
