import Mathlib.NumberTheory.ArithmeticFunction.Misc
import Mathlib.Basic.Real.Basic
import Mathlib.Algebra.Order.BigOperators.Group.List

/-!
# Σ (−a)^{Ω(n)}: the 2-adic identity

`WN a N = Σ_{n=1}^{N} (−a)^{Ω(n)}`. Writing `n = 2^j m` with `m` odd,
`WN a N = Σ_{m ≤ N odd} (−a)^{Ω(m)} · Σ_{j : m 2^j ≤ N} (−a)^j` (`WN_eq`).
Also `3^{Ω(m)} ≤ m` and `2^{3Ω(m)} ≤ m²` for odd `m`.
-/

open Finset
open scoped ArithmeticFunction.Omega

namespace OmegaLean

/-- `WN a N = Σ_{n=1}^{N} (−a)^{Ω(n)}`. -/
noncomputable def WN (a : ℝ) (N : ℕ) : ℝ := ∑ n ∈ Icc 1 N, (-a) ^ (Ω n)

/-- Odd part: `ON a N = Σ_{n ≤ N, n odd} (−a)^{Ω(n)}`. -/
noncomputable def ON (a : ℝ) (N : ℕ) : ℝ := ∑ n ∈ Icc 1 N with Odd n, (-a) ^ (Ω n)

/-- Inner geometric sum `E a N m = Σ_{j < N, m 2^j ≤ N} (−a)^j`. -/
noncomputable def E (a : ℝ) (N m : ℕ) : ℝ := ∑ j ∈ range N with m * 2 ^ j ≤ N, (-a) ^ j

lemma cardFactors_two_mul {n : ℕ} (hn : n ≠ 0) : Ω (2 * n) = Ω n + 1 := by
  rw [ArithmeticFunction.cardFactors_mul two_ne_zero hn,
    ArithmeticFunction.cardFactors_apply_prime Nat.prime_two, add_comm]

lemma sum_even (a : ℝ) (N : ℕ) :
    ∑ n ∈ Icc 1 N with ¬ Odd n, (-a) ^ (Ω n) = (-a) * WN a (N / 2) := by
  rw [WN, mul_sum]
  symm
  refine sum_nbij' (fun n => 2 * n) (fun n => n / 2) ?_ ?_ ?_ ?_ ?_
  · intro n hn
    simp only [mem_Icc, mem_filter, Nat.not_odd_iff_even] at hn ⊢
    exact ⟨⟨by omega, by omega⟩, even_two_mul n⟩
  · intro n hn
    simp only [mem_Icc, mem_filter, Nat.not_odd_iff_even] at hn ⊢
    obtain ⟨⟨h1, h2⟩, ⟨r, hr⟩⟩ := hn
    omega
  · intro n _
    omega
  · intro n hn
    simp only [mem_Icc, mem_filter, Nat.not_odd_iff_even] at hn
    obtain ⟨_, ⟨r, hr⟩⟩ := hn
    omega
  · intro n hn
    have h1 : 1 ≤ n := (mem_Icc.mp hn).1
    rw [cardFactors_two_mul (by omega), pow_succ]
    ring

lemma WN_rec (a : ℝ) (N : ℕ) : WN a N = ON a N + (-a) * WN a (N / 2) := by
  rw [← sum_even, ON]
  unfold WN
  exact (sum_filter_add_sum_filter_not _ _ _).symm

lemma WN_zero (a : ℝ) : WN a 0 = 0 := by simp [WN]

lemma WN_iter (a : ℝ) (N : ℕ) : ∀ J : ℕ,
    WN a N = ∑ j ∈ range J, (-a) ^ j * ON a (N / 2 ^ j) + (-a) ^ J * WN a (N / 2 ^ J)
  | 0 => by simp
  | J + 1 => by
    have h2 : N / 2 ^ (J + 1) = N / 2 ^ J / 2 := by rw [Nat.div_div_eq_div_mul, pow_succ]
    rw [WN_iter a N J, sum_range_succ, WN_rec a (N / 2 ^ J), h2, pow_succ]
    ring

lemma WN_eq_sum_ON (a : ℝ) (N : ℕ) :
    WN a N = ∑ j ∈ range N, (-a) ^ j * ON a (N / 2 ^ j) := by
  rw [WN_iter a N N, Nat.div_eq_of_lt Nat.lt_two_pow_self, WN_zero, mul_zero, add_zero]

/-- The 2-adic identity. -/
theorem WN_eq (a : ℝ) (N : ℕ) :
    WN a N = ∑ m ∈ Icc 1 N with Odd m, (-a) ^ (Ω m) * E a N m := by
  have hmem : ∀ (j m : ℕ), j ∈ range N ∧ m ∈ (Icc 1 (N / 2 ^ j)).filter Odd ↔
      j ∈ (range N).filter (fun j => m * 2 ^ j ≤ N) ∧ m ∈ (Icc 1 N).filter Odd := by
    intro j m
    have key : m ≤ N / 2 ^ j ↔ m * 2 ^ j ≤ N := Nat.le_div_iff_mul_le (by positivity)
    have hle : m * 2 ^ j ≤ N → m ≤ N := fun h =>
      le_trans (Nat.le_mul_of_pos_right m (by positivity)) h
    simp only [mem_range, mem_filter, mem_Icc, key]
    constructor
    · rintro ⟨hj, ⟨h1, h2⟩, hodd⟩
      exact ⟨⟨hj, h2⟩, ⟨h1, hle h2⟩, hodd⟩
    · rintro ⟨⟨hj, h2⟩, ⟨h1, -⟩, hodd⟩
      exact ⟨hj, ⟨h1, h2⟩, hodd⟩
  rw [WN_eq_sum_ON]
  simp only [ON, E, mul_sum]
  rw [sum_comm' hmem]
  exact sum_congr rfl fun m _ => sum_congr rfl fun j _ => mul_comm _ _

/-- For odd `m`, `3^{Ω(m)} ≤ m`. -/
lemma three_pow_cardFactors_le {m : ℕ} (hm : Odd m) : 3 ^ (Ω m) ≤ m := by
  have hm0 : m ≠ 0 := by
    obtain ⟨r, hr⟩ := hm
    omega
  rw [ArithmeticFunction.cardFactors_apply]
  conv_rhs => rw [← Nat.prod_primeFactorsList hm0]
  apply List.pow_length_le_prod
  intro p hp
  have hpr := Nat.prime_of_mem_primeFactorsList hp
  have hpd := Nat.dvd_of_mem_primeFactorsList hp
  have hp2 : p ≠ 2 := by
    rintro rfl
    exact (Nat.not_even_iff_odd.mpr hm) (even_iff_two_dvd.mpr hpd)
  have := hpr.two_le
  omega

/-- For odd `m`, `2^{3Ω(m)} ≤ m²` (from `8 ≤ 9`). -/
lemma two_pow_three_mul_le {m : ℕ} (hm : Odd m) : 2 ^ (3 * Ω m) ≤ m ^ 2 := by
  have h := three_pow_cardFactors_le hm
  calc 2 ^ (3 * Ω m) = 8 ^ (Ω m) := by rw [pow_mul]; norm_num
    _ ≤ 9 ^ (Ω m) := Nat.pow_le_pow_left (by norm_num) _
    _ = (3 ^ (Ω m)) ^ 2 := by rw [← pow_mul, mul_comm, pow_mul]; norm_num
    _ ≤ m ^ 2 := Nat.pow_le_pow_left h 2

end OmegaLean
