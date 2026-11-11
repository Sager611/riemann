import Riemann.Deep.ConjFals.Lean.Identity
import Mathlib.Algebra.Field.GeomSum
import Mathlib.Data.Nat.Log
import Mathlib.Tactic.LinearCombination

/-!
# The 2-adic identity in closed form (step (1) of the brief)

For `a ≠ −1` and every `N`:
`WN a N = (1/(1+a)) Σ_{m ≤ N odd} (−a)^{Ω(m)} + (a/(1+a)) Σ_{m ≤ N odd} (−a)^{Ω(m)} (−a)^{K(m)}`
with `K(m) = Nat.log 2 (N / m)`, the largest `j` with `m 2^j ≤ N` (`mul_two_pow_le_iff`), i.e.
`K(m) = ⌊log₂(x/m)⌋` for `N = ⌊x⌋`. Not needed for the main theorem.
-/

open Finset
open scoped ArithmeticFunction.Omega

namespace OmegaLean

/-- `Nat.log 2 (N/m)` is the largest `j` with `m 2^j ≤ N`. -/
lemma mul_two_pow_le_iff {N m j : ℕ} (hm : 1 ≤ m) (hmN : m ≤ N) :
    m * 2 ^ j ≤ N ↔ j ≤ Nat.log 2 (N / m) := by
  have hpos : 0 < N / m := Nat.div_pos hmN (by omega)
  rw [Nat.le_log_iff_pow_le (by norm_num) hpos.ne', Nat.le_div_iff_mul_le (by omega), mul_comm]

lemma E_eq_range (a : ℝ) {N m : ℕ} (hm : 1 ≤ m) (hmN : m ≤ N) :
    E a N m = ∑ j ∈ range (Nat.log 2 (N / m) + 1), (-a) ^ j := by
  unfold E
  congr 1
  ext j
  simp only [mem_filter, mem_range]
  rw [mul_two_pow_le_iff hm hmN]
  constructor
  · rintro ⟨_, h⟩
    omega
  · intro h
    refine ⟨?_, by omega⟩
    have h1 : m * 2 ^ j ≤ N := (mul_two_pow_le_iff hm hmN).mpr (by omega)
    have h2 : j < 2 ^ j := Nat.lt_two_pow_self
    have h3 : 2 ^ j ≤ m * 2 ^ j := Nat.le_mul_of_pos_left (2 ^ j) (by omega)
    omega

lemma E_closed {a : ℝ} (ha : a ≠ -1) {N m : ℕ} (hm : 1 ≤ m) (hmN : m ≤ N) :
    E a N m = 1 / (1 + a) + a / (1 + a) * (-a) ^ (Nat.log 2 (N / m)) := by
  have h1 : (1 + a) ≠ 0 := by
    intro h; apply ha; linarith
  have h2 : (-a) - 1 ≠ 0 := by
    intro h; apply ha; linarith
  rw [E_eq_range a hm hmN, geom_sum_eq (by intro h; apply ha; linarith), pow_succ,
    div_eq_iff h2]
  have hw : (1 + a) * (1 + a)⁻¹ = 1 := mul_inv_cancel₀ h1
  linear_combination (1 + a * (-a) ^ (Nat.log 2 (N / m))) * hw

/-- The 2-adic identity in closed form. -/
theorem WN_two_adic {a : ℝ} (ha : a ≠ -1) (N : ℕ) :
    WN a N = 1 / (1 + a) * ON a N +
      a / (1 + a) * ∑ m ∈ Icc 1 N with Odd m, (-a) ^ (Ω m) * (-a) ^ (Nat.log 2 (N / m)) := by
  rw [WN_eq, ON, mul_sum, mul_sum, ← sum_add_distrib]
  apply sum_congr rfl
  intro m hm
  have hm' := mem_Icc.mp (mem_filter.mp hm).1
  rw [E_closed ha hm'.1 hm'.2]
  ring

end OmegaLean
