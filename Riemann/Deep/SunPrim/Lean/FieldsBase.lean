import Mathlib.Data.ZMod.Basic
import Mathlib.Algebra.Field.ZMod
import Mathlib.Data.Set.Card
import Mathlib.Algebra.BigOperators.Group.Finset.Basic
import Mathlib.Order.Interval.Finset.Nat

/-!
# Weighted sums of distinct elements: shared definitions for 1811.10503 Conj. 4.9 and MathOverflow 316143

Both conjectures are stated for every field `F` (with `p(F)` as in the sources) and refuted in
`F = ℤ/11`, `ℤ/13`: the sum set is contained in an explicit finite set that is too small (`decide`).
Sources and fidelity: `../LEAN.md`.
-/

namespace SunLean

open Finset

/-- `p(F)`: the characteristic if it is a prime `p`, and `+∞` if the characteristic is zero. -/
noncomputable def pF (F : Type) [Field F] : ℕ∞ := if ringChar F = 0 then ⊤ else (ringChar F : ℕ∞)

theorem pF_zmod (p : ℕ) [Fact p.Prime] : pF (ZMod p) = (p : ℕ∞) := by
  unfold pF
  simp [ZMod.ringChar_zmod_n, (Fact.out : p.Prime).ne_zero]

/-- `{∑_{k=1}^n w_k a_k : a_1, …, a_n distinct elements of A}`. -/
def WSumSet {F : Type} [Field F] (n : ℕ) (w : ℕ → ℕ) (A : Finset F) : Set F :=
  {s | ∃ a : ℕ → F, (∀ k ∈ Icc 1 n, a k ∈ A) ∧ (∀ i ∈ Icc 1 n, ∀ j ∈ Icc 1 n, a i = a j → i = j) ∧
    s = ∑ k ∈ Icc 1 n, (w k : F) * a k}

theorem sum_Icc_one_three {M : Type} [AddCommMonoid M] (f : ℕ → M) :
    ∑ k ∈ Icc 1 3, f k = f 1 + f 2 + f 3 := by
  rw [show Icc 1 3 = {1, 2, 3} by decide, sum_insert (by decide), sum_insert (by decide),
    sum_singleton, add_assoc]

theorem sum_Icc_one_four {M : Type} [AddCommMonoid M] (f : ℕ → M) :
    ∑ k ∈ Icc 1 4, f k = f 1 + f 2 + f 3 + f 4 := by
  rw [show Icc 1 4 = {1, 2, 3, 4} by decide, sum_insert (by decide), sum_insert (by decide),
    sum_insert (by decide), sum_singleton, add_assoc, add_assoc]

theorem WSumSet_three_sub {F : Type} [Field F] (w : ℕ → ℕ) (A T : Finset F)
    (hT : ∀ x ∈ A, ∀ y ∈ A, ∀ z ∈ A, x ≠ y → x ≠ z → y ≠ z →
      (w 1 : F) * x + (w 2 : F) * y + (w 3 : F) * z ∈ T) :
    WSumSet 3 w A ⊆ T := by
  rintro s ⟨a, hA, hinj, rfl⟩
  rw [sum_Icc_one_three]
  have h1 : 1 ∈ Icc 1 3 := by decide
  have h2 : 2 ∈ Icc 1 3 := by decide
  have h3 : 3 ∈ Icc 1 3 := by decide
  exact hT _ (hA 1 h1) _ (hA 2 h2) _ (hA 3 h3) (fun e => absurd (hinj 1 h1 2 h2 e) (by decide))
    (fun e => absurd (hinj 1 h1 3 h3 e) (by decide)) (fun e => absurd (hinj 2 h2 3 h3 e) (by decide))

theorem WSumSet_four_sub {F : Type} [Field F] (w : ℕ → ℕ) (A T : Finset F)
    (hT : ∀ x ∈ A, ∀ y ∈ A, ∀ z ∈ A, ∀ u ∈ A, x ≠ y → x ≠ z → x ≠ u → y ≠ z → y ≠ u → z ≠ u →
      (w 1 : F) * x + (w 2 : F) * y + (w 3 : F) * z + (w 4 : F) * u ∈ T) :
    WSumSet 4 w A ⊆ T := by
  rintro s ⟨a, hA, hinj, rfl⟩
  rw [sum_Icc_one_four]
  have h1 : 1 ∈ Icc 1 4 := by decide
  have h2 : 2 ∈ Icc 1 4 := by decide
  have h3 : 3 ∈ Icc 1 4 := by decide
  have h4 : 4 ∈ Icc 1 4 := by decide
  exact hT _ (hA 1 h1) _ (hA 2 h2) _ (hA 3 h3) _ (hA 4 h4)
    (fun e => absurd (hinj 1 h1 2 h2 e) (by decide)) (fun e => absurd (hinj 1 h1 3 h3 e) (by decide))
    (fun e => absurd (hinj 1 h1 4 h4 e) (by decide)) (fun e => absurd (hinj 2 h2 3 h3 e) (by decide))
    (fun e => absurd (hinj 2 h2 4 h4 e) (by decide)) (fun e => absurd (hinj 3 h3 4 h4 e) (by decide))

instance fact_prime_11 : Fact (Nat.Prime 11) := ⟨by decide⟩
instance fact_prime_13 : Fact (Nat.Prime 13) := ⟨by decide⟩

end SunLean
