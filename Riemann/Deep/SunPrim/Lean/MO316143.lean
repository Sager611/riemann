import Riemann.Deep.SunPrim.Lean.FieldsBase

/-! Sun, MathOverflow question 316143: refuted in `ℤ/11` with `λ = (1, 1, 2)`. -/

namespace SunLean

open Finset

/-- `∑_{k=1}^{⌊n/2⌋} (n + 1 - 2k)(λ_{n+1-k} - λ_k)` and the bound of MathOverflow 316143, with
`m = |A|`: `(λ_1 + … + λ_n)(m - n) + ∑_{k=1}^{⌊n/2⌋} (n + 1 - 2k)(λ_{n+1-k} - λ_k) + 1` (in `ℤ`). -/
def bound316143 (n : ℕ) (lam : ℕ → ℕ) (m : ℕ) : ℤ :=
  (∑ k ∈ Icc 1 n, (lam k : ℤ)) * ((m : ℤ) - n) +
    (∑ k ∈ Icc 1 (n / 2), ((n : ℤ) + 1 - 2 * k) * ((lam (n + 1 - k) : ℤ) - lam k)) + 1

/-- Sanity checks: the instance's bound is 7; for `λ_k = k` the bound is the one of 1811.10503 Conj. 4.9
(`(|A| - n) n(n+1)/2 + n(n²-1)/6 + 1`: 11 for `n = |A| = 4`, 51 for `n = 5`, `|A| = 7`), as Sun remarks. -/
theorem bound316143_check : bound316143 3 (fun k => if k = 3 then 2 else 1) 4 = 7 ∧
    bound316143 4 (fun k => k) 4 = 11 ∧ bound316143 5 (fun k => k) 7 = 51 := by decide

/-- Z.-W. Sun, MathOverflow question 316143 (2018-11-25), Conjecture: `λ_1, …, λ_n` (`n ≥ 3`)
positive integers with `gcd(λ_1, …, λ_n) = 1` (only the divisor 1 is common) and
`λ_{k+1} - λ_k ∈ {0, 1}` (`0 < k < n`); `F` a field with `p(F) > n + 1`; then for every finite `A ⊆ F`
with `|A| ≥ n + δ_{n,3}`: `|{∑ λ_k a_k : a_1, …, a_n distinct elements of A}| ≥ min{p(F), bound}`. -/
def Conj_MO316143 : Prop :=
  ∀ (n : ℕ) (lam : ℕ → ℕ), 3 ≤ n → (∀ k ∈ Icc 1 n, 0 < lam k) →
    (∀ d : ℕ, (∀ k ∈ Icc 1 n, d ∣ lam k) → d = 1) →
    (∀ k ∈ Icc 1 (n - 1), (lam (k + 1) : ℤ) - lam k = 0 ∨ (lam (k + 1) : ℤ) - lam k = 1) →
    ∀ (F : Type) [Field F], ((n + 1 : ℕ) : ℕ∞) < pF F →
    ∀ A : Finset F, n + (if n = 3 then 1 else 0) ≤ A.card →
      min (pF F) ((bound316143 n lam A.card).toNat : ℕ∞) ≤ (WSumSet n lam A).encard

/-- `λ = (1, 1, 2)`, `A = {0, 1, 2, 3}` in `ℤ/11`: the sums are `{3, 4, 5, 7, 8, 9}` (no wrap-around, the
same six sums as over `ℚ`), `6 < min{11, 7}`. -/
theorem not_conj_MO316143 : ¬ Conj_MO316143 := by
  intro h
  have hS := h 3 (fun k => if k = 3 then 2 else 1) le_rfl (by decide)
    (fun d hd => Nat.dvd_one.1 (hd 1 (by decide))) (by decide) (ZMod 11)
    (by rw [pF_zmod]; exact_mod_cast (by decide : 4 < 11)) ({0, 1, 2, 3} : Finset (ZMod 11)) (by decide)
  have hsub : WSumSet 3 (fun k => if k = 3 then 2 else 1) ({0, 1, 2, 3} : Finset (ZMod 11)) ⊆
      (({3, 4, 5, 7, 8, 9} : Finset (ZMod 11)) : Set (ZMod 11)) :=
    WSumSet_three_sub _ _ _ (by decide)
  have h1 := (hS.trans (Set.encard_le_encard hsub)).trans_eq (Set.encard_coe_eq_coe_finsetCard _)
  rw [pF_zmod] at h1
  revert h1
  decide

end SunLean

#print axioms SunLean.not_conj_MO316143
