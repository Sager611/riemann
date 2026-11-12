import Riemann.Deep.SunPrim.Lean.FieldsBase

/-! Sun, arXiv 1811.10503, Conjecture 4.9: refuted in `ℤ/11` (`n = 3, 4`) and `ℤ/13` (`n = 4`). -/

namespace SunLean

open Finset

/-- Z.-W. Sun, arXiv 1811.10503 v6, Conjecture 4.9 (TeX l. 767-770): `n ≥ 1`, `F` a field with
`p(F) > n + 1`, `A ⊆ F` finite with `|A| ≥ n + δ_{n,3}`; then
`S(A) = {∑_{k=1}^n k a_k : a_1, …, a_n distinct elements of A}` has
`|S(A)| ≥ min{p(F), (|A| - n) n(n+1)/2 + n(n²-1)/6 + 1}` (all divisions exact, `|A| ≥ n`). -/
def Conj_4_9 : Prop :=
  ∀ n : ℕ, 0 < n → ∀ (F : Type) [Field F], ((n + 1 : ℕ) : ℕ∞) < pF F →
    ∀ A : Finset F, n + (if n = 3 then 1 else 0) ≤ A.card →
      min (pF F) (((A.card - n) * (n * (n + 1) / 2) + n * (n ^ 2 - 1) / 6 + 1 : ℕ) : ℕ∞) ≤
        (WSumSet n (fun k => k) A).encard

/-- `n = 3`, `F = ℤ/11`, `A = {0, 1, 2, 5}`: `S(A) = {0, 1, 2, 4, 5, 6, 7, 8, 9}`, `9 < min{11, 11}`. -/
theorem not_conj_4_9_n3 : ¬ Conj_4_9 := by
  intro h
  have hS := h 3 (by decide) (ZMod 11) (by rw [pF_zmod]; exact_mod_cast (by decide : 4 < 11))
    ({0, 1, 2, 5} : Finset (ZMod 11)) (by decide)
  have hsub : WSumSet 3 (fun k => k) ({0, 1, 2, 5} : Finset (ZMod 11)) ⊆
      (({0, 1, 2, 4, 5, 6, 7, 8, 9} : Finset (ZMod 11)) : Set (ZMod 11)) :=
    WSumSet_three_sub _ _ _ (by decide)
  have h1 := (hS.trans (Set.encard_le_encard hsub)).trans_eq (Set.encard_coe_eq_coe_finsetCard _)
  rw [pF_zmod] at h1
  revert h1
  decide

theorem not_conj_4_9 : ¬ Conj_4_9 := not_conj_4_9_n3

end SunLean

#print axioms SunLean.not_conj_4_9_n3
