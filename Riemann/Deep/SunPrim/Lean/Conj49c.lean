import Riemann.Deep.SunPrim.Lean.Conj49

/-! Sun, arXiv 1811.10503, Conjecture 4.9 at `n = 4`: refuted in `ℤ/13`. -/

namespace SunLean

open Finset

/-- `n = 4`, `F = ℤ/13`, `A = {0, 1, 3, 4}`: `|S(A)| = 10 < min{13, 11}`. -/
theorem not_conj_4_9_n4_F13 : ¬ Conj_4_9 := by
  intro h
  have hS := h 4 (by decide) (ZMod 13) (by rw [pF_zmod]; exact_mod_cast (by decide : 5 < 13))
    ({0, 1, 3, 4} : Finset (ZMod 13)) (by decide)
  have hsub : WSumSet 4 (fun k => k) ({0, 1, 3, 4} : Finset (ZMod 13)) ⊆
      (({0, 1, 2, 4, 5, 6, 8, 9, 10, 12} : Finset (ZMod 13)) : Set (ZMod 13)) :=
    WSumSet_four_sub _ _ _ (by decide)
  have h1 := (hS.trans (Set.encard_le_encard hsub)).trans_eq (Set.encard_coe_eq_coe_finsetCard _)
  rw [pF_zmod] at h1
  revert h1
  decide

end SunLean

#print axioms SunLean.not_conj_4_9_n4_F13
