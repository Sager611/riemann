import Riemann.Deep.SunPrim.Lean.Conj49

/-! Sun, arXiv 1811.10503, Conjecture 4.9 at `n = 4`: refuted in `ℤ/11`. -/

namespace SunLean

open Finset

/-- `n = 4`, `F = ℤ/11`, `A = {0, 1, 2, 5}`: `|S(A)| = 9 < min{11, 11}`. -/
theorem not_conj_4_9_n4_F11 : ¬ Conj_4_9 := by
  intro h
  have hS := h 4 (by decide) (ZMod 11) (by rw [pF_zmod]; exact_mod_cast (by decide : 5 < 11))
    ({0, 1, 2, 5} : Finset (ZMod 11)) (by decide)
  have hsub : WSumSet 4 (fun k => k) ({0, 1, 2, 5} : Finset (ZMod 11)) ⊆
      (({1, 2, 3, 4, 5, 6, 8, 9, 10} : Finset (ZMod 11)) : Set (ZMod 11)) :=
    WSumSet_four_sub _ _ _ (by decide)
  have h1 := (hS.trans (Set.encard_le_encard hsub)).trans_eq (Set.encard_coe_eq_coe_finsetCard _)
  rw [pF_zmod] at h1
  revert h1
  decide

end SunLean

#print axioms SunLean.not_conj_4_9_n4_F11
