import Mathlib.RingTheory.IntegralDomain
import Mathlib.GroupTheory.SpecificGroups.Cyclic
import Mathlib.Algebra.GroupWithZero.Units.Fintype

/-!
# sunlean · FFCoreC: `F` has exactly `φ(q − 1)` elements of order `q − 1` (Mathlib, light imports)
-/

namespace SunLean.FFCore

variable (F : Type) [Field F] [Fintype F]

/-- `F` has exactly `φ(q − 1)` elements of order `q − 1` (`Fˣ` cyclic of order `q − 1`). -/
theorem card_prim [DecidableEq F] :
    (Finset.univ.filter (fun x : F => orderOf x = Fintype.card F - 1)).card =
      Nat.totient (Fintype.card F - 1) := by
  have h2 : 2 ≤ Fintype.card F := Fintype.one_lt_card
  have hU : Fintype.card Fˣ = Fintype.card F - 1 := Fintype.card_units F
  have h1 := IsCyclic.card_orderOf_eq_totient (α := Fˣ) (d := Fintype.card Fˣ) (dvd_refl _)
  rw [hU] at h1
  rw [← h1]
  have e : Finset.univ.filter (fun x : F => orderOf x = Fintype.card F - 1) =
      (Finset.univ.filter (fun u : Fˣ => orderOf u = Fintype.card F - 1)).map
        ⟨Units.val, Units.val_injective⟩ := by
    ext x
    simp only [Finset.mem_map, Finset.mem_filter, Finset.mem_univ, true_and,
      Function.Embedding.coeFn_mk]
    constructor
    · intro hx
      have hx0 : x ≠ 0 := by
        intro h0
        rw [h0] at hx
        have h3 := pow_orderOf_eq_one (0 : F)
        rw [hx, zero_pow (by omega)] at h3
        exact zero_ne_one h3
      refine ⟨Units.mk0 x hx0, ?_, rfl⟩
      rw [← orderOf_units]
      exact hx
    · rintro ⟨u, hu, rfl⟩
      rw [orderOf_units]
      exact hu
  rw [e, Finset.card_map]

end SunLean.FFCore

#print axioms SunLean.FFCore.card_prim
