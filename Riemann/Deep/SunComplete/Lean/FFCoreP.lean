import Mathlib.FieldTheory.Finite.Basic

/-!
# sunlean · FFCoreP: `q` is a prime power; `2 = 0` for even `q` (Mathlib, `FieldTheory.Finite.Basic` only)
-/

namespace SunLean.FFCore

variable (F : Type) [Field F] [Fintype F]

/-- `q = card F` is a prime power (`FiniteField.isPrimePow_card`). -/
theorem primePow : ∃ p k : ℕ, p.Prime ∧ 0 < k ∧ p ^ k = Fintype.card F :=
  (isPrimePow_nat_iff _).mp (FiniteField.isPrimePow_card F)

/-- For even `q`, `2 = 0` in `F` (`FiniteField.even_card_iff_char_two`). -/
theorem two_eq_zero (h : Fintype.card F % 2 = 0) : (2 : F) = 0 := by
  have hc : ringChar F = 2 := FiniteField.even_card_iff_char_two.mpr h
  have h2 : ((2 : ℕ) : F) = 0 := (ringChar.spec F 2).mpr (by rw [hc])
  exact_mod_cast h2

end SunLean.FFCore
