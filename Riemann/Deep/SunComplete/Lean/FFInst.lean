import Riemann.Deep.SunComplete.Lean.Sieve
import Mathlib.FieldTheory.Finite.Basic

/-!
# sunlean · FFInst: `FFacts F` for every finite field, from Mathlib — NOT COMPILED

Compiling this file needs `Sieve` and `Mathlib.FieldTheory.Finite.Basic` together: killed at 312 MB
(cap 300 MB). `card_prim` is compiled separately in `FFCoreC.lean` (284 MB); `FFCoreP.lean`
(`primePow`, `two_eq_zero`, `FieldTheory.Finite.Basic` alone) was killed at 309 MB.

`primePow` ← `FiniteField.isPrimePow_card`; `card_prim` ← `IsCyclic.card_orderOf_eq_totient` on
`Fˣ` (cyclic: `RingTheory.IntegralDomain`) with `Fintype.card_units` and `orderOf_units`;
`two_eq_zero` ← `FiniteField.even_card_iff_char_two`, `ringChar.spec`. Also: `phi = Nat.totient`,
`PrimePow ↔ IsPrimePow`, `primeDivisors m = Nat.primeFactors m` (`m ≠ 0`).
-/

namespace SunLean

theorem phi_eq_totient (m : ℕ) : phi m = Nat.totient m := rfl

theorem primePow_iff (q : ℕ) : PrimePow q ↔ IsPrimePow q := (isPrimePow_nat_iff q).symm

theorem primeDivisors_eq (m : ℕ) (hm : m ≠ 0) : primeDivisors m = m.primeFactors := by
  ext p
  simp only [primeDivisors, Finset.mem_filter, Finset.mem_range, Nat.mem_primeFactors]
  constructor
  · rintro ⟨_, hp, hd⟩
    exact ⟨hp, hd, hm⟩
  · rintro ⟨hp, hd, -⟩
    exact ⟨Nat.lt_succ_of_le (Nat.le_of_dvd (Nat.pos_of_ne_zero hm) hd), hp, hd⟩

/-- Every finite field satisfies `FFacts`. -/
theorem ffacts (F : Type) [Field F] [Fintype F] : FFacts F where
  primePow := (primePow_iff _).mpr (FiniteField.isPrimePow_card F)
  card_prim := by
    classical
    have h2 : 2 ≤ Fintype.card F := Fintype.one_lt_card
    have hU : Fintype.card Fˣ = Fintype.card F - 1 := Fintype.card_units
    have h1 := IsCyclic.card_orderOf_eq_totient (α := Fˣ) (d := Fintype.card Fˣ) (dvd_refl _)
    rw [hU] at h1
    rw [phi_eq_totient, ← h1]
    have e : primSet F = (Finset.univ.filter (fun u : Fˣ => orderOf u = Fintype.card F - 1)).map
        ⟨Units.val, Units.val_injective⟩ := by
      ext x
      simp only [mem_primSet, Finset.mem_map, Finset.mem_filter, Finset.mem_univ, true_and,
        Function.Embedding.coeFn_mk]
      constructor
      · intro hx
        have hx0 : x ≠ 0 := by
          intro h0
          rw [h0] at hx
          exact not_isPrim_zero h2 hx
        refine ⟨Units.mk0 x hx0, ?_, rfl⟩
        rw [← orderOf_units]
        exact hx
      · rintro ⟨u, hu, rfl⟩
        unfold IsPrim
        rw [orderOf_units]
        exact hu
    rw [e, Finset.card_map]
  two_eq_zero := fun h => by
    have hc : ringChar F = 2 := FiniteField.even_card_iff_char_two.mpr h
    have h2 : ((2 : ℕ) : F) = 0 := (ringChar.spec F 2).mpr (by rw [hc])
    exact_mod_cast h2

end SunLean
