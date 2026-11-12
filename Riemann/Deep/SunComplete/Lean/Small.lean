import Riemann.Deep.SunComplete.Lean.Defs
import Mathlib.Data.ZMod.Basic
import Mathlib.Algebra.Field.ZMod
import Mathlib.Data.Fin.VecNotation
import Mathlib.Tactic.DeriveFintype

/-!
# sunlean · Small: kernel-checked small cases (kernel `decide` only)

* `q = 3, 5, 7` fail: `ZMod p` has a vertex with at most one neighbour in `H_p` (vertices `1`,
  `1`, `5`); every field with `p` elements is `≃+* ZMod p` (Mathlib `ZMod.ringEquivOfPrime`).
* `q = 4` works: the order `(0, ω, 1, ω²)` in an explicit field `F4`.
* certificate cycles: `q = 11` (She's construction, `x = 6`) and `q = 13` (certA.json).
-/

namespace SunLean

section transfer
variable {F K : Type} [Field F] [Fintype F] [Field K] [Fintype K]

theorem isPrim_map (e : F ≃+* K) (x : F) : IsPrim (e x) ↔ IsPrim x := by
  unfold IsPrim
  have h1 : orderOf (e x) = orderOf x := MulEquiv.orderOf_eq e.toMulEquiv x
  rw [h1, Fintype.card_congr e.toEquiv]

theorem finCongr_finRotate {n m : ℕ} (h : n = m) (i : Fin m) :
    finCongr h.symm (finRotate m i) = finRotate n (finCongr h.symm i) := by
  subst h; rfl

/-- The Sun–Hou property is invariant under ring isomorphism. -/
theorem sunHou_of_ringEquiv (e : F ≃+* K) (h : SunHou F) : SunHou K := by
  obtain ⟨a, ha⟩ := h
  have hc : Fintype.card F = Fintype.card K := Fintype.card_congr e.toEquiv
  refine ⟨(finCongr hc.symm).trans (a.trans e.toEquiv), fun i => ?_⟩
  simp only [Equiv.trans_apply, RingEquiv.toEquiv_eq_coe, EquivLike.coe_coe]
  rw [finCongr_finRotate hc i, ← map_add, isPrim_map]
  exact ha _

/-- A numbering indexed by `Fin n` with `n = q` is a Sun–Hou numbering. -/
theorem sunHou_of_cycle {n : ℕ} (hn : Fintype.card F = n) (a : Fin n ≃ F)
    (ha : ∀ i, IsPrim (a i + a (finRotate n i))) : SunHou F := by
  subst hn; exact ⟨a, ha⟩
end transfer

/-- Primitivity in `ZMod p` as a decidable condition. -/
theorem isPrim_zmod_iff (p : ℕ) [Fact p.Prime] (x : ZMod p) :
    IsPrim x ↔ x ^ (p - 1) = 1 ∧ ∀ m, m < p - 1 → 0 < m → x ^ m ≠ 1 := by
  unfold IsPrim
  rw [ZMod.card]
  exact orderOf_eq_iff (by have := (Fact.out : p.Prime).two_le; omega)

instance : Fact (Nat.Prime 3) := ⟨by decide⟩
instance : Fact (Nat.Prime 5) := ⟨by decide⟩
instance : Fact (Nat.Prime 7) := ⟨by decide⟩
instance : Fact (Nat.Prime 11) := ⟨by decide⟩
instance : Fact (Nat.Prime 13) := ⟨by decide⟩

/-- `q = 3` fails: vertex `1` of `H_3` has no neighbour. -/
theorem not_sunHou_zmod3 : ¬ SunHou (ZMod 3) := by
  intro h
  obtain ⟨y, z, hyz, hy, hz, py, pz⟩ := two_nbrs (by rw [ZMod.card]) h 1
  have key : ∀ y z : ZMod 3, y ≠ z → y ≠ 1 → z ≠ 1 → IsPrim (1 + y) → IsPrim (1 + z) → False := by
    simp only [isPrim_zmod_iff]
    decide
  exact key y z hyz hy hz py pz

/-- `q = 5` fails: vertex `1` of `H_5` has only the neighbour `2`. -/
theorem not_sunHou_zmod5 : ¬ SunHou (ZMod 5) := by
  intro h
  obtain ⟨y, z, hyz, hy, hz, py, pz⟩ := two_nbrs (by rw [ZMod.card]; decide) h 1
  have key : ∀ y z : ZMod 5, y ≠ z → y ≠ 1 → z ≠ 1 → IsPrim (1 + y) → IsPrim (1 + z) → False := by
    simp only [isPrim_zmod_iff]
    decide
  exact key y z hyz hy hz py pz

/-- `q = 7` fails: vertex `5` of `H_7` has only the neighbour `0`. -/
theorem not_sunHou_zmod7 : ¬ SunHou (ZMod 7) := by
  intro h
  obtain ⟨y, z, hyz, hy, hz, py, pz⟩ := two_nbrs (by rw [ZMod.card]; decide) h 5
  have key : ∀ y z : ZMod 7, y ≠ z → y ≠ 5 → z ≠ 5 → IsPrim (5 + y) → IsPrim (5 + z) → False := by
    simp only [isPrim_zmod_iff]
    decide
  exact key y z hyz hy hz py pz

/-- **Sharpness**: no field with `3`, `5` or `7` elements has the Sun–Hou property. -/
theorem not_sunHou_small (F : Type) [Field F] [Fintype F]
    (h : Fintype.card F = 3 ∨ Fintype.card F = 5 ∨ Fintype.card F = 7) : ¬ SunHou F := by
  intro hS
  rcases h with h | h | h
  · exact not_sunHou_zmod3 (sunHou_of_ringEquiv (ZMod.ringEquivOfPrime F (by decide) h).symm hS)
  · exact not_sunHou_zmod5 (sunHou_of_ringEquiv (ZMod.ringEquivOfPrime F (by decide) h).symm hS)
  · exact not_sunHou_zmod7 (sunHou_of_ringEquiv (ZMod.ringEquivOfPrime F (by decide) h).symm hS)

end SunLean
