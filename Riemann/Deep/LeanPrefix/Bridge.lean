import Riemann.Deep.LeanPrefix.Prim0
import Mathlib.GroupTheory.GroupAction.Primitive

/-!
# Bridge to Mathlib's blocks and primitivity

For a subgroup H of Perm (Fin n) acting on Fin n (Mathlib's subgroup action, h • x = h x):
IsBlockOf H B ↔ MulAction.IsBlock H B, and IsPrimitive H ↔ MulAction.IsPreprimitive H (Fin n).
-/

namespace Riemann.Deep.LeanPrefix

open Equiv MulAction
open scoped Pointwise

theorem smul_eq_apply {n : ℕ} {H : Subgroup (Perm (Fin n))} (h : H) (x : Fin n) :
    h • x = (h : Perm (Fin n)) x := rfl

theorem smul_set_eq_image {n : ℕ} {H : Subgroup (Perm (Fin n))} (h : H) (B : Set (Fin n)) :
    h • B = ⇑(h : Perm (Fin n)) '' B := rfl

theorem isBlock_iff_isBlockOf {n : ℕ} (H : Subgroup (Perm (Fin n))) (B : Set (Fin n)) :
    MulAction.IsBlock H B ↔ IsBlockOf H B := by
  rw [MulAction.isBlock_iff_smul_eq_or_disjoint]
  constructor
  · intro h g hg
    have := h ⟨g, hg⟩
    rwa [smul_set_eq_image] at this
  · intro h g
    rw [smul_set_eq_image]
    exact h g g.2

theorem isPretransitive_iff {n : ℕ} (H : Subgroup (Perm (Fin n))) :
    MulAction.IsPretransitive H (Fin n) ↔ ∀ x y : Fin n, ∃ h ∈ H, h x = y := by
  constructor
  · intro ht x y
    obtain ⟨g, hg⟩ := ht.exists_smul_eq x y
    exact ⟨g, g.2, hg⟩
  · intro h
    refine ⟨fun x y => ?_⟩
    obtain ⟨g, hg, hgx⟩ := h x y
    exact ⟨⟨g, hg⟩, hgx⟩

theorem isPreprimitive_iff_isPrimitive {n : ℕ} (H : Subgroup (Perm (Fin n))) :
    MulAction.IsPreprimitive H (Fin n) ↔ IsPrimitive H := by
  constructor
  · intro hp
    refine ⟨(isPretransitive_iff H).mp hp.toIsPretransitive, fun B hB => ?_⟩
    exact hp.isTrivialBlock_of_isBlock ((isBlock_iff_isBlockOf H B).mpr hB)
  · rintro ⟨ht, hb⟩
    have := (isPretransitive_iff H).mpr ht
    exact { isTrivialBlock_of_isBlock := fun {B} hB => hb B ((isBlock_iff_isBlockOf H B).mp hB) }

end Riemann.Deep.LeanPrefix

