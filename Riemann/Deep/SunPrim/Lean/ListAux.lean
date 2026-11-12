import Mathlib.Algebra.Order.Field.Basic
import Mathlib.Algebra.Order.Ring.Cast
import Mathlib.Data.List.Permutation
import Mathlib.Data.List.Rotate
import Mathlib.Data.Finset.Card
import Mathlib.Data.Int.Cast.Lemmas

/-! List helpers for the refutations over fields (casts `ℤ → R`). -/

namespace SunLean

/-- If `l` is a permutation of `L.map f`, then `l = l₀.map f` for a permutation `l₀` of `L`. -/
theorem exists_of_perm_map {α β : Type} (f : α → β) {l : List β} {L : List α}
    (h : l.Perm (L.map f)) : ∃ l₀ ∈ L.permutations', l = l₀.map f := by
  have hm : l ∈ (L.map f).permutations' := List.mem_permutations'.2 h
  rw [← List.map_permutations'] at hm
  obtain ⟨l₀, h₀, rfl⟩ := List.mem_map.1 hm
  exact ⟨l₀, h₀, rfl⟩

/-- Circular adjacent operations commute with a homomorphism `f`. -/
theorem zipWith_rotate_map {α β : Type} (f : α → β) (g : β → β → β) (g₀ : α → α → α)
    (hg : ∀ x y, g (f x) (f y) = f (g₀ x y)) (l : List α) :
    List.zipWith g (l.map f) ((l.map f).rotate 1) = (List.zipWith g₀ l (l.rotate 1)).map f := by
  rw [← List.map_rotate, List.zipWith_map, List.map_zipWith]
  simp only [hg]

end SunLean
