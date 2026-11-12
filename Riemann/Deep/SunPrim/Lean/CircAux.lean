import Mathlib.Data.List.Rotate
import Mathlib.Data.Fintype.Card

/-! List lemmas for circular orders: rotation, duplicate-free prefixes, covering lists. -/

namespace SunLean

theorem nodup_rev_take {α : Type} {l : List α} (h : l.Nodup) (k : ℕ) : (l.take k).reverse.Nodup :=
  List.nodup_reverse.2 (h.sublist (List.take_sublist k l))

/-- A list of length `card α` containing every element is duplicate-free. -/
theorem nodup_of_cover {α : Type} [Fintype α] [DecidableEq α] (s : List α)
    (hlen : s.length = Fintype.card α) (hcov : ∀ y, y ∈ s) : s.Nodup := by
  have h1 : s.toFinset = Finset.univ :=
    Finset.eq_univ_iff_forall.2 (fun y => List.mem_toFinset.2 (hcov y))
  have h2 : s.dedup.length = s.length := by
    rw [hlen, ← Finset.card_univ, ← h1, List.card_toFinset]
  have h3 : s.dedup = s := (List.dedup_sublist s).eq_of_length h2
  rw [← h3]
  exact List.nodup_dedup s

/-- Rotating a circular order rotates its list of adjacent values. -/
theorem sums_rotate {α β : Type} (f : α → α → β) (l : List α) (k : ℕ) :
    List.zipWith f (l.rotate k) ((l.rotate k).rotate 1) = (List.zipWith f l (l.rotate 1)).rotate k := by
  rw [List.zipWith_rotate_distrib f l (l.rotate 1) k (List.length_rotate l 1).symm,
    List.rotate_rotate, List.rotate_rotate, Nat.add_comm]

theorem exists_rotate_head {α : Type} [DecidableEq α] {l : List α} {x : α} (hx : x ∈ l) :
    ∃ k, (l.rotate k).head? = some x := by
  refine ⟨l.idxOf x, ?_⟩
  rw [List.head?_rotate (List.idxOf_lt_length_of_mem hx), List.getElem?_idxOf hx]

theorem eq_of_length_eight {α : Type} (l : List α) (h : l.length = 8) :
    ∃ x0 x1 x2 x3 x4 x5 x6 x7, l = [x0, x1, x2, x3, x4, x5, x6, x7] := by
  match l, h with
  | [x0, x1, x2, x3, x4, x5, x6, x7], _ => exact ⟨x0, x1, x2, x3, x4, x5, x6, x7, rfl⟩

theorem eq_of_length_nine {α : Type} (l : List α) (h : l.length = 9) :
    ∃ x0 x1 x2 x3 x4 x5 x6 x7 x8, l = [x0, x1, x2, x3, x4, x5, x6, x7, x8] := by
  match l, h with
  | [x0, x1, x2, x3, x4, x5, x6, x7, x8], _ => exact ⟨x0, x1, x2, x3, x4, x5, x6, x7, x8, rfl⟩

end SunLean
