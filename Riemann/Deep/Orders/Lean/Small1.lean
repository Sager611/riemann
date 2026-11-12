import Riemann.Deep.Orders.Lean.Search

/-! Exact sets of perfect N below N0 (working orders) or below 6 (failing orders):
witnesses by the checker, non-existence by the sound search noPerfect. -/
namespace Orders

theorem small_pxmd : ∀ N, N < 9 → ((∃ ms, Perfect (⟨.add, .mul, .sub, .div⟩ : Order) N ms) ↔ N ∈ []) := by
  intro N hN
  rcases (by omega : N = 0 ∨ N = 1 ∨ N = 2 ∨ N = 3 ∨ N = 4 ∨ N = 5 ∨ N = 6 ∨ N = 7 ∨ N = 8) with rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl
  · exact ⟨fun ⟨_, h⟩ => absurd h.1 (by omega), fun h => absurd h (by decide)⟩
  · exact ⟨fun ⟨ms, h⟩ => absurd h (noPerfect_sound (by decide +kernel) ms), fun h => absurd h (by decide)⟩
  · exact ⟨fun ⟨ms, h⟩ => absurd h (noPerfect_sound (by decide +kernel) ms), fun h => absurd h (by decide)⟩
  · exact ⟨fun ⟨ms, h⟩ => absurd h (noPerfect_sound (by decide +kernel) ms), fun h => absurd h (by decide)⟩
  · exact ⟨fun ⟨ms, h⟩ => absurd h (noPerfect_sound (by decide +kernel) ms), fun h => absurd h (by decide)⟩
  · exact ⟨fun ⟨ms, h⟩ => absurd h (noPerfect_sound (by decide +kernel) ms), fun h => absurd h (by decide)⟩
  · exact ⟨fun ⟨ms, h⟩ => absurd h (noPerfect_sound (by decide +kernel) ms), fun h => absurd h (by decide)⟩
  · exact ⟨fun ⟨ms, h⟩ => absurd h (noPerfect_sound (by decide +kernel) ms), fun h => absurd h (by decide)⟩
  · exact ⟨fun ⟨ms, h⟩ => absurd h (noPerfect_sound (by decide +kernel) ms), fun h => absurd h (by decide)⟩

theorem small_pxdm : ∀ N, N < 10 → ((∃ ms, Perfect (⟨.add, .mul, .div, .sub⟩ : Order) N ms) ↔ N ∈ []) := by
  intro N hN
  rcases (by omega : N = 0 ∨ N = 1 ∨ N = 2 ∨ N = 3 ∨ N = 4 ∨ N = 5 ∨ N = 6 ∨ N = 7 ∨ N = 8 ∨ N = 9) with rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl
  · exact ⟨fun ⟨_, h⟩ => absurd h.1 (by omega), fun h => absurd h (by decide)⟩
  · exact ⟨fun ⟨ms, h⟩ => absurd h (noPerfect_sound (by decide +kernel) ms), fun h => absurd h (by decide)⟩
  · exact ⟨fun ⟨ms, h⟩ => absurd h (noPerfect_sound (by decide +kernel) ms), fun h => absurd h (by decide)⟩
  · exact ⟨fun ⟨ms, h⟩ => absurd h (noPerfect_sound (by decide +kernel) ms), fun h => absurd h (by decide)⟩
  · exact ⟨fun ⟨ms, h⟩ => absurd h (noPerfect_sound (by decide +kernel) ms), fun h => absurd h (by decide)⟩
  · exact ⟨fun ⟨ms, h⟩ => absurd h (noPerfect_sound (by decide +kernel) ms), fun h => absurd h (by decide)⟩
  · exact ⟨fun ⟨ms, h⟩ => absurd h (noPerfect_sound (by decide +kernel) ms), fun h => absurd h (by decide)⟩
  · exact ⟨fun ⟨ms, h⟩ => absurd h (noPerfect_sound (by decide +kernel) ms), fun h => absurd h (by decide)⟩
  · exact ⟨fun ⟨ms, h⟩ => absurd h (noPerfect_sound (by decide +kernel) ms), fun h => absurd h (by decide)⟩
  · exact ⟨fun ⟨ms, h⟩ => absurd h (noPerfect_sound (by decide +kernel) ms), fun h => absurd h (by decide)⟩

theorem small_pmxd : ∀ N, N < 6 → ((∃ ms, Perfect (⟨.add, .sub, .mul, .div⟩ : Order) N ms) ↔ N ∈ []) := by
  intro N hN
  rcases (by omega : N = 0 ∨ N = 1 ∨ N = 2 ∨ N = 3 ∨ N = 4 ∨ N = 5) with rfl | rfl | rfl | rfl | rfl | rfl
  · exact ⟨fun ⟨_, h⟩ => absurd h.1 (by omega), fun h => absurd h (by decide)⟩
  · exact ⟨fun ⟨ms, h⟩ => absurd h (noPerfect_sound (by decide +kernel) ms), fun h => absurd h (by decide)⟩
  · exact ⟨fun ⟨ms, h⟩ => absurd h (noPerfect_sound (by decide +kernel) ms), fun h => absurd h (by decide)⟩
  · exact ⟨fun ⟨ms, h⟩ => absurd h (noPerfect_sound (by decide +kernel) ms), fun h => absurd h (by decide)⟩
  · exact ⟨fun ⟨ms, h⟩ => absurd h (noPerfect_sound (by decide +kernel) ms), fun h => absurd h (by decide)⟩
  · exact ⟨fun ⟨ms, h⟩ => absurd h (noPerfect_sound (by decide +kernel) ms), fun h => absurd h (by decide)⟩

theorem small_pmdx : ∀ N, N < 6 → ((∃ ms, Perfect (⟨.add, .sub, .div, .mul⟩ : Order) N ms) ↔ N ∈ []) := by
  intro N hN
  rcases (by omega : N = 0 ∨ N = 1 ∨ N = 2 ∨ N = 3 ∨ N = 4 ∨ N = 5) with rfl | rfl | rfl | rfl | rfl | rfl
  · exact ⟨fun ⟨_, h⟩ => absurd h.1 (by omega), fun h => absurd h (by decide)⟩
  · exact ⟨fun ⟨ms, h⟩ => absurd h (noPerfect_sound (by decide +kernel) ms), fun h => absurd h (by decide)⟩
  · exact ⟨fun ⟨ms, h⟩ => absurd h (noPerfect_sound (by decide +kernel) ms), fun h => absurd h (by decide)⟩
  · exact ⟨fun ⟨ms, h⟩ => absurd h (noPerfect_sound (by decide +kernel) ms), fun h => absurd h (by decide)⟩
  · exact ⟨fun ⟨ms, h⟩ => absurd h (noPerfect_sound (by decide +kernel) ms), fun h => absurd h (by decide)⟩
  · exact ⟨fun ⟨ms, h⟩ => absurd h (noPerfect_sound (by decide +kernel) ms), fun h => absurd h (by decide)⟩

end Orders
