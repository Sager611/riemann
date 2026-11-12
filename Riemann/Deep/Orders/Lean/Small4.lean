import Riemann.Deep.Orders.Lean.Search

/-! Exact sets of perfect N below N0 (working orders) or below 6 (failing orders):
witnesses by the checker, non-existence by the sound search noPerfect. -/
namespace Orders

theorem small_mpxd : ∀ N, N < 6 → ((∃ ms, Perfect (⟨.sub, .add, .mul, .div⟩ : Order) N ms) ↔ N ∈ []) := by
  intro N hN
  rcases (by omega : N = 0 ∨ N = 1 ∨ N = 2 ∨ N = 3 ∨ N = 4 ∨ N = 5) with rfl | rfl | rfl | rfl | rfl | rfl
  · exact ⟨fun ⟨_, h⟩ => absurd h.1 (by omega), fun h => absurd h (by decide)⟩
  · exact ⟨fun ⟨ms, h⟩ => absurd h (noPerfect_sound (by decide +kernel) ms), fun h => absurd h (by decide)⟩
  · exact ⟨fun ⟨ms, h⟩ => absurd h (noPerfect_sound (by decide +kernel) ms), fun h => absurd h (by decide)⟩
  · exact ⟨fun ⟨ms, h⟩ => absurd h (noPerfect_sound (by decide +kernel) ms), fun h => absurd h (by decide)⟩
  · exact ⟨fun ⟨ms, h⟩ => absurd h (noPerfect_sound (by decide +kernel) ms), fun h => absurd h (by decide)⟩
  · exact ⟨fun ⟨ms, h⟩ => absurd h (noPerfect_sound (by decide +kernel) ms), fun h => absurd h (by decide)⟩

theorem small_mpdx : ∀ N, N < 6 → ((∃ ms, Perfect (⟨.sub, .add, .div, .mul⟩ : Order) N ms) ↔ N ∈ []) := by
  intro N hN
  rcases (by omega : N = 0 ∨ N = 1 ∨ N = 2 ∨ N = 3 ∨ N = 4 ∨ N = 5) with rfl | rfl | rfl | rfl | rfl | rfl
  · exact ⟨fun ⟨_, h⟩ => absurd h.1 (by omega), fun h => absurd h (by decide)⟩
  · exact ⟨fun ⟨ms, h⟩ => absurd h (noPerfect_sound (by decide +kernel) ms), fun h => absurd h (by decide)⟩
  · exact ⟨fun ⟨ms, h⟩ => absurd h (noPerfect_sound (by decide +kernel) ms), fun h => absurd h (by decide)⟩
  · exact ⟨fun ⟨ms, h⟩ => absurd h (noPerfect_sound (by decide +kernel) ms), fun h => absurd h (by decide)⟩
  · exact ⟨fun ⟨ms, h⟩ => absurd h (noPerfect_sound (by decide +kernel) ms), fun h => absurd h (by decide)⟩
  · exact ⟨fun ⟨ms, h⟩ => absurd h (noPerfect_sound (by decide +kernel) ms), fun h => absurd h (by decide)⟩

theorem small_mxpd : ∀ N, N < 6 → ((∃ ms, Perfect (⟨.sub, .mul, .add, .div⟩ : Order) N ms) ↔ N ∈ []) := by
  intro N hN
  rcases (by omega : N = 0 ∨ N = 1 ∨ N = 2 ∨ N = 3 ∨ N = 4 ∨ N = 5) with rfl | rfl | rfl | rfl | rfl | rfl
  · exact ⟨fun ⟨_, h⟩ => absurd h.1 (by omega), fun h => absurd h (by decide)⟩
  · exact ⟨fun ⟨ms, h⟩ => absurd h (noPerfect_sound (by decide +kernel) ms), fun h => absurd h (by decide)⟩
  · exact ⟨fun ⟨ms, h⟩ => absurd h (noPerfect_sound (by decide +kernel) ms), fun h => absurd h (by decide)⟩
  · exact ⟨fun ⟨ms, h⟩ => absurd h (noPerfect_sound (by decide +kernel) ms), fun h => absurd h (by decide)⟩
  · exact ⟨fun ⟨ms, h⟩ => absurd h (noPerfect_sound (by decide +kernel) ms), fun h => absurd h (by decide)⟩
  · exact ⟨fun ⟨ms, h⟩ => absurd h (noPerfect_sound (by decide +kernel) ms), fun h => absurd h (by decide)⟩

theorem small_mxdp : ∀ N, N < 6 → ((∃ ms, Perfect (⟨.sub, .mul, .div, .add⟩ : Order) N ms) ↔ N ∈ []) := by
  intro N hN
  rcases (by omega : N = 0 ∨ N = 1 ∨ N = 2 ∨ N = 3 ∨ N = 4 ∨ N = 5) with rfl | rfl | rfl | rfl | rfl | rfl
  · exact ⟨fun ⟨_, h⟩ => absurd h.1 (by omega), fun h => absurd h (by decide)⟩
  · exact ⟨fun ⟨ms, h⟩ => absurd h (noPerfect_sound (by decide +kernel) ms), fun h => absurd h (by decide)⟩
  · exact ⟨fun ⟨ms, h⟩ => absurd h (noPerfect_sound (by decide +kernel) ms), fun h => absurd h (by decide)⟩
  · exact ⟨fun ⟨ms, h⟩ => absurd h (noPerfect_sound (by decide +kernel) ms), fun h => absurd h (by decide)⟩
  · exact ⟨fun ⟨ms, h⟩ => absurd h (noPerfect_sound (by decide +kernel) ms), fun h => absurd h (by decide)⟩
  · exact ⟨fun ⟨ms, h⟩ => absurd h (noPerfect_sound (by decide +kernel) ms), fun h => absurd h (by decide)⟩

end Orders
