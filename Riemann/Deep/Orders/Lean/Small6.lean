import Riemann.Deep.Orders.Lean.Search

/-! Exact sets of perfect N below N0 (working orders) or below 6 (failing orders):
witnesses by the checker, non-existence by the sound search noPerfect. -/
namespace Orders

theorem small_dxpm : ∀ N, N < 12 → ((∃ ms, Perfect (⟨.div, .mul, .add, .sub⟩ : Order) N ms) ↔ N ∈ [2, 4, 5, 6, 8, 9, 10]) := by
  intro N hN
  rcases (by omega : N = 0 ∨ N = 1 ∨ N = 2 ∨ N = 3 ∨ N = 4 ∨ N = 5 ∨ N = 6 ∨ N = 7 ∨ N = 8 ∨ N = 9 ∨ N = 10 ∨ N = 11) with rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl
  · exact ⟨fun ⟨_, h⟩ => absurd h.1 (by omega), fun h => absurd h (by decide)⟩
  · exact ⟨fun ⟨ms, h⟩ => absurd h (noPerfect_sound (by decide +kernel) ms), fun h => absurd h (by decide)⟩
  · exact ⟨fun _ => by decide, fun _ => ⟨[⟨2, 2, 1⟩], (check_iff _ _ _).1 (by decide +kernel)⟩⟩
  · exact ⟨fun ⟨ms, h⟩ => absurd h (noPerfect_sound (by decide +kernel) ms), fun h => absurd h (by decide)⟩
  · exact ⟨fun _ => by decide, fun _ => ⟨[⟨2, 2, 1⟩, ⟨2, 2, 4⟩, ⟨2, 1, 3⟩], (check_iff _ _ _).1 (by decide +kernel)⟩⟩
  · exact ⟨fun _ => by decide, fun _ => ⟨[⟨2, 2, 1⟩, ⟨2, 2, 4⟩, ⟨1, 4, 5⟩, ⟨4, 1, 3⟩], (check_iff _ _ _).1 (by decide +kernel)⟩⟩
  · exact ⟨fun _ => by decide, fun _ => ⟨[⟨2, 2, 1⟩, ⟨2, 2, 4⟩, ⟨2, 4, 6⟩, ⟨6, 1, 5⟩, ⟨6, 2, 3⟩], (check_iff _ _ _).1 (by decide +kernel)⟩⟩
  · exact ⟨fun ⟨ms, h⟩ => absurd h (noPerfect_sound (by decide +kernel) ms), fun h => absurd h (by decide)⟩
  · exact ⟨fun _ => by decide, fun _ => ⟨[⟨2, 2, 1⟩, ⟨2, 2, 4⟩, ⟨2, 4, 6⟩, ⟨6, 1, 5⟩, ⟨6, 2, 3⟩, ⟨2, 4, 8⟩, ⟨2, 5, 7⟩], (check_iff _ _ _).1 (by decide +kernel)⟩⟩
  · exact ⟨fun _ => by decide, fun _ => ⟨[⟨2, 2, 1⟩, ⟨2, 2, 4⟩, ⟨4, 4, 8⟩, ⟨8, 2, 6⟩, ⟨6, 2, 3⟩, ⟨3, 3, 9⟩, ⟨1, 6, 7⟩, ⟨8, 3, 5⟩], (check_iff _ _ _).1 (by decide +kernel)⟩⟩
  · exact ⟨fun _ => by decide, fun _ => ⟨[⟨2, 2, 1⟩, ⟨2, 2, 4⟩, ⟨4, 4, 8⟩, ⟨8, 2, 6⟩, ⟨6, 2, 3⟩, ⟨3, 3, 9⟩, ⟨2, 8, 10⟩, ⟨8, 1, 7⟩, ⟨10, 2, 5⟩], (check_iff _ _ _).1 (by decide +kernel)⟩⟩
  · exact ⟨fun ⟨ms, h⟩ => absurd h (noPerfect_sound (by decide +kernel) ms), fun h => absurd h (by decide)⟩

theorem small_dxmp : ∀ N, N < 6 → ((∃ ms, Perfect (⟨.div, .mul, .sub, .add⟩ : Order) N ms) ↔ N ∈ [2, 4, 5]) := by
  intro N hN
  rcases (by omega : N = 0 ∨ N = 1 ∨ N = 2 ∨ N = 3 ∨ N = 4 ∨ N = 5) with rfl | rfl | rfl | rfl | rfl | rfl
  · exact ⟨fun ⟨_, h⟩ => absurd h.1 (by omega), fun h => absurd h (by decide)⟩
  · exact ⟨fun ⟨ms, h⟩ => absurd h (noPerfect_sound (by decide +kernel) ms), fun h => absurd h (by decide)⟩
  · exact ⟨fun _ => by decide, fun _ => ⟨[⟨2, 2, 1⟩], (check_iff _ _ _).1 (by decide +kernel)⟩⟩
  · exact ⟨fun ⟨ms, h⟩ => absurd h (noPerfect_sound (by decide +kernel) ms), fun h => absurd h (by decide)⟩
  · exact ⟨fun _ => by decide, fun _ => ⟨[⟨2, 2, 1⟩, ⟨2, 2, 4⟩, ⟨4, 1, 3⟩], (check_iff _ _ _).1 (by decide +kernel)⟩⟩
  · exact ⟨fun _ => by decide, fun _ => ⟨[⟨2, 2, 1⟩, ⟨2, 2, 4⟩, ⟨4, 1, 3⟩, ⟨2, 3, 5⟩], (check_iff _ _ _).1 (by decide +kernel)⟩⟩

theorem small_dmpx : ∀ N, N < 6 → ((∃ ms, Perfect (⟨.div, .sub, .add, .mul⟩ : Order) N ms) ↔ N ∈ [2]) := by
  intro N hN
  rcases (by omega : N = 0 ∨ N = 1 ∨ N = 2 ∨ N = 3 ∨ N = 4 ∨ N = 5) with rfl | rfl | rfl | rfl | rfl | rfl
  · exact ⟨fun ⟨_, h⟩ => absurd h.1 (by omega), fun h => absurd h (by decide)⟩
  · exact ⟨fun ⟨ms, h⟩ => absurd h (noPerfect_sound (by decide +kernel) ms), fun h => absurd h (by decide)⟩
  · exact ⟨fun _ => by decide, fun _ => ⟨[⟨2, 2, 1⟩], (check_iff _ _ _).1 (by decide +kernel)⟩⟩
  · exact ⟨fun ⟨ms, h⟩ => absurd h (noPerfect_sound (by decide +kernel) ms), fun h => absurd h (by decide)⟩
  · exact ⟨fun ⟨ms, h⟩ => absurd h (noPerfect_sound (by decide +kernel) ms), fun h => absurd h (by decide)⟩
  · exact ⟨fun ⟨ms, h⟩ => absurd h (noPerfect_sound (by decide +kernel) ms), fun h => absurd h (by decide)⟩

theorem small_dmxp : ∀ N, N < 6 → ((∃ ms, Perfect (⟨.div, .sub, .mul, .add⟩ : Order) N ms) ↔ N ∈ [2]) := by
  intro N hN
  rcases (by omega : N = 0 ∨ N = 1 ∨ N = 2 ∨ N = 3 ∨ N = 4 ∨ N = 5) with rfl | rfl | rfl | rfl | rfl | rfl
  · exact ⟨fun ⟨_, h⟩ => absurd h.1 (by omega), fun h => absurd h (by decide)⟩
  · exact ⟨fun ⟨ms, h⟩ => absurd h (noPerfect_sound (by decide +kernel) ms), fun h => absurd h (by decide)⟩
  · exact ⟨fun _ => by decide, fun _ => ⟨[⟨2, 2, 1⟩], (check_iff _ _ _).1 (by decide +kernel)⟩⟩
  · exact ⟨fun ⟨ms, h⟩ => absurd h (noPerfect_sound (by decide +kernel) ms), fun h => absurd h (by decide)⟩
  · exact ⟨fun ⟨ms, h⟩ => absurd h (noPerfect_sound (by decide +kernel) ms), fun h => absurd h (by decide)⟩
  · exact ⟨fun ⟨ms, h⟩ => absurd h (noPerfect_sound (by decide +kernel) ms), fun h => absurd h (by decide)⟩

end Orders
