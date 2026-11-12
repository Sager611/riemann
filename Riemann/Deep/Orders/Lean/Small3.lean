import Riemann.Deep.Orders.Lean.Search

/-! Exact sets of perfect N below N0 (working orders) or below 6 (failing orders):
witnesses by the checker, non-existence by the sound search noPerfect. -/
namespace Orders

theorem small_xmpd : ∀ N, N < 6 → ((∃ ms, Perfect (⟨.mul, .sub, .add, .div⟩ : Order) N ms) ↔ N ∈ []) := by
  intro N hN
  rcases (by omega : N = 0 ∨ N = 1 ∨ N = 2 ∨ N = 3 ∨ N = 4 ∨ N = 5) with rfl | rfl | rfl | rfl | rfl | rfl
  · exact ⟨fun ⟨_, h⟩ => absurd h.1 (by omega), fun h => absurd h (by decide)⟩
  · exact ⟨fun ⟨ms, h⟩ => absurd h (noPerfect_sound (by decide +kernel) ms), fun h => absurd h (by decide)⟩
  · exact ⟨fun ⟨ms, h⟩ => absurd h (noPerfect_sound (by decide +kernel) ms), fun h => absurd h (by decide)⟩
  · exact ⟨fun ⟨ms, h⟩ => absurd h (noPerfect_sound (by decide +kernel) ms), fun h => absurd h (by decide)⟩
  · exact ⟨fun ⟨ms, h⟩ => absurd h (noPerfect_sound (by decide +kernel) ms), fun h => absurd h (by decide)⟩
  · exact ⟨fun ⟨ms, h⟩ => absurd h (noPerfect_sound (by decide +kernel) ms), fun h => absurd h (by decide)⟩

theorem small_xmdp : ∀ N, N < 6 → ((∃ ms, Perfect (⟨.mul, .sub, .div, .add⟩ : Order) N ms) ↔ N ∈ []) := by
  intro N hN
  rcases (by omega : N = 0 ∨ N = 1 ∨ N = 2 ∨ N = 3 ∨ N = 4 ∨ N = 5) with rfl | rfl | rfl | rfl | rfl | rfl
  · exact ⟨fun ⟨_, h⟩ => absurd h.1 (by omega), fun h => absurd h (by decide)⟩
  · exact ⟨fun ⟨ms, h⟩ => absurd h (noPerfect_sound (by decide +kernel) ms), fun h => absurd h (by decide)⟩
  · exact ⟨fun ⟨ms, h⟩ => absurd h (noPerfect_sound (by decide +kernel) ms), fun h => absurd h (by decide)⟩
  · exact ⟨fun ⟨ms, h⟩ => absurd h (noPerfect_sound (by decide +kernel) ms), fun h => absurd h (by decide)⟩
  · exact ⟨fun ⟨ms, h⟩ => absurd h (noPerfect_sound (by decide +kernel) ms), fun h => absurd h (by decide)⟩
  · exact ⟨fun ⟨ms, h⟩ => absurd h (noPerfect_sound (by decide +kernel) ms), fun h => absurd h (by decide)⟩

theorem small_xdpm : ∀ N, N < 12 → ((∃ ms, Perfect (⟨.mul, .div, .add, .sub⟩ : Order) N ms) ↔ N ∈ [4, 5, 6, 8, 9, 10]) := by
  intro N hN
  rcases (by omega : N = 0 ∨ N = 1 ∨ N = 2 ∨ N = 3 ∨ N = 4 ∨ N = 5 ∨ N = 6 ∨ N = 7 ∨ N = 8 ∨ N = 9 ∨ N = 10 ∨ N = 11) with rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl
  · exact ⟨fun ⟨_, h⟩ => absurd h.1 (by omega), fun h => absurd h (by decide)⟩
  · exact ⟨fun ⟨ms, h⟩ => absurd h (noPerfect_sound (by decide +kernel) ms), fun h => absurd h (by decide)⟩
  · exact ⟨fun ⟨ms, h⟩ => absurd h (noPerfect_sound (by decide +kernel) ms), fun h => absurd h (by decide)⟩
  · exact ⟨fun ⟨ms, h⟩ => absurd h (noPerfect_sound (by decide +kernel) ms), fun h => absurd h (by decide)⟩
  · exact ⟨fun _ => by decide, fun _ => ⟨[⟨2, 2, 4⟩, ⟨2, 2, 1⟩, ⟨2, 1, 3⟩], (check_iff _ _ _).1 (by decide +kernel)⟩⟩
  · exact ⟨fun _ => by decide, fun _ => ⟨[⟨2, 2, 4⟩, ⟨2, 2, 1⟩, ⟨4, 1, 5⟩, ⟨4, 1, 3⟩], (check_iff _ _ _).1 (by decide +kernel)⟩⟩
  · exact ⟨fun _ => by decide, fun _ => ⟨[⟨2, 2, 4⟩, ⟨2, 2, 1⟩, ⟨4, 1, 5⟩, ⟨4, 1, 3⟩, ⟨2, 3, 6⟩], (check_iff _ _ _).1 (by decide +kernel)⟩⟩
  · exact ⟨fun ⟨ms, h⟩ => absurd h (noPerfect_sound (by decide +kernel) ms), fun h => absurd h (by decide)⟩
  · exact ⟨fun _ => by decide, fun _ => ⟨[⟨2, 2, 4⟩, ⟨2, 2, 1⟩, ⟨2, 4, 6⟩, ⟨6, 1, 5⟩, ⟨2, 4, 8⟩, ⟨6, 2, 3⟩, ⟨2, 5, 7⟩], (check_iff _ _ _).1 (by decide +kernel)⟩⟩
  · exact ⟨fun _ => by decide, fun _ => ⟨[⟨2, 2, 4⟩, ⟨2, 2, 1⟩, ⟨2, 4, 6⟩, ⟨6, 1, 5⟩, ⟨2, 4, 8⟩, ⟨6, 2, 3⟩, ⟨4, 5, 9⟩, ⟨8, 1, 7⟩], (check_iff _ _ _).1 (by decide +kernel)⟩⟩
  · exact ⟨fun _ => by decide, fun _ => ⟨[⟨2, 2, 4⟩, ⟨2, 2, 1⟩, ⟨2, 4, 6⟩, ⟨6, 1, 5⟩, ⟨2, 5, 10⟩, ⟨6, 2, 3⟩, ⟨4, 5, 9⟩, ⟨10, 3, 7⟩, ⟨2, 4, 8⟩], (check_iff _ _ _).1 (by decide +kernel)⟩⟩
  · exact ⟨fun ⟨ms, h⟩ => absurd h (noPerfect_sound (by decide +kernel) ms), fun h => absurd h (by decide)⟩

theorem small_xdmp : ∀ N, N < 12 → ((∃ ms, Perfect (⟨.mul, .div, .sub, .add⟩ : Order) N ms) ↔ N ∈ [4, 5, 6]) := by
  intro N hN
  rcases (by omega : N = 0 ∨ N = 1 ∨ N = 2 ∨ N = 3 ∨ N = 4 ∨ N = 5 ∨ N = 6 ∨ N = 7 ∨ N = 8 ∨ N = 9 ∨ N = 10 ∨ N = 11) with rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl
  · exact ⟨fun ⟨_, h⟩ => absurd h.1 (by omega), fun h => absurd h (by decide)⟩
  · exact ⟨fun ⟨ms, h⟩ => absurd h (noPerfect_sound (by decide +kernel) ms), fun h => absurd h (by decide)⟩
  · exact ⟨fun ⟨ms, h⟩ => absurd h (noPerfect_sound (by decide +kernel) ms), fun h => absurd h (by decide)⟩
  · exact ⟨fun ⟨ms, h⟩ => absurd h (noPerfect_sound (by decide +kernel) ms), fun h => absurd h (by decide)⟩
  · exact ⟨fun _ => by decide, fun _ => ⟨[⟨2, 2, 4⟩, ⟨2, 2, 1⟩, ⟨4, 1, 3⟩], (check_iff _ _ _).1 (by decide +kernel)⟩⟩
  · exact ⟨fun _ => by decide, fun _ => ⟨[⟨2, 2, 4⟩, ⟨2, 2, 1⟩, ⟨4, 1, 3⟩, ⟨2, 3, 5⟩], (check_iff _ _ _).1 (by decide +kernel)⟩⟩
  · exact ⟨fun _ => by decide, fun _ => ⟨[⟨2, 2, 4⟩, ⟨2, 2, 1⟩, ⟨4, 1, 3⟩, ⟨2, 3, 5⟩, ⟨2, 3, 6⟩], (check_iff _ _ _).1 (by decide +kernel)⟩⟩
  · exact ⟨fun ⟨ms, h⟩ => absurd h (noPerfect_sound (by decide +kernel) ms), fun h => absurd h (by decide)⟩
  · exact ⟨fun ⟨ms, h⟩ => absurd h (noPerfect_sound (by decide +kernel) ms), fun h => absurd h (by decide)⟩
  · exact ⟨fun ⟨ms, h⟩ => absurd h (noPerfect_sound (by decide +kernel) ms), fun h => absurd h (by decide)⟩
  · exact ⟨fun ⟨ms, h⟩ => absurd h (noPerfect_sound (by decide +kernel) ms), fun h => absurd h (by decide)⟩
  · exact ⟨fun ⟨ms, h⟩ => absurd h (noPerfect_sound (by decide +kernel) ms), fun h => absurd h (by decide)⟩

end Orders
