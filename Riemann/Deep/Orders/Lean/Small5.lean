import Riemann.Deep.Orders.Lean.Search

/-! Exact sets of perfect N below N0 (working orders) or below 6 (failing orders):
witnesses by the checker, non-existence by the sound search noPerfect. -/
namespace Orders

theorem small_mdpx : ∀ N, N < 6 → ((∃ ms, Perfect (⟨.sub, .div, .add, .mul⟩ : Order) N ms) ↔ N ∈ []) := by
  intro N hN
  rcases (by omega : N = 0 ∨ N = 1 ∨ N = 2 ∨ N = 3 ∨ N = 4 ∨ N = 5) with rfl | rfl | rfl | rfl | rfl | rfl
  · exact ⟨fun ⟨_, h⟩ => absurd h.1 (by omega), fun h => absurd h (by decide)⟩
  · exact ⟨fun ⟨ms, h⟩ => absurd h (noPerfect_sound (by decide +kernel) ms), fun h => absurd h (by decide)⟩
  · exact ⟨fun ⟨ms, h⟩ => absurd h (noPerfect_sound (by decide +kernel) ms), fun h => absurd h (by decide)⟩
  · exact ⟨fun ⟨ms, h⟩ => absurd h (noPerfect_sound (by decide +kernel) ms), fun h => absurd h (by decide)⟩
  · exact ⟨fun ⟨ms, h⟩ => absurd h (noPerfect_sound (by decide +kernel) ms), fun h => absurd h (by decide)⟩
  · exact ⟨fun ⟨ms, h⟩ => absurd h (noPerfect_sound (by decide +kernel) ms), fun h => absurd h (by decide)⟩

theorem small_mdxp : ∀ N, N < 6 → ((∃ ms, Perfect (⟨.sub, .div, .mul, .add⟩ : Order) N ms) ↔ N ∈ []) := by
  intro N hN
  rcases (by omega : N = 0 ∨ N = 1 ∨ N = 2 ∨ N = 3 ∨ N = 4 ∨ N = 5) with rfl | rfl | rfl | rfl | rfl | rfl
  · exact ⟨fun ⟨_, h⟩ => absurd h.1 (by omega), fun h => absurd h (by decide)⟩
  · exact ⟨fun ⟨ms, h⟩ => absurd h (noPerfect_sound (by decide +kernel) ms), fun h => absurd h (by decide)⟩
  · exact ⟨fun ⟨ms, h⟩ => absurd h (noPerfect_sound (by decide +kernel) ms), fun h => absurd h (by decide)⟩
  · exact ⟨fun ⟨ms, h⟩ => absurd h (noPerfect_sound (by decide +kernel) ms), fun h => absurd h (by decide)⟩
  · exact ⟨fun ⟨ms, h⟩ => absurd h (noPerfect_sound (by decide +kernel) ms), fun h => absurd h (by decide)⟩
  · exact ⟨fun ⟨ms, h⟩ => absurd h (noPerfect_sound (by decide +kernel) ms), fun h => absurd h (by decide)⟩

theorem small_dpxm : ∀ N, N < 9 → ((∃ ms, Perfect (⟨.div, .add, .mul, .sub⟩ : Order) N ms) ↔ N ∈ [2, 3, 4]) := by
  intro N hN
  rcases (by omega : N = 0 ∨ N = 1 ∨ N = 2 ∨ N = 3 ∨ N = 4 ∨ N = 5 ∨ N = 6 ∨ N = 7 ∨ N = 8) with rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl
  · exact ⟨fun ⟨_, h⟩ => absurd h.1 (by omega), fun h => absurd h (by decide)⟩
  · exact ⟨fun ⟨ms, h⟩ => absurd h (noPerfect_sound (by decide +kernel) ms), fun h => absurd h (by decide)⟩
  · exact ⟨fun _ => by decide, fun _ => ⟨[⟨2, 2, 1⟩], (check_iff _ _ _).1 (by decide +kernel)⟩⟩
  · exact ⟨fun _ => by decide, fun _ => ⟨[⟨2, 2, 1⟩, ⟨2, 1, 3⟩], (check_iff _ _ _).1 (by decide +kernel)⟩⟩
  · exact ⟨fun _ => by decide, fun _ => ⟨[⟨2, 2, 1⟩, ⟨2, 1, 3⟩, ⟨2, 2, 4⟩], (check_iff _ _ _).1 (by decide +kernel)⟩⟩
  · exact ⟨fun ⟨ms, h⟩ => absurd h (noPerfect_sound (by decide +kernel) ms), fun h => absurd h (by decide)⟩
  · exact ⟨fun ⟨ms, h⟩ => absurd h (noPerfect_sound (by decide +kernel) ms), fun h => absurd h (by decide)⟩
  · exact ⟨fun ⟨ms, h⟩ => absurd h (noPerfect_sound (by decide +kernel) ms), fun h => absurd h (by decide)⟩
  · exact ⟨fun ⟨ms, h⟩ => absurd h (noPerfect_sound (by decide +kernel) ms), fun h => absurd h (by decide)⟩

theorem small_dpmx : ∀ N, N < 12 → ((∃ ms, Perfect (⟨.div, .add, .sub, .mul⟩ : Order) N ms) ↔ N ∈ [2, 3, 4]) := by
  intro N hN
  rcases (by omega : N = 0 ∨ N = 1 ∨ N = 2 ∨ N = 3 ∨ N = 4 ∨ N = 5 ∨ N = 6 ∨ N = 7 ∨ N = 8 ∨ N = 9 ∨ N = 10 ∨ N = 11) with rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl
  · exact ⟨fun ⟨_, h⟩ => absurd h.1 (by omega), fun h => absurd h (by decide)⟩
  · exact ⟨fun ⟨ms, h⟩ => absurd h (noPerfect_sound (by decide +kernel) ms), fun h => absurd h (by decide)⟩
  · exact ⟨fun _ => by decide, fun _ => ⟨[⟨2, 2, 1⟩], (check_iff _ _ _).1 (by decide +kernel)⟩⟩
  · exact ⟨fun _ => by decide, fun _ => ⟨[⟨2, 2, 1⟩, ⟨2, 1, 3⟩], (check_iff _ _ _).1 (by decide +kernel)⟩⟩
  · exact ⟨fun _ => by decide, fun _ => ⟨[⟨2, 2, 1⟩, ⟨2, 2, 4⟩, ⟨4, 1, 3⟩], (check_iff _ _ _).1 (by decide +kernel)⟩⟩
  · exact ⟨fun ⟨ms, h⟩ => absurd h (noPerfect_sound (by decide +kernel) ms), fun h => absurd h (by decide)⟩
  · exact ⟨fun ⟨ms, h⟩ => absurd h (noPerfect_sound (by decide +kernel) ms), fun h => absurd h (by decide)⟩
  · exact ⟨fun ⟨ms, h⟩ => absurd h (noPerfect_sound (by decide +kernel) ms), fun h => absurd h (by decide)⟩
  · exact ⟨fun ⟨ms, h⟩ => absurd h (noPerfect_sound (by decide +kernel) ms), fun h => absurd h (by decide)⟩
  · exact ⟨fun ⟨ms, h⟩ => absurd h (noPerfect_sound (by decide +kernel) ms), fun h => absurd h (by decide)⟩
  · exact ⟨fun ⟨ms, h⟩ => absurd h (noPerfect_sound (by decide +kernel) ms), fun h => absurd h (by decide)⟩
  · exact ⟨fun ⟨ms, h⟩ => absurd h (noPerfect_sound (by decide +kernel) ms), fun h => absurd h (by decide)⟩

end Orders
