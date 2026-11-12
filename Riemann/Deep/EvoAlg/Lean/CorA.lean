import Riemann.Deep.EvoAlg.Lean.TheoremA
import Mathlib.Algebra.Field.Rat

/-!
# Corollary A (NOTE §0): over `ℚ` the survey's statement fails; it holds over every field of
characteristic 0 in which every positive rational number is a cube (e.g. `ℝ`, `ℂ`)

`K_{1,2} = P₃` is Mathlib's `completeBipartiteGraph (Fin 1) (Fin 2)`: biregular with degrees
`2` (centre) and `1` (leaves); `2` and `1/2` are not rational cubes, so `A(K_{1,2}) ≇ A_RW(K_{1,2})`
over `ℚ` by Theorem A.
-/

open Finset

namespace EvoAlg

instance (V W : Type*) : DecidableRel (completeBipartiteGraph V W).Adj := fun v w =>
  inferInstanceAs (Decidable (v.isLeft ∧ w.isRight ∨ v.isRight ∧ w.isLeft))

/-- `K_{1,2} = P₃`. -/
abbrev K12 : SimpleGraph (Fin 1 ⊕ Fin 2) := completeBipartiteGraph (Fin 1) (Fin 2)

theorem K12_deg_inl : K12.degree (Sum.inl 0) = 2 := by decide

theorem K12_deg_inr : ∀ j : Fin 2, K12.degree (Sum.inr j) = 1 := by decide

theorem K12_connected : K12.Connected := by
  have h : ∀ v, K12.Reachable (Sum.inl 0) v := by
    intro v
    rcases v with a | b
    · obtain rfl : a = 0 := Subsingleton.elim a 0
      exact SimpleGraph.Reachable.refl _
    · exact SimpleGraph.Adj.reachable (Or.inl ⟨rfl, rfl⟩)
  exact (SimpleGraph.connected_iff_exists_forall_reachable _).mpr ⟨Sum.inl 0, h⟩

/-- `K_{1,2}` is biregular with sides `{centre}`, `{leaves}` and degrees `2`, `1`. -/
theorem K12_biregular : IsBiregular K12 (Set.range Sum.inl) (Set.range Sum.inr) 2 1 where
  bip := SimpleGraph.IsBipartiteWith.completeBipartiteGraph (Fin 1) (Fin 2)
  cover v := by
    rcases v with a | b
    · exact Or.inl ⟨a, rfl⟩
    · exact Or.inr ⟨b, rfl⟩
  deg₁ v hv := by
    obtain ⟨a, rfl⟩ := hv
    obtain rfl : a = 0 := Subsingleton.elim a 0
    exact K12_deg_inl
  deg₂ v hv := by
    obtain ⟨b, rfl⟩ := hv
    exact K12_deg_inr b

theorem nat_cube_ne_two (m : ℕ) : m ^ 3 ≠ 2 := by
  intro h
  match m, h with
  | 0, h => simp at h
  | 1, h => simp at h
  | m + 2, h =>
    have h8 : 2 ^ 3 ≤ (m + 2) ^ 3 := Nat.pow_le_pow_left (by omega) 3
    rw [h] at h8
    exact absurd h8 (by decide)

/-- `2` is not a rational cube. -/
theorem not_cube_two : ¬ ∃ t : ℚ, t ^ 3 = 2 := by
  rintro ⟨t, ht⟩
  have hnum : t.num ^ 3 = 2 := by
    have := congrArg Rat.num ht
    rw [Rat.num_pow] at this
    exact this
  have h3 : t.num.natAbs ^ 3 = 2 := by
    have := congrArg Int.natAbs hnum
    rwa [Int.natAbs_pow] at this
  exact nat_cube_ne_two _ h3

/-- **Corollary A over `ℚ`**: `A(K_{1,2}) ≇ A_RW(K_{1,2})` over `ℚ`. -/
theorem corA_K12 : ¬ AIsoRW K12 ℚ := by
  intro h
  rcases (theoremA K12 K12_connected (by decide)).mp h with
    ⟨d, hd⟩ | ⟨s, t, d₁, d₂, hb, c, hc, hcube⟩
  · have h1 := hd (Sum.inl 0)
    have h2 := hd (Sum.inr 0)
    rw [K12_deg_inl] at h1
    rw [K12_deg_inr] at h2
    omega
  · have hadj : K12.Adj (Sum.inl 0) (Sum.inr 0) := Or.inl ⟨rfl, rfl⟩
    rcases hb.bip.mem_of_adj hadj with ⟨h1, h2⟩ | ⟨h1, h2⟩
    · have e1 := hb.deg₁ _ h1
      have e2 := hb.deg₂ _ h2
      rw [K12_deg_inl] at e1
      rw [K12_deg_inr] at e2
      subst e1 e2
      exact not_cube_two ⟨c, by rw [← hcube]; norm_num⟩
    · have e1 := hb.deg₂ _ h1
      have e2 := hb.deg₁ _ h2
      rw [K12_deg_inl] at e1
      rw [K12_deg_inr] at e2
      subst e1 e2
      exact not_cube_two ⟨c⁻¹, by rw [inv_pow, ← hcube]; norm_num⟩

/-- The survey's statement (arXiv 2610.22248, Question "Graphs versus random walks", over an
arbitrary field) fails over `ℚ`: `K_{1,2}` is connected, biregular, and `A(K_{1,2}) ≇ A_RW(K_{1,2})`. -/
theorem survey_false_over_Q : K12.Connected ∧ 2 ≤ Fintype.card (Fin 1 ⊕ Fin 2) ∧
    (∃ (s t : Set (Fin 1 ⊕ Fin 2)) (d₁ d₂ : ℕ), IsBiregular K12 s t d₁ d₂) ∧ ¬ AIsoRW K12 ℚ :=
  ⟨K12_connected, by decide, ⟨_, _, _, _, K12_biregular⟩, corA_K12⟩

variable {V : Type*} [Fintype V] [DecidableEq V] {𝕜 : Type*} [Field 𝕜]
variable (G : SimpleGraph V) [DecidableRel G.Adj]

/-- **Corollary A**, positive part: over a field of characteristic 0 in which every positive
rational number is a cube (e.g. `ℝ`, `ℂ`), `A(G) ≅ A_RW(G)` iff `G` is regular or biregular. -/
theorem corA_cubes [CharZero 𝕜]
    (hcub : ∀ a b : ℕ, 0 < a → 0 < b → ∃ t : 𝕜, t ^ 3 = (a : 𝕜) / b)
    (hG : G.Connected) (hn : 2 ≤ Fintype.card V) :
    AIsoRW G 𝕜 ↔ IsRegular G ∨ ∃ (s t : Set V) (d₁ d₂ : ℕ), IsBiregular G s t d₁ d₂ := by
  rw [theoremA G hG hn]
  refine or_congr Iff.rfl ⟨fun ⟨s, t, d₁, d₂, hb, _⟩ => ⟨s, t, d₁, d₂, hb⟩, ?_⟩
  rintro ⟨s, t, d₁, d₂, hb⟩
  obtain ⟨v₀⟩ := hG.nonempty
  obtain ⟨w₀, hw₀⟩ := exists_adj_of_connected G hG hn v₀
  have hpos : ∀ v, 0 < G.degree v := fun v => by
    obtain ⟨w, hw⟩ := exists_adj_of_connected G hG hn v
    exact (G.degree_pos_iff_exists_adj v).mpr ⟨w, hw⟩
  have hd : 0 < d₁ ∧ 0 < d₂ := by
    rcases hb.bip.mem_of_adj hw₀ with ⟨h1, h2⟩ | ⟨h1, h2⟩
    · exact ⟨hb.deg₁ _ h1 ▸ hpos v₀, hb.deg₂ _ h2 ▸ hpos w₀⟩
    · exact ⟨hb.deg₁ _ h2 ▸ hpos w₀, hb.deg₂ _ h1 ▸ hpos v₀⟩
  obtain ⟨c, hc⟩ := hcub d₁ d₂ hd.1 hd.2
  have hc0 : c ≠ 0 := by
    rintro rfl
    have h0 : (d₁ : 𝕜) / d₂ ≠ 0 := div_ne_zero (Nat.cast_ne_zero.mpr hd.1.ne')
      (Nat.cast_ne_zero.mpr hd.2.ne')
    exact h0 (by rw [← hc]; ring)
  exact ⟨s, t, d₁, d₂, hb, c, hc0, hc.symm⟩

end EvoAlg

