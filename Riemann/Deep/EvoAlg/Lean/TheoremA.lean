import Riemann.Deep.EvoAlg.Lean.Cor41

/-!
# §5 and Theorem A (NOTE §0, §5)
-/

open Finset

namespace EvoAlg

variable {V : Type*} [Fintype V] [DecidableEq V] {𝕜 : Type*} [Field 𝕜]
variable (G : SimpleGraph V) [DecidableRel G.Adj]

/-- NOTE §5: in characteristic 0, on a connected graph without isolated vertices, a solution of
the scaling system (i), (ii) forces `G` regular, or bipartite biregular with `d₁/d₂ = (c_Y/c_X)³`.
(`c` alternates between `x = c_{v₀}` and `y = c_{u₀}` along walks; `x = y` gives `d ≡ x`, `x ≠ y`
gives the bipartition `{c = x} ⊔ {c = y}` with `d = y²/x` and `d = x²/y` on its sides.) -/
theorem sec5 [CharZero 𝕜] (hG : G.Connected) (hnb : ∀ v, ∃ w, G.Adj v w) (c : V → 𝕜)
    (hc0 : ∀ v, c v ≠ 0) (hi : ∀ u v v', G.Adj u v → G.Adj u v' → c v = c v')
    (hii : ∀ u v, G.Adj u v → c v ^ 2 = c u * G.degree u) :
    IsRegular G ∨ BiregularCube G 𝕜 := by
  obtain ⟨v₀⟩ := hG.nonempty
  obtain ⟨u₀, hu₀⟩ := hnb v₀
  have hP : ∀ w, (c w = c v₀ ∧ ∀ z, G.Adj w z → c z = c u₀) ∨
      (c w = c u₀ ∧ ∀ z, G.Adj w z → c z = c v₀) := by
    intro w
    have hr : Relation.ReflTransGen G.Adj v₀ w :=
      (SimpleGraph.reachable_iff_reflTransGen v₀ w).mp (hG.preconnected v₀ w)
    induction hr with
    | refl => exact Or.inl ⟨rfl, fun z hz => hi v₀ z u₀ hz hu₀⟩
    | tail _ hbc ih =>
      rcases ih with ⟨h1, h2⟩ | ⟨h1, h2⟩
      · exact Or.inr ⟨h2 _ hbc, fun z hz => (hi _ z _ hz hbc.symm).trans h1⟩
      · exact Or.inl ⟨h2 _ hbc, fun z hz => (hi _ z _ hz hbc.symm).trans h1⟩
  have hdpos : ∀ w, G.degree w ≠ 0 := fun w => by
    obtain ⟨z, hz⟩ := hnb w
    exact Nat.pos_iff_ne_zero.mp ((G.degree_pos_iff_exists_adj w).mpr ⟨z, hz⟩)
  have hdw : ∀ w z, G.Adj w z → (G.degree w : 𝕜) * c w = c z ^ 2 := fun w z hz => by
    rw [hii w z hz, mul_comm]
  by_cases hxy : c v₀ = c u₀
  · left
    have hconst : ∀ w, c w = c v₀ := fun w => by
      rcases hP w with ⟨h, _⟩ | ⟨h, _⟩
      · exact h
      · exact h.trans hxy.symm
    have hd : ∀ w, (G.degree w : 𝕜) = c v₀ := fun w => by
      obtain ⟨z, hz⟩ := hnb w
      have h := hdw w z hz
      rw [hconst w, hconst z, sq] at h
      exact mul_right_cancel₀ (hc0 v₀) h
    refine ⟨G.degree v₀, fun w => ?_⟩
    exact_mod_cast (hd w).trans (hd v₀).symm
  · right
    have hS : ∀ w, c w = c v₀ → ∀ z, G.Adj w z → c z = c u₀ := fun w hw z hz => by
      rcases hP w with ⟨_, h⟩ | ⟨h1, _⟩
      · exact h z hz
      · exact absurd (hw.symm.trans h1) hxy
    have hT : ∀ w, c w = c u₀ → ∀ z, G.Adj w z → c z = c v₀ := fun w hw z hz => by
      rcases hP w with ⟨h1, _⟩ | ⟨_, h⟩
      · exact absurd (h1.symm.trans hw) hxy
      · exact h z hz
    have hdS : ∀ w, c w = c v₀ → (G.degree w : 𝕜) * c v₀ = c u₀ ^ 2 := fun w hw => by
      obtain ⟨z, hz⟩ := hnb w
      rw [← hw, ← hS w hw z hz]
      exact hdw w z hz
    have hdT : ∀ w, c w = c u₀ → (G.degree w : 𝕜) * c u₀ = c v₀ ^ 2 := fun w hw => by
      obtain ⟨z, hz⟩ := hnb w
      rw [← hw, ← hT w hw z hz]
      exact hdw w z hz
    have hd₂ne : (G.degree u₀ : 𝕜) ≠ 0 := Nat.cast_ne_zero.mpr (hdpos u₀)
    have hx0 := hc0 v₀
    have hy0 := hc0 u₀
    refine ⟨{w | c w = c v₀}, {w | c w = c u₀}, G.degree v₀, G.degree u₀,
      ⟨⟨?_, ?_⟩, ?_, ?_, ?_⟩, c u₀ / c v₀, div_ne_zero hy0 hx0, ?_⟩
    · exact Set.disjoint_left.mpr
        (fun w (h1 : c w = c v₀) (h2 : c w = c u₀) => hxy (h1.symm.trans h2))
    · intro v w hvw
      rcases hP v with ⟨h1, h2⟩ | ⟨h1, h2⟩
      · exact Or.inl ⟨h1, h2 w hvw⟩
      · exact Or.inr ⟨h1, h2 w hvw⟩
    · intro w
      rcases hP w with ⟨h1, _⟩ | ⟨h1, _⟩
      · exact Or.inl h1
      · exact Or.inr h1
    · intro w (hw : c w = c v₀)
      exact_mod_cast mul_right_cancel₀ hx0 ((hdS w hw).trans (hdS v₀ rfl).symm)
    · intro w (hw : c w = c u₀)
      exact_mod_cast mul_right_cancel₀ hy0 ((hdT w hw).trans (hdT u₀ rfl).symm)
    · have h1 := hdS v₀ rfl
      have h2 := hdT u₀ rfl
      rw [div_pow, div_eq_div_iff hd₂ne (pow_ne_zero 3 hx0)]
      linear_combination (c v₀ ^ 2) * h1 - (c u₀ ^ 2) * h2

/-- **Theorem A** (NOTE §0, §5). Let `char 𝕜 = 0` and `G` a finite connected graph with `n ≥ 2`
vertices. Then `A(G) ≅ A_RW(G)` over `𝕜` iff `G` is regular, or `G` is bipartite with every vertex
of side `V_i` of degree `d_i` and `d₁/d₂ ∈ (𝕜^×)³`. -/
theorem theoremA [CharZero 𝕜] (hG : G.Connected) (hn : 2 ≤ Fintype.card V) :
    AIsoRW G 𝕜 ↔ IsRegular G ∨ BiregularCube G 𝕜 := by
  refine ⟨fun h => ?_, theoremA_if G hG hn⟩
  have hnb := exists_adj_of_connected G hG hn
  have hdeg : ∀ v, (G.degree v : 𝕜) ≠ 0 := fun v => by
    obtain ⟨w, hw⟩ := hnb v
    exact Nat.cast_ne_zero.mpr (Nat.pos_iff_ne_zero.mp ((G.degree_pos_iff_exists_adj v).mpr ⟨w, hw⟩))
  have hk : ∀ u, (#(twin G u) : 𝕜) ≠ 0 := fun u =>
    Nat.cast_ne_zero.mpr (Finset.card_ne_zero_of_mem ((mem_twin G).mpr rfl))
  obtain ⟨c, hc0, hi, hii⟩ := (cor41 G hnb hdeg hk).mp h
  exact sec5 G hG hnb c hc0 hi hii

end EvoAlg

