import Riemann.Deep.EvoAlg.Lean.Examples

/-!
# Corollary A (NOTE §0), complete in characteristic 0: the survey's statement holds over `𝕜` for
all finite connected graphs iff every positive rational number is a cube in `𝕜`; `K_{2,16}` example
-/

namespace EvoAlg

theorem cbg_neighborFinset_inl (m n : ℕ) (x : Fin m) :
    (completeBipartiteGraph (Fin m) (Fin n)).neighborFinset (Sum.inl x)
      = Finset.univ.map Function.Embedding.inr := by
  ext w
  rcases w with y | y <;> simp

theorem cbg_neighborFinset_inr (m n : ℕ) (y : Fin n) :
    (completeBipartiteGraph (Fin m) (Fin n)).neighborFinset (Sum.inr y)
      = Finset.univ.map Function.Embedding.inl := by
  ext w
  rcases w with x | x <;> simp

theorem cbg_degree_inl (m n : ℕ) (x : Fin m) :
    (completeBipartiteGraph (Fin m) (Fin n)).degree (Sum.inl x) = n := by
  rw [← SimpleGraph.card_neighborFinset_eq_degree, cbg_neighborFinset_inl, Finset.card_map,
    Finset.card_univ, Fintype.card_fin]

theorem cbg_degree_inr (m n : ℕ) (y : Fin n) :
    (completeBipartiteGraph (Fin m) (Fin n)).degree (Sum.inr y) = m := by
  rw [← SimpleGraph.card_neighborFinset_eq_degree, cbg_neighborFinset_inr, Finset.card_map,
    Finset.card_univ, Fintype.card_fin]

/-- `K_{m,n}` is biregular: degree `n` on the `m` left vertices, degree `m` on the `n` right ones. -/
theorem cbg_biregular (m n : ℕ) :
    IsBiregular (completeBipartiteGraph (Fin m) (Fin n)) (Set.range Sum.inl) (Set.range Sum.inr)
      n m where
  bip := SimpleGraph.IsBipartiteWith.completeBipartiteGraph (Fin m) (Fin n)
  cover v := by
    rcases v with a | b
    · exact Or.inl ⟨a, rfl⟩
    · exact Or.inr ⟨b, rfl⟩
  deg₁ v hv := by
    obtain ⟨a, rfl⟩ := hv
    exact cbg_degree_inl m n a
  deg₂ v hv := by
    obtain ⟨b, rfl⟩ := hv
    exact cbg_degree_inr m n b

/-- `A(K_{2,16}) ≅ A_RW(K_{2,16})` over `ℚ`: ratio `16/2 = 2³`. -/
theorem K2_16_iso : AIsoRW (completeBipartiteGraph (Fin 2) (Fin 16)) ℚ :=
  aIsoRW_of_biregular _ (cbg_biregular 2 16) (c := (2 : ℚ)) two_ne_zero (by norm_num)

variable {𝕜 : Type*} [Field 𝕜]

/-- The survey's statement over `𝕜`, for all finite connected graphs with `n ≥ 2` vertices
(vertex types in `Type`). -/
def SurveyHolds (𝕜 : Type*) [Field 𝕜] : Prop :=
  ∀ (V : Type) [Fintype V] [DecidableEq V] (G : SimpleGraph V) [DecidableRel G.Adj],
    G.Connected → 2 ≤ Fintype.card V →
      (AIsoRW G 𝕜 ↔ IsRegular G ∨ ∃ (s t : Set V) (d₁ d₂ : ℕ), IsBiregular G s t d₁ d₂)

/-- **Corollary A**, necessity: if the survey's statement holds over `𝕜` (characteristic 0), every
positive rational number is a cube in `𝕜` (apply it to `K_{b,a}`, ratio `a/b`). -/
theorem corA_necessity [CharZero 𝕜] (h : SurveyHolds 𝕜) :
    ∀ a b : ℕ, 0 < a → 0 < b → ∃ t : 𝕜, t ^ 3 = (a : 𝕜) / b := by
  intro a b ha hb
  have : NeZero a := ⟨ha.ne'⟩
  have : NeZero b := ⟨hb.ne'⟩
  have hcard : 2 ≤ Fintype.card (Fin b ⊕ Fin a) := by
    rw [Fintype.card_sum, Fintype.card_fin, Fintype.card_fin]
    omega
  have hiso : AIsoRW (completeBipartiteGraph (Fin b) (Fin a)) 𝕜 :=
    (h _ _ (cbg_connected b a) hcard).mpr (Or.inr ⟨_, _, _, _, cbg_biregular b a⟩)
  have ha0 : (a : 𝕜) ≠ 0 := Nat.cast_ne_zero.mpr ha.ne'
  have hb0 : (b : 𝕜) ≠ 0 := Nat.cast_ne_zero.mpr hb.ne'
  rcases (theoremA _ (cbg_connected b a) hcard).mp hiso with
    ⟨d, hd⟩ | ⟨s, t, d₁, d₂, hbi, c, hc, hcube⟩
  · have e1 := hd (Sum.inl 0)
    have e2 := hd (Sum.inr 0)
    rw [cbg_degree_inl] at e1
    rw [cbg_degree_inr] at e2
    refine ⟨1, ?_⟩
    rw [e1, ← e2, one_pow, div_self hb0]
  · have hadj : (completeBipartiteGraph (Fin b) (Fin a)).Adj (Sum.inl 0) (Sum.inr 0) :=
      Or.inl ⟨rfl, rfl⟩
    rcases hbi.bip.mem_of_adj hadj with ⟨e1, e2⟩ | ⟨e1, e2⟩
    · have f1 := hbi.deg₁ _ e1
      have f2 := hbi.deg₂ _ e2
      rw [cbg_degree_inl] at f1
      rw [cbg_degree_inr] at f2
      subst f1 f2
      exact ⟨c, hcube.symm⟩
    · have f1 := hbi.deg₂ _ e1
      have f2 := hbi.deg₁ _ e2
      rw [cbg_degree_inl] at f1
      rw [cbg_degree_inr] at f2
      subst f1 f2
      refine ⟨c⁻¹, ?_⟩
      rw [inv_pow, ← hcube, inv_div]

/-- **Corollary A** (NOTE §0, second item): in characteristic 0 the survey's statement holds over
`𝕜` for all finite connected graphs iff every positive rational number is a cube in `𝕜`. -/
theorem corA_iff [CharZero 𝕜] :
    SurveyHolds 𝕜 ↔ ∀ a b : ℕ, 0 < a → 0 < b → ∃ t : 𝕜, t ^ 3 = (a : 𝕜) / b :=
  ⟨corA_necessity, fun hcub _ _ _ G _ hG hn => corA_cubes G hcub hG hn⟩

/-- **Corollary A** (first item), abstractly: over a field of characteristic 0 in which every
element is a cube (`ℝ`, `ℂ`, every algebraically closed field) the survey's statement holds. -/
theorem corA_allCubes [CharZero 𝕜] (h3 : ∀ x : 𝕜, ∃ z, z ^ 3 = x) : SurveyHolds 𝕜 :=
  corA_iff.mpr (fun _ _ _ _ => h3 _)

/-- The survey's statement is false over `ℚ`. -/
theorem not_surveyHolds_rat : ¬ SurveyHolds ℚ := fun h =>
  corA_K12 ((h _ K12 K12_connected (by decide)).mpr (Or.inr ⟨_, _, _, _, K12_biregular⟩))

end EvoAlg

