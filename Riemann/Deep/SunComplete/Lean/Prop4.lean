import Riemann.Deep.SunComplete.Lean.Graph
import Riemann.Deep.SunComplete.Lean.Arith

/-!
# sunlean · Prop4: Proposition 4 and Lemma 5 of NOTE §3, for an abstract graph

Input (`L3Data`): a pair count `e(X,Y)` obeying Lemma 3 in squared, √-free form
(`L ≤ τE√Z ⟺ L ≤ 0 ∨ L² ≤ τ²E²Z`), the diagonal bound `e(I,I) ≤ ε|I|` on independent sets,
and `e(A,B) = 0` across an empty cut. Output: a Hamiltonian cycle, via Chvátal–Erdős.
All arithmetic is over `ℚ`.
-/

namespace SunLean

variable {V : Type} [Fintype V] [DecidableEq V]

/-- The data of NOTE §§2–3 for a graph `G` on `n = card V` vertices (Lemma 3 in squared form). -/
structure L3Data (G : SimpleGraph V) where
  e : Finset V → Finset V → ℕ
  τ : ℚ
  E : ℚ
  Nt : ℚ
  ε : ℚ
  hτ : 0 < τ
  hE : 0 ≤ E
  hε : 0 ≤ ε
  hdiag : ∀ I : Finset V, IndepSet G I → (e I I : ℚ) ≤ ε * I.card
  hcross : ∀ A B : Finset V, Disjoint A B → (∀ x ∈ A, ∀ y ∈ B, ¬ G.Adj x y) → e A B = 0
  hL3 : ∀ X Y : Finset V,
    τ * (Nt * X.card * Y.card) - (Fintype.card V : ℚ) * (e X Y : ℚ) ≤ 0 ∨
    (τ * (Nt * X.card * Y.card) - (Fintype.card V : ℚ) * (e X Y : ℚ)) ^ 2 ≤
      (τ * E) ^ 2 * ((X.card * ((Fintype.card V : ℚ) - X.card)) *
        (Y.card * ((Fintype.card V : ℚ) - Y.card)))

/-- Degree of `v`. -/
def deg (G : SimpleGraph V) [DecidableRel G.Adj] (v : V) : ℕ := (Finset.univ.filter (G.Adj v)).card


/-! ### The graph argument -/

namespace L3Data

variable {G : SimpleGraph V}

/-- Proposition 4 (i), multiplied out: `|I|·τ(Ñ + E) ≤ n(ε + τE)` for independent `I`. -/
theorem indep_bound (D : L3Data G) {I : Finset V} (hI : IndepSet G I) :
    (I.card : ℚ) * (D.τ * (D.Nt + D.E)) ≤ (Fintype.card V : ℚ) * (D.ε + D.τ * D.E) :=
  indep_arith D.hτ D.hE D.hε (Nat.cast_nonneg _) (by exact_mod_cast Finset.card_le_univ I)
    (D.hdiag I hI) (D.hL3 I I)

/-- Proposition 4 (i): `α ≤ α* = n(ε/τ + E)/(Ñ + E)`. -/
theorem alpha_le (D : L3Data G) (hNE : 0 < D.Nt + D.E) {I : Finset V} (hI : IndepSet G I) :
    (I.card : ℚ) ≤ (Fintype.card V : ℚ) * (D.ε / D.τ + D.E) / (D.Nt + D.E) := by
  rw [le_div_iff₀ hNE]
  have hτ := D.hτ
  have hτ0 : D.τ ≠ 0 := hτ.ne'
  have hτε : D.τ * (D.ε / D.τ) = D.ε := by field_simp
  have e1 : D.τ * ((Fintype.card V : ℚ) * (D.ε / D.τ + D.E)) =
      (Fintype.card V : ℚ) * (D.ε + D.τ * D.E) := by
    rw [show D.τ * ((Fintype.card V : ℚ) * (D.ε / D.τ + D.E)) =
      (Fintype.card V : ℚ) * (D.τ * (D.ε / D.τ) + D.τ * D.E) by ring, hτε]
  have h1 := D.indep_bound hI
  have h3 : D.τ * ((I.card : ℚ) * (D.Nt + D.E)) ≤
      D.τ * ((Fintype.card V : ℚ) * (D.ε / D.τ + D.E)) := by
    rw [e1]; linarith [h1]
  exact le_of_mul_le_mul_left h3 hτ

variable [DecidableRel G.Adj]

/-- Proposition 4 (ii): for an empty cut `(A, B)` of `G − S` with `|A| ≤ |B|`:
`δ₀ ≤ a − 1 + c`, `2a + c ≤ n` and `Φ(a, c) ≤ 0`. -/
theorem cut_bound (D : L3Data G) {δ₀ : ℚ} (hdeg : ∀ v, δ₀ ≤ (deg G v : ℚ)) (hNt : 0 < D.Nt)
    {A B S : Finset V} (hAB : Disjoint A B) (hAS : Disjoint A S) (hBS : Disjoint B S)
    (hcov : ∀ w, w ∈ A ∨ w ∈ B ∨ w ∈ S) (hA : A.Nonempty) (hB : B.Nonempty)
    (hno : ∀ x ∈ A, ∀ y ∈ B, ¬ G.Adj x y) (hle : A.card ≤ B.card) :
    δ₀ ≤ (A.card : ℚ) - 1 + S.card ∧ 2 * (A.card : ℚ) + S.card ≤ Fintype.card V ∧
      D.Nt ^ 2 * A.card * ((Fintype.card V : ℚ) - A.card - S.card) ≤
        D.E ^ 2 * (A.card + S.card) * ((Fintype.card V : ℚ) - A.card) := by
  have hcard : A.card + B.card + S.card = Fintype.card V := by
    have hU : A ∪ B ∪ S = Finset.univ := by
      ext w
      simp only [Finset.mem_union, Finset.mem_univ, iff_true]
      rcases hcov w with h | h | h
      · exact Or.inl (Or.inl h)
      · exact Or.inl (Or.inr h)
      · exact Or.inr h
    rw [← Finset.card_univ, ← hU,
      Finset.card_union_of_disjoint (Finset.disjoint_union_left.mpr ⟨hAS, hBS⟩),
      Finset.card_union_of_disjoint hAB]
  obtain ⟨v, hv⟩ := hA
  have hsub : Finset.univ.filter (G.Adj v) ⊆ A.erase v ∪ S := by
    intro w hw
    simp only [Finset.mem_filter, Finset.mem_univ, true_and] at hw
    rcases hcov w with h | h | h
    · exact Finset.mem_union_left _ (Finset.mem_erase.mpr ⟨fun hwv => hw.ne hwv.symm, h⟩)
    · exact absurd hw (hno v hv w h)
    · exact Finset.mem_union_right _ h
  have hdeg' : deg G v + 1 ≤ A.card + S.card := by
    have h1 := Finset.card_le_card hsub
    have h2 := Finset.card_union_le (A.erase v) S
    have h3 := Finset.card_erase_add_one hv
    unfold deg
    omega
  have hd1 : δ₀ ≤ (A.card : ℚ) - 1 + S.card := by
    have h0 := hdeg v
    have h' : (deg G v : ℚ) + 1 ≤ A.card + S.card := by exact_mod_cast hdeg'
    linarith
  have hcardQ : (A.card : ℚ) + B.card + S.card = Fintype.card V := by exact_mod_cast hcard
  have hleQ : (A.card : ℚ) ≤ B.card := by exact_mod_cast hle
  refine ⟨hd1, by linarith, ?_⟩
  have he : (D.e A B : ℚ) = 0 := by rw [D.hcross A B hAB hno]; simp
  have hA1 : (1 : ℚ) ≤ A.card := by exact_mod_cast Finset.card_pos.mpr ⟨v, hv⟩
  have hB1 : (1 : ℚ) ≤ B.card := by exact_mod_cast Finset.card_pos.mpr hB
  exact cut_arith D.hτ hNt hA1 hB1 hcardQ.symm he (D.hL3 A B)

/-- The common skeleton: Chvátal–Erdős with `k`, given `α ≤ k` and `Φ(a, c) > 0` for every
admissible `(a, c)` with `c < k` (so no vertex cut of size `< k` exists). -/
theorem ham_of_bounds (hCE : ChvatalErdos) (D : L3Data G) {δ₀ : ℚ}
    (hdeg : ∀ v, δ₀ ≤ (deg G v : ℚ)) (hNt : 0 < D.Nt) (hn : 3 ≤ Fintype.card V) (k : ℕ)
    (hk : k ≤ Fintype.card V) (hindep : ∀ I : Finset V, IndepSet G I → I.card ≤ k)
    (hcut : ∀ a c : ℕ, c < k → 1 ≤ a → δ₀ ≤ (a : ℚ) - 1 + c → 2 * (a : ℚ) + c ≤ Fintype.card V →
      D.E ^ 2 * (a + c) * ((Fintype.card V : ℚ) - a) <
        D.Nt ^ 2 * a * ((Fintype.card V : ℚ) - a - c)) :
    HasHamCycle G := by
  refine hCE V G k hn hindep (fun S hS => ⟨?_, ?_⟩)
  · rw [← Finset.card_pos, Finset.card_sdiff_of_subset (Finset.subset_univ S), Finset.card_univ]
    omega
  · intro A hA hAne hBne
    by_contra hcon
    push Not at hcon
    have hAB : Disjoint A ((Finset.univ \ S) \ A) := Finset.disjoint_sdiff
    have hAS : Disjoint A S := Finset.disjoint_of_subset_left hA Finset.sdiff_disjoint
    have hBS : Disjoint ((Finset.univ \ S) \ A) S :=
      Finset.disjoint_of_subset_left Finset.sdiff_subset Finset.sdiff_disjoint
    have hcov : ∀ w, w ∈ A ∨ w ∈ (Finset.univ \ S) \ A ∨ w ∈ S := by
      intro w
      by_cases h1 : w ∈ S
      · exact Or.inr (Or.inr h1)
      · by_cases h2 : w ∈ A
        · exact Or.inl h2
        · exact Or.inr (Or.inl (by simp [h1, h2]))
    rcases le_total A.card ((Finset.univ \ S) \ A).card with hle | hle
    · obtain ⟨h1, h2, h3⟩ := D.cut_bound hdeg hNt hAB hAS hBS hcov hAne hBne hcon hle
      have h4 := hcut A.card S.card hS (Finset.card_pos.mpr hAne) h1 h2
      linarith
    · have hcov' : ∀ w, w ∈ (Finset.univ \ S) \ A ∨ w ∈ A ∨ w ∈ S := by
        intro w
        rcases hcov w with h | h | h
        · exact Or.inr (Or.inl h)
        · exact Or.inl h
        · exact Or.inr (Or.inr h)
      obtain ⟨h1, h2, h3⟩ := D.cut_bound hdeg hNt hAB.symm hBS hAS hcov' hBne hAne
        (fun x hx y hy hxy => hcon y hy x hx hxy.symm) hle
      have h4 := hcut _ S.card hS (Finset.card_pos.mpr hBne) h1 h2
      linarith

/-- A vertex has degree at most `n − 1`. -/
theorem deg_lt (v : V) : deg G v + 1 ≤ Fintype.card V := by
  have h1 : Finset.univ.filter (G.Adj v) ⊆ Finset.univ.erase v := by
    intro w hw
    simp only [Finset.mem_filter, Finset.mem_univ, true_and] at hw
    exact Finset.mem_erase.mpr ⟨fun h => hw.ne h.symm, Finset.mem_univ w⟩
  have h2 := Finset.card_le_card h1
  rw [Finset.card_erase_of_mem (Finset.mem_univ v), Finset.card_univ] at h2
  have h3 : 0 < Fintype.card V := Fintype.card_pos_iff.mpr ⟨v⟩
  unfold deg
  omega

/-- **Proposition 4** (NOTE §3). With `α* = n(ε/τ + E)/(Ñ + E)`, any natural `C` with
`α* < C + 2` (NOTE takes `C = ⌊α*⌋ − 1`), `a₁ = δ₀ − C + 1 ≥ 1` and either `a₁ > n/2` or
`F(a₁) > 0` and `F(n/2) > 0`: the graph is Hamiltonian. -/
theorem prop4 (hCE : ChvatalErdos) (D : L3Data G) {δ₀ : ℚ} (hdeg : ∀ v, δ₀ ≤ (deg G v : ℚ))
    (hn : 3 ≤ Fintype.card V) (hNE : D.E < D.Nt) (C : ℕ)
    (hC : (Fintype.card V : ℚ) * (D.ε / D.τ + D.E) / (D.Nt + D.E) < C + 2)
    (ha1 : 1 ≤ δ₀ - C + 1)
    (hF : (Fintype.card V : ℚ) / 2 < δ₀ - C + 1 ∨
      (0 < Fq D.Nt D.E (Fintype.card V) C (δ₀ - C + 1) ∧
        0 < Fq D.Nt D.E (Fintype.card V) C ((Fintype.card V : ℚ) / 2))) :
    HasHamCycle G := by
  have hNt : 0 < D.Nt := lt_of_le_of_lt D.hE hNE
  have hNE' : 0 < D.Nt + D.E := by linarith [D.hE]
  obtain ⟨v0⟩ : Nonempty V := Fintype.card_pos_iff.mp (by omega)
  have hdv : (deg G v0 : ℚ) + 1 ≤ Fintype.card V := by exact_mod_cast deg_lt v0
  have hCk : C + 1 ≤ Fintype.card V := by
    have : (C : ℚ) + 1 ≤ Fintype.card V := by linarith [hdeg v0]
    exact_mod_cast this
  refine D.ham_of_bounds hCE hdeg hNt hn (C + 1) hCk ?_ ?_
  · intro I hI
    have h1 : (I.card : ℚ) < (C : ℚ) + 2 := lt_of_le_of_lt (D.alpha_le hNE' hI) hC
    have h2 : I.card < C + 2 := by exact_mod_cast h1
    omega
  · intro a c hc ha hda h2a
    have hcC : (c : ℚ) ≤ C := by exact_mod_cast Nat.lt_succ_iff.mp hc
    exact prop4_arith D.hE hNE hcC hda h2a (Nat.cast_nonneg c) (by exact_mod_cast ha) hF

/-- **Lemma 5** (NOTE §3, clean form (CS)). For any `A' ≥ α*` and `D ≤ δ₀ + 1` with `D > A'`,
`2A' < n` and `(D − A')·Ñ²·(1 − 2A'/n) > E²·D`: the graph is Hamiltonian. -/
theorem lemma5 (hCE : ChvatalErdos) (D : L3Data G) {δ₀ : ℚ} (hdeg : ∀ v, δ₀ ≤ (deg G v : ℚ))
    (hn : 3 ≤ Fintype.card V) (hNt : 0 ≤ D.Nt) {A' Dd : ℚ}
    (hA' : (Fintype.card V : ℚ) * (D.ε / D.τ + D.E) / (D.Nt + D.E) ≤ A')
    (hD : Dd ≤ δ₀ + 1) (hDA : A' < Dd) (h2A : 2 * A' < Fintype.card V)
    (hCS : D.E ^ 2 * Dd < (Dd - A') * D.Nt ^ 2 * (1 - 2 * A' / Fintype.card V)) :
    HasHamCycle G := by
  have hn0 : (0 : ℚ) < Fintype.card V := by exact_mod_cast (show 0 < Fintype.card V by omega)
  have hα0 : 0 ≤ (Fintype.card V : ℚ) * (D.ε / D.τ + D.E) / (D.Nt + D.E) :=
    div_nonneg (mul_nonneg hn0.le (add_nonneg (div_nonneg D.hε D.hτ.le) D.hE))
      (by linarith [D.hE])
  have hA0 : 0 ≤ A' := le_trans hα0 hA'
  have hNtpos : 0 < D.Nt := by
    rcases hNt.lt_or_eq with h | h
    · exact h
    · exfalso
      rw [← h] at hCS
      have : 0 ≤ D.E ^ 2 * Dd := mul_nonneg (sq_nonneg _) (by linarith)
      simp at hCS
      linarith
  have hNE' : 0 < D.Nt + D.E := by linarith [D.hE]
  refine D.ham_of_bounds hCE hdeg hNtpos hn ⌈A'⌉₊ (Nat.ceil_le.mpr (by linarith)) ?_ ?_
  · intro I hI
    have h2 : (I.card : ℚ) ≤ (⌈A'⌉₊ : ℚ) :=
      ((D.alpha_le hNE' hI).trans hA').trans (Nat.le_ceil A')
    exact_mod_cast h2
  · intro a c hc ha hda h2a
    have hcA : (c : ℚ) < A' := Nat.lt_ceil.mp hc
    exact lemma5_arith D.hE hn0 (Nat.cast_nonneg c) hcA hD hDA h2A hCS hda h2a
      (by exact_mod_cast ha)

end L3Data

end SunLean
