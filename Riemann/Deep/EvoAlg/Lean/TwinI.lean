import Riemann.Deep.EvoAlg.Lean.Basic

/-!
# Theorem I (NOTE §3), vertex-level lemmas

Setting: `φ : A(G) → A_RW(G)` an algebra isomorphism, `P a w = (φ e_a)_w`, every vertex has a
neighbour (so both algebras are non-degenerate) and every degree is non-zero in `𝕜`.
Twin classes: `twin G u = {u' : N(u') = N(u)}`. Lemma 1 (Basic.lean) gives Theorem N's
consequence: the support of each column `φ(e_a)` lies in one twin class (`L1a`); `TW` says that
the induced class map `τ` is well defined and injective; `adj_iff` says `τ ∈ Aut(H)`.
-/

open Finset

namespace EvoAlg

variable {V : Type*} [Fintype V] [DecidableEq V] {𝕜 : Type*} [Field 𝕜]
variable (G : SimpleGraph V) [DecidableRel G.Adj]

/-- The twin class of `u`: the vertices with the same neighbourhood as `u`. -/
def twin (u : V) : Finset V := univ.filter (fun u' => G.neighborFinset u' = G.neighborFinset u)

theorem mem_twin {u u' : V} : u' ∈ twin G u ↔ G.neighborFinset u' = G.neighborFinset u := by
  simp [twin]

theorem adjM_eq_ind (a : V) : adjM G 𝕜 a = ind (G.neighborFinset a) := by
  funext w
  simp [adjM, ind]

theorem rwM_eq_ind (i : V) :
    rwM G 𝕜 i = ((G.degree i : 𝕜))⁻¹ • (ind (G.neighborFinset i) : V → 𝕜) := by
  funext w
  by_cases h : G.Adj i w <;> simp [rwM, ind, h]

omit [Fintype V] [DecidableEq V] in
/-- `N(i)` is a union of twin classes. -/
theorem adj_of_twin [Fintype V] {i w w' : V} (hw : G.neighborFinset w = G.neighborFinset w')
    (h : G.Adj i w') : G.Adj i w := by
  have : i ∈ G.neighborFinset w' := (G.mem_neighborFinset w' i).mpr h.symm
  rw [← hw] at this
  exact ((G.mem_neighborFinset w i).mp this).symm

omit [DecidableEq V] in
theorem degree_eq_of_twin {u u' : V} (h : G.neighborFinset u = G.neighborFinset u') :
    G.degree u = G.degree u' := by
  show #(G.neighborFinset u) = #(G.neighborFinset u')
  rw [h]

section iso

variable {G} (φ : (V → 𝕜) ≃ₗ[𝕜] (V → 𝕜))
  (hφ : ∀ x y, φ (evoMul (adjM G 𝕜) x y) = evoMul (rwM G 𝕜) (φ x) (φ y))
  (hnb : ∀ v, ∃ w, G.Adj v w) (hdeg : ∀ v, (G.degree v : 𝕜) ≠ 0)

include hφ hnb hdeg

/-- Lemma 1 for `A(G) → A_RW(G)`: if `(φ e_a)_i ≠ 0` then `φ(1_{N(a)}) = μ 1_{N(i)}`, `μ ≠ 0`. -/
theorem K1 (a i : V) (h : φ (Pi.single a 1) i ≠ 0) :
    ∃ μ : 𝕜, μ ≠ 0 ∧ φ (ind (G.neighborFinset a)) = μ • ind (G.neighborFinset i) := by
  have hL := lemma1 φ hφ a i h
  rw [rwM_eq_ind, adjM_eq_ind] at hL
  generalize (φ (Pi.single a 1) i)⁻¹ * φ.symm (Pi.single i 1) a = κ at hL
  obtain ⟨j, hj⟩ := hnb i
  have hκ : κ ≠ 0 := by
    intro hκ0
    have := congrFun hL j
    rw [hκ0, zero_smul] at this
    simp [ind, hj, hdeg i] at this
  refine ⟨κ⁻¹ * ((G.degree i : 𝕜))⁻¹, mul_ne_zero (inv_ne_zero hκ) (inv_ne_zero (hdeg i)), ?_⟩
  rw [mul_smul, hL, smul_smul, inv_mul_cancel₀ hκ, one_smul]

/-- The support of a column lies in one twin class. -/
theorem L1a {a i j : V} (hi : φ (Pi.single a 1) i ≠ 0) (hj : φ (Pi.single a 1) j ≠ 0) :
    G.neighborFinset i = G.neighborFinset j := by
  obtain ⟨μ, hμ, h1⟩ := K1 φ hφ hnb hdeg a i hi
  obtain ⟨μ', hμ', h2⟩ := K1 φ hφ hnb hdeg a j hj
  exact eq_of_smul_ind hμ hμ' (h1.symm.trans h2)

/-- Twin columns are supported in the same twin class (`τ` is well defined on classes). -/
theorem L1b {a b i j : V} (hi : φ (Pi.single a 1) i ≠ 0) (hj : φ (Pi.single b 1) j ≠ 0)
    (hab : G.neighborFinset a = G.neighborFinset b) : G.neighborFinset i = G.neighborFinset j := by
  obtain ⟨μ, hμ, h1⟩ := K1 φ hφ hnb hdeg a i hi
  obtain ⟨μ', hμ', h2⟩ := K1 φ hφ hnb hdeg b j hj
  rw [hab] at h1
  exact eq_of_smul_ind hμ hμ' (h1.symm.trans h2)

/-- Columns supported in the same twin class are twins (`τ` is injective on classes). -/
theorem L1d {a b i j : V} (hi : φ (Pi.single a 1) i ≠ 0) (hj : φ (Pi.single b 1) j ≠ 0)
    (hij : G.neighborFinset i = G.neighborFinset j) : G.neighborFinset a = G.neighborFinset b := by
  obtain ⟨μ, hμ, h1⟩ := K1 φ hφ hnb hdeg a i hi
  obtain ⟨μ', hμ', h2⟩ := K1 φ hφ hnb hdeg b j hj
  rw [← hij] at h2
  have h3 : φ (μ' • ind (G.neighborFinset a)) = φ (μ • ind (G.neighborFinset b)) := by
    rw [map_smul, map_smul, h1, h2, smul_smul, smul_smul, mul_comm]
  exact eq_of_smul_ind hμ' hμ (φ.injective h3)

theorem TW {a b i j : V} (hi : φ (Pi.single a 1) i ≠ 0) (hj : φ (Pi.single b 1) j ≠ 0) :
    G.neighborFinset a = G.neighborFinset b ↔ G.neighborFinset i = G.neighborFinset j :=
  ⟨L1b φ hφ hnb hdeg hi hj, L1d φ hφ hnb hdeg hi hj⟩

/-- Definition of `c` (NOTE §3, Step 2): `φ(e_a)² = φ(1_{N(a)}) = c_a 1'_{N(i)}` for every `i`
in the support of `φ(e_a)`. -/
theorem exists_c (a : V) : ∃ μ : 𝕜, μ ≠ 0 ∧ ∀ i, φ (Pi.single a 1) i ≠ 0 →
    φ (ind (G.neighborFinset a)) = μ • ind (G.neighborFinset i) := by
  obtain ⟨i₀, hi₀⟩ := col_ne_zero φ a
  obtain ⟨μ, hμ, h⟩ := K1 φ hφ hnb hdeg a i₀ hi₀
  refine ⟨μ, hμ, fun i hi => ?_⟩
  rw [h, L1a φ hφ hnb hdeg hi₀ hi]

/-- Twin classes of `A_RW(G)` have the same size as their preimage classes (trace trick). -/
theorem card_twin_eq {a i : V} (hi : φ (Pi.single a 1) i ≠ 0) :
    (#(twin G a) : 𝕜) = #(twin G i) := by
  refine card_eq_of_support φ _ _ (fun b j hj => ?_)
  rw [mem_twin, mem_twin]
  exact TW φ hφ hnb hdeg hj hi

section c

variable (c : V → 𝕜) (hc0 : ∀ a, c a ≠ 0)
  (hc : ∀ a i, φ (Pi.single a 1) i ≠ 0 →
    φ (ind (G.neighborFinset a)) = c a • ind (G.neighborFinset i))

include hc0 hc

/-- `c` is constant on twin classes. -/
theorem c_twin {a b : V} (hab : G.neighborFinset a = G.neighborFinset b) : c a = c b := by
  obtain ⟨i, hi⟩ := col_ne_zero φ a
  obtain ⟨j, hj⟩ := col_ne_zero φ b
  have hij := L1b φ hφ hnb hdeg hi hj hab
  have h1 := hc a i hi
  have h2 := hc b j hj
  rw [hab, h2, ← hij] at h1
  obtain ⟨k, hk⟩ := hnb i
  exact (smul_eq_of_smul_ind ((G.mem_neighborFinset i k).mpr hk) h1).symm

/-- `τ` maps `H`-edges to `H`-edges (NOTE §3, Step 3, first part). -/
theorem adj_fwd {v u w i : V} (hvu : G.Adj v u) (hw : φ (Pi.single u 1) w ≠ 0)
    (hi : φ (Pi.single v 1) i ≠ 0) : G.Adj i w := by
  have hne : φ (ind (twin G u)) ≠ 0 := by
    intro h0
    exact ind_ne_zero (𝕜 := 𝕜) ((mem_twin G).mpr rfl) (φ.injective (h0.trans (map_zero φ).symm))
  obtain ⟨w', hw'⟩ : ∃ w', φ (ind (twin G u)) w' ≠ 0 := by
    by_contra h
    exact hne (funext fun x => by
      by_contra hx
      exact h ⟨x, hx⟩)
  rw [map_ind_apply] at hw'
  obtain ⟨u'', hu''S, hu''⟩ := Finset.exists_ne_zero_of_sum_ne_zero hw'
  have hu''N : G.neighborFinset u'' = G.neighborFinset u := (mem_twin G).mp hu''S
  have hsub : twin G u ⊆ G.neighborFinset v := by
    intro x hx
    rw [mem_twin] at hx
    exact (G.mem_neighborFinset v x).mpr (adj_of_twin G hx hvu)
  have hsum : ∑ a ∈ G.neighborFinset v, φ (Pi.single a 1) w'
      = ∑ a ∈ twin G u, φ (Pi.single a 1) w' := by
    refine (Finset.sum_subset hsub (fun x _ hxS => ?_)).symm
    by_contra hx
    exact hxS ((mem_twin G).mpr ((L1d φ hφ hnb hdeg hx hu'' rfl).trans hu''N))
  have hval : φ (ind (G.neighborFinset v)) w' ≠ 0 := by
    rw [map_ind_apply, hsum]
    exact hw'
  rw [hc v i hi] at hval
  have hw'i : G.Adj i w' := by
    by_contra hadj
    exact hval (by simp [ind, hadj])
  have hww' : G.neighborFinset w' = G.neighborFinset w := L1b φ hφ hnb hdeg hu'' hw hu''N
  exact adj_of_twin G hww'.symm hw'i

/-- `τ⁻¹` maps `H`-edges to `H`-edges (NOTE §3, Step 3, second part). -/
theorem adj_bwd {a b i j : V} (hij : G.Adj i j) (hi : φ (Pi.single a 1) i ≠ 0)
    (hj : φ (Pi.single b 1) j ≠ 0) : G.Adj a b := by
  have h := congrFun (hc a i hi) j
  rw [map_ind_apply] at h
  have hne : ∑ x ∈ G.neighborFinset a, φ (Pi.single x 1) j ≠ 0 := by
    rw [h]
    simp [ind, hij, hc0 a]
  obtain ⟨a', ha'N, ha'⟩ := Finset.exists_ne_zero_of_sum_ne_zero hne
  have htw : G.neighborFinset a' = G.neighborFinset b := L1d φ hφ hnb hdeg ha' hj rfl
  exact adj_of_twin G htw.symm ((G.mem_neighborFinset a a').mp ha'N)

/-- `τ ∈ Aut(H)`. -/
theorem adj_iff {a b i j : V} (hi : φ (Pi.single a 1) i ≠ 0) (hj : φ (Pi.single b 1) j ≠ 0) :
    G.Adj a b ↔ G.Adj i j :=
  ⟨fun h => adj_fwd φ hφ hnb hdeg c hc0 hc h hj hi, fun h => adj_bwd φ hφ hnb hdeg c hc0 hc h hi hj⟩

/-- `d_{τS} = d_S` (trace trick on `N(a) × N(i)`). -/
theorem deg_eq {a i : V} (hi : φ (Pi.single a 1) i ≠ 0) : (G.degree a : 𝕜) = G.degree i := by
  show (#(G.neighborFinset a) : 𝕜) = #(G.neighborFinset i)
  refine card_eq_of_support φ _ _ (fun b j hj => ?_)
  rw [G.mem_neighborFinset, G.mem_neighborFinset]
  exact adj_iff φ hφ hnb hdeg c hc0 hc hi hj

/-- (a) `P_SᵀP_S = c_S d_S I` (in fact `Σ_w (φ e_a)_w (φ e_b)_w = δ_{ab} c_a d_a` for all `a, b`). -/
theorem blockA (a b : V) : ∑ w, φ (Pi.single a 1) w * φ (Pi.single b 1) w
    = if a = b then c a * G.degree a else 0 := by
  obtain ⟨i, hi⟩ := col_ne_zero φ a
  obtain ⟨j, hj⟩ := hnb i
  have hsq := congrFun (sq_expand φ hφ a b) j
  rw [evoMul_apply] at hsq
  have hterm : ∀ w, φ (Pi.single a 1) w * φ (Pi.single b 1) w * rwM G 𝕜 w j
      = (φ (Pi.single a 1) w * φ (Pi.single b 1) w) * ((G.degree i : 𝕜))⁻¹ := by
    intro w
    by_cases hw : φ (Pi.single a 1) w = 0
    · simp [hw]
    · have hN : G.neighborFinset w = G.neighborFinset i := L1a φ hφ hnb hdeg hw hi
      have hadj : G.Adj w j := by
        have : j ∈ G.neighborFinset w := by
          rw [hN]
          exact (G.mem_neighborFinset i j).mpr hj
        exact (G.mem_neighborFinset w j).mp this
      simp only [rwM, hadj, ↓reduceIte, degree_eq_of_twin G hN]
  rw [Finset.sum_congr rfl (fun w _ => hterm w), ← Finset.sum_mul] at hsq
  rw [deg_eq φ hφ hnb hdeg c hc0 hc hi]
  have hdi := hdeg i
  by_cases hab : a = b
  · subst hab
    rw [ite_eq_left rfl] at hsq ⊢
    rw [adjM_eq_ind, hc a i hi] at hsq
    simp only [ind, Pi.smul_apply, smul_eq_mul, (G.mem_neighborFinset i j).mpr hj, ↓reduceIte,
      mul_one] at hsq
    rw [← hsq, mul_assoc, inv_mul_cancel₀ hdi, mul_one]
  · rw [ite_eq_right hab] at hsq ⊢
    rw [Pi.zero_apply] at hsq
    exact (mul_eq_zero.mp hsq).resolve_right (inv_ne_zero hdi)

/-- (b) `P_S𝟙 = c_T𝟙` for `T ∼ S`: for `v ∼ u` and `w` in the class `τ(twin u)`,
`Σ_{u' ∈ twin u} (φ e_{u'})_w = c_v`. -/
theorem blockB {v u w₀ w : V} (hvu : G.Adj v u) (hw₀ : φ (Pi.single u 1) w₀ ≠ 0)
    (hw : G.neighborFinset w = G.neighborFinset w₀) :
    ∑ u' ∈ twin G u, φ (Pi.single u' 1) w = c v := by
  obtain ⟨i, hi⟩ := col_ne_zero φ v
  have hiw₀ : G.Adj i w₀ := adj_fwd φ hφ hnb hdeg c hc0 hc hvu hw₀ hi
  have hiw : G.Adj i w := adj_of_twin G hw hiw₀
  have h := congrFun (hc v i hi) w
  rw [map_ind_apply] at h
  simp only [Pi.smul_apply, smul_eq_mul, ind, (G.mem_neighborFinset i w).mpr hiw, ↓reduceIte,
    mul_one] at h
  rw [← h]
  refine Finset.sum_subset (fun x hx => ?_) (fun x _ hxS => ?_)
  · rw [mem_twin] at hx
    exact (G.mem_neighborFinset v x).mpr (adj_of_twin G hx hvu)
  · by_contra hx
    exact hxS ((mem_twin G).mpr (L1d φ hφ hnb hdeg hx hw₀ hw))

end c

end iso

end EvoAlg

