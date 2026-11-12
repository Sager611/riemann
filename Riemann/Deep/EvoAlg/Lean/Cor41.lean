import Riemann.Deep.EvoAlg.Lean.TwinI
import Riemann.Deep.EvoAlg.Lean.Lemma2
import Riemann.Deep.EvoAlg.Lean.Diag

/-!
# Theorem I (NOTE §3) and Corollary 4.1 (NOTE §4): the scaling system

Vertex-level form. A function `c : V → 𝕜` satisfying (i) is constant on twin classes, so it is
a function on `V(H)`; (i) for vertices is (i) for classes (two classes have a common
`H`-neighbour iff two of their vertices have a common neighbour), and (ii) for an edge `u ∼ v`
is (ii) for the edge `S ∼ T` of `H` (`u ∈ S`, `v ∈ T`, `d_S = deg u`).
-/

open Finset

namespace EvoAlg

variable {V : Type*} [Fintype V] [DecidableEq V] {𝕜 : Type*} [Field 𝕜]

section iso

variable {G : SimpleGraph V} [DecidableRel G.Adj] (φ : (V → 𝕜) ≃ₗ[𝕜] (V → 𝕜))
  (hφ : ∀ x y, φ (evoMul (adjM G 𝕜) x y) = evoMul (rwM G 𝕜) (φ x) (φ y))
  (hnb : ∀ v, ∃ w, G.Adj v w) (hdeg : ∀ v, (G.degree v : 𝕜) ≠ 0)
  (c : V → 𝕜) (hc0 : ∀ a, c a ≠ 0)
  (hc : ∀ a i, φ (Pi.single a 1) i ≠ 0 →
    φ (ind (G.neighborFinset a)) = c a • ind (G.neighborFinset i))

include hφ hnb hdeg hc0 hc

/-- (ii) from Lemma 2: for `u ∼ v`, `c_v² = c_u d_u` (block `P_S`, `S = twin u`, rows `τS`, with
`α = c_S d_S` by (a), `γ = c_T` by (b), and `|S| = |τS| ≠ 0`). -/
theorem scaling_ii (hk : ∀ u, (#(twin G u) : 𝕜) ≠ 0) {u v : V} (huv : G.Adj u v) :
    c v ^ 2 = c u * G.degree u := by
  obtain ⟨w₀, hw₀⟩ := col_ne_zero φ u
  have hA : ∀ a ∈ twin G u, ∀ b ∈ twin G u,
      ∑ w ∈ twin G w₀, φ (Pi.single a 1) w * φ (Pi.single b 1) w
        = if a = b then c u * G.degree u else 0 := by
    intro a ha b _
    have haN := (mem_twin G).mp ha
    rw [← c_twin φ hφ hnb hdeg c hc0 hc haN, ← degree_eq_of_twin G haN,
      ← blockA φ hφ hnb hdeg c hc0 hc a b]
    refine Finset.sum_subset (Finset.subset_univ _) (fun w _ hw => ?_)
    by_contra hne
    exact hw ((mem_twin G).mpr (L1b φ hφ hnb hdeg (left_ne_zero_of_mul hne) hw₀ haN))
  have hB : ∀ w ∈ twin G w₀, ∑ a ∈ twin G u, φ (Pi.single a 1) w = c v :=
    fun w hw => blockB φ hφ hnb hdeg c hc0 hc huv.symm hw₀ ((mem_twin G).mp hw)
  have hL := lemma2_blocks (twin G u) (twin G w₀) (fun a w => φ (Pi.single a 1) w) _ _ hA hB
  rw [← card_twin_eq φ hφ hnb hdeg hw₀] at hL
  exact (mul_left_cancel₀ (hk u) hL).symm

end iso

variable (G : SimpleGraph V) [DecidableRel G.Adj]

/-- **Theorem I (1)** (NOTE §3), vertex-level form, over any field: for an isomorphism
`φ : A(G) → A_RW(G)` with matrix `P a w = (φ e_a)_w` there is `c : V → 𝕜^×` with
* `φ(e_a)² = c_a 1'_{N(i)}` for `i` in the support of `φ(e_a)` (definition of `c_S`);
* every column and every row of `P` is non-zero, and for `P a i ≠ 0 ≠ P b j`:
  `a ≡ b ↔ i ≡ j` (the class map `τ` is a well-defined bijection) and `a ∼ b ↔ i ∼ j` (`τ ∈ Aut H`);
* `P a i ≠ 0` gives `|twin a| = |twin i|` and `d_a = d_i` in `𝕜` (`k_{τS} = k_S`, `d_{τS} = d_S`);
* (a) `P_SᵀP_S = c_S d_S I`: `Σ_w P a w P b w = δ_{ab} c_a d_a` for twins `a, b`;
* (b) `P_S𝟙 = c_T𝟙` for `T ∼ S`: `Σ_{u' ≡ u} P u' w = c_v` for `v ∼ u` and `w ∈ τ(twin u)`. -/
theorem theoremI (hnb : ∀ v, ∃ w, G.Adj v w) (hdeg : ∀ v, (G.degree v : 𝕜) ≠ 0)
    (φ : (V → 𝕜) ≃ₗ[𝕜] (V → 𝕜))
    (hφ : ∀ x y, φ (evoMul (adjM G 𝕜) x y) = evoMul (rwM G 𝕜) (φ x) (φ y)) :
    ∃ c : V → 𝕜, (∀ a, c a ≠ 0) ∧
      (∀ a i, φ (Pi.single a 1) i ≠ 0 → evoMul (rwM G 𝕜) (φ (Pi.single a 1)) (φ (Pi.single a 1))
        = c a • ind (G.neighborFinset i)) ∧
      (∀ a, ∃ i, φ (Pi.single a 1) i ≠ 0) ∧ (∀ i, ∃ a, φ (Pi.single a 1) i ≠ 0) ∧
      (∀ a b i j, φ (Pi.single a 1) i ≠ 0 → φ (Pi.single b 1) j ≠ 0 →
        (G.neighborFinset a = G.neighborFinset b ↔ G.neighborFinset i = G.neighborFinset j) ∧
        (G.Adj a b ↔ G.Adj i j)) ∧
      (∀ a i, φ (Pi.single a 1) i ≠ 0 →
        (#(twin G a) : 𝕜) = #(twin G i) ∧ (G.degree a : 𝕜) = G.degree i) ∧
      (∀ a b, G.neighborFinset a = G.neighborFinset b →
        ∑ w, φ (Pi.single a 1) w * φ (Pi.single b 1) w = if a = b then c a * G.degree a else 0) ∧
      (∀ v u w₀ w, G.Adj v u → φ (Pi.single u 1) w₀ ≠ 0 →
        G.neighborFinset w = G.neighborFinset w₀ → ∑ u' ∈ twin G u, φ (Pi.single u' 1) w = c v) := by
  choose c hc0 hc using exists_c φ hφ hnb hdeg
  refine ⟨c, hc0, fun a i hi => ?_, col_ne_zero φ, row_ne_zero φ,
    fun a b i j hi hj => ⟨TW φ hφ hnb hdeg hi hj, adj_iff φ hφ hnb hdeg c hc0 hc hi hj⟩,
    fun a i hi => ⟨card_twin_eq φ hφ hnb hdeg hi, deg_eq φ hφ hnb hdeg c hc0 hc hi⟩,
    fun a b _ => blockA φ hφ hnb hdeg c hc0 hc a b,
    fun v u w₀ w hvu hw₀ hw => blockB φ hφ hnb hdeg c hc0 hc hvu hw₀ hw⟩
  rw [sq_expand φ hφ a a, ite_eq_left rfl, adjM_eq_ind, hc a i hi]

/-- **Corollary 4.1** (NOTE §4), for every twin-class size non-zero in `𝕜` (always true in
characteristic 0): `A(G) ≅ A_RW(G)` iff the scaling system has a solution
`c : V → 𝕜^×`: (i) `c_v = c_{v'}` when `v, v'` have a common neighbour `u`, and
(ii) `c_v² = c_u d_u` for every edge `u ∼ v`. -/
theorem cor41 (hnb : ∀ v, ∃ w, G.Adj v w) (hdeg : ∀ v, (G.degree v : 𝕜) ≠ 0)
    (hk : ∀ u, (#(twin G u) : 𝕜) ≠ 0) :
    AIsoRW G 𝕜 ↔ ∃ c : V → 𝕜, (∀ v, c v ≠ 0) ∧
      (∀ u v v', G.Adj u v → G.Adj u v' → c v = c v') ∧
      (∀ u v, G.Adj u v → c v ^ 2 = c u * G.degree u) := by
  constructor
  · intro h
    obtain ⟨φ, hφ⟩ := h
    choose c hc0 hc using exists_c φ hφ hnb hdeg
    refine ⟨c, hc0, fun u v v' huv huv' => ?_,
      fun u v huv => scaling_ii φ hφ hnb hdeg c hc0 hc hk huv⟩
    obtain ⟨w₀, hw₀⟩ := col_ne_zero φ u
    exact (blockB φ hφ hnb hdeg c hc0 hc huv.symm hw₀ rfl).symm.trans
      (blockB φ hφ hnb hdeg c hc0 hc huv'.symm hw₀ rfl)
  · rintro ⟨c, hc0, hi, hii⟩
    choose nb hnb' using hnb
    refine aIsoRW_of_diag G (fun u => c (nb u)) (fun u => hc0 _) (fun u w huw => ?_)
    show c (nb u) ^ 2 * ((G.degree u : 𝕜))⁻¹ = c (nb w)
    rw [hii u (nb u) (hnb' u), mul_assoc, mul_inv_cancel₀ (hdeg u), mul_one]
    exact hi w u (nb w) huw.symm (hnb' w)

end EvoAlg

