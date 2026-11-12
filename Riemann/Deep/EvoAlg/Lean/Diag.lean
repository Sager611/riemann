import Riemann.Deep.EvoAlg.Lean.Defs
import Mathlib.Combinatorics.SimpleGraph.Connectivity.Connected
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.LinearCombination

/-!
# Diagonal isomorphisms and the "if" direction of Theorem A (NOTE §4 Cor. 4.1 "⇐", §5 "Conversely")
-/

open Finset

namespace EvoAlg

variable {V : Type*} [Fintype V] [DecidableEq V] {𝕜 : Type*} [Field 𝕜]

/-- The diagonal linear equivalence `x ↦ (γ_w x_w)_w` of `𝕜^V` (all `γ_w ≠ 0`). -/
def diagEquiv (γ : V → 𝕜) (hγ : ∀ v, γ v ≠ 0) : (V → 𝕜) ≃ₗ[𝕜] (V → 𝕜) where
  toFun x := fun w => γ w * x w
  invFun x := fun w => (γ w)⁻¹ * x w
  map_add' x y := by funext w; simp [mul_add]
  map_smul' a x := by funext w; simp only [Pi.smul_apply, smul_eq_mul, RingHom.id_apply]; ring
  left_inv x := by funext w; simp [← mul_assoc, inv_mul_cancel₀ (hγ w)]
  right_inv x := by funext w; simp [← mul_assoc, mul_inv_cancel₀ (hγ w)]

omit [DecidableEq V] in
/-- A diagonal map `e_v ↦ γ_v e'_v` is an isomorphism as soon as `γ_v² M' v w = M v w γ_w`. -/
theorem evoIso_of_diag (M M' : V → V → 𝕜) (γ : V → 𝕜) (hγ : ∀ v, γ v ≠ 0)
    (h : ∀ v w, γ v ^ 2 * M' v w = M v w * γ w) : EvoIso M M' := by
  refine ⟨diagEquiv γ hγ, fun x y => ?_⟩
  funext w
  show γ w * evoMul M x y w = evoMul M' (fun u => γ u * x u) (fun u => γ u * y u) w
  simp only [evoMul, Finset.mul_sum]
  refine Finset.sum_congr rfl (fun v _ => ?_)
  linear_combination (-(x v * y v)) * h v w

variable (G : SimpleGraph V) [DecidableRel G.Adj]

omit [DecidableEq V] in
/-- Diagonal criterion for `A(G) → A_RW(G)`: `γ_v²/d_v = γ_w` for every edge `v ∼ w`. -/
theorem aIsoRW_of_diag (γ : V → 𝕜) (hγ : ∀ v, γ v ≠ 0)
    (h : ∀ v w, G.Adj v w → γ v ^ 2 * ((G.degree v : 𝕜))⁻¹ = γ w) : AIsoRW G 𝕜 := by
  refine evoIso_of_diag _ _ γ hγ (fun v w => ?_)
  by_cases hvw : G.Adj v w
  · simp only [adjM, rwM, hvw, ↓reduceIte, one_mul]
    exact h v w hvw
  · simp [adjM, rwM, hvw]

/-- Regular case: `G` `d`-regular with `d ≠ 0` in `𝕜` gives `φ = d · id`. -/
theorem aIsoRW_of_regular {d : ℕ} (hd : (d : 𝕜) ≠ 0) (hreg : G.IsRegularOfDegree d) :
    AIsoRW G 𝕜 := by
  refine aIsoRW_of_diag G (fun _ => (d : 𝕜)) (fun _ => hd) (fun v w _ => ?_)
  rw [hreg.degree_eq v, sq, mul_assoc, mul_inv_cancel₀ hd, mul_one]

/-- Biregular case (NOTE §5, "Conversely"): if `d₁/d₂ = c³` with `c ≠ 0`, then
`φ(e_v) = c²d₂ e'_v` (`v ∈ s`), `φ(e_v) = c d₂ e'_v` (`v ∈ t`) is an isomorphism. -/
theorem aIsoRW_of_biregular {s t : Set V} {d₁ d₂ : ℕ} (hb : IsBiregular G s t d₁ d₂)
    {c : 𝕜} (hc : c ≠ 0) (hcube : (d₁ : 𝕜) / d₂ = c ^ 3) : AIsoRW G 𝕜 := by
  classical
  have hc3 : c ^ 3 ≠ 0 := pow_ne_zero 3 hc
  have hd₂ : (d₂ : 𝕜) ≠ 0 := by
    intro h; rw [h, div_zero] at hcube; exact hc3 hcube.symm
  have hd₁ : (d₁ : 𝕜) = c ^ 3 * d₂ := by rw [← hcube, div_mul_cancel₀ _ hd₂]
  have hd₁0 : (d₁ : 𝕜) ≠ 0 := by rw [hd₁]; exact mul_ne_zero hc3 hd₂
  refine aIsoRW_of_diag G (fun v => if v ∈ s then c ^ 2 * d₂ else c * d₂) ?_ ?_
  · intro v
    by_cases hv : v ∈ s
    · simp only [hv, ↓reduceIte]; exact mul_ne_zero (pow_ne_zero 2 hc) hd₂
    · simp only [hv, ↓reduceIte]; exact mul_ne_zero hc hd₂
  · intro v w hvw
    rcases hb.bip.mem_of_adj hvw with ⟨hv, hw⟩ | ⟨hv, hw⟩
    · have hws : w ∉ s := fun h => Set.disjoint_left.mp hb.bip.disjoint h hw
      simp only [hv, hws, ↓reduceIte]
      rw [hb.deg₁ v hv, hd₁]
      field_simp
    · have hvs : v ∉ s := fun h => Set.disjoint_left.mp hb.bip.disjoint h hv
      simp only [hvs, hw, ↓reduceIte]
      rw [hb.deg₂ v hv]
      field_simp

omit [DecidableEq V] [DecidableRel G.Adj] in
/-- In a connected graph with at least two vertices every vertex has a neighbour. -/
theorem exists_adj_of_connected (hG : G.Connected) (hn : 2 ≤ Fintype.card V) (v : V) :
    ∃ w, G.Adj v w := by
  obtain ⟨u, hu⟩ := Fintype.exists_ne_of_one_lt_card (by omega) v
  obtain ⟨w, hw⟩ := (hG.preconnected v u).nonempty_neighborSet_left (Ne.symm hu)
  exact ⟨w, hw⟩

/-- **Theorem A, "if"** (NOTE §5): for `char 𝕜 = 0`, `G` connected with `n ≥ 2` vertices, regular or
biregular with `d₁/d₂ ∈ (𝕜^×)³` implies `A(G) ≅ A_RW(G)`. -/
theorem theoremA_if [CharZero 𝕜] (hG : G.Connected) (hn : 2 ≤ Fintype.card V)
    (h : IsRegular G ∨ BiregularCube G 𝕜) : AIsoRW G 𝕜 := by
  rcases h with ⟨d, hd⟩ | ⟨s, t, d₁, d₂, hb, c, hc, hcube⟩
  · obtain ⟨v⟩ := hG.nonempty
    obtain ⟨w, hw⟩ := exists_adj_of_connected G hG hn v
    have hpos : 0 < G.degree v := (G.degree_pos_iff_exists_adj v).mpr ⟨w, hw⟩
    rw [hd.degree_eq v] at hpos
    exact aIsoRW_of_regular G (Nat.cast_ne_zero.mpr (by omega)) hd
  · exact aIsoRW_of_biregular G hb hc hcube

end EvoAlg

