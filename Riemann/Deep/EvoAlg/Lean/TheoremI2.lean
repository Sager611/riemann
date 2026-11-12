import Riemann.Deep.EvoAlg.Lean.TwinI
import Riemann.Deep.EvoAlg.Lean.TheoremN
import Mathlib.LinearAlgebra.Matrix.SemiringInverse

/-!
# Theorem I (2) (NOTE §3): the converse

Any `c : V → 𝕜^×` and block matrix `P` (`P a w` = coefficient of `e'_w` in `φ(e_a)`, supported in
the twin class of `a`, i.e. `τ = id`) satisfying (a) `P_SᵀP_S = c_S d_S I` and (b) `P_S𝟙 = c_T𝟙`
(`T ∼ S`) define an isomorphism `A(G) → A_RW(G)`. Bijectivity: `P` has the left inverse
`D⁻¹Pᵀ` by (a); a one-sided inverse of a square matrix over a commutative ring is two-sided (`mul_eq_one_comm`,
with Mathlib's `IsStablyFiniteRing` instance for commutative semirings, `Matrix.SemiringInverse`).
-/

open Finset

namespace EvoAlg

variable {V : Type*} [Fintype V] [DecidableEq V] {𝕜 : Type*} [Field 𝕜]

omit [DecidableEq V] in
theorem evoMul_eq_sum (M : V → V → 𝕜) (x y : V → 𝕜) :
    evoMul M x y = ∑ a, (x a * y a) • M a := by
  funext w
  simp [evoMul, Finset.sum_apply]

theorem map_eq_sum_smul (φ : (V → 𝕜) →ₗ[𝕜] (V → 𝕜)) (x : V → 𝕜) :
    φ x = ∑ a, x a • φ (Pi.single a 1) := by
  have hx : x = ∑ a, x a • (Pi.single a 1 : V → 𝕜) := by
    funext u
    simp [Finset.sum_apply, Pi.single_apply]
  conv_lhs => rw [hx]
  simp [map_sum, map_smul]

/-- A linear map is multiplicative as soon as it is on the natural basis. -/
theorem map_mul_of_basis (M M' : V → V → 𝕜) (φ : (V → 𝕜) →ₗ[𝕜] (V → 𝕜))
    (h : ∀ a b, evoMul M' (φ (Pi.single a 1)) (φ (Pi.single b 1)) = if a = b then φ (M a) else 0)
    (x y : V → 𝕜) : φ (evoMul M x y) = evoMul M' (φ x) (φ y) := by
  rw [map_eq_sum_smul φ x, map_eq_sum_smul φ y, evoMul_sum_smul_right, evoMul_eq_sum, map_sum]
  refine Finset.sum_congr rfl (fun b _ => ?_)
  rw [evoMul_comm, evoMul_sum_smul_right, Finset.sum_eq_single b
    (fun a _ hab => by rw [evoMul_comm, h a b, ite_eq_right hab, smul_zero])
    (fun hb => absurd (Finset.mem_univ b) hb)]
  rw [evoMul_comm, h b b, ite_eq_left rfl, map_smul, smul_smul, mul_comm]

/-- The linear map with matrix `P`: `e_a ↦ Σ_w P a w e_w`. -/
def matMap (P : V → V → 𝕜) : (V → 𝕜) →ₗ[𝕜] (V → 𝕜) where
  toFun x w := ∑ a, x a * P a w
  map_add' x y := by
    funext w
    simp [add_mul, Finset.sum_add_distrib]
  map_smul' r x := by
    funext w
    simp [Finset.mul_sum, mul_assoc]

omit [DecidableEq V] in
theorem matMap_apply (P : V → V → 𝕜) (x : V → 𝕜) (w : V) :
    matMap P x w = ∑ a, x a * P a w := rfl

theorem matMap_single (P : V → V → 𝕜) (a : V) : matMap P (Pi.single a 1) = P a := by
  funext w
  simp [matMap_apply, Pi.single_apply]

omit [DecidableEq V] in
theorem matMap_comp (P R : V → V → 𝕜) (x : V → 𝕜) :
    matMap R (matMap P x) = matMap (fun a w => ∑ v, P a v * R v w) x := by
  funext w
  simp only [matMap_apply, Finset.mul_sum, Finset.sum_mul]
  conv_lhs => rw [Finset.sum_comm]
  exact Finset.sum_congr rfl (fun a _ => Finset.sum_congr rfl (fun v _ => by ring))

theorem matMap_one (P : V → V → 𝕜) (hP : ∀ a w, P a w = if a = w then 1 else 0) (x : V → 𝕜) :
    matMap P x = x := by
  funext w
  simp [matMap_apply, hP]

variable (G : SimpleGraph V) [DecidableRel G.Adj]

omit [DecidableEq V] in
theorem rwM_twin {w a : V} (h : G.neighborFinset w = G.neighborFinset a) (j : V) :
    rwM G 𝕜 w j = rwM G 𝕜 a j := by
  have hadj : G.Adj w j ↔ G.Adj a j := by
    rw [← G.mem_neighborFinset, ← G.mem_neighborFinset, h]
  simp only [rwM, hadj, degree_eq_of_twin G h]

/-- **Theorem I (2)** (NOTE §3, `τ = id`). -/
theorem theoremI_converse (hdeg : ∀ v, (G.degree v : 𝕜) ≠ 0) (c : V → 𝕜) (hc0 : ∀ v, c v ≠ 0)
    (P : V → V → 𝕜)
    (hblock : ∀ a w, P a w ≠ 0 → G.neighborFinset w = G.neighborFinset a)
    (ha : ∀ a b, G.neighborFinset a = G.neighborFinset b →
      ∑ w, P a w * P b w = if a = b then c a * G.degree a else 0)
    (hb : ∀ v u w, G.Adj v u → G.neighborFinset w = G.neighborFinset u →
      ∑ u' ∈ twin G u, P u' w = c v) :
    ∃ φ : (V → 𝕜) ≃ₗ[𝕜] (V → 𝕜), (∀ a, φ (Pi.single a 1) = P a) ∧
      ∀ x y, φ (evoMul (adjM G 𝕜) x y) = evoMul (rwM G 𝕜) (φ x) (φ y) := by
  -- (a) for all pairs
  have hall : ∀ a b, ∑ w, P a w * P b w = if a = b then c a * G.degree a else 0 := by
    intro a b
    by_cases hab : G.neighborFinset a = G.neighborFinset b
    · exact ha a b hab
    · have hne : a ≠ b := fun h => hab (by rw [h])
      rw [ite_eq_right hne]
      refine Finset.sum_eq_zero (fun w _ => ?_)
      by_contra h
      exact hab ((hblock a w (left_ne_zero_of_mul h)).symm.trans
        (hblock b w (right_ne_zero_of_mul h)))
  have hcd : ∀ a, c a * (G.degree a : 𝕜) ≠ 0 := fun a => mul_ne_zero (hc0 a) (hdeg a)
  -- the left inverse R = D⁻¹ Pᵀ
  set R : V → V → 𝕜 := fun v w => P w v / (c w * G.degree w) with hR
  have hPR : ∀ a w, ∑ v, P a v * R v w = if a = w then 1 else 0 := by
    intro a w
    simp only [hR, div_eq_mul_inv, ← mul_assoc]
    rw [← Finset.sum_mul, hall a w]
    by_cases h : a = w
    · subst h
      rw [ite_eq_left rfl, ite_eq_left rfl, mul_inv_cancel₀ (hcd a)]
    · rw [ite_eq_right h, ite_eq_right h, zero_mul]
  have hAB : (Matrix.of fun a v => P a v) * (Matrix.of fun v w => R v w) = 1 := by
    ext a w
    rw [Matrix.mul_apply, Matrix.one_apply]
    simpa using hPR a w
  have hBA : (Matrix.of fun v w => R v w) * (Matrix.of fun a v => P a v) = 1 :=
    mul_eq_one_comm.mp hAB
  have hRP : ∀ a w, ∑ v, R a v * P v w = if a = w then 1 else 0 := by
    intro a w
    have := congrFun (congrFun hBA a) w
    rw [Matrix.mul_apply, Matrix.one_apply] at this
    simpa using this
  let φ : (V → 𝕜) ≃ₗ[𝕜] (V → 𝕜) :=
    { matMap P with
      invFun := matMap R
      left_inv := fun x => by
        show matMap R (matMap P x) = x
        rw [matMap_comp]
        exact matMap_one _ hPR x
      right_inv := fun y => by
        show matMap P (matMap R y) = y
        rw [matMap_comp]
        exact matMap_one _ hRP y }
  have hφ : ∀ x, φ x = matMap P x := fun _ => rfl
  refine ⟨φ, fun a => by rw [hφ, matMap_single], ?_⟩
  -- φ(e_a²) = c_a 1'_{N(a)}
  have hsq : ∀ a j, matMap P (adjM G 𝕜 a) j = if G.Adj a j then c a else 0 := by
    intro a j
    rw [matMap_apply]
    have h1 : ∑ u, adjM G 𝕜 a u * P u j = ∑ u ∈ G.neighborFinset a, P u j := by
      rw [← Finset.sum_filter_add_sum_filter_not Finset.univ (fun u => G.Adj a u)]
      have h2 : ∑ u ∈ Finset.univ.filter (fun u => ¬ G.Adj a u), adjM G 𝕜 a u * P u j = 0 :=
        Finset.sum_eq_zero (fun u hu => by
          rw [Finset.mem_filter] at hu
          simp [adjM, hu.2])
      rw [h2, add_zero]
      have h3 : Finset.univ.filter (fun u => G.Adj a u) = G.neighborFinset a := by
        ext u
        simp
      rw [h3]
      exact Finset.sum_congr rfl (fun u hu => by
        rw [G.mem_neighborFinset] at hu
        simp [adjM, hu])
    rw [h1]
    by_cases haj : G.Adj a j
    · rw [ite_eq_left haj, ← hb a j j haj rfl]
      refine (Finset.sum_subset (fun u hu => ?_) (fun u _ hu => ?_)).symm
      · rw [mem_twin] at hu
        exact (G.mem_neighborFinset a u).mpr (adj_of_twin G hu haj)
      · by_contra h
        exact hu ((mem_twin G).mpr (hblock u j h).symm)
    · rw [ite_eq_right haj]
      refine Finset.sum_eq_zero (fun u hu => ?_)
      by_contra h
      rw [G.mem_neighborFinset] at hu
      exact haj (adj_of_twin G (hblock u j h) hu)
  refine map_mul_of_basis _ _ φ.toLinearMap (fun a b => ?_)
  show evoMul (rwM G 𝕜) (matMap P (Pi.single a 1)) (matMap P (Pi.single b 1))
    = if a = b then matMap P (adjM G 𝕜 a) else 0
  rw [matMap_single, matMap_single]
  funext j
  rw [evoMul_apply]
  by_cases hab : G.neighborFinset a = G.neighborFinset b
  · have hterm : ∀ w, P a w * P b w * rwM G 𝕜 w j = P a w * P b w * rwM G 𝕜 a j := by
      intro w
      by_cases hw : P a w = 0
      · simp [hw]
      · rw [rwM_twin G (hblock a w hw) j]
    rw [Finset.sum_congr rfl (fun w _ => hterm w), ← Finset.sum_mul, ha a b hab]
    by_cases he : a = b
    · subst he
      rw [ite_eq_left rfl, ite_eq_left rfl, hsq a j]
      by_cases haj : G.Adj a j
      · simp only [rwM, haj, ↓reduceIte]
        rw [mul_assoc, mul_inv_cancel₀ (hdeg a), mul_one]
      · simp [rwM, haj]
    · rw [ite_eq_right he, ite_eq_right he, zero_mul, Pi.zero_apply]
  · have hne : a ≠ b := fun h => hab (by rw [h])
    rw [ite_eq_right hne, Pi.zero_apply]
    refine Finset.sum_eq_zero (fun w _ => ?_)
    by_cases h : P a w * P b w = 0
    · rw [h, zero_mul]
    · exact absurd ((hblock a w (left_ne_zero_of_mul h)).symm.trans
        (hblock b w (right_ne_zero_of_mul h))) hab

end EvoAlg

