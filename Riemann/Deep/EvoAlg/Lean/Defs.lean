import Mathlib.Combinatorics.SimpleGraph.Bipartite
import Mathlib.Algebra.Module.Equiv.Defs
import Mathlib.Algebra.Module.Pi
import Mathlib.Algebra.Field.Basic
import Mathlib.Algebra.BigOperators.Pi
import Mathlib.Tactic.Ring

/-!
# Evolution algebras A(G) and A_RW(G): definitions (NOTE §0–§1)

Structure-constant model. An evolution algebra over `𝕜` with natural basis `(e_v)_{v ∈ V}` and
`e_v² = Σ_w M v w e_w` (and `e_v e_w = 0` for `v ≠ w`) is, in coordinates, `𝕜^V = (V → 𝕜)` with
the commutative bilinear product `(x ⋆ y)_w = Σ_v x_v y_v M v w` and `e_v = Pi.single v 1`
(`evoMul`, lemmas `evoMul_single_single`, `evoMul_single_self`, `evoMul_comm`, bilinearity).
An algebra isomorphism is a `𝕜`-linear equivalence preserving the product (`EvoIso`).
-/

open Finset

namespace EvoAlg

variable {V : Type*} [Fintype V] [DecidableEq V]

section alg

variable {𝕜 : Type*} [Field 𝕜]

/-- The evolution product on `𝕜^V` with structure matrix `M`. -/
def evoMul (M : V → V → 𝕜) (x y : V → 𝕜) : V → 𝕜 := fun w => ∑ v, x v * y v * M v w

/-- `EvoIso M M'`: the evolution algebras with structure matrices `M`, `M'` are isomorphic, i.e.
there is a `𝕜`-linear equivalence `φ` of `𝕜^V` with `φ (x ⋆ y) = φ x ⋆' φ y`. -/
def EvoIso (M M' : V → V → 𝕜) : Prop :=
  ∃ φ : (V → 𝕜) ≃ₗ[𝕜] (V → 𝕜), ∀ x y, φ (evoMul M x y) = evoMul M' (φ x) (φ y)

section basic

variable (M : V → V → 𝕜)

omit [DecidableEq V] in
theorem evoMul_apply (x y : V → 𝕜) (w : V) : evoMul M x y w = ∑ v, x v * y v * M v w := rfl

theorem evoMul_comm (x y : V → 𝕜) : evoMul M x y = evoMul M y x := by
  funext w
  simp only [evoMul]
  exact Finset.sum_congr rfl (fun v _ => by ring)

theorem evoMul_add_left (x x' y : V → 𝕜) :
    evoMul M (x + x') y = evoMul M x y + evoMul M x' y := by
  funext w
  simp only [evoMul, Pi.add_apply, ← Finset.sum_add_distrib]
  exact Finset.sum_congr rfl (fun v _ => by ring)

theorem evoMul_smul_left (a : 𝕜) (x y : V → 𝕜) : evoMul M (a • x) y = a • evoMul M x y := by
  funext w
  simp only [evoMul, Pi.smul_apply, smul_eq_mul, Finset.mul_sum]
  exact Finset.sum_congr rfl (fun v _ => by ring)

theorem evoMul_single_left (v : V) (y : V → 𝕜) : evoMul M (Pi.single v 1) y = y v • M v := by
  funext w
  simp [evoMul, Pi.single_apply]

theorem evoMul_single_right (x : V → 𝕜) (v : V) : evoMul M x (Pi.single v 1) = x v • M v := by
  rw [evoMul_comm, evoMul_single_left]

/-- natural basis: `e_v e_u = 0` for `v ≠ u` and `e_v e_v = M v`. -/
theorem evoMul_single_single (v u : V) :
    evoMul M (Pi.single v 1) (Pi.single u 1) = if v = u then M v else 0 := by
  rw [evoMul_single_left]
  by_cases h : v = u
  · subst h; simp
  · simp [h]

/-- `e_v² = Σ_w M v w e_w`. -/
theorem evoMul_single_self (v : V) :
    evoMul M (Pi.single v 1) (Pi.single v 1) = ∑ w, M v w • (Pi.single w 1 : V → 𝕜) := by
  rw [evoMul_single_single, ite_eq_left rfl]
  funext u
  simp [Finset.sum_apply, Pi.single_apply]

end basic

end alg

section graph

variable (G : SimpleGraph V) [DecidableRel G.Adj] (𝕜 : Type*) [Field 𝕜]

/-- Structure matrix of `A(G)`: `e_v² = Σ_{w ∼ v} e_w`. -/
def adjM : V → V → 𝕜 := fun v w => if G.Adj v w then 1 else 0

/-- Structure matrix of `A_RW(G)`: `e'_v² = (deg v)⁻¹ Σ_{w ∼ v} e'_w`. -/
def rwM : V → V → 𝕜 := fun v w => if G.Adj v w then ((G.degree v : 𝕜))⁻¹ else 0

/-- `A(G) ≅ A_RW(G)` over `𝕜`. -/
def AIsoRW : Prop := EvoIso (adjM G 𝕜) (rwM G 𝕜)

/-- `G` is regular. -/
def IsRegular : Prop := ∃ d : ℕ, G.IsRegularOfDegree d

/-- `G` is biregular with sides `s`, `t` and degrees `d₁`, `d₂` (NOTE §0; CRR l. 159): `G` is
bipartite with sides `s`, `t` (a partition of `V`), every vertex of `s` has degree `d₁` and every
vertex of `t` has degree `d₂`. -/
structure IsBiregular (s t : Set V) (d₁ d₂ : ℕ) : Prop where
  bip : G.IsBipartiteWith s t
  cover : ∀ v, v ∈ s ∨ v ∈ t
  deg₁ : ∀ v ∈ s, G.degree v = d₁
  deg₂ : ∀ v ∈ t, G.degree v = d₂

/-- Condition 2 of Theorem A: `G` is biregular and `d₁/d₂ ∈ (𝕜^×)³`. -/
def BiregularCube : Prop :=
  ∃ (s t : Set V) (d₁ d₂ : ℕ), IsBiregular G s t d₁ d₂ ∧ ∃ c : 𝕜, c ≠ 0 ∧ (d₁ : 𝕜) / d₂ = c ^ 3

variable {𝕜}

/-- In `A(G)`: `e_v² = Σ_{w ∈ N(v)} e_w`. -/
theorem adj_sq (v : V) : evoMul (adjM G 𝕜) (Pi.single v 1) (Pi.single v 1)
    = ∑ w ∈ G.neighborFinset v, (Pi.single w 1 : V → 𝕜) := by
  rw [evoMul_single_single, ite_eq_left rfl]
  funext u
  simp [adjM, Finset.sum_apply, Pi.single_apply]

/-- In `A_RW(G)`: `e'_v² = (deg v)⁻¹ • Σ_{w ∈ N(v)} e'_w`. -/
theorem rw_sq (v : V) : evoMul (rwM G 𝕜) (Pi.single v 1) (Pi.single v 1)
    = ((G.degree v : 𝕜))⁻¹ • ∑ w ∈ G.neighborFinset v, (Pi.single w 1 : V → 𝕜) := by
  rw [evoMul_single_single, ite_eq_left rfl]
  funext u
  simp [rwM, Finset.sum_apply, Pi.single_apply]

end graph

end EvoAlg

