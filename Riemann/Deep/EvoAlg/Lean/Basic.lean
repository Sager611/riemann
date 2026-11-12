import Riemann.Deep.EvoAlg.Lean.Defs

/-!
# Linear-algebra facts used in NOTE §2–§3

Indicator vectors `1_S`, the coordinates of a linear equivalence `φ` of `𝕜^V` (`P a w = (φ e_a)_w`,
`Q w a = (φ⁻¹ e_w)_a`, so `QP = PQ = I`), the trace trick that replaces the counting of Theorem N
(`|A| = |B|` in `𝕜` when the support of `P` respects `A × B`), and Lemma 1 (BCS22 Thm 2.4(i) "⇒").
-/

open Finset

namespace EvoAlg

variable {V : Type*} [Fintype V] [DecidableEq V] {𝕜 : Type*} [Field 𝕜]

/-- The indicator vector `1_S = Σ_{u ∈ S} e_u`. -/
def ind (S : Finset V) : V → 𝕜 := fun w => if w ∈ S then 1 else 0

theorem ind_eq_sum (S : Finset V) : (ind S : V → 𝕜) = ∑ u ∈ S, (Pi.single u 1 : V → 𝕜) := by
  funext w
  simp [ind, Finset.sum_apply, Pi.single_apply]

theorem ind_ne_zero {S : Finset V} {u : V} (hu : u ∈ S) : (ind S : V → 𝕜) ≠ 0 := by
  intro h
  have := congrFun h u
  simp [ind, hu] at this

theorem eq_of_smul_ind {μ μ' : 𝕜} (hμ : μ ≠ 0) (hμ' : μ' ≠ 0) {A B : Finset V}
    (h : μ • (ind A : V → 𝕜) = μ' • ind B) : A = B := by
  ext w
  have hw := congrFun h w
  simp only [ind, Pi.smul_apply, smul_eq_mul] at hw
  by_cases hA : w ∈ A <;> by_cases hB : w ∈ B <;> simp_all

theorem smul_eq_of_smul_ind {μ μ' : 𝕜} {A : Finset V} {u : V} (hu : u ∈ A)
    (h : μ • (ind A : V → 𝕜) = μ' • ind A) : μ = μ' := by
  have hw := congrFun h u
  simpa [ind, hu] using hw

theorem sum_indicator (A : Finset V) : ∑ a, (if a ∈ A then (1 : 𝕜) else 0) = #A := by
  rw [Finset.sum_ite_mem, Finset.univ_inter]
  simp

section linear

variable (φ : (V → 𝕜) ≃ₗ[𝕜] (V → 𝕜))

/-- Coordinates: `(φ x)_w = Σ_a x_a (φ e_a)_w`. -/
theorem map_apply_eq_sum (x : V → 𝕜) (w : V) :
    φ x w = ∑ a, x a * φ (Pi.single a 1) w := by
  have hx : x = ∑ a, x a • (Pi.single a 1 : V → 𝕜) := by
    funext u
    simp [Finset.sum_apply, Pi.single_apply]
  calc φ x w = φ (∑ a, x a • (Pi.single a 1 : V → 𝕜)) w := by rw [← hx]
    _ = ∑ a, x a * φ (Pi.single a 1) w := by simp [map_sum, map_smul, Finset.sum_apply]

theorem map_ind_apply (S : Finset V) (w : V) :
    φ (ind S) w = ∑ a ∈ S, φ (Pi.single a 1) w := by
  rw [ind_eq_sum, map_sum, Finset.sum_apply]

/-- `(QP)_{aa} = 1`. -/
theorem sum_P_Q (a : V) : ∑ w, φ (Pi.single a 1) w * φ.symm (Pi.single w 1) a = 1 := by
  have h := map_apply_eq_sum φ.symm (φ (Pi.single a 1)) a
  rw [LinearEquiv.symm_apply_apply, Pi.single_eq_same] at h
  exact h.symm

/-- `(PQ)_{ww} = 1`. -/
theorem sum_Q_P (w : V) : ∑ a, φ.symm (Pi.single w 1) a * φ (Pi.single a 1) w = 1 := by
  have h := map_apply_eq_sum φ (φ.symm (Pi.single w 1)) w
  rw [LinearEquiv.apply_symm_apply, Pi.single_eq_same] at h
  exact h.symm

/-- Every column of `P` is non-zero (`φ` injective). -/
theorem col_ne_zero (a : V) : ∃ i, φ (Pi.single a 1) i ≠ 0 := by
  by_contra h
  have h' : φ (Pi.single a 1) = 0 := funext fun i => by
    by_contra hi
    exact h ⟨i, hi⟩
  have h0 : (Pi.single a 1 : V → 𝕜) = 0 := φ.injective (h'.trans (map_zero φ).symm)
  have := congrFun h0 a
  simp at this

/-- Every row of `P` is non-zero (`φ` surjective). -/
theorem row_ne_zero (i : V) : ∃ a, φ (Pi.single a 1) i ≠ 0 := by
  by_contra h
  have h' : ∀ a, φ (Pi.single a 1) i = 0 := fun a => by
    by_contra ha
    exact h ⟨a, ha⟩
  have := sum_Q_P φ i
  simp [h'] at this

/-- Trace trick (replaces the counting in Theorem N): if `(φ e_a)_w ≠ 0` always gives
`(a ∈ A ↔ w ∈ B)`, then `|A| = |B|` in `𝕜` (`|A| = Σ_{a∈A} (QP)_{aa} = Σ_{w∈B} (PQ)_{ww} = |B|`). -/
theorem card_eq_of_support (A B : Finset V)
    (h : ∀ a w, φ (Pi.single a 1) w ≠ 0 → (a ∈ A ↔ w ∈ B)) : (#A : 𝕜) = #B := by
  have key : ∀ a w, (if a ∈ A then (1 : 𝕜) else 0) * (φ (Pi.single a 1) w * φ.symm (Pi.single w 1) a)
      = (if w ∈ B then (1 : 𝕜) else 0) * (φ (Pi.single a 1) w * φ.symm (Pi.single w 1) a) := by
    intro a w
    by_cases hp : φ (Pi.single a 1) w = 0
    · simp [hp]
    · by_cases ha : a ∈ A
      · have hw : w ∈ B := (h a w hp).mp ha
        simp [ha, hw]
      · have hw : w ∉ B := fun hw => ha ((h a w hp).mpr hw)
        simp [ha, hw]
  calc (#A : 𝕜)
      = ∑ a, (if a ∈ A then (1 : 𝕜) else 0) * ∑ w, φ (Pi.single a 1) w * φ.symm (Pi.single w 1) a := by
        simp only [sum_P_Q, mul_one, sum_indicator]
    _ = ∑ a, ∑ w, (if w ∈ B then (1 : 𝕜) else 0) * (φ (Pi.single a 1) w * φ.symm (Pi.single w 1) a) := by
        simp only [Finset.mul_sum, key]
    _ = ∑ w, ∑ a, (if w ∈ B then (1 : 𝕜) else 0) * (φ.symm (Pi.single w 1) a * φ (Pi.single a 1) w) := by
        rw [Finset.sum_comm]
        exact Finset.sum_congr rfl (fun w _ => Finset.sum_congr rfl (fun a _ => by ring))
    _ = ∑ w, (if w ∈ B then (1 : 𝕜) else 0) * ∑ a, φ.symm (Pi.single w 1) a * φ (Pi.single a 1) w := by
        simp only [Finset.mul_sum]
    _ = #B := by simp only [sum_Q_P, mul_one, sum_indicator]

variable {M M' : V → V → 𝕜}

/-- `φ(e_a) φ(e_b) = φ(e_a e_b)` for an algebra isomorphism. -/
theorem sq_expand (hφ : ∀ x y, φ (evoMul M x y) = evoMul M' (φ x) (φ y)) (a b : V) :
    evoMul M' (φ (Pi.single a 1)) (φ (Pi.single b 1)) = if a = b then φ (M a) else 0 := by
  rw [← hφ, evoMul_single_single]
  split_ifs <;> simp

/-- **Lemma 1** (NOTE §2.2 = BCS22 Thm 2.4(i) "⇒", in coordinates): for the natural basis
`f_a = φ(e_a)` of the target, if the `e'_i`-coordinate of `f_a` is non-zero, then
`e'_i² ∈ 𝕜 · φ(e_a²)`. (Proof: `f_a e'_i = (f_a)_i e'_i²` and
`f_a e'_i = φ(e_a · φ⁻¹ e'_i) = (φ⁻¹ e'_i)_a φ(e_a²)`.) -/
theorem lemma1 (hφ : ∀ x y, φ (evoMul M x y) = evoMul M' (φ x) (φ y)) (a i : V)
    (h : φ (Pi.single a 1) i ≠ 0) :
    M' i = ((φ (Pi.single a 1) i)⁻¹ * φ.symm (Pi.single i 1) a) • φ (M a) := by
  have h1 : evoMul M' (φ (Pi.single a 1)) (Pi.single i 1) = φ (Pi.single a 1) i • M' i :=
    evoMul_single_right M' _ i
  have h2 : evoMul M' (φ (Pi.single a 1)) (Pi.single i 1)
      = φ.symm (Pi.single i 1) a • φ (M a) := by
    conv_lhs => rw [← φ.apply_symm_apply (Pi.single i 1)]
    rw [← hφ, evoMul_single_left, map_smul]
  rw [mul_smul, ← h2, h1, smul_smul, inv_mul_cancel₀ h, one_smul]

end linear

end EvoAlg

