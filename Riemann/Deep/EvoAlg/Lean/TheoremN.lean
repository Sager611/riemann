import Riemann.Deep.EvoAlg.Lean.Basic

/-!
# Theorem N (NOTE §2.2) for an arbitrary non-degenerate evolution algebra, in coordinates

The evolution algebra is `(𝕜^V, evoMul M)` with natural basis `(e_i)`; non-degenerate means
`e_i² = M i ≠ 0` for all `i`. A natural basis `B'` (indexed by `V`, as `|B'| = |V|`) is
`f_a = ψ e_a` for a linear equivalence `ψ` with `f_a f_b = 0` (`a ≠ b`). The class of `i` is
`{j : e_j² ∈ 𝕜 e_i²}`; on a class `T` with `e_i² = λ_i u_T`, `b_T(x, y) = Σ_{i ∈ T} λ_i x_i y_i`.
-/

open Finset

namespace EvoAlg

variable {V : Type*} [Fintype V] [DecidableEq V] {𝕜 : Type*} [Field 𝕜]

omit [DecidableEq V] in
theorem evoMul_sum_smul_right (M : V → V → 𝕜) (x : V → 𝕜) (β : V → 𝕜) (f : V → V → 𝕜) :
    evoMul M x (∑ b, β b • f b) = ∑ b, β b • evoMul M x (f b) := by
  funext w
  simp only [evoMul, Finset.sum_apply, Pi.smul_apply, smul_eq_mul, Finset.mul_sum, Finset.sum_mul]
  conv_lhs => rw [Finset.sum_comm]
  exact Finset.sum_congr rfl (fun b _ => Finset.sum_congr rfl (fun v _ => by ring))

/-- **Lemma 1** (BCS22 Thm 2.4(i) "⇒") for an arbitrary natural basis `f_a = ψ e_a`:
`(f_a)_i ≠ 0 ⟹ e_i² ∈ 𝕜 f_a²`. -/
theorem lemma1_nat (M : V → V → 𝕜) (ψ : (V → 𝕜) ≃ₗ[𝕜] (V → 𝕜))
    (horth : ∀ a b, a ≠ b → evoMul M (ψ (Pi.single a 1)) (ψ (Pi.single b 1)) = 0)
    (a i : V) (h : ψ (Pi.single a 1) i ≠ 0) :
    M i = ((ψ (Pi.single a 1) i)⁻¹ * ψ.symm (Pi.single i 1) a) •
      evoMul M (ψ (Pi.single a 1)) (ψ (Pi.single a 1)) := by
  have h1 : evoMul M (ψ (Pi.single a 1)) (Pi.single i 1) = ψ (Pi.single a 1) i • M i :=
    evoMul_single_right M _ i
  have hexp : (Pi.single i 1 : V → 𝕜) = ∑ b, ψ.symm (Pi.single i 1) b • ψ (Pi.single b 1) := by
    funext w
    rw [Finset.sum_apply]
    simp only [Pi.smul_apply, smul_eq_mul]
    rw [← map_apply_eq_sum ψ (ψ.symm (Pi.single i 1)) w, LinearEquiv.apply_symm_apply]
  have h2 : evoMul M (ψ (Pi.single a 1)) (Pi.single i 1)
      = ψ.symm (Pi.single i 1) a • evoMul M (ψ (Pi.single a 1)) (ψ (Pi.single a 1)) := by
    conv_lhs => rw [hexp]
    rw [evoMul_sum_smul_right, Finset.sum_eq_single a
      (fun b _ hb => by rw [horth a b (Ne.symm hb), smul_zero])
      (fun h => absurd (Finset.mem_univ a) h)]
  rw [mul_smul, ← h2, h1, smul_smul, inv_mul_cancel₀ h, one_smul]

open scoped Classical in
/-- **Theorem N, "⇒"** (NOTE §2.2). In a non-degenerate evolution algebra every natural basis
`(f_a)` is a union of `b_T`-orthogonal bases of the blocks `E_T`:
1. the support of each `f_a` lies in one class of proportional squares (so `f_a ∈ E_T`);
2. `f_a² ≠ 0`;
3. if `e_i² = λ_i u` on `T` (`u ≠ 0`) and `f_a ∈ E_T`, then `b_T(f_a, f_b) = 0` for `b ≠ a`;
4. for `T` a union of classes, the number of `f_a` in `E_T` is `|T|` (in `𝕜`; in characteristic 0
   this is `|B'_T| = |T|`). -/
theorem theoremN_mp (M : V → V → 𝕜) (hnd : ∀ i, M i ≠ 0) (ψ : (V → 𝕜) ≃ₗ[𝕜] (V → 𝕜))
    (horth : ∀ a b, a ≠ b → evoMul M (ψ (Pi.single a 1)) (ψ (Pi.single b 1)) = 0) :
    (∀ a i j, ψ (Pi.single a 1) i ≠ 0 → ψ (Pi.single a 1) j ≠ 0 →
      ∃ μ : 𝕜, μ ≠ 0 ∧ M i = μ • M j) ∧
    (∀ a, evoMul M (ψ (Pi.single a 1)) (ψ (Pi.single a 1)) ≠ 0) ∧
    (∀ (T : Finset V) (lam u : V → 𝕜), u ≠ 0 → (∀ i ∈ T, M i = lam i • u) →
      ∀ a b, a ≠ b → (∀ i, ψ (Pi.single a 1) i ≠ 0 → i ∈ T) →
        ∑ i ∈ T, lam i * ψ (Pi.single a 1) i * ψ (Pi.single b 1) i = 0) ∧
    (∀ T : Finset V, (∀ i j, i ∈ T → (∃ μ : 𝕜, M j = μ • M i) → j ∈ T) →
      ((univ.filter (fun a => ∃ i ∈ T, ψ (Pi.single a 1) i ≠ 0)).card : 𝕜) = #T) := by
  -- the coefficient in Lemma 1 is non-zero
  have hκ : ∀ a i, ψ (Pi.single a 1) i ≠ 0 →
      ∃ κ : 𝕜, κ ≠ 0 ∧ M i = κ • evoMul M (ψ (Pi.single a 1)) (ψ (Pi.single a 1)) := by
    intro a i h
    refine ⟨_, fun h0 => hnd i ?_, lemma1_nat M ψ horth a i h⟩
    rw [lemma1_nat M ψ horth a i h, h0, zero_smul]
  have hN1 : ∀ a i j, ψ (Pi.single a 1) i ≠ 0 → ψ (Pi.single a 1) j ≠ 0 →
      ∃ μ : 𝕜, μ ≠ 0 ∧ M i = μ • M j := by
    intro a i j hi hj
    obtain ⟨κ, hκ0, hκi⟩ := hκ a i hi
    obtain ⟨κ', hκ0', hκj⟩ := hκ a j hj
    refine ⟨κ * κ'⁻¹, mul_ne_zero hκ0 (inv_ne_zero hκ0'), ?_⟩
    rw [hκi, hκj, smul_smul, mul_assoc, inv_mul_cancel₀ hκ0', mul_one]
  refine ⟨hN1, fun a h0 => ?_, ?_, ?_⟩
  · obtain ⟨i, hi⟩ := col_ne_zero ψ a
    obtain ⟨κ, _, hκi⟩ := hκ a i hi
    exact hnd i (by rw [hκi, h0, smul_zero])
  · intro T lam u hu hM a b hab hsupp
    obtain ⟨w, hw⟩ : ∃ w, u w ≠ 0 := by
      by_contra h
      exact hu (funext fun x => by
        by_contra hx
        exact h ⟨x, hx⟩)
    have h0 := congrFun (horth a b hab) w
    rw [evoMul_apply, Pi.zero_apply] at h0
    rw [← Finset.sum_subset (Finset.subset_univ T) (fun i _ hi => by
      have : ψ (Pi.single a 1) i = 0 := by
        by_contra h
        exact hi (hsupp i h)
      rw [this, zero_mul, zero_mul])] at h0
    have h1 : ∑ i ∈ T, ψ (Pi.single a 1) i * ψ (Pi.single b 1) i * M i w
        = (∑ i ∈ T, lam i * ψ (Pi.single a 1) i * ψ (Pi.single b 1) i) * u w := by
      rw [Finset.sum_mul]
      refine Finset.sum_congr rfl (fun i hi => ?_)
      rw [hM i hi, Pi.smul_apply, smul_eq_mul]
      ring
    rw [h1] at h0
    exact (mul_eq_zero.mp h0).resolve_right hw
  · intro T hT
    refine card_eq_of_support ψ _ _ (fun a w hw => ?_)
    rw [Finset.mem_filter]
    constructor
    · rintro ⟨_, i, hiT, hi⟩
      obtain ⟨μ, _, hμ⟩ := hN1 a w i hw hi
      exact hT i w hiT ⟨μ, hμ⟩
    · intro hwT
      exact ⟨Finset.mem_univ a, w, hwT, hw⟩

/-- **Theorem N, "⇐"** (NOTE §2.2): if `e_i² = λ_i u_i` with `u` constant on the support of each
`f_a` (so `f_a ∈ E_T`) and `b_T(f_a, f_b) = 0` for `a ≠ b`, then `(f_a)` is a natural basis. -/
theorem theoremN_mpr (M : V → V → 𝕜) (ψ : (V → 𝕜) ≃ₗ[𝕜] (V → 𝕜)) (lam : V → 𝕜)
    (u : V → V → 𝕜) (hM : ∀ i, M i = lam i • u i)
    (hblock : ∀ a i j, ψ (Pi.single a 1) i ≠ 0 → ψ (Pi.single a 1) j ≠ 0 → u i = u j)
    (horth : ∀ a b, a ≠ b → ∀ i, ψ (Pi.single a 1) i ≠ 0 →
      ∑ j, lam j * ψ (Pi.single a 1) j * ψ (Pi.single b 1) j = 0) :
    ∀ a b, a ≠ b → evoMul M (ψ (Pi.single a 1)) (ψ (Pi.single b 1)) = 0 := by
  intro a b hab
  obtain ⟨i, hi⟩ := col_ne_zero ψ a
  funext w
  rw [evoMul_apply, Pi.zero_apply]
  have h1 : ∀ j, ψ (Pi.single a 1) j * ψ (Pi.single b 1) j * M j w
      = (lam j * ψ (Pi.single a 1) j * ψ (Pi.single b 1) j) * u i w := by
    intro j
    by_cases hj : ψ (Pi.single a 1) j = 0
    · simp [hj]
    · rw [hM j, Pi.smul_apply, smul_eq_mul, hblock a j i hj hi]
      ring
  rw [Finset.sum_congr rfl (fun j _ => h1 j), ← Finset.sum_mul, horth a b hab i hi, zero_mul]

end EvoAlg

