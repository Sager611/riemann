import Riemann.Deep.LeanNBRate.Dual1

/-!
# Dual problem: the representer of evaluation at 0 in the span of a configuration
-/

open MeasureTheory Filter Topology Set Complex ComplexConjugate
open scoped InnerProductSpace

noncomputable section

namespace Riemann.Deep.LeanNBRate

namespace ExpSum

/-- Same nodes and multiplicities, new coefficients. -/
def withCoef (P : ExpSum) (c : Fin P.J → ℕ → ℂ) : ExpSum := ⟨P.J, P.node, P.mult, c⟩

theorem eval_withCoef_add (P : ExpSum) (c d : Fin P.J → ℕ → ℂ) (u : ℝ) :
    (P.withCoef (c + d)).eval u = (P.withCoef c).eval u + (P.withCoef d).eval u := by
  unfold eval withCoef
  simp only [Pi.add_apply, add_mul, Finset.sum_add_distrib]

theorem eval_withCoef_smul (P : ExpSum) (a : ℂ) (c : Fin P.J → ℕ → ℂ) (u : ℝ) :
    (P.withCoef (a • c)).eval u = a * (P.withCoef c).eval u := by
  unfold eval withCoef
  simp only [Pi.smul_apply, smul_eq_mul, Finset.mul_sum]
  refine Finset.sum_congr rfl fun j _ => Finset.sum_congr rfl fun k _ => by ring

theorem dmod_eval_withCoef_add (P : ExpSum) (c d : Fin P.J → ℕ → ℂ) (u : ℝ) :
    ((P.withCoef (c + d)).dmod 0).eval u = ((P.withCoef c).dmod 0).eval u + ((P.withCoef d).dmod 0).eval u := by
  have h1 := (P.withCoef (c + d)).hasDerivAt_eval u
  have h2 := ((P.withCoef c).hasDerivAt_eval u).add ((P.withCoef d).hasDerivAt_eval u)
  have e : (P.withCoef (c + d)).eval = fun v => (P.withCoef c).eval v + (P.withCoef d).eval v :=
    funext fun v => eval_withCoef_add P c d v
  rw [e] at h1
  exact h1.unique h2

theorem dmod_eval_withCoef_smul (P : ExpSum) (a : ℂ) (c : Fin P.J → ℕ → ℂ) (u : ℝ) :
    ((P.withCoef (a • c)).dmod 0).eval u = a * ((P.withCoef c).dmod 0).eval u := by
  have h1 := (P.withCoef (a • c)).hasDerivAt_eval u
  have h2 := ((P.withCoef c).hasDerivAt_eval u).const_mul a
  have e : (P.withCoef (a • c)).eval = fun v => a * (P.withCoef c).eval v :=
    funext fun v => eval_withCoef_smul P a c v
  rw [e] at h1
  exact h1.unique h2

/-- The finite index set of coefficients. -/
abbrev Idx (P : ExpSum) : Type := Σ j : Fin P.J, Fin (P.mult j)

/-- Coefficients from the finite index set. -/
def icoef (P : ExpSum) (c : P.Idx → ℂ) : Fin P.J → ℕ → ℂ :=
  fun j k => if h : k < P.mult j then c ⟨j, ⟨k, h⟩⟩ else 0

/-- Truncation to the finite index set. -/
def tcoef (P : ExpSum) (d : Fin P.J → ℕ → ℂ) : P.Idx → ℂ := fun i => d i.1 i.2

theorem icoef_add (P : ExpSum) (c d : P.Idx → ℂ) : P.icoef (c + d) = P.icoef c + P.icoef d := by
  funext j k; simp only [icoef, Pi.add_apply]
  by_cases h : k < P.mult j
  · simp [h]
  · simp [h]

theorem icoef_smul (P : ExpSum) (a : ℂ) (c : P.Idx → ℂ) : P.icoef (a • c) = a • P.icoef c := by
  funext j k; simp only [icoef, Pi.smul_apply]
  by_cases h : k < P.mult j
  · simp [h]
  · simp [h]

theorem eval_icoef_tcoef (P : ExpSum) (d : Fin P.J → ℕ → ℂ) (u : ℝ) :
    (P.withCoef (P.icoef (P.tcoef d))).eval u = (P.withCoef d).eval u := by
  unfold eval
  refine Finset.sum_congr rfl fun j _ => Finset.sum_congr rfl fun k hk => ?_
  have hk' : k < P.mult j := Finset.mem_range.mp hk
  simp [withCoef, icoef, tcoef, hk']

theorem dmod_eval_icoef_tcoef (P : ExpSum) (d : Fin P.J → ℕ → ℂ) (u : ℝ) :
    ((P.withCoef (P.icoef (P.tcoef d))).dmod 0).eval u = ((P.withCoef d).dmod 0).eval u := by
  have h1 := (P.withCoef (P.icoef (P.tcoef d))).hasDerivAt_eval u
  have h2 := (P.withCoef d).hasDerivAt_eval u
  have e : (P.withCoef (P.icoef (P.tcoef d))).eval = (P.withCoef d).eval :=
    funext fun v => eval_icoef_tcoef P d v
  rw [e] at h1
  exact h1.unique h2

end ExpSum

theorem Temb_congr {ℓ : ℝ} {F G F' G' : ℝ → ℂ} (hFG : F = G) (hFG' : F' = G') (hF : Continuous F)
    (hF' : Continuous F') (hG : Continuous G) (hG' : Continuous G') :
    Temb ℓ F F' hF hF' = Temb ℓ G G' hG hG' := by
  subst hFG; subst hFG'; rfl

theorem Temb_add {ℓ : ℝ} (F F' G G' : ℝ → ℂ) (hF : Continuous F) (hF' : Continuous F')
    (hG : Continuous G) (hG' : Continuous G') :
    Temb ℓ (fun u => F u + G u) (fun u => F' u + G' u) (hF.add hG) (hF'.add hG') =
      Temb ℓ F F' hF hF' + Temb ℓ G G' hG hG' := by
  unfold Temb
  rw [← WithLp.toLp_add, Prod.mk_add_mk, ← WithLp.toLp_add, Prod.mk_add_mk, ← WithLp.toLp_add,
    Prod.mk_add_mk, ← MemLp.toLp_add, ← MemLp.toLp_add]
  congr 1
  refine Prod.ext ?_ ?_
  · dsimp only
    congr 1
    refine Prod.ext ?_ ?_
    · rfl
    · exact MemLp.toLp_congr (memLp_muL_of_continuous (by fun_prop)) (memLp_muL_of_continuous (by fun_prop))
        (Filter.Eventually.of_forall fun u => by simp only [Pi.add_apply]; ring)
  · dsimp only
    congr 1
    refine Prod.ext ?_ ?_ <;> (dsimp only; ring)

theorem Temb_smul {ℓ : ℝ} (a : ℂ) (F F' : ℝ → ℂ) (hF : Continuous F) (hF' : Continuous F') :
    Temb ℓ (fun u => a * F u) (fun u => a * F' u) (continuous_const.mul hF) (continuous_const.mul hF') =
      a • Temb ℓ F F' hF hF' := by
  unfold Temb
  rw [← WithLp.toLp_smul, Prod.smul_mk, ← WithLp.toLp_smul, Prod.smul_mk, ← WithLp.toLp_smul,
    Prod.smul_mk, ← MemLp.toLp_const_smul, ← MemLp.toLp_const_smul]
  congr 1
  refine Prod.ext ?_ ?_
  · dsimp only
    congr 1
    refine Prod.ext ?_ ?_
    · rfl
    · exact MemLp.toLp_congr (memLp_muL_of_continuous (by fun_prop))
        (memLp_muL_of_continuous (continuous_const.smul (by fun_prop)))
        (Filter.Eventually.of_forall fun u => by
          show a * F u * (1 / 2) = a * (F u * (1 / 2)); ring)
  · dsimp only
    congr 1
    refine Prod.ext ?_ ?_ <;> (dsimp only; rw [smul_eq_mul]; ring)

/-- The embedding of a configuration with given finite coefficients, as a linear map. -/
def Lam (ℓ : ℝ) (P : ExpSum) : (P.Idx → ℂ) →ₗ[ℂ] WSp ℓ where
  toFun c := Temb ℓ (P.withCoef (P.icoef c)).eval ((P.withCoef (P.icoef c)).dmod 0).eval
    (P.withCoef (P.icoef c)).continuous_eval ((P.withCoef (P.icoef c)).dmod 0).continuous_eval
  map_add' c d := by
    rw [← Temb_add]
    refine Temb_congr ?_ ?_ _ _ _ _
    · funext u; rw [P.icoef_add, ExpSum.eval_withCoef_add]
    · funext u; rw [P.icoef_add, ExpSum.dmod_eval_withCoef_add]
  map_smul' a c := by
    simp only [RingHom.id_apply]
    rw [← Temb_smul]
    refine Temb_congr ?_ ?_ _ _ _ _
    · funext u; rw [P.icoef_smul, ExpSum.eval_withCoef_smul]
    · funext u; rw [P.icoef_smul, ExpSum.dmod_eval_withCoef_smul]

/-- **The Riesz representer** of f ↦ f(0) for the OU form on the span of a configuration. -/
theorem exists_representer (P : ExpSum) {ℓ : ℝ} (hℓ : 0 ≤ ℓ) :
    ∃ c : Fin P.J → ℕ → ℂ, ∀ d : Fin P.J → ℕ → ℂ,
      NForm ℓ (P.withCoef d).eval ((P.withCoef d).dmod 0).eval (P.withCoef c).eval
        ((P.withCoef c).dmod 0).eval = (P.withCoef d).eval 0 := by
  set U := LinearMap.range (Lam ℓ P)
  have hv : U.starProjection (evec ℓ) ∈ U := U.starProjection_apply_mem (evec ℓ)
  obtain ⟨c₀, hc₀⟩ := LinearMap.mem_range.mp hv
  refine ⟨P.icoef c₀, fun d => ?_⟩
  set d₀ := P.tcoef d
  have hd : (P.withCoef d).eval = (P.withCoef (P.icoef d₀)).eval := funext fun u => (P.eval_icoef_tcoef d u).symm
  have hd' : ((P.withCoef d).dmod 0).eval = ((P.withCoef (P.icoef d₀)).dmod 0).eval :=
    funext fun u => (P.dmod_eval_icoef_tcoef d u).symm
  rw [hd, hd']
  rw [← inner_Temb hℓ _ _ _ _ (P.withCoef (P.icoef d₀)).continuous_eval
    ((P.withCoef (P.icoef d₀)).dmod 0).continuous_eval (P.withCoef (P.icoef c₀)).continuous_eval
    ((P.withCoef (P.icoef c₀)).dmod 0).continuous_eval]
  have e1 : Temb ℓ (P.withCoef (P.icoef c₀)).eval ((P.withCoef (P.icoef c₀)).dmod 0).eval
      (P.withCoef (P.icoef c₀)).continuous_eval ((P.withCoef (P.icoef c₀)).dmod 0).continuous_eval =
      U.starProjection (evec ℓ) := hc₀
  have e2 : Temb ℓ (P.withCoef (P.icoef d₀)).eval ((P.withCoef (P.icoef d₀)).dmod 0).eval
      (P.withCoef (P.icoef d₀)).continuous_eval ((P.withCoef (P.icoef d₀)).dmod 0).continuous_eval =
      Lam ℓ P d₀ := rfl
  rw [e1, e2]
  have hmem : Lam ℓ P d₀ ∈ U := LinearMap.mem_range_self _ _
  have h0 := U.starProjection_inner_eq_zero (evec ℓ) (Lam ℓ P d₀) hmem
  rw [inner_sub_left, sub_eq_zero] at h0
  rw [← h0, ← e2]
  exact inner_evec _ _ _ _

end Riemann.Deep.LeanNBRate

