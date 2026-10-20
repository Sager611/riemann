import Riemann.Fayed.Defs
import Mathlib.Analysis.Calculus.FDeriv.Partial

/-!
# Fayed.Inputs.leibniz44 (Leibniz rule for F) from Mathlib

Substitute u = eʳ: ψ(eʳ, y) = A(y) e^{rW(y)} ζ(W(y)) is jointly continuous with jointly continuous
r-derivative on ℝ × ℝ. P(r, w) = ∫₀^w … has continuous partials (dominated differentiation,
FTC, continuity of parametric primitives), so it is strictly differentiable
(hasStrictFDerivAt_uncurry_coprod); compose with r ↦ (r, ±(eʳ)⁸) and with log.
-/

open Complex Set Filter Topology MeasureTheory intervalIntegral

namespace Fayed.InputsProof

noncomputable section

variable (η σ₀ t₀ : ℝ)

/-- W(y) = 1 - η + i(t₀ + y). -/
def Wl (y : ℝ) : ℂ := ((1 - η : ℝ) : ℂ) + (t₀ + y) * I

/-- A(y) = 1/(1-σ₀-η+iy) - 1/(σ₀-η+iy). -/
def Al (y : ℝ) : ℂ := 1 / (((1 - σ₀ - η : ℝ) : ℂ) + y * I) - 1 / (((σ₀ - η : ℝ) : ℂ) + y * I)

/-- integrand in the variable r = log u. -/
def fl (r y : ℝ) : ℂ := Al η σ₀ y * cexp (r * Wl η t₀ y) * riemannZeta (Wl η t₀ y)

/-- its r-derivative. -/
def gl (r y : ℝ) : ℂ := Al η σ₀ y * (Wl η t₀ y * cexp (r * Wl η t₀ y)) * riemannZeta (Wl η t₀ y)

/-- P(r, w) = ∫₀^w fl r y dy. -/
def Pl (r w : ℝ) : ℂ := ∫ y in (0 : ℝ)..w, fl η σ₀ t₀ r y

variable {η σ₀ t₀}

lemma continuous_Wl : Continuous (Wl η t₀) := by unfold Wl; fun_prop

lemma ne_zero_of_re_pos {a : ℝ} (ha : 0 < a) (y : ℝ) : ((a : ℂ) + y * I) ≠ 0 := by
  intro h
  have := congrArg Complex.re h
  simp at this
  linarith

lemma continuous_Al (h1 : η < σ₀) (h2 : η < 1 - σ₀) : Continuous (Al η σ₀) := by
  unfold Al
  refine Continuous.sub ?_ ?_
  · exact continuous_const.div (by fun_prop) (ne_zero_of_re_pos (by linarith))
  · exact continuous_const.div (by fun_prop) (ne_zero_of_re_pos (by linarith))

lemma continuous_Zl (hη : 0 < η) : Continuous fun y => riemannZeta (Wl η t₀ y) := by
  rw [continuous_iff_continuousAt]
  intro y
  have hne : Wl η t₀ y ≠ 1 := by
    intro h
    have := congrArg Complex.re h
    simp [Wl] at this
    linarith
  exact (differentiableAt_riemannZeta hne).continuousAt.comp continuous_Wl.continuousAt

lemma continuous_fl (hη : 0 < η) (h1 : η < σ₀) (h2 : η < 1 - σ₀) :
    Continuous (Function.uncurry (fl η σ₀ t₀)) := by
  have hA := continuous_Al h1 h2
  have hW := continuous_Wl (η := η) (t₀ := t₀)
  have hZ := continuous_Zl (t₀ := t₀) hη
  unfold fl
  exact ((hA.comp continuous_snd).mul (Complex.continuous_exp.comp
    ((continuous_ofReal.comp continuous_fst).mul (hW.comp continuous_snd)))).mul
    (hZ.comp continuous_snd)

lemma continuous_gl (hη : 0 < η) (h1 : η < σ₀) (h2 : η < 1 - σ₀) :
    Continuous (Function.uncurry (gl η σ₀ t₀)) := by
  have hA := continuous_Al h1 h2
  have hW := continuous_Wl (η := η) (t₀ := t₀)
  have hZ := continuous_Zl (t₀ := t₀) hη
  unfold gl
  exact ((hA.comp continuous_snd).mul ((hW.comp continuous_snd).mul (Complex.continuous_exp.comp
    ((continuous_ofReal.comp continuous_fst).mul (hW.comp continuous_snd))))).mul
    (hZ.comp continuous_snd)

lemma cont_slice {f : ℝ → ℝ → ℂ} (hf : Continuous (Function.uncurry f)) (x : ℝ) :
    Continuous (f x) :=
  hf.comp (continuous_const.prodMk continuous_id)

lemma fl_hasDerivAt (r y : ℝ) : HasDerivAt (fun x => fl η σ₀ t₀ x y) (gl η σ₀ t₀ r y) r := by
  have h := (((hasDerivAt_id' r).ofReal_comp.mul_const (Wl η t₀ y)).cexp.const_mul
    (Al η σ₀ y)).mul_const (riemannZeta (Wl η t₀ y))
  refine h.congr_deriv ?_
  simp only [gl, ofReal_one, one_mul]
  ring

lemma Pl_hasDerivAt_fst (hη : 0 < η) (h1 : η < σ₀) (h2 : η < 1 - σ₀) (r w : ℝ) :
    HasDerivAt (fun x => Pl η σ₀ t₀ x w) (∫ y in (0 : ℝ)..w, gl η σ₀ t₀ r y) r := by
  have hf := continuous_fl (t₀ := t₀) hη h1 h2
  have hg := continuous_gl (t₀ := t₀) hη h1 h2
  have hK : IsCompact (Metric.closedBall r 1 ×ˢ uIcc (0 : ℝ) w) :=
    (isCompact_closedBall r 1).prod isCompact_uIcc
  obtain ⟨M, hM⟩ := hK.exists_bound_of_continuousOn hg.continuousOn
  have hcy : ∀ x : ℝ, Continuous (fl η σ₀ t₀ x) := cont_slice hf
  have hcg : Continuous (gl η σ₀ t₀ r) := cont_slice hg r
  have key := hasDerivAt_integral_of_dominated_loc_of_deriv_le (μ := volume) (a := 0) (b := w)
    (F := fun x y => fl η σ₀ t₀ x y) (F' := fun x y => gl η σ₀ t₀ x y) (x₀ := r)
    (bound := fun _ => M) (s := Metric.ball r 1) (Metric.ball_mem_nhds r one_pos)
    (Eventually.of_forall fun x => (hcy x).aestronglyMeasurable)
    ((hcy r).intervalIntegrable _ _) hcg.aestronglyMeasurable
    (Eventually.of_forall fun t ht x hx =>
      hM (x, t) ⟨Metric.ball_subset_closedBall hx, uIoc_subset_uIcc ht⟩)
    intervalIntegrable_const
    (Eventually.of_forall fun t _ x _ => fl_hasDerivAt x t)
  unfold Pl
  exact key.2

lemma Pl_hasDerivAt_snd (hη : 0 < η) (h1 : η < σ₀) (h2 : η < 1 - σ₀) (r w : ℝ) :
    HasDerivAt (fun x => Pl η σ₀ t₀ r x) (fl η σ₀ t₀ r w) w := by
  have hcy : Continuous (fl η σ₀ t₀ r) := cont_slice (continuous_fl (t₀ := t₀) hη h1 h2) r
  unfold Pl
  exact integral_hasDerivAt_right (f := fl η σ₀ t₀ r) (a := 0) (b := w)
    (hcy.intervalIntegrable _ _) (hcy.stronglyMeasurableAtFilter _ _) hcy.continuousAt

lemma smulRight_continuous {X : Type*} [TopologicalSpace X] {F : X → ℂ} (hF : Continuous F) :
    Continuous fun x => (1 : ℝ →L[ℝ] ℝ).smulRight (F x) :=
  (ContinuousLinearMap.smulRightL ℝ ℝ ℂ (1 : ℝ →L[ℝ] ℝ)).continuous.comp hF

lemma Pl_hasStrictFDerivAt (hη : 0 < η) (h1 : η < σ₀) (h2 : η < 1 - σ₀) (p : ℝ × ℝ) :
    HasStrictFDerivAt (Function.uncurry (Pl η σ₀ t₀))
      (((1 : ℝ →L[ℝ] ℝ).smulRight (∫ y in (0 : ℝ)..p.2, gl η σ₀ t₀ p.1 y)).coprod
        ((1 : ℝ →L[ℝ] ℝ).smulRight (fl η σ₀ t₀ p.1 p.2))) p := by
  have hf := continuous_fl (t₀ := t₀) hη h1 h2
  have hg := continuous_gl (t₀ := t₀) hη h1 h2
  have c1 : Continuous fun q : ℝ × ℝ =>
      (1 : ℝ →L[ℝ] ℝ).smulRight (∫ y in (0 : ℝ)..q.2, gl η σ₀ t₀ q.1 y) :=
    smulRight_continuous (intervalIntegral.continuous_parametric_primitive_of_continuous
      (a₀ := 0) hg)
  have c2 : Continuous fun q : ℝ × ℝ => (1 : ℝ →L[ℝ] ℝ).smulRight (fl η σ₀ t₀ q.1 q.2) :=
    smulRight_continuous hf
  exact hasStrictFDerivAt_uncurry_coprod
    (f₁ := fun r w => (1 : ℝ →L[ℝ] ℝ).smulRight (∫ y in (0 : ℝ)..w, gl η σ₀ t₀ r y))
    (f₂ := fun r w => (1 : ℝ →L[ℝ] ℝ).smulRight (fl η σ₀ t₀ r w))
    (Eventually.of_forall fun v => (Pl_hasDerivAt_fst hη h1 h2 v.1 v.2).hasFDerivAt)
    (Eventually.of_forall fun v => (Pl_hasDerivAt_snd hη h1 h2 v.1 v.2).hasFDerivAt)
    c1.continuousAt c2.continuousAt

/-- r ↦ P(r, c·(eʳ)⁸) has derivative ∫₀^{c e^{8r}} g + 8c e^{8r} f(r, c e^{8r}). -/
lemma Pl_curve (hη : 0 < η) (h1 : η < σ₀) (h2 : η < 1 - σ₀) (c r : ℝ) :
    HasDerivAt (fun x => Pl η σ₀ t₀ x (c * Real.exp x ^ 8))
      ((∫ y in (0 : ℝ)..(c * Real.exp r ^ 8), gl η σ₀ t₀ r y) +
        ((8 * c * Real.exp r ^ 8 : ℝ) : ℂ) * fl η σ₀ t₀ r (c * Real.exp r ^ 8)) r := by
  have hc : HasDerivAt (fun x => (x, c * Real.exp x ^ 8)) ((1 : ℝ), 8 * c * Real.exp r ^ 8) r := by
    refine (hasDerivAt_id' r).prodMk ?_
    have := ((Real.hasDerivAt_exp r).pow 8).const_mul c
    refine this.congr_deriv ?_
    push_cast; ring
  have := (Pl_hasStrictFDerivAt (t₀ := t₀) hη h1 h2 (r, c * Real.exp r ^ 8)).hasFDerivAt.comp_hasDerivAt
    r hc
  refine this.congr_deriv ?_
  rw [ContinuousLinearMap.coprod_apply, ContinuousLinearMap.smulRight_apply,
    ContinuousLinearMap.smulRight_apply, Complex.real_smul, Complex.real_smul]
  simp

lemma fl_log {v : ℝ} (hv : 0 < v) (y : ℝ) : fl η σ₀ t₀ (Real.log v) y = ψ η σ₀ t₀ v y := by
  have hv0 : (v : ℂ) ≠ 0 := by exact_mod_cast hv.ne'
  simp only [fl, ψ, Al, Wl]
  rw [cpow_def_of_ne_zero hv0, ← ofReal_log hv.le]

lemma Ffun_eq (hη : 0 < η) (h1 : η < σ₀) (h2 : η < 1 - σ₀) {v : ℝ} (hv : 0 < v) :
    Ffun η σ₀ t₀ v = ((1 / (2 * Real.pi) : ℝ) : ℂ) *
      (Pl η σ₀ t₀ (Real.log v) (1 * Real.exp (Real.log v) ^ 8) -
        Pl η σ₀ t₀ (Real.log v) (-1 * Real.exp (Real.log v) ^ 8)) := by
  have hcy : Continuous (fl η σ₀ t₀ (Real.log v)) :=
    cont_slice (continuous_fl (t₀ := t₀) hη h1 h2) _
  rw [Real.exp_log hv, one_mul, neg_one_mul, Ffun]
  unfold Pl
  rw [integral_interval_sub_left (hcy.intervalIntegrable _ _) (hcy.intervalIntegrable _ _)]
  congr 1
  exact intervalIntegral.integral_congr fun y _ => (fl_log hv y).symm

lemma gl_log {u : ℝ} (hu : 0 < u) (y : ℝ) :
    gl η σ₀ t₀ (Real.log u) y * ((u⁻¹ : ℝ) : ℂ) =
      (1 / (((1 - σ₀ - η : ℝ) : ℂ) + y * I) - 1 / (((σ₀ - η : ℝ) : ℂ) + y * I)) *
      ((((1 - η : ℝ) : ℂ) + (t₀ + y) * I) * (u : ℂ) ^ (((-η : ℝ) : ℂ) + (t₀ + y) * I)) *
      riemannZeta (((1 - η : ℝ) : ℂ) + (t₀ + y) * I) := by
  have hu0 : (u : ℂ) ≠ 0 := by exact_mod_cast hu.ne'
  have hpow : (u : ℂ) ^ (((-η : ℝ) : ℂ) + (t₀ + y) * I) =
      (u : ℂ) ^ (((1 - η : ℝ) : ℂ) + (t₀ + y) * I) * (u : ℂ)⁻¹ := by
    rw [show ((-η : ℝ) : ℂ) + (t₀ + y) * I = (((1 - η : ℝ) : ℂ) + (t₀ + y) * I) - 1 by
      push_cast; ring, cpow_sub _ _ hu0, cpow_one, div_eq_mul_inv]
  have hcp : cexp ((Real.log u : ℂ) * Wl η t₀ y) =
      (u : ℂ) ^ (((1 - η : ℝ) : ℂ) + (t₀ + y) * I) := by
    rw [cpow_def_of_ne_zero hu0, ← ofReal_log hu.le]
    rfl
  simp only [gl, Al]
  rw [hcp, hpow]
  simp only [Wl]
  push_cast
  ring

lemma alg_final (c G1 G2 f1 f2 : ℂ) {u : ℝ} (hu : u ≠ 0) :
    c * ((u⁻¹ : ℝ) • ((G1 + ((8 * 1 * u ^ 8 : ℝ) : ℂ) * f1) -
      (G2 + ((8 * (-1) * u ^ 8 : ℝ) : ℂ) * f2))) =
    c * ((G1 - G2) * ((u⁻¹ : ℝ) : ℂ) + 8 * (u : ℂ) ^ 7 * (f1 + f2)) := by
  have hu' : (u : ℂ) ≠ 0 := by exact_mod_cast hu
  rw [Complex.real_smul]
  push_cast
  field_simp
  ring

theorem leibniz44_holds : ∀ η σ₀ t₀ : ℝ, 0 < η → η < σ₀ → η < 1 - σ₀ → ∀ u : ℝ, 0 < u →
    HasDerivAt (Ffun η σ₀ t₀) (Fder η σ₀ t₀ u) u := by
  intro η σ₀ t₀ hη h1 h2 u hu
  have hg := continuous_gl (t₀ := t₀) hη h1 h2
  have hcg : Continuous (gl η σ₀ t₀ (Real.log u)) := cont_slice hg _
  have hD := (((Pl_curve (t₀ := t₀) hη h1 h2 1 (Real.log u)).sub
    (Pl_curve (t₀ := t₀) hη h1 h2 (-1) (Real.log u))).scomp u
    (Real.hasDerivAt_log hu.ne')).const_mul (((1 / (2 * Real.pi) : ℝ) : ℂ))
  have hev : Ffun η σ₀ t₀ =ᶠ[𝓝 u] fun v => ((1 / (2 * Real.pi) : ℝ) : ℂ) *
      (Pl η σ₀ t₀ (Real.log v) (1 * Real.exp (Real.log v) ^ 8) -
        Pl η σ₀ t₀ (Real.log v) (-1 * Real.exp (Real.log v) ^ 8)) := by
    filter_upwards [Ioi_mem_nhds hu] with v hv
    exact Ffun_eq hη h1 h2 hv
  refine (hD.congr_of_eventuallyEq hev).congr_deriv ?_
  have hI : (∫ y in (-(u ^ 8))..(u ^ 8), gl η σ₀ t₀ (Real.log u) y) * ((u⁻¹ : ℝ) : ℂ) =
      ∫ y in (-(u ^ 8))..(u ^ 8),
        (1 / (((1 - σ₀ - η : ℝ) : ℂ) + y * I) - 1 / (((σ₀ - η : ℝ) : ℂ) + y * I)) *
        ((((1 - η : ℝ) : ℂ) + (t₀ + y) * I) * (u : ℂ) ^ (((-η : ℝ) : ℂ) + (t₀ + y) * I)) *
        riemannZeta (((1 - η : ℝ) : ℂ) + (t₀ + y) * I) := by
    rw [← intervalIntegral.integral_mul_const]
    exact intervalIntegral.integral_congr fun y _ => gl_log hu y
  rw [Fder, ← hI, ← integral_interval_sub_left (hcg.intervalIntegrable 0 (u ^ 8))
      (hcg.intervalIntegrable 0 (-(u ^ 8))), ← fl_log hu (u ^ 8), ← fl_log hu (-(u ^ 8)),
    Real.exp_log hu]
  have e1 : (1 : ℝ) * u ^ 8 = u ^ 8 := one_mul _
  have e2 : (-1 : ℝ) * u ^ 8 = -(u ^ 8) := neg_one_mul _
  rw [e1, e2]
  exact alg_final _ _ _ _ _ hu.ne'

end

end Fayed.InputsProof
