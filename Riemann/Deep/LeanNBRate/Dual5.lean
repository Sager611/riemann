import Riemann.Deep.LeanNBRate.Dual4

/-!
# Dual problem: E_low is entire and interpolates 1 at the configuration's zeros
-/

open MeasureTheory Filter Topology Set Complex ComplexConjugate
open scoped InnerProductSpace

noncomputable section

namespace Riemann.Deep.LeanNBRate

theorem hasDerivAt_eS (s : ℂ) (u : ℝ) : HasDerivAt (eS s) (eS' s u) u := by
  unfold eS eS'
  have h := ((hasDerivAt_id (u : ℂ)).const_mul (-(s - 1 / 2))).cexp
  have h2 := h.comp_ofReal (z := u)
  simpa [mul_comm] using h2

/-- The Laplace-type integral ∫_0^ℓ e^{-(s-1/2)u} q(u) du is entire. -/
theorem differentiable_laplace {ℓ : ℝ} (hℓ : 0 ≤ ℓ) {q : ℝ → ℂ} (hq : Continuous q) :
    Differentiable ℂ fun s : ℂ => ∫ u in (0 : ℝ)..ℓ, cexp (-(s - 1 / 2) * u) * q u := by
  intro s₀
  obtain ⟨M, hM⟩ := (isCompact_Icc (a := (0 : ℝ)) (b := ℓ)).exists_bound_of_continuousOn hq.continuousOn
  set B : ℝ := ℓ * Real.exp ((‖s₀‖ + 2) * ℓ) * max M 0
  have key := intervalIntegral.hasDerivAt_integral_of_dominated_loc_of_deriv_le (μ := volume) (a := 0) (b := ℓ)
    (F := fun s u => cexp (-(s - 1 / 2) * u) * q u)
    (F' := fun s u => (-(u : ℂ)) * cexp (-(s - 1 / 2) * u) * q u) (x₀ := s₀)
    (s := Metric.ball s₀ 1) (bound := fun _ => B) (Metric.ball_mem_nhds s₀ one_pos)
    (Filter.Eventually.of_forall fun s => (by fun_prop : Continuous fun u : ℝ =>
      cexp (-(s - 1 / 2) * u) * q u).aestronglyMeasurable)
    ((by fun_prop : Continuous fun u : ℝ => cexp (-(s₀ - 1 / 2) * u) * q u).intervalIntegrable _ _)
    ((by fun_prop : Continuous fun u : ℝ => (-(u : ℂ)) * cexp (-(s₀ - 1 / 2) * u) * q u).aestronglyMeasurable)
    ?_ intervalIntegrable_const ?_
  · exact key.2.differentiableAt
  · refine Filter.Eventually.of_forall fun u hu s hs => ?_
    rw [Set.uIoc_of_le hℓ] at hu
    have hu0 : 0 ≤ u := hu.1.le
    have hul : u ≤ ℓ := hu.2
    have hs' : ‖s‖ ≤ ‖s₀‖ + 1 := by
      have := Metric.mem_ball.mp hs
      have h2 := norm_le_norm_add_norm_sub' s s₀
      rw [dist_eq_norm] at this
      linarith
    rw [norm_mul, norm_mul, norm_neg, Complex.norm_real, Real.norm_eq_abs, abs_of_nonneg hu0, Complex.norm_exp]
    have hre : (-(s - 1 / 2) * (u : ℂ)).re ≤ (‖s₀‖ + 2) * ℓ := by
      have : (-(s - 1 / 2) * (u : ℂ)).re = -(s.re - 1 / 2) * u := by simp
      rw [this]
      have h1 : |s.re| ≤ ‖s‖ := Complex.abs_re_le_norm s
      have h2 : -(s.re - 1 / 2) ≤ ‖s₀‖ + 2 := by have := neg_abs_le s.re; linarith
      calc -(s.re - 1 / 2) * u ≤ (‖s₀‖ + 2) * u := mul_le_mul_of_nonneg_right h2 hu0
        _ ≤ (‖s₀‖ + 2) * ℓ := mul_le_mul_of_nonneg_left hul (by positivity)
    have hqu : ‖q u‖ ≤ max M 0 := (hM u ⟨hu0, hul⟩).trans (le_max_left _ _)
    have hfin : u * Real.exp (-(s - 1 / 2) * (u : ℂ)).re * ‖q u‖ ≤
        ℓ * Real.exp ((‖s₀‖ + 2) * ℓ) * max M 0 := by
      gcongr
    exact hfin
  · refine Filter.Eventually.of_forall fun u _ s _ => ?_
    have h1 : HasDerivAt (fun y : ℂ => -(y - 1 / 2) * (u : ℂ)) (-(u : ℂ)) s := by
      have := ((hasDerivAt_id s).sub_const (1 / 2 : ℂ)).const_mul (-(u : ℂ))
      convert this using 1
      · ext y; simp only [id]; ring
      · ring
    have h := (h1.cexp).mul_const (q u)
    convert h using 1
    ring

/-- E_low in Green form. -/
theorem Elow_green {ℓ : ℝ} (p p' p'' : ℝ → ℂ) (hp' : ∀ x, HasDerivAt p' (p'' x) x) (hcp : Continuous p)
    (hcp' : Continuous p') (hcp'' : Continuous p'') (s : ℂ) :
    Elow ℓ p p' s = (∫ u in (0 : ℝ)..ℓ, cexp (-(s - 1 / 2) * u) * conj (gDens p p' p'' u)) +
      conj (gA p p') + cexp (-(s - 1 / 2) * ℓ) * conj (gB ℓ p p') := by
  unfold Elow
  rw [NForm_green (eS s) (eS' s) p p' p'' (hasDerivAt_eS s) hp' (continuous_eS s) (continuous_eS' s) hcp hcp' hcp'']
  unfold eS
  simp

theorem differentiable_Elow {ℓ : ℝ} (hℓ : 0 ≤ ℓ) (p p' p'' : ℝ → ℂ) (hp' : ∀ x, HasDerivAt p' (p'' x) x)
    (hcp : Continuous p) (hcp' : Continuous p') (hcp'' : Continuous p'') :
    Differentiable ℂ (Elow ℓ p p') := by
  have e : Elow ℓ p p' = fun s => (∫ u in (0 : ℝ)..ℓ, cexp (-(s - 1 / 2) * u) * conj (gDens p p' p'' u)) +
      conj (gA p p') + cexp (-(s - 1 / 2) * ℓ) * conj (gB ℓ p p') :=
    funext (Elow_green p p' p'' hp' hcp hcp' hcp'')
  rw [e]
  have hq : Continuous fun u => conj (gDens p p' p'' u) := by unfold gDens; fun_prop
  exact ((differentiable_laplace hℓ hq).add_const _).add (by fun_prop)

end Riemann.Deep.LeanNBRate

