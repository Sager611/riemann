import Riemann.Deep.JerbyFals.Lean.Defs

/-!
# θ is the continuous branch of arg Γ(1/4 + it/2) − (t/2) log π with θ(0) = 0

* `gamma_phase`: Γ(1/4 + it/2) = |Γ(1/4 + it/2)| e^{i(θ(t) + (t/2) log π)} for every real t.
* `theta_unique`: every continuous φ with φ(0) = 0 and the same phase property is θ.

So `theta` is the Riemann–Siegel θ (Edwards §6.5), and the function certify_cx.py evaluates
as Im lgamma(1/4 + it/2) − (t/2) log π (Arb's principal log Γ: continuous on Re z > 0, real at
t = 0, and its imaginary part is an argument of Γ).
-/

open Complex Set
open scoped Real

namespace JerbyLean

noncomputable section

lemma quarter_line_re_pos (t : ℝ) : 0 < (1 / 4 + t / 2 * I : ℂ).re := by
  rw [re_quarter_line]; norm_num

lemma continuous_digamma_line : Continuous (fun u : ℝ => digamma (1 / 4 + u / 2 * I)) :=
  continuousOn_digamma_rhp.comp_continuous (by fun_prop)
    (fun u => by show 0 < (1 / 4 + (u : ℂ) / 2 * I).re; exact quarter_line_re_pos u)

lemma hasDerivAt_quarter_line (t : ℝ) :
    HasDerivAt (fun u : ℝ => (1 / 4 + u / 2 * I : ℂ)) (I / 2) t :=
  ((((hasDerivAt_id' (x := t)).ofReal_comp).div_const 2).mul_const I |>.const_add
    (1 / 4 : ℂ)).congr_deriv (by push_cast; ring)

lemma deriv_Gamma_eq {s : ℂ} (hs : 0 < s.re) : deriv Gamma s = digamma s * Gamma s := by
  rw [digamma_def, logDeriv_apply, div_mul_cancel₀ _ (Gamma_ne_zero_of_re_pos hs)]

lemma hasDerivAt_Gamma_line (t : ℝ) :
    HasDerivAt (fun u : ℝ => Gamma (1 / 4 + u / 2 * I))
      (digamma (1 / 4 + t / 2 * I) * Gamma (1 / 4 + t / 2 * I) * (I / 2)) t := by
  have hpos := quarter_line_re_pos t
  have hd := (differentiableAt_Gamma _ (ne_neg_nat_of_re_pos hpos)).hasDerivAt
  rw [deriv_Gamma_eq hpos] at hd
  exact hd.comp t (hasDerivAt_quarter_line t)

/-- θ(t) + (t/2) log π, the claimed continuous argument of Γ(1/4 + it/2). -/
def argGamma (t : ℝ) : ℝ := theta t + t / 2 * Real.log π

lemma hasDerivAt_argGamma (t : ℝ) :
    HasDerivAt argGamma ((digamma (1 / 4 + t / 2 * I)).re / 2) t := by
  have h := (hasDerivAt_theta t).add
    (((hasDerivAt_id' (x := t)).div_const 2).mul_const (Real.log π))
  refine h.congr_deriv ?_
  unfold thetaDeriv; ring

/-- −½ Im ψ(1/4 + it/2) = d/dt log |Γ(1/4 + it/2)|. -/
def phaseRate (u : ℝ) : ℝ := -(digamma (1 / 4 + u / 2 * I)).im / 2

lemma continuous_phaseRate : Continuous phaseRate :=
  ((continuous_im.comp continuous_digamma_line).neg).div_const 2

def phaseLog (t : ℝ) : ℝ := ∫ u in (0 : ℝ)..t, phaseRate u

lemma hasDerivAt_phaseLog (t : ℝ) : HasDerivAt phaseLog (phaseRate t) t :=
  intervalIntegral.integral_hasDerivAt_right (continuous_phaseRate.intervalIntegrable _ _)
    continuous_phaseRate.aestronglyMeasurable.stronglyMeasurableAtFilter
    continuous_phaseRate.continuousAt

/-- Γ(1/4 + it/2) e^{−i(θ(t) + (t/2) log π)} e^{−phaseLog t}; its derivative vanishes. -/
def phaseK (t : ℝ) : ℂ :=
  Gamma (1 / 4 + t / 2 * I) * exp (-(I * argGamma t)) * exp (-(phaseLog t : ℂ))

lemma phase_alg (ψ G e1 e2 : ℂ) :
    (ψ * G * (I / 2) * e1 + G * (e1 * -(I * ((ψ.re / 2 : ℝ) : ℂ)))) * e2 +
      G * e1 * (e2 * -((-ψ.im / 2 : ℝ) : ℂ)) = 0 := by
  have hψ : ψ = (ψ.re : ℂ) + ψ.im * I := (re_add_im ψ).symm
  push_cast
  linear_combination (G * e1 * e2 * (I / 2)) * hψ + (G * e1 * e2 * ((ψ.im : ℂ) / 2)) * I_sq

lemma hasDerivAt_phaseK (x : ℝ) : HasDerivAt phaseK 0 x := by
  have h1 := hasDerivAt_Gamma_line x
  have h2 := (((hasDerivAt_argGamma x).ofReal_comp.const_mul I).neg).cexp
  have h3 := ((hasDerivAt_phaseLog x).ofReal_comp.neg).cexp
  exact ((h1.mul h2).mul h3).congr_deriv (phase_alg _ _ _ _)

lemma phaseK_const (t : ℝ) : phaseK t = phaseK 0 :=
  is_const_of_deriv_eq_zero (fun x => (hasDerivAt_phaseK x).differentiableAt)
    (fun x => (hasDerivAt_phaseK x).deriv) t 0

lemma phaseK_zero : phaseK 0 = (Real.Gamma (1 / 4) : ℂ) := by
  simp only [phaseK, argGamma, theta_zero, phaseLog, intervalIntegral.integral_same]
  rw [← Gamma_ofReal]
  push_cast
  simp

/-- **Phase identity**: θ(t) + (t/2) log π is an argument of Γ(1/4 + it/2). -/
theorem gamma_phase (t : ℝ) :
    Gamma (1 / 4 + t / 2 * I) =
      ‖Gamma (1 / 4 + t / 2 * I)‖ * exp (I * ((theta t + t / 2 * Real.log π : ℝ) : ℂ)) := by
  have hK := phaseK_const t
  rw [phaseK_zero] at hK
  unfold phaseK at hK
  have e1 : exp (-(I * (argGamma t : ℂ))) * exp (I * (argGamma t : ℂ)) = 1 := by
    rw [← Complex.exp_add, neg_add_cancel, Complex.exp_zero]
  have e2 : exp (-(phaseLog t : ℂ)) * exp (phaseLog t : ℂ) = 1 := by
    rw [← Complex.exp_add, neg_add_cancel, Complex.exp_zero]
  have hG : Gamma (1 / 4 + t / 2 * I) =
      ((Real.Gamma (1 / 4) * Real.exp (phaseLog t) : ℝ) : ℂ) * exp (I * (argGamma t : ℂ)) := by
    push_cast
    linear_combination (exp (I * (argGamma t : ℂ)) * exp (phaseLog t : ℂ)) * hK -
      (Gamma (1 / 4 + t / 2 * I) * exp (-(phaseLog t : ℂ)) * exp (phaseLog t : ℂ)) * e1 -
      Gamma (1 / 4 + t / 2 * I) * e2
  have hρ : 0 < Real.Gamma (1 / 4) * Real.exp (phaseLog t) :=
    mul_pos (Real.Gamma_pos_of_pos (by norm_num)) (Real.exp_pos _)
  have hu : ‖exp (I * (argGamma t : ℂ))‖ = 1 := by
    rw [mul_comm]; exact norm_exp_ofReal_mul_I _
  have hn : ‖Gamma (1 / 4 + t / 2 * I)‖ = Real.Gamma (1 / 4) * Real.exp (phaseLog t) := by
    rw [hG, norm_mul, hu, mul_one, norm_real, Real.norm_of_nonneg hρ.le]
  rw [hn]
  exact hG

/-- **θ is the continuous branch**: a continuous φ with φ(0) = 0 whose φ(t) + (t/2) log π is an
argument of Γ(1/4 + it/2) for every t is θ. -/
theorem theta_unique (φ : ℝ → ℝ) (hφ : Continuous φ) (h0 : φ 0 = 0)
    (hphase : ∀ t : ℝ, Gamma (1 / 4 + t / 2 * I) =
      ‖Gamma (1 / 4 + t / 2 * I)‖ * exp (I * ((φ t + t / 2 * Real.log π : ℝ) : ℂ))) :
    φ = theta := by
  have hint : ∀ t, ∃ k : ℤ, φ t - theta t = 2 * π * k := by
    intro t
    have hG0 : Gamma (1 / 4 + t / 2 * I) ≠ 0 := Gamma_ne_zero_of_re_pos (quarter_line_re_pos t)
    have hn : ((‖Gamma (1 / 4 + t / 2 * I)‖ : ℝ) : ℂ) ≠ 0 := by
      exact_mod_cast norm_ne_zero_iff.2 hG0
    have heq : exp (I * ((φ t + t / 2 * Real.log π : ℝ) : ℂ)) =
        exp (I * ((theta t + t / 2 * Real.log π : ℝ) : ℂ)) :=
      mul_left_cancel₀ hn ((hphase t).symm.trans (gamma_phase t))
    obtain ⟨k, hk⟩ := Complex.exp_eq_exp_iff_exists_int.1 heq
    refine ⟨k, ?_⟩
    have := congrArg Complex.im hk
    simp at this
    linarith
  have hf : Continuous (fun u => φ u - theta u) := hφ.sub continuous_theta
  funext t
  obtain ⟨k, hk⟩ := hint t
  by_contra hne
  have hk0 : k ≠ 0 := by
    rintro rfl
    apply hne
    simp at hk
    linarith
  have hπ := Real.pi_pos
  have h00 : φ 0 - theta 0 = 0 := by rw [h0, theta_zero]; ring
  -- the value v = ±π lies between f 0 = 0 and f t = 2πk, so f c = v = 2πm for some c, m
  have key : ∀ v : ℝ, v ∈ uIcc (φ 0 - theta 0) (φ t - theta t) →
      ∃ m : ℤ, v = 2 * π * m := by
    intro v hv
    obtain ⟨c, -, hc⟩ := intermediate_value_uIcc hf.continuousOn hv
    obtain ⟨m, hm⟩ := hint c
    exact ⟨m, hc.symm.trans hm⟩
  rcases lt_or_gt_of_ne hk0 with hneg | hpos
  · have hk1 : (k : ℝ) ≤ -1 := by exact_mod_cast (show k ≤ -1 by omega)
    obtain ⟨m, hm⟩ := key (-π) (by
      rw [h00, hk]; exact mem_uIcc_of_ge (by nlinarith) (by linarith))
    have h1 : (2 * m : ℝ) = -1 := by
      have := mul_left_cancel₀ hπ.ne' (show π * (2 * m) = π * (-1) by linarith)
      linarith
    have : (2 * m : ℤ) = -1 := by exact_mod_cast h1
    omega
  · have hk1 : (1 : ℝ) ≤ k := by exact_mod_cast (show 1 ≤ k by omega)
    obtain ⟨m, hm⟩ := key π (by
      rw [h00, hk]; exact mem_uIcc_of_le hπ.le (by nlinarith))
    have h1 : (2 * m : ℝ) = 1 := by
      have := mul_left_cancel₀ hπ.ne' (show π * (2 * m) = π * 1 by linarith)
      linarith
    have : (2 * m : ℤ) = 1 := by exact_mod_cast h1
    omega

end

end JerbyLean
