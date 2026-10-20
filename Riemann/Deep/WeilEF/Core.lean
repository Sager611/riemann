/-
Riemann/Deep/WeilEF/Core.lean — the analytic core of the sufficiency half of Weil's criterion.

Abstract setting (no zeta): weights `w i` (absolutely summable), exponents `p i` with `|Re p i| ≤ 1/2`,
the exponential sum `D τ = Σ w i e^{p i τ}`.  If `D` is bounded on `τ ≥ 0` then every exponent
`p i0` with `Re p i0 > 0` that is "isolated from the right" (no exponent with larger real part
within ordinate distance `η`) carries zero weight.

Proof: Laplace transform `L s = ∫_0^∞ D τ e^{-sτ} dτ` is holomorphic on `Re s > 0` and bounded by
`M / Re s`; on `Re s > 1/2` it equals `R s = Σ w i / (s - p i)` (Fubini); `R` is holomorphic on the
rectangle `V = {Re s > Re p i0, |Im s - Im p i0| < η}`; identity theorem on the convex `V`; then
`x · L (p i0 + x) → 0` while `x · R (p i0 + x) → w i0` as `x → 0⁺` (dominated convergence for sums).
-/
import Mathlib.Analysis.SpecialFunctions.ImproperIntegrals
import Mathlib.Analysis.Complex.LocallyUniformLimit
import Mathlib.Analysis.Complex.CauchyIntegral
import Mathlib.Analysis.Analytic.Uniqueness
import Mathlib.Analysis.Calculus.ParametricIntegral
import Mathlib.Analysis.Normed.Group.Tannery
import Mathlib.MeasureTheory.Integral.DominatedConvergence
import Mathlib.Analysis.Complex.Convex

noncomputable section

open Complex MeasureTheory Filter Topology Set

namespace Riemann.Deep.WeilEF.Core

variable {ι : Type*} (p w : ι → ℂ)

/-- The exponential sum `D τ = Σ w i e^{p i τ}`. -/
def expSum (τ : ℝ) : ℂ := ∑' i, w i * cexp (p i * τ)

/-- Laplace transform on `(0, ∞)`. -/
def lap (F : ℝ → ℂ) (s : ℂ) : ℂ := ∫ τ in Ioi (0 : ℝ), F τ * cexp (-s * τ)

/-- The partial-fraction sum `Σ w i / (s - p i)`. -/
def resSum (s : ℂ) : ℂ := ∑' i, w i / (s - p i)

lemma norm_cexp_mul_ofReal (z : ℂ) (τ : ℝ) : ‖cexp (z * τ)‖ = Real.exp (z.re * τ) := by
  rw [Complex.norm_exp]; simp

variable {p w}

lemma term_norm_le (hp : ∀ i, |(p i).re| ≤ 1 / 2) (i : ι) (τ : ℝ) :
    ‖w i * cexp (p i * τ)‖ ≤ ‖w i‖ * Real.exp |τ| := by
  rw [norm_mul, norm_cexp_mul_ofReal]
  refine mul_le_mul_of_nonneg_left (Real.exp_le_exp.mpr ?_) (norm_nonneg _)
  have h1 : (p i).re * τ ≤ |(p i).re| * |τ| := by rw [← abs_mul]; exact le_abs_self _
  have h2 : |(p i).re| * |τ| ≤ (1 / 2) * |τ| := mul_le_mul_of_nonneg_right (hp i) (abs_nonneg _)
  linarith [abs_nonneg τ]

lemma continuous_expSum (hw : Summable (fun i => ‖w i‖)) (hp : ∀ i, |(p i).re| ≤ 1 / 2) :
    Continuous (expSum p w) := by
  rw [continuous_iff_continuousAt]
  intro τ0
  have hon : ContinuousOn (fun τ : ℝ => ∑' i, w i * cexp (p i * τ)) (Ioo (τ0 - 1) (τ0 + 1)) := by
    refine continuousOn_tsum (u := fun i => ‖w i‖ * Real.exp (|τ0| + 1))
      (fun i => (by fun_prop : Continuous fun τ : ℝ => w i * cexp (p i * τ)).continuousOn)
      (hw.mul_right _) (fun i τ hτ => (term_norm_le hp i τ).trans ?_)
    refine mul_le_mul_of_nonneg_left (Real.exp_le_exp.mpr ?_) (norm_nonneg _)
    have h1 := abs_sub_abs_le_abs_sub τ τ0
    have h2 : |τ - τ0| < 1 := abs_sub_lt_iff.mpr ⟨by linarith [hτ.2], by linarith [hτ.1]⟩
    linarith
  exact hon.continuousAt (Ioo_mem_nhds (by linarith) (by linarith))

lemma norm_lap_le {F : ℝ → ℂ} {M : ℝ} (hM : ∀ τ : ℝ, 0 ≤ τ → ‖F τ‖ ≤ M)
    {s : ℂ} (hs : 0 < s.re) : ‖lap F s‖ ≤ M / s.re := by
  have hg : Integrable (fun τ : ℝ => M * Real.exp (-s.re * τ)) (volume.restrict (Ioi 0)) :=
    (integrableOn_exp_mul_Ioi (neg_lt_zero.mpr hs) 0).const_mul M
  have hb : ∀ᵐ τ ∂(volume.restrict (Ioi (0 : ℝ))),
      ‖F τ * cexp (-s * τ)‖ ≤ M * Real.exp (-s.re * τ) := by
    refine (ae_restrict_iff' measurableSet_Ioi).mpr (Eventually.of_forall fun τ hτ => ?_)
    rw [norm_mul, norm_cexp_mul_ofReal, neg_re]
    exact mul_le_mul_of_nonneg_right (hM τ (le_of_lt hτ)) (Real.exp_pos _).le
  refine (norm_integral_le_of_norm_le hg hb).trans (le_of_eq ?_)
  rw [integral_const_mul, integral_exp_mul_Ioi (neg_lt_zero.mpr hs)]
  field_simp
  simp

lemma lap_eq_resSum [Countable ι] (hw : Summable (fun i => ‖w i‖))
    (hp : ∀ i, |(p i).re| ≤ 1 / 2) {s : ℂ} (hs : 1 / 2 < s.re) :
    lap (expSum p w) s = resSum p w s := by
  have hre : ∀ i, (p i - s).re < 0 := fun i => by
    rw [sub_re]; linarith [(abs_le.mp (hp i)).2]
  set F : ι → ℝ → ℂ := fun i τ => w i * cexp ((p i - s) * τ) with hF
  have hint : ∀ i, Integrable (F i) (volume.restrict (Ioi 0)) := fun i =>
    (integrableOn_exp_mul_complex_Ioi (hre i) 0).const_mul (w i)
  have hnorm : ∀ i, ∫ τ in Ioi (0 : ℝ), ‖F i τ‖ = ‖w i‖ / (s.re - (p i).re) := by
    intro i
    have hlt : (p i).re - s.re < 0 := by have := hre i; rwa [sub_re] at this
    have : (fun τ : ℝ => ‖F i τ‖) = fun τ => ‖w i‖ * Real.exp (((p i).re - s.re) * τ) := by
      funext τ; rw [hF, norm_mul, norm_cexp_mul_ofReal, sub_re]
    rw [this, integral_const_mul, integral_exp_mul_Ioi hlt, mul_zero, Real.exp_zero,
      show (p i).re - s.re = -(s.re - (p i).re) by ring, div_neg, neg_div, neg_neg, mul_one_div]
  have hsum : Summable (fun i => ∫ τ in Ioi (0 : ℝ), ‖F i τ‖) := by
    simp_rw [hnorm]
    refine Summable.of_nonneg_of_le (fun i => div_nonneg (norm_nonneg _) ?_) (fun i => ?_)
      (hw.div_const (s.re - 1 / 2))
    · linarith [(abs_le.mp (hp i)).2]
    · exact div_le_div_of_nonneg_left (norm_nonneg _) (by linarith)
        (by linarith [(abs_le.mp (hp i)).2])
  have key := integral_tsum_of_summable_integral_norm hint hsum
  unfold lap resSum
  have hfun : (fun τ : ℝ => expSum p w τ * cexp (-s * τ)) = fun τ => ∑' i, F i τ := by
    funext τ
    rw [expSum, ← tsum_mul_right]
    congr 1; funext i
    rw [hF]
    simp only
    rw [mul_assoc, ← Complex.exp_add]
    congr 2; ring
  rw [hfun, ← key]
  congr 1; funext i
  rw [hF]
  simp only
  rw [integral_const_mul, integral_exp_mul_complex_Ioi (hre i)]
  have hne : p i - s ≠ 0 := fun h => by have := hre i; rw [h] at this; simp at this
  have hne' : s - p i ≠ 0 := fun h => hne (by rw [← neg_sub, h, neg_zero])
  simp only [ofReal_zero, mul_zero, Complex.exp_zero]
  field_simp
  ring

lemma differentiableOn_lap {F : ℝ → ℂ} (hF : Continuous F) {M : ℝ}
    (hM : ∀ τ : ℝ, 0 ≤ τ → ‖F τ‖ ≤ M) : DifferentiableOn ℂ (lap F) {s : ℂ | 0 < s.re} := by
  intro s0 hs0
  have ha : 0 < s0.re := hs0
  set a := s0.re with ha_def
  refine DifferentiableAt.differentiableWithinAt ?_
  have hM0 : 0 ≤ M := (norm_nonneg _).trans (hM 0 le_rfl)
  have key := hasDerivAt_integral_of_dominated_loc_of_deriv_le (μ := volume.restrict (Ioi (0 : ℝ)))
    (F := fun (s : ℂ) (τ : ℝ) => F τ * cexp (-s * τ))
    (F' := fun (s : ℂ) (τ : ℝ) => F τ * (cexp (-s * τ) * -(τ : ℂ)))
    (x₀ := s0) (s := Metric.ball s0 (a / 2))
    (bound := fun τ => M * (4 / a) * Real.exp (-(a / 4) * τ))
    (Metric.ball_mem_nhds _ (by positivity))
    (Eventually.of_forall fun s => (hF.mul (by fun_prop)).aestronglyMeasurable)
    ?_ ((hF.mul (by fun_prop)).aestronglyMeasurable) ?_
    ((integrableOn_exp_mul_Ioi (by linarith : -(a / 4) < 0) 0).const_mul (M * (4 / a))) ?_
  · exact key.2.differentiableAt
  · -- integrability at s0
    refine Integrable.mono' ((integrableOn_exp_mul_Ioi (neg_lt_zero.mpr ha) 0).const_mul M)
      ((hF.mul (by fun_prop)).aestronglyMeasurable) ?_
    refine (ae_restrict_iff' measurableSet_Ioi).mpr (Eventually.of_forall fun τ hτ => ?_)
    show ‖F τ * cexp (-s0 * τ)‖ ≤ M * Real.exp (-a * τ)
    rw [norm_mul, norm_cexp_mul_ofReal, neg_re]
    exact mul_le_mul_of_nonneg_right (hM τ (le_of_lt hτ)) (Real.exp_pos _).le
  · refine (ae_restrict_iff' measurableSet_Ioi).mpr (Eventually.of_forall fun τ hτ s hs => ?_)
    have hτ0 : (0 : ℝ) < τ := hτ
    have hsre : a / 2 < s.re := by
      have h1 : |s.re - s0.re| ≤ ‖s - s0‖ := by rw [← sub_re]; exact abs_re_le_norm _
      have h2 : ‖s - s0‖ < a / 2 := by rw [← dist_eq_norm]; exact hs
      have := (abs_lt.mp (lt_of_le_of_lt h1 h2)).1
      linarith
    show ‖F τ * (cexp (-s * τ) * -(τ : ℂ))‖ ≤ M * (4 / a) * Real.exp (-(a / 4) * τ)
    rw [norm_mul, norm_mul, norm_cexp_mul_ofReal, neg_re, norm_neg, Complex.norm_real,
      Real.norm_eq_abs, abs_of_pos hτ0]
    have hexp1 : Real.exp (-s.re * τ) ≤ Real.exp (-(a / 2) * τ) :=
      Real.exp_le_exp.mpr (by nlinarith)
    have hτle : τ ≤ (4 / a) * Real.exp ((a / 4) * τ) := by
      have := Real.add_one_le_exp ((a / 4) * τ)
      rw [div_mul_eq_mul_div, le_div_iff₀ ha]
      nlinarith
    have hprod : Real.exp ((a / 4) * τ) * Real.exp (-(a / 2) * τ) = Real.exp (-(a / 4) * τ) := by
      rw [← Real.exp_add]; congr 1; ring
    calc ‖F τ‖ * (Real.exp (-s.re * τ) * τ)
        ≤ M * (Real.exp (-(a / 2) * τ) * ((4 / a) * Real.exp ((a / 4) * τ))) := by
          apply mul_le_mul (hM τ hτ0.le) _ (by positivity) hM0
          exact mul_le_mul hexp1 hτle hτ0.le (Real.exp_pos _).le
      _ = M * (4 / a) * Real.exp (-(a / 4) * τ) := by
          rw [← hprod]; ring
  · refine Eventually.of_forall fun τ s _ => ?_
    have h1 : HasDerivAt (fun s : ℂ => -s * (τ : ℂ)) (-(τ : ℂ)) s := by
      simpa using ((hasDerivAt_id s).neg.mul_const (τ : ℂ))
    exact h1.cexp.const_mul (F τ)

/-- Distance from the exponents, on a shrunken rectangle. -/
lemma le_norm_sub {σ0 t0 η ε : ℝ} (hsep : ∀ i, σ0 < (p i).re → η ≤ |(p i).im - t0|) {s : ℂ}
    (hre : σ0 + ε < s.re) (him1 : t0 - η + ε < s.im) (him2 : s.im < t0 + η - ε) (i : ι) :
    ε ≤ ‖s - p i‖ := by
  by_cases h : σ0 < (p i).re
  · have h1 := hsep i h
    have h2 : ε ≤ |(s - p i).im| := by
      rw [sub_im]
      have h4 : |s.im - t0| < η - ε := abs_sub_lt_iff.mpr ⟨by linarith, by linarith⟩
      have h5 := abs_sub_le ((p i).im) (s.im) t0
      have h6 : |(p i).im - s.im| = |s.im - (p i).im| := abs_sub_comm _ _
      linarith
    exact h2.trans (abs_im_le_norm _)
  · have h2 : ε ≤ (s - p i).re := by rw [sub_re]; linarith
    exact h2.trans ((le_abs_self _).trans (abs_re_le_norm _))

lemma differentiableOn_resSum (hw : Summable (fun i => ‖w i‖)) {σ0 t0 η : ℝ}
    (hsep : ∀ i, σ0 < (p i).re → η ≤ |(p i).im - t0|) :
    DifferentiableOn ℂ (resSum p w) {s : ℂ | σ0 < s.re ∧ t0 - η < s.im ∧ s.im < t0 + η} := by
  rintro s ⟨h1, h2, h3⟩
  set ε := min (s.re - σ0) (min (s.im - (t0 - η)) (t0 + η - s.im)) / 2 with hε
  have hε0 : 0 < ε := by
    rw [hε]; apply half_pos; exact lt_min (by linarith) (lt_min (by linarith) (by linarith))
  have hεa : ε < s.re - σ0 := by
    have := min_le_left (s.re - σ0) (min (s.im - (t0 - η)) (t0 + η - s.im)); linarith
  have hεb : ε < s.im - (t0 - η) := by
    have := min_le_right (s.re - σ0) (min (s.im - (t0 - η)) (t0 + η - s.im))
    have := min_le_left (s.im - (t0 - η)) (t0 + η - s.im); linarith
  have hεc : ε < t0 + η - s.im := by
    have := min_le_right (s.re - σ0) (min (s.im - (t0 - η)) (t0 + η - s.im))
    have := min_le_right (s.im - (t0 - η)) (t0 + η - s.im); linarith
  set U := {z : ℂ | σ0 + ε < z.re ∧ t0 - η + ε < z.im ∧ z.im < t0 + η - ε} with hU
  have hUo : IsOpen U :=
    (isOpen_lt continuous_const continuous_re).inter
      ((isOpen_lt continuous_const continuous_im).inter (isOpen_lt continuous_im continuous_const))
  have hsU : s ∈ U := ⟨by linarith, by linarith, by linarith⟩
  have hd : DifferentiableOn ℂ (resSum p w) U := by
    refine differentiableOn_tsum_of_summable_norm (u := fun i => ‖w i‖ / ε) (hw.div_const ε)
      (fun i => ?_) hUo (fun i z hz => ?_)
    · refine (differentiableOn_const _).div (differentiableOn_id.sub_const _) (fun z hz => ?_)
      have := le_norm_sub hsep hz.1 hz.2.1 hz.2.2 i
      intro h0; rw [h0, norm_zero] at this; linarith
    · rw [norm_div]
      exact div_le_div_of_nonneg_left (norm_nonneg _) hε0 (le_norm_sub hsep hz.1 hz.2.1 hz.2.2 i)
  exact (hd.differentiableAt (hUo.mem_nhds hsU)).differentiableWithinAt

/-- **Core.** Bounded exponential sums carry no weight at right-isolated exponents with `Re > 0`. -/
theorem weight_eq_zero [Countable ι] (hw : Summable (fun i => ‖w i‖))
    (hp : ∀ i, |(p i).re| ≤ 1 / 2) {M : ℝ} (hD : ∀ τ : ℝ, 0 ≤ τ → ‖expSum p w τ‖ ≤ M)
    (i0 : ι) (hσ : 0 < (p i0).re) {η : ℝ} (hη : 0 < η)
    (hsep : ∀ i, (p i0).re < (p i).re → η ≤ |(p i).im - (p i0).im|)
    (hinj : ∀ i, p i = p i0 → i = i0) : w i0 = 0 := by
  classical
  set σ0 := (p i0).re with hσ0
  set t0 := (p i0).im with ht0
  set s0 := p i0 with hs0
  set V := {s : ℂ | σ0 < s.re ∧ t0 - η < s.im ∧ s.im < t0 + η} with hV
  have hVo : IsOpen V :=
    (isOpen_lt continuous_const continuous_re).inter
      ((isOpen_lt continuous_const continuous_im).inter (isOpen_lt continuous_im continuous_const))
  have hVc : Convex ℝ V := by
    have : V = {s : ℂ | σ0 < s.re} ∩ ({s : ℂ | t0 - η < s.im} ∩ {s : ℂ | s.im < t0 + η}) := by
      ext s; simp [hV]
    rw [this]
    exact (convex_halfSpace_re_gt σ0).inter ((convex_halfSpace_im_gt _).inter
      (convex_halfSpace_im_lt _))
  have hM0 : 0 ≤ M := (norm_nonneg _).trans (hD 0 le_rfl)
  have hcont := continuous_expSum hw hp
  have hLd : DifferentiableOn ℂ (lap (expSum p w)) V :=
    (differentiableOn_lap hcont hD).mono fun s hs => lt_trans hσ hs.1
  have hRd : DifferentiableOn ℂ (resSum p w) V := differentiableOn_resSum hw hsep
  -- identity theorem on V
  set z1 : ℂ := ⟨σ0 + 1, t0⟩ with hz1
  have hz1V : z1 ∈ V := ⟨by simp [hz1], by simp [hz1]; linarith, by simp [hz1]; linarith⟩
  have hev : lap (expSum p w) =ᶠ[𝓝 z1] resSum p w := by
    have hop : IsOpen {s : ℂ | 1 / 2 < s.re} := isOpen_lt continuous_const continuous_re
    have hz : z1 ∈ {s : ℂ | 1 / 2 < s.re} := by
      show (1 : ℝ) / 2 < z1.re; simp [hz1]; linarith
    exact Filter.eventually_of_mem (hop.mem_nhds hz) fun s hs => lap_eq_resSum hw hp hs
  have hEq : EqOn (lap (expSum p w)) (resSum p w) V :=
    (hLd.analyticOnNhd hVo).eqOn_of_preconnected_of_eventuallyEq (hRd.analyticOnNhd hVo)
      hVc.isPreconnected hz1V hev
  -- the approach s0 + x, x → 0⁺
  have hmem : ∀ x : ℝ, 0 < x → s0 + (x : ℂ) ∈ V := fun x hx =>
    ⟨by simp [hσ0, hs0]; linarith, by simp [ht0, hs0]; linarith, by simp [ht0, hs0]; linarith⟩
  have hLlim : Tendsto (fun x : ℝ => (x : ℂ) * lap (expSum p w) (s0 + x)) (𝓝[>] 0) (𝓝 0) := by
    have hlim0 : Tendsto (fun x : ℝ => x * (M / σ0)) (𝓝[>] 0) (𝓝 0) := by
      have : Tendsto (fun x : ℝ => x * (M / σ0)) (𝓝 0) (𝓝 (0 * (M / σ0))) :=
        (continuous_id.mul continuous_const).tendsto 0
      rw [zero_mul] at this
      exact this.mono_left nhdsWithin_le_nhds
    refine squeeze_zero_norm' ?_ hlim0
    filter_upwards [self_mem_nhdsWithin] with x hx
    have hx0 : (0 : ℝ) < x := hx
    have hre : 0 < (s0 + (x : ℂ)).re := by simp [hs0]; linarith
    rw [norm_mul, Complex.norm_real, Real.norm_eq_abs, abs_of_pos hx0]
    refine mul_le_mul_of_nonneg_left ((norm_lap_le hD hre).trans ?_) hx0.le
    have : (s0 + (x : ℂ)).re = σ0 + x := by simp [hs0, hσ0]
    rw [this]
    exact div_le_div_of_nonneg_left hM0 hσ (by linarith)
  have hRlim : Tendsto (fun x : ℝ => (x : ℂ) * resSum p w (s0 + x)) (𝓝[>] 0) (𝓝 (w i0)) := by
    have hrw : (fun x : ℝ => (x : ℂ) * resSum p w (s0 + x))
        = fun x : ℝ => ∑' i, w i * ((x : ℂ) / (s0 + x - p i)) := by
      funext x
      rw [resSum, ← tsum_mul_left]
      congr 1; funext i; ring
    rw [hrw]
    have hg : (∑' i, if i = i0 then w i0 else 0) = w i0 := tsum_ite_eq i0 (fun _ => w i0)
    rw [← hg]
    refine tendsto_tsum_of_dominated_convergence (bound := fun i => ‖w i‖ * (1 + 1 / η))
      (hw.mul_right _) (fun i => ?_) ?_
    · by_cases hi : i = i0
      · subst hi
        simp only [ite_true]
        refine tendsto_const_nhds.congr' ?_
        filter_upwards [self_mem_nhdsWithin] with x hx
        have hx0 : (x : ℂ) ≠ 0 := ofReal_ne_zero.mpr (ne_of_gt hx)
        rw [show s0 + (x : ℂ) - p i = x by rw [hs0]; ring, div_self hx0, mul_one]
      · simp only [hi, ite_false]
        have hne : s0 - p i ≠ 0 := fun h => hi (hinj i (by rw [hs0] at h; linear_combination -h))
        have hc : ContinuousAt (fun x : ℝ => w i * ((x : ℂ) / (s0 + x - p i))) 0 := by
          refine continuousAt_const.mul (ContinuousAt.div (by fun_prop) (by fun_prop) ?_)
          simpa using hne
        have := hc.tendsto
        simp only [ofReal_zero, zero_div, mul_zero] at this
        exact this.mono_left nhdsWithin_le_nhds
    · filter_upwards [Ioo_mem_nhdsGT (show (0 : ℝ) < 1 by norm_num)] with x hx i
      have hx0 : 0 < x := hx.1
      rw [norm_mul, norm_div, Complex.norm_real, Real.norm_eq_abs, abs_of_pos hx0]
      refine mul_le_mul_of_nonneg_left ?_ (norm_nonneg _)
      by_cases h : σ0 < (p i).re
      · have h1 := hsep i h
        have h2 : η ≤ ‖s0 + (x : ℂ) - p i‖ := by
          refine h1.trans ((le_of_eq ?_).trans (abs_im_le_norm _))
          simp [hs0, ht0, abs_sub_comm]
        have h3 : x / ‖s0 + (x : ℂ) - p i‖ ≤ x / η := div_le_div_of_nonneg_left hx0.le hη h2
        have h4 : x / η ≤ 1 / η := div_le_div_of_nonneg_right hx.2.le hη.le
        have h5 : 0 ≤ 1 / η := by positivity
        linarith
      · have h2 : x ≤ ‖s0 + (x : ℂ) - p i‖ := by
          refine le_trans ?_ ((le_abs_self _).trans (abs_re_le_norm _))
          simp [hs0]; push Not at h; linarith
        have h3 : x / ‖s0 + (x : ℂ) - p i‖ ≤ 1 := by
          rw [div_le_one (lt_of_lt_of_le hx0 h2)]; exact h2
        have h5 : 0 ≤ 1 / η := by positivity
        linarith
  have hLR : (fun x : ℝ => (x : ℂ) * lap (expSum p w) (s0 + x))
      =ᶠ[𝓝[>] 0] fun x : ℝ => (x : ℂ) * resSum p w (s0 + x) := by
    filter_upwards [self_mem_nhdsWithin] with x hx
    rw [hEq (hmem x hx)]
  exact (tendsto_nhds_unique (hLlim.congr' hLR) hRlim).symm

end Riemann.Deep.WeilEF.Core
