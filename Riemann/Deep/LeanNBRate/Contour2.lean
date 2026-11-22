import Riemann.Deep.LeanNBRate.Contour1

/-!
# U1 contour argument, part 2: bounds on the right, Fubini for G, the residue at s
-/

open Complex Filter Topology Set MeasureTheory intervalIntegral
open scoped Real Interval

noncomputable section

namespace Riemann.Deep.LeanNBRate

variable {H : ℂ → ℂ} {M β D : ℝ} {N : ℕ}

theorem continuousAt_Fq (hH : HBound H M β D N) {z : ℂ} (hz1 : z ≠ 1) (hz0 : riemannZeta z ≠ 0) :
    ContinuousAt (Fq H) z :=
  (hH.diff z).continuousAt.div (differentiableAt_riemannZeta hz1).continuousAt hz0

theorem differentiableAt_Fq (hH : HBound H M β D N) {z : ℂ} (hz1 : z ≠ 1) (hz0 : riemannZeta z ≠ 0) :
    DifferentiableAt ℂ (Fq H) z :=
  (hH.diff z).div (differentiableAt_riemannZeta hz1) hz0

theorem ne_one_of_two_le {z : ℂ} (hz : 2 ≤ z.re) : z ≠ 1 := fun h => by rw [h] at hz; norm_num at hz

theorem one_add_norm_le (σ y : ℝ) : 1 + ‖(σ : ℂ) + y * I‖ ≤ (1 + |σ|) * (1 + |y|) := by
  have : ‖(σ : ℂ) + y * I‖ ≤ |σ| + |y| := by
    calc ‖(σ : ℂ) + y * I‖ ≤ ‖(σ : ℂ)‖ + ‖(y : ℂ) * I‖ := norm_add_le _ _
      _ = |σ| + |y| := by simp
  nlinarith [abs_nonneg σ, abs_nonneg y]

/-- Right region bound for F = H/zeta. -/
theorem norm_Fq_right (hH : HBound H M β D N) {z : ℂ} (hz : 2 ≤ z.re) :
    ‖Fq H z‖ ≤ 3 * M * (D / (D + ‖z - 1 / 2‖)) ^ N := by
  rw [Fq, div_eq_mul_inv, norm_mul]
  have h1 := hH.bound z
  have hmax : max 0 (1 / 2 - z.re) = 0 := max_eq_left (by linarith)
  rw [hmax, mul_zero, Real.exp_zero, mul_one] at h1
  have h2 := norm_inv_zeta_le_three hz
  have hM := hH.M_nonneg
  have hD := hH.D_pos
  calc ‖H z‖ * ‖(riemannZeta z)⁻¹‖ ≤ (M * (D / (D + ‖z - 1 / 2‖)) ^ N) * 3 :=
        mul_le_mul h1 h2 (norm_nonneg _) (by positivity)
    _ = _ := by ring

/-- In the right region: |F z| <= C/(1+|Im z|)^2. -/
theorem norm_Fq_right_le (hH : HBound H M β D N) {z : ℂ} (hz : 2 ≤ z.re) :
    ‖Fq H z‖ ≤ 3 * M * (D / min D 1) ^ 2 * ((1 + |z.im|) ^ 2)⁻¹ := by
  have h := norm_Fq_right hH hz
  have hw := weight_mono hH.D_pos N (z := z)
  have hwl := weight_le_inv_sq hH.D_pos (by linarith [hH.N_ge] : 2 ≤ N) z.im
  have hM := hH.M_nonneg
  calc ‖Fq H z‖ ≤ 3 * M * (D / (D + ‖z - 1 / 2‖)) ^ N := h
    _ ≤ 3 * M * ((D / min D 1) ^ 2 * ((1 + |z.im|) ^ 2)⁻¹) := by gcongr; exact hw.trans hwl
    _ = _ := by ring

theorem norm_Fq_line_right (hH : HBound H M β D N) {σ : ℝ} (hσ : 2 ≤ σ) (y : ℝ) :
    ‖Fq H (σ + y * I)‖ ≤ 3 * M * (D / min D 1) ^ 2 * ((1 + |y|) ^ 2)⁻¹ := by
  have hz : 2 ≤ ((σ : ℂ) + y * I).re := by simp; exact hσ
  have h := norm_Fq_right_le hH hz
  have him : ((σ : ℂ) + y * I).im = y := by simp
  rwa [him] at h

theorem continuous_Fq_line (hH : HBound H M β D N) {σ : ℝ} (hσ : 2 ≤ σ) :
    Continuous fun y : ℝ => Fq H (σ + y * I) := by
  rw [continuous_iff_continuousAt]
  intro y
  have hz : 2 ≤ ((σ : ℂ) + y * I).re := by simp; exact hσ
  exact ContinuousAt.comp (f := fun y : ℝ => (σ : ℂ) + y * I) (x := y)
    (continuousAt_Fq hH (ne_one_of_two_le hz) (riemannZeta_ne_zero_of_one_lt_re (by linarith)))
    (by fun_prop : Continuous fun y : ℝ => (σ : ℂ) + y * I).continuousAt

theorem integrable_Fq_line (hH : HBound H M β D N) {σ : ℝ} (hσ : 2 ≤ σ) :
    Integrable fun y : ℝ => Fq H (σ + y * I) :=
  integrable_of_le_inv_sq (continuous_Fq_line hH hσ).aestronglyMeasurable (norm_Fq_line_right hH hσ)

theorem norm_kap_le_one {L : ℝ} (hL : 0 ≤ L) {x : ℂ} (hx : 0 ≤ x.re) : ‖kap L x‖ ≤ 1 := by
  rw [norm_kap, Real.exp_le_one_iff]; nlinarith

/-- **Fubini**: G(s) = (2 pi)^{-1} (int Psi0 - int Psi) on Re z = 3, for Re s > 3. -/
theorem GF_eq (hH : HBound H M β D N) {L : ℝ} (hL : 0 ≤ L) {s : ℂ} (hs : 3 < s.re) :
    GF H L s = (2 * (π : ℂ))⁻¹ * ((∫ y : ℝ, Psi0 H s (3 + y * I)) - ∫ y : ℝ, Psi H L s (3 + y * I)) := by
  have hFi : Integrable (fun y : ℝ => Fq H (3 + y * I)) := by
    have := integrable_Fq_line hH (σ := 3) (by norm_num); simpa using this
  have hFc : Continuous (fun y : ℝ => Fq H (3 + y * I)) := by
    have := continuous_Fq_line hH (σ := 3) (by norm_num); simpa using this
  have hsz : ∀ y : ℝ, s - (3 + y * I) ≠ 0 := fun y h => by
    have := congrArg Complex.re h; simp at this; linarith
  have hsz' : ∀ y : ℝ, 3 - s.re ≤ -(s - (3 + y * I)).re := fun y => by simp
  -- the product integrand
  set φ : ℝ → ℝ → ℂ := fun u y => Fq H (3 + y * I) * cexp (u * (3 + y * I)) * cexp (-(u * s)) with hφ
  have hφc : Continuous (Function.uncurry φ) := by
    simp only [hφ, Function.uncurry]
    exact ((hFc.comp continuous_snd).mul (by fun_prop)).mul (by fun_prop)
  haveI : IsFiniteMeasure (volume.restrict (Set.uIoc (0 : ℝ) L)) := by
    rw [Set.uIoc_of_le hL]; infer_instance
  have hφi : Integrable (Function.uncurry φ) ((volume.restrict (Set.uIoc 0 L)).prod volume) := by
    have hdom : Integrable (fun p : ℝ × ℝ => (fun _ : ℝ => Real.exp (3 * L)) p.1 * (fun y => ‖Fq H (3 + y * I)‖) p.2)
        ((volume.restrict (Set.uIoc 0 L)).prod volume) :=
      Integrable.mul_prod (integrable_const _) hFi.norm
    refine hdom.mono' hφc.aestronglyMeasurable ?_
    rw [Measure.restrict_prod_eq_prod_univ]
    filter_upwards [ae_restrict_mem (measurableSet_uIoc.prod MeasurableSet.univ)] with p hp
    obtain ⟨hp1, -⟩ := hp
    rw [Set.uIoc_of_le hL] at hp1
    simp only [hφ, Function.uncurry, norm_mul, Complex.norm_exp]
    have e1 : ((p.1 : ℂ) * (3 + p.2 * I)).re = 3 * p.1 := by simp; ring
    have e2 : (-((p.1 : ℂ) * s)).re = -(p.1 * s.re) := by simp
    rw [e1, e2]
    have h3 : Real.exp (3 * p.1) * Real.exp (-(p.1 * s.re)) ≤ Real.exp (3 * L) := by
      rw [← Real.exp_add]; apply Real.exp_le_exp.mpr; nlinarith [hp1.1, hp1.2]
    calc ‖Fq H (3 + p.2 * I)‖ * Real.exp (3 * p.1) * Real.exp (-(p.1 * s.re))
        = ‖Fq H (3 + p.2 * I)‖ * (Real.exp (3 * p.1) * Real.exp (-(p.1 * s.re))) := by ring
      _ ≤ ‖Fq H (3 + p.2 * I)‖ * Real.exp (3 * L) := mul_le_mul_of_nonneg_left h3 (norm_nonneg _)
      _ = _ := by ring
  have hP0i : Integrable (fun y : ℝ => Psi0 H s (3 + y * I)) := by
    refine (hFi.norm.const_mul (1 / (s.re - 3))).mono' ?_ (Filter.Eventually.of_forall fun y => ?_)
    · exact (hFc.div (by fun_prop) hsz).aestronglyMeasurable
    · rw [Psi0, norm_div, div_eq_mul_inv, mul_comm (1 / (s.re - 3))]
      refine mul_le_mul_of_nonneg_left ?_ (norm_nonneg _)
      rw [one_div]
      apply inv_anti₀ (by linarith)
      have := Complex.abs_re_le_norm (s - (3 + y * I))
      simp at this
      have h' : s.re - 3 ≤ |s.re - 3| := le_abs_self _
      linarith
  have hPi : Integrable (fun y : ℝ => Psi H L s (3 + y * I)) := by
    refine hP0i.norm.mono' ?_ (Filter.Eventually.of_forall fun y => ?_)
    · exact ((hFc.mul (by unfold kap; fun_prop)).div (by fun_prop) hsz).aestronglyMeasurable
    · rw [Psi, Psi0, norm_div, norm_div, norm_mul]
      have hk : ‖kap L (s - (3 + y * I))‖ ≤ 1 := norm_kap_le_one hL (by simp; linarith)
      rw [div_le_div_iff_of_pos_right (norm_pos_iff.mpr (hsz y))]
      exact mul_le_of_le_one_right (norm_nonneg _) hk
  unfold GF gF
  have e1 : ∀ u : ℝ, ((2 * (π : ℂ))⁻¹ * ∫ y : ℝ, Fq H (3 + y * I) * cexp (u * (3 + y * I))) * cexp (-(u * s)) =
      (2 * (π : ℂ))⁻¹ * ∫ y : ℝ, φ u y := by
    intro u
    simp only [hφ]
    rw [mul_assoc, ← MeasureTheory.integral_mul_const]
  simp_rw [e1]
  rw [intervalIntegral.integral_const_mul, intervalIntegral_integral_swap hφi, ← integral_sub hP0i hPi]
  congr 1
  refine integral_congr_ae (Filter.Eventually.of_forall fun y => ?_)
  simp only [hφ]
  have hc : (3 + (y : ℂ) * I - s) ≠ 0 := fun h => hsz y (by rw [← neg_sub, h, neg_zero])
  have e2 : ∀ u : ℝ, Fq H (3 + y * I) * cexp (u * (3 + y * I)) * cexp (-(u * s)) =
      Fq H (3 + y * I) * cexp ((3 + y * I - s) * u) := by
    intro u; rw [mul_assoc, ← Complex.exp_add]; congr 2; ring
  simp_rw [e2]
  rw [intervalIntegral.integral_const_mul, integral_exp_mul_complex hc]
  simp only [Psi, Psi0, kap, Complex.ofReal_zero, mul_zero, Complex.exp_zero]
  have hsz0 := hsz y
  field_simp
  ring_nf

end Riemann.Deep.LeanNBRate

