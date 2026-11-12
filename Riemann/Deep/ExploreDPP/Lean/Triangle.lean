import Riemann.Deep.ExploreDPP.Lean.TheoremA

/-!
# The periodized triangle: Fourier series (makes Theorem A unconditional)

𝓕Δ = sinc² (direct integration, FTC), the Fourier transform of a real even function is real, and
Mathlib's Poisson summation formula (with shift) give TriangleSeries a for every a > 0:
a + 2a Σ_{k≥1} sinc²(ka) cos(2π α k a) = Σ_m Δ(α - m/a). Hence theoremA' (no hypothesis).
-/

open Real MeasureTheory FourierTransform Filter Asymptotics
open scoped Real

noncomputable section

namespace ExploreDPP

lemma tri_neg (x : ℝ) : tri (-x) = tri x := by simp [tri, abs_neg]

lemma continuous_tri : Continuous tri := by unfold tri; fun_prop

lemma tri_eq_zero {x : ℝ} (hx : 1 ≤ |x|) : tri x = 0 := by
  unfold tri; exact max_eq_left (by linarith)

lemma hasCompactSupport_tri : HasCompactSupport tri := by
  apply HasCompactSupport.intro (isCompact_closedBall (0:ℝ) 1)
  intro x hx
  simp only [Metric.mem_closedBall, dist_zero_right, Real.norm_eq_abs, not_le] at hx
  exact tri_eq_zero hx.le

/-- The Fourier transform of a real even integrable function is real. -/
lemma fourier_real_of_even {f : ℝ → ℝ}
    (hev : ∀ x, f (-x) = f x) (ξ : ℝ) : 𝓕 (fun x => (f x : ℂ)) ξ = (FT f ξ : ℂ) := by
  have key : ∀ v : ℝ, (starRingEnd ℂ) (Complex.exp (↑(-2 * π * v * ξ) * Complex.I) • (f v : ℂ)) =
      Complex.exp (↑(-2 * π * (-v) * ξ) * Complex.I) • (f (-v) : ℂ) := by
    intro v
    rw [hev v, smul_eq_mul, smul_eq_mul, map_mul, Complex.conj_ofReal, ← Complex.exp_conj, map_mul,
      Complex.conj_ofReal, Complex.conj_I]
    congr 2; push_cast; ring
  have hconj : (starRingEnd ℂ) (𝓕 (fun x => (f x : ℂ)) ξ) = 𝓕 (fun x => (f x : ℂ)) ξ := by
    calc (starRingEnd ℂ) (𝓕 (fun x => (f x : ℂ)) ξ)
        = ∫ v, (starRingEnd ℂ) (Complex.exp (↑(-2 * π * v * ξ) * Complex.I) • (f v : ℂ)) := by
          rw [Real.fourier_real_eq_integral_exp_smul, integral_conj]
      _ = ∫ v, Complex.exp (↑(-2 * π * (-v) * ξ) * Complex.I) • (f (-v) : ℂ) := by
          congr 1; funext v; exact key v
      _ = ∫ v, Complex.exp (↑(-2 * π * v * ξ) * Complex.I) • (f v : ℂ) :=
          integral_neg_eq_self (fun v => Complex.exp (↑(-2 * π * v * ξ) * Complex.I) • (f v : ℂ)) volume
      _ = 𝓕 (fun x => (f x : ℂ)) ξ := (Real.fourier_real_eq_integral_exp_smul _ _).symm
  exact (Complex.conj_eq_iff_re.mp hconj).symm

lemma FT_tri (ξ : ℝ) : FT tri ξ = Real.sinc (π * ξ) ^ 2 := by
  rw [FT_eq_integral_cos continuous_tri hasCompactSupport_tri]
  have hcont : Continuous (fun x => tri x * Real.cos (2 * π * x * ξ)) := by
    have := continuous_tri; fun_prop
  have h1 : ∫ x, tri x * Real.cos (2 * π * x * ξ) =
      ∫ x in (-1:ℝ)..1, tri x * Real.cos (2 * π * x * ξ) := by
    rw [intervalIntegral.integral_of_le (by norm_num)]
    symm
    apply setIntegral_eq_integral_of_forall_compl_eq_zero
    intro x hx
    simp only [Set.mem_Ioc, not_and_or, not_lt, not_le] at hx
    rw [tri_eq_zero, zero_mul]
    rcases hx with hx | hx
    · rw [abs_of_nonpos (by linarith)]; linarith
    · rw [abs_of_pos (by linarith)]; linarith
  have hadd := intervalIntegral.integral_add_adjacent_intervals (μ := volume)
    (hcont.intervalIntegrable (-1) 0) (hcont.intervalIntegrable 0 1)
  have hsym : ∫ x in (-1:ℝ)..0, tri x * Real.cos (2 * π * x * ξ) =
      ∫ x in (0:ℝ)..1, tri x * Real.cos (2 * π * x * ξ) := by
    have := intervalIntegral.integral_comp_neg (a := 0) (b := 1)
      (fun x => tri x * Real.cos (2 * π * x * ξ))
    simp only [neg_zero, tri_neg, mul_neg, neg_mul, Real.cos_neg] at this
    exact this.symm
  have h01 : ∫ x in (0:ℝ)..1, tri x * Real.cos (2 * π * x * ξ) =
      ∫ x in (0:ℝ)..1, (1 - x) * Real.cos (2 * π * x * ξ) := by
    apply intervalIntegral.integral_congr
    intro x hx
    rw [Set.uIcc_of_le zero_le_one] at hx
    simp only [tri, abs_of_nonneg hx.1]
    rw [max_eq_right (by linarith [hx.2])]
  rw [h1, ← hadd, hsym, h01]
  rcases eq_or_ne ξ 0 with hξ | hξ
  · subst hξ
    simp only [mul_zero, Real.cos_zero, mul_one, Real.sinc_zero, one_pow]
    have h12 : ∫ x in (0:ℝ)..1, (1 - x) = 1 / 2 := by
      rw [intervalIntegral.integral_sub (f := fun _ => (1:ℝ)) (g := fun x => x)
        (continuous_const.intervalIntegrable 0 1) (continuous_id.intervalIntegrable 0 1)]
      simp
      norm_num
    rw [h12]; norm_num
  · have hπ := Real.pi_pos.ne'
    have hd : ∀ x : ℝ, HasDerivAt (fun x => (1 - x) * Real.sin (2 * π * x * ξ) / (2 * π * ξ) -
        Real.cos (2 * π * x * ξ) / (2 * π * ξ) ^ 2) ((1 - x) * Real.cos (2 * π * x * ξ)) x := by
      intro x
      have hlin : HasDerivAt (fun x : ℝ => 2 * π * x * ξ) (2 * π * ξ) x := by
        have := ((hasDerivAt_id x).const_mul (2 * π)).mul_const ξ
        simpa using this
      have hs := (Real.hasDerivAt_sin (2 * π * x * ξ)).comp x hlin
      have hc := (Real.hasDerivAt_cos (2 * π * x * ξ)).comp x hlin
      have h1' : HasDerivAt (fun x : ℝ => 1 - x) (-1) x := by
        simpa using (hasDerivAt_id x).const_sub 1
      have := ((h1'.mul hs).div_const (2 * π * ξ)).sub (hc.div_const ((2 * π * ξ) ^ 2))
      convert this using 1
      · ext y; simp [Pi.mul_apply, Function.comp_apply]
      · simp only [Function.comp_apply]; field_simp; ring
    rw [intervalIntegral.integral_eq_sub_of_hasDerivAt (fun x _ => hd x)
      ((by fun_prop : Continuous (fun x : ℝ => (1 - x) * Real.cos (2 * π * x * ξ))).intervalIntegrable 0 1)]
    have hc2 : Real.cos (2 * π * ξ) = 1 - 2 * Real.sin (π * ξ) ^ 2 := by
      rw [show 2 * π * ξ = 2 * (π * ξ) by ring, Real.cos_two_mul]
      linarith [Real.sin_sq_add_cos_sq (π * ξ)]
    rw [Real.sinc_of_ne_zero (mul_ne_zero hπ hξ)]
    simp only [sub_self, zero_mul, zero_div, mul_one, mul_zero, Real.sin_zero, Real.cos_zero,
      sub_zero, zero_sub]
    rw [hc2]
    field_simp
    ring


lemma fourier_comp_div (g : ℝ → ℂ) {h : ℝ} (hh : 0 < h) (ξ : ℝ) :
    𝓕 (fun x => g (x / h)) ξ = (h : ℂ) * 𝓕 g (h * ξ) := by
  rw [Real.fourier_real_eq, Real.fourier_real_eq]
  have := MeasureTheory.Measure.integral_comp_div (fun y : ℝ => 𝐞 (-((h * y) * ξ)) • g y) h
  rw [abs_of_pos hh] at this
  have e1 : (fun v : ℝ => 𝐞 (-(v * ξ)) • g (v / h)) =
      fun v => 𝐞 (-((h * (v / h)) * ξ)) • g (v / h) := by
    funext v; congr 3; field_simp
  rw [e1, this, Complex.real_smul]
  congr 1
  congr 1; funext y; congr 3; ring

/-- The Fourier series of the periodized triangle, from Mathlib's Poisson summation and
𝓕Δ = sinc²: TriangleSeries holds for every a > 0. -/
theorem triangleSeries {a : ℝ} (ha : 0 < a) : TriangleSeries a := by
  intro α
  have hTF : ∀ ξ : ℝ, 𝓕 (fun x => (tri x : ℂ)) ξ = ((Real.sinc (π * ξ) ^ 2 : ℝ) : ℂ) := fun ξ => by
    rw [fourier_real_of_even tri_neg, FT_tri]
  set G : ℝ → ℂ := fun t => ((tri (t / a) : ℝ) : ℂ) with hG
  have hFG : ∀ ξ : ℝ, 𝓕 G ξ = ((a * Real.sinc (π * (a * ξ)) ^ 2 : ℝ) : ℂ) := by
    intro ξ
    rw [hG, fourier_comp_div (fun x => (tri x : ℂ)) ha ξ, hTF]
    push_cast; ring
  have hGc : Continuous G :=
    (Complex.continuous_ofReal.comp continuous_tri).comp (continuous_id.div_const a)
  have hGO : G =O[cocompact ℝ] (fun x : ℝ => |x| ^ (-(2:ℝ))) := by
    have hev : G =ᶠ[cocompact ℝ] 0 := by
      have hmem := (isCompact_closedBall (0:ℝ) a).compl_mem_cocompact
      filter_upwards [hmem] with x hx
      simp only [Set.mem_compl_iff, Metric.mem_closedBall, dist_zero_right, Real.norm_eq_abs,
        not_le] at hx
      simp only [hG, Pi.zero_apply, Complex.ofReal_eq_zero]
      apply tri_eq_zero
      rw [abs_div, abs_of_pos ha, le_div_iff₀ ha]; linarith
    exact hev.isBigO.trans (isBigO_zero _ _)
  have hFGO : (𝓕 G) =O[cocompact ℝ] (fun x : ℝ => |x| ^ (-(2:ℝ))) := by
    refine IsBigO.of_bound (1 / (π ^ 2 * a)) ?_
    have hmem : ({0}ᶜ : Set ℝ) ∈ cocompact ℝ := isCompact_singleton.compl_mem_cocompact
    filter_upwards [hmem] with x hx
    have hx0 : x ≠ 0 := hx
    have hy : π * (a * x) ≠ 0 := mul_ne_zero Real.pi_pos.ne' (mul_ne_zero ha.ne' hx0)
    rw [hFG x, Complex.norm_real, Real.norm_eq_abs, Real.norm_eq_abs,
      abs_of_nonneg (by positivity : (0:ℝ) ≤ a * Real.sinc (π * (a * x)) ^ 2),
      abs_of_nonneg (by positivity : (0:ℝ) ≤ |x| ^ (-(2:ℝ))),
      Real.rpow_neg (abs_nonneg x), Real.rpow_two, sq_abs]
    calc a * Real.sinc (π * (a * x)) ^ 2 ≤ a * (1 / (π * (a * x)) ^ 2) := by
          gcongr; exact sinc_sq_le hy
      _ = 1 / (π ^ 2 * a) * (x ^ 2)⁻¹ := by
          have := Real.pi_pos.ne'
          field_simp
  have hsFGn : Summable (fun n : ℤ => ‖𝓕 G n‖) :=
    summable_of_isBigO (Real.summable_abs_int_rpow one_lt_two)
      ((hFGO.norm_left).comp_tendsto Int.tendsto_coe_cofinite)
  have hpois := Real.tsum_eq_tsum_fourier_of_rpow_decay hGc one_lt_two hGO hFGO (-(a * α))
  have hL : ∀ n : ℤ, G (-(a * α) + n) = ((tri (α - n / a) : ℝ) : ℂ) := by
    intro n
    simp only [hG]
    congr 1
    rw [← tri_neg]
    congr 1
    field_simp
    ring
  have hR : ∀ n : ℤ, 𝓕 G n * fourier n ((-(a * α) : ℝ) : UnitAddCircle) =
      ((a * Real.sinc (π * (a * n)) ^ 2 : ℝ) : ℂ) *
        Complex.exp (((-(2 * π * n * (a * α))) : ℝ) * Complex.I) := by
    intro n
    rw [hFG, fourier_coe_apply]
    congr 2
    push_cast
    ring
  simp only [hL, hR] at hpois
  have hsumR : Summable (fun n : ℤ => ((a * Real.sinc (π * (a * n)) ^ 2 : ℝ) : ℂ) *
      Complex.exp (((-(2 * π * n * (a * α))) : ℝ) * Complex.I)) := by
    refine Summable.of_norm_bounded hsFGn (fun n => ?_)
    rw [norm_mul, Complex.norm_exp_ofReal_mul_I, mul_one, hFG]
  have e3 := congrArg Complex.re hpois
  rw [← Complex.ofReal_tsum, Complex.ofReal_re, Complex.re_tsum hsumR] at e3
  simp only [Complex.re_ofReal_mul, Complex.exp_ofReal_mul_I_re, Real.cos_neg] at e3
  set c : ℤ → ℝ := fun k => a * Real.sinc (π * (a * k)) ^ 2 * Real.cos (2 * π * k * (a * α))
    with hc
  have hcev : ∀ k : ℤ, c (-k) = c k := by
    intro k
    simp only [hc, Int.cast_neg, mul_neg, Real.sinc_neg, neg_mul, Real.cos_neg]
  have hcs : Summable c := by
    have := Complex.reCLM.summable hsumR
    refine this.congr (fun n => ?_)
    simp only [Complex.reCLM_apply, Complex.re_ofReal_mul, Complex.exp_ofReal_mul_I_re,
      Real.cos_neg, hc]
  have hfold := tsum_int_fold hcev hcs
  have h0 : c 0 = a := by simp [hc]
  have hn : ∀ n : ℕ, c ((n : ℤ) + 1) = a * (Real.sinc (π * ((n + 1) * a)) ^ 2 *
      Real.cos (2 * π * α * ((n + 1) * a))) := by
    intro n; simp only [hc]; push_cast; ring_nf
  rw [e3, show tsum c = ∑' n : ℤ, c n from rfl, hfold, h0, tsum_congr hn, tsum_mul_left]
  ring

/-- Theorem A, unconditional: for 0 < θ ≤ 1 and λ ≤ 1/(1+θ) every admissible f has
J_θ(f) ≤ 1 = f^(0). -/
theorem theoremA' {θ lam : ℝ} (hθ0 : 0 < θ) (hθ1 : θ ≤ 1) (hlam : lam ≤ 1 / (1 + θ))
    {f : ℝ → ℝ} (hf : Admissible lam f) : J θ f ≤ 1 :=
  theoremA hθ0 hθ1 (triangleSeries (by positivity)) hlam hf


/-- Corollary in Lamzouri's range (1/2 < θ < 4/7): no λ ≤ 7/11 meets the hypotheses of Thm Main2. -/
theorem not_main2Hyp_of_le_7_11 {θ lam : ℝ} (hlam : lam ≤ 7 / 11) (f : ℝ → ℝ) :
    ¬ Main2Hyp θ lam f := by
  rintro ⟨hθ1, hθ2, -, -, hf, hJ⟩
  have h711 : (7:ℝ) / 11 ≤ 1 / (1 + θ) := by
    rw [div_le_div_iff₀ (by norm_num) (by linarith)]; linarith
  have := theoremA' (by linarith) (by linarith) (hlam.trans h711) hf
  linarith

/-- Corollary at θ = 1 (the RH / Montgomery range): λ ≤ 1/2 forces J_1(f) ≤ 1. -/
theorem theoremA_RH {lam : ℝ} (hlam : lam ≤ 1 / 2) {f : ℝ → ℝ} (hf : Admissible lam f) :
    J 1 f ≤ 1 :=
  theoremA' one_pos le_rfl (by rw [show (1:ℝ) / (1 + 1) = 1 / 2 by norm_num]; exact hlam) hf

end ExploreDPP

