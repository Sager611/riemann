import Riemann.Deep.Hardy.Repr

/-!
# Hardy's theorem, part 3: contour shift and exponential decay of Λ on the critical line

* integral_shift: for f holomorphic on a closed horizontal strip 0 ≤ Im ≤ c (or c ≤ Im ≤ 0)
  with ‖f z‖ ≤ C e^{-a|Re z|}, ∫ f(x) dx = ∫ f(x + c i) dx (rectangle Cauchy + limits).
* Lam_mul_exp: Λ(1/2+it) e^{αt} = (1/2) 𝓕 (Hk(· - 2αi)) (t/4π) for |α| < π/4.
* Lam_decay: ‖Λ(1/2+it)‖ ≤ C e^{-γ|t|} for 0 ≤ γ < π/4.
-/

open Complex Real MeasureTheory Set Filter Topology
open scoped ComplexConjugate FourierTransform Interval

namespace Fayed.Hardy

lemma integrable_exp_neg_abs {a : ℝ} (ha : 0 < a) :
    Integrable (fun x : ℝ => rexp (-a * |x|)) := by
  rw [← integrableOn_univ, ← Iic_union_Ioi (a := (0 : ℝ)), integrableOn_union]
  constructor
  · refine (integrableOn_exp_mul_Iic ha 0).congr_fun (fun x hx => ?_) measurableSet_Iic
    rw [abs_of_nonpos (mem_Iic.mp hx)]; ring_nf
  · refine (exp_neg_integrableOn_Ioi 0 ha).congr_fun (fun x hx => ?_) measurableSet_Ioi
    rw [abs_of_pos (mem_Ioi.mp hx)]

/-- Horizontal contour shift. -/
theorem integral_shift {f : ℂ → ℂ} {c a C : ℝ} (ha : 0 < a)
    (hd : DifferentiableOn ℂ f {z | z.im ∈ [[0, c]]})
    (hb : ∀ z : ℂ, z.im ∈ [[0, c]] → ‖f z‖ ≤ C * rexp (-a * |z.re|)) :
    ∫ x : ℝ, f x = ∫ x : ℝ, f (x + c * I) := by
  have hcont : ∀ y ∈ [[0, c]], Continuous (fun x : ℝ => f (x + y * I)) := by
    intro y hy
    refine hd.continuousOn.comp_continuous (by fun_prop) (fun x => ?_)
    show (↑x + ↑y * I).im ∈ [[0, c]]
    simpa using hy
  have hint : ∀ y ∈ [[0, c]], Integrable (fun x : ℝ => f (x + y * I)) := by
    intro y hy
    refine Integrable.mono' ((integrable_exp_neg_abs ha).const_mul C)
      (hcont y hy).aestronglyMeasurable (Eventually.of_forall fun x => ?_)
    have := hb (x + y * I) (by simpa using hy)
    simpa using this
  have h0 : (0 : ℝ) ∈ [[0, c]] := left_mem_uIcc
  have hc : c ∈ [[0, c]] := right_mem_uIcc
  have hint0 : Integrable (fun x : ℝ => f x) := by simpa using hint 0 h0
  -- rectangle identity
  have hrect : ∀ R : ℝ, ((∫ x in -R..R, f x) - ∫ x in -R..R, f (x + c * I)) +
      I • (∫ y in (0 : ℝ)..c, f (R + y * I)) - I • (∫ y in (0 : ℝ)..c, f (-R + y * I)) = 0 := by
    intro R
    have := integral_boundary_rect_eq_zero_of_differentiableOn f (-R : ℂ) (R + c * I)
      (hd.mono (fun z hz => by
        show z.im ∈ [[0, c]]
        have := (mem_reProdIm.mp hz).2
        simpa using this))
    simpa using this
  -- limits
  have hA : Tendsto (fun R : ℝ => ∫ x in -R..R, f x) atTop (𝓝 (∫ x : ℝ, f x)) :=
    intervalIntegral_tendsto_integral hint0 tendsto_neg_atTop_atBot tendsto_id
  have hB : Tendsto (fun R : ℝ => ∫ x in -R..R, f (x + c * I)) atTop
      (𝓝 (∫ x : ℝ, f (x + c * I))) :=
    intervalIntegral_tendsto_integral (hint c hc) tendsto_neg_atTop_atBot tendsto_id
  have hV : ∀ s : ℝ, s = 1 ∨ s = -1 →
      Tendsto (fun R : ℝ => ∫ y in (0 : ℝ)..c, f (s * R + y * I)) atTop (𝓝 0) := by
    intro s hs
    have hs1 : |s| = 1 := by rcases hs with rfl | rfl <;> simp
    have hbound : ∀ R : ℝ, ‖∫ y in (0 : ℝ)..c, f (s * R + y * I)‖ ≤
        C * rexp (-a * |R|) * |c - 0| := by
      intro R
      refine intervalIntegral.norm_integral_le_of_norm_le_const (fun y hy => ?_)
      have hy' : y ∈ [[0, c]] := uIoc_subset_uIcc hy
      have := hb (s * R + y * I) (by simpa using hy')
      simpa [abs_mul, hs1] using this
    have hlim : Tendsto (fun R : ℝ => C * rexp (-a * |R|) * |c - 0|) atTop (𝓝 0) := by
      have h1 : Tendsto (fun R : ℝ => -a * |R|) atTop atBot := by
        refine (tendsto_abs_atTop_atTop.comp tendsto_id).const_mul_atTop_of_neg ?_ |>.congr
          (fun R => rfl)
        linarith
      have h2 := (Real.tendsto_exp_atBot.comp h1).const_mul C |>.mul_const |c - 0|
      simpa using h2
    exact squeeze_zero_norm hbound hlim
  have hV1 := (hV 1 (Or.inl rfl)).const_smul I
  have hV2 := (hV (-1) (Or.inr rfl)).const_smul I
  simp only [smul_zero] at hV1 hV2
  have hlim := ((hA.sub hB).add hV1).sub hV2
  rw [add_zero, sub_zero] at hlim
  have hzero : Tendsto (fun R : ℝ => ((∫ x in -R..R, f x) - ∫ x in -R..R, f (x + c * I)) +
      I • (∫ y in (0 : ℝ)..c, f (R + y * I)) - I • (∫ y in (0 : ℝ)..c, f (-R + y * I)))
      atTop (𝓝 0) := by simp only [hrect]; exact tendsto_const_nhds
  exact sub_eq_zero.mp (tendsto_nhds_unique hlim (hzero.congr (fun R => by simp)))

lemma mem_strip_of {z : ℂ} {b : ℝ} (hb : b < π / 2) (h : |z.im| ≤ b) : z ∈ strip :=
  show |z.im| < π / 2 from lt_of_le_of_lt h hb

lemma uIcc_abs_le {y c : ℝ} (hy : y ∈ [[0, c]]) : |y| ≤ |c| := by
  rcases le_total 0 c with h | h
  · rw [uIcc_of_le h] at hy
    rw [abs_of_nonneg hy.1, abs_of_nonneg h]; exact hy.2
  · rw [uIcc_of_ge h] at hy
    rw [abs_of_nonpos hy.2, abs_of_nonpos h]; linarith [hy.1]

/-- Shifted kernel lines are continuous. -/
lemma continuous_Hk_line {c : ℝ} (hc : |c| < π / 2) :
    Continuous (fun x : ℝ => Hk (x + c * I)) := by
  refine differentiableOn_Hk.continuousOn.comp_continuous (by fun_prop) (fun x => ?_)
  show |(↑x + ↑c * I).im| < π / 2
  simpa using hc

/-- Shifted kernel lines are integrable. -/
lemma integrable_Hk_line {c : ℝ} (hc : |c| < π / 2) :
    Integrable (fun x : ℝ => Hk (x + c * I)) := by
  obtain ⟨C, hC⟩ := Hk_bound (abs_nonneg c) hc
  refine Integrable.mono' ((integrable_exp_neg_abs (by norm_num : (0 : ℝ) < 1 / 4)).const_mul C)
    (continuous_Hk_line hc).aestronglyMeasurable (Eventually.of_forall fun x => ?_)
  have := hC (x + c * I) (by simp)
  simp only [add_re, ofReal_re, mul_re, ofReal_im, I_re, mul_zero, I_im, mul_one, sub_self,
    add_zero] at this
  rw [show -|x| / 4 = -(1 / 4) * |x| by ring] at this
  exact this

lemma norm_fourier_le {h : ℝ → ℂ} (w : ℝ) : ‖𝓕 h w‖ ≤ ∫ v, ‖h v‖ := by
  rw [Real.fourier_real_eq]
  refine (norm_integral_le_integral_norm _).trans (le_of_eq ?_)
  congr 1; ext v; rw [Circle.norm_smul]

/-- Shift of the Fourier representation. -/
theorem Lam_mul_exp {α : ℝ} (hα : |α| < π / 4) (t : ℝ) :
    completedRiemannZeta (1 / 2 + t * I) * rexp (α * t) =
      (1 / 2) * 𝓕 (fun v : ℝ => Hk (v + (-2 * α : ℝ) * I)) (t / (4 * π)) := by
  set w := t / (4 * π) with hw
  set c : ℝ := -2 * α with hc_def
  have hc : |c| < π / 2 := by
    rw [hc_def, abs_mul]; norm_num; linarith
  set f : ℂ → ℂ := fun z => cexp (-(2 * π * w : ℝ) * z * I) * Hk z with hf
  obtain ⟨C, hC⟩ := Hk_bound (abs_nonneg c) hc
  have hshift : ∫ x : ℝ, f x = ∫ x : ℝ, f (x + c * I) := by
    refine integral_shift (a := 1 / 4) (C := rexp (2 * π * |w| * |c|) * C) (by norm_num) ?_ ?_
    · intro z hz
      have hz' : z ∈ strip := mem_strip_of hc (uIcc_abs_le hz)
      have hd1 : Differentiable ℂ fun z : ℂ => cexp (-(2 * π * w : ℝ) * z * I) := by fun_prop
      exact ((hd1 z).mul ((differentiableAt_Rk hz').sub (differentiable_Ck _))).differentiableWithinAt
    · intro z hz
      have hz1 := uIcc_abs_le hz
      rw [hf]; simp only
      rw [norm_mul, Complex.norm_exp]
      have him : (-(2 * π * w : ℝ) * z * I : ℂ).re = 2 * π * w * z.im := by
        simp [mul_re, mul_im]
      rw [him]
      have hk := hC z hz1
      have he : rexp (2 * π * w * z.im) ≤ rexp (2 * π * |w| * |c|) := by
        apply Real.exp_le_exp.mpr
        have h1 : w * z.im ≤ |w| * |c| := (le_abs_self _).trans (by
          rw [abs_mul]; exact mul_le_mul_of_nonneg_left hz1 (abs_nonneg _))
        nlinarith [Real.pi_pos]
      calc rexp (2 * π * w * z.im) * ‖Hk z‖ ≤ rexp (2 * π * |w| * |c|) * (C * rexp (-|z.re| / 4)) :=
            mul_le_mul he hk (norm_nonneg _) (Real.exp_pos _).le
        _ = rexp (2 * π * |w| * |c|) * C * rexp (-(1 / 4) * |z.re|) := by ring_nf
  rw [Lam_eq_fourier, Real.fourier_real_eq_integral_exp_smul, Real.fourier_real_eq_integral_exp_smul]
  have hL : ∫ v : ℝ, cexp (↑(-2 * π * v * w) * I) • Hk ↑v = ∫ x : ℝ, f x := by
    congr 1; ext v; rw [hf, smul_eq_mul]; congr 2; push_cast; ring
  have hR : ∫ x : ℝ, f (x + c * I) =
      (rexp (-(α * t)) : ℂ) * ∫ v : ℝ, cexp (↑(-2 * π * v * w) * I) • Hk (↑v + ↑c * I) := by
    rw [← integral_const_mul]
    congr 1; ext v
    rw [hf, smul_eq_mul]; simp only
    rw [← mul_assoc]; congr 1
    rw [ofReal_exp, ← Complex.exp_add]; congr 1
    have hreal : 2 * π * w * c = -(α * t) := by
      rw [hw, hc_def]; field_simp; ring
    rw [← hreal]; push_cast
    linear_combination (-(2 * π * w * c : ℂ)) * I_sq
  rw [hL, hshift, hR]
  have hexp : (rexp (α * t) : ℂ) * (rexp (-(α * t)) : ℂ) = 1 := by
    rw [← ofReal_mul, ← Real.exp_add, add_neg_cancel, Real.exp_zero, ofReal_one]
  calc (1 / 2 : ℂ) * ((rexp (-(α * t)) : ℂ) *
        ∫ v : ℝ, cexp (↑(-2 * π * v * w) * I) • Hk (↑v + ↑c * I)) * (rexp (α * t) : ℂ)
      = (1 / 2) * ((rexp (α * t) : ℂ) * (rexp (-(α * t)) : ℂ)) *
        ∫ v : ℝ, cexp (↑(-2 * π * v * w) * I) • Hk (↑v + ↑c * I) := by ring
    _ = _ := by rw [hexp, mul_one]

/-- Exponential decay of Λ on the critical line. -/
theorem Lam_decay {γ : ℝ} (hγ0 : 0 ≤ γ) (hγ : γ < π / 4) :
    ∃ C, ∀ t : ℝ, ‖completedRiemannZeta (1 / 2 + t * I)‖ ≤ C * rexp (-γ * |t|) := by
  have h1 : |γ| < π / 4 := by rwa [abs_of_nonneg hγ0]
  have h2 : |-γ| < π / 4 := by rwa [abs_neg, abs_of_nonneg hγ0]
  set N₁ := ∫ v : ℝ, ‖Hk (v + (-2 * γ : ℝ) * I)‖
  set N₂ := ∫ v : ℝ, ‖Hk (v + (-2 * -γ : ℝ) * I)‖
  refine ⟨(1 / 2) * max N₁ N₂, fun t => ?_⟩
  have key : ∀ α : ℝ, |α| < π / 4 → ‖completedRiemannZeta (1 / 2 + t * I)‖ * rexp (α * t) ≤
      (1 / 2) * ∫ v : ℝ, ‖Hk (v + (-2 * α : ℝ) * I)‖ := by
    intro α hα
    have := congrArg norm (Lam_mul_exp hα t)
    rw [norm_mul, norm_mul, Complex.norm_real, Real.norm_eq_abs, abs_of_pos (Real.exp_pos _)]
      at this
    rw [this]
    simp only [one_div, norm_inv, Complex.norm_ofNat]
    gcongr
    exact norm_fourier_le _
  have hN₁ : 0 ≤ N₁ := integral_nonneg (fun _ => norm_nonneg _)
  have hN₂ : 0 ≤ N₂ := integral_nonneg (fun _ => norm_nonneg _)
  rcases le_total 0 t with ht | ht
  · have := key γ h1
    rw [abs_of_nonneg ht]
    have hpos := Real.exp_pos (γ * t)
    rw [show -γ * t = -(γ * t) by ring, Real.exp_neg, ← div_eq_mul_inv, le_div_iff₀ hpos]
    refine this.trans ?_
    gcongr; exact le_max_left _ _
  · have := key (-γ) h2
    rw [abs_of_nonpos ht]
    have hpos := Real.exp_pos (-γ * t)
    rw [show -γ * -t = -(-γ * t) by ring, Real.exp_neg, ← div_eq_mul_inv, le_div_iff₀ hpos]
    refine this.trans ?_
    gcongr; exact le_max_right _ _

end Fayed.Hardy
