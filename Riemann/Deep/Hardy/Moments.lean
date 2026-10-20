import Riemann.Deep.Hardy.Shift

/-!
# Hardy's theorem, part 4: moments of Λ on the critical line

J m α = ∫ Λ(1/2+it) t^m e^{αt} dt. For |α| < π/4:
* J_zero: J 0 α = 2π Hk(-2αi)   (Fourier inversion at 0 of the shifted kernel line);
* J_hasDerivAt: d/dα J m = J (m+1)   (differentiation under the integral, via Lam_decay);
* J_eq: J m α = 2π (-2i)^m Hk^{(m)}(-2αi).
-/

open Complex Real MeasureTheory Set Filter Topology
open scoped ComplexConjugate FourierTransform Interval Nat

namespace Fayed.Hardy

/-- The moments. -/
noncomputable def J (m : ℕ) (α : ℝ) : ℂ :=
  ∫ t : ℝ, completedRiemannZeta (1 / 2 + t * I) * (t : ℂ) ^ m * (rexp (α * t) : ℂ)

lemma continuous_Lam : Continuous (fun t : ℝ => completedRiemannZeta (1 / 2 + t * I)) := by
  rw [continuous_iff_continuousAt]
  intro t
  have h0 : (1 / 2 + t * I : ℂ) ≠ 0 := by
    intro h; have := congrArg Complex.re h; simp at this
  have h1 : (1 / 2 + t * I : ℂ) ≠ 1 := by
    intro h; have := congrArg Complex.re h; simp at this
  exact ContinuousAt.comp (g := completedRiemannZeta)
    (differentiableAt_completedZeta h0 h1).continuousAt (by fun_prop)

lemma continuous_moment (m : ℕ) (α : ℝ) : Continuous (fun t : ℝ =>
    completedRiemannZeta (1 / 2 + t * I) * (t : ℂ) ^ m * (rexp (α * t) : ℂ)) := by
  have := continuous_Lam; fun_prop

/-- |t|^k e^{-b|t|} ≤ k! (2/b)^k e^{-(b/2)|t|}. -/
lemma pow_abs_mul_exp_le {b : ℝ} (hb : 0 < b) (k : ℕ) (t : ℝ) :
    |t| ^ k * rexp (-b * |t|) ≤ (k ! * (2 / b) ^ k) * rexp (-(b / 2) * |t|) := by
  have hx : 0 ≤ b / 2 * |t| := by positivity
  have h := Real.pow_div_factorial_le_exp _ hx k
  have hk : (0 : ℝ) < k ! := by exact_mod_cast Nat.factorial_pos k
  rw [div_le_iff₀ hk] at h
  have h3 : (b / 2 * |t|) ^ k * rexp (-(b / 2 * |t|)) ≤ k ! := by
    rw [Real.exp_neg, ← div_eq_mul_inv, div_le_iff₀ (Real.exp_pos _)]; linarith
  have h2 : |t| ^ k = (2 / b) ^ k * (b / 2 * |t|) ^ k := by
    rw [← mul_pow, show 2 / b * (b / 2 * |t|) = |t| by field_simp]
  have he : rexp (-b * |t|) = rexp (-(b / 2 * |t|)) * rexp (-(b / 2) * |t|) := by
    rw [← Real.exp_add]; congr 1; ring
  have hpos : 0 ≤ (2 / b) ^ k := by positivity
  calc |t| ^ k * rexp (-b * |t|)
      = (2 / b) ^ k * ((b / 2 * |t|) ^ k * rexp (-(b / 2 * |t|))) * rexp (-(b / 2) * |t|) := by
        rw [h2, he]; ring
    _ ≤ (2 / b) ^ k * k ! * rexp (-(b / 2) * |t|) := by gcongr
    _ = (k ! * (2 / b) ^ k) * rexp (-(b / 2) * |t|) := by ring

lemma integrable_pow_abs_mul_exp {b : ℝ} (hb : 0 < b) (k : ℕ) :
    Integrable (fun t : ℝ => |t| ^ k * rexp (-b * |t|)) := by
  refine Integrable.mono' ((integrable_exp_neg_abs (half_pos hb)).const_mul (k ! * (2 / b) ^ k))
    (by fun_prop) (Eventually.of_forall fun t => ?_)
  rw [Real.norm_eq_abs, abs_of_nonneg (by positivity)]
  exact pow_abs_mul_exp_le hb k t

/-- Pointwise bound for the moment integrand. -/
lemma norm_moment_le {C γ α : ℝ} (hC : ∀ t : ℝ, ‖completedRiemannZeta (1 / 2 + t * I)‖ ≤
    C * rexp (-γ * |t|)) (m : ℕ) (t : ℝ) (hα : |α| ≤ γ) :
    ‖completedRiemannZeta (1 / 2 + t * I) * (t : ℂ) ^ m * (rexp (α * t) : ℂ)‖ ≤
      C * (|t| ^ m * rexp (-(γ - |α|) * |t|)) := by
  rw [norm_mul, norm_mul, norm_pow, Complex.norm_real, Complex.norm_real, Real.norm_eq_abs,
    Real.norm_eq_abs, abs_of_pos (Real.exp_pos _)]
  have h1 : rexp (α * t) ≤ rexp (|α| * |t|) := Real.exp_le_exp.mpr (by
    rw [← abs_mul]; exact le_abs_self _)
  have h2 : rexp (-γ * |t|) * rexp (|α| * |t|) = rexp (-(γ - |α|) * |t|) := by
    rw [← Real.exp_add]; congr 1; ring
  have hC0 : 0 ≤ C * rexp (-γ * |t|) := (norm_nonneg _).trans (hC t)
  calc ‖completedRiemannZeta (1 / 2 + t * I)‖ * |t| ^ m * rexp (α * t)
      ≤ (C * rexp (-γ * |t|)) * |t| ^ m * rexp (|α| * |t|) := by gcongr; exact hC t
    _ = C * (|t| ^ m * (rexp (-γ * |t|) * rexp (|α| * |t|))) := by ring
    _ = _ := by rw [h2]

lemma integrable_moment {α : ℝ} (hα : |α| < π / 4) (m : ℕ) :
    Integrable (fun t : ℝ => completedRiemannZeta (1 / 2 + t * I) * (t : ℂ) ^ m *
      (rexp (α * t) : ℂ)) := by
  set γ := (|α| + π / 4) / 2
  have hγ0 : 0 ≤ γ := by positivity
  have hγ : γ < π / 4 := by simp only [γ]; linarith
  have hαγ : |α| < γ := by simp only [γ]; linarith
  obtain ⟨C, hC⟩ := Lam_decay hγ0 hγ
  refine Integrable.mono' ((integrable_pow_abs_mul_exp (sub_pos.mpr hαγ) m).const_mul C)
    (continuous_moment m α).aestronglyMeasurable (Eventually.of_forall fun t => ?_)
  exact norm_moment_le hC m t hαγ.le

/-- Fourier inversion: the zeroth moment. -/
theorem J_zero {α : ℝ} (hα : |α| < π / 4) : J 0 α = 2 * π * Hk (-2 * α * I) := by
  have hc : |(-2 * α : ℝ)| < π / 2 := by rw [abs_mul]; norm_num; linarith
  set h : ℝ → ℂ := fun v => Hk (v + (-2 * α : ℝ) * I) with hh
  have hF : 𝓕 h = fun w => 2 * completedRiemannZeta (1 / 2 + ((4 * π * w : ℝ) : ℂ) * I) *
      (rexp (α * (4 * π * w)) : ℂ) := by
    ext w
    have := Lam_mul_exp hα (4 * π * w)
    rw [show 4 * π * w / (4 * π) = w by field_simp] at this
    calc 𝓕 h w = 2 * (1 / 2 * 𝓕 h w) := by ring
      _ = _ := by rw [← this]; ring
  have hFint : Integrable (𝓕 h) := by
    rw [hF]
    have := ((integrable_moment hα 0).comp_mul_left' (by positivity : (4 * π : ℝ) ≠ 0)).const_mul 2
    refine this.congr (Eventually.of_forall fun w => ?_)
    simp only [pow_zero, mul_one]; ring
  have hinv := (integrable_Hk_line hc).fourierInv_fourier_eq hFint
    (continuous_Hk_line hc).continuousAt (v := 0)
  rw [Real.fourierInv_eq] at hinv
  simp only [inner_zero_right, AddChar.map_zero_eq_one, one_smul] at hinv
  have hcomp : ∫ w : ℝ, 𝓕 h w = 2 * ((1 / (4 * π) : ℝ) : ℂ) * J 0 α := by
    rw [hF]
    have := Measure.integral_comp_mul_left (fun t : ℝ => 2 * completedRiemannZeta (1 / 2 + (t : ℂ) * I) *
      (rexp (α * t) : ℂ)) (4 * π)
    rw [this, abs_of_pos (by positivity : (0 : ℝ) < (4 * π)⁻¹), J]
    simp only [pow_zero, mul_one, mul_assoc (2 : ℂ)]
    rw [integral_const_mul, Complex.real_smul]
    push_cast; ring
  rw [hcomp] at hinv
  have hv : Hk (((0 : ℝ) : ℂ) + ((-2 * α : ℝ) : ℂ) * I) = Hk (-2 * α * I) := by
    push_cast; ring_nf
  try rw [hv] at hinv
  have hπ : (π : ℂ) ≠ 0 := ofReal_ne_zero.mpr Real.pi_pos.ne'
  rw [← hinv]; push_cast; field_simp; ring

/-- Differentiation under the integral sign. -/
theorem J_hasDerivAt {α₀ : ℝ} (hα₀ : |α₀| < π / 4) (m : ℕ) :
    HasDerivAt (J m) (J (m + 1) α₀) α₀ := by
  set δ := π / 4 - |α₀| with hδ
  have hδ0 : 0 < δ := by simp only [δ]; linarith
  set γ := |α₀| + 2 * δ / 3
  have hγ0 : 0 ≤ γ := by positivity
  have hγ : γ < π / 4 := by simp only [γ, δ]; linarith
  obtain ⟨C, hC⟩ := Lam_decay hγ0 hγ
  have hball : ∀ α ∈ Metric.ball α₀ (δ / 3), |α| ≤ |α₀| + δ / 3 := by
    intro α hα
    rw [Metric.mem_ball, Real.dist_eq] at hα
    have := abs_sub_abs_le_abs_sub α α₀
    linarith
  have hball' : ∀ α ∈ Metric.ball α₀ (δ / 3), |α| < π / 4 := by
    intro α hα; have := hball α hα; simp only [δ] at this; linarith
  have key := hasDerivAt_integral_of_dominated_loc_of_deriv_le (μ := volume)
    (F := fun α t => completedRiemannZeta (1 / 2 + t * I) * (t : ℂ) ^ m * (rexp (α * t) : ℂ))
    (F' := fun α t => completedRiemannZeta (1 / 2 + t * I) * (t : ℂ) ^ (m + 1) *
      (rexp (α * t) : ℂ))
    (bound := fun t => C * (|t| ^ (m + 1) * rexp (-(δ / 3) * |t|)))
    (x₀ := α₀) (Metric.ball_mem_nhds α₀ (by positivity : 0 < δ / 3))
    (Eventually.of_forall fun α => (continuous_moment m α).aestronglyMeasurable)
    (integrable_moment hα₀ m)
    ((continuous_moment (m + 1) α₀).aestronglyMeasurable)
    (Eventually.of_forall fun t α hα => by
      have h1 := norm_moment_le hC (m + 1) t (show |α| ≤ γ by
        have := hball α hα; simp only [γ]; linarith)
      have h2 : -(γ - |α|) ≤ -(δ / 3) := by have := hball α hα; simp only [γ]; linarith
      refine h1.trans ?_
      have hC0 : 0 ≤ C := by
        have := (norm_nonneg _).trans (hC 0); simpa using this
      gcongr)
    ((integrable_pow_abs_mul_exp (by positivity) (m + 1)).const_mul C)
    (Eventually.of_forall fun t α _ => by
      have h1 : HasDerivAt (fun α : ℝ => rexp (α * t)) (rexp (α * t) * t) α := by
        have := ((hasDerivAt_id α).mul_const t).exp
        simpa using this
      have h2 := (h1.ofReal_comp).const_mul (completedRiemannZeta (1 / 2 + t * I) * (t : ℂ) ^ m)
      convert h2 using 1
      push_cast; ring)
  exact key.2

lemma analyticOnNhd_Hk : AnalyticOnNhd ℂ Hk strip :=
  differentiableOn_Hk.analyticOnNhd isOpen_strip

lemma differentiableAt_iteratedDeriv_Hk (m : ℕ) {z : ℂ} (hz : z ∈ strip) :
    DifferentiableAt ℂ (iteratedDeriv m Hk) z := by
  rw [iteratedDeriv_eq_iterate]
  exact (analyticOnNhd_Hk.iterated_deriv m z hz).differentiableAt

/-- All moments in terms of derivatives of the kernel. -/
theorem J_eq (m : ℕ) : ∀ {α : ℝ}, |α| < π / 4 →
    J m α = 2 * π * (-2 * I) ^ m * iteratedDeriv m Hk (-2 * α * I) := by
  induction m with
  | zero => intro α hα; rw [J_zero hα]; simp
  | succ m ih =>
    intro α hα
    have hmem : ∀ β : ℝ, |β| < π / 4 → (-2 * (β : ℂ) * I) ∈ strip := by
      intro β hβ
      show |(-2 * (β : ℂ) * I).im| < π / 2
      simp [abs_mul]; linarith
    have hD := (differentiableAt_iteratedDeriv_Hk m (hmem α hα)).hasDerivAt
    rw [← iteratedDeriv_succ] at hD
    have hlin : HasDerivAt (fun z : ℂ => -2 * z * I) (-2 * I) (α : ℂ) := by
      simpa using ((hasDerivAt_id (α : ℂ)).const_mul (-2)).mul_const I
    have hcomp := (hD.comp (α : ℂ) hlin).comp_ofReal
    have hg := hcomp.const_mul (2 * π * (-2 * I) ^ m)
    have hopen : IsOpen {β : ℝ | |β| < π / 4} := isOpen_lt continuous_abs continuous_const
    have heq : (fun β : ℝ => 2 * π * (-2 * I) ^ m * iteratedDeriv m Hk (-2 * (β : ℂ) * I))
        =ᶠ[𝓝 α] J m := by
      filter_upwards [hopen.mem_nhds hα] with β hβ
      exact (ih hβ).symm
    have hJ := J_hasDerivAt hα m
    have := (hg.congr_of_eventuallyEq heq.symm).unique hJ
    rw [← this]
    ring

end Fayed.Hardy
