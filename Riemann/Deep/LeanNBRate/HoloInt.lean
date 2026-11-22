import Riemann.Deep.LeanNBRate.Contour8
import Mathlib.Analysis.Calculus.ParametricIntegral

/-!
# Holomorphy of parametric integrals (Cauchy estimates for the derivative bound)
-/

open Complex Filter Topology Set MeasureTheory intervalIntegral Metric
open scoped Real Interval

noncomputable section

namespace Riemann.Deep.LeanNBRate

/-- **Holomorphic parametric integrals**: if F(., w) is holomorphic on a closed ball of radius 2r, F(s, .)
is measurable there, and |F| <= g with g integrable, then s -> int F(s, w) dmu(w) is complex differentiable
at the centre. -/
theorem differentiableAt_integral_of_holo {α : Type*} [MeasurableSpace α] {μ : Measure α}
    {F : ℂ → α → ℂ} {s₀ : ℂ} {r : ℝ} (hr : 0 < r) {g : α → ℝ} (hg : Integrable g μ)
    (hd : ∀ᵐ w ∂μ, DifferentiableOn ℂ (fun s => F s w) (closedBall s₀ (2 * r)))
    (hm : ∀ s ∈ closedBall s₀ (2 * r), AEStronglyMeasurable (F s) μ)
    (hb : ∀ᵐ w ∂μ, ∀ s ∈ closedBall s₀ (2 * r), ‖F s w‖ ≤ g w) :
    DifferentiableAt ℂ (fun s => ∫ w, F s w ∂μ) s₀ := by
  set F' : ℂ → α → ℂ := fun s w => deriv (fun s => F s w) s with hF'
  have hsub : ∀ s ∈ ball s₀ r, closedBall s r ⊆ closedBall s₀ (2 * r) := by
    intro s hs x hx
    rw [mem_closedBall] at hx ⊢; rw [mem_ball] at hs
    calc dist x s₀ ≤ dist x s + dist s s₀ := dist_triangle _ _ _
      _ ≤ r + r := by linarith
      _ = 2 * r := by ring
  have hnhds : ∀ s ∈ ball s₀ r, closedBall s₀ (2 * r) ∈ 𝓝 s := by
    intro s hs
    refine mem_of_superset (ball_mem_nhds s hr) ((ball_subset_closedBall).trans (hsub s hs))
  have hdiff : ∀ᵐ w ∂μ, ∀ s ∈ ball s₀ r, HasDerivAt (fun s => F s w) (F' s w) s := by
    filter_upwards [hd] with w hw
    intro s hs
    exact (hw.differentiableAt (hnhds s hs)).hasDerivAt
  have hbound : ∀ᵐ w ∂μ, ∀ s ∈ ball s₀ r, ‖F' s w‖ ≤ g w / r := by
    filter_upwards [hd, hb] with w hw hbw
    intro s hs
    have hdc : DiffContOnCl ℂ (fun s => F s w) (ball s r) := by
      refine DifferentiableOn.diffContOnCl ?_
      rw [closure_ball s hr.ne']
      exact hw.mono (hsub s hs)
    refine Complex.norm_deriv_le_of_forall_mem_sphere_norm_le hr hdc fun z hz => ?_
    exact hbw z (hsub s hs (sphere_subset_closedBall hz))
  have hs₀ : s₀ ∈ closedBall s₀ (2 * r) := mem_closedBall_self (by linarith)
  have hF_int : Integrable (F s₀) μ :=
    hg.mono' (hm s₀ hs₀) (hb.mono fun w hw => hw s₀ hs₀)
  have hF_meas : ∀ᶠ s in 𝓝 s₀, AEStronglyMeasurable (F s) μ := by
    filter_upwards [ball_mem_nhds s₀ (by linarith : 0 < 2 * r)] with s hs
    exact hm s (ball_subset_closedBall hs)
  have hF'_meas : AEStronglyMeasurable (F' s₀) μ := by
    set h : ℕ → ℂ := fun n => ((r / (n + 1) : ℝ) : ℂ) with hh
    have hh0 : ∀ n, h n ≠ 0 := fun n => by
      simp only [hh]; exact_mod_cast (by positivity : r / (n + 1) ≠ 0)
    have hhmem : ∀ n, s₀ + h n ∈ closedBall s₀ (2 * r) := by
      intro n
      rw [mem_closedBall, dist_eq_norm, add_sub_cancel_left, hh, Complex.norm_real, Real.norm_eq_abs,
        abs_of_pos (by positivity)]
      have : r / (n + 1) ≤ r := div_le_self hr.le (by linarith [(Nat.cast_nonneg n : (0:ℝ) ≤ n)])
      linarith
    refine aestronglyMeasurable_of_tendsto_ae atTop
      (f := fun n w => (h n)⁻¹ • (F (s₀ + h n) w - F s₀ w)) (fun n => ?_) ?_
    · exact ((hm _ (hhmem n)).sub (hm s₀ hs₀)).const_smul (h n)⁻¹
    · filter_upwards [hdiff] with w hw
      have hD := (hw s₀ (mem_ball_self hr)).tendsto_slope_zero
      have hT : Tendsto h atTop (𝓝[≠] 0) := by
        refine tendsto_nhdsWithin_iff.mpr ⟨?_, Filter.Eventually.of_forall fun n => hh0 n⟩
        have : Tendsto (fun n : ℕ => r / ((n : ℝ) + 1)) atTop (𝓝 0) := by
          have := (tendsto_one_div_add_atTop_nhds_zero_nat (𝕜 := ℝ)).const_mul r
          simpa [div_eq_mul_inv] using this
        have h1 := (Complex.continuous_ofReal.tendsto 0).comp this
        rw [Complex.ofReal_zero] at h1
        exact h1
      exact hD.comp hT
  have h := hasDerivAt_integral_of_dominated_loc_of_deriv_le (μ := μ) (F := F) (F' := F') (x₀ := s₀)
    (bound := fun w => g w / r) (ball_mem_nhds s₀ hr) hF_meas hF_int hF'_meas
    hbound (hg.div_const r) hdiff
  exact h.2.differentiableAt

end Riemann.Deep.LeanNBRate

