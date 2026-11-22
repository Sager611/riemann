import Riemann.Deep.LeanNBRate.HoloInt

/-!
# The divided-difference correction K (analytic off s = 1)

K(s) = (2 pi i)^{-1} (UpperU - LowerU)(w -> H(w) e^{-L(s-w)} zeta[s,w]/zeta(w)), zeta[s,w] = dslope zeta w s.
-/

open Complex Filter Topology Set MeasureTheory intervalIntegral Metric
open scoped Real Interval

noncomputable section

namespace Riemann.Deep.LeanNBRate

variable {H : ℂ → ℂ} {M β D : ℝ} {N : ℕ}

def kint (H : ℂ → ℂ) (L : ℝ) (s w : ℂ) : ℂ := H w * kap L (s - w) * dslope riemannZeta w s / riemannZeta w

def Kf (H : ℂ → ℂ) (L T₂ : ℝ) (s : ℂ) : ℂ :=
  (2 * π * I)⁻¹ * (UpperU (kint H L s) (-1) 2 T₂ - LowerU (kint H L s) (-1) 2 T₂)

theorem dslope_comm' (f : ℂ → ℂ) (a b : ℂ) : dslope f a b = dslope f b a := by
  by_cases h : a = b
  · rw [h]
  · rw [dslope_of_ne _ (Ne.symm h), dslope_of_ne _ h, slope_comm]

theorem differentiableOn_dslope_zeta {w : ℂ} (hw : w ≠ 1) :
    DifferentiableOn ℂ (dslope riemannZeta w) {z | z ≠ 1} := by
  have hU : {z : ℂ | z ≠ 1} ∈ 𝓝 w := isOpen_ne.mem_nhds hw
  exact (differentiableOn_dslope hU).mpr fun z hz => (differentiableAt_riemannZeta hz).differentiableWithinAt

theorem differentiableOn_deriv_zeta : DifferentiableOn ℂ (deriv riemannZeta) {z | z ≠ 1} := by
  have h : DifferentiableOn ℂ riemannZeta {z | z ≠ 1} := fun z hz =>
    (differentiableAt_riemannZeta hz).differentiableWithinAt
  exact h.deriv isOpen_ne

/-- Analyticity of one contour piece of K. -/
theorem differentiableAt_piece {L : ℝ} (hL : 0 ≤ L) {S : Set ℝ} (hS : MeasurableSet S) {γ : ℝ → ℂ}
    (hγ : Continuous γ) {B : ℝ} (hB : 0 ≤ B) (hγ1 : ∀ y ∈ S, γ y ≠ 1) (hγz : ∀ y ∈ S, riemannZeta (γ y) ≠ 0)
    (hγB : ∀ y ∈ S, ‖(riemannZeta (γ y))⁻¹‖ ≤ B) (hγre : ∀ y ∈ S, (γ y).re ≤ 2)
    (hHc : Continuous H) (hHi : IntegrableOn (fun y => ‖H (γ y)‖) S) {s₀ : ℂ} (hs₀ : s₀ ≠ 1) :
    DifferentiableAt ℂ (fun s => ∫ y in S, kint H L s (γ y)) s₀ := by
  set r := ‖s₀ - 1‖ / 4 with hr
  have hpos := norm_pos_iff.mpr (sub_ne_zero.mpr hs₀)
  have hr0 : 0 < r := by rw [hr]; positivity
  have hball : ∀ z ∈ closedBall s₀ (3 * r), z ≠ 1 := by
    intro z hz h
    rw [h, mem_closedBall, dist_comm, dist_eq_norm] at hz
    rw [hr] at hz
    linarith
  have hball2 : closedBall s₀ (2 * r) ⊆ closedBall s₀ (3 * r) := closedBall_subset_closedBall (by linarith)
  obtain ⟨Cz, hCz⟩ := (isCompact_closedBall s₀ (2 * r)).exists_bound_of_continuousOn
    (f := riemannZeta) fun z hz => (differentiableAt_riemannZeta (hball z (hball2 hz))).continuousAt.continuousWithinAt
  obtain ⟨Cd, hCd⟩ := (isCompact_closedBall s₀ (3 * r)).exists_bound_of_continuousOn
    (f := deriv riemannZeta) fun z hz =>
      (differentiableOn_deriv_zeta.differentiableAt (isOpen_ne.mem_nhds (hball z hz))).continuousAt.continuousWithinAt
  have hCd0 : 0 ≤ Cd := (norm_nonneg _).trans (hCd s₀ (mem_closedBall_self (by linarith)))
  have hCz0 : 0 ≤ Cz := (norm_nonneg _).trans (hCz s₀ (mem_closedBall_self (by linarith)))
  set E₀ := Real.exp (L * (2 + ‖s₀‖ + 2 * r)) with hE₀
  set C := max ((Cz * B + 1) / r) (Cd * B) with hC
  -- the pointwise bound
  have hkb : ∀ y ∈ S, ∀ s ∈ closedBall s₀ (2 * r), ‖kint H L s (γ y)‖ ≤ ‖H (γ y)‖ * E₀ * C := by
    intro y hy s hs
    set w := γ y with hw
    have hk : ‖kap L (s - w)‖ ≤ E₀ := by
      rw [norm_kap, hE₀]
      apply Real.exp_le_exp.mpr
      have h1 : (s - w).re = s.re - w.re := by simp
      rw [h1]
      have h2 : -s.re ≤ ‖s₀‖ + 2 * r := by
        have := neg_abs_le s.re
        have h3 : |s.re| ≤ ‖s‖ := Complex.abs_re_le_norm s
        have h4 : ‖s‖ ≤ ‖s₀‖ + 2 * r := by
          rw [mem_closedBall, dist_eq_norm] at hs
          calc ‖s‖ = ‖s₀ + (s - s₀)‖ := by ring_nf
            _ ≤ ‖s₀‖ + ‖s - s₀‖ := norm_add_le _ _
            _ ≤ ‖s₀‖ + 2 * r := by linarith
        linarith
      have := hγre y hy
      nlinarith
    have hds : ‖dslope riemannZeta w s‖ * ‖(riemannZeta w)⁻¹‖ ≤ C := by
      by_cases hsw : r ≤ ‖s - w‖
      · have hne : s ≠ w := fun h => by rw [h, sub_self, norm_zero] at hsw; linarith
        rw [dslope_of_ne _ hne, slope_def_field, norm_div, div_mul_eq_mul_div, ← norm_mul]
        have e : (riemannZeta s - riemannZeta w) * (riemannZeta w)⁻¹ = riemannZeta s * (riemannZeta w)⁻¹ - 1 := by
          rw [sub_mul, mul_inv_cancel₀ (hγz y hy)]
        rw [e]
        refine le_trans ?_ (le_max_left _ _)
        rw [div_le_div_iff₀ (lt_of_lt_of_le hr0 hsw) hr0]
        have h1 : ‖riemannZeta s * (riemannZeta w)⁻¹ - 1‖ ≤ Cz * B + 1 := by
          calc ‖riemannZeta s * (riemannZeta w)⁻¹ - 1‖ ≤ ‖riemannZeta s * (riemannZeta w)⁻¹‖ + ‖(1 : ℂ)‖ :=
                norm_sub_le _ _
            _ ≤ Cz * B + 1 := by
                rw [norm_mul, norm_one]
                have := hCz s hs
                have := hγB y hy
                gcongr
        nlinarith [norm_nonneg (riemannZeta s * (riemannZeta w)⁻¹ - 1)]
      · push Not at hsw
        have hw3 : w ∈ closedBall s₀ (3 * r) := by
          rw [mem_closedBall] at hs ⊢
          calc dist w s₀ ≤ dist w s + dist s s₀ := dist_triangle _ _ _
            _ ≤ r + 2 * r := by rw [dist_comm, dist_eq_norm]; linarith
            _ = 3 * r := by ring
        have hdl : ‖dslope riemannZeta w s‖ ≤ Cd := by
          by_cases hne : s = w
          · rw [hne, dslope_same]; exact hCd w hw3
          · rw [dslope_of_ne _ hne, slope_def_field, norm_div]
            rw [div_le_iff₀ (norm_pos_iff.mpr (sub_ne_zero.mpr hne))]
            exact Convex.norm_image_sub_le_of_norm_deriv_le
              (fun z hz => differentiableAt_riemannZeta (hball z hz)) hCd (convex_closedBall s₀ (3 * r))
              hw3 (hball2 hs)
        refine le_trans ?_ (le_max_right _ _)
        exact mul_le_mul hdl (hγB y hy) (norm_nonneg _) hCd0
    rw [kint, norm_div, norm_mul, norm_mul, div_eq_mul_inv, ← norm_inv]
    calc ‖H w‖ * ‖kap L (s - w)‖ * ‖dslope riemannZeta w s‖ * ‖(riemannZeta w)⁻¹‖
        = ‖H w‖ * ‖kap L (s - w)‖ * (‖dslope riemannZeta w s‖ * ‖(riemannZeta w)⁻¹‖) := by ring
      _ ≤ ‖H w‖ * E₀ * C := by gcongr
  refine differentiableAt_integral_of_holo (μ := volume.restrict S) (F := fun s y => kint H L s (γ y))
    (g := fun y => ‖H (γ y)‖ * E₀ * C) hr0 ((hHi.mul_const E₀).mul_const C) ?_ ?_ ?_
  · filter_upwards [ae_restrict_mem hS] with y hy
    intro s hs
    have hw1 := hγ1 y hy
    have hd : DifferentiableAt ℂ (fun s => kint H L s (γ y)) s := by
      unfold kint kap
      have h1 : DifferentiableAt ℂ (dslope riemannZeta (γ y)) s :=
        (differentiableOn_dslope_zeta hw1).differentiableAt (isOpen_ne.mem_nhds (hball s (hball2 hs)))
      exact (((differentiableAt_const _).mul (by fun_prop)).mul h1).div_const _
    exact hd.differentiableWithinAt
  · intro s hs
    refine ContinuousOn.aestronglyMeasurable (fun y hy => ?_) hS
    have hs1 : s ≠ 1 := hball s (hball2 hs)
    have hdsl : ContinuousAt (fun y => dslope riemannZeta (γ y) s) y := by
      have e : (fun y => dslope riemannZeta (γ y) s) = fun y => dslope riemannZeta s (γ y) :=
        funext fun y => dslope_comm' _ _ _
      rw [e]
      refine ContinuousAt.comp (f := γ) (g := dslope riemannZeta s) ?_ hγ.continuousAt
      by_cases h : γ y = s
      · rw [h]; exact continuousAt_dslope_same.mpr (differentiableAt_riemannZeta hs1)
      · exact (continuousAt_dslope_of_ne h).mpr (differentiableAt_riemannZeta (hγ1 y hy)).continuousAt
    have hz : ContinuousAt (fun y => riemannZeta (γ y)) y :=
      ContinuousAt.comp (f := γ) (differentiableAt_riemannZeta (hγ1 y hy)).continuousAt hγ.continuousAt
    refine ContinuousAt.continuousWithinAt ?_
    unfold kint kap
    exact (((hHc.comp hγ).continuousAt.mul (by fun_prop)).mul hdsl).div hz (hγz y hy)
  · filter_upwards [ae_restrict_mem hS] with y hy
    exact hkb y hy

end Riemann.Deep.LeanNBRate
