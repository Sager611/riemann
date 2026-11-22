import Riemann.Deep.LeanNBRate.KJ4
import Riemann.Deep.LeanNBRate.Dual5

/-!
# Analyticity of R_tr and G, and the identity theorem: the formula on the critical line
-/

open Complex Filter Topology Set MeasureTheory intervalIntegral Metric
open scoped Real Interval Nat

noncomputable section

namespace Riemann.Deep.LeanNBRate

variable {H : ℂ → ℂ} {M β D : ℝ} {N : ℕ}

theorem norm_rk_le (hH : HBound H M β D N) (j : ℕ) :
    ‖rk H j‖ ≤ M * Real.exp (β * (2 * j + 5 / 2)) * (6 * (2 * π) ^ (2 * j + 2) / (2 * j + 2).factorial) := by
  have hM := hH.M_nonneg
  have hD := hH.D_pos
  have htz : (tz j).re = -2 * j - 2 := tz_re j
  have hj0 : (0 : ℝ) ≤ j := Nat.cast_nonneg j
  have hHb := hH.bound (tz j)
  have hw : (D / (D + ‖tz j - 1 / 2‖)) ^ N ≤ 1 :=
    pow_le_one₀ (by positivity) (div_le_one_of_le₀ (by linarith [norm_nonneg (tz j - 1 / 2)]) (by positivity))
  have hmax : max 0 (1 / 2 - (tz j).re) = 2 * j + 5 / 2 := by
    rw [htz, max_eq_right (by linarith)]; ring
  rw [hmax] at hHb
  have hd := norm_deriv_zeta_neg_two_ge j
  have hfpos : (0 : ℝ) < ((2 * j + 2).factorial : ℝ) := by exact_mod_cast Nat.factorial_pos _
  have hdpos : (0 : ℝ) < ((2 * j + 2).factorial : ℝ) / (6 * (2 * π) ^ (2 * j + 2)) := by positivity
  rw [rk, norm_div]
  have h1 : ‖H (-2 * ((j : ℂ) + 1))‖ ≤ M * Real.exp (β * (2 * j + 5 / 2)) := by
    have := hHb; rw [tz] at this
    calc ‖H (-2 * ((j : ℂ) + 1))‖ ≤ M * Real.exp (β * (2 * j + 5 / 2)) * (D / (D + ‖tz j - 1 / 2‖)) ^ N := this
      _ ≤ M * Real.exp (β * (2 * j + 5 / 2)) * 1 := by gcongr
      _ = _ := mul_one _
  rw [div_le_iff₀ (lt_of_lt_of_le hdpos hd)]
  calc ‖H (-2 * ((j : ℂ) + 1))‖ ≤ M * Real.exp (β * (2 * j + 5 / 2)) := h1
    _ = M * Real.exp (β * (2 * j + 5 / 2)) * (6 * (2 * π) ^ (2 * j + 2) / (2 * j + 2).factorial) *
        (((2 * j + 2).factorial : ℝ) / (6 * (2 * π) ^ (2 * j + 2))) := by field_simp
    _ ≤ _ := by gcongr

/-- Uniform bound on Re s > -1. -/
theorem norm_Ares_le_strip (hH : HBound H M β D N) {L : ℝ} (hL : 0 ≤ L) {s : ℂ} (hs : -1 < s.re) (j : ℕ) :
    ‖Ares H L s j‖ ≤ M * Real.exp (β * (2 * j + 5 / 2)) * (6 * (2 * π) ^ (2 * j + 2) / j !) := by
  have hj0 : (0 : ℝ) ≤ j := Nat.cast_nonneg j
  have htz : (tz j).re = -2 * j - 2 := tz_re j
  have hsz : 1 ≤ ‖s - tz j‖ := by
    have h1 := Complex.abs_re_le_norm (s - tz j)
    rw [sub_re, htz] at h1
    have h2 : (0 : ℝ) < s.re - (-2 * j - 2) := by linarith
    rw [abs_of_pos h2] at h1
    linarith
  have hk : ‖kap L (s - tz j)‖ ≤ 1 := norm_kap_le_one hL (by rw [sub_re, htz]; linarith)
  have hM := hH.M_nonneg
  have hfac : ((j ! : ℕ) : ℝ) ≤ ((2 * j + 2).factorial : ℝ) := by exact_mod_cast Nat.factorial_le (by omega)
  have hjpos : (0 : ℝ) < (j ! : ℝ) := by exact_mod_cast Nat.factorial_pos _
  rw [Ares, norm_div, norm_mul]
  calc ‖rk H j‖ * ‖kap L (s - tz j)‖ / ‖s - tz j‖ ≤ ‖rk H j‖ * ‖kap L (s - tz j)‖ := div_le_self (by positivity) hsz
    _ ≤ ‖rk H j‖ * 1 := mul_le_mul_of_nonneg_left hk (norm_nonneg _)
    _ ≤ M * Real.exp (β * (2 * j + 5 / 2)) * (6 * (2 * π) ^ (2 * j + 2) / (2 * j + 2).factorial) := by
        rw [mul_one]; exact norm_rk_le hH j
    _ ≤ _ := by gcongr

theorem summable_Ares_bound (hH : HBound H M β D N) :
    Summable fun j : ℕ => M * Real.exp (β * (2 * j + 5 / 2)) * (6 * (2 * π) ^ (2 * j + 2) / j !) := by
  set r := Real.exp (2 * β) * (2 * π) ^ 2 with hr
  have h := (Real.summable_pow_div_factorial r).mul_left (M * Real.exp (β * (5 / 2)) * 6 * (2 * π) ^ 2)
  refine h.congr fun j => ?_
  have e1 : r ^ j = Real.exp (j * (2 * β)) * (2 * π) ^ (2 * j) := by
    rw [hr, mul_pow (Real.exp (2 * β)) ((2 * π) ^ 2) j, ← Real.exp_nat_mul, ← pow_mul]
  have e2 : Real.exp (β * (2 * j + 5 / 2)) = Real.exp (β * (5 / 2)) * Real.exp (j * (2 * β)) := by
    rw [← Real.exp_add]; congr 1; ring
  have e3 : (2 * π) ^ (2 * j + 2) = (2 * π) ^ 2 * (2 * π) ^ (2 * j) := by
    rw [← pow_add]; ring_nf
  rw [e1, e2, e3]
  ring

theorem differentiableOn_Rtr (hH : HBound H M β D N) {L : ℝ} (hL : 0 ≤ L) :
    DifferentiableOn ℂ (Rtr H L) {s : ℂ | -1 < s.re} := by
  have e : Rtr H L = fun s => ∑' j : ℕ, Ares H L s j := funext fun s => Rtr_eq H L s
  rw [e]
  refine differentiableOn_tsum_of_summable_norm (summable_Ares_bound hH) (fun j => ?_)
    (isOpen_lt continuous_const Complex.continuous_re) (fun j s hs => norm_Ares_le_strip hH hL hs j)
  intro s hs
  have hne : s - tz j ≠ 0 := by
    intro h; have := congrArg Complex.re h; rw [sub_re, tz_re, Complex.zero_re] at this
    have hs' : -1 < s.re := hs
    linarith [(Nat.cast_nonneg j : (0:ℝ) ≤ j)]
  have hk : DifferentiableAt ℂ (fun s : ℂ => cexp (-(L : ℂ) * (s - tz j))) s := by fun_prop
  have hd : DifferentiableAt ℂ (fun s : ℂ => rk H j * cexp (-(L : ℂ) * (s - tz j)) / (s - tz j)) s :=
    ((differentiableAt_const _).mul hk).div ((differentiableAt_id).sub (differentiableAt_const _)) hne
  exact hd.differentiableWithinAt

theorem continuous_gF (hH : HBound H M β D N) : Continuous (gF H) := by
  have hFi : Integrable (fun y : ℝ => Fq H (3 + y * I)) := by
    have := integrable_Fq_line hH (σ := 3) (by norm_num); simpa using this
  have hFc : Continuous (fun y : ℝ => Fq H (3 + y * I)) := by
    have := continuous_Fq_line hH (σ := 3) (by norm_num); simpa using this
  rw [continuous_iff_continuousAt]
  intro u₀
  unfold gF
  refine ContinuousAt.mul continuousAt_const ?_
  refine continuousAt_of_dominated (bound := fun y => ‖Fq H (3 + y * I)‖ * Real.exp (3 * (|u₀| + 1))) ?_ ?_
    (hFi.norm.mul_const _) ?_
  · exact Filter.Eventually.of_forall fun u => (hFc.mul (by fun_prop)).aestronglyMeasurable
  · filter_upwards [Metric.ball_mem_nhds u₀ one_pos] with u hu
    refine Filter.Eventually.of_forall fun y => ?_
    rw [norm_mul, Complex.norm_exp]
    refine mul_le_mul_of_nonneg_left (Real.exp_le_exp.mpr ?_) (norm_nonneg _)
    have e : ((u : ℂ) * (3 + y * I)).re = 3 * u := by simp; ring
    rw [e]
    rw [Metric.mem_ball, Real.dist_eq] at hu
    have := abs_sub_abs_le_abs_sub u u₀
    have := le_abs_self u
    linarith
  · exact Filter.Eventually.of_forall fun y => by fun_prop

theorem differentiable_GF (hH : HBound H M β D N) {L : ℝ} (hL : 0 ≤ L) : Differentiable ℂ (GF H L) := by
  have hq : Continuous fun u : ℝ => gF H u * cexp (-(u : ℂ) / 2) := (continuous_gF hH).mul (by fun_prop)
  have h := differentiable_laplace hL hq
  have e : GF H L = fun s : ℂ => ∫ u in (0 : ℝ)..L, cexp (-(s - 1 / 2) * u) * (gF H u * cexp (-(u : ℂ) / 2)) := by
    funext s
    unfold GF
    refine intervalIntegral.integral_congr fun u _ => ?_
    show gF H u * cexp (-(u * s)) = cexp (-(s - 1 / 2) * u) * (gF H u * cexp (-(u : ℂ) / 2))
    rw [show gF H u * cexp (-(u * s)) = gF H u * (cexp (-(s - 1 / 2) * u) * cexp (-(u : ℂ) / 2)) by
      rw [← Complex.exp_add]; congr 2; ring]
    ring
  rw [e]
  exact h

/-- **The formula on the critical line** (identity theorem): zeta G = H - K - zeta R_tr. -/
theorem line_formula (hH : HBound H M β D N) {T₂ : ℝ} (hZ : HZeros H T₂) (hT1 : 1 ≤ T₂) {L : ℝ}
    (hL : 0 ≤ L) (hLβ : β + Real.log (2 * π) < L) {s : ℂ} (hs1 : -1 < s.re) (hs2 : s.re < 1) :
    riemannZeta s * GF H L s = H s - Kf H L T₂ s - riemannZeta s * Rtr H L s := by
  set Φ : ℂ → ℂ := fun s => riemannZeta s * GF H L s + riemannZeta s * Rtr H L s + Kf H L T₂ s - H s with hΦ
  have hΦd : DifferentiableOn ℂ Φ {s : ℂ | -1 < s.re ∧ s ≠ 1} := by
    intro z hz
    have hz1 : z ≠ 1 := hz.2
    refine DifferentiableAt.differentiableWithinAt ?_
    have hR : DifferentiableAt ℂ (Rtr H L) z :=
      (differentiableOn_Rtr hH hL).differentiableAt ((isOpen_lt continuous_const Complex.continuous_re).mem_nhds hz.1)
    exact ((((differentiableAt_riemannZeta hz1).mul (differentiable_GF hH hL z)).add
      ((differentiableAt_riemannZeta hz1).mul hR)).add (differentiableAt_Kf hH hZ hT1 hL hz1)).sub (hH.diff z)
  have hΦ3 : ∀ z : ℂ, 3 < z.re → Φ z = 0 := by
    intro z hz
    have hG := GF_eq_main hH hZ hL hLβ hz
    have hK := Kf_eq_zeta_Jf hH hZ hT1 hL hz
    have hz0 : riemannZeta z ≠ 0 := riemannZeta_ne_zero_of_one_lt_re (by linarith)
    have hF : riemannZeta z * Fq H z = H z := by rw [Fq]; field_simp
    simp only [hΦ, hG, hK]
    linear_combination hF
  have hanal : ∀ U : Set ℂ, IsOpen U → U ⊆ {s : ℂ | -1 < s.re ∧ s ≠ 1} → AnalyticOnNhd ℂ Φ U :=
    fun U hU hsub => (hΦd.mono hsub).analyticOnNhd hU
  have hzero : AnalyticOnNhd ℂ (fun _ : ℂ => (0 : ℂ)) univ := fun z _ => analyticAt_const
  -- U1 = {1 < re}
  set U₁ : Set ℂ := {s : ℂ | 1 < s.re}
  have hU₁o : IsOpen U₁ := isOpen_lt continuous_const Complex.continuous_re
  have hU₁s : U₁ ⊆ {s : ℂ | -1 < s.re ∧ s ≠ 1} := fun z hz =>
    ⟨by have : 1 < z.re := hz; linarith, fun h => by have : 1 < z.re := hz; rw [h] at this; simp at this⟩
  have hU₁c : Convex ℝ U₁ := convex_halfSpace_re_gt 1
  have h1 : EqOn Φ (fun _ => 0) U₁ := by
    refine (hanal U₁ hU₁o hU₁s).eqOn_of_preconnected_of_eventuallyEq (hzero.mono (subset_univ _))
      hU₁c.isPreconnected (z₀ := 4) (by show (1 : ℝ) < (4 : ℂ).re; norm_num) ?_
    filter_upwards [(isOpen_lt continuous_const Complex.continuous_re).mem_nhds
      (show (3 : ℝ) < (4 : ℂ).re by norm_num)] with z hz
    exact hΦ3 z hz
  -- U2 = {-1 < re, 0 < im}
  set U₂ : Set ℂ := {s : ℂ | -1 < s.re ∧ 0 < s.im}
  have hU₂o : IsOpen U₂ := (isOpen_lt continuous_const Complex.continuous_re).inter
    (isOpen_lt continuous_const Complex.continuous_im)
  have hU₂s : U₂ ⊆ {s : ℂ | -1 < s.re ∧ s ≠ 1} := fun z hz =>
    ⟨hz.1, fun h => by have := hz.2; rw [h] at this; simp at this⟩
  have hU₂c : Convex ℝ U₂ := (convex_halfSpace_re_gt (-1)).inter (convex_halfSpace_im_gt 0)
  have h2 : EqOn Φ (fun _ => 0) U₂ := by
    refine (hanal U₂ hU₂o hU₂s).eqOn_of_preconnected_of_eventuallyEq (hzero.mono (subset_univ _))
      hU₂c.isPreconnected (z₀ := 2 + I) ⟨by simp; norm_num, by simp⟩ ?_
    filter_upwards [hU₁o.mem_nhds (show (1 : ℝ) < (2 + I : ℂ).re by norm_num)] with z hz
    exact h1 hz
  -- U4 = {-1 < re < 1}
  set U₄ : Set ℂ := {s : ℂ | -1 < s.re ∧ s.re < 1}
  have hU₄o : IsOpen U₄ := (isOpen_lt continuous_const Complex.continuous_re).inter
    (isOpen_lt Complex.continuous_re continuous_const)
  have hU₄s : U₄ ⊆ {s : ℂ | -1 < s.re ∧ s ≠ 1} := fun z hz =>
    ⟨hz.1, fun h => by have := hz.2; rw [h] at this; simp at this⟩
  have hU₄c : Convex ℝ U₄ := (convex_halfSpace_re_gt (-1)).inter (convex_halfSpace_re_lt 1)
  have h4 : EqOn Φ (fun _ => 0) U₄ := by
    refine (hanal U₄ hU₄o hU₄s).eqOn_of_preconnected_of_eventuallyEq (hzero.mono (subset_univ _))
      hU₄c.isPreconnected (z₀ := 1 / 2 + I) ⟨by simp; norm_num, by simp; norm_num⟩ ?_
    filter_upwards [hU₂o.mem_nhds (show (1 / 2 + I : ℂ) ∈ U₂ from ⟨by simp; norm_num, by simp⟩)] with z hz
    exact h2 hz
  have := h4 ⟨hs1, hs2⟩
  simp only [hΦ] at this
  linear_combination this

end Riemann.Deep.LeanNBRate

