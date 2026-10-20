import Mathlib

/-!
# ζ(s) = s/(s-1) - s ∫₁^∞ {t} t^{-s-1} dt on Re s > 0, Im s > 0

Partial summation (Mathlib `LSeries_eq_mul_integral'`) for Re s > 1, then analytic continuation
(identity theorem on the convex set {Re s > 0, Im s > 0}). Consequence: crude polynomial bound
‖ζ s‖ ≤ ‖s‖/‖s-1‖ + ‖s‖/Re s.
-/

open Complex Filter Topology Set MeasureTheory Asymptotics
open scoped Real ComplexConjugate

namespace Fayed.Convexity

noncomputable section

/-- Fractional part via `Nat.floor`. -/
def fr (t : ℝ) : ℝ := t - (⌊t⌋₊ : ℝ)

lemma fr_nonneg {t : ℝ} (ht : 0 ≤ t) : 0 ≤ fr t := by
  unfold fr; linarith [Nat.floor_le ht]

lemma fr_lt_one (t : ℝ) : fr t < 1 := by
  unfold fr; linarith [Nat.lt_floor_add_one t]

lemma measurable_fr : Measurable fr :=
  measurable_id.sub (measurable_from_nat.comp measurable_id.nat_floor)

/-- `{t}` on `(1,∞)`, zero elsewhere. -/
def F : ℝ → ℂ := (Ioi (1 : ℝ)).indicator (fun t => ((fr t : ℝ) : ℂ))

lemma norm_F_le (t : ℝ) : ‖F t‖ ≤ 1 := by
  unfold F indicator
  split_ifs with h
  · have h0 : 0 ≤ fr t := fr_nonneg (le_of_lt (lt_trans zero_lt_one h))
    rw [Complex.norm_real, Real.norm_of_nonneg h0]; exact (fr_lt_one t).le
  · simp

lemma measurable_F : Measurable F :=
  (Complex.measurable_ofReal.comp measurable_fr).indicator measurableSet_Ioi

/-- I(s) = ∫₁^∞ {t} t^{-s-1} dt as a Mellin transform. -/
def Iz (s : ℂ) : ℂ := mellin F (-s)

lemma mellin_F (z : ℂ) : mellin F z = ∫ t in Ioi (1 : ℝ), (t : ℂ) ^ (z - 1) * ((fr t : ℝ) : ℂ) := by
  rw [mellin]
  have : (fun t : ℝ => (t : ℂ) ^ (z - 1) • F t) =
      (Ioi (1 : ℝ)).indicator (fun t => (t : ℂ) ^ (z - 1) * ((fr t : ℝ) : ℂ)) := by
    ext t; simp only [F, smul_eq_mul, indicator]; split_ifs <;> simp
  rw [this, setIntegral_indicator measurableSet_Ioi,
    inter_eq_right.mpr (Ioi_subset_Ioi zero_le_one)]

lemma Iz_eq (s : ℂ) : Iz s = ∫ t in Ioi (1 : ℝ), (t : ℂ) ^ (-s - 1) * ((fr t : ℝ) : ℂ) :=
  mellin_F (-s)

lemma kernel_norm_le {s : ℂ} {t : ℝ} (ht : 1 < t) :
    ‖(t : ℂ) ^ (-s - 1) * ((fr t : ℝ) : ℂ)‖ ≤ t ^ (-s.re - 1) := by
  have ht0 : 0 < t := lt_trans zero_lt_one ht
  rw [norm_mul, Complex.norm_cpow_eq_rpow_re_of_pos ht0, Complex.norm_real,
    Real.norm_of_nonneg (fr_nonneg ht0.le)]
  have h1 : (-s - 1).re = -s.re - 1 := by simp
  rw [h1]
  exact mul_le_of_le_one_right (Real.rpow_nonneg ht0.le _) (fr_lt_one t).le

lemma kernel_integrable {s : ℂ} (hs : 0 < s.re) :
    IntegrableOn (fun t : ℝ => (t : ℂ) ^ (-s - 1) * ((fr t : ℝ) : ℂ)) (Ioi 1) := by
  refine Integrable.mono' (integrableOn_Ioi_rpow_of_lt (a := -s.re - 1) (by linarith) one_pos)
    ?_ ?_
  · exact ((Complex.measurable_ofReal.pow_const _).mul
      (Complex.measurable_ofReal.comp measurable_fr)).aestronglyMeasurable
  · exact (ae_restrict_iff' measurableSet_Ioi).mpr (Eventually.of_forall fun t ht =>
      kernel_norm_le ht)

lemma norm_Iz_le {s : ℂ} (hs : 0 < s.re) : ‖Iz s‖ ≤ 1 / s.re := by
  rw [Iz_eq]
  refine (norm_integral_le_of_norm_le
    (integrableOn_Ioi_rpow_of_lt (a := -s.re - 1) (by linarith) one_pos)
    ((ae_restrict_iff' measurableSet_Ioi).mpr (Eventually.of_forall fun t ht =>
      kernel_norm_le ht))).trans (le_of_eq ?_)
  rw [integral_Ioi_rpow_of_lt (by linarith) one_pos, Real.one_rpow]
  have : -s.re - 1 + 1 = -s.re := by ring
  rw [this, neg_div_neg_eq]

lemma locallyIntegrable_F : LocallyIntegrable F := by
  intro x
  refine ⟨Ioo (x - 1) (x + 1), Ioo_mem_nhds (by linarith) (by linarith), ?_⟩
  exact Measure.integrableOn_of_bounded (M := 1) measure_Ioo_lt_top.ne
    measurable_F.aestronglyMeasurable (ae_of_all _ norm_F_le)

lemma differentiableAt_Iz {s : ℂ} (hs : 0 < s.re) : DifferentiableAt ℂ Iz s := by
  have hm : DifferentiableAt ℂ (mellin F) (-s) := by
    refine mellin_differentiableAt_of_isBigO_rpow (a := 0) (b := (-s).re - 1)
      (locallyIntegrable_F.locallyIntegrableOn _) ?_ (by simp; linarith) ?_ (by linarith)
    · refine IsBigO.of_bound 1 ?_
      filter_upwards [eventually_gt_atTop (0 : ℝ)] with t ht
      rw [neg_zero, Real.rpow_zero, norm_one, one_mul]
      exact norm_F_le t
    · have hz : F =ᶠ[𝓝[>] 0] fun _ => (0 : ℂ) := by
        filter_upwards [nhdsWithin_le_nhds (Iio_mem_nhds (zero_lt_one : (0 : ℝ) < 1))] with t ht
        simp only [F, indicator, mem_Ioi]
        rw [ite_eq_right_iff.mpr (fun h => absurd (lt_trans h ht) (lt_irrefl 1))]
      exact (isBigO_zero _ _).congr' hz.symm EventuallyEq.rfl
  exact hm.comp s differentiableAt_id.neg

/-- Partial summation for Re s > 1. -/
lemma zeta_eq_of_one_lt {s : ℂ} (hs : 1 < s.re) :
    riemannZeta s = s / (s - 1) - s * Iz s := by
  have hO : (fun n : ℕ => ∑ k ∈ Finset.Icc 1 n, ‖(1 : ℕ → ℂ) k‖) =O[atTop]
      fun n : ℕ => (n : ℝ) ^ (1 : ℝ) := by
    simp only [Pi.one_apply, norm_one, Finset.sum_const, Nat.card_Icc, nsmul_eq_mul, mul_one,
      Real.rpow_one, add_tsub_cancel_right]
    exact isBigO_refl _ _
  have h1 := LSeries_eq_mul_integral' (1 : ℕ → ℂ) zero_le_one (by simpa using hs) hO
  rw [LSeries_one_eq_riemannZeta hs] at h1
  rw [h1]
  have hs1 : -s + 1 ≠ 0 := by
    intro h; have := congrArg Complex.re h; simp at this; linarith
  have heq : EqOn (fun t : ℝ => (∑ k ∈ Finset.Icc 1 ⌊t⌋₊, (1 : ℕ → ℂ) k) * (t : ℂ) ^ (-(s + 1)))
      (fun t : ℝ => (t : ℂ) ^ (-s) - (t : ℂ) ^ (-s - 1) * ((fr t : ℝ) : ℂ)) (Ioi 1) := by
    intro t ht
    have ht0 : (t : ℂ) ≠ 0 := by exact_mod_cast (lt_trans zero_lt_one ht).ne'
    simp only [Pi.one_apply, Finset.sum_const, Nat.card_Icc, add_tsub_cancel_right, nsmul_eq_mul,
      mul_one]
    have e : (t : ℂ) ^ (-s) = (t : ℂ) * (t : ℂ) ^ (-s - 1) := by
      rw [show -s = 1 + (-s - 1) by ring, cpow_add _ _ ht0, cpow_one]
      ring_nf
    rw [e, fr, show -(s + 1) = -s - 1 by ring]
    push_cast; ring
  rw [setIntegral_congr_fun measurableSet_Ioi heq,
    integral_sub (integrableOn_Ioi_cpow_of_lt (by simp; linarith) zero_lt_one)
      (kernel_integrable (by linarith)),
    integral_Ioi_cpow_of_lt (by simp; linarith) zero_lt_one, ← Iz_eq]
  simp only [ofReal_one, one_cpow]
  have hs2 : s - 1 ≠ 0 := by
    intro h; apply hs1; linear_combination -h
  field_simp
  ring

/-- Upper-right open quadrant of the right half-plane. -/
def U : Set ℂ := {s | 0 < s.re} ∩ {s | 0 < s.im}

lemma isOpen_U : IsOpen U :=
  (isOpen_lt continuous_const Complex.continuous_re).inter
    (isOpen_lt continuous_const Complex.continuous_im)

/-- Analytic continuation of the partial-summation formula to `U`. -/
lemma zeta_eq_on_U {s : ℂ} (hs : s ∈ U) : riemannZeta s = s / (s - 1) - s * Iz s := by
  have hne : ∀ z ∈ U, z ≠ 1 := by
    intro z hz h; have := hz.2; rw [h] at this; simp at this
  have hf : AnalyticOnNhd ℂ riemannZeta U := by
    refine DifferentiableOn.analyticOnNhd (fun z hz => ?_) isOpen_U
    exact (differentiableAt_riemannZeta (hne z hz)).differentiableWithinAt
  have hg : AnalyticOnNhd ℂ (fun s => s / (s - 1) - s * Iz s) U := by
    refine DifferentiableOn.analyticOnNhd (fun z hz => ?_) isOpen_U
    have h1 : z - 1 ≠ 0 := sub_ne_zero.mpr (hne z hz)
    have e1 : DifferentiableAt ℂ (fun s : ℂ => s / (s - 1)) z := by
      fun_prop (disch := exact h1)
    have e2 : DifferentiableAt ℂ (fun s : ℂ => s * Iz s) z :=
      DifferentiableAt.mul differentiableAt_id (differentiableAt_Iz hz.1)
    exact (e1.sub e2).differentiableWithinAt
  have hU : IsPreconnected U :=
    ((convex_halfSpace_re_gt 0).inter (convex_halfSpace_im_gt 0)).isPreconnected
  have h0 : (2 + I : ℂ) ∈ U := ⟨by simp, by simp⟩
  have hfg : riemannZeta =ᶠ[𝓝 (2 + I)] fun s => s / (s - 1) - s * Iz s := by
    have : {z : ℂ | 1 < z.re} ∈ 𝓝 (2 + I) :=
      (isOpen_lt continuous_const Complex.continuous_re).mem_nhds (by norm_num)
    filter_upwards [this] with z hz using zeta_eq_of_one_lt hz
  exact hf.eqOn_of_preconnected_of_eventuallyEq hg hU h0 hfg hs

/-- Crude bound in the quadrant `U`. -/
lemma norm_zeta_le_on_U {s : ℂ} (hs : s ∈ U) :
    ‖riemannZeta s‖ ≤ ‖s‖ / ‖s - 1‖ + ‖s‖ * (1 / s.re) := by
  rw [zeta_eq_on_U hs]
  refine (norm_sub_le _ _).trans ?_
  rw [norm_div, norm_mul]
  gcongr
  exact norm_Iz_le hs.1

end

end Fayed.Convexity
