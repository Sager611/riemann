import Mathlib.NumberTheory.LSeries.SumCoeff
import Mathlib.NumberTheory.LSeries.RiemannZeta
import Mathlib.NumberTheory.LSeries.Dirichlet
import Mathlib.Analysis.MellinTransform
import Mathlib.Analysis.Complex.CauchyIntegral
import Mathlib.Analysis.Complex.Convex
import Mathlib.Analysis.SpecialFunctions.Integrals.Basic
import Mathlib.MeasureTheory.Function.Floor

/-!
# NBLp in Lean, part 1: Mellin transforms of the fractional part (copied from NymanBeurling/Suff.lean)

Sections A-F of Riemann/Deep/NymanBeurling/Suff.lean (namespace Fayed.NB) verbatim, and
isPreconnected_U from Riemann/Deep/Mertens/ZeroFree.lean, so that this file needs Mathlib only.
Main facts: zeta_eq (zeta(s) = s/(s-1) - s F(s), Re s > 0, s /= 1), M_eq_zero, hasMellin_rho,
integral_rho_eq_zero (zeta(s) = 0, 0 < Re s < 1 => int_0^1 rho_theta x^(s-1) dx = 0).
-/

open Complex Asymptotics Filter Topology MeasureTheory Set

namespace NBLp

lemma isPreconnected_U (θ : ℝ) : IsPreconnected ({s : ℂ | θ < s.re} \ {1}) := by
  let c : ℂ := ((max θ 1 + 1 : ℝ) : ℂ)
  have hθ1 := le_max_left θ 1
  have h11 := le_max_right θ 1
  have cP1 : IsPreconnected ({s : ℂ | θ < s.re} ∩ {s | 0 < s.im}) :=
    ((convex_halfSpace_re_gt θ).inter (convex_halfSpace_im_gt 0)).isPreconnected
  have cP2 : IsPreconnected ({s : ℂ | θ < s.re} ∩ {s | s.im < 0}) :=
    ((convex_halfSpace_re_gt θ).inter (convex_halfSpace_im_lt 0)).isPreconnected
  have cP3 : IsPreconnected {s : ℂ | max θ 1 < s.re} :=
    (convex_halfSpace_re_gt _).isPreconnected
  have cP4 : IsPreconnected ({s : ℂ | θ < s.re} ∩ {s | s.re < 1}) :=
    ((convex_halfSpace_re_gt θ).inter (convex_halfSpace_re_lt 1)).isPreconnected
  have sP1 : {s : ℂ | θ < s.re} ∩ {s | 0 < s.im} ⊆ {s : ℂ | θ < s.re} \ {1} := by
    rintro s ⟨h1, h2⟩
    refine ⟨h1, fun hs => ?_⟩
    rw [mem_singleton_iff] at hs; subst hs; simp at h2
  have sP2 : {s : ℂ | θ < s.re} ∩ {s | s.im < 0} ⊆ {s : ℂ | θ < s.re} \ {1} := by
    rintro s ⟨h1, h2⟩
    refine ⟨h1, fun hs => ?_⟩
    rw [mem_singleton_iff] at hs; subst hs; simp at h2
  have sP3 : {s : ℂ | max θ 1 < s.re} ⊆ {s : ℂ | θ < s.re} \ {1} := by
    intro s (h1 : max θ 1 < s.re)
    refine ⟨show θ < s.re by linarith, fun hs => ?_⟩
    rw [mem_singleton_iff] at hs; subst hs; simp at h1
  have sP4 : {s : ℂ | θ < s.re} ∩ {s | s.re < 1} ⊆ {s : ℂ | θ < s.re} \ {1} := by
    rintro s ⟨h1, h2⟩
    refine ⟨h1, fun hs => ?_⟩
    rw [mem_singleton_iff] at hs; subst hs; simp at h2
  have hc3 : c ∈ {s : ℂ | max θ 1 < s.re} := by
    show max θ 1 < c.re
    simp only [c, ofReal_re]; linarith
  have hcI3 : c + I ∈ {s : ℂ | max θ 1 < s.re} := by
    show max θ 1 < (c + I).re
    simp only [c, add_re, ofReal_re, I_re]; linarith
  have hcI1 : c + I ∈ {s : ℂ | θ < s.re} ∩ {s | 0 < s.im} := by
    refine ⟨show θ < (c + I).re by simp only [c, add_re, ofReal_re, I_re]; linarith,
      show 0 < (c + I).im by simp [c]⟩
  have hcJ3 : c - I ∈ {s : ℂ | max θ 1 < s.re} := by
    show max θ 1 < (c - I).re
    simp only [c, sub_re, ofReal_re, I_re]; linarith
  have hcJ2 : c - I ∈ {s : ℂ | θ < s.re} ∩ {s | s.im < 0} := by
    refine ⟨show θ < (c - I).re by simp only [c, sub_re, ofReal_re, I_re]; linarith,
      show (c - I).im < 0 by simp [c]⟩
  refine isPreconnected_of_forall c fun y hy => ?_
  obtain ⟨hyθ, hy1⟩ := hy
  have hyθ : θ < y.re := hyθ
  by_cases h3 : max θ 1 < y.re
  · exact ⟨_, sP3, hc3, h3, cP3⟩
  by_cases hpos : 0 < y.im
  · exact ⟨_, union_subset sP1 sP3, Or.inr hc3, Or.inl ⟨hyθ, hpos⟩,
      cP1.union (c + I) hcI1 hcI3 cP3⟩
  by_cases hneg : y.im < 0
  · exact ⟨_, union_subset sP2 sP3, Or.inr hc3, Or.inl ⟨hyθ, hneg⟩,
      cP2.union (c - I) hcJ2 hcJ3 cP3⟩
  have him : y.im = 0 := by linarith [not_lt.mp hpos, not_lt.mp hneg]
  have hre : y.re < 1 := by
    by_contra hge
    push Not at hge h3
    have : y.re = 1 := by
      rcases max_cases θ 1 with ⟨hm, _⟩ | ⟨hm, _⟩ <;> rw [hm] at h3 <;> linarith
    exact hy1 (Complex.ext (by simpa using this) (by simpa using him))
  let w : ℂ := ⟨y.re, 1⟩
  have hw4 : w ∈ {s : ℂ | θ < s.re} ∩ {s | s.re < 1} := ⟨hyθ, hre⟩
  have hw1 : w ∈ {s : ℂ | θ < s.re} ∩ {s | 0 < s.im} := ⟨hyθ, show (0 : ℝ) < 1 from one_pos⟩
  exact ⟨_, union_subset (union_subset sP4 sP1) sP3, Or.inr hc3, Or.inl (Or.inl ⟨hyθ, hre⟩),
    (cP4.union w hw4 hw1 cP1).union (c + I) (Or.inr hcI1) hcI3 cP3⟩

noncomputable def rhoNB (θ x : ℝ) : ℝ := Int.fract (θ / x) - θ * Int.fract (1 / x)


/-! ## A. F(s) = ∫_1^∞ {t} t^{-s-1} dt -/

/-- {t} on [1, ∞), zero below. -/
noncomputable def Q (t : ℝ) : ℂ := if 1 ≤ t then ((Int.fract t : ℝ) : ℂ) else 0

noncomputable def F (s : ℂ) : ℂ := mellin Q (-s)

lemma norm_fract_le (t : ℝ) : ‖((Int.fract t : ℝ) : ℂ)‖ ≤ 1 := by
  rw [Complex.norm_real, Real.norm_of_nonneg (Int.fract_nonneg t)]
  exact (Int.fract_lt_one t).le

lemma norm_Q_le (t : ℝ) : ‖Q t‖ ≤ 1 := by
  unfold Q; split_ifs
  · exact norm_fract_le t
  · simp

lemma Q_eq_zero {t : ℝ} (ht : t < 1) : Q t = 0 := by simp [Q, not_le.mpr ht]

lemma measurable_Q : Measurable Q :=
  Measurable.ite measurableSet_Ici (Complex.measurable_ofReal.comp measurable_fract)
    measurable_const

lemma locallyIntegrableOn_Q : LocallyIntegrableOn Q (Ioi 0) := by
  refine ((continuous_const (y := (1 : ℝ))).locallyIntegrable.locallyIntegrableOn (Ioi 0)).mono
    measurable_Q.aestronglyMeasurable (Eventually.of_forall fun t => ?_)
  simpa using norm_Q_le t

lemma isBigO_Q_atTop : Q =O[atTop] (fun t : ℝ => t ^ (-(0:ℝ))) := by
  refine IsBigO.of_bound 1 (Eventually.of_forall fun t => ?_)
  simpa using norm_Q_le t

lemma isBigO_Q_zero (b : ℝ) : Q =O[𝓝[>] 0] (fun t : ℝ => t ^ (-b)) := by
  have : Q =ᶠ[𝓝[>] (0 : ℝ)] (fun _ => (0 : ℂ)) := by
    filter_upwards [Ioo_mem_nhdsGT (zero_lt_one' ℝ)] with t ht using Q_eq_zero ht.2
  exact (isBigO_zero _ _).congr' this.symm EventuallyEq.rfl

lemma mellinConvergent_Q {s : ℂ} (hs : 0 < s.re) : MellinConvergent Q (-s) :=
  mellinConvergent_of_isBigO_rpow (a := 0) (b := -s.re - 1) locallyIntegrableOn_Q isBigO_Q_atTop
    (by rw [neg_re]; linarith) (isBigO_Q_zero _) (by rw [neg_re]; linarith)

lemma differentiableAt_F {s : ℂ} (hs : 0 < s.re) : DifferentiableAt ℂ F s := by
  have hm : DifferentiableAt ℂ (mellin Q) (-s) :=
    mellin_differentiableAt_of_isBigO_rpow (a := 0) (b := -s.re - 1) locallyIntegrableOn_Q
      isBigO_Q_atTop (by rw [neg_re]; linarith) (isBigO_Q_zero _) (by rw [neg_re]; linarith)
  show DifferentiableAt ℂ (fun z => mellin Q (-z)) s
  exact hm.comp s differentiableAt_id.neg

lemma F_eq_integral (s : ℂ) :
    F s = ∫ t in Ici (1 : ℝ), (t : ℂ) ^ (-s - 1) * ((Int.fract t : ℝ) : ℂ) := by
  rw [F, mellin, setIntegral_eq_of_subset_of_ae_sdiff_eq_zero measurableSet_Ioi.nullMeasurableSet
      (fun t (ht : 1 ≤ t) => (zero_lt_one.trans_le ht : (0 : ℝ) < t))
      (Eventually.of_forall fun t ht => ?_)]
  · refine setIntegral_congr_fun measurableSet_Ici fun t (ht : 1 ≤ t) => ?_
    simp [Q, ht]
  · rw [Q_eq_zero (not_le.mp ht.2), smul_zero]

lemma integrableOn_hi {s : ℂ} (hs : 0 < s.re) :
    IntegrableOn (fun t : ℝ => (t : ℂ) ^ (-s - 1) * ((Int.fract t : ℝ) : ℂ)) (Ici 1) := by
  refine (IntegrableOn.mono_set (mellinConvergent_Q hs)
    (fun t (ht : 1 ≤ t) => (zero_lt_one.trans_le ht : (0 : ℝ) < t))).congr_fun
    (fun t (ht : 1 ≤ t) => ?_) measurableSet_Ici
  simp [Q, ht]

/-! ## B. ζ(s) = s/(s-1) - s F(s) for Re s > 1 -/

lemma zeta_eq_of_one_lt {s : ℂ} (hs : 1 < s.re) :
    riemannZeta s = s / (s - 1) - s * F s := by
  have hO : (fun n : ℕ => ∑ k ∈ Finset.Icc 1 n, ‖(1 : ℕ → ℂ) k‖) =O[atTop]
      (fun n : ℕ => (n : ℝ) ^ (1 : ℝ)) := by
    refine IsBigO.of_bound 1 (Eventually.of_forall fun n => ?_)
    simp
  rw [← LSeries_one_eq_riemannZeta hs, LSeries_eq_mul_integral' _ zero_le_one (by linarith) hO]
  have h1 : ∫ t in Ioi (1:ℝ), (∑ k ∈ Finset.Icc 1 ⌊t⌋₊, (1 : ℕ → ℂ) k) * (t : ℂ) ^ (-(s + 1)) =
      ∫ t in Ioi (1:ℝ), ((t : ℂ) ^ (-s) - (t : ℂ) ^ (-s - 1) * ((Int.fract t : ℝ) : ℂ)) := by
    refine setIntegral_congr_fun measurableSet_Ioi fun t (ht : 1 < t) => ?_
    have h0 : (0:ℝ) ≤ t := zero_le_one.trans ht.le
    have ht0 : (t : ℂ) ≠ 0 := by exact_mod_cast (zero_lt_one.trans ht).ne'
    have hfl : ((⌊t⌋₊ : ℕ) : ℂ) = ((⌊t⌋ : ℤ) : ℂ) := by
      rw [← Int.natCast_floor_eq_floor h0, Int.cast_natCast]
    have hfr : ((Int.fract t : ℝ) : ℂ) = (t : ℂ) - ((⌊t⌋ : ℤ) : ℂ) := by
      rw [Int.fract]; push_cast; ring
    simp only [Pi.one_apply, Finset.sum_const, Nat.card_Icc, nsmul_eq_mul, mul_one,
      Nat.add_sub_cancel]
    rw [hfl, hfr, show -(s + 1) = -s - 1 by ring, cpow_sub _ _ ht0, cpow_one]
    field_simp
    ring
  have hi1 : IntegrableOn (fun t : ℝ => (t : ℂ) ^ (-s)) (Ioi 1) :=
    integrableOn_Ioi_cpow_of_lt (by rw [neg_re]; linarith) zero_lt_one
  have hi2 : IntegrableOn (fun t : ℝ => (t : ℂ) ^ (-s - 1) * ((Int.fract t : ℝ) : ℂ)) (Ioi 1) :=
    (integrableOn_hi (by linarith)).mono_set Ioi_subset_Ici_self
  have h3 : ∫ t in Ioi (1:ℝ), (t : ℂ) ^ (-s - 1) * ((Int.fract t : ℝ) : ℂ) = F s := by
    rw [F_eq_integral, integral_Ici_eq_integral_Ioi]
  have hs1 : s - 1 ≠ 0 := by
    intro h; have := congrArg Complex.re h; simp at this; linarith
  have hs1' : -s + 1 ≠ 0 := by
    intro h; apply hs1; linear_combination -h
  have h4 : ∫ t in Ioi (1:ℝ), (t : ℂ) ^ (-s) = 1 / (s - 1) := by
    rw [integral_Ioi_cpow_of_lt (by rw [neg_re]; linarith) zero_lt_one, ofReal_one, one_cpow]
    field_simp
    ring
  rw [h1, integral_sub hi1 hi2, h3, h4]
  ring

/-! ## C. Analytic continuation to {Re s > 0} \ {1} -/

lemma zeta_eq {s : ℂ} (hs : 0 < s.re) (hs1 : s ≠ 1) :
    riemannZeta s = s / (s - 1) - s * F s := by
  have hUo : IsOpen ({s : ℂ | 0 < s.re} \ {1}) :=
    (isOpen_lt continuous_const Complex.continuous_re).sdiff isClosed_singleton
  have hd1 : DifferentiableOn ℂ riemannZeta ({s : ℂ | 0 < s.re} \ {1}) :=
    fun z hz => (differentiableAt_riemannZeta hz.2).differentiableWithinAt
  have hd2 : DifferentiableOn ℂ (fun z => z / (z - 1) - z * F z) ({s : ℂ | 0 < s.re} \ {1}) := by
    intro z hz
    have hz1 : z - 1 ≠ 0 := sub_ne_zero.mpr hz.2
    have hF := differentiableAt_F hz.1
    have : DifferentiableAt ℂ (fun z => z / (z - 1) - z * F z) z := by
      fun_prop (disch := exact hz1)
    exact this.differentiableWithinAt
  have hc : (2 : ℂ) ∈ {s : ℂ | 0 < s.re} \ {1} :=
    ⟨show (0:ℝ) < (2:ℂ).re by norm_num, fun h => by norm_num at h⟩
  have hev : riemannZeta =ᶠ[𝓝 (2:ℂ)] (fun z => z / (z - 1) - z * F z) := by
    have hmem : {z : ℂ | 1 < z.re} ∈ 𝓝 (2:ℂ) :=
      (isOpen_lt continuous_const Complex.continuous_re).mem_nhds
        (show (1:ℝ) < (2:ℂ).re by norm_num)
    filter_upwards [hmem] with z hz using zeta_eq_of_one_lt hz
  exact (hd1.analyticOnNhd hUo).eqOn_of_preconnected_of_eventuallyEq (hd2.analyticOnNhd hUo)
    (isPreconnected_U 0) hc hev ⟨hs, hs1⟩

/-! ## D. M(s) = ∫_0^∞ {1/x} x^{s-1} dx = 1/(1-s) + F(s) on the strip -/

noncomputable def M (s : ℂ) : ℂ := mellin (fun x : ℝ => ((Int.fract (1 / x) : ℝ) : ℂ)) s

lemma lo_eq (s : ℂ) {t : ℝ} (ht : t ∈ Ioo (0:ℝ) 1) :
    (t : ℂ) ^ (-s - 1) * ((Int.fract t : ℝ) : ℂ) = (t : ℂ) ^ (-s) := by
  have ht0 : (t : ℂ) ≠ 0 := by exact_mod_cast ht.1.ne'
  rw [Int.fract_eq_self.mpr ⟨ht.1.le, ht.2⟩, cpow_sub _ _ ht0, cpow_one]
  field_simp

lemma integrableOn_lo {s : ℂ} (hs : s.re < 1) :
    IntegrableOn (fun t : ℝ => (t : ℂ) ^ (-s - 1) * ((Int.fract t : ℝ) : ℂ)) (Ioo 0 1) := by
  have h : IntegrableOn (fun t : ℝ => (t : ℂ) ^ (-s)) (Ioo 0 1) := by
    have := intervalIntegral.intervalIntegrable_cpow' (a := 0) (b := 1) (r := -s)
      (by rw [neg_re]; linarith)
    rw [intervalIntegrable_iff_integrableOn_Ioc_of_le zero_le_one] at this
    exact this.mono_set Ioo_subset_Ioc_self
  exact h.congr_fun (fun t ht => (lo_eq s ht).symm) measurableSet_Ioo

lemma disj : Disjoint (Ioo (0:ℝ) 1) (Ici 1) :=
  Set.disjoint_left.mpr fun _ ht ht' => absurd ht.2 (not_lt.mpr ht')

lemma mellinConvergent_fract {s : ℂ} (hs0 : 0 < s.re) (hs1 : s.re < 1) :
    MellinConvergent (fun t : ℝ => ((Int.fract t : ℝ) : ℂ)) (-s) := by
  unfold MellinConvergent
  rw [← Ioo_union_Ici_eq_Ioi zero_lt_one]
  simp_rw [smul_eq_mul]
  exact (integrableOn_lo hs1).union (integrableOn_hi hs0)

lemma mellin_fract {s : ℂ} (hs0 : 0 < s.re) (hs1 : s.re < 1) :
    mellin (fun t : ℝ => ((Int.fract t : ℝ) : ℂ)) (-s) = 1 / (1 - s) + F s := by
  have hs : -s + 1 ≠ 0 := by
    intro h; have := congrArg Complex.re h; simp at this; linarith
  rw [mellin, ← Ioo_union_Ici_eq_Ioi zero_lt_one]
  simp_rw [smul_eq_mul]
  rw [setIntegral_union disj measurableSet_Ici (integrableOn_lo hs1) (integrableOn_hi hs0),
    F_eq_integral]
  congr 1
  rw [setIntegral_congr_fun measurableSet_Ioo (fun t ht => lo_eq s ht),
    ← integral_Ioc_eq_integral_Ioo, ← intervalIntegral.integral_of_le zero_le_one,
    integral_cpow (Or.inl (by rw [neg_re]; linarith)), ofReal_one, ofReal_zero, one_cpow,
    zero_cpow hs, sub_zero, neg_add_eq_sub]

lemma M_eq {s : ℂ} (hs0 : 0 < s.re) (hs1 : s.re < 1) : M s = 1 / (1 - s) + F s := by
  rw [M, ← mellin_fract hs0 hs1, ← mellin_comp_inv]
  simp_rw [one_div]

lemma mellinConvergent_M {s : ℂ} (hs0 : 0 < s.re) (hs1 : s.re < 1) :
    MellinConvergent (fun x : ℝ => ((Int.fract (1 / x) : ℝ) : ℂ)) s := by
  have h := (MellinConvergent.comp_rpow (f := fun t : ℝ => ((Int.fract t : ℝ) : ℂ)) (s := s)
    (a := -1) (by norm_num)).mpr (by rw [ofReal_neg, ofReal_one, div_neg, div_one]; exact mellinConvergent_fract hs0 hs1)
  have e : (fun x : ℝ => ((Int.fract (1 / x) : ℝ) : ℂ)) =
      fun t => ((Int.fract (t ^ (-1:ℝ)) : ℝ) : ℂ) := by
    ext x; rw [Real.rpow_neg_one, one_div]
  rw [e]; exact h

lemma M_eq_zero {s : ℂ} (hs0 : 0 < s.re) (hs1 : s.re < 1) (hz : riemannZeta s = 0) : M s = 0 := by
  have hs : s ≠ 0 := by rintro rfl; simp at hs0
  have hne : s ≠ 1 := by rintro rfl; simp at hs1
  have hs1' : s - 1 ≠ 0 := sub_ne_zero.mpr hne
  have h := zeta_eq hs0 hne
  rw [hz] at h
  have h2 : s * (F s - 1 / (s - 1)) = 0 := by
    rw [mul_sub, mul_one_div]; linear_combination h
  have hF : F s = 1 / (s - 1) := sub_eq_zero.mp ((mul_eq_zero.mp h2).resolve_left hs)
  rw [M_eq hs0 hs1, hF, show (1:ℂ) - s = -(s - 1) by ring, one_div_neg_eq_neg_one_div]
  ring

/-! ## E/F. ∫_0^1 rho_θ x^{s-1} = θ^s M(s) - θ M(s) -/

lemma fract_div_eq (θ x : ℝ) : Int.fract (θ / x) = Int.fract (1 / (θ⁻¹ * x)) := by
  rw [one_div, mul_inv, inv_inv, div_eq_mul_inv]

lemma mellin_theta {θ : ℝ} (hθ : 0 < θ) (s : ℂ) :
    mellin (fun x : ℝ => ((Int.fract (θ / x) : ℝ) : ℂ)) s = ((θ⁻¹ : ℝ) : ℂ) ^ (-s) • M s := by
  simp_rw [fract_div_eq θ]
  exact mellin_comp_mul_left (fun x : ℝ => ((Int.fract (1 / x) : ℝ) : ℂ)) s (inv_pos.mpr hθ)

lemma mellinConvergent_theta {θ : ℝ} (hθ : 0 < θ) {s : ℂ} (hs0 : 0 < s.re) (hs1 : s.re < 1) :
    MellinConvergent (fun x : ℝ => ((Int.fract (θ / x) : ℝ) : ℂ)) s := by
  simp_rw [fract_div_eq θ]
  exact (MellinConvergent.comp_mul_left (f := fun x : ℝ => ((Int.fract (1 / x) : ℝ) : ℂ))
    (inv_pos.mpr hθ)).mpr (mellinConvergent_M hs0 hs1)

lemma hasMellin_rho {θ : ℝ} (hθ : 0 < θ) {s : ℂ} (hs0 : 0 < s.re) (hs1 : s.re < 1) :
    HasMellin (fun x : ℝ => ((rhoNB θ x : ℝ) : ℂ)) s
      (((θ⁻¹ : ℝ) : ℂ) ^ (-s) • M s - (θ : ℂ) • M s) := by
  have h1 := hasMellin_const_smul (mellinConvergent_M hs0 hs1) (θ : ℂ)
  have h := hasMellin_sub (mellinConvergent_theta hθ hs0 hs1) h1.1
  rw [h1.2, mellin_theta hθ] at h
  have e : (fun x : ℝ => ((rhoNB θ x : ℝ) : ℂ)) = fun x =>
      ((Int.fract (θ / x) : ℝ) : ℂ) - (θ : ℂ) • ((Int.fract (1 / x) : ℝ) : ℂ) := by
    ext x; simp [rhoNB]
  rw [e]; exact h

lemma rhoNB_eq_zero {θ x : ℝ} (hθ : 0 < θ) (hθ1 : θ ≤ 1) (hx : 1 < x) : rhoNB θ x = 0 := by
  have hx0 : 0 < x := zero_lt_one.trans hx
  unfold rhoNB
  rw [Int.fract_eq_self.mpr ⟨(div_pos hθ hx0).le, (div_lt_one hx0).mpr (hθ1.trans_lt hx)⟩,
    Int.fract_eq_self.mpr ⟨(div_pos one_pos hx0).le, (div_lt_one hx0).mpr hx⟩]
  ring

lemma integral_rho_eq_zero {θ : ℝ} (hθ : 0 < θ) (hθ1 : θ ≤ 1) {s : ℂ} (hs0 : 0 < s.re)
    (hs1 : s.re < 1) (hz : riemannZeta s = 0) :
    ∫ x in Ioc (0:ℝ) 1, ((rhoNB θ x : ℝ) : ℂ) * (x : ℂ) ^ (s - 1) = 0 := by
  have h := hasMellin_rho hθ hs0 hs1
  rw [M_eq_zero hs0 hs1 hz, smul_zero, smul_zero, sub_zero] at h
  rw [← h.2, mellin, setIntegral_eq_of_subset_of_ae_sdiff_eq_zero
    measurableSet_Ioi.nullMeasurableSet Ioc_subset_Ioi_self
    (Eventually.of_forall fun x hx => ?_)]
  · refine setIntegral_congr_fun measurableSet_Ioc fun x _ => ?_
    rw [smul_eq_mul, mul_comm]
  · rw [rhoNB_eq_zero hθ hθ1 (not_le.mp fun h => hx.2 ⟨hx.1, h⟩), ofReal_zero, smul_zero]

lemma integrableOn_rho {θ : ℝ} (hθ : 0 < θ) {s : ℂ} (hs0 : 0 < s.re) (hs1 : s.re < 1) :
    IntegrableOn (fun x : ℝ => ((rhoNB θ x : ℝ) : ℂ) * (x : ℂ) ^ (s - 1)) (Ioc 0 1) := by
  refine (IntegrableOn.mono_set (hasMellin_rho hθ hs0 hs1).1 Ioc_subset_Ioi_self).congr_fun
    (fun x _ => ?_) measurableSet_Ioc
  dsimp only
  rw [smul_eq_mul, mul_comm]

lemma measurable_rhoNB (θ : ℝ) : Measurable (rhoNB θ) :=
  (measurable_fract.comp (measurable_const.div measurable_id)).sub
    (measurable_const.mul (measurable_fract.comp (measurable_const.div measurable_id)))

lemma abs_rhoNB_le {θ : ℝ} (hθ : 0 < θ) (x : ℝ) : |rhoNB θ x| ≤ 1 + θ := by
  unfold rhoNB
  have h1 := Int.fract_nonneg (θ / x)
  have h2 := Int.fract_lt_one (θ / x)
  have h3 := Int.fract_nonneg (1 / x)
  have h4 := Int.fract_lt_one (1 / x)
  rw [abs_le]; constructor <;> nlinarith

end NBLp
