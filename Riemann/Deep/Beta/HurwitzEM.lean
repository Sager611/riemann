import Riemann.Deep.FirstZero.EM
import Mathlib.NumberTheory.LSeries.HurwitzZeta
import Mathlib.Analysis.MellinTransform

/-!
# Order-8 Euler–Maclaurin for the Hurwitz zeta function

For 0 < a, J : ℕ, c = J + a, 0 < Re s, s ≠ 1:
ζ(s, a) = Σ_{j<J} (j+a)^{-s} - c^{1-s}/(1-s) + c^{-s}/2 + s·∫_c^∞ saw_a(y) y^{-s-1} dy,
saw_a(y) = ⌊y - a⌋ + a + 1/2 - y (first-order formula: Re s > 1 by telescoping, then the identity
theorem; the integral is holomorphic as a Mellin transform). The integral is then expanded with the
order-8 machinery of Riemann/Deep/FirstZero/EM.lean on the intervals [c + i, c + i + 1].
-/

open Complex Set Filter Topology MeasureTheory Interval

namespace Fayed.Beta

open Fayed.FirstZero HurwitzZeta

/-- Sawtooth with jumps at a + ℤ. -/
noncomputable def saw (a y : ℝ) : ℝ := (⌊y - a⌋ : ℝ) + a + 1 / 2 - y

lemma abs_saw_le (a y : ℝ) : |saw a y| ≤ 1 / 2 := by
  unfold saw
  have h1 := Int.floor_le (y - a)
  have h2 := Int.lt_floor_add_one (y - a)
  rw [abs_le]; constructor <;> linarith

lemma measurable_saw (a : ℝ) : Measurable (saw a) := by
  have h : Measurable (fun y : ℝ => ((⌊y - a⌋ : ℤ) : ℝ)) :=
    (measurable_from_top (f := fun z : ℤ => (z : ℝ))).comp
      (Int.measurable_floor.comp (measurable_id.sub_const a))
  exact ((h.add_const a).add_const _).sub measurable_id

lemma integrableOn_saw (a : ℝ) {c : ℝ} (hc : 0 < c) {s : ℂ} (hs : 0 < s.re) :
    IntegrableOn (fun y : ℝ => (saw a y : ℂ) * (y : ℂ) ^ (-(s + 1))) (Ioi c) := by
  refine Integrable.bdd_mul (c := 1 / 2) ?_ ?_ ?_
  · exact (integrableOn_Ioi_cpow_iff hc).mpr (by simp only [neg_re, add_re, one_re]; linarith)
  · exact (Complex.measurable_ofReal.comp (measurable_saw a)).aestronglyMeasurable
  · exact ae_of_all _ fun y => by rw [Complex.norm_real, Real.norm_eq_abs]; exact abs_saw_le a y

lemma hasDerivAt_ofReal_cpow {y : ℝ} (hy : 0 < y) (c : ℂ) :
    HasDerivAt (fun x : ℝ => (x : ℂ) ^ c) (c * (y : ℂ) ^ (c - 1)) y :=
  ((Complex.hasStrictDerivAt_cpow_const (c := c) (ofReal_mem_slitPlane.mpr hy)).hasDerivAt).comp_ofReal

lemma ofReal_cpow_eq_mul {y : ℝ} (hy : 0 < y) (c : ℂ) :
    (y : ℂ) ^ c = (y : ℂ) * (y : ℂ) ^ (c - 1) := by
  have hy0 : (y : ℂ) ≠ 0 := by exact_mod_cast hy.ne'
  conv_lhs => rw [show c = 1 + (c - 1) by ring]
  rw [cpow_add _ _ hy0, cpow_one]

/-- One unit interval of the first-order Euler–Maclaurin formula. -/
lemma piece {s : ℂ} (hs0 : s ≠ 0) (hs1 : s ≠ 1) {n : ℝ} (hn : 0 < n) :
    ∫ y in n..n + 1, ((n + 1 / 2 - y : ℝ) : ℂ) * (y : ℂ) ^ (-(s + 1)) =
      ((n : ℂ) ^ (-s) / 2 + ((n + 1 : ℝ) : ℂ) ^ (-s) / 2
        - ((((n + 1 : ℝ) : ℂ) ^ (1 - s) - (n : ℂ) ^ (1 - s)) / (1 - s))) / s := by
  have hpos : ∀ y ∈ [[n, n + 1]], 0 < y := by
    intro y hy; rw [uIcc_of_le (by linarith)] at hy; linarith [hy.1]
  have h1s : (1 : ℂ) - s ≠ 0 := sub_ne_zero.mpr (Ne.symm hs1)
  have hderiv : ∀ y ∈ [[n, n + 1]], HasDerivAt
      (fun y : ℝ => ((n + 1 / 2 : ℝ) : ℂ) * (y : ℂ) ^ (-s) / (-s) - (y : ℂ) ^ (1 - s) / (1 - s))
      (((n + 1 / 2 - y : ℝ) : ℂ) * (y : ℂ) ^ (-(s + 1))) y := by
    intro y hy
    have hy0 := hpos y hy
    have d1 := ((hasDerivAt_ofReal_cpow hy0 (-s)).const_mul ((n + 1 / 2 : ℝ) : ℂ)).div_const (-s)
    have d2 := (hasDerivAt_ofReal_cpow hy0 (1 - s)).div_const (1 - s)
    convert d1.sub d2 using 1
    rw [show (1 : ℂ) - s - 1 = -s by ring, show -s - 1 = -(s + 1) by ring,
      ofReal_cpow_eq_mul hy0 (-s), show -s - 1 = -(s + 1) by ring]
    have hs0' : -s ≠ 0 := neg_ne_zero.mpr hs0
    field_simp
    push_cast
    ring
  have hcont : ContinuousOn (fun y : ℝ => ((n + 1 / 2 - y : ℝ) : ℂ) * (y : ℂ) ^ (-(s + 1)))
      [[n, n + 1]] := by
    intro y hy
    exact ((continuous_ofReal.comp (continuous_const.sub continuous_id)).continuousAt.mul
      ((hasDerivAt_ofReal_cpow (hpos y hy) _).continuousAt)).continuousWithinAt
  rw [intervalIntegral.integral_eq_sub_of_hasDerivAt hderiv hcont.intervalIntegrable]
  have hn' : (0 : ℝ) < n + 1 := by linarith
  rw [ofReal_cpow_eq_mul hn (1 - s), ofReal_cpow_eq_mul hn' (1 - s), show (1 : ℂ) - s - 1 = -s by ring]
  field_simp
  push_cast
  ring

/-- On [c + i, c + i + 1) with c = J + a the sawtooth is c + i + 1/2 - y. -/
lemma saw_piece (a : ℝ) (J i : ℕ) (s : ℂ) :
    ∫ y in ((J : ℝ) + a + i)..((J : ℝ) + a + i) + 1, (saw a y : ℂ) * (y : ℂ) ^ (-(s + 1)) =
      ∫ y in ((J : ℝ) + a + i)..((J : ℝ) + a + i) + 1,
        ((((J : ℝ) + a + i) + 1 / 2 - y : ℝ) : ℂ) * (y : ℂ) ^ (-(s + 1)) := by
  refine intervalIntegral.integral_congr_ae ?_
  have hae : ∀ᵐ y : ℝ, y ≠ ((J : ℝ) + a + i) + 1 := by
    rw [ae_iff]; simp [measure_singleton]
  filter_upwards [hae] with y hy hmem
  rw [uIoc_of_le (by linarith)] at hmem
  have hfl : ⌊y - a⌋ = ((J + i : ℕ) : ℤ) := by
    rw [Int.floor_eq_iff]; push_cast
    constructor
    · linarith [hmem.1]
    · have := lt_of_le_of_ne hmem.2 hy; linarith
  simp only [saw, hfl]
  push_cast
  ring_nf

lemma intervalIntegrable_saw (a : ℝ) {c : ℝ} (hc : 0 < c) {s : ℂ} (hs : 0 < s.re) {u v : ℝ}
    (hu : c ≤ u) (huv : u ≤ v) :
    IntervalIntegrable (fun y : ℝ => (saw a y : ℂ) * (y : ℂ) ^ (-(s + 1))) volume u v := by
  rw [intervalIntegrable_iff_integrableOn_Ioc_of_le huv]
  exact (integrableOn_saw a hc hs).mono_set (Ioc_subset_Ioi_self.trans (Ioi_subset_Ioi hu))

/-- Finite first-order Euler–Maclaurin over [c, c + K], c = J + a. -/
lemma partial_sum {a : ℝ} {J : ℕ} (hc : 0 < (J : ℝ) + a) {s : ℂ} (hs : 0 < s.re) (hs1 : s ≠ 1)
    (K : ℕ) :
    ∑ i ∈ Finset.range K, ((((J : ℝ) + a + i) : ℝ) : ℂ) ^ (-s) =
      (((((J : ℝ) + a + K) : ℝ) : ℂ) ^ (1 - s) - (((J : ℝ) + a : ℝ) : ℂ) ^ (1 - s)) / (1 - s)
        + (((J : ℝ) + a : ℝ) : ℂ) ^ (-s) / 2 - (((((J : ℝ) + a + K) : ℝ) : ℂ) ^ (-s)) / 2
        + s * ∫ y in ((J : ℝ) + a)..((J : ℝ) + a + K), (saw a y : ℂ) * (y : ℂ) ^ (-(s + 1)) := by
  have hs0 : s ≠ 0 := by rintro rfl; simp at hs
  have h1s : (1 : ℂ) - s ≠ 0 := sub_ne_zero.mpr (Ne.symm hs1)
  induction K with
  | zero => simp
  | succ K ih =>
    have hle : (J : ℝ) + a ≤ (J : ℝ) + a + K := by
      have : (0 : ℝ) ≤ K := Nat.cast_nonneg K; linarith
    have hK : ((J : ℝ) + a + ((K + 1 : ℕ) : ℝ)) = ((J : ℝ) + a + K) + 1 := by push_cast; ring
    rw [Finset.sum_range_succ, ih, hK, ← intervalIntegral.integral_add_adjacent_intervals
      (intervalIntegrable_saw a hc hs le_rfl hle)
      (intervalIntegrable_saw a hc hs hle (v := (J : ℝ) + a + K + 1) (by linarith)), saw_piece a J K s,
      piece hs0 hs1 (by linarith)]
    field_simp
    ring

lemma tendsto_cpow_shift {c : ℝ} (hc : 0 < c) {w : ℂ} (hw : w.re < 0) :
    Tendsto (fun K : ℕ => (((c + K : ℝ)) : ℂ) ^ w) atTop (𝓝 0) := by
  rw [tendsto_zero_iff_norm_tendsto_zero]
  have hb : Tendsto (fun K : ℕ => c + (K : ℝ)) atTop atTop :=
    tendsto_atTop_add_const_left _ _ tendsto_natCast_atTop_atTop
  have h1 : Tendsto (fun y : ℝ => y ^ (-(-w.re))) atTop (𝓝 0) :=
    tendsto_rpow_neg_atTop (by linarith)
  refine (h1.comp hb).congr' ?_
  filter_upwards with K
  have hp : 0 < c + (K : ℝ) := by have : (0 : ℝ) ≤ K := Nat.cast_nonneg K; linarith
  simp only [Function.comp, norm_cpow_eq_rpow_re_of_pos hp, neg_neg]

/-- First-order formula for Re s > 1. -/
lemma hurwitz_first_order_of_one_lt {a : ℝ} (ha : a ∈ Icc (0 : ℝ) 1) {J : ℕ}
    (hc : 0 < (J : ℝ) + a) {s : ℂ} (hs : 1 < s.re) :
    hurwitzZeta a s = ∑ j ∈ Finset.range J, (((j : ℝ) + a : ℝ) : ℂ) ^ (-s)
      + (-(((J : ℝ) + a : ℝ) : ℂ) ^ (1 - s) / (1 - s) + (((J : ℝ) + a : ℝ) : ℂ) ^ (-s) / 2
        + s * ∫ y in Ioi ((J : ℝ) + a), (saw a y : ℂ) * (y : ℂ) ^ (-(s + 1))) := by
  have hs0 : 0 < s.re := by linarith
  have hs1 : s ≠ 1 := by rintro rfl; simp at hs
  have hH := hasSum_hurwitzZeta_of_one_lt_re ha hs
  have hterm : ∀ n : ℕ, 1 / ((n : ℂ) + a) ^ s = (((n : ℝ) + a : ℝ) : ℂ) ^ (-s) := by
    intro n; rw [cpow_neg, one_div]; push_cast; rfl
  simp only [hterm] at hH
  have hT := ((hasSum_nat_add_iff' J).mpr hH).tendsto_sum_nat
  have hlim : Tendsto (fun K : ℕ => ∑ i ∈ Finset.range K, ((((J : ℝ) + a + i) : ℝ) : ℂ) ^ (-s))
      atTop (𝓝 (-(((J : ℝ) + a : ℝ) : ℂ) ^ (1 - s) / (1 - s) + (((J : ℝ) + a : ℝ) : ℂ) ^ (-s) / 2
        + s * ∫ y in Ioi ((J : ℝ) + a), (saw a y : ℂ) * (y : ℂ) ^ (-(s + 1)))) := by
    have hb : Tendsto (fun K : ℕ => (J : ℝ) + a + K) atTop atTop :=
      tendsto_atTop_add_const_left _ _ tendsto_natCast_atTop_atTop
    have h1 := tendsto_cpow_shift hc (w := 1 - s) (by simp only [sub_re, one_re]; linarith)
    have h2 := tendsto_cpow_shift hc (w := -s) (by simp only [neg_re]; linarith)
    have h3 := intervalIntegral_tendsto_integral_Ioi ((J : ℝ) + a) (integrableOn_saw a hc hs0) hb
    have := ((((h1.sub_const ((((J : ℝ) + a : ℝ) : ℂ) ^ (1 - s))).div_const (1 - s)).add_const
      ((((J : ℝ) + a : ℝ) : ℂ) ^ (-s) / 2)).sub (h2.div_const 2)).add (h3.const_mul s)
    have e : (-(((J : ℝ) + a : ℝ) : ℂ) ^ (1 - s) / (1 - s) + (((J : ℝ) + a : ℝ) : ℂ) ^ (-s) / 2
        + s * ∫ y in Ioi ((J : ℝ) + a), (saw a y : ℂ) * (y : ℂ) ^ (-(s + 1))) =
        (0 - (((J : ℝ) + a : ℝ) : ℂ) ^ (1 - s)) / (1 - s) + (((J : ℝ) + a : ℝ) : ℂ) ^ (-s) / 2
          - 0 / 2 + s * ∫ y in Ioi ((J : ℝ) + a), (saw a y : ℂ) * (y : ℂ) ^ (-(s + 1)) := by ring
    rw [e]
    exact this.congr' (Eventually.of_forall fun K => (partial_sum hc hs0 hs1 K).symm)
  have heq := tendsto_nhds_unique hT (hlim.congr' (Eventually.of_forall fun K => by
    refine Finset.sum_congr rfl fun i _ => ?_
    push_cast; ring_nf))
  rw [← heq]
  ring

/-- The sawtooth tail, cut off below c, as a function for the Mellin transform. -/
noncomputable def sawInd (a c : ℝ) : ℝ → ℂ := (Ioi c).indicator (fun y => (saw a y : ℂ))

lemma sawIntegral_eq_mellin (a : ℝ) {c : ℝ} (hc : 0 < c) (s : ℂ) :
    ∫ y in Ioi c, (saw a y : ℂ) * (y : ℂ) ^ (-(s + 1)) = mellin (sawInd a c) (-s) := by
  rw [mellin]
  have h : ∀ t : ℝ, (t : ℂ) ^ (-s - 1) • sawInd a c t =
      (Ioi c).indicator (fun y => (saw a y : ℂ) * (y : ℂ) ^ (-(s + 1))) t := by
    intro t
    by_cases ht : t ∈ Ioi c
    · rw [sawInd, indicator_of_mem ht, indicator_of_mem ht, smul_eq_mul, mul_comm,
        show -s - 1 = -(s + 1) by ring]
    · simp [sawInd, ht]
  simp_rw [h]
  rw [setIntegral_indicator measurableSet_Ioi, Ioi_inter_Ioi, sup_of_le_right hc.le]

lemma norm_sawInd_le (a c y : ℝ) : ‖sawInd a c y‖ ≤ 1 / 2 :=
  le_trans (norm_indicator_le_norm_self _ _)
    (by rw [Complex.norm_real, Real.norm_eq_abs]; exact abs_saw_le a y)

lemma differentiableAt_sawIntegral (a : ℝ) {c : ℝ} (hc : 0 < c) {s : ℂ} (hs : 0 < s.re) :
    DifferentiableAt ℂ (fun z : ℂ => ∫ y in Ioi c, (saw a y : ℂ) * (y : ℂ) ^ (-(z + 1))) s := by
  simp_rw [sawIntegral_eq_mellin a hc]
  have hm : DifferentiableAt ℂ (mellin (sawInd a c)) (-s) := by
    refine mellin_differentiableAt_of_isBigO_rpow (a := 0) (b := (-s).re - 1) ?_ ?_ ?_ ?_ ?_
    · refine (locallyIntegrable_iff.mpr fun k hk => ?_).locallyIntegrableOn _
      refine Measure.integrableOn_of_bounded (M := 1 / 2) hk.measure_lt_top.ne ?_ ?_
      · exact ((Complex.measurable_ofReal.comp (measurable_saw a)).indicator
          measurableSet_Ioi).aestronglyMeasurable
      · exact ae_of_all _ fun y => norm_sawInd_le a c y
    · refine Asymptotics.IsBigO.of_bound (1 / 2) (Eventually.of_forall fun y => ?_)
      rw [neg_zero, Real.rpow_zero, norm_one, mul_one]
      exact norm_sawInd_le a c y
    · simp only [neg_re]; linarith
    · refine Asymptotics.IsBigO.of_bound 0 ?_
      filter_upwards [Ioo_mem_nhdsGT hc] with y hy
      have : sawInd a c y = 0 :=
        indicator_of_notMem (by simp only [mem_Ioi, not_lt]; exact hy.2.le) _
      rw [this, norm_zero, zero_mul]
    · linarith
  exact hm.comp s differentiableAt_id.neg

/-- The first-order Euler–Maclaurin expression for ζ(s, a). -/
noncomputable def hurwitzFO (a : ℝ) (J : ℕ) (s : ℂ) : ℂ :=
  ∑ j ∈ Finset.range J, (((j : ℝ) + a : ℝ) : ℂ) ^ (-s)
    + (-(((J : ℝ) + a : ℝ) : ℂ) ^ (1 - s) / (1 - s) + (((J : ℝ) + a : ℝ) : ℂ) ^ (-s) / 2
      + s * ∫ y in Ioi ((J : ℝ) + a), (saw a y : ℂ) * (y : ℂ) ^ (-(s + 1)))

lemma differentiableOn_hurwitzFO {a : ℝ} (ha : 0 < a) (J : ℕ) :
    DifferentiableOn ℂ (hurwitzFO a J) {s : ℂ | s ≠ 1 ∧ 0 < s.re} := by
  intro s ⟨hs1, hs⟩
  have hc : 0 < (J : ℝ) + a := by positivity
  refine DifferentiableAt.differentiableWithinAt ?_
  have hpow : ∀ x : ℝ, 0 < x → ∀ f : ℂ → ℂ, DifferentiableAt ℂ f s →
      DifferentiableAt ℂ (fun z => ((x : ℝ) : ℂ) ^ f z) s :=
    fun x hx f hf => hf.const_cpow (Or.inl (by exact_mod_cast hx.ne'))
  have h1 : DifferentiableAt ℂ (fun z : ℂ => ∑ j ∈ Finset.range J, (((j : ℝ) + a : ℝ) : ℂ) ^ (-z)) s :=
    DifferentiableAt.fun_sum fun j _ => hpow _ (by positivity) _ differentiableAt_id.neg
  have h2 : DifferentiableAt ℂ (fun z : ℂ => -(((J : ℝ) + a : ℝ) : ℂ) ^ (1 - z) / (1 - z)) s :=
    (hpow _ hc _ ((differentiableAt_const _).sub differentiableAt_id)).neg.div
      ((differentiableAt_const _).sub differentiableAt_id) (sub_ne_zero.mpr (Ne.symm hs1))
  have h3 : DifferentiableAt ℂ (fun z : ℂ => (((J : ℝ) + a : ℝ) : ℂ) ^ (-z) / 2) s :=
    (hpow _ hc _ differentiableAt_id.neg).div_const 2
  have h4 := differentiableAt_id.mul (differentiableAt_sawIntegral a hc hs)
  exact h1.add ((h2.add h3).add h4)

/-- **First-order formula for the Hurwitz zeta function on Re s > 0.** -/
theorem hurwitz_first_order {a : ℝ} (ha : a ∈ Ioc (0 : ℝ) 1) (J : ℕ) {s : ℂ} (hs : 0 < s.re)
    (hs1 : s ≠ 1) : hurwitzZeta a s = hurwitzFO a J s := by
  have hc : 0 < (J : ℝ) + a := by have := ha.1; positivity
  set U := {z : ℂ | z ≠ 1 ∧ 0 < z.re} with hU
  have f_an : AnalyticOnNhd ℂ (hurwitzZeta (a : UnitAddCircle)) U :=
    DifferentiableOn.analyticOnNhd
      (fun z hz => (differentiableAt_hurwitzZeta _ hz.1).differentiableWithinAt) PNT.isOpen_aux
  have g_an : AnalyticOnNhd ℂ (hurwitzFO a J) U :=
    (differentiableOn_hurwitzFO ha.1 J).analyticOnNhd PNT.isOpen_aux
  have preconU : IsPreconnected U :=
    ((IsOpen.isConnected_iff_isPathConnected PNT.isOpen_aux).mpr PNT.isPathConnected_aux).isPreconnected
  have h2 : (2 : ℂ) ∈ U := ⟨by norm_num, by simp⟩
  refine AnalyticOnNhd.eqOn_of_preconnected_of_eventuallyEq f_an g_an preconU h2 ?_ ⟨hs1, hs⟩
  have u_mem : {z : ℂ | 1 < z.re} ∈ 𝓝 (2 : ℂ) :=
    (isOpen_lt continuous_const continuous_re).mem_nhds (by simp)
  filter_upwards [u_mem] with z hz
  exact hurwitz_first_order_of_one_lt (Ioc_subset_Icc_self ha) hc hz

lemma saw_piece_BP (a : ℝ) (J i : ℕ) (s : ℂ) :
    ∫ y in ((J : ℝ) + a + i)..((J : ℝ) + a + i) + 1, (saw a y : ℂ) * (y : ℂ) ^ (-(s + 1)) =
      -∫ y in ((J : ℝ) + a + i)..((J : ℝ) + a + i) + 1,
        (BP 1 (y - ((J : ℝ) + a + i)) : ℂ) * gf s 0 y := by
  rw [saw_piece, ← intervalIntegral.integral_neg]
  congr 1; ext y
  have hexp : (y : ℂ) ^ (-s - 1 - ((0 : ℕ) : ℂ)) = (y : ℂ) ^ (-(s + 1)) := by
    congr 1; push_cast; ring
  simp only [gf, hexp, gc, BP]
  push_cast
  ring

/-- Order-8 Euler–Maclaurin for the shifted sawtooth integral. -/
theorem saw_integral_EM8 {a : ℝ} {J : ℕ} (hc : 0 < (J : ℝ) + a) {s : ℂ} (hs : 0 < s.re) :
    ‖(∫ y in Ioi ((J : ℝ) + a), (saw a y : ℂ) * (y : ℂ) ^ (-(s + 1))) - EMPhi s ((J : ℝ) + a)‖ ≤
      ‖gc s 7‖ / 1000000 * (((J : ℝ) + a) ^ (-(s.re + 7)) / (s.re + 7)) := by
  have hpiece := saw_piece_BP a J
  set c := (J : ℝ) + a with hcdef
  set f : ℝ → ℂ := fun y => (saw a y : ℂ) * (y : ℂ) ^ (-(s + 1)) with hf
  have hfi : IntegrableOn f (Ioi c) := integrableOn_saw a hc hs
  set R : ℕ → ℂ := fun i => ∫ x in (c + i)..(c + i) + 1, (BP 8 (x - (c + i)) : ℂ) * gf s 7 x with hR
  have hii : ∀ u v : ℝ, c ≤ u → u ≤ v → IntervalIntegrable f volume u v :=
    fun u v hu huv => intervalIntegrable_saw a hc hs hu huv
  have key : ∀ m : ℕ, ∫ x in c..(c + m), f x =
      EMPhi s c - EMPhi s (c + m) + ∑ i ∈ Finset.range m, R i := by
    intro m
    induction m with
    | zero => simp
    | succ m ih =>
      have hcm : c + ((m + 1 : ℕ) : ℝ) = (c + m) + 1 := by push_cast; ring
      have hle : c ≤ c + m := by have : (0 : ℝ) ≤ m := Nat.cast_nonneg m; linarith
      rw [hcm, ← intervalIntegral.integral_add_adjacent_intervals (hii _ _ le_rfl hle)
        (hii _ _ hle (by linarith)), ih, Finset.sum_range_succ]
      have hp : 0 < c + m := by linarith
      have h1 := hpiece m s
      have h2 := interval_EM8 s hp
      simp only [hf] at h1 ⊢
      rw [h1, h2]
      simp only [hR]
      ring
  have hbound : ∀ m : ℕ, ‖∑ i ∈ Finset.range m, R i‖ ≤
      ‖gc s 7‖ / 1000000 * (c ^ (-(s.re + 7)) / (s.re + 7)) := by
    intro m
    set h : ℕ → ℝ := fun i => (c + i) ^ (-(s.re + 7)) / (s.re + 7) with hh
    have hpiece : ∀ i, ‖R i‖ ≤ ‖gc s 7‖ / 1000000 * (h i - h (i + 1)) := by
      intro i
      have hp : 0 < c + i := by have : (0 : ℝ) ≤ i := Nat.cast_nonneg i; linarith
      have := remainder_piece s hs hp
      have hc1 : c + ((i + 1 : ℕ) : ℝ) = (c + i) + 1 := by push_cast; ring
      simp only [hh, hc1]
      rw [← sub_div]
      exact this
    calc ‖∑ i ∈ Finset.range m, R i‖ ≤ ∑ i ∈ Finset.range m, ‖R i‖ := norm_sum_le _ _
      _ ≤ ∑ i ∈ Finset.range m, ‖gc s 7‖ / 1000000 * (h i - h (i + 1)) :=
          Finset.sum_le_sum fun i _ => hpiece i
      _ = ‖gc s 7‖ / 1000000 * (h 0 - h m) := by
          rw [← Finset.mul_sum, Finset.sum_range_sub']
      _ ≤ ‖gc s 7‖ / 1000000 * h 0 := by
          apply mul_le_mul_of_nonneg_left _ (by positivity)
          have : 0 ≤ h m := by
            simp only [hh]
            apply div_nonneg (Real.rpow_nonneg (by positivity) _) (by linarith)
          linarith
      _ = _ := by simp [hh]
  have hb : Tendsto (fun m : ℕ => c + (m : ℝ)) atTop atTop :=
    tendsto_atTop_add_const_left _ _ tendsto_natCast_atTop_atTop
  have hlim1 : Tendsto (fun m : ℕ => ∫ x in c..(c + m), f x) atTop (𝓝 (∫ x in Ioi c, f x)) :=
    intervalIntegral_tendsto_integral_Ioi _ hfi hb
  have hshift : ∀ j : ℕ, Tendsto (fun m : ℕ => gf s j (c + m)) atTop (𝓝 0) := by
    intro j
    have hj : (0 : ℝ) ≤ j := Nat.cast_nonneg j
    have := (tendsto_cpow_shift hc (w := -s - 1 - (j : ℂ))
      (by simp only [sub_re, neg_re, one_re, natCast_re]; linarith)).const_mul (gc s j)
    rw [mul_zero] at this
    exact this
  have hlim2 : Tendsto (fun m : ℕ => EMPhi s (c + m)) atTop (𝓝 0) := by
    have := ((((hshift 0).const_mul (1 / 12 : ℂ)).sub ((hshift 2).const_mul (1 / 720 : ℂ))).add
      ((hshift 4).const_mul (1 / 30240 : ℂ))).sub ((hshift 6).const_mul (1 / 1209600 : ℂ))
    simp only [mul_zero, sub_zero, add_zero] at this
    exact this
  have hlim : Tendsto (fun m : ℕ => ∑ i ∈ Finset.range m, R i) atTop
      (𝓝 ((∫ x in Ioi c, f x) - EMPhi s c + 0)) := by
    refine ((hlim1.sub_const (EMPhi s c)).add hlim2).congr' ?_
    filter_upwards with m
    rw [key m]; ring
  rw [add_zero] at hlim
  exact le_of_tendsto' hlim.norm hbound

/-- **Order-8 Euler–Maclaurin for ζ(s, a)** on Re s > 0, s ≠ 1, with c = J + a. -/
theorem hurwitz_EM8 {a : ℝ} (ha : a ∈ Ioc (0 : ℝ) 1) (J : ℕ) {s : ℂ} (hs : 0 < s.re)
    (hs1 : s ≠ 1) :
    ‖hurwitzZeta a s - (∑ j ∈ Finset.range J, (((j : ℝ) + a : ℝ) : ℂ) ^ (-s)
      + (-(((J : ℝ) + a : ℝ) : ℂ) ^ (1 - s) / (1 - s) + (((J : ℝ) + a : ℝ) : ℂ) ^ (-s) / 2
        + s * EMPhi s ((J : ℝ) + a)))‖ ≤
      ‖s‖ * (‖gc s 7‖ / 1000000 * (((J : ℝ) + a) ^ (-(s.re + 7)) / (s.re + 7))) := by
  have hc : 0 < (J : ℝ) + a := by have := ha.1; positivity
  rw [hurwitz_first_order ha J hs hs1, hurwitzFO]
  have := saw_integral_EM8 hc hs
  set I := ∫ y in Ioi ((J : ℝ) + a), (saw a y : ℂ) * (y : ℂ) ^ (-(s + 1))
  have e : ∑ j ∈ Finset.range J, (((j : ℝ) + a : ℝ) : ℂ) ^ (-s)
      + (-(((J : ℝ) + a : ℝ) : ℂ) ^ (1 - s) / (1 - s) + (((J : ℝ) + a : ℝ) : ℂ) ^ (-s) / 2 + s * I)
      - (∑ j ∈ Finset.range J, (((j : ℝ) + a : ℝ) : ℂ) ^ (-s)
      + (-(((J : ℝ) + a : ℝ) : ℂ) ^ (1 - s) / (1 - s) + (((J : ℝ) + a : ℝ) : ℂ) ^ (-s) / 2
        + s * EMPhi s ((J : ℝ) + a))) = s * (I - EMPhi s ((J : ℝ) + a)) := by ring
  rw [e, norm_mul]
  exact mul_le_mul_of_nonneg_left this (norm_nonneg _)

end Fayed.Beta
