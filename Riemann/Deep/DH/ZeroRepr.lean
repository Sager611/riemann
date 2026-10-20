import Riemann.Deep.DH.FE
import Riemann.Deep.FirstZero.NextLib

/-!
# Davenport–Heilbronn zero, part 1: a convergent representation of DH on Re s > 0

For `M : ℕ`, Re s > 0, s ≠ 1:

  DH s = Σ_{m<M} W s m + 5^{-s} · (R(1/5) − R(4/5) + κ (R(2/5) − R(3/5))),

`W s m = (5m+1)^{-s} − (5m+4)^{-s} + κ ((5m+2)^{-s} − (5m+3)^{-s})`,
`R(a) = A(M+a) + s · J_a(s)`, `A(x) = x^{1-s}/(s-1) + x^{-s}/2`,
`J_a(s) = ∫_{M+a}^∞ (⌊y-a⌋ + 1/2 - (y-a)) y^{-s-1} dy`  (first-order Euler–Maclaurin, Hurwitz form).

Proof: equality for Re s > 1 from the Dirichlet series (blocks of 5 terms, per-class telescoping),
then the identity theorem on {Re s > 0, s ≠ 1}; `J_a` is analytic there because it is a Mellin
transform of a bounded function (`mellin_differentiableAt_of_isBigO_rpow`).
-/

open Complex Set Filter Topology MeasureTheory Interval Asymptotics

namespace Fayed.DH.Zero

open Fayed.DH
open Fayed.FirstZero (nps)

/-! ### The shifted sawtooth and its Mellin integral -/

/-- `saw a y = ⌊y - a⌋ + 1/2 - (y - a)`. -/
noncomputable def saw (a y : ℝ) : ℝ := ⌊y - a⌋ + 1 / 2 - (y - a)

lemma abs_saw_le (a y : ℝ) : |saw a y| ≤ 1 / 2 := by
  unfold saw
  have h1 := Int.floor_le (y - a)
  have h2 := Int.lt_floor_add_one (y - a)
  rw [abs_le]; constructor <;> linarith

lemma measurable_saw (a : ℝ) : Measurable (saw a) := by
  unfold saw
  have h : Measurable (fun y : ℝ => ((⌊y - a⌋ : ℤ) : ℝ)) :=
    (measurable_of_countable (fun z : ℤ => (z : ℝ))).comp
      (Int.measurable_floor.comp (measurable_id.sub_const a))
  exact (h.add_const _).sub (measurable_id.sub_const a)

lemma integrableOn_saw {a x0 : ℝ} (hx0 : 0 < x0) {s : ℂ} (hs : 0 < s.re) :
    IntegrableOn (fun y : ℝ => (saw a y : ℂ) * (y : ℂ) ^ (-(s + 1))) (Ioi x0) := by
  refine Integrable.bdd_mul (c := 1 / 2) ?_ ?_ ?_
  · exact (integrableOn_Ioi_cpow_iff hx0).mpr (by simp only [neg_re, add_re, one_re]; linarith)
  · exact (Complex.measurable_ofReal.comp (measurable_saw a)).aestronglyMeasurable
  · exact ae_of_all _ fun y => by rw [Complex.norm_real, Real.norm_eq_abs]; exact abs_saw_le a y

/-- `J a x0 s = ∫_{x0}^∞ saw a y · y^{-(s+1)} dy`. -/
noncomputable def J (a x0 : ℝ) (s : ℂ) : ℂ :=
  ∫ y in Ioi x0, (saw a y : ℂ) * (y : ℂ) ^ (-(s + 1))

noncomputable def sawInd (a x0 : ℝ) : ℝ → ℂ := (Ioi x0).indicator (fun y => (saw a y : ℂ))

lemma J_eq_mellin (a : ℝ) {x0 : ℝ} (hx0 : 0 ≤ x0) (s : ℂ) :
    J a x0 s = mellin (sawInd a x0) (-s) := by
  rw [mellin, J, sawInd]
  have h : (fun t : ℝ => (t : ℂ) ^ (-s - 1) • (Ioi x0).indicator (fun y => (saw a y : ℂ)) t)
      = (Ioi x0).indicator (fun t : ℝ => (saw a t : ℂ) * (t : ℂ) ^ (-(s + 1))) := by
    funext t
    by_cases ht : t ∈ Ioi x0
    · simp only [indicator_of_mem ht, smul_eq_mul]
      rw [show -s - 1 = -(s + 1) by ring, mul_comm]
    · simp [indicator_of_notMem ht]
  rw [h, setIntegral_indicator measurableSet_Ioi, Ioi_inter_Ioi, sup_of_le_right hx0]

lemma norm_sawInd_le (a x0 t : ℝ) : ‖sawInd a x0 t‖ ≤ 1 / 2 := by
  unfold sawInd
  by_cases ht : t ∈ Ioi x0
  · rw [indicator_of_mem ht, Complex.norm_real, Real.norm_eq_abs]; exact abs_saw_le a t
  · rw [indicator_of_notMem ht, norm_zero]; norm_num

lemma differentiableAt_J (a : ℝ) {x0 : ℝ} (hx0 : 0 < x0) {s : ℂ} (hs : 0 < s.re) :
    DifferentiableAt ℂ (J a x0) s := by
  have hfun : J a x0 = fun s => mellin (sawInd a x0) (-s) := funext (J_eq_mellin a hx0.le)
  rw [hfun]
  have hmeas : Measurable (sawInd a x0) :=
    (Complex.measurable_ofReal.comp (measurable_saw a)).indicator measurableSet_Ioi
  have hloc : LocallyIntegrableOn (sawInd a x0) (Ioi 0) := by
    refine (locallyIntegrable_iff.mpr fun K hK => ?_).locallyIntegrableOn _
    exact Measure.integrableOn_of_bounded hK.measure_lt_top.ne hmeas.aestronglyMeasurable
      (M := 1 / 2) (ae_of_all _ (norm_sawInd_le a x0))
  have htop : sawInd a x0 =O[atTop] (fun t : ℝ => t ^ (-(0 : ℝ))) := by
    refine IsBigO.of_bound 1 (Eventually.of_forall fun t => ?_)
    simp only [neg_zero, Real.rpow_zero, norm_one, mul_one]
    linarith [norm_sawInd_le a x0 t]
  have hbot : sawInd a x0 =O[𝓝[>] 0] (fun t : ℝ => t ^ (-(-s.re - 1))) := by
    refine IsBigO.of_bound 0 ?_
    filter_upwards [Ioo_mem_nhdsGT hx0] with t ht
    have : t ∉ Ioi x0 := fun h => by simp only [mem_Ioi] at h; linarith [ht.2]
    simp [sawInd, indicator_of_notMem this]
  have hD := mellin_differentiableAt_of_isBigO_rpow (s := -s) hloc htop
    (by simp only [neg_re]; linarith) hbot (by simp only [neg_re]; linarith)
  exact hD.comp s differentiableAt_id.neg

/-! ### Per-class telescoping -/

/-- `Af s x = x^{1-s}/(s-1) + x^{-s}/2`. -/
noncomputable def Af (s : ℂ) (x : ℝ) : ℂ := (x : ℂ) ^ (1 - s) / (s - 1) + (x : ℂ) ^ (-s) / 2

lemma hasDerivAt_ofReal_cpow {c : ℂ} {y : ℝ} (hy : 0 < y) :
    HasDerivAt (fun x : ℝ => (x : ℂ) ^ c) (c * (y : ℂ) ^ (c - 1)) y :=
  ((Complex.hasStrictDerivAt_cpow_const (c := c) (ofReal_mem_slitPlane.mpr hy)).hasDerivAt).comp_ofReal

lemma cpow_one_sub' {x : ℝ} (hx : 0 < x) (s : ℂ) : (x : ℂ) ^ (1 - s) = x * (x : ℂ) ^ (-s) := by
  rw [show (1 - s) = 1 + (-s) by ring, cpow_add _ _ (by exact_mod_cast hx.ne'), cpow_one]

lemma piece_eq {s : ℂ} (hs1 : s ≠ 1) {n : ℝ} (hn : 0 < n) :
    s * ∫ y in n..n + 1, ((n : ℂ) + 1 / 2 - y) * (y : ℂ) ^ (-(s + 1)) =
      (n : ℂ) ^ (-s) - (Af s n - Af s (n + 1)) := by
  have hs1' : s - 1 ≠ 0 := sub_ne_zero.mpr hs1
  have hpos : ∀ y ∈ [[n, n + 1]], 0 < y := by
    intro y hy; rw [uIcc_of_le (by linarith)] at hy; linarith [hy.1]
  set F : ℝ → ℂ := fun y => -((n : ℂ) + 1 / 2) * (y : ℂ) ^ (-s) + s / (s - 1) * (y : ℂ) ^ (1 - s)
    with hFdef
  have hF : ∀ y ∈ [[n, n + 1]],
      HasDerivAt F (s * (((n : ℂ) + 1 / 2 - y) * (y : ℂ) ^ (-(s + 1)))) y := by
    intro y hy
    have hy0 := hpos y hy
    have hy' : (y : ℂ) ≠ 0 := by exact_mod_cast hy0.ne'
    have h1 := (hasDerivAt_ofReal_cpow (c := -s) hy0).const_mul (-((n : ℂ) + 1 / 2))
    have h2 := (hasDerivAt_ofReal_cpow (c := 1 - s) hy0).const_mul (s / (s - 1))
    convert h1.add h2 using 1
    have e1 : (y : ℂ) ^ (-s - 1) = (y : ℂ) ^ (-(s + 1)) := by ring_nf
    have e2 : (y : ℂ) ^ (1 - s - 1) = (y : ℂ) * (y : ℂ) ^ (-(s + 1)) := by
      rw [show (1 - s - 1 : ℂ) = 1 + (-(s + 1)) by ring, cpow_add _ _ hy', cpow_one]
    rw [e1, e2]
    field_simp
    ring
  have hint : IntervalIntegrable
      (fun y : ℝ => s * (((n : ℂ) + 1 / 2 - y) * (y : ℂ) ^ (-(s + 1)))) volume n (n + 1) := by
    apply ContinuousOn.intervalIntegrable
    intro y hy
    apply ContinuousAt.continuousWithinAt
    have hc := (hasDerivAt_ofReal_cpow (c := -(s + 1)) (hpos y hy)).continuousAt
    exact continuousAt_const.mul
      ((continuousAt_const.sub Complex.continuous_ofReal.continuousAt).mul hc)
  rw [← intervalIntegral.integral_const_mul, intervalIntegral.integral_eq_sub_of_hasDerivAt hF hint]
  simp only [hFdef, Af]
  have hn1 : (0 : ℝ) < n + 1 := by linarith
  rw [cpow_one_sub' hn, cpow_one_sub' hn1]
  push_cast
  field_simp
  ring

lemma saw_piece (a : ℝ) (k : ℕ) (s : ℂ) :
    ∫ y in ((k : ℝ) + a)..((k : ℝ) + a) + 1, (saw a y : ℂ) * (y : ℂ) ^ (-(s + 1)) =
      ∫ y in ((k : ℝ) + a)..((k : ℝ) + a) + 1,
        ((((k : ℝ) + a : ℝ) : ℂ) + 1 / 2 - y) * (y : ℂ) ^ (-(s + 1)) := by
  refine intervalIntegral.integral_congr_ae ?_
  have hae : ∀ᵐ x : ℝ, x ≠ ((k : ℝ) + a) + 1 := by
    rw [ae_iff]; simp [measure_singleton]
  filter_upwards [hae] with y hy hmem
  rw [uIoc_of_le (by linarith)] at hmem
  have hfl : ⌊y - a⌋ = (k : ℤ) := by
    rw [Int.floor_eq_iff]; push_cast
    constructor
    · linarith [hmem.1]
    · have := lt_of_le_of_ne hmem.2 hy; linarith
  simp only [saw, hfl]
  push_cast
  ring

lemma class_finite {s : ℂ} (hs : 0 < s.re) (hs1 : s ≠ 1) {M : ℕ} {a : ℝ}
    (hx0 : 0 < (M : ℝ) + a) (K : ℕ) :
    ∑ i ∈ Finset.range K, ((((M + i : ℕ) : ℝ) + a : ℝ) : ℂ) ^ (-s) =
      Af s ((M : ℝ) + a) - Af s (((M + K : ℕ) : ℝ) + a) +
        s * ∫ y in ((M : ℝ) + a)..(((M + K : ℕ) : ℝ) + a),
          (saw a y : ℂ) * (y : ℂ) ^ (-(s + 1)) := by
  have hii : ∀ p q : ℝ, (M : ℝ) + a ≤ p → p ≤ q →
      IntervalIntegrable (fun y : ℝ => (saw a y : ℂ) * (y : ℂ) ^ (-(s + 1))) volume p q := by
    intro p q hp hpq
    rw [intervalIntegrable_iff_integrableOn_Ioc_of_le hpq]
    exact (integrableOn_saw hx0 hs).mono_set (Ioc_subset_Ioi_self.trans (Ioi_subset_Ioi hp))
  induction K with
  | zero => simp
  | succ K ih =>
    rw [Finset.sum_range_succ, ih]
    have hc : (((M + (K + 1) : ℕ) : ℝ) + a) = (((M + K : ℕ) : ℝ) + a) + 1 := by push_cast; ring
    have hle : (M : ℝ) + a ≤ ((M + K : ℕ) : ℝ) + a := by
      push_cast; linarith [(Nat.cast_nonneg K : (0 : ℝ) ≤ K)]
    have hpos : 0 < ((M + K : ℕ) : ℝ) + a := by linarith
    rw [hc, ← intervalIntegral.integral_add_adjacent_intervals (hii _ _ le_rfl hle)
      (hii _ ((((M + K : ℕ) : ℝ) + a) + 1) hle (by linarith)), saw_piece a (M + K) s]
    have hp := piece_eq hs1 hpos
    rw [mul_add, hp]
    ring

lemma cpow_tendsto_zero {c : ℂ} (hc : c.re < 0) :
    Tendsto (fun x : ℝ => (x : ℂ) ^ c) atTop (𝓝 0) := by
  rw [tendsto_zero_iff_norm_tendsto_zero]
  refine (tendsto_rpow_neg_atTop (y := -c.re) (by linarith)).congr' ?_
  filter_upwards [eventually_gt_atTop 0] with x hx
  rw [norm_cpow_eq_rpow_re_of_pos hx, neg_neg]

lemma Af_tendsto {s : ℂ} (hs : 1 < s.re) : Tendsto (Af s) atTop (𝓝 0) := by
  have h1 := (cpow_tendsto_zero (c := 1 - s) (by simp only [sub_re, one_re]; linarith)).div_const
    (s - 1)
  have h2 := (cpow_tendsto_zero (c := -s) (by simp only [neg_re]; linarith)).div_const 2
  have e : Af s = fun x : ℝ => (x : ℂ) ^ (1 - s) / (s - 1) + (x : ℂ) ^ (-s) / 2 := rfl
  rw [e]
  simpa using h1.add h2

lemma class_tendsto {s : ℂ} (hs : 1 < s.re) {M : ℕ} {a : ℝ} (hx0 : 0 < (M : ℝ) + a) :
    Tendsto (fun K => ∑ i ∈ Finset.range K, ((((M + i : ℕ) : ℝ) + a : ℝ) : ℂ) ^ (-s)) atTop
      (𝓝 (Af s ((M : ℝ) + a) + s * J a ((M : ℝ) + a) s)) := by
  have hs0 : 0 < s.re := by linarith
  have hs1 : s ≠ 1 := by intro h; rw [h] at hs; simp at hs
  have hb : Tendsto (fun K : ℕ => ((M + K : ℕ) : ℝ) + a) atTop atTop := by
    refine tendsto_atTop_add_const_right _ a ?_
    exact tendsto_natCast_atTop_atTop.comp
      (tendsto_atTop_mono (fun K => Nat.le_add_left K M) tendsto_id)
  have hA : Tendsto (fun K : ℕ => Af s (((M + K : ℕ) : ℝ) + a)) atTop (𝓝 0) :=
    (Af_tendsto hs).comp hb
  have hI := intervalIntegral_tendsto_integral_Ioi _ (integrableOn_saw (a := a) hx0 hs0) hb
  have h := ((tendsto_const_nhds (x := Af s ((M : ℝ) + a))).sub hA).add (hI.const_mul s)
  rw [sub_zero] at h
  exact h.congr' (Eventually.of_forall fun K => (class_finite hs0 hs1 hx0 K).symm)

/-! ### Blocks of five Dirichlet terms -/

/-- One block `5m+1, …, 5m+4` of the DH Dirichlet series. -/
noncomputable def W (s : ℂ) (m : ℕ) : ℂ :=
  nps s (5 * m + 1) - nps s (5 * m + 4) + (κ : ℂ) * (nps s (5 * m + 2) - nps s (5 * m + 3))

noncomputable def Gmain (M : ℕ) (s : ℂ) : ℂ := ∑ m ∈ Finset.range M, W s m

lemma coeff_shift (K j : ℕ) : coeff (5 * K + j) = coeff j := by
  have h : ((5 * K + j : ℕ) : ZMod 5) = (j : ZMod 5) := by
    push_cast; rw [show (5 : ZMod 5) = 0 from rfl]; ring
  simp only [coeff, h]

lemma coeff_zero' : coeff 0 = 0 := by
  rcases coeff_values 0 with ⟨-, h⟩ | ⟨h1, -⟩ | ⟨h1, -⟩ | ⟨h1, -⟩ | ⟨h1, -⟩
  · exact h
  all_goals exact absurd h1 (by decide)

lemma coeff_one' : coeff 1 = 1 := by
  rcases coeff_values 1 with ⟨h1, -⟩ | ⟨-, h⟩ | ⟨h1, -⟩ | ⟨h1, -⟩ | ⟨h1, -⟩
  · exact absurd h1 (by decide)
  · exact h
  all_goals exact absurd h1 (by decide)

lemma coeff_two' : coeff 2 = κ := by
  rcases coeff_values 2 with ⟨h1, -⟩ | ⟨h1, -⟩ | ⟨-, h⟩ | ⟨h1, -⟩ | ⟨h1, -⟩
  · exact absurd h1 (by decide)
  · exact absurd h1 (by decide)
  · exact h
  all_goals exact absurd h1 (by decide)

lemma coeff_three' : coeff 3 = -κ := by
  rcases coeff_values 3 with ⟨h1, -⟩ | ⟨h1, -⟩ | ⟨h1, -⟩ | ⟨-, h⟩ | ⟨h1, -⟩
  · exact absurd h1 (by decide)
  · exact absurd h1 (by decide)
  · exact absurd h1 (by decide)
  · exact h
  · exact absurd h1 (by decide)

lemma coeff_four' : coeff 4 = -1 := by
  rcases coeff_values 4 with ⟨h1, -⟩ | ⟨h1, -⟩ | ⟨h1, -⟩ | ⟨h1, -⟩ | ⟨-, h⟩
  all_goals first | exact h | exact absurd h1 (by decide)

lemma term_eq (s : ℂ) {n : ℕ} (hn : n ≠ 0) : LSeries.term coeff s n = coeff n * nps s n := by
  rw [LSeries.term_of_ne_zero hn, nps, cpow_neg, div_eq_mul_inv]

lemma term_block (s : ℂ) (K : ℕ) :
    ∑ n ∈ Finset.range (5 * K), LSeries.term coeff s n = ∑ m ∈ Finset.range K, W s m := by
  induction K with
  | zero => simp
  | succ K ih =>
    rw [show 5 * (K + 1) = 5 * K + 1 + 1 + 1 + 1 + 1 by ring]
    simp only [Finset.sum_range_succ]
    rw [ih]
    have t0 : LSeries.term coeff s (5 * K) = 0 := by
      rcases Nat.eq_zero_or_pos K with hK | hK
      · subst hK; simp
      · rw [term_eq s (by omega), show 5 * K = 5 * K + 0 by ring, coeff_shift, coeff_zero',
          zero_mul]
    have t1 : LSeries.term coeff s (5 * K + 1) = nps s (5 * K + 1) := by
      rw [term_eq s (by omega), coeff_shift, coeff_one', one_mul]
    have t2 : LSeries.term coeff s (5 * K + 1 + 1) = κ * nps s (5 * K + 2) := by
      rw [show 5 * K + 1 + 1 = 5 * K + 2 by ring, term_eq s (by omega), coeff_shift, coeff_two']
    have t3 : LSeries.term coeff s (5 * K + 1 + 1 + 1) = -κ * nps s (5 * K + 3) := by
      rw [show 5 * K + 1 + 1 + 1 = 5 * K + 3 by ring, term_eq s (by omega), coeff_shift,
        coeff_three']
    have t4 : LSeries.term coeff s (5 * K + 1 + 1 + 1 + 1) = -nps s (5 * K + 4) := by
      rw [show 5 * K + 1 + 1 + 1 + 1 = 5 * K + 4 by ring, term_eq s (by omega), coeff_shift,
        coeff_four']
      ring
    rw [t0, t1, t2, t3, t4, W]
    ring

lemma nps_shift (s : ℂ) (M i j : ℕ) :
    nps s (5 * (M + i) + j) = (5 : ℂ) ^ (-s) * ((((M + i : ℕ) : ℝ) + j / 5 : ℝ) : ℂ) ^ (-s) := by
  rw [nps, show ((5 * (M + i) + j : ℕ) : ℂ) = ((5 : ℝ) : ℂ) * ((((M + i : ℕ) : ℝ) + j / 5 : ℝ) : ℂ)
    by push_cast; ring, mul_cpow_ofReal_nonneg (by norm_num) (by positivity)]
  norm_cast

/-- `R a = A(M+a) + s·J_a(s)`. -/
noncomputable def Rj (M : ℕ) (a : ℝ) (s : ℂ) : ℂ := Af s ((M : ℝ) + a) + s * J a ((M : ℝ) + a) s

/-- The representation of DH on Re s > 0, s ≠ 1. -/
noncomputable def Grep (M : ℕ) (s : ℂ) : ℂ :=
  Gmain M s + (5 : ℂ) ^ (-s) *
    (Rj M (1 / 5) s - Rj M (4 / 5) s + (κ : ℂ) * (Rj M (2 / 5) s - Rj M (3 / 5) s))

lemma tail_block (s : ℂ) (M K : ℕ) :
    ∑ i ∈ Finset.range K, W s (M + i) = (5 : ℂ) ^ (-s) *
      ((∑ i ∈ Finset.range K, ((((M + i : ℕ) : ℝ) + (1 : ℕ) / 5 : ℝ) : ℂ) ^ (-s))
        - (∑ i ∈ Finset.range K, ((((M + i : ℕ) : ℝ) + (4 : ℕ) / 5 : ℝ) : ℂ) ^ (-s))
        + (κ : ℂ) * ((∑ i ∈ Finset.range K, ((((M + i : ℕ) : ℝ) + (2 : ℕ) / 5 : ℝ) : ℂ) ^ (-s))
          - (∑ i ∈ Finset.range K, ((((M + i : ℕ) : ℝ) + (3 : ℕ) / 5 : ℝ) : ℂ) ^ (-s)))) := by
  induction K with
  | zero => simp
  | succ K ih =>
    simp only [Finset.sum_range_succ]
    rw [ih, W, nps_shift s M K 1, nps_shift s M K 2, nps_shift s M K 3, nps_shift s M K 4]
    ring

lemma DH_eq_Grep_of_one_lt (M : ℕ) {s : ℂ} (hs : 1 < s.re) : DH s = Grep M s := by
  have hsum : LSeriesSummable coeff s := by
    have h := ((DirichletCharacter.LSeriesSummable_of_one_lt_re χ hs).smul a).add
      ((DirichletCharacter.LSeriesSummable_of_one_lt_re χ⁻¹ hs).smul b)
    convert h using 1
    funext n; simp [coeff]
  have hL := hsum.hasSum.tendsto_sum_nat
  have h5 : Tendsto (fun K : ℕ => 5 * (M + K)) atTop atTop :=
    tendsto_atTop_mono (fun K => by simp only [id]; omega) tendsto_id
  have hL' := hL.comp h5
  have hpos : ∀ j : ℕ, 0 < j → 0 < (M : ℝ) + (j : ℕ) / 5 := fun j hj => by
    have : (0 : ℝ) < j := by exact_mod_cast hj
    positivity
  have hT := ((((class_tendsto hs (hpos 1 (by norm_num))).sub
    (class_tendsto hs (hpos 4 (by norm_num)))).add
    (((class_tendsto hs (hpos 2 (by norm_num))).sub
      (class_tendsto hs (hpos 3 (by norm_num)))).const_mul (κ : ℂ))).const_mul
      ((5 : ℂ) ^ (-s))).const_add (Gmain M s)
  have heq : (fun K : ℕ => ∑ n ∈ Finset.range (5 * (M + K)), LSeries.term coeff s n) =
      fun K => Gmain M s + (5 : ℂ) ^ (-s) *
      ((∑ i ∈ Finset.range K, ((((M + i : ℕ) : ℝ) + (1 : ℕ) / 5 : ℝ) : ℂ) ^ (-s))
        - (∑ i ∈ Finset.range K, ((((M + i : ℕ) : ℝ) + (4 : ℕ) / 5 : ℝ) : ℂ) ^ (-s))
        + (κ : ℂ) * ((∑ i ∈ Finset.range K, ((((M + i : ℕ) : ℝ) + (2 : ℕ) / 5 : ℝ) : ℂ) ^ (-s))
          - (∑ i ∈ Finset.range K, ((((M + i : ℕ) : ℝ) + (3 : ℕ) / 5 : ℝ) : ℂ) ^ (-s)))) := by
    funext K
    rw [term_block, Finset.sum_range_add, tail_block, Gmain]
  rw [Function.comp_def, heq] at hL'
  rw [DH_eq_LSeries hs, LSeries, tendsto_nhds_unique hL' hT, Grep, Rj, Rj, Rj, Rj]
  push_cast
  ring

/-! ### Analyticity and the identity theorem -/

lemma differentiable_nps (n : ℕ) (hn : n ≠ 0) : Differentiable ℂ (fun s => nps s n) := fun s =>
  (differentiableAt_id.neg).const_cpow (Or.inl (by exact_mod_cast hn))

lemma differentiable_W (m : ℕ) : Differentiable ℂ (fun s => W s m) := by
  have h := fun j : ℕ => differentiable_nps (5 * m + j + 1) (by omega)
  unfold W
  exact ((h 0).sub (h 3)).add ((h 1).sub (h 2) |>.const_mul _)

lemma differentiable_Gmain (M : ℕ) : Differentiable ℂ (Gmain M) := by
  unfold Gmain
  induction M with
  | zero => simp only [Finset.range_zero, Finset.sum_empty]; exact differentiable_const 0
  | succ M ih =>
    simp only [Finset.sum_range_succ]
    exact ih.add (differentiable_W M)

lemma differentiableAt_Af {x : ℝ} (hx : 0 < x) {s : ℂ} (hs1 : s ≠ 1) :
    DifferentiableAt ℂ (fun s => Af s x) s := by
  have hx' : (x : ℂ) ≠ 0 := by exact_mod_cast hx.ne'
  unfold Af
  refine (((differentiableAt_const _).sub differentiableAt_id).const_cpow (Or.inl hx')).div
    (differentiableAt_id.sub (differentiableAt_const _)) (sub_ne_zero.mpr hs1) |>.add ?_
  exact ((differentiableAt_id.neg).const_cpow (Or.inl hx')).div_const 2

lemma differentiableAt_Rj (M : ℕ) {a : ℝ} (ha : 0 < (M : ℝ) + a) {s : ℂ} (hs : 0 < s.re)
    (hs1 : s ≠ 1) : DifferentiableAt ℂ (Rj M a) s :=
  (differentiableAt_Af ha hs1).add (differentiableAt_id.mul (differentiableAt_J a ha hs))

lemma differentiableAt_Grep (M : ℕ) {s : ℂ} (hs : 0 < s.re) (hs1 : s ≠ 1) :
    DifferentiableAt ℂ (Grep M) s := by
  have hpos : ∀ j : ℝ, 0 < j → 0 < (M : ℝ) + j := fun j hj => by positivity
  have r := fun j (hj : (0 : ℝ) < j) => differentiableAt_Rj M (hpos j hj) hs hs1
  unfold Grep
  refine (differentiable_Gmain M s).add (DifferentiableAt.mul ?_ ?_)
  · exact (differentiableAt_id.neg).const_cpow (Or.inl (by norm_num))
  · exact ((r _ (by norm_num)).sub (r _ (by norm_num))).add
      (((r _ (by norm_num)).sub (r _ (by norm_num))).const_mul _)

lemma differentiable_DH : Differentiable ℂ DH := by
  have h1 := DirichletCharacter.differentiable_LFunction χ_ne_one
  have h2 := DirichletCharacter.differentiable_LFunction (inv_ne_one.mpr χ_ne_one)
  exact (h1.const_mul a).add (h2.const_mul b)

/-- **Representation of DH on Re s > 0, s ≠ 1.** -/
theorem DH_eq_Grep (M : ℕ) {s : ℂ} (hs : 0 < s.re) (hs1 : s ≠ 1) : DH s = Grep M s := by
  let U := {z : ℂ | z ≠ 1 ∧ 0 < z.re}
  have f_an : AnalyticOnNhd ℂ DH U :=
    differentiable_DH.differentiableOn.analyticOnNhd Fayed.FirstZero.PNT.isOpen_aux
  have g_an : AnalyticOnNhd ℂ (Grep M) U :=
    DifferentiableOn.analyticOnNhd
      (fun z hz => (differentiableAt_Grep M hz.2 hz.1).differentiableWithinAt)
      Fayed.FirstZero.PNT.isOpen_aux
  have preconU : IsPreconnected U :=
    ((IsOpen.isConnected_iff_isPathConnected Fayed.FirstZero.PNT.isOpen_aux).mpr
      Fayed.FirstZero.PNT.isPathConnected_aux).isPreconnected
  have h2 : (2 : ℂ) ∈ U := by
    refine ⟨?_, by simp⟩
    intro h; have := congrArg Complex.re h; norm_num at this
  refine AnalyticOnNhd.eqOn_of_preconnected_of_eventuallyEq f_an g_an preconU h2 ?_ ⟨hs1, hs⟩
  have u_mem : {z : ℂ | 1 < z.re} ∈ 𝓝 (2 : ℂ) :=
    (isOpen_lt continuous_const Complex.continuous_re).mem_nhds (by simp)
  filter_upwards [u_mem] with z hz
  exact DH_eq_Grep_of_one_lt M hz

end Fayed.DH.Zero
