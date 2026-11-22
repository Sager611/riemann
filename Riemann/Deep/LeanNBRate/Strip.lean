import Riemann.Deep.LeanNBRate.Rect
import Mathlib.MeasureTheory.Integral.IntegralEqImproper

/-!
# Vertical lines, half-strip boundary integrals, and contour shifts across strips
-/

open Complex Filter Topology Set MeasureTheory intervalIntegral Real
open scoped Interval

noncomputable section

namespace Riemann.Deep.LeanNBRate

/-- The vertical line integral of f along Re = sigma (dy, no factor I). -/
def VI (f : ℂ → ℂ) (σ : ℝ) : ℂ := ∫ y : ℝ, f (σ + y * I)

/-- Boundary integral of the upper half-strip [a,b]×[T,∞), positively oriented. -/
def UpperU (f : ℂ → ℂ) (a b T : ℝ) : ℂ :=
  (∫ x : ℝ in a..b, f (x + T * I)) + I * (∫ y in Ioi T, f (b + y * I)) - I * ∫ y in Ioi T, f (a + y * I)

/-- Boundary integral of the lower half-strip [a,b]×(-∞,-T], negatively oriented. -/
def LowerU (f : ℂ → ℂ) (a b T : ℝ) : ℂ :=
  (∫ x : ℝ in a..b, f (x + ((-T : ℝ) : ℂ) * I)) - I * (∫ y in Iic (-T), f (b + y * I)) +
    I * ∫ y in Iic (-T), f (a + y * I)

theorem integral_split3 {g : ℝ → ℂ} (hg : Integrable g) (T : ℝ) :
    ∫ y, g y = (∫ y in Iic (-T), g y) + (∫ y in (-T)..T, g y) + ∫ y in Ioi T, g y := by
  rw [← intervalIntegral.integral_Iic_sub_Iic hg.integrableOn hg.integrableOn,
    ← integral_Iic_add_Ioi (b := T) hg.integrableOn hg.integrableOn]
  ring

/-- Decomposition (as in PNT+'s DiffVertRect_eq_UpperLowerUs). -/
theorem VI_sub_VI_eq {f : ℂ → ℂ} {a b : ℝ} (T : ℝ) (ha : Integrable fun y : ℝ => f (a + y * I))
    (hb : Integrable fun y : ℝ => f (b + y * I)) :
    I * VI f b - I * VI f a = RectI f ⟨a, -T⟩ ⟨b, T⟩ + UpperU f a b T - LowerU f a b T := by
  simp only [VI, RectI, UpperU, LowerU]
  rw [integral_split3 ha T, integral_split3 hb T]
  ring

theorem tendsto_integral_Ioi_zero {g : ℝ → ℂ} (hg : Integrable g) :
    Tendsto (fun T : ℝ => ∫ y in Ioi T, g y) atTop (𝓝 0) := by
  have h : Tendsto (fun T : ℝ => ∫ y, (Ioi T).indicator g y) atTop (𝓝 (∫ _y : ℝ, (0 : ℂ))) := by
    refine tendsto_integral_filter_of_dominated_convergence (fun y => ‖g y‖) ?_ ?_ hg.norm ?_
    · exact Filter.Eventually.of_forall fun T => (hg.indicator measurableSet_Ioi).aestronglyMeasurable
    · refine Filter.Eventually.of_forall fun T => Filter.Eventually.of_forall fun y => ?_
      by_cases hy : y ∈ Ioi T
      · rw [indicator_of_mem hy]
      · rw [indicator_of_notMem hy, norm_zero]; exact norm_nonneg _
    · refine Filter.Eventually.of_forall fun y => ?_
      refine tendsto_const_nhds.congr' ?_
      filter_upwards [Filter.eventually_ge_atTop y] with T hT
      rw [indicator_of_notMem (show y ∉ Ioi T from fun h => by simp at h; linarith)]
  simp only [integral_zero] at h
  refine h.congr fun T => ?_
  rw [MeasureTheory.integral_indicator measurableSet_Ioi]

theorem tendsto_integral_Iic_zero {g : ℝ → ℂ} (hg : Integrable g) :
    Tendsto (fun T : ℝ => ∫ y in Iic (-T), g y) atTop (𝓝 0) := by
  have h : Tendsto (fun T : ℝ => ∫ y, (Iic (-T)).indicator g y) atTop (𝓝 (∫ _y : ℝ, (0 : ℂ))) := by
    refine tendsto_integral_filter_of_dominated_convergence (fun y => ‖g y‖) ?_ ?_ hg.norm ?_
    · exact Filter.Eventually.of_forall fun T => (hg.indicator measurableSet_Iic).aestronglyMeasurable
    · refine Filter.Eventually.of_forall fun T => Filter.Eventually.of_forall fun y => ?_
      by_cases hy : y ∈ Iic (-T)
      · rw [indicator_of_mem hy]
      · rw [indicator_of_notMem hy, norm_zero]; exact norm_nonneg _
    · refine Filter.Eventually.of_forall fun y => ?_
      refine tendsto_const_nhds.congr' ?_
      filter_upwards [Filter.eventually_ge_atTop (-y + 1)] with T hT
      rw [indicator_of_notMem (show y ∉ Iic (-T) from fun h => by simp at h; linarith)]
  simp only [integral_zero] at h
  refine h.congr fun T => ?_
  rw [MeasureTheory.integral_indicator measurableSet_Iic]

theorem norm_hint_le {f : ℂ → ℂ} {a b ε : ℝ} (hab : a ≤ b) {y : ℝ}
    (h : ∀ x ∈ Icc a b, ‖f (x + y * I)‖ ≤ ε) : ‖∫ x : ℝ in a..b, f (x + y * I)‖ ≤ ε * (b - a) := by
  have := intervalIntegral.norm_integral_le_of_norm_le_const (a := a) (b := b)
    (f := fun x : ℝ => f (x + y * I)) (C := ε) (fun x hx => h x (by
      rw [Set.uIoc_of_le hab] at hx; exact ⟨hx.1.le, hx.2⟩))
  rwa [abs_of_nonneg (by linarith)] at this

theorem tendsto_hint_zero {f : ℂ → ℂ} {a b : ℝ} (hab : a ≤ b) {Y : ℝ → ℝ}
    (hdec : ∀ ε > 0, ∀ᶠ T : ℝ in atTop, ∀ x ∈ Icc a b, ‖f (x + Y T * I)‖ ≤ ε) :
    Tendsto (fun T => ∫ x : ℝ in a..b, f (x + Y T * I)) atTop (𝓝 0) := by
  rw [Metric.tendsto_nhds]
  intro ε hε
  have hε' : 0 < ε / (b - a + 1) := by positivity
  filter_upwards [hdec _ hε'] with T hT
  rw [dist_zero_right]
  calc ‖∫ x : ℝ in a..b, f (x + Y T * I)‖ ≤ ε / (b - a + 1) * (b - a) := norm_hint_le hab hT
    _ < ε := by
        rw [div_mul_eq_mul_div, div_lt_iff₀ (by linarith)]
        nlinarith

/-- **Contour shift across a strip with finitely many simple poles.** -/
theorem strip_residues {f : ℂ → ℂ} {a b : ℝ} (hab : a < b) (P : Finset ℂ) (A : ℂ → ℂ)
    (hP : ∀ p ∈ P, a < p.re ∧ p.re < b)
    (hc : ∀ x : ℂ, a ≤ x.re → x.re ≤ b → x ∉ P → ContinuousAt f x)
    (hd : ∀ x : ℂ, a < x.re → x.re < b → x ∉ P → DifferentiableAt ℂ f x)
    (hpole : ∀ p ∈ P, ∃ g : ℂ → ℂ, ContinuousAt g p ∧ ∀ᶠ s in 𝓝[≠] p, f s = A p / (s - p) + g s)
    (ha : Integrable fun y : ℝ => f (a + y * I)) (hb : Integrable fun y : ℝ => f (b + y * I))
    (hdec : ∀ ε > 0, ∀ᶠ T : ℝ in atTop, ∀ x ∈ Icc a b, ‖f (x + T * I)‖ ≤ ε ∧
      ‖f (x + ((-T : ℝ) : ℂ) * I)‖ ≤ ε) :
    I * VI f b - I * VI f a = 2 * π * I * ∑ p ∈ P, A p := by
  obtain ⟨T₁, hT₁⟩ : ∃ T₁ : ℝ, 0 < T₁ ∧ ∀ p ∈ P, |p.im| < T₁ := by
    refine ⟨1 + ∑ p ∈ P, |p.im|, by positivity, fun p hp => ?_⟩
    have := Finset.single_le_sum (f := fun p : ℂ => |p.im|) (fun q _ => abs_nonneg _) hp
    linarith
  have hrect : ∀ T ≥ T₁, RectI f ⟨a, -T⟩ ⟨b, T⟩ = 2 * π * I * ∑ p ∈ P, A p := by
    intro T hT
    have hT0 : 0 < T := lt_of_lt_of_le hT₁.1 hT
    refine RectI_eq_sum_residues (z := ⟨a, -T⟩) (w := ⟨b, T⟩) hab (by simp; linarith) P A ?_ ?_ ?_ hpole
    · intro p hp
      obtain ⟨h1, h2⟩ := hP p hp
      have h3 := hT₁.2 p hp
      rw [abs_lt] at h3
      exact ⟨h1, h2, by simp; linarith, by simp; linarith⟩
    · intro x hx hxP
      rw [Complex.mem_reProdIm] at hx
      simp only [Set.uIcc_of_le hab.le] at hx
      exact hc x hx.1.1 hx.1.2 hxP
    · intro x hx hxP
      rw [Complex.mem_reProdIm] at hx
      exact hd x hx.1.1 hx.1.2 hxP
  have hU : Tendsto (fun T => UpperU f a b T) atTop (𝓝 0) := by
    have h1 : Tendsto (fun T : ℝ => ∫ x : ℝ in a..b, f (x + T * I)) atTop (𝓝 0) :=
      tendsto_hint_zero (Y := fun T => T) hab.le fun ε hε => (hdec ε hε).mono fun T hT x hx => (hT x hx).1
    have := (h1.add ((tendsto_integral_Ioi_zero hb).const_mul I)).sub
      ((tendsto_integral_Ioi_zero ha).const_mul I)
    simpa [UpperU] using this
  have hL : Tendsto (fun T => LowerU f a b T) atTop (𝓝 0) := by
    have h1 : Tendsto (fun T => ∫ x : ℝ in a..b, f (x + ((-T : ℝ) : ℂ) * I)) atTop (𝓝 0) :=
      tendsto_hint_zero (Y := fun T => -T) hab.le fun ε hε => (hdec ε hε).mono fun T hT x hx => (hT x hx).2
    have := (h1.sub ((tendsto_integral_Iic_zero hb).const_mul I)).add
      ((tendsto_integral_Iic_zero ha).const_mul I)
    simpa [LowerU] using this
  have hlim : Tendsto (fun T => I * VI f b - I * VI f a - (UpperU f a b T - LowerU f a b T)) atTop
      (𝓝 (I * VI f b - I * VI f a)) := by
    have := tendsto_const_nhds (x := I * VI f b - I * VI f a) |>.sub (hU.sub hL)
    simpa using this
  have hev : ∀ᶠ T in atTop, I * VI f b - I * VI f a - (UpperU f a b T - LowerU f a b T) =
      2 * π * I * ∑ p ∈ P, A p := by
    filter_upwards [Filter.eventually_ge_atTop T₁] with T hT
    rw [VI_sub_VI_eq T ha hb, hrect T hT]
    ring
  exact tendsto_nhds_unique hlim (tendsto_const_nhds.congr' (hev.mono fun T hT => hT.symm))

/-- **Cauchy for the upper half-strip** (no poles): the boundary integral vanishes. -/
theorem upperU_eq_zero {f : ℂ → ℂ} {a b T : ℝ} (hab : a < b)
    (hc : ∀ x : ℂ, a ≤ x.re → x.re ≤ b → T ≤ x.im → ContinuousAt f x)
    (hd : ∀ x : ℂ, a < x.re → x.re < b → T < x.im → DifferentiableAt ℂ f x)
    (ha : IntegrableOn (fun y : ℝ => f (a + y * I)) (Ioi T))
    (hb : IntegrableOn (fun y : ℝ => f (b + y * I)) (Ioi T))
    (hdec : ∀ ε > 0, ∀ᶠ R : ℝ in atTop, ∀ x ∈ Icc a b, ‖f (x + R * I)‖ ≤ ε) :
    UpperU f a b T = 0 := by
  have hrect : ∀ R > T, RectI f ⟨a, T⟩ ⟨b, R⟩ = 0 := by
    intro R hR
    refine RectI_eq_zero_of_off_countable ∅ Set.countable_empty ?_ ?_
    · intro x hx
      rw [Complex.mem_reProdIm] at hx
      simp only [Set.uIcc_of_le hab.le, Set.uIcc_of_le hR.le] at hx
      exact (hc x hx.1.1 hx.1.2 hx.2.1).continuousWithinAt
    · intro x hx
      rw [Set.sdiff_empty, Complex.mem_reProdIm] at hx
      simp only [min_eq_left hab.le, max_eq_right hab.le, min_eq_left hR.le, max_eq_right hR.le] at hx
      exact hd x hx.1.1 hx.1.2 hx.2.1
  have h1 : Tendsto (fun R : ℝ => ∫ x : ℝ in a..b, f (x + R * I)) atTop (𝓝 0) :=
    tendsto_hint_zero (Y := fun R => R) hab.le hdec
  have h2 := intervalIntegral_tendsto_integral_Ioi T hb tendsto_id
  have h3 := intervalIntegral_tendsto_integral_Ioi T ha tendsto_id
  have hlim : Tendsto (fun R => RectI f ⟨a, T⟩ ⟨b, R⟩) atTop (𝓝 (UpperU f a b T)) := by
    have := ((tendsto_const_nhds (x := ∫ x : ℝ in a..b, f (x + T * I))).sub h1).add
      ((h2.const_mul I).sub (h3.const_mul I))
    simp only [sub_zero] at this
    have e : UpperU f a b T = (∫ x : ℝ in a..b, f (x + T * I)) +
        (I * (∫ y in Ioi T, f (b + y * I)) - I * ∫ y in Ioi T, f (a + y * I)) := by
      simp only [UpperU]; ring
    rw [e]
    refine this.congr fun R => ?_
    simp only [RectI, id]
    ring
  have hlim0 : Tendsto (fun R => RectI f ⟨a, T⟩ ⟨b, R⟩) atTop (𝓝 0) :=
    tendsto_const_nhds.congr' ((Filter.eventually_gt_atTop T).mono fun R hR => (hrect R hR).symm)
  exact tendsto_nhds_unique hlim hlim0

/-- **Cauchy for the lower half-strip** (no poles). -/
theorem lowerU_eq_zero {f : ℂ → ℂ} {a b T : ℝ} (hab : a < b)
    (hc : ∀ x : ℂ, a ≤ x.re → x.re ≤ b → x.im ≤ -T → ContinuousAt f x)
    (hd : ∀ x : ℂ, a < x.re → x.re < b → x.im < -T → DifferentiableAt ℂ f x)
    (ha : IntegrableOn (fun y : ℝ => f (a + y * I)) (Iic (-T)))
    (hb : IntegrableOn (fun y : ℝ => f (b + y * I)) (Iic (-T)))
    (hdec : ∀ ε > 0, ∀ᶠ R : ℝ in atTop, ∀ x ∈ Icc a b, ‖f (x + ((-R : ℝ) : ℂ) * I)‖ ≤ ε) :
    LowerU f a b T = 0 := by
  have hrect : ∀ R > T, RectI f ⟨a, -R⟩ ⟨b, -T⟩ = 0 := by
    intro R hR
    have hR' : -R ≤ -T := by linarith
    refine RectI_eq_zero_of_off_countable ∅ Set.countable_empty ?_ ?_
    · intro x hx
      rw [Complex.mem_reProdIm] at hx
      simp only [Set.uIcc_of_le hab.le, Set.uIcc_of_le hR'] at hx
      exact (hc x hx.1.1 hx.1.2 hx.2.2).continuousWithinAt
    · intro x hx
      rw [Set.sdiff_empty, Complex.mem_reProdIm] at hx
      simp only [min_eq_left hab.le, max_eq_right hab.le, min_eq_left hR', max_eq_right hR'] at hx
      exact hd x hx.1.1 hx.1.2 hx.2.2
  have h1 : Tendsto (fun R => ∫ x : ℝ in a..b, f (x + ((-R : ℝ) : ℂ) * I)) atTop (𝓝 0) :=
    tendsto_hint_zero (Y := fun R => -R) hab.le hdec
  have h2 := intervalIntegral_tendsto_integral_Iic (-T) hb (tendsto_neg_atTop_atBot)
  have h3 := intervalIntegral_tendsto_integral_Iic (-T) ha (tendsto_neg_atTop_atBot)
  have hlim : Tendsto (fun R => RectI f ⟨a, -R⟩ ⟨b, -T⟩) atTop (𝓝 (-LowerU f a b T)) := by
    have := (h1.sub (tendsto_const_nhds (x := ∫ x : ℝ in a..b, f (x + ((-T : ℝ) : ℂ) * I)))).add
      ((h2.const_mul I).sub (h3.const_mul I))
    have e : -LowerU f a b T = (0 - ∫ x : ℝ in a..b, f (x + ((-T : ℝ) : ℂ) * I)) +
        (I * (∫ y in Iic (-T), f (b + y * I)) - I * ∫ y in Iic (-T), f (a + y * I)) := by
      simp only [LowerU]; ring
    rw [e]
    refine this.congr fun R => ?_
    simp only [RectI]
    ring
  have hlim0 : Tendsto (fun R => RectI f ⟨a, -R⟩ ⟨b, -T⟩) atTop (𝓝 0) :=
    tendsto_const_nhds.congr' ((Filter.eventually_gt_atTop T).mono fun R hR => (hrect R hR).symm)
  have := tendsto_nhds_unique hlim hlim0
  exact neg_eq_zero.mp this

end Riemann.Deep.LeanNBRate

