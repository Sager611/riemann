module

public import Riemann.Deep.LeanLi.SmoothSum

/-!
# LeanLi: integration by parts against `S(t) dt = dS₁(t)`, `S₁(t) = ∫_0^t S`
-/

@[expose] public section

open Real Filter Set Finset MeasureTheory Topology

namespace LeanLi

/-- `S₁(t) = ∫_0^t S(u) du`. -/
noncomputable def S1 (t : ℝ) : ℝ := ∫ u in (0 : ℝ)..t, zetaS u

lemma continuous_S1 : Continuous S1 :=
  intervalIntegral.continuous_primitive (fun a b => intervalIntegrable_zetaS a b) 0

/-- Off the ordinates, the averaged count is locally constant. -/
lemma zeroCountAvg_eventuallyEq {x : ℝ} (hx : 0 < x) (hno : ∀ ρ ∈ NontrivialZeros, ρ.im ≠ x) :
    ∀ᶠ y in 𝓝 x, zeroCountAvg y = zeroCountAvg x := by
  set O : Finset ℝ := (zset (x + 1)).image Complex.im with hO
  have hxO : x ∉ O := by
    intro h; rw [hO, Finset.mem_image] at h
    obtain ⟨ρ, hρ, hρx⟩ := h
    rw [mem_zset] at hρ
    exact hno ρ hρ.1 hρx
  have hopen : IsOpen ((↑O : Set ℝ)ᶜ ∩ Ioo 0 (x + 1)) :=
    O.finite_toSet.isClosed.isOpen_compl.inter isOpen_Ioo
  have hmem : x ∈ (↑O : Set ℝ)ᶜ ∩ Ioo 0 (x + 1) := ⟨hxO, hx, by linarith⟩
  obtain ⟨ε, hε, hball⟩ := Metric.isOpen_iff.1 hopen x hmem
  have hgood : ∀ z, |z - x| < ε → ∀ ρ ∈ NontrivialZeros, ρ.im ≠ z := by
    intro z hz ρ hρ heq
    have hzU := hball (show z ∈ Metric.ball x ε by rw [Metric.mem_ball, Real.dist_eq]; exact hz)
    apply hzU.1
    rw [hO, Finset.coe_image, Set.mem_image]
    refine ⟨ρ, ?_, heq⟩
    rw [Finset.mem_coe, mem_zset]
    exact ⟨hρ, by rw [heq]; exact hzU.2.1, by rw [heq]; exact hzU.2.2.le⟩
  filter_upwards [Metric.ball_mem_nhds x hε] with y hy
  rw [Metric.mem_ball, Real.dist_eq] at hy
  rw [avg_eq_of_not_ordinate (hgood y hy), avg_eq_of_not_ordinate hno, zeroCount_eq, zeroCount_eq]
  have hzs : zset y = zset x := by
    ext ρ
    rw [mem_zset, mem_zset]
    simp only [zerosUpTo, Set.mem_setOf_eq]
    constructor
    · rintro ⟨h1, h2, h3⟩
      refine ⟨h1, h2, ?_⟩
      by_contra h4; push Not at h4
      exact hgood ρ.im (by rw [abs_lt]; constructor <;> [linarith; linarith [abs_lt.1 hy]]) ρ h1 rfl
    · rintro ⟨h1, h2, h3⟩
      refine ⟨h1, h2, ?_⟩
      by_contra h4; push Not at h4
      exact hgood ρ.im (by
        rw [abs_lt]; constructor
        · have := (abs_lt.1 hy).1; linarith
        · linarith) ρ h1 rfl
  rw [hzs]

lemma measurable_zetaS : Measurable zetaS := by
  have h1 : Measurable zeroCountAvg := zeroCountAvg_mono.measurable
  have h2 : Continuous Afun := continuous_iff_continuousAt.2 fun t => (hasDerivAt_Afun t).continuousAt
  have : zetaS = fun t => zeroCountAvg t - Afun t := funext zetaS_eq
  rw [this]
  exact h1.sub h2.measurable

lemma hasDerivAt_S1 {x : ℝ} (hx : 0 < x) (hno : ∀ ρ ∈ NontrivialZeros, ρ.im ≠ x) :
    HasDerivAt S1 (zetaS x) x := by
  have hc : ContinuousAt zetaS x := by
    have h2 : Continuous Afun := continuous_iff_continuousAt.2 fun t => (hasDerivAt_Afun t).continuousAt
    have he : zeroCountAvg =ᶠ[𝓝 x] fun _ => zeroCountAvg x := zeroCountAvg_eventuallyEq hx hno
    have h1 : ContinuousAt zeroCountAvg x := (continuousAt_const (y := zeroCountAvg x)).congr he.symm
    have : zetaS = fun t => zeroCountAvg t - Afun t := funext zetaS_eq
    rw [this]
    exact h1.sub h2.continuousAt
  exact intervalIntegral.integral_hasDerivAt_right (intervalIntegrable_zetaS 0 x)
    measurable_zetaS.aestronglyMeasurable.stronglyMeasurableAtFilter hc

/-- **Integration by parts**: `∫_a^b S G = S₁(b)G(b) - S₁(a)G(a) - ∫_a^b S₁ G'`. -/
theorem ibp_S {a b : ℝ} (ha : 0 ≤ a) (hab : a ≤ b) (G G' : ℝ → ℝ) (hG : ∀ t, HasDerivAt G (G' t) t)
    (hG' : Continuous G') :
    ∫ t in a..b, zetaS t * G t = S1 b * G b - S1 a * G a - ∫ t in a..b, S1 t * G' t := by
  have hGc : Continuous G := continuous_iff_continuousAt.2 fun t => (hG t).continuousAt
  set s : Set ℝ := ↑((zset (b + 1)).image Complex.im) with hs
  have hsc : s.Countable := (Finset.finite_toSet _).countable
  have hi1 : IntervalIntegrable (fun t => zetaS t * G t) volume a b :=
    (intervalIntegrable_zetaS a b).mul_continuousOn hGc.continuousOn
  have hi2 : IntervalIntegrable (fun t => S1 t * G' t) volume a b :=
    (continuous_S1.mul hG').intervalIntegrable a b
  have key := integral_eq_of_hasDerivAt_off_countable_of_le (fun t => S1 t * G t)
    (fun t => zetaS t * G t + S1 t * G' t) hab hsc ((continuous_S1.mul hGc).continuousOn)
    (fun x hx => by
      obtain ⟨hx1, hx2⟩ := hx
      have hx0 : 0 < x := lt_of_le_of_lt ha hx1.1
      have hno : ∀ ρ ∈ NontrivialZeros, ρ.im ≠ x := by
        intro ρ hρ heq
        apply hx2
        rw [hs, Finset.coe_image, Set.mem_image]
        refine ⟨ρ, ?_, heq⟩
        rw [Finset.mem_coe, mem_zset]
        exact ⟨hρ, by rw [heq]; exact hx0, by rw [heq]; linarith [hx1.2]⟩
      exact (hasDerivAt_S1 hx0 hno).mul (hG x))
    (hi1.add hi2)
  rw [intervalIntegral.integral_add hi1 hi2] at key
  linarith

end LeanLi
