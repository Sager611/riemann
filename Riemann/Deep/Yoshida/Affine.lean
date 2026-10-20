/-
Riemann/Deep/Yoshida/Affine.lean — affine cells and the "wider support" form of the moment certificate.

Two changes against `Moment.lean`:
1. **Affine cells** `(lo, hi, A, B)`: the deficit `t` is `A + B|r|` on `lo ≤ |r| < hi`. With the chord bounds
   of `Chord.lean` the loss per cell is second order in the width.
2. **Wider support.** For `f` supported in `[−a, a]` with `2a < log 2` and any `a₁ ≥ a`, the pole term is
   `(1/2π)∫|F|² J(a₁, r) dr` (the kernel `2cosh(x/2)` may be cut anywhere outside the support of the
   autocorrelation), and the Taylor and Bessel steps hold with `a₁`. So one certificate at a rational
   `a₁ ≥ log 2 / 2` gives Weil positivity on **every** window `L < log 2`.
-/
import Riemann.Deep.Yoshida.Chord
import Riemann.Deep.Yoshida.Moment

noncomputable section

open Complex MeasureTheory Set Filter Finset
open scoped ComplexConjugate

namespace Fayed.Yoshida

open Riemann.Deep.Weil Riemann.Deep.WeilEF Zeta23 Zeta23.EF

/-! ### Wider support -/

/-- Below `log 2`, with the pole kernel cut at `2a₁ ≥ 2a`. -/
theorem weilQ_re_eq_W' {a a1 : ℝ} (ha : 0 ≤ a) (ha1 : a ≤ a1) (hlog : 2 * a < Real.log 2)
    {f : ℝ → ℂ} (hf : TestFn a f) :
    Integrable (fun r : ℝ => ‖paperFT f r‖ ^ 2 * W a1 r) ∧
    (weilQ f).re = 1 / (2 * Real.pi) * ∫ r : ℝ, ‖paperFT f r‖ ^ 2 * W a1 r := by
  have h0 : 0 ≤ a1 := ha.trans ha1
  have hk := hf.autocorr_contDiff
  have hkc := hf.autocorr_hasCompactSupport
  have hFk := integrable_fourier_of_contDiff_two hk hkc
  have hφi : Integrable (fun r : ℝ => ‖paperFT f r‖ ^ 2) := by
    refine (integrable_paperFT_ofReal hFk).re.congr (Eventually.of_forall fun r => ?_)
    simp only [paperFT_autocorr_real hf, RCLike.re_to_complex, Complex.ofReal_re]
  obtain ⟨hμi, harch⟩ := arch_re hf
  have hJi : Integrable (fun r : ℝ => ‖paperFT f r‖ ^ 2 * J a1 r) :=
    hφi.mul_bdd (continuous_J a1).aestronglyMeasurable
      (Eventually.of_forall fun r => by rw [Real.norm_eq_abs]; exact abs_J_le h0 r)
  have hμi' : Integrable (fun r : ℝ => ‖paperFT f r‖ ^ 2 * (2 * Real.pi * mu r)) := by
    refine (hμi.const_mul (2 * Real.pi)).congr (Eventually.of_forall fun r => ?_)
    simp only
    ring
  have hsum : (fun r : ℝ => ‖paperFT f r‖ ^ 2 * W a1 r)
      = fun r : ℝ => ‖paperFT f r‖ ^ 2 * (2 * Real.pi * mu r) + ‖paperFT f r‖ ^ 2 * J a1 r := by
    funext r
    unfold W
    ring
  refine ⟨by rw [hsum]; exact hμi'.add hJi, ?_⟩
  rw [weilQ_eq_fix]
  unfold weilQFix
  simp only [Complex.add_re, Complex.sub_re]
  rw [primeTerm_eq_zero hf.2 hlog, harch, pole_re_eq_J h0 (hf.mono ha1), Complex.zero_re, sub_zero,
    hsum, integral_add hμi' hJi]
  have hπ := Real.pi_pos
  have e : ∫ r : ℝ, ‖paperFT f r‖ ^ 2 * (2 * Real.pi * mu r)
      = 2 * Real.pi * ∫ r : ℝ, ‖paperFT f r‖ ^ 2 * mu r := by
    rw [← integral_const_mul]
    congr 1
    funext r
    ring
  rw [e]
  field_simp
  ring

/-! ### Affine cells -/

/-- An affine cell `(lo, hi, A, B)`. -/
abbrev Cell4 := ℝ × ℝ × ℝ × ℝ

/-- Weight `A + B|r|` on `lo ≤ |r| < hi` (both signs; `r = 0` may count twice). -/
def affFn (c : Cell4) (r : ℝ) : ℝ :=
  (Ico c.1 c.2.1).indicator (fun r => c.2.2.1 + c.2.2.2 * r) r
    + (Ioc (-c.2.1) (-c.1)).indicator (fun r => c.2.2.1 - c.2.2.2 * r) r

def affStep : List Cell4 → ℝ → ℝ
  | [], _ => 0
  | c :: cs, r => affFn c r + affStep cs r

/-- Consecutive cells from `x` to `R`. -/
def Chain4 : ℝ → List Cell4 → ℝ → Prop
  | x, [], R => x = R
  | x, c :: cs, R => c.1 = x ∧ Chain4 c.2.1 cs R

/-- A certified affine cell: `W ≥ s − (A + B|r|)` on `lo ≤ |r| < hi`, weight `≥ 0` there. -/
def AffOK (a s : ℝ) (c : Cell4) : Prop :=
  0 ≤ c.1 ∧ c.1 ≤ c.2.1 ∧ 0 ≤ c.2.2.1 + c.2.2.2 * c.1 ∧ 0 ≤ c.2.2.1 + c.2.2.2 * c.2.1 ∧
    ∀ r : ℝ, c.1 ≤ |r| → |r| < c.2.1 → s - (c.2.2.1 + c.2.2.2 * |r|) ≤ W a r

lemma aff_nonneg {A B lo hi r : ℝ} (h1 : 0 ≤ A + B * lo) (h2 : 0 ≤ A + B * hi) (hr1 : lo ≤ r)
    (hr2 : r ≤ hi) : 0 ≤ A + B * r := by
  rcases eq_or_lt_of_le (hr1.trans hr2) with h | h
  · have : r = lo := le_antisymm (h ▸ hr2) hr1
    rw [this]; exact h1
  · have key : (A + B * r) * (hi - lo) = (A + B * lo) * (hi - r) + (A + B * hi) * (r - lo) := by ring
    have : 0 ≤ (A + B * r) * (hi - lo) := by
      rw [key]
      have := mul_nonneg h1 (sub_nonneg.2 hr2)
      have := mul_nonneg h2 (sub_nonneg.2 hr1)
      linarith
    exact nonneg_of_mul_nonneg_left this (by linarith)

lemma affFn_nonneg {c : Cell4} (_h0 : 0 ≤ c.1) (h1 : 0 ≤ c.2.2.1 + c.2.2.2 * c.1)
    (h2 : 0 ≤ c.2.2.1 + c.2.2.2 * c.2.1) (r : ℝ) : 0 ≤ affFn c r := by
  unfold affFn
  refine add_nonneg (Set.indicator_nonneg (fun r hr => aff_nonneg h1 h2 hr.1 hr.2.le) r)
    (Set.indicator_nonneg (fun r hr => ?_) r)
  have := aff_nonneg h1 h2 (r := -r) (by linarith [hr.2]) (by linarith [hr.1])
  linarith

lemma integrable_affFn (c : Cell4) : Integrable (affFn c) := by
  unfold affFn
  refine Integrable.add ?_ ?_
  · exact ((Continuous.integrableOn_Icc (by fun_prop)).mono_set Ico_subset_Icc_self).integrable_indicator
      measurableSet_Ico
  · exact ((Continuous.integrableOn_Icc (by fun_prop)).mono_set Ioc_subset_Icc_self).integrable_indicator
      measurableSet_Ioc

lemma le_affFn {c : Cell4} (_h0 : 0 ≤ c.1) (h1 : 0 ≤ c.2.2.1 + c.2.2.2 * c.1)
    (h2 : 0 ≤ c.2.2.1 + c.2.2.2 * c.2.1) {r : ℝ} (hr1 : c.1 ≤ |r|) (hr2 : |r| < c.2.1) :
    c.2.2.1 + c.2.2.2 * |r| ≤ affFn c r := by
  unfold affFn
  have i1 : 0 ≤ (Ico c.1 c.2.1).indicator (fun r => c.2.2.1 + c.2.2.2 * r) r :=
    Set.indicator_nonneg (fun r hr => aff_nonneg h1 h2 hr.1 hr.2.le) r
  have i2 : 0 ≤ (Ioc (-c.2.1) (-c.1)).indicator (fun r => c.2.2.1 - c.2.2.2 * r) r := by
    refine Set.indicator_nonneg (fun r hr => ?_) r
    have := aff_nonneg h1 h2 (r := -r) (by linarith [hr.2]) (by linarith [hr.1])
    linarith
  rcases le_total 0 r with hr | hr
  · rw [abs_of_nonneg hr] at hr1 hr2 ⊢
    rw [Set.indicator_of_mem (show r ∈ Ico c.1 c.2.1 from ⟨hr1, hr2⟩)]
    linarith
  · rw [abs_of_nonpos hr] at hr1 hr2 ⊢
    rw [Set.indicator_of_mem (show r ∈ Ioc (-c.2.1) (-c.1) from ⟨by linarith, by linarith⟩) (f := fun r => c.2.2.1 - c.2.2.2 * r)]
    linarith

lemma affStep_nonneg : ∀ (cs : List Cell4), (∀ c ∈ cs, AffOK 0 0 c ∨ True) →
    (∀ c ∈ cs, 0 ≤ c.1 ∧ 0 ≤ c.2.2.1 + c.2.2.2 * c.1 ∧ 0 ≤ c.2.2.1 + c.2.2.2 * c.2.1) →
    ∀ r, 0 ≤ affStep cs r
  | [], _, _, _ => le_refl _
  | c :: cs, h', h, r => by
    simp only [List.mem_cons, forall_eq_or_imp] at h h'
    exact add_nonneg (affFn_nonneg h.1.1 h.1.2.1 h.1.2.2 r) (affStep_nonneg cs h'.2 h.2 r)

lemma integrable_affStep : ∀ (cs : List Cell4), Integrable (affStep cs)
  | [] => by
    have : affStep [] = fun _ => (0 : ℝ) := rfl
    rw [this]
    exact integrable_zero _ _ _
  | c :: cs => (integrable_affFn c).add (integrable_affStep cs)

lemma le_affStep : ∀ (cs : List Cell4),
    (∀ c ∈ cs, 0 ≤ c.1 ∧ 0 ≤ c.2.2.1 + c.2.2.2 * c.1 ∧ 0 ≤ c.2.2.1 + c.2.2.2 * c.2.1) →
    ∀ c ∈ cs, ∀ r : ℝ, c.1 ≤ |r| → |r| < c.2.1 → c.2.2.1 + c.2.2.2 * |r| ≤ affStep cs r
  | [], _, c, hc, _, _, _ => by simp at hc
  | d :: cs, h, c, hc, r, h1, h2 => by
    have h' := h
    simp only [List.mem_cons, forall_eq_or_imp] at h
    have hs : affStep (d :: cs) r = affFn d r + affStep cs r := rfl
    rw [hs]
    rcases List.mem_cons.1 hc with rfl | hc'
    · have := affStep_nonneg cs (fun _ _ => Or.inr trivial) h.2 r
      linarith [le_affFn h.1.1 h.1.2.1 h.1.2.2 h1 h2]
    · have := affFn_nonneg h.1.1 h.1.2.1 h.1.2.2 r
      linarith [le_affStep cs h.2 c hc' r h1 h2]

lemma chain4_cover : ∀ (cs : List Cell4) (x R : ℝ), Chain4 x cs R → ∀ y, x ≤ y → y < R →
    ∃ c ∈ cs, c.1 ≤ y ∧ y < c.2.1
  | [], x, R, h, y, h1, h2 => by
    simp only [Chain4] at h
    linarith
  | c :: cs, x, R, h, y, h1, h2 => by
    simp only [Chain4] at h
    rcases lt_or_ge y c.2.1 with hy | hy
    · exact ⟨c, List.mem_cons_self .., by rw [h.1]; exact h1, hy⟩
    · obtain ⟨d, hd, hd'⟩ := chain4_cover cs c.2.1 R h.2 y hy h2
      exact ⟨d, List.mem_cons_of_mem _ hd, hd'⟩

lemma chain4_append : ∀ (cs ds : List Cell4) (x y R : ℝ), Chain4 x cs y → Chain4 y ds R →
    Chain4 x (cs ++ ds) R
  | [], ds, x, y, R, h1, h2 => by
    simp only [Chain4] at h1
    rw [h1]
    exact h2
  | c :: cs, ds, x, y, R, h1, h2 => by
    simp only [Chain4] at h1
    simp only [List.cons_append, Chain4]
    exact ⟨h1.1, chain4_append cs ds c.2.1 y R h1.2 h2⟩

/-! ### Moments of affine cells -/

/-- `∫_lo^hi r^p (A + B r) dr`. -/
def gm (A B lo hi : ℝ) (p : ℕ) : ℝ :=
  A * ((hi ^ (p + 1) - lo ^ (p + 1)) / ((p : ℝ) + 1)) + B * ((hi ^ (p + 2) - lo ^ (p + 2)) / ((p : ℝ) + 2))

lemma integral_aff_pow (A B lo hi : ℝ) (p : ℕ) :
    ∫ r in lo..hi, r ^ p * (A + B * r) = gm A B lo hi p := by
  have e : (fun r : ℝ => r ^ p * (A + B * r)) = fun r => A * r ^ p + B * r ^ (p + 1) := by
    funext r; ring
  rw [e, intervalIntegral.integral_add, intervalIntegral.integral_const_mul,
    intervalIntegral.integral_const_mul, integral_pow, integral_pow]
  · unfold gm
    push_cast
    ring
  all_goals exact Continuous.intervalIntegrable (by fun_prop) _ _

lemma integral_poly_sq_aff (γ : ℕ → ℝ) (K : ℕ) (A B lo hi : ℝ) :
    ∫ r in lo..hi, (∑ m ∈ range K, γ m * r ^ m) ^ 2 * (A + B * r)
      = ∑ m ∈ range K, ∑ k ∈ range K, γ m * γ k * gm A B lo hi (m + k) := by
  have e : ∀ r : ℝ, (∑ m ∈ range K, γ m * r ^ m) ^ 2 * (A + B * r)
      = ∑ m ∈ range K, ∑ k ∈ range K, γ m * γ k * (r ^ (m + k) * (A + B * r)) := by
    intro r
    rw [sq, Finset.sum_mul_sum, Finset.sum_mul]
    refine Finset.sum_congr rfl fun m _ => ?_
    rw [Finset.sum_mul]
    refine Finset.sum_congr rfl fun k _ => ?_
    ring
  simp_rw [e]
  rw [intervalIntegral.integral_finsetSum (fun m _ => ?_)]
  · refine Finset.sum_congr rfl fun m _ => ?_
    rw [intervalIntegral.integral_finsetSum (fun k _ => ?_)]
    · refine Finset.sum_congr rfl fun k _ => ?_
      rw [intervalIntegral.integral_const_mul, integral_aff_pow]
    · exact Continuous.intervalIntegrable (by fun_prop) _ _
  · exact Continuous.intervalIntegrable (by fun_prop) _ _

lemma integral_phiMaj_aff (f : ℝ → ℂ) (K : ℕ) (ε Cβ A B lo hi : ℝ) :
    ∫ r in lo..hi, phiMaj f K ε Cβ r * (A + B * r)
      = (1 + ε) * ∑ m ∈ range K, ∑ k ∈ range K, ((momB f m).re * (momB f k).re
          + (momB f m).im * (momB f k).im) * gm A B lo hi (m + k)
        + (1 + 1 / ε) * Cβ * gm A B lo hi (2 * K) := by
  have e : (fun r : ℝ => phiMaj f K ε Cβ r * (A + B * r))
      = fun r => (1 + ε) * ((∑ m ∈ range K, (momB f m).re * r ^ m) ^ 2 * (A + B * r))
        + (1 + ε) * ((∑ m ∈ range K, (momB f m).im * r ^ m) ^ 2 * (A + B * r))
        + ((1 + 1 / ε) * Cβ) * (r ^ (2 * K) * (A + B * r)) := by
    funext r
    unfold phiMaj
    ring
  rw [e, intervalIntegral.integral_add, intervalIntegral.integral_add, intervalIntegral.integral_const_mul,
    intervalIntegral.integral_const_mul, intervalIntegral.integral_const_mul, integral_poly_sq_aff,
    integral_poly_sq_aff, integral_aff_pow]
  · simp only [add_mul, Finset.sum_add_distrib, mul_add]
  all_goals exact Continuous.intervalIntegrable (by fun_prop) _ _

/-- Moment of an affine cell (both signs). -/
def affMom (c : Cell4) (p : ℕ) : ℝ :=
  gm c.2.2.1 c.2.2.2 c.1 c.2.1 p + gm c.2.2.1 (-c.2.2.2) (-c.2.1) (-c.1) p

def affSum : List Cell4 → ℕ → ℝ
  | [], _ => 0
  | c :: cs, p => affMom c p + affSum cs p

lemma affMom_even_nonneg {c : Cell4} (_h0 : 0 ≤ c.1) (hlh : c.1 ≤ c.2.1) (h1 : 0 ≤ c.2.2.1 + c.2.2.2 * c.1)
    (h2 : 0 ≤ c.2.2.1 + c.2.2.2 * c.2.1) (K : ℕ) : 0 ≤ affMom c (2 * K) := by
  unfold affMom
  rw [← integral_aff_pow, ← integral_aff_pow]
  have i1 : 0 ≤ ∫ r in c.1..c.2.1, r ^ (2 * K) * (c.2.2.1 + c.2.2.2 * r) :=
    intervalIntegral.integral_nonneg hlh fun r hr =>
      mul_nonneg (by rw [pow_mul]; positivity) (aff_nonneg h1 h2 hr.1 hr.2)
  have i2 : 0 ≤ ∫ r in (-c.2.1)..(-c.1), r ^ (2 * K) * (c.2.2.1 + -c.2.2.2 * r) :=
    intervalIntegral.integral_nonneg (by linarith) fun r hr => by
      refine mul_nonneg (by rw [pow_mul]; positivity) ?_
      have := aff_nonneg h1 h2 (r := -r) (by linarith [hr.2]) (by linarith [hr.1])
      linarith
  linarith

lemma affSum_even_nonneg : ∀ cs : List Cell4,
    (∀ c ∈ cs, 0 ≤ c.1 ∧ c.1 ≤ c.2.1 ∧ 0 ≤ c.2.2.1 + c.2.2.2 * c.1 ∧ 0 ≤ c.2.2.1 + c.2.2.2 * c.2.1) →
    ∀ K : ℕ, 0 ≤ affSum cs (2 * K)
  | [], _, _ => by simp [affSum]
  | c :: cs, h, K => by
    simp only [List.mem_cons, forall_eq_or_imp] at h
    simp only [affSum]
    have h1 := affMom_even_nonneg h.1.1 h.1.2.1 h.1.2.2.1 h.1.2.2.2 K
    have h2 := affSum_even_nonneg cs h.2 K
    linarith

/-- One affine cell: `∫ φ · affFn ≤` (moment bound). -/
lemma aff_cell_bound {a : ℝ} (ha : 0 < a) {f : ℝ → ℂ} (hf : TestFn a f)
    (hφi : Integrable (fun r : ℝ => ‖paperFT f r‖ ^ 2)) (K : ℕ) {ε : ℝ} (hε : 0 < ε)
    {c : Cell4} (_h0 : 0 ≤ c.1) (h1 : c.1 ≤ c.2.1) (hA : 0 ≤ c.2.2.1 + c.2.2.2 * c.1)
    (hB : 0 ≤ c.2.2.1 + c.2.2.2 * c.2.1)
    (hK : (c.2.2.1 = 0 ∧ c.2.2.2 = 0) ∨ c.2.1 * a ≤ (K + 1) / 2) :
    Integrable (fun r : ℝ => ‖paperFT f r‖ ^ 2 * affFn c r) ∧
    ∫ r : ℝ, ‖paperFT f r‖ ^ 2 * affFn c r ≤ (1 + ε) * ∑ m ∈ range K, ∑ k ∈ range K,
        ((momB f m).re * (momB f k).re + (momB f m).im * (momB f k).im) * affMom c (m + k)
      + (1 + 1 / ε) * (4 * a ^ (2 * K) * (∫ u, ‖f u‖) ^ 2 / (K.factorial : ℝ) ^ 2)
        * affMom c (2 * K) := by
  obtain ⟨lo, hi, A, B⟩ := c
  simp only at h1 hA hB hK ⊢
  set φ := fun r : ℝ => ‖paperFT f r‖ ^ 2 with hφ
  have hφ0 : ∀ r, 0 ≤ φ r := fun r => by positivity
  have e : (fun r => φ r * affFn (lo, hi, A, B) r)
      = fun r => (Ico lo hi).indicator (fun r => φ r * (A + B * r)) r
        + (Ioc (-hi) (-lo)).indicator (fun r => φ r * (A - B * r)) r := by
    funext r
    unfold affFn
    simp only [Set.indicator]
    split_ifs <;> ring
  have hc1 : Continuous fun r : ℝ => A + B * r := by fun_prop
  have hc2 : Continuous fun r : ℝ => A - B * r := by fun_prop
  have i1 : Integrable fun r => (Ico lo hi).indicator (fun r => φ r * (A + B * r)) r := by
    refine IntegrableOn.integrable_indicator ?_ measurableSet_Ico
    refine (hφi.integrableOn.mul_continuousOn (hc1.continuousOn) isCompact_Icc).mono_set Ico_subset_Icc_self
  have i2 : Integrable fun r => (Ioc (-hi) (-lo)).indicator (fun r => φ r * (A - B * r)) r := by
    refine IntegrableOn.integrable_indicator ?_ measurableSet_Ioc
    exact (hφi.integrableOn.mul_continuousOn (hc2.continuousOn) isCompact_Icc).mono_set Ioc_subset_Icc_self
  refine ⟨by rw [e]; exact i1.add i2, ?_⟩
  rw [e, integral_add i1 i2, integral_indicator measurableSet_Ico,
    integral_indicator measurableSet_Ioc]
  rcases hK with ⟨hA0, hB0⟩ | hK
  · subst hA0 hB0
    simp [affMom, gm]
  set Cβ := 4 * a ^ (2 * K) * (∫ u, ‖f u‖) ^ 2 / (K.factorial : ℝ) ^ 2
  have hle : ∀ r : ℝ, |r| ≤ hi → φ r ≤ phiMaj f K ε Cβ r := by
    intro r hr
    have hr' : |r| * a ≤ (K + 1) / 2 :=
      le_trans (mul_le_mul_of_nonneg_right hr ha.le) hK
    exact phi_le_taylor ha hf K hε hr'
  have hP := continuous_phiMaj f K ε Cβ
  have j1 : ∫ r in Ico lo hi, φ r * (A + B * r) ≤ ∫ r in lo..hi, phiMaj f K ε Cβ r * (A + B * r) := by
    rw [integral_Ico_eq_integral_Ioc, ← intervalIntegral.integral_of_le h1]
    refine intervalIntegral.integral_mono_on h1
      ((intervalIntegrable_iff_integrableOn_Ioc_of_le h1).2
        ((hφi.integrableOn.mul_continuousOn (hc1.continuousOn) isCompact_Icc).mono_set Ioc_subset_Icc_self))
      ((hP.mul hc1).intervalIntegrable _ _) fun r hr => ?_
    have hw := aff_nonneg hA hB hr.1 hr.2
    refine mul_le_mul_of_nonneg_right (hle r ?_) hw
    rw [abs_of_nonneg (by linarith [hr.1])]
    exact hr.2
  have j2 : ∫ r in Ioc (-hi) (-lo), φ r * (A - B * r)
      ≤ ∫ r in (-hi)..(-lo), phiMaj f K ε Cβ r * (A + -B * r) := by
    rw [← intervalIntegral.integral_of_le (by linarith)]
    refine intervalIntegral.integral_mono_on (by linarith)
      ((intervalIntegrable_iff_integrableOn_Ioc_of_le (by linarith)).2
        ((hφi.integrableOn.mul_continuousOn (hc2.continuousOn) isCompact_Icc).mono_set Ioc_subset_Icc_self))
      ((hP.mul (by fun_prop : Continuous fun r : ℝ => A + -B * r)).intervalIntegrable _ _) fun r hr => ?_
    have hw : 0 ≤ A + -B * r := by
      have := aff_nonneg hA hB (r := -r) (by linarith [hr.2]) (by linarith [hr.1])
      linarith
    have e3 : A - B * r = A + -B * r := by ring
    rw [e3]
    refine mul_le_mul_of_nonneg_right (hle r ?_) hw
    rw [abs_of_nonpos (by linarith [hr.2])]
    linarith [hr.1]
  refine (add_le_add j1 j2).trans (le_of_eq ?_)
  rw [integral_phiMaj_aff, integral_phiMaj_aff]
  unfold affMom
  simp only [mul_add, Finset.sum_add_distrib]
  ring

/-- Bound for the whole affine step function. -/
lemma aff_step_bound {a : ℝ} (ha : 0 < a) {f : ℝ → ℂ} (hf : TestFn a f)
    (hφi : Integrable (fun r : ℝ => ‖paperFT f r‖ ^ 2)) (K : ℕ) {ε : ℝ} (hε : 0 < ε) :
    ∀ cs : List Cell4,
    (∀ c ∈ cs, 0 ≤ c.1 ∧ c.1 ≤ c.2.1 ∧ 0 ≤ c.2.2.1 + c.2.2.2 * c.1 ∧ 0 ≤ c.2.2.1 + c.2.2.2 * c.2.1) →
    (∀ c ∈ cs, (c.2.2.1 = 0 ∧ c.2.2.2 = 0) ∨ c.2.1 * a ≤ (K + 1) / 2) →
    Integrable (fun r : ℝ => ‖paperFT f r‖ ^ 2 * affStep cs r) ∧
    ∫ r : ℝ, ‖paperFT f r‖ ^ 2 * affStep cs r ≤ (1 + ε) * ∑ m ∈ range K, ∑ k ∈ range K,
        ((momB f m).re * (momB f k).re + (momB f m).im * (momB f k).im) * affSum cs (m + k)
      + (1 + 1 / ε) * (4 * a ^ (2 * K) * (∫ u, ‖f u‖) ^ 2 / (K.factorial : ℝ) ^ 2)
        * affSum cs (2 * K)
  | [], _, _ => by
    have e : (fun r : ℝ => ‖paperFT f r‖ ^ 2 * affStep [] r) = fun _ => (0 : ℝ) := by
      funext r; simp [affStep]
    rw [e]
    simp [affSum]
  | c :: cs, h, hK => by
    simp only [List.mem_cons, forall_eq_or_imp] at h hK
    obtain ⟨hc1, hc2⟩ := aff_cell_bound ha hf hφi K hε h.1.1 h.1.2.1 h.1.2.2.1 h.1.2.2.2 hK.1
    obtain ⟨hs1, hs2⟩ := aff_step_bound ha hf hφi K hε cs h.2 hK.2
    have e : (fun r : ℝ => ‖paperFT f r‖ ^ 2 * affStep (c :: cs) r)
        = fun r : ℝ => ‖paperFT f r‖ ^ 2 * affFn c r + ‖paperFT f r‖ ^ 2 * affStep cs r := by
      funext r
      simp only [affStep]
      ring
    refine ⟨by rw [e]; exact hc1.add hs1, ?_⟩
    rw [e, integral_add hc1 hs1]
    refine (add_le_add hc2 hs2).trans (le_of_eq ?_)
    simp only [affSum, mul_add, Finset.sum_add_distrib, Finset.mul_sum]
    ring_nf

/-! ### A certified affine cell from end values -/

lemma two_pi_mu_abs (r : ℝ) : 2 * Real.pi * mu |r| = 2 * Real.pi * mu r := by
  have e : (|r| / 2) ^ 2 = (r / 2) ^ 2 := by rw [div_pow, div_pow, sq_abs]
  rw [two_pi_mu_series, two_pi_mu_series, e]

/-- **Affine cell from rational end data.** `mL, mH` bound `2πμ` and `jL, jH` bound `J` at the ends;
`δ = (mH − mL)(hi − lo)/(4(hi + lo)) + E_j (hi − lo)²` is the second-order loss. -/
lemma affOK_of {a s lo hi A B mL mH jL jH Ej : ℝ} (ha : 0 ≤ a) (h0 : 0 ≤ lo) (hlh : lo < hi)
    (hmL : mL ≤ 2 * Real.pi * mu lo) (hmH : mH ≤ 2 * Real.pi * mu hi) (hjL : jL ≤ J a lo)
    (hjH : jH ≤ J a hi) (hE : 4 / 3 * a ^ 3 * Real.cosh a ≤ Ej) (hΔ : mL ≤ mH)
    (hn1 : 0 ≤ A + B * lo) (hn2 : 0 ≤ A + B * hi)
    (hc1 : s - (A + B * lo) ≤ mL + jL - ((mH - mL) * ((hi - lo) / (4 * (hi + lo))) + Ej * (hi - lo) ^ 2))
    (hc2 : s - (A + B * hi) ≤ mH + jH - ((mH - mL) * ((hi - lo) / (4 * (hi + lo))) + Ej * (hi - lo) ^ 2)) :
    AffOK a s (lo, hi, A, B) := by
  refine ⟨h0, hlh.le, hn1, hn2, fun r h1 h2 => ?_⟩
  simp only at h1 h2 ⊢
  set x := |r| with hxdef
  have hJr : J a r = J a x := by
    rcases le_total 0 r with h | h
    · rw [hxdef, abs_of_nonneg h]
    · rw [hxdef, abs_of_nonpos h, J_even]
  have hW : W a r = 2 * Real.pi * mu x + J a x := by
    unfold W
    rw [hJr, hxdef, two_pi_mu_abs]
  rw [hW]
  have hw : 0 < hi - lo := by linarith
  have hD : 0 < hi ^ 2 - lo ^ 2 := by nlinarith
  set l := (x - lo) / (hi - lo) with hl
  have hl0 : 0 ≤ l := div_nonneg (by linarith) hw.le
  have hl1 : l ≤ 1 := by rw [hl, div_le_one hw]; linarith
  have hxl : x = (1 - l) * lo + l * hi := by
    rw [hl]; field_simp; ring
  set θ := (x ^ 2 - lo ^ 2) / (hi ^ 2 - lo ^ 2) with hθ
  have hθ0 : 0 ≤ θ := div_nonneg (by nlinarith) hD.le
  have hθ1 : θ ≤ 1 := by rw [hθ, div_le_one hD]; nlinarith
  have hx2 : (x / 2) ^ 2 = (1 - θ) * (lo / 2) ^ 2 + θ * (hi / 2) ^ 2 := by
    rw [hθ]; field_simp; ring
  -- μ part
  have hmu := mu_chord hθ0 hθ1 hx2
  have k1 : (1 - θ) * mL + θ * mH ≤ 2 * Real.pi * mu x := by
    have := mul_le_mul_of_nonneg_left hmL (by linarith : (0 : ℝ) ≤ 1 - θ)
    have := mul_le_mul_of_nonneg_left hmH hθ0
    linarith
  set q := (hi - lo) / (4 * (hi + lo)) with hq
  have hθl : l - q ≤ θ := by
    have hs : hi + lo ≠ 0 := by linarith
    have hw' : hi - lo ≠ 0 := by linarith
    have hD' : hi ^ 2 - lo ^ 2 ≠ 0 := by linarith
    have e : (θ - (l - q)) * (hi ^ 2 - lo ^ 2) = (x - (hi + lo) / 2) ^ 2 := by
      rw [hθ, hl, hq]
      field_simp
      ring
    have : 0 ≤ (θ - (l - q)) * (hi ^ 2 - lo ^ 2) := by rw [e]; positivity
    have := nonneg_of_mul_nonneg_left this hD
    linarith
  have k2 : (1 - l) * mL + l * mH - (mH - mL) * q ≤ 2 * Real.pi * mu x := by
    have h3 := mul_le_mul_of_nonneg_left hθl (sub_nonneg.2 hΔ)
    have e1 : (1 - θ) * mL + θ * mH = mL + (mH - mL) * θ := by ring
    have e2 : (1 - l) * mL + l * mH - (mH - mL) * q = mL + (mH - mL) * (l - q) := by ring
    rw [e1] at k1
    rw [e2]
    linarith only [k1, h3]
  -- J part
  have hJ := J_chord (lo := lo) (hi := hi) ha hl0 hl1
  rw [← hxl] at hJ
  have hll : l * (1 - l) ≤ 1 / 4 := by nlinarith [sq_nonneg (l - 1 / 2)]
  have hc0 : 0 ≤ a ^ 3 * Real.cosh a := by have := Real.cosh_pos a; positivity
  have k3 : 16 / 3 * a ^ 3 * Real.cosh a * (l * (1 - l) * (hi - lo) ^ 2) ≤ Ej * (hi - lo) ^ 2 := by
    have hw2 : 0 ≤ (hi - lo) ^ 2 := sq_nonneg _
    have h16 : 0 ≤ 16 / 3 * a ^ 3 * Real.cosh a := by positivity
    have t1 := mul_le_mul_of_nonneg_left hll h16
    have t2 : 16 / 3 * a ^ 3 * Real.cosh a * (l * (1 - l)) ≤ Ej := by linarith only [t1, hE]
    have t3 := mul_le_mul_of_nonneg_right t2 hw2
    have e4 : 16 / 3 * a ^ 3 * Real.cosh a * (l * (1 - l) * (hi - lo) ^ 2)
        = 16 / 3 * a ^ 3 * Real.cosh a * (l * (1 - l)) * (hi - lo) ^ 2 := by ring
    rw [e4]
    exact t3
  have k4 : (1 - l) * jL + l * jH ≤ (1 - l) * J a lo + l * J a hi := by
    have := mul_le_mul_of_nonneg_left hjL (by linarith : (0 : ℝ) ≤ 1 - l)
    have := mul_le_mul_of_nonneg_left hjH hl0
    linarith
  -- combine
  have k5 := mul_le_mul_of_nonneg_left hc1 (by linarith : (0 : ℝ) ≤ 1 - l)
  have k6 := mul_le_mul_of_nonneg_left hc2 hl0
  have e : s - (A + B * x) = (1 - l) * (s - (A + B * lo)) + l * (s - (A + B * hi)) := by
    rw [hxl]; ring
  have e3 : (1 - l) * (mL + jL - ((mH - mL) * q + Ej * (hi - lo) ^ 2))
      + l * (mH + jH - ((mH - mL) * q + Ej * (hi - lo) ^ 2))
      = ((1 - l) * mL + l * mH - (mH - mL) * q) + ((1 - l) * jL + l * jH) - Ej * (hi - lo) ^ 2 := by
    ring
  rw [e]
  linarith only [k2, hJ, k3, k4, k5, k6, e3]

/-! ### The certificate theorem -/

/-- **Affine moment certificate at `a₁` ⇒ Weil positivity on every window `L < log 2`, `L ≤ 2a₁`.** -/
theorem weilPos_of_affine {a1 s R ε α β : ℝ} (ha1 : 0 < a1)
    (cs : List Cell4) (hchain : Chain4 0 cs R) (hok : ∀ c ∈ cs, AffOK a1 s c)
    (htail : ∀ r : ℝ, R ≤ |r| → s ≤ W a1 r) (K : ℕ)
    (hK : ∀ c ∈ cs, (c.2.2.1 = 0 ∧ c.2.2.2 = 0) ∨ c.2.1 * a1 ≤ (K + 1) / 2) (Λ : ℕ → ℕ → ℝ)
    (hpsd : ∀ x : ℕ → ℝ, ∑ m ∈ range K, ∑ k ∈ range K, x m * x k * affSum cs (m + k)
      ≤ α * (2 * ∑ m ∈ range K, (∑ i ∈ range K, Λ m i * x i) * x m
        - ∑ m ∈ range K, ∑ k ∈ range K,
            ((∑ i ∈ range K, Λ m i * x i) * (∑ i ∈ range K, Λ k i * x i)) * Ghat a1 m k))
    (hε : 0 < ε) (hα : 0 ≤ α)
    (hβ : 8 * a1 ^ (2 * K + 1) * affSum cs (2 * K) / (K.factorial : ℝ) ^ 2 ≤ β)
    (hnum : (1 + ε) * α + (1 + 1 / ε) * β ≤ 2 * Real.pi * s)
    {L : ℝ} (hL : L < Real.log 2) (hL1 : L ≤ 2 * a1) : WeilPos L := by
  have key : ∀ L : ℝ, 0 ≤ L → L < Real.log 2 → L ≤ 2 * a1 → WeilPos L := by
    intro L hL0 hL hL1 f hf
    set a := L / 2 with hadef
    have ha0 : 0 ≤ a := by rw [hadef]; linarith
    have haa : a ≤ a1 := by rw [hadef]; linarith
    have hlog : 2 * a < Real.log 2 := by rw [hadef]; linarith
    have hf1 : TestFn a1 f := hf.mono haa
    obtain ⟨hφw, hre⟩ := weilQ_re_eq_W' ha0 haa hlog hf
    have hk := hf1.autocorr_contDiff
    have hkc := hf1.autocorr_hasCompactSupport
    have hFk := integrable_fourier_of_contDiff_two hk hkc
    set N := ∫ u, ‖f u‖ ^ 2 with hN
    have hN0 : 0 ≤ N := integral_nonneg fun u => by positivity
    have hφi : Integrable (fun r : ℝ => ‖paperFT f r‖ ^ 2) := by
      refine (integrable_paperFT_ofReal hFk).re.congr (Eventually.of_forall fun r => ?_)
      simp only [paperFT_autocorr_real hf1, RCLike.re_to_complex, Complex.ofReal_re]
    have hcs : ∀ c ∈ cs, 0 ≤ c.1 ∧ c.1 ≤ c.2.1 ∧ 0 ≤ c.2.2.1 + c.2.2.2 * c.1
        ∧ 0 ≤ c.2.2.1 + c.2.2.2 * c.2.1 :=
      fun c hc => ⟨(hok c hc).1, (hok c hc).2.1, (hok c hc).2.2.1, (hok c hc).2.2.2.1⟩
    have hcs' : ∀ c ∈ cs, 0 ≤ c.1 ∧ 0 ≤ c.2.2.1 + c.2.2.2 * c.1 ∧ 0 ≤ c.2.2.1 + c.2.2.2 * c.2.1 :=
      fun c hc => ⟨(hcs c hc).1, (hcs c hc).2.2.1, (hcs c hc).2.2.2⟩
    have hw : ∀ r : ℝ, s - affStep cs r ≤ W a1 r := by
      intro r
      rcases le_or_gt R |r| with hR | hR
      · have := affStep_nonneg cs (fun _ _ => Or.inr trivial) hcs' r
        linarith [htail r hR]
      · obtain ⟨c, hc, h1, h2⟩ := chain4_cover cs 0 R hchain |r| (abs_nonneg r) hR
        have := le_affStep cs hcs' c hc r h1 h2
        linarith [(hok c hc).2.2.2.2 r h1 h2]
    obtain ⟨hsi, hsb⟩ := aff_step_bound ha1 hf1 hφi K hε cs hcs hK
    have hlow : s * (2 * Real.pi * N) - ∫ r : ℝ, ‖paperFT f r‖ ^ 2 * affStep cs r
        ≤ ∫ r : ℝ, ‖paperFT f r‖ ^ 2 * W a1 r := by
      have hint : Integrable (fun r : ℝ => s * ‖paperFT f r‖ ^ 2 - ‖paperFT f r‖ ^ 2 * affStep cs r) :=
        (hφi.const_mul s).sub hsi
      have h := integral_mono hint hφw (fun r => by
        have h1 := mul_le_mul_of_nonneg_left (hw r) (by positivity : (0 : ℝ) ≤ ‖paperFT f r‖ ^ 2)
        simp only
        nlinarith)
      rw [integral_sub (hφi.const_mul s) hsi, integral_const_mul, plancherel hf1] at h
      exact h
    set x : ℕ → ℝ := fun m => (momB f m).re
    set y : ℕ → ℝ := fun m => (momB f m).im
    have hX : ∑ m ∈ range K, ∑ k ∈ range K, ((momB f m).re * (momB f k).re
        + (momB f m).im * (momB f k).im) * affSum cs (m + k) ≤ α * N := by
      have hx := hpsd x
      have hy := hpsd y
      have hb := bessel ha1 hf1 K (fun m => ∑ i ∈ range K, Λ m i * x i)
        (fun m => ∑ i ∈ range K, Λ m i * y i)
      have e1 : ∑ m ∈ range K, ∑ k ∈ range K, ((momB f m).re * (momB f k).re
          + (momB f m).im * (momB f k).im) * affSum cs (m + k)
          = ∑ m ∈ range K, ∑ k ∈ range K, x m * x k * affSum cs (m + k)
            + ∑ m ∈ range K, ∑ k ∈ range K, y m * y k * affSum cs (m + k) := by
        rw [← Finset.sum_add_distrib]
        refine Finset.sum_congr rfl fun m _ => ?_
        rw [← Finset.sum_add_distrib]
        refine Finset.sum_congr rfl fun k _ => ?_
        ring
      have e2 : 2 * ∑ m ∈ range K, ((∑ i ∈ range K, Λ m i * x i) * (momB f m).re
          + (∑ i ∈ range K, Λ m i * y i) * (momB f m).im)
          - ∑ m ∈ range K, ∑ k ∈ range K, ((∑ i ∈ range K, Λ m i * x i) * (∑ i ∈ range K, Λ k i * x i)
            + (∑ i ∈ range K, Λ m i * y i) * (∑ i ∈ range K, Λ k i * y i)) * Ghat a1 m k
          = (2 * ∑ m ∈ range K, (∑ i ∈ range K, Λ m i * x i) * x m
            - ∑ m ∈ range K, ∑ k ∈ range K,
              ((∑ i ∈ range K, Λ m i * x i) * (∑ i ∈ range K, Λ k i * x i)) * Ghat a1 m k)
          + (2 * ∑ m ∈ range K, (∑ i ∈ range K, Λ m i * y i) * y m
            - ∑ m ∈ range K, ∑ k ∈ range K,
              ((∑ i ∈ range K, Λ m i * y i) * (∑ i ∈ range K, Λ k i * y i)) * Ghat a1 m k) := by
        simp only [x, y, add_mul, Finset.sum_add_distrib, mul_add]
        ring
      rw [e1]
      rw [e2] at hb
      nlinarith
    have hA : (∫ u, ‖f u‖) ^ 2 ≤ 2 * a1 * N := sq_integral_norm_le ha1 hf1
    have hM0 : 0 ≤ affSum cs (2 * K) := affSum_even_nonneg cs hcs K
    set A := ∫ u, ‖f u‖
    have hβN : (1 + 1 / ε) * (4 * a1 ^ (2 * K) * A ^ 2 / (K.factorial : ℝ) ^ 2) * affSum cs (2 * K)
        ≤ (1 + 1 / ε) * β * N := by
      have hc : 0 ≤ 1 + 1 / ε := by positivity
      rw [mul_assoc, mul_assoc (1 + 1 / ε)]
      refine mul_le_mul_of_nonneg_left ?_ hc
      have h1 : 4 * a1 ^ (2 * K) * A ^ 2 / (K.factorial : ℝ) ^ 2 * affSum cs (2 * K)
          ≤ 4 * a1 ^ (2 * K) * (2 * a1 * N) / (K.factorial : ℝ) ^ 2 * affSum cs (2 * K) := by
        gcongr
      have h2 : 4 * a1 ^ (2 * K) * (2 * a1 * N) / (K.factorial : ℝ) ^ 2 * affSum cs (2 * K)
          = (8 * a1 ^ (2 * K + 1) * affSum cs (2 * K) / (K.factorial : ℝ) ^ 2) * N := by
        ring
      rw [h2] at h1
      exact h1.trans (mul_le_mul_of_nonneg_right hβ hN0)
    have hstep : ∫ r : ℝ, ‖paperFT f r‖ ^ 2 * affStep cs r ≤ (1 + ε) * α * N + (1 + 1 / ε) * β * N := by
      have h1 := mul_le_mul_of_nonneg_left hX (by linarith : (0 : ℝ) ≤ 1 + ε)
      linarith
    have hfin : 0 ≤ ∫ r : ℝ, ‖paperFT f r‖ ^ 2 * W a1 r := by
      have h3 : ((1 + ε) * α + (1 + 1 / ε) * β) * N ≤ (2 * Real.pi * s) * N :=
        mul_le_mul_of_nonneg_right hnum hN0
      nlinarith
    rw [hre]
    have hπ := Real.pi_pos
    positivity
  rcases le_or_gt 0 L with h | h
  · exact key L h hL hL1
  · exact WeilPos.mono h.le (key 0 le_rfl (Real.log_pos (by norm_num)) (by linarith))

end Fayed.Yoshida
