import Riemann.Deep.FirstZero.MinModel
import Riemann.Deep.WeilEF.Z23.RvM.LocalCount
import Riemann.Deep.WeilEF.Z23.Analytic.RectangleLogDeriv
import Riemann.Deep.WeilEF.Z23.ZetaReflect

/-!
# Minimality of the first zero, part 2: the argument principle on [σL, 2] × [−T, T]

winding: for the rectangle R = [σL, 2] × [−T, T] (0 < σL < 1) with ζ ≠ 0 on the upper half of
its border, Σ_{zeros ρ ∈ R} ord_ρ ζ − 1 = (1/π)·(arg change of ζ along 2 → 2+iT → σL+iT → σL),
written as −Im ∫_top + Re ∫_right − Re ∫_left of ζ'/ζ. (Weighted argument principle of zeta23
with the simple pole at 1, folded by ζ(s̄) = conj ζ(s).)

TopOK / LeftOK package nonvanishing plus an upper bound for one piece; topOK_of_slit,
leftOK_of_slit compute a piece from an anchor c with Re(ζ·c̄) > 0 along it (log of ζ·c̄).
-/

open Complex Set Filter Topology MeasureTheory

namespace Fayed.FirstZero.Min
open Fayed.FirstZero

/-! ## Conjugation and continuity -/

lemma deriv_zeta_conj {s : ℂ} (hs : s ≠ 1) :
    deriv riemannZeta ((starRingEnd ℂ) s) = (starRingEnd ℂ) (deriv riemannZeta s) := by
  have h2 := (differentiableAt_riemannZeta hs).hasDerivAt.conj_conj
  have hf : (⇑(starRingEnd ℂ) ∘ riemannZeta ∘ ⇑(starRingEnd ℂ)) = riemannZeta := by
    funext z; simp [Function.comp, riemannZeta_conj]
  rw [hf] at h2
  exact h2.deriv

lemma logDeriv_zeta_conj {s : ℂ} (hs : s ≠ 1) :
    logDeriv riemannZeta ((starRingEnd ℂ) s) = (starRingEnd ℂ) (logDeriv riemannZeta s) := by
  simp only [logDeriv_apply, deriv_zeta_conj hs, riemannZeta_conj, map_div₀]

lemma continuousAt_logDeriv_zeta {s : ℂ} (hs : s ≠ 1) (h0 : riemannZeta s ≠ 0) :
    ContinuousAt (logDeriv riemannZeta) s := by
  have han : AnalyticAt ℂ riemannZeta s := Zeta23.riemannZeta_analyticOnNhd_compl_one s hs
  show ContinuousAt (fun x => deriv riemannZeta x / riemannZeta x) s
  exact han.deriv.continuousAt.div han.continuousAt h0

lemma ne_one_of_im {x T : ℝ} (hT : T ≠ 0) : (x : ℂ) + T * I ≠ 1 := by
  intro h; have := congrArg Complex.im h; simp at this; exact hT this

lemma ne_one_of_re {σ y : ℝ} (hσ : σ ≠ 1) : (σ : ℂ) + y * I ≠ 1 := by
  intro h; have := congrArg Complex.re h; simp at this; exact hσ this

lemma intInt_top {T a b : ℝ} (hT : T ≠ 0) (hab : a ≤ b)
    (h0 : ∀ x ∈ Icc a b, riemannZeta (x + T * I) ≠ 0) :
    IntervalIntegrable (fun x : ℝ => logDeriv riemannZeta (x + T * I)) volume a b := by
  apply ContinuousOn.intervalIntegrable
  rw [uIcc_of_le hab]
  intro x hx
  have h := continuousAt_logDeriv_zeta (ne_one_of_im hT) (h0 x hx)
  exact (ContinuousAt.comp (g := logDeriv riemannZeta) (f := fun x : ℝ => (x : ℂ) + T * I) (x := x)
    h (by fun_prop)).continuousWithinAt

lemma intInt_vert {σ a b : ℝ} (hσ : σ ≠ 1) (hab : a ≤ b)
    (h0 : ∀ y ∈ Icc a b, riemannZeta (σ + y * I) ≠ 0) :
    IntervalIntegrable (fun y : ℝ => logDeriv riemannZeta (σ + y * I)) volume a b := by
  apply ContinuousOn.intervalIntegrable
  rw [uIcc_of_le hab]
  intro y hy
  have h := continuousAt_logDeriv_zeta (ne_one_of_re hσ) (h0 y hy)
  exact (ContinuousAt.comp (g := logDeriv riemannZeta) (f := fun y : ℝ => (σ : ℂ) + y * I) (x := y)
    h (by fun_prop)).continuousWithinAt

/-! ## One piece: ∫ ζ'/ζ = log(ζ c̄)(end) − log(ζ c̄)(start) -/

lemma top_ftc {T a b : ℝ} (hT : T ≠ 0) (hab : a ≤ b) {c : ℂ}
    (hslit : ∀ x ∈ Icc a b, 0 < (riemannZeta (x + T * I) * (starRingEnd ℂ) c).re) :
    ∫ x in a..b, logDeriv riemannZeta (x + T * I) =
      log (riemannZeta (b + T * I) * (starRingEnd ℂ) c)
        - log (riemannZeta (a + T * I) * (starRingEnd ℂ) c) := by
  have hnz : ∀ x ∈ Icc a b, riemannZeta (x + T * I) ≠ 0 ∧ (starRingEnd ℂ) c ≠ 0 := by
    intro x hx
    have h := hslit x hx
    constructor <;> intro h0 <;> rw [h0] at h <;> simp at h
  refine intervalIntegral.integral_eq_sub_of_hasDerivAt
    (f := fun x : ℝ => log (riemannZeta (x + T * I) * (starRingEnd ℂ) c)) ?_
    (intInt_top hT hab fun x hx => (hnz x hx).1)
  intro x hx
  rw [uIcc_of_le hab] at hx
  have hz : HasDerivAt riemannZeta (deriv riemannZeta (x + T * I)) (x + T * I) :=
    (differentiableAt_riemannZeta (ne_one_of_im hT)).hasDerivAt
  have hγ : HasDerivAt (fun y : ℝ => (y : ℂ) + T * I) 1 x := by
    simpa using ((hasDerivAt_id x).ofReal_comp).add_const ((T : ℂ) * I)
  have hφ := (hz.comp x hγ).mul_const ((starRingEnd ℂ) c)
  have hlog := hφ.clog_real (mem_slitPlane_iff.mpr (Or.inl (hslit x hx)))
  refine hlog.congr_deriv ?_
  obtain ⟨h1, h2⟩ := hnz x hx
  simp only [Function.comp, logDeriv_apply]
  field_simp

lemma vert_ftc {σ a b : ℝ} (hσ : σ ≠ 1) (hab : a ≤ b) {c : ℂ}
    (hslit : ∀ y ∈ Icc a b, 0 < (riemannZeta (σ + y * I) * (starRingEnd ℂ) c).re) :
    (∫ y in a..b, logDeriv riemannZeta (σ + y * I)) * I =
      log (riemannZeta (σ + b * I) * (starRingEnd ℂ) c)
        - log (riemannZeta (σ + a * I) * (starRingEnd ℂ) c) := by
  have hnz : ∀ y ∈ Icc a b, riemannZeta (σ + y * I) ≠ 0 ∧ (starRingEnd ℂ) c ≠ 0 := by
    intro y hy
    have h := hslit y hy
    constructor <;> intro h0 <;> rw [h0] at h <;> simp at h
  rw [← intervalIntegral.integral_mul_const]
  refine intervalIntegral.integral_eq_sub_of_hasDerivAt
    (f := fun y : ℝ => log (riemannZeta (σ + y * I) * (starRingEnd ℂ) c)) ?_
    ((intInt_vert hσ hab fun y hy => (hnz y hy).1).mul_const I)
  intro y hy
  rw [uIcc_of_le hab] at hy
  have hz : HasDerivAt riemannZeta (deriv riemannZeta (σ + y * I)) (σ + y * I) :=
    (differentiableAt_riemannZeta (ne_one_of_re hσ)).hasDerivAt
  have hγ : HasDerivAt (fun y : ℝ => (σ : ℂ) + y * I) I y := by
    simpa using (((hasDerivAt_id y).ofReal_comp).mul_const I).const_add (σ : ℂ)
  have hφ := (hz.comp y hγ).mul_const ((starRingEnd ℂ) c)
  have hlog := hφ.clog_real (mem_slitPlane_iff.mpr (Or.inl (hslit y hy)))
  refine hlog.congr_deriv ?_
  obtain ⟨h1, h2⟩ := hnz y hy
  simp only [Function.comp, logDeriv_apply]
  field_simp

/-! ## Pieces with bounds -/

/-- Top piece on [a, b]: ζ ≠ 0 and −Im ∫_a^b ζ'/ζ(x + iT) dx ≤ U. -/
def TopOK (T a b U : ℝ) : Prop :=
  (∀ x ∈ Icc a b, riemannZeta (x + T * I) ≠ 0) ∧
    -(∫ x in a..b, logDeriv riemannZeta (x + T * I)).im ≤ U

/-- Vertical piece on [a, b]: ζ ≠ 0 and −Re ∫_a^b ζ'/ζ(σ + iy) dy ≤ U. -/
def LeftOK (σ a b U : ℝ) : Prop :=
  (∀ y ∈ Icc a b, riemannZeta (σ + y * I) ≠ 0) ∧
    -(∫ y in a..b, logDeriv riemannZeta (σ + y * I)).re ≤ U

lemma TopOK.trans {T a b c U1 U2 : ℝ} (hT : T ≠ 0) (hab : a ≤ b) (hbc : b ≤ c)
    (h1 : TopOK T a b U1) (h2 : TopOK T b c U2) : TopOK T a c (U1 + U2) := by
  refine ⟨fun x hx => ?_, ?_⟩
  · rcases le_total x b with h | h
    · exact h1.1 x ⟨hx.1, h⟩
    · exact h2.1 x ⟨h, hx.2⟩
  · rw [← intervalIntegral.integral_add_adjacent_intervals (intInt_top hT hab h1.1)
      (intInt_top hT hbc h2.1), add_im]
    linarith [h1.2, h2.2]

lemma LeftOK.trans {σ a b c U1 U2 : ℝ} (hσ : σ ≠ 1) (hab : a ≤ b) (hbc : b ≤ c)
    (h1 : LeftOK σ a b U1) (h2 : LeftOK σ b c U2) : LeftOK σ a c (U1 + U2) := by
  refine ⟨fun x hx => ?_, ?_⟩
  · rcases le_total x b with h | h
    · exact h1.1 x ⟨hx.1, h⟩
    · exact h2.1 x ⟨h, hx.2⟩
  · rw [← intervalIntegral.integral_add_adjacent_intervals (intInt_vert hσ hab h1.1)
      (intInt_vert hσ hbc h2.1), add_re]
    linarith [h1.2, h2.2]

lemma topOK_of_slit {T a b : ℝ} (hT : T ≠ 0) (hab : a ≤ b) {c : ℂ}
    (hslit : ∀ x ∈ Icc a b, 0 < (riemannZeta (x + T * I) * (starRingEnd ℂ) c).re) {Ua Lb U : ℝ}
    (hUa : arg (riemannZeta (a + T * I) * (starRingEnd ℂ) c) ≤ Ua)
    (hLb : Lb ≤ arg (riemannZeta (b + T * I) * (starRingEnd ℂ) c)) (hU : Ua - Lb ≤ U) :
    TopOK T a b U := by
  refine ⟨fun x hx h0 => ?_, ?_⟩
  · have := hslit x hx; rw [h0] at this; simp at this
  · rw [top_ftc hT hab hslit, sub_im, log_im, log_im]
    linarith

lemma leftOK_of_slit {σ a b : ℝ} (hσ : σ ≠ 1) (hab : a ≤ b) {c : ℂ}
    (hslit : ∀ y ∈ Icc a b, 0 < (riemannZeta (σ + y * I) * (starRingEnd ℂ) c).re) {Ua Lb U : ℝ}
    (hUa : arg (riemannZeta (σ + a * I) * (starRingEnd ℂ) c) ≤ Ua)
    (hLb : Lb ≤ arg (riemannZeta (σ + b * I) * (starRingEnd ℂ) c)) (hU : Ua - Lb ≤ U) :
    LeftOK σ a b U := by
  refine ⟨fun x hx h0 => ?_, ?_⟩
  · have := hslit x hx; rw [h0] at this; simp at this
  · have h := congrArg Complex.im (vert_ftc hσ hab hslit)
    rw [mul_I_im, sub_im, log_im, log_im] at h
    rw [h]
    linarith

/-- The right edge Re s = 2: Re ζ > 0, so its arg change is arg ζ(2 + iT). -/
lemma right_edge {T : ℝ} (hT : 0 ≤ T) :
    (∫ y in (0 : ℝ)..T, logDeriv riemannZeta ((2 : ℝ) + y * I)).re =
      arg (riemannZeta ((2 : ℝ) + T * I)) := by
  have hslit : ∀ y ∈ Icc (0 : ℝ) T,
      0 < (riemannZeta ((2 : ℝ) + y * I) * (starRingEnd ℂ) 1).re := by
    intro y _
    rw [map_one, mul_one]
    have hre : (2 : ℝ) ≤ ((2 : ℝ) + y * I : ℂ).re := by simp
    have h := Zeta23.RvM.norm_riemannZeta_sub_one_le hre
    have hpi : Real.pi ^ 2 / 6 - 1 < 1 := by
      have := Real.pi_lt_d2; nlinarith [Real.pi_pos]
    have := Complex.abs_re_le_norm (riemannZeta ((2 : ℝ) + y * I) - 1)
    rw [abs_le, sub_re, one_re] at this
    linarith [this.1]
  have h := congrArg Complex.im (vert_ftc (by norm_num) hT hslit)
  rw [mul_I_im, sub_im, log_im, log_im, map_one, mul_one, mul_one] at h
  rw [h]
  have h2 : riemannZeta ((2 : ℝ) + (0 : ℝ) * I) = ((Real.pi ^ 2 / 6 : ℝ) : ℂ) := by
    simp only [ofReal_zero, zero_mul, add_zero]
    rw [show ((2 : ℝ) : ℂ) = 2 by norm_num, riemannZeta_two]
    push_cast; ring
  rw [h2, arg_ofReal_of_nonneg (by positivity), sub_zero]

/-! ## The folded argument principle -/

lemma intervalIntegral_conj {f : ℝ → ℂ} {a b : ℝ} (hf : IntervalIntegrable f volume a b) :
    ∫ x in a..b, (starRingEnd ℂ) (f x) = (starRingEnd ℂ) (∫ x in a..b, f x) :=
  Complex.conjCLE.toContinuousLinearMap.intervalIntegral_comp_comm hf

lemma conj_top (x T : ℝ) : (starRingEnd ℂ) ((x : ℂ) + T * I) = (x : ℂ) + ((-T : ℝ) : ℂ) * I := by
  simp [Complex.conj_ofReal]

lemma vert_full {σ T : ℝ} (hσ : σ ≠ 1) (hT : 0 ≤ T)
    (h0 : ∀ y ∈ Icc 0 T, riemannZeta (σ + y * I) ≠ 0) :
    IntervalIntegrable (fun y : ℝ => logDeriv riemannZeta (σ + y * I)) volume (-T) 0 ∧
    IntervalIntegrable (fun y : ℝ => logDeriv riemannZeta (σ + y * I)) volume 0 T ∧
    ∫ y in (-T)..T, logDeriv riemannZeta (σ + y * I) =
      (starRingEnd ℂ) (∫ y in (0 : ℝ)..T, logDeriv riemannZeta (σ + y * I))
        + ∫ y in (0 : ℝ)..T, logDeriv riemannZeta (σ + y * I) := by
  have hc : ∀ y : ℝ, logDeriv riemannZeta (σ + ((-y : ℝ) : ℂ) * I) =
      (starRingEnd ℂ) (logDeriv riemannZeta (σ + y * I)) := by
    intro y
    rw [← logDeriv_zeta_conj (ne_one_of_re hσ)]
    congr 1
    simp [Complex.conj_ofReal]
  have h0' : ∀ y ∈ Icc (-T) 0, riemannZeta (σ + y * I) ≠ 0 := by
    intro y hy hz
    have h1 := h0 (-y) ⟨by linarith [hy.2], by linarith [hy.1]⟩
    have e : (σ : ℂ) + ((-y : ℝ) : ℂ) * I = (starRingEnd ℂ) ((σ : ℂ) + y * I) := by
      simp [Complex.conj_ofReal]
    rw [e, riemannZeta_conj, hz, map_zero] at h1
    exact h1 rfl
  have i1 := intInt_vert hσ (by linarith) h0'
  have i2 := intInt_vert hσ hT h0
  refine ⟨i1, i2, ?_⟩
  rw [← intervalIntegral.integral_add_adjacent_intervals i1 i2]
  congr 1
  have e := intervalIntegral.integral_comp_neg (a := 0) (b := T)
    (fun y : ℝ => logDeriv riemannZeta (σ + y * I))
  try simp only [neg_zero] at e
  rw [← e]
  simp only [hc]
  exact intervalIntegral_conj i2

theorem winding {σL T : ℝ} (hσ0 : 0 < σL) (hσ1 : σL < 1) (hT : 0 < T)
    (hleft : ∀ y ∈ Icc 0 T, riemannZeta (σL + y * I) ≠ 0)
    (htop : ∀ x ∈ Icc σL 2, riemannZeta (x + T * I) ≠ 0) :
    ∃ Z : Finset ℂ,
      (∀ s ∈ Rectangle ((σL : ℂ) - T * I) (2 + T * I), s ≠ 1 → riemannZeta s = 0 → s ∈ Z) ∧
      ((∑ ρ ∈ Z, (analyticOrderNatAt riemannZeta ρ : ℝ)) - 1) * Real.pi =
        -(∫ x in σL..2, logDeriv riemannZeta (x + T * I)).im
          + (∫ y in (0 : ℝ)..T, logDeriv riemannZeta ((2 : ℝ) + y * I)).re
          - (∫ y in (0 : ℝ)..T, logDeriv riemannZeta (σL + y * I)).re := by
  classical
  set z : ℂ := (σL : ℂ) - T * I with hz
  set w : ℂ := 2 + T * I with hw
  have hzre : z.re = σL := by simp [z]
  have hzim : z.im = -T := by simp [z]
  have hwre : w.re = 2 := by simp [w]
  have hwim : w.im = T := by simp [w]
  have hcomp : IsCompact (Rectangle z w) := isCompact_uIcc.reProdIm isCompact_uIcc
  have hfin := Zeta23.RvM.riemannZeta_zeros_finite_of_isCompact hcomp
  set Z := hfin.toFinset with hZdef
  have hZmem : ∀ s, s ∈ Z ↔ s ∈ Rectangle z w ∧ s ≠ 1 ∧ riemannZeta s = 0 := by
    intro s; rw [hZdef, Set.Finite.mem_toFinset]; rfl
  refine ⟨Z, fun s hs hs1 h0 => (hZmem s).mpr ⟨hs, hs1, h0⟩, ?_⟩
  -- nonvanishing on the whole border
  have hstd : ∀ s : ℂ, (s.re : ℂ) + s.im * I = s := fun s => re_add_im s
  have hborder : ∀ s ∈ RectangleBorder z w, riemannZeta s ≠ 0 := by
    intro s hs
    simp only [RectangleBorder, mem_union, mem_reProdIm, mem_singleton_iff, hzre, hzim, hwre,
      hwim, uIcc_of_le (show σL ≤ 2 by linarith), uIcc_of_le (show -T ≤ T by linarith),
      mem_Icc] at hs
    rcases hs with ((⟨hr, hi⟩ | ⟨hr, hi⟩) | ⟨hr, hi⟩) | ⟨hr, _⟩
    · -- bottom
      have h1 := htop s.re hr
      have e : (s.re : ℂ) + T * I = (starRingEnd ℂ) s := by
        apply Complex.ext <;> simp [hi]
      rw [e, riemannZeta_conj] at h1
      intro h0; rw [h0, map_zero] at h1; exact h1 rfl
    · -- left
      rcases le_total 0 s.im with hp | hn
      · have h1 := hleft s.im ⟨hp, hi.2⟩
        rw [← hr, hstd] at h1; exact h1
      · have h1 := hleft (-s.im) ⟨by linarith, by linarith [hi.1]⟩
        have e : (σL : ℂ) + ((-s.im : ℝ) : ℂ) * I = (starRingEnd ℂ) s := by
          apply Complex.ext <;> simp [hr]
        rw [e, riemannZeta_conj] at h1
        intro h0; rw [h0, map_zero] at h1; exact h1 rfl
    · -- top
      have h1 := htop s.re hr
      rw [← hi, hstd] at h1; exact h1
    · -- right
      exact riemannZeta_ne_zero_of_one_le_re (by rw [hr]; norm_num)
  have hPint : ∀ p ∈ ({1} : Finset ℂ), Rectangle z w ∈ 𝓝 p := by
    intro p hp
    rw [Finset.mem_singleton] at hp
    rw [hp, rectangle_mem_nhds_iff, mem_reProdIm, hzre, hzim, hwre, hwim]
    simp only [one_re, one_im, uIoo, mem_Ioo]
    constructor
    · constructor
      · exact lt_of_le_of_lt (min_le_left _ _) hσ1
      · exact lt_of_lt_of_le (by norm_num) (le_max_right _ _)
    · constructor
      · exact lt_of_le_of_lt (min_le_left _ _) (by linarith)
      · exact lt_of_lt_of_le hT (le_max_right _ _)
  have hf : AnalyticOnNhd ℂ riemannZeta (Rectangle z w \ (({1} : Finset ℂ) : Set ℂ)) := by
    intro s hs
    exact Zeta23.riemannZeta_analyticOnNhd_compl_one s (by simpa using hs.2)
  have hg : AnalyticOnNhd ℂ (fun _ : ℂ => (1 : ℂ)) (Rectangle z w) := fun _ _ => analyticAt_const
  have hZ : ∀ s ∈ Rectangle z w \ (({1} : Finset ℂ) : Set ℂ), riemannZeta s = 0 ↔ s ∈ Z := by
    intro s hs
    rw [hZmem]
    have hs1 : s ≠ 1 := by simpa using hs.2
    exact ⟨fun h => ⟨hs.1, hs1, h⟩, fun h => h.2.2⟩
  have hZsub : (Z : Set ℂ) ⊆ Rectangle z w := fun s hs => ((hZmem s).mp hs).1
  have hdisj : Disjoint Z ({1} : Finset ℂ) := by
    rw [Finset.disjoint_singleton_right]
    intro h; exact ((hZmem 1).mp h).2.1 rfl
  have hpole : ∀ p ∈ ({1} : Finset ℂ), ∃ c : ℂ, c ≠ 0 ∧
      Tendsto (fun s => (s - p) ^ (fun _ : ℂ => 1) p * riemannZeta s) (𝓝[≠] p) (𝓝 c) := by
    intro p hp
    rw [Finset.mem_singleton] at hp
    subst hp
    exact ⟨1, one_ne_zero, by simpa using riemannZeta_residue_one⟩
  have key := Zeta23.Analytic.rectangleIntegral'_mul_logDeriv_of_poles
    (f := riemannZeta) (g := fun _ => 1) (by rw [hzre, hwre]; linarith)
    (by rw [hzim, hwim]; linarith) Z {1} hdisj hPint hf hg hborder hZ hZsub (fun _ => 1) hpole
  simp only [one_mul, mul_one, Finset.sum_singleton, Nat.cast_one] at key
  -- evaluate the rectangle integral by the conjugation symmetry
  have htopI := intInt_top (ne_of_gt hT) (show σL ≤ 2 by linarith) htop
  have hbot : (∫ x in σL..2, logDeriv riemannZeta (x + ((-T : ℝ) : ℂ) * I)) =
      (starRingEnd ℂ) (∫ x in σL..2, logDeriv riemannZeta (x + T * I)) := by
    rw [← intervalIntegral_conj htopI]
    congr 1; funext x
    rw [← conj_top, logDeriv_zeta_conj (ne_one_of_im (ne_of_gt hT))]
  have hright := vert_full (σ := 2) (by norm_num) hT.le (fun y _ =>
    riemannZeta_ne_zero_of_one_le_re (by simp))
  have hleftf := vert_full (σ := σL) (by linarith) hT.le hleft
  have hRI : RectangleIntegral (logDeriv riemannZeta) z w =
      2 * I * ((-(∫ x in σL..2, logDeriv riemannZeta (x + T * I)).im
          + (∫ y in (0 : ℝ)..T, logDeriv riemannZeta ((2 : ℝ) + y * I)).re
          - (∫ y in (0 : ℝ)..T, logDeriv riemannZeta (σL + y * I)).re : ℝ) : ℂ) := by
    simp only [RectangleIntegral, HIntegral, VIntegral, hzre, hzim, hwre, hwim]
    rw [hbot, hright.2.2, hleftf.2.2]
    apply Complex.ext <;> simp <;> ring
  rw [show (fun s => logDeriv riemannZeta s) = logDeriv riemannZeta from rfl] at key
  simp only [RectangleIntegral', smul_eq_mul, hRI] at key
  set Q : ℝ := -(∫ x in σL..2, logDeriv riemannZeta (x + T * I)).im
          + (∫ y in (0 : ℝ)..T, logDeriv riemannZeta ((2 : ℝ) + y * I)).re
          - (∫ y in (0 : ℝ)..T, logDeriv riemannZeta (σL + y * I)).re with hQ
  have hpi : (Real.pi : ℂ) ≠ 0 := by exact_mod_cast Real.pi_ne_zero
  have hk : ((Q / Real.pi : ℝ) : ℂ) =
      (((∑ ρ ∈ Z, (analyticOrderNatAt riemannZeta ρ : ℝ)) - 1 : ℝ) : ℂ) := by
    have e : (((∑ ρ ∈ Z, (analyticOrderNatAt riemannZeta ρ : ℝ)) - 1 : ℝ) : ℂ) =
        ∑ x ∈ Z, (analyticOrderNatAt riemannZeta x : ℂ) - 1 := by push_cast; rfl
    rw [e, ← key]
    push_cast
    field_simp
  have hk' := Complex.ofReal_injective hk
  rw [← hk']
  field_simp

end Fayed.FirstZero.Min
