import Riemann.Deep.FirstZero.MinLib

/-!
# Simplicity of the 2nd–5th zeros, part 1: the argument principle on an off-axis box

For a box B = [σL, 2] × [a, b] with 0 < a (so the pole s = 1 is outside), the argument principle
(zeta23 `rectangleIntegral'_mul_logDeriv`) gives
  2π · Σ_{ρ ∈ B} ord_ρ ζ = Im ∫_bottom − Im ∫_top + Re ∫_right − Re ∫_left   (of ζ'/ζ).
The bottom edge needs an upper bound for +Im ∫ (`BotOK`, `botOK_seg`); top and left reuse
`TopOK` / `LeftOK` of MinArg; the right edge Re s = 2 is free (Re ζ > 0 there).
A certified total < 4π forces Σ ord < 2: a zero 1/2 + it in B is simple and the only zero in B.

Also: log n to ~3·10⁻⁹ for 9 ≤ n ≤ 20 (Euler–Maclaurin main terms with N ≤ 20).
-/

open Complex Set Filter Topology MeasureTheory

set_option linter.unusedSimpArgs false

namespace Fayed.FirstZero.SimpleNext
open Fayed.FirstZero Fayed.FirstZero.Min

/-! ## Logarithms -/

lemma log17_near : |Real.log 17 - 28332133441 / 10 ^ 10| ≤ 3 / 10 ^ 9 := by
  have h := Real.abs_log_sub_add_sum_range_le (x := (-1 / 16 : ℝ)) (by norm_num [abs_of_neg]) 20
  have e : Real.log 17 = 4 * Real.log 2 + Real.log (1 - (-1 / 16)) := by
    rw [show (1 : ℝ) - (-1 / 16) = 17 / 2 ^ 4 by norm_num, Real.log_div (by norm_num) (by norm_num),
      Real.log_pow]
    push_cast; ring
  have h2 := log2_near
  rw [abs_le] at h h2 ⊢
  norm_num [Finset.sum_range_succ, abs_of_neg] at h h2 ⊢
  rw [e]
  norm_num
  constructor <;> linarith [h.1, h.2, h2.1, h2.2]

lemma log19_near : |Real.log 19 - 29444389792 / 10 ^ 10| ≤ 3 / 10 ^ 9 := by
  have h := Real.abs_log_sub_add_sum_range_le (x := (-1 / 18 : ℝ)) (by norm_num [abs_of_neg]) 20
  have e : Real.log 19 = Real.log 2 + 2 * Real.log 3 + Real.log (1 - (-1 / 18)) := by
    rw [show (1 : ℝ) - (-1 / 18) = 19 / (2 * 3 ^ 2) by norm_num,
      Real.log_div (by norm_num) (by norm_num), Real.log_mul (by norm_num) (by norm_num),
      Real.log_pow]
    push_cast; ring
  have h2 := log2_near
  have h3 := log3_near
  rw [abs_le] at h h2 h3 ⊢
  norm_num [Finset.sum_range_succ, abs_of_neg] at h h2 h3 ⊢
  rw [e]
  norm_num
  constructor <;> linarith [h.1, h.2, h2.1, h2.2, h3.1, h3.2]

lemma log11_near' : |Real.log ((11 : ℕ) : ℝ) - 23978952728 / 10 ^ 10| ≤ 3 / 10 ^ 9 := by
  simpa using log11_near
lemma log13_near' : |Real.log ((13 : ℕ) : ℝ) - 25649493575 / 10 ^ 10| ≤ 3 / 10 ^ 9 := by
  simpa using log13_near
lemma log17_near' : |Real.log ((17 : ℕ) : ℝ) - 28332133441 / 10 ^ 10| ≤ 3 / 10 ^ 9 := by
  simpa using log17_near
lemma log19_near' : |Real.log ((19 : ℕ) : ℝ) - 29444389792 / 10 ^ 10| ≤ 3 / 10 ^ 9 := by
  simpa using log19_near

lemma log9_near' : |Real.log ((9 : ℕ) : ℝ) - (10986122887 / 10 ^ 10 + 10986122887 / 10 ^ 10)| ≤ 1 / 10 ^ 9 + 1 / 10 ^ 9 :=
  log_near_mul (m := 3) (n := 3) (by norm_num) (by norm_num) log3_near' log3_near'
lemma log10_near' : |Real.log ((10 : ℕ) : ℝ) - (69314718055 / 10 ^ 11 + 16094379124 / 10 ^ 10)| ≤ 3 / 10 ^ 10 + 1 / 10 ^ 9 :=
  log_near_mul (m := 2) (n := 5) (by norm_num) (by norm_num) log2_near' log5_near'
lemma log12_near' : |Real.log ((12 : ℕ) : ℝ) - (69314718055 / 10 ^ 11 + (69314718055 / 10 ^ 11 + 10986122887 / 10 ^ 10))| ≤ 3 / 10 ^ 10 + (3 / 10 ^ 10 + 1 / 10 ^ 9) :=
  log_near_mul (m := 2) (n := 6) (by norm_num) (by norm_num) log2_near' log6_near'
lemma log14_near' : |Real.log ((14 : ℕ) : ℝ) - (69314718055 / 10 ^ 11 + 19459101491 / 10 ^ 10)| ≤ 3 / 10 ^ 10 + 3 / 10 ^ 9 :=
  log_near_mul (m := 2) (n := 7) (by norm_num) (by norm_num) log2_near' log7_near'
lemma log15_near' : |Real.log ((15 : ℕ) : ℝ) - (10986122887 / 10 ^ 10 + 16094379124 / 10 ^ 10)| ≤ 1 / 10 ^ 9 + 1 / 10 ^ 9 :=
  log_near_mul (m := 3) (n := 5) (by norm_num) (by norm_num) log3_near' log5_near'
lemma log16_near' : |Real.log ((16 : ℕ) : ℝ) - (69314718055 / 10 ^ 11 + (69314718055 / 10 ^ 11 + (69314718055 / 10 ^ 11 + 69314718055 / 10 ^ 11)))| ≤ 3 / 10 ^ 10 + (3 / 10 ^ 10 + (3 / 10 ^ 10 + 3 / 10 ^ 10)) :=
  log_near_mul (m := 2) (n := 8) (by norm_num) (by norm_num) log2_near' log8_near'
lemma log18_near' : |Real.log ((18 : ℕ) : ℝ) - (69314718055 / 10 ^ 11 + (10986122887 / 10 ^ 10 + 10986122887 / 10 ^ 10))| ≤ 3 / 10 ^ 10 + (1 / 10 ^ 9 + 1 / 10 ^ 9) :=
  log_near_mul (m := 2) (n := 9) (by norm_num) (by norm_num) log2_near' log9_near'
lemma log20_near' : |Real.log ((20 : ℕ) : ℝ) - (69314718055 / 10 ^ 11 + (69314718055 / 10 ^ 11 + 16094379124 / 10 ^ 10))| ≤ 3 / 10 ^ 10 + (3 / 10 ^ 10 + 1 / 10 ^ 9) :=
  log_near_mul (m := 2) (n := 10) (by norm_num) (by norm_num) log2_near' log10_near'

/-! ## Bottom edge: upper bound for +Im ∫ -/

/-- Bottom piece on [a, b] at height T: ζ ≠ 0 and Im ∫_a^b ζ'/ζ(x + iT) dx ≤ U. -/
def BotOK (T a b U : ℝ) : Prop :=
  (∀ x ∈ Icc a b, riemannZeta (x + T * I) ≠ 0) ∧
    (∫ x in a..b, logDeriv riemannZeta (x + T * I)).im ≤ U

lemma BotOK.trans {T a b c U1 U2 : ℝ} (hT : T ≠ 0) (hab : a ≤ b) (hbc : b ≤ c)
    (h1 : BotOK T a b U1) (h2 : BotOK T b c U2) : BotOK T a c (U1 + U2) := by
  refine ⟨fun x hx => ?_, ?_⟩
  · rcases le_total x b with h | h
    · exact h1.1 x ⟨hx.1, h⟩
    · exact h2.1 x ⟨h, hx.2⟩
  · rw [← intervalIntegral.integral_add_adjacent_intervals (intInt_top hT hab h1.1)
      (intInt_top hT hbc h2.1), add_im]
    linarith [h1.2, h2.2]

lemma botOK_of_slit {T a b : ℝ} (hT : T ≠ 0) (hab : a ≤ b) {c : ℂ}
    (hslit : ∀ x ∈ Icc a b, 0 < (riemannZeta (x + T * I) * (starRingEnd ℂ) c).re) {La Ub U : ℝ}
    (hLa : La ≤ arg (riemannZeta (a + T * I) * (starRingEnd ℂ) c))
    (hUb : arg (riemannZeta (b + T * I) * (starRingEnd ℂ) c) ≤ Ub) (hU : Ub - La ≤ U) :
    BotOK T a b U := by
  refine ⟨fun x hx h0 => ?_, ?_⟩
  · have := hslit x hx; rw [h0] at this; simp at this
  · rw [top_ftc hT hab hslit, sub_im, log_im, log_im]
    linarith

/-- One certified bottom piece (as `topOK_seg`, with the arg bounds swapped). -/
lemma botOK_seg {T a b : ℝ} (hT : 0 < T) (ha : 0 < a) (hab : a < b) {N : ℕ}
    (hN : 2 ≤ N) {Ar Ai Br Bi rA rB : ℝ}
    (hA : ‖Mf N ((a : ℂ) + T * I) - (Ar + Ai * I)‖ ≤ rA)
    (hB : ‖M1f N ((a : ℂ) + T * I) - (Br + Bi * I)‖ ≤ rB)
    {R v v2 L L1 C ε : ℝ} (hR0 : 0 ≤ R) (hR : b ^ 2 + T ^ 2 ≤ R ^ 2)
    (hv : (N : ℝ) ^ (-a) ≤ v) (hv2 : (2 : ℝ) ^ (-a) ≤ v2) (hL : Real.log N ≤ L)
    (hL1 : Real.log ((N - 1 : ℕ) : ℝ) ≤ L1)
    (hC : ((N - 2 : ℕ) : ℝ) * (L1 ^ 2 * v2) +
      v * (K2b N R T + 2 * L * K1b N R T + L ^ 2 * K0b N R T) ≤ C)
    (cj : ℕ → ℝ) (hc : ∀ j : ℕ, j < 7 → 0 ≤ cj j ∧ (b + 1 + j) ^ 2 + T ^ 2 ≤ cj j ^ 2)
    (hε : R * ((∏ j ∈ Finset.range 7, cj j) / 1000000 * (v / (N : ℝ) ^ 7 / (a + 7))) ≤ ε)
    {cr ci cn E : ℝ} (hcn0 : 0 ≤ cn) (hcn : cr ^ 2 + ci ^ 2 ≤ cn ^ 2)
    (hE : ε + rA + rB * (b - a) + C * (b - a) ^ 2 ≤ E)
    (h0 : cn * E < Ar * cr + Ai * ci)
    (h1 : cn * E < (Ar + (b - a) * Br) * cr + (Ai + (b - a) * Bi) * ci)
    {La Ub U : ℝ}
    (hLa : La ≤ arg (riemannZeta ((a : ℂ) + T * I) * (starRingEnd ℂ) ((cr : ℂ) + ci * I)))
    (hUb : arg (riemannZeta ((b : ℂ) + T * I) * (starRingEnd ℂ) ((cr : ℂ) + ci * I)) ≤ Ub)
    (hU : Ub - La ≤ U) : BotOK T a b U := by
  set h := b - a with hh
  have hh0 : 0 < h := by rw [hh]; linarith
  set p : ℂ := (a : ℂ) + T * I with hp
  have hpu : ∀ u : ℝ, p + (u : ℂ) * 1 = ((a + u : ℝ) : ℂ) + (T : ℂ) * I := by
    intro u; rw [hp]; push_cast; ring
  have hs1 : ∀ u ∈ Icc (0 : ℝ) h, p + (u : ℂ) * 1 ≠ 1 := by
    intro u _; rw [hpu]; exact ne_one_of_im hT.ne'
  have hreu : ∀ u : ℝ, (p + (u : ℂ) * 1).re = a + u := by intro u; rw [hpu]; simp
  have himu : ∀ u : ℝ, (p + (u : ℂ) * 1).im = T := by intro u; rw [hpu]; simp
  have hnorm : ∀ u ∈ Icc (0 : ℝ) h, ‖p + (u : ℂ) * 1‖ ≤ R := by
    intro u hu
    refine norm_le_of_sq_le hR0 ?_
    rw [hreu, himu]
    have : (a + u) ^ 2 ≤ b ^ 2 := pow_le_pow_left₀ (by linarith [hu.1]) (by linarith [hu.2]) 2
    linarith
  have hC' : ∀ u ∈ Icc (0 : ℝ) h, ‖M2f N (p + (u : ℂ) * 1)‖ ≤ C := by
    intro u hu
    refine le_trans (norm_M2f_le' hN (hnorm u hu) hT ?_ ha.le (by rw [hreu]; linarith [hu.1])
      hv hv2 hL hL1) hC
    have := Complex.abs_im_le_norm (1 - (p + (u : ℂ) * 1))
    rw [sub_im, one_im, himu, zero_sub, abs_neg, abs_of_pos hT] at this
    exact this
  have hε' : ∀ u ∈ Icc (0 : ℝ) h, ‖riemannZeta (p + (u : ℂ) * 1) - Mf N (p + (u : ℂ) * 1)‖ ≤ ε := by
    intro u hu
    refine em_err_le hN cj ha (by rw [hreu]; linarith [hu.1]) (hs1 u hu) (hnorm u hu) ?_ hv hε
    intro j hj
    obtain ⟨hc0, hc1⟩ := hc j hj
    have := norm_shift_le (s := p + (u : ℂ) * 1) (x1 := b) (Y := T) (a := 1 + j) (c := cj j)
      (by positivity) (by rw [hreu]; linarith [hu.1]) (by rw [hreu]; linarith [hu.2])
      (by rw [himu, abs_of_pos hT]) hc0 (by linarith)
    simpa [add_assoc] using this
  have hmod := seg_model hN (p := p) (d := 1) (by simp) hs1 hA hB hC' hε'
  have hrB : 0 ≤ rB := (norm_nonneg _).trans hB
  have hC0 : 0 ≤ C := (norm_nonneg _).trans (hC' 0 ⟨le_rfl, hh0.le⟩)
  have hmod' : ∀ u ∈ Icc (0 : ℝ) h, ‖riemannZeta (p + (u : ℂ) * 1) - (((Ar : ℂ) + Ai * I) +
      (u : ℂ) * 1 * ((Br : ℂ) + Bi * I))‖ ≤ E := by
    intro u hu
    refine (hmod u hu).trans (le_trans ?_ hE)
    have : rB * u ≤ rB * h := mul_le_mul_of_nonneg_left hu.2 hrB
    have : C * u ^ 2 ≤ C * h ^ 2 :=
      mul_le_mul_of_nonneg_left (pow_le_pow_left₀ hu.1 hu.2 2) hC0
    linarith
  have hslit := re_pos_of_model hh0 hmod' (norm_c_le hcn0 hcn)
    (by rw [conjc_re]; exact h0)
    (by
      have e : (((Ar : ℂ) + Ai * I) + (h : ℂ) * 1 * ((Br : ℂ) + Bi * I)) =
          ((Ar + h * Br : ℝ) : ℂ) + ((Ai + h * Bi : ℝ) : ℂ) * I := by
        push_cast; ring
      rw [e, conjc_re]; exact h1)
  refine botOK_of_slit hT.ne' hab.le (fun x hx => ?_) hLa hUb hU
  have := hslit (x - a) ⟨by linarith [hx.1], by rw [hh]; linarith [hx.2]⟩
  rwa [hpu, show a + (x - a) = x by ring] at this

/-- arg ζ from a ball (anchor 1), the remaining sign cases. -/
lemma arg_lo_one {z : ℂ} {a b r L : ℝ} (hz : ‖z - (a + b * I)‖ ≤ r) (h1 : r < a)
    (h2 : 0 ≤ r - b) (h3 : L ≤ (b - r) / (a - r)) : L ≤ arg z := le_arg_of_ball_pos hz h1 h2 h3

lemma arg_lo_one_neg {z : ℂ} {a b r M L : ℝ} (hz : ‖z - (a + b * I)‖ ≤ r) (h1 : r < a)
    (h2 : r - b < 0) (hM : 0 ≤ M - r) (hM2 : a ^ 2 + b ^ 2 ≤ (M - r) ^ 2)
    (h3 : L ≤ (b - r) / M) : L ≤ arg z := le_arg_of_ball_neg hz h1 h2 hM hM2 h3

lemma arg_up_one_neg {z : ℂ} {a b r M U : ℝ} (hz : ‖z - (a + b * I)‖ ≤ r) (h1 : r < a)
    (h2 : b + r < 0) (hM : 0 ≤ M - r) (hM2 : a ^ 2 + b ^ 2 ≤ (M - r) ^ 2)
    (h3 : (b + r) / M ≤ U) : arg z ≤ U := arg_le_of_ball_neg hz h1 h2 hM hM2 h3

/-! ## Right edge and the box count -/

lemma right_edge_ab {a b : ℝ} (hab : a ≤ b) :
    (∫ y in a..b, logDeriv riemannZeta ((2 : ℝ) + y * I)).re =
      arg (riemannZeta ((2 : ℝ) + b * I)) - arg (riemannZeta ((2 : ℝ) + a * I)) := by
  have hslit : ∀ y ∈ Icc a b,
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
  have h := congrArg Complex.im (vert_ftc (by norm_num) hab hslit)
  rw [mul_I_im, sub_im, log_im, log_im, map_one, mul_one, mul_one] at h
  exact h

lemma mem_box {σL a b : ℝ} (hσ : σL ≤ 2) (hab : a ≤ b) {x : ℂ} :
    x ∈ Rectangle ((σL : ℂ) + a * I) (2 + b * I) ↔
      σL ≤ x.re ∧ x.re ≤ 2 ∧ a ≤ x.im ∧ x.im ≤ b := by
  simp only [Rectangle, mem_reProdIm, add_re, ofReal_re, mul_re, I_re, I_im, ofReal_im,
    add_im, mul_im, re_ofNat, im_ofNat]
  rw [uIcc_of_le (by norm_num; linarith), uIcc_of_le (by norm_num; linarith)]
  simp only [mem_Icc, mul_zero, sub_zero, mul_one, zero_add, add_zero]
  tauto

/-- **Box count**: with the certified arg change < 4π, the zeros of ζ in [σL, 2] × [a, b] have
total order < 2. -/
theorem box_count {σL a b Ubot Utop Uleft Lra Urb : ℝ} (hσ0 : 0 < σL) (hσ1 : σL < 1)
    (ha : 0 < a) (hab : a < b)
    (hbot : BotOK a σL 2 Ubot) (htop : TopOK b σL 2 Utop) (hleft : LeftOK σL a b Uleft)
    (hra : Lra ≤ arg (riemannZeta ((2 : ℝ) + a * I)))
    (hrb : arg (riemannZeta ((2 : ℝ) + b * I)) ≤ Urb)
    (hsum : Ubot + Utop + Uleft + (Urb - Lra) < 4 * Real.pi) :
    ∃ Z : Finset ℂ,
      (∀ s ∈ Rectangle ((σL : ℂ) + a * I) (2 + b * I), riemannZeta s = 0 → s ∈ Z) ∧
      (∑ ρ ∈ Z, (analyticOrderNatAt riemannZeta ρ : ℝ)) < 2 := by
  classical
  set z : ℂ := (σL : ℂ) + a * I with hz
  set w : ℂ := 2 + b * I with hw
  have hzre : z.re = σL := by simp [z]
  have hzim : z.im = a := by simp [z]
  have hwre : w.re = 2 := by simp [w]
  have hwim : w.im = b := by simp [w]
  have hσ2 : σL ≤ 2 := by linarith
  have hcomp : IsCompact (Rectangle z w) := isCompact_uIcc.reProdIm isCompact_uIcc
  have hfin := Zeta23.RvM.riemannZeta_zeros_finite_of_isCompact hcomp
  set Z := hfin.toFinset with hZdef
  have hZmem : ∀ s, s ∈ Z ↔ s ∈ Rectangle z w ∧ s ≠ 1 ∧ riemannZeta s = 0 := by
    intro s; rw [hZdef, Set.Finite.mem_toFinset]; rfl
  have hne1 : ∀ s ∈ Rectangle z w, s ≠ 1 := by
    intro s hs h
    have h3 := ((mem_box hσ2 hab.le).mp hs).2.2.1
    rw [h, one_im] at h3; linarith
  have hZ : ∀ s ∈ Rectangle z w, riemannZeta s = 0 ↔ s ∈ Z := by
    intro s hs
    rw [hZmem]
    exact ⟨fun h => ⟨hs, hne1 s hs, h⟩, fun h => h.2.2⟩
  refine ⟨Z, fun s hs h0 => (hZ s hs).mp h0, ?_⟩
  have hstd : ∀ s : ℂ, (s.re : ℂ) + s.im * I = s := fun s => re_add_im s
  have hborder : ∀ s ∈ RectangleBorder z w, riemannZeta s ≠ 0 := by
    intro s hs
    simp only [RectangleBorder, mem_union, mem_reProdIm, mem_singleton_iff, hzre, hzim, hwre,
      hwim, uIcc_of_le hσ2, uIcc_of_le hab.le, mem_Icc] at hs
    rcases hs with ((⟨hr, hi⟩ | ⟨hr, hi⟩) | ⟨hr, hi⟩) | ⟨hr, _⟩
    · have h1 := hbot.1 s.re hr
      rw [← hi, hstd] at h1; exact h1
    · have h1 := hleft.1 s.im hi
      rw [← hr, hstd] at h1; exact h1
    · have h1 := htop.1 s.re hr
      rw [← hi, hstd] at h1; exact h1
    · exact riemannZeta_ne_zero_of_one_le_re (by rw [hr]; norm_num)
  have hf : AnalyticOnNhd ℂ riemannZeta (Rectangle z w) := fun s hs =>
    Zeta23.riemannZeta_analyticOnNhd_compl_one s (hne1 s hs)
  have hg : AnalyticOnNhd ℂ (fun _ : ℂ => (1 : ℂ)) (Rectangle z w) := fun _ _ => analyticAt_const
  have hZsub : (Z : Set ℂ) ⊆ Rectangle z w := fun s hs => ((hZmem s).mp hs).1
  have key := Zeta23.Analytic.rectangleIntegral'_mul_logDeriv (f := riemannZeta)
    (g := fun _ => 1) (by rw [hzre, hwre]; exact hσ2) (by rw [hzim, hwim]; exact hab.le)
    hf hg hborder Z hZ hZsub
  simp only [one_mul, mul_one] at key
  set S : ℝ := ∑ ρ ∈ Z, (analyticOrderNatAt riemannZeta ρ : ℝ) with hS
  have hSc : ((S : ℝ) : ℂ) = ∑ ρ ∈ Z, (analyticOrderNatAt riemannZeta ρ : ℂ) := by
    rw [hS]; push_cast; rfl
  have hpi : (2 * (Real.pi : ℂ) * I) ≠ 0 := by
    simp [Real.pi_ne_zero, I_ne_zero]
  have hRI : RectangleIntegral (fun s => logDeriv riemannZeta s) z w =
      2 * (Real.pi : ℂ) * I * ((S : ℝ) : ℂ) := by
    rw [hSc, ← key, RectangleIntegral', smul_eq_mul]
    field_simp
  have him := congrArg Complex.im hRI
  simp only [RectangleIntegral, HIntegral, VIntegral, hzre, hzim, hwre, hwim, smul_eq_mul,
    sub_im, add_im, mul_im, I_re, I_im] at him
  simp only [mul_im, mul_re, ofReal_re, ofReal_im, I_re, I_im, re_ofNat, im_ofNat] at him
  have hR := right_edge_ab hab.le
  have hb := hbot.2
  have ht := htop.2
  have hl := hleft.2
  have hS2 : 2 * Real.pi * S < 4 * Real.pi := by
    nlinarith [hR, hb, ht, hl, him]
  have := Real.pi_pos
  nlinarith

/-- **A certified box makes a critical zero simple and unique.** If the arg change around
[σL, 2] × [a, b] is certified < 4π and ζ(1/2 + it) = 0 with a < t < b, then that zero has order 1
and is the only zero of ζ in the critical strip with a ≤ Im s ≤ b. -/
theorem simple_of_box {σL a b Ubot Utop Uleft Lra Urb : ℝ} (hσ0 : 0 < σL) (hσh : σL < 1 / 2)
    (ha : 0 < a) (hab : a < b)
    (hbot : BotOK a σL 2 Ubot) (htop : TopOK b σL 2 Utop) (hleft : LeftOK σL a b Uleft)
    (hra : Lra ≤ arg (riemannZeta ((2 : ℝ) + a * I)))
    (hrb : arg (riemannZeta ((2 : ℝ) + b * I)) ≤ Urb)
    (hsum : Ubot + Utop + Uleft + (Urb - Lra) < 4 * Real.pi)
    (t : ℝ) (hta : a < t) (htb : t < b) (hz : riemannZeta (1 / 2 + t * I) = 0) :
    analyticOrderAt riemannZeta (1 / 2 + t * I) = 1 ∧
      ∀ s : ℂ, riemannZeta s = 0 → 0 < s.re → s.re < 1 → a ≤ s.im → s.im ≤ b →
        s = 1 / 2 + t * I := by
  classical
  obtain ⟨Z, hZ, hlt⟩ := box_count hσ0 (by linarith) ha hab hbot htop hleft hra hrb hsum
  have hσ2 : σL ≤ 2 := by linarith
  set c : ℂ := 1 / 2 + (t : ℂ) * I with hc
  have hcre : c.re = 1 / 2 := by simp [c]
  have hcim : c.im = t := by simp [c]
  have hc1 : c ≠ 1 := by
    intro h; have := congrArg Complex.im h; rw [hcim] at this; simp at this; linarith
  have hcZ : c ∈ Z := hZ c ((mem_box hσ2 hab.le).mpr
    ⟨by rw [hcre]; linarith, by rw [hcre]; norm_num, by rw [hcim]; linarith,
      by rw [hcim]; linarith⟩) hz
  have hge : 1 ≤ analyticOrderNatAt riemannZeta c :=
    one_le_ord_of_zero hz (by rw [hcre]; norm_num) (by rw [hcre]; norm_num)
  have hsub := Finset.single_le_sum (f := fun ρ => (analyticOrderNatAt riemannZeta ρ : ℝ))
    (fun _ _ => Nat.cast_nonneg _) hcZ
  have hnat : analyticOrderNatAt riemannZeta c = 1 := by
    have : (analyticOrderNatAt riemannZeta c : ℝ) < 2 := lt_of_le_of_lt hsub hlt
    have : analyticOrderNatAt riemannZeta c < 2 := by exact_mod_cast this
    omega
  refine ⟨?_, ?_⟩
  · have hnetop := Zeta23.analyticOrderAt_riemannZeta_ne_top hc1
    rw [← Nat.cast_analyticOrderNatAt hnetop, hnat, Nat.cast_one]
  · -- uniqueness: a second zero in the box (or its reflection) would push the count to 2
    have key : ∀ x : ℂ, riemannZeta x = 0 → σL ≤ x.re → x.re < 1 → a ≤ x.im → x.im ≤ b →
        x = c := by
      intro x hx h1 h2 h3 h4
      by_contra hne
      have hxZ : x ∈ Z := hZ x ((mem_box hσ2 hab.le).mpr ⟨h1, by linarith, h3, h4⟩) hx
      have hx1 : 1 ≤ analyticOrderNatAt riemannZeta x := one_le_ord_of_zero hx (by linarith) h2
      have hpair := Finset.sum_le_sum_of_subset_of_nonneg
        (f := fun ρ => (analyticOrderNatAt riemannZeta ρ : ℝ))
        (show ({x, c} : Finset ℂ) ⊆ Z by
          intro y hy
          simp only [Finset.mem_insert, Finset.mem_singleton] at hy
          rcases hy with rfl | rfl
          exacts [hxZ, hcZ]) (fun _ _ _ => Nat.cast_nonneg _)
      rw [Finset.sum_insert (by simpa using hne), Finset.sum_singleton] at hpair
      have e1 : (1 : ℝ) ≤ analyticOrderNatAt riemannZeta x := by exact_mod_cast hx1
      have e2 : (1 : ℝ) ≤ analyticOrderNatAt riemannZeta c := by exact_mod_cast hge
      linarith
    intro s hs h0 h1 h2 h3
    rcases le_or_gt σL s.re with hσ | hσ
    · exact key s hs hσ h1 h2 h3
    · exfalso
      have hr := (Zeta23.zeta_reflect_zero s ⟨hs, h0, h1⟩)
      have hre : (Zeta23.reflect s).re = 1 - s.re := by simp [Zeta23.reflect]
      have hrim : (Zeta23.reflect s).im = s.im := by simp [Zeta23.reflect]
      have e := key _ hr.1 (by rw [hre]; linarith) (by rw [hre]; linarith) (by rw [hrim]; exact h2)
        (by rw [hrim]; exact h3)
      have := congrArg Complex.re e
      rw [hre, hcre] at this
      linarith

end Fayed.FirstZero.SimpleNext
