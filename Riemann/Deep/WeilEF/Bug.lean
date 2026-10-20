/-
Riemann/Deep/WeilEF/Bug.lean — the ORIGINAL (unparenthesised) `WeilPos L` of
`Riemann/Deep/Weil/Defs.lean` was FALSE for every `L > 0`. Defs.lean is now fixed; this file keeps a
frozen copy (`archTermOld`, `weilQOld`, `WeilPosOld`) of the buggy definitions.

Cause (`archTermOld_parse`): `∫ r, h r * Re ψ − g 0 * log π` puts the constant inside the
integrand, which is then not integrable (`g 0 = ‖f‖² > 0`), so `archTerm = 0` (Bochner junk value).
Witness: `f = φ − φ(· − r)` (`φ` a small bump): no prime enters, the archimedean term is `0`, and the
pole term `2 F(i/2) F(−i/2) < 0` because `F(i/2) = A(1 − e^{−r/2}) > 0`, `F(−i/2) = B(1 − e^{r/2}) < 0`.
Consequently, before the fix, every cited window Prop in Defs.lean (Yoshida, Zhu, Liu, Mori,
CertWindow) was false as stated in Lean, and `WeilCriterion` was equivalent to `¬ RH`.
-/
import Riemann.Deep.WeilEF.Suff

noncomputable section

open Complex MeasureTheory Set Filter
open scoped ComplexConjugate

namespace Fayed.WeilSuff.Bug

open Riemann.Deep.Weil Riemann.Deep.WeilEF Zeta23 Zeta23.EF Fayed.WeilSuff

/-- Frozen copy of the original, buggy `archTerm` (binder swallows `- g 0 * log π`). -/
def archTermOld (g : ℝ → ℂ) : ℂ :=
  (1 / (2 * Real.pi) : ℂ) *
      ∫ r : ℝ, ft g r * ((logDeriv Complex.Gamma (1 / 4 + I * r / 2)).re : ℂ) -
    g 0 * (Real.log Real.pi : ℂ)

/-- How `archTermOld` parses (`rfl`). -/
lemma archTermOld_parse (g : ℝ → ℂ) :
    archTermOld g = (1 / (2 * Real.pi) : ℂ) * ∫ r : ℝ,
      (ft g r * ((logDeriv Complex.Gamma (1 / 4 + I * r / 2)).re : ℂ) -
        g 0 * (Real.log Real.pi : ℂ)) := rfl

/-- The original Weil functional. -/
def weilQOld (f : ℝ → ℂ) : ℂ :=
  poleTerm (autocorr f) + archTermOld (autocorr f) - primeTerm (autocorr f)

/-- The original window positivity. -/
def WeilPosOld (L : ℝ) : Prop :=
  ∀ f : ℝ → ℂ, TestFn (L / 2) f → 0 ≤ (weilQOld f).re

/-- Bump of radius `r`. -/
def bumpR {r : ℝ} (hr : 0 < r) : ContDiffBump (0 : ℝ) := ⟨r / 2, r, half_pos hr, half_lt_self hr⟩

/-- `φ(u) = bump_r(u)` as a complex function. -/
def phi {r : ℝ} (hr : 0 < r) : ℝ → ℂ := fun u => (bumpR hr u : ℂ)

lemma testFn_phi {r : ℝ} (hr : 0 < r) : TestFn r (phi hr) := by
  refine ⟨ofRealCLM.contDiff.comp (bumpR hr).contDiff, ?_⟩
  intro u hu
  have hb : bumpR hr u ≠ 0 := by
    intro h0; apply hu; simp [phi, h0]
  have : u ∈ Metric.ball (0 : ℝ) r := by
    have := (bumpR hr).support_eq ▸ Function.mem_support.mpr hb
    simpa [bumpR] using this
  rw [Metric.mem_ball, dist_zero_right, Real.norm_eq_abs] at this
  exact ⟨by linarith [neg_abs_le u], by linarith [le_abs_self u]⟩

lemma paperFT_phi {r : ℝ} (hr : 0 < r) (b : ℝ) :
    paperFT (phi hr) (-(I * b)) = ((∫ u : ℝ, bumpR hr u * Real.exp (b * u) : ℝ) : ℂ) := by
  rw [← integral_complex_ofReal]
  unfold paperFT phi
  congr 1; funext u
  push_cast
  congr 2
  ring_nf
  rw [I_sq]
  ring

lemma integral_pos {r : ℝ} (hr : 0 < r) (b : ℝ) :
    0 < ∫ u : ℝ, bumpR hr u * Real.exp (b * u) := by
  refine Continuous.integral_pos_of_hasCompactSupport_nonneg_nonzero (x := 0)
    ((bumpR hr).continuous.mul (by fun_prop)) (bumpR hr).hasCompactSupport.mul_right
    (fun u => mul_nonneg (bumpR hr).nonneg (Real.exp_pos _).le) ?_
  rw [(bumpR hr).one_of_mem_closedBall (by simp [bumpR]; linarith)]
  simp

/-- The witness `f = φ − φ(· − r)`. -/
def wit {r : ℝ} (hr : 0 < r) : ℝ → ℂ := shiftComb (phi hr) (-1) r

lemma testFn_wit {r : ℝ} (hr : 0 < r) : TestFn (2 * r) (wit hr) := by
  have := testFn_shiftComb (testFn_phi hr) (-1) r
  rwa [abs_of_pos hr, show r + r = 2 * r by ring] at this

lemma wit_zero {r : ℝ} (hr : 0 < r) : wit hr 0 = 1 := by
  have h0 : bumpR hr 0 = 1 := (bumpR hr).one_of_mem_closedBall (by simp [bumpR]; linarith)
  have h1 : bumpR hr (-r) = 0 := by
    rw [← Function.notMem_support, (bumpR hr).support_eq]
    simp [bumpR, abs_of_pos hr]
  simp [wit, shiftComb, phi, h0, h1]

/-- `F(−(i·b))` for the witness, as a real number. -/
lemma paperFT_wit {r : ℝ} (hr : 0 < r) (b : ℝ) :
    paperFT (wit hr) (-(I * b)) =
      (((∫ u : ℝ, bumpR hr u * Real.exp (b * u)) * (1 - Real.exp (b * r)) : ℝ) : ℂ) := by
  rw [wit, paperFT_shiftComb (testFn_phi hr), paperFT_phi]
  have : cexp (I * -(I * (b : ℂ)) * (r : ℂ)) = ((Real.exp (b * r) : ℝ) : ℂ) := by
    push_cast; congr 1; ring_nf; rw [I_sq]; ring
  rw [this]; push_cast; ring

lemma autocorr_wit_zero_ne {r : ℝ} (hr : 0 < r) : autocorr (wit hr) 0 ≠ 0 := by
  have hf := testFn_wit hr
  have hpt : (fun y : ℝ => wit hr y * conj (wit hr (y - 0)))
      = fun y : ℝ => ((Complex.normSq (wit hr y) : ℝ) : ℂ) := by
    funext y; rw [sub_zero, Complex.mul_conj]
  unfold autocorr
  rw [hpt, integral_complex_ofReal]
  refine ofReal_ne_zero.mpr (ne_of_gt ?_)
  refine Continuous.integral_pos_of_hasCompactSupport_nonneg_nonzero (x := 0)
    (Complex.continuous_normSq.comp hf.continuous) ?_ (fun y => Complex.normSq_nonneg _) ?_
  · exact hf.hasCompactSupport.comp_left (g := Complex.normSq) (map_zero _)
  · simp [wit_zero hr]

/-- The archimedean term of the witness is the junk value `0`. -/
lemma archTerm_wit {r : ℝ} (hr : 0 < r) : archTermOld (autocorr (wit hr)) = 0 := by
  have hf := testFn_wit hr
  have hk := hf.autocorr_contDiff
  have hkc := hf.autocorr_hasCompactSupport
  rw [archTermOld_parse]
  suffices hni : ¬ Integrable (fun r : ℝ => ft (autocorr (wit hr)) r *
      ((logDeriv Complex.Gamma (1 / 4 + I * r / 2)).re : ℂ) - (autocorr (wit hr)) 0 * (Real.log Real.pi : ℂ)) by
    rw [integral_undef hni, mul_zero]
  intro hint
  -- the h·Re ψ part is integrable (as in `Bridge.arch_eq`)
  have hFk := integrable_fourier_of_contDiff_two hk hkc
  have hI : Integrable (fun r : ℝ => paperFT (autocorr (wit hr)) r) := integrable_paperFT_ofReal hFk
  have hmu := integrable_paperFT_mul_mu hk hkc gammaFacts
  have hIg : Integrable (fun r : ℝ => paperFT (autocorr (wit hr)) r * (gammaBracket r : ℂ)) := by
    have := hmu.const_mul ((2 * Real.pi : ℝ) : ℂ)
    refine this.congr (Filter.Eventually.of_forall fun r => ?_)
    simp only [mu, gammaBracket]
    push_cast
    field_simp
  have hpsi : Integrable (fun r : ℝ => ft (autocorr (wit hr)) r *
      ((logDeriv Complex.Gamma (1 / 4 + I * r / 2)).re : ℂ)) := by
    refine (hIg.add (hI.const_mul (Real.log Real.pi : ℂ))).congr
      (Filter.Eventually.of_forall fun r => ?_)
    simp only [Pi.add_apply, ft_eq_paperFT, gammaBracket, Complex.digamma_def]
    push_cast
    ring
  have hconst : Integrable (fun _ : ℝ => (autocorr (wit hr)) 0 * (Real.log Real.pi : ℂ)) := by
    refine (hpsi.sub hint).congr (Filter.Eventually.of_forall fun r => ?_)
    simp
  rcases integrable_const_iff.mp hconst with h0 | hfin
  · have hlog : (Real.log Real.pi : ℂ) ≠ 0 :=
      ofReal_ne_zero.mpr (Real.log_pos (by linarith [Real.pi_gt_three])).ne'
    exact (mul_ne_zero (autocorr_wit_zero_ne hr) hlog) h0
  · have := hfin.measure_univ_lt_top
    rw [Real.volume_univ] at this
    exact lt_irrefl _ this

/-- **The bug**: the original `WeilPos L` (`WeilPosOld`) fails for every `L > 0`. -/
theorem not_weilPos {L : ℝ} (hL : 0 < L) : ¬ WeilPosOld L := by
  intro hW
  set r := min (L / 4) (1 / 8) with hr_def
  have hr : 0 < r := lt_min (by linarith) (by norm_num)
  have hrL : 2 * r ≤ L / 2 := by have := min_le_left (L / 4) (1 / 8); linarith
  have hr8 : r ≤ 1 / 8 := min_le_right _ _
  have hf := testFn_wit hr
  have hpos := hW (wit hr) (hf.mono hrL)
  -- prime term vanishes (support ⊆ [−2r, 2r], 4r < log 2)
  have hprime : primeTerm (autocorr (wit hr)) = 0 :=
    primeTerm_eq_zero hf.2 (by linarith [Real.log_two_gt_d9])
  -- pole term: 2 F(i/2) F(−i/2) < 0
  set A := ∫ u : ℝ, bumpR hr u * Real.exp ((-(1 / 2) : ℝ) * u) with hA
  set B := ∫ u : ℝ, bumpR hr u * Real.exp ((1 / 2 : ℝ) * u) with hB
  have hA0 : 0 < A := integral_pos hr _
  have hB0 : 0 < B := integral_pos hr _
  have hP : paperFT (wit hr) (I / 2) = ((A * (1 - Real.exp ((-(1 / 2) : ℝ) * r)) : ℝ) : ℂ) := by
    rw [show I / 2 = -(I * (((-(1 / 2) : ℝ)) : ℂ)) by push_cast; ring, paperFT_wit]
  have hN : paperFT (wit hr) (-(I / 2)) = ((B * (1 - Real.exp ((1 / 2 : ℝ) * r)) : ℝ) : ℂ) := by
    rw [show -(I / 2) = -(I * (((1 / 2 : ℝ)) : ℂ)) by push_cast; ring, paperFT_wit]
  have hcI : conj (I / 2) = -(I / 2) := by rw [map_div₀, Complex.conj_I, map_ofNat]; ring
  have hcI' : conj (-(I / 2)) = I / 2 := by rw [map_neg, hcI, neg_neg]
  have hpole : poleTerm (autocorr (wit hr)) =
      ((2 * (A * (1 - Real.exp ((-(1 / 2) : ℝ) * r))) * (B * (1 - Real.exp ((1 / 2 : ℝ) * r))) : ℝ)
        : ℂ) := by
    unfold poleTerm
    rw [ft_eq_paperFT, ft_eq_paperFT, paperFT_autocorr hf, paperFT_autocorr hf, hcI, hcI', hP, hN,
      Complex.conj_ofReal, Complex.conj_ofReal]
    push_cast; ring
  have hval : weilQOld (wit hr) = ((2 * (A * (1 - Real.exp ((-(1 / 2) : ℝ) * r))) *
      (B * (1 - Real.exp ((1 / 2 : ℝ) * r))) : ℝ) : ℂ) := by
    rw [weilQOld, hpole, archTerm_wit hr, hprime]; ring
  rw [hval, ofReal_re] at hpos
  have h1 : 0 < 1 - Real.exp ((-(1 / 2) : ℝ) * r) := by
    have : Real.exp ((-(1 / 2) : ℝ) * r) < 1 := Real.exp_lt_one_iff.mpr (by nlinarith)
    linarith
  have h2 : 1 - Real.exp ((1 / 2 : ℝ) * r) < 0 := by
    have : 1 < Real.exp ((1 / 2 : ℝ) * r) := Real.one_lt_exp_iff.mpr (by positivity)
    linarith
  have : 2 * (A * (1 - Real.exp ((-(1 / 2) : ℝ) * r))) * (B * (1 - Real.exp ((1 / 2 : ℝ) * r))) < 0 :=
    mul_neg_of_pos_of_neg (by positivity) (mul_neg_of_pos_of_neg hB0 h2)
  linarith

/-- Hence `∀ L > 0, WeilPos L` (the hypothesis of the brief's theorem) is refutable outright. -/
theorem not_forall_weilPos : ¬ ∀ L : ℝ, 0 < L → WeilPosOld L := fun h => not_weilPos one_pos (h 1 one_pos)

end Fayed.WeilSuff.Bug
