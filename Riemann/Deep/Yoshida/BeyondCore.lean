/-
Riemann/Deep/Yoshida/BeyondCore.lean — Weil positivity beyond `log 2`: the Fourier form with the prime 2,
free cosine terms, chord bounds, and the affine moment certificate theorem.

For `f` supported in `[−a, a]` with `2a < log 3` only `n = 2` enters the prime sum, and
`g(log 2) + g(−log 2) = (1/π) ∫ |F|² cos(r log 2)` (Fourier inversion). Every `cos(r x)` with `|x| > 2a`
integrates to zero against `|F|²` (the autocorrelation vanishes there). So for any `b ≥ a` and free terms
`(d_k, x_k)` with `|x_k| > 2a`:

  `Re Q(f) = (1/2π) ∫ |F(r)|² W3(r) dr`,  `W3 = 2πμ + J(b, ·) − √2 log 2 · cos(r log 2) + Σ d_k cos(r x_k)`.

The free terms do not change `Q`; they reshape the weight so that the moment certificate (`Affine.lean`)
works on windows beyond `log 2`. The chord bounds add `|d| x²` per cosine term to the second-order loss.
-/
import Riemann.Deep.Yoshida.Affine

noncomputable section

open Complex MeasureTheory Set Filter Finset
open scoped ComplexConjugate

namespace Fayed.Yoshida.Beyond

open Riemann.Deep.Weil Riemann.Deep.WeilEF Zeta23 Zeta23.EF Fayed.Yoshida

/-! ### The weight -/

/-- `c₂ = √2 log 2`: the prime `2` enters the weight as `−c₂ cos(r log 2)`. -/
def cP : ℝ := Real.sqrt 2 * Real.log 2

lemma cP_nonneg : 0 ≤ cP := mul_nonneg (Real.sqrt_nonneg _) (Real.log_nonneg (by norm_num))

/-- Free cosine terms `Σ d cos(r x)` over pairs `(d, x)`. -/
def trigF : List (ℝ × ℝ) → ℝ → ℝ
  | [], _ => 0
  | p :: ps, r => p.1 * Real.cos (r * p.2) + trigF ps r

/-- The trigonometric part of the weight. -/
def trigP (fs : List (ℝ × ℝ)) (r : ℝ) : ℝ := -cP * Real.cos (r * Real.log 2) + trigF fs r

/-- The weight `W3 = 2πμ + J(b, ·) + trigP`. -/
def W3 (b : ℝ) (fs : List (ℝ × ℝ)) (r : ℝ) : ℝ := W b r + trigP fs r

/-- `Σ |d| x²` (curvature of the free terms). -/
def absSum : List (ℝ × ℝ) → ℝ
  | [] => 0
  | p :: ps => |p.1| * p.2 ^ 2 + absSum ps

/-- `Σ |d|`. -/
def absSum0 : List (ℝ × ℝ) → ℝ
  | [] => 0
  | p :: ps => |p.1| + absSum0 ps

/-- Curvature bound of the trigonometric part. -/
def curvP (fs : List (ℝ × ℝ)) : ℝ := cP * Real.log 2 ^ 2 + absSum fs

lemma absSum_nonneg : ∀ fs : List (ℝ × ℝ), 0 ≤ absSum fs
  | [] => le_refl _
  | p :: ps => by
    simp only [absSum]
    have := absSum_nonneg ps
    positivity

lemma curvP_nonneg (fs : List (ℝ × ℝ)) : 0 ≤ curvP fs := by
  unfold curvP
  have := cP_nonneg
  have := absSum_nonneg fs
  positivity

lemma trigF_neg : ∀ (fs : List (ℝ × ℝ)) (r : ℝ), trigF fs (-r) = trigF fs r
  | [], _ => rfl
  | p :: ps, r => by
    simp only [trigF, neg_mul, Real.cos_neg, trigF_neg ps r]

lemma trigP_neg (fs : List (ℝ × ℝ)) (r : ℝ) : trigP fs (-r) = trigP fs r := by
  unfold trigP
  rw [trigF_neg, neg_mul r, Real.cos_neg]

lemma trigP_abs (fs : List (ℝ × ℝ)) (r : ℝ) : trigP fs |r| = trigP fs r := by
  rcases le_total 0 r with h | h
  · rw [abs_of_nonneg h]
  · rw [abs_of_nonpos h, trigP_neg]

lemma abs_trigF_le : ∀ (fs : List (ℝ × ℝ)) (r : ℝ), |trigF fs r| ≤ absSum0 fs
  | [], _ => by simp [trigF, absSum0]
  | p :: ps, r => by
    simp only [trigF, absSum0]
    have h1 : |p.1 * Real.cos (r * p.2)| ≤ |p.1| := by
      rw [abs_mul]
      exact mul_le_of_le_one_right (abs_nonneg _) (Real.abs_cos_le_one _)
    have h2 := abs_trigF_le ps r
    exact (abs_add_le _ _).trans (by linarith)

lemma continuous_trigF : ∀ fs : List (ℝ × ℝ), Continuous (trigF fs)
  | [] => continuous_const
  | p :: ps => by
    have h : trigF (p :: ps) = fun r => p.1 * Real.cos (r * p.2) + trigF ps r := rfl
    rw [h]
    exact (continuous_const.mul (Real.continuous_cos.comp (continuous_id.mul continuous_const))).add
      (continuous_trigF ps)

/-! ### The Fourier form beyond `log 2` -/

lemma integrable_phi {a : ℝ} {f : ℝ → ℂ} (hf : TestFn a f) :
    Integrable (fun r : ℝ => ‖paperFT f r‖ ^ 2) := by
  have hk := hf.autocorr_contDiff
  have hkc := hf.autocorr_hasCompactSupport
  have hFk := integrable_fourier_of_contDiff_two hk hkc
  refine (integrable_paperFT_ofReal hFk).re.congr (Eventually.of_forall fun r => ?_)
  simp only [paperFT_autocorr_real hf, RCLike.re_to_complex, Complex.ofReal_re]

/-- `g(y) + g(−y) = (1/π) ∫ |F|² cos(r y)`. -/
lemma autocorr_add_neg {a : ℝ} {f : ℝ → ℂ} (hf : TestFn a f) (y : ℝ) :
    autocorr f y + autocorr f (-y)
      = (((1 / Real.pi) * ∫ r : ℝ, ‖paperFT f r‖ ^ 2 * Real.cos (r * y) : ℝ) : ℂ) := by
  have hk := hf.autocorr_contDiff
  have hkc := hf.autocorr_hasCompactSupport
  have hFk := integrable_fourier_of_contDiff_two hk hkc
  rw [k_add_k_neg hk.continuous (hk.continuous.integrable_of_hasCompactSupport hkc) hFk y]
  have e : (fun τ : ℝ => paperFT (autocorr f) τ * ((Real.cos (τ * y) : ℝ) : ℂ))
      = fun τ : ℝ => ((‖paperFT f τ‖ ^ 2 * Real.cos (τ * y) : ℝ) : ℂ) := by
    funext τ
    rw [paperFT_autocorr_real hf]
    push_cast
    ring
  rw [e, integral_complex_ofReal]
  push_cast
  ring

/-- Cosines beyond the support integrate to zero. -/
lemma integral_cos_eq_zero {a : ℝ} {f : ℝ → ℂ} (hf : TestFn a f) {x : ℝ} (hx : 2 * a < |x|) :
    ∫ r : ℝ, ‖paperFT f r‖ ^ 2 * Real.cos (r * x) = 0 := by
  have h := autocorr_add_neg hf x
  rw [autocorr_eq_zero hf.2 hx, autocorr_eq_zero hf.2 (by rwa [abs_neg]), add_zero] at h
  have h2 : (1 / Real.pi) * ∫ r : ℝ, ‖paperFT f r‖ ^ 2 * Real.cos (r * x) = 0 := by
    exact_mod_cast h.symm
  have hπ : (1 / Real.pi) ≠ 0 := by have := Real.pi_pos; positivity
  exact (mul_eq_zero.1 h2).resolve_left hπ

lemma integrable_phi_cos {a : ℝ} {f : ℝ → ℂ} (hf : TestFn a f) (x : ℝ) :
    Integrable (fun r : ℝ => ‖paperFT f r‖ ^ 2 * Real.cos (r * x)) :=
  (integrable_phi hf).mul_bdd (by fun_prop : Continuous fun r : ℝ => Real.cos (r * x)).aestronglyMeasurable
    (Eventually.of_forall fun r => by rw [Real.norm_eq_abs]; exact Real.abs_cos_le_one _)

lemma trigF_integral {a : ℝ} {f : ℝ → ℂ} (hf : TestFn a f) :
    ∀ fs : List (ℝ × ℝ), (∀ p ∈ fs, 2 * a < |p.2|) →
    Integrable (fun r : ℝ => ‖paperFT f r‖ ^ 2 * trigF fs r) ∧
    ∫ r : ℝ, ‖paperFT f r‖ ^ 2 * trigF fs r = 0
  | [], _ => by simp [trigF]
  | p :: ps, h => by
    simp only [List.mem_cons, forall_eq_or_imp] at h
    obtain ⟨ih1, ih2⟩ := trigF_integral hf ps h.2
    have hc : Integrable (fun r : ℝ => ‖paperFT f r‖ ^ 2 * (p.1 * Real.cos (r * p.2))) := by
      have := (integrable_phi_cos hf p.2).const_mul p.1
      refine this.congr (Eventually.of_forall fun r => ?_)
      simp only
      ring
    have e : (fun r : ℝ => ‖paperFT f r‖ ^ 2 * trigF (p :: ps) r)
        = fun r : ℝ => ‖paperFT f r‖ ^ 2 * (p.1 * Real.cos (r * p.2)) + ‖paperFT f r‖ ^ 2 * trigF ps r := by
      funext r
      simp only [trigF]
      ring
    refine ⟨by rw [e]; exact hc.add ih1, ?_⟩
    rw [e, integral_add hc ih1, ih2, add_zero]
    have e2 : (fun r : ℝ => ‖paperFT f r‖ ^ 2 * (p.1 * Real.cos (r * p.2)))
        = fun r : ℝ => p.1 * (‖paperFT f r‖ ^ 2 * Real.cos (r * p.2)) := by
      funext r
      ring
    rw [e2, integral_const_mul, integral_cos_eq_zero hf h.1, mul_zero]

/-- Below `log 3` only `n = 2` enters the prime sum. -/
lemma primeTerm_eq_two {a : ℝ} {f : ℝ → ℂ} (hf : Function.support f ⊆ Set.Icc (-a) a)
    (h3 : 2 * a < Real.log 3) :
    primeTerm (autocorr f) = ((Real.log 2 / Real.sqrt 2 : ℝ) : ℂ)
      * (autocorr f (Real.log 2) + autocorr f (-Real.log 2)) := by
  unfold primeTerm
  rw [tsum_eq_single 2]
  · rw [ArithmeticFunction.vonMangoldt_apply_prime Nat.prime_two]
    norm_num
  · intro n hn
    rcases Nat.lt_or_ge n 2 with h | h
    · interval_cases n <;> simp [ArithmeticFunction.vonMangoldt_apply_one]
    · have h3n : 3 ≤ n := by omega
      have hlog : Real.log 3 ≤ Real.log n :=
        Real.log_le_log (by norm_num) (by exact_mod_cast h3n)
      have hl0 : 0 ≤ Real.log n := le_trans (Real.log_nonneg (by norm_num)) hlog
      rw [autocorr_eq_zero hf (by rw [abs_of_nonneg hl0]; linarith),
        autocorr_eq_zero hf (by rw [abs_neg, abs_of_nonneg hl0]; linarith)]
      simp

/-- **The Fourier form beyond `log 2`** (any window `2a < log 3`, pole kernel cut at `b ≥ a`). -/
theorem weilQ_re_eq_W3 {a b : ℝ} (ha : 0 ≤ a) (hab : a ≤ b) (h3 : 2 * a < Real.log 3)
    (fs : List (ℝ × ℝ)) (hfs : ∀ p ∈ fs, 2 * a < |p.2|) {f : ℝ → ℂ} (hf : TestFn a f) :
    Integrable (fun r : ℝ => ‖paperFT f r‖ ^ 2 * W3 b fs r) ∧
    (weilQ f).re = 1 / (2 * Real.pi) * ∫ r : ℝ, ‖paperFT f r‖ ^ 2 * W3 b fs r := by
  have hb0 : 0 ≤ b := ha.trans hab
  have hφi := integrable_phi hf
  obtain ⟨hμi, harch⟩ := arch_re hf
  have hJi : Integrable (fun r : ℝ => ‖paperFT f r‖ ^ 2 * J b r) :=
    hφi.mul_bdd (continuous_J b).aestronglyMeasurable
      (Eventually.of_forall fun r => by rw [Real.norm_eq_abs]; exact abs_J_le hb0 r)
  have hμi' : Integrable (fun r : ℝ => ‖paperFT f r‖ ^ 2 * (2 * Real.pi * mu r)) := by
    refine (hμi.const_mul (2 * Real.pi)).congr (Eventually.of_forall fun r => ?_)
    simp only
    ring
  have hci := (integrable_phi_cos hf (Real.log 2)).const_mul (-cP)
  obtain ⟨hti, ht0⟩ := trigF_integral hf fs hfs
  have hsum : (fun r : ℝ => ‖paperFT f r‖ ^ 2 * W3 b fs r)
      = fun r : ℝ => ‖paperFT f r‖ ^ 2 * (2 * Real.pi * mu r) + ‖paperFT f r‖ ^ 2 * J b r
        + (-cP) * (‖paperFT f r‖ ^ 2 * Real.cos (r * Real.log 2)) + ‖paperFT f r‖ ^ 2 * trigF fs r := by
    funext r
    unfold W3 trigP W
    ring
  have i1 : Integrable (fun r : ℝ => ‖paperFT f r‖ ^ 2 * (2 * Real.pi * mu r) + ‖paperFT f r‖ ^ 2 * J b r) :=
    hμi'.add hJi
  have i2 : Integrable (fun r : ℝ => ‖paperFT f r‖ ^ 2 * (2 * Real.pi * mu r) + ‖paperFT f r‖ ^ 2 * J b r
      + (-cP) * (‖paperFT f r‖ ^ 2 * Real.cos (r * Real.log 2))) := i1.add hci
  refine ⟨by rw [hsum]; exact i2.add hti, ?_⟩
  have hp := primeTerm_eq_two hf.2 h3
  rw [autocorr_add_neg hf, ← Complex.ofReal_mul] at hp
  rw [weilQ_eq_fix]
  unfold weilQFix
  simp only [Complex.add_re, Complex.sub_re]
  rw [hp, harch, pole_re_eq_J hb0 (hf.mono hab), Complex.ofReal_re, hsum,
    integral_add i2 hti, integral_add i1 hci, integral_add hμi' hJi,
    integral_const_mul, ht0]
  have hπ := Real.pi_pos
  have e : ∫ r : ℝ, ‖paperFT f r‖ ^ 2 * (2 * Real.pi * mu r)
      = 2 * Real.pi * ∫ r : ℝ, ‖paperFT f r‖ ^ 2 * mu r := by
    rw [← integral_const_mul]
    congr 1
    funext r
    ring
  rw [e]
  have hs2 : Real.sqrt 2 ^ 2 = 2 := Real.sq_sqrt (by norm_num)
  have hs0 : 0 < Real.sqrt 2 := by positivity
  have hinv : Real.log 2 / Real.sqrt 2 = Real.sqrt 2 * Real.log 2 / 2 := by
    rw [div_eq_div_iff hs0.ne' (by norm_num)]
    linear_combination (-Real.log 2) * hs2
  rw [hinv]
  unfold cP
  field_simp
  ring

/-! ### Chord bounds for the trigonometric part -/

/-- `y ↦ −cos(yx) − x²y²/2` is concave. -/
lemma concave_negcos_sub (x : ℝ) :
    ConcaveOn ℝ Set.univ (fun y : ℝ => -Real.cos (y * x) - x ^ 2 * y ^ 2 / 2) := by
  have hd : ∀ y, HasDerivAt (fun y : ℝ => -Real.cos (y * x) - x ^ 2 * y ^ 2 / 2)
      (Real.sin (y * x) * x - x ^ 2 * y) y := by
    intro y
    have h1 : HasDerivAt (fun y : ℝ => -Real.cos (y * x)) (-(-Real.sin (y * x) * x)) y :=
      ((hasDerivAt_mul_const x).cos).neg
    have h2 : HasDerivAt (fun y : ℝ => x ^ 2 * y ^ 2 / 2) (x ^ 2 * y) y := by
      have := ((hasDerivAt_pow 2 y).const_mul (x ^ 2)).div_const 2
      convert this using 1
      push_cast
      ring
    convert h1.sub h2 using 1
    ring
  have hd2 : ∀ y, HasDerivAt (fun y : ℝ => Real.sin (y * x) * x - x ^ 2 * y)
      (Real.cos (y * x) * x * x - x ^ 2) y := by
    intro y
    have h1 : HasDerivAt (fun y : ℝ => Real.sin (y * x) * x) (Real.cos (y * x) * x * x) y :=
      (hasDerivAt_mul_const x).sin.mul_const x
    have h2 : HasDerivAt (fun y : ℝ => x ^ 2 * y) (x ^ 2) y := by
      simpa using (hasDerivAt_id y).const_mul (x ^ 2)
    exact h1.sub h2
  have hD : deriv (fun y : ℝ => -Real.cos (y * x) - x ^ 2 * y ^ 2 / 2)
      = fun y => Real.sin (y * x) * x - x ^ 2 * y := funext fun y => (hd y).deriv
  refine concaveOn_univ_of_deriv2_nonpos (fun y => (hd y).differentiableAt) ?_ ?_
  · rw [hD]; exact fun y => (hd2 y).differentiableAt
  · intro y
    rw [Function.iterate_succ, Function.iterate_one, Function.comp_apply, hD, (hd2 y).deriv]
    have := Real.cos_le_one (y * x)
    nlinarith [sq_nonneg x]

lemma negcos_chord (x lo hi θ : ℝ) (hθ0 : 0 ≤ θ) (hθ1 : θ ≤ 1) :
    (1 - θ) * (-Real.cos (lo * x)) + θ * (-Real.cos (hi * x))
        - x ^ 2 * (θ * (1 - θ) * (hi - lo) ^ 2) / 2
      ≤ -Real.cos (((1 - θ) * lo + θ * hi) * x) := by
  have h := (concave_negcos_sub x).2 (Set.mem_univ lo) (Set.mem_univ hi) (by linarith : 0 ≤ 1 - θ) hθ0
    (by ring)
  simp only [smul_eq_mul] at h
  nlinarith

/-- Chord bound for `d cos(rx)`, any sign of `d`. -/
lemma dcos_chord (d x lo hi θ : ℝ) (hθ0 : 0 ≤ θ) (hθ1 : θ ≤ 1) :
    (1 - θ) * (d * Real.cos (lo * x)) + θ * (d * Real.cos (hi * x))
        - |d| * x ^ 2 * (θ * (1 - θ) * (hi - lo) ^ 2) / 2
      ≤ d * Real.cos (((1 - θ) * lo + θ * hi) * x) := by
  rcases le_total 0 d with hd | hd
  · have h := cos_chord x lo hi θ hθ0 hθ1
    rw [abs_of_nonneg hd]
    have := mul_le_mul_of_nonneg_left h hd
    nlinarith
  · have h := negcos_chord x lo hi θ hθ0 hθ1
    rw [abs_of_nonpos hd]
    have := mul_le_mul_of_nonneg_left h (neg_nonneg.2 hd)
    nlinarith

lemma trigF_chord (lo hi θ : ℝ) (hθ0 : 0 ≤ θ) (hθ1 : θ ≤ 1) : ∀ fs : List (ℝ × ℝ),
    (1 - θ) * trigF fs lo + θ * trigF fs hi - absSum fs * (θ * (1 - θ) * (hi - lo) ^ 2) / 2
      ≤ trigF fs ((1 - θ) * lo + θ * hi)
  | [] => by simp [trigF, absSum]
  | p :: ps => by
    have h1 := dcos_chord p.1 p.2 lo hi θ hθ0 hθ1
    have h2 := trigF_chord lo hi θ hθ0 hθ1 ps
    simp only [trigF, absSum]
    nlinarith

lemma trigP_chord (fs : List (ℝ × ℝ)) (lo hi θ : ℝ) (hθ0 : 0 ≤ θ) (hθ1 : θ ≤ 1) :
    (1 - θ) * trigP fs lo + θ * trigP fs hi - curvP fs * (θ * (1 - θ) * (hi - lo) ^ 2) / 2
      ≤ trigP fs ((1 - θ) * lo + θ * hi) := by
  have h1 := dcos_chord (-cP) (Real.log 2) lo hi θ hθ0 hθ1
  have h2 := trigF_chord lo hi θ hθ0 hθ1 fs
  rw [abs_neg, abs_of_nonneg cP_nonneg] at h1
  unfold trigP curvP
  nlinarith

/-! ### Certified affine cells for `W3` -/

/-- A certified affine cell: `W3 ≥ s − (A + B|r|)` on `lo ≤ |r| < hi`, weight `≥ 0` there. -/
def AffOK3 (b s : ℝ) (fs : List (ℝ × ℝ)) (c : Cell4) : Prop :=
  0 ≤ c.1 ∧ c.1 ≤ c.2.1 ∧ 0 ≤ c.2.2.1 + c.2.2.2 * c.1 ∧ 0 ≤ c.2.2.1 + c.2.2.2 * c.2.1 ∧
    ∀ r : ℝ, c.1 ≤ |r| → |r| < c.2.1 → s - (c.2.2.1 + c.2.2.2 * |r|) ≤ W3 b fs r

/-- **Affine cell for `W3` from rational end data** (as `affOK_of`, with the trigonometric part). -/
lemma affOK3_of {b s lo hi A B mL mH jL jH Ej : ℝ} {fs : List (ℝ × ℝ)} (hb : 0 ≤ b) (h0 : 0 ≤ lo)
    (hlh : lo < hi) (hmL : mL ≤ 2 * Real.pi * mu lo) (hmH : mH ≤ 2 * Real.pi * mu hi)
    (hjL : jL ≤ J b lo + trigP fs lo) (hjH : jH ≤ J b hi + trigP fs hi)
    (hE : 4 / 3 * b ^ 3 * Real.cosh b + curvP fs / 8 ≤ Ej) (hΔ : mL ≤ mH)
    (hn1 : 0 ≤ A + B * lo) (hn2 : 0 ≤ A + B * hi)
    (hc1 : s - (A + B * lo) ≤ mL + jL - ((mH - mL) * ((hi - lo) / (4 * (hi + lo))) + Ej * (hi - lo) ^ 2))
    (hc2 : s - (A + B * hi) ≤ mH + jH - ((mH - mL) * ((hi - lo) / (4 * (hi + lo))) + Ej * (hi - lo) ^ 2)) :
    AffOK3 b s fs (lo, hi, A, B) := by
  refine ⟨h0, hlh.le, hn1, hn2, fun r h1 h2 => ?_⟩
  simp only at h1 h2 ⊢
  set x := |r| with hxdef
  have hJr : J b r = J b x := by
    rcases le_total 0 r with h | h
    · rw [hxdef, abs_of_nonneg h]
    · rw [hxdef, abs_of_nonpos h, J_even]
  have hW : W3 b fs r = 2 * Real.pi * mu x + (J b x + trigP fs x) := by
    unfold W3 W
    rw [hJr, hxdef, two_pi_mu_abs, trigP_abs]
    ring
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
  -- J and trigonometric part
  have hJ := J_chord (lo := lo) (hi := hi) hb hl0 hl1
  have hT := trigP_chord fs lo hi l hl0 hl1
  rw [← hxl] at hJ hT
  have hll : l * (1 - l) ≤ 1 / 4 := by nlinarith only [sq_nonneg (l - 1 / 2)]
  have hc0 : 0 ≤ b ^ 3 * Real.cosh b := by have := Real.cosh_pos b; positivity
  have hcv := curvP_nonneg fs
  have k3 : (16 / 3 * b ^ 3 * Real.cosh b + curvP fs / 2) * (l * (1 - l) * (hi - lo) ^ 2)
      ≤ Ej * (hi - lo) ^ 2 := by
    have hw2 : 0 ≤ (hi - lo) ^ 2 := sq_nonneg _
    have h16 : 0 ≤ 16 / 3 * b ^ 3 * Real.cosh b + curvP fs / 2 := by positivity
    have t1 := mul_le_mul_of_nonneg_left hll h16
    have t2 : (16 / 3 * b ^ 3 * Real.cosh b + curvP fs / 2) * (l * (1 - l)) ≤ Ej := by
      linarith only [t1, hE]
    have t3 := mul_le_mul_of_nonneg_right t2 hw2
    have e4 : (16 / 3 * b ^ 3 * Real.cosh b + curvP fs / 2) * (l * (1 - l) * (hi - lo) ^ 2)
        = (16 / 3 * b ^ 3 * Real.cosh b + curvP fs / 2) * (l * (1 - l)) * (hi - lo) ^ 2 := by ring
    rw [e4]
    exact t3
  have k4 : (1 - l) * jL + l * jH ≤ (1 - l) * (J b lo + trigP fs lo) + l * (J b hi + trigP fs hi) := by
    have := mul_le_mul_of_nonneg_left hjL (by linarith : (0 : ℝ) ≤ 1 - l)
    have := mul_le_mul_of_nonneg_left hjH hl0
    linarith
  have k5 := mul_le_mul_of_nonneg_left hc1 (by linarith : (0 : ℝ) ≤ 1 - l)
  have k6 := mul_le_mul_of_nonneg_left hc2 hl0
  have e : s - (A + B * x) = (1 - l) * (s - (A + B * lo)) + l * (s - (A + B * hi)) := by
    rw [hxl]; ring
  have e3 : (1 - l) * (mL + jL - ((mH - mL) * q + Ej * (hi - lo) ^ 2))
      + l * (mH + jH - ((mH - mL) * q + Ej * (hi - lo) ^ 2))
      = ((1 - l) * mL + l * mH - (mH - mL) * q) + ((1 - l) * jL + l * jH) - Ej * (hi - lo) ^ 2 := by
    ring
  have e5 : (16 / 3 * b ^ 3 * Real.cosh b + curvP fs / 2) * (l * (1 - l) * (hi - lo) ^ 2)
      = 16 / 3 * b ^ 3 * Real.cosh b * (l * (1 - l) * (hi - lo) ^ 2)
        + curvP fs * (l * (1 - l) * (hi - lo) ^ 2) / 2 := by ring
  rw [e]
  linarith only [k2, hJ, hT, k3, k4, k5, k6, e3, e5]

/-! ### The certificate theorem beyond `log 2` -/

/-- **Affine moment certificate for `W3` ⇒ Weil positivity on every window `L ≤ 2a₁`**
(`2a₁ < log 3`; pole kernel cut at `b ≥ a₁`; free frequencies `|x| > 2a₁`). -/
theorem weilPos_of_affine3 {a1 b s R ε α β : ℝ} (ha1 : 0 < a1) (hb : a1 ≤ b) (h3 : 2 * a1 < Real.log 3)
    (fs : List (ℝ × ℝ)) (hfs : ∀ p ∈ fs, 2 * a1 < |p.2|)
    (cs : List Cell4) (hchain : Chain4 0 cs R) (hok : ∀ c ∈ cs, AffOK3 b s fs c)
    (htail : ∀ r : ℝ, R ≤ |r| → s ≤ W3 b fs r) (K : ℕ)
    (hK : ∀ c ∈ cs, (c.2.2.1 = 0 ∧ c.2.2.2 = 0) ∨ c.2.1 * a1 ≤ (K + 1) / 2) (Λ : ℕ → ℕ → ℝ)
    (hpsd : ∀ x : ℕ → ℝ, ∑ m ∈ range K, ∑ k ∈ range K, x m * x k * affSum cs (m + k)
      ≤ α * (2 * ∑ m ∈ range K, (∑ i ∈ range K, Λ m i * x i) * x m
        - ∑ m ∈ range K, ∑ k ∈ range K,
            ((∑ i ∈ range K, Λ m i * x i) * (∑ i ∈ range K, Λ k i * x i)) * Ghat a1 m k))
    (hε : 0 < ε) (hα : 0 ≤ α)
    (hβ : 8 * a1 ^ (2 * K + 1) * affSum cs (2 * K) / (K.factorial : ℝ) ^ 2 ≤ β)
    (hnum : (1 + ε) * α + (1 + 1 / ε) * β ≤ 2 * Real.pi * s)
    {L : ℝ} (hL1 : L ≤ 2 * a1) : WeilPos L := by
  have key : ∀ L : ℝ, 0 ≤ L → L ≤ 2 * a1 → WeilPos L := by
    intro L hL0 hL1 f hf
    set a := L / 2 with hadef
    have ha0 : 0 ≤ a := by rw [hadef]; linarith
    have haa : a ≤ a1 := by rw [hadef]; linarith
    have hf1 : TestFn a1 f := hf.mono haa
    obtain ⟨hφw, hre⟩ := weilQ_re_eq_W3 ha0 (haa.trans hb) (by linarith) fs
      (fun p hp => by linarith [hfs p hp]) hf
    set N := ∫ u, ‖f u‖ ^ 2 with hN
    have hN0 : 0 ≤ N := integral_nonneg fun u => by positivity
    have hφi := integrable_phi hf1
    have hcs : ∀ c ∈ cs, 0 ≤ c.1 ∧ c.1 ≤ c.2.1 ∧ 0 ≤ c.2.2.1 + c.2.2.2 * c.1
        ∧ 0 ≤ c.2.2.1 + c.2.2.2 * c.2.1 :=
      fun c hc => ⟨(hok c hc).1, (hok c hc).2.1, (hok c hc).2.2.1, (hok c hc).2.2.2.1⟩
    have hcs' : ∀ c ∈ cs, 0 ≤ c.1 ∧ 0 ≤ c.2.2.1 + c.2.2.2 * c.1 ∧ 0 ≤ c.2.2.1 + c.2.2.2 * c.2.1 :=
      fun c hc => ⟨(hcs c hc).1, (hcs c hc).2.2.1, (hcs c hc).2.2.2⟩
    have hw : ∀ r : ℝ, s - affStep cs r ≤ W3 b fs r := by
      intro r
      rcases le_or_gt R |r| with hR | hR
      · have := affStep_nonneg cs (fun _ _ => Or.inr trivial) hcs' r
        linarith [htail r hR]
      · obtain ⟨c, hc, h1, h2⟩ := chain4_cover cs 0 R hchain |r| (abs_nonneg r) hR
        have := le_affStep cs hcs' c hc r h1 h2
        linarith [(hok c hc).2.2.2.2 r h1 h2]
    obtain ⟨hsi, hsb⟩ := aff_step_bound ha1 hf1 hφi K hε cs hcs hK
    have hlow : s * (2 * Real.pi * N) - ∫ r : ℝ, ‖paperFT f r‖ ^ 2 * affStep cs r
        ≤ ∫ r : ℝ, ‖paperFT f r‖ ^ 2 * W3 b fs r := by
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
      have hb' := bessel ha1 hf1 K (fun m => ∑ i ∈ range K, Λ m i * x i)
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
      rw [e2] at hb'
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
    have hfin : 0 ≤ ∫ r : ℝ, ‖paperFT f r‖ ^ 2 * W3 b fs r := by
      have h3 : ((1 + ε) * α + (1 + 1 / ε) * β) * N ≤ (2 * Real.pi * s) * N :=
        mul_le_mul_of_nonneg_right hnum hN0
      nlinarith
    rw [hre]
    have hπ := Real.pi_pos
    positivity
  rcases le_or_gt 0 L with h | h
  · exact key L h hL1
  · exact WeilPos.mono h.le (key 0 le_rfl (by linarith))

end Fayed.Yoshida.Beyond
