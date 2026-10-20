/-
Riemann/Deep/Yoshida/Small.lean — Weil positivity on a small window, sorry-free.

For `f` supported in `[−a, a]` with `2a < log 2` no prime enters, and (`F = paperFT f`, `N = ∫|f|²`)

  Re Q(f) = 2 Re(F(i/2) · conj F(−i/2)) + ∫ |F(r)|² μ(r) dr,   μ = (Re ψ(1/4 + ir/2) − log π)/2π.

Proof of `Q ≥ 0` on the window `1/16` ("bathtub" / layer-cake argument):
* pole term `≥ −‖F(i/2) − F(−i/2)‖² ≥ −(2a ∫|f|)² ≥ −8a³N`;
* Plancherel `∫ |F|² = 2πN` (Fourier inversion of the autocorrelation at 0);
* `|F(r)|² ≤ (∫|f|)² ≤ 2aN` (Cauchy–Schwarz via AM–GM), so `∫_{|r|<R} |F|² ≤ 4aRN`;
* `μ` is even and increasing in `|r|` with `μ ≥ μ(0) > −1`, and `μ(R) ≥ (log(R/2π) − 20/R²)/2π`
  (vertical Stirling), proved in the ported zeta23 `GammaFacts`;
* a 4-step lower bound for `μ` then gives `∫|F|²μ ≥ v₄·2πN − 2aN·2Σ(v_k − v_{k−1})R_k`.
-/
import Riemann.Deep.WeilEF.Suff

noncomputable section

open Complex MeasureTheory Set Filter
open scoped ComplexConjugate

namespace Fayed.Yoshida

open Riemann.Deep.Weil Riemann.Deep.WeilEF Zeta23 Zeta23.EF

/-! ### Fourier side -/

lemma norm_paperFT_le (f : ℝ → ℂ) (r : ℝ) : ‖paperFT f r‖ ≤ ∫ u, ‖f u‖ := by
  unfold paperFT
  refine (norm_integral_le_integral_norm _).trans (le_of_eq ?_)
  congr 1
  funext u
  rw [norm_mul, show I * (r : ℂ) * (u : ℂ) = ((r * u : ℝ) : ℂ) * I by push_cast; ring,
    Complex.norm_exp_ofReal_mul_I, mul_one]

/-- `ĝ(r) = |F(r)|²` on the real line. -/
lemma paperFT_autocorr_real {a : ℝ} {f : ℝ → ℂ} (hf : TestFn a f) (r : ℝ) :
    paperFT (autocorr f) r = ((‖paperFT f r‖ ^ 2 : ℝ) : ℂ) := by
  rw [paperFT_autocorr hf, Complex.conj_ofReal, Complex.mul_conj, Complex.normSq_eq_norm_sq]

lemma autocorr_zero (f : ℝ → ℂ) : autocorr f 0 = ((∫ u, ‖f u‖ ^ 2 : ℝ) : ℂ) := by
  rw [← integral_complex_ofReal]
  unfold autocorr
  congr 1
  funext y
  rw [sub_zero, Complex.mul_conj, Complex.normSq_eq_norm_sq]

/-- Plancherel: `∫ |F|² = 2π ∫ |f|²`. -/
lemma plancherel {a : ℝ} {f : ℝ → ℂ} (hf : TestFn a f) :
    ∫ r : ℝ, ‖paperFT f r‖ ^ 2 = 2 * Real.pi * ∫ u, ‖f u‖ ^ 2 := by
  have hk := hf.autocorr_contDiff
  have hkc := hf.autocorr_hasCompactSupport
  have hFk := integrable_fourier_of_contDiff_two hk hkc
  have hinv := paper_inversion hk.continuous (hk.continuous.integrable_of_hasCompactSupport hkc)
    hFk 0
  simp only [Complex.ofReal_zero, mul_zero, Complex.exp_zero, mul_one] at hinv
  simp only [paperFT_autocorr_real hf, autocorr_zero, integral_complex_ofReal] at hinv
  set Φ := ∫ r : ℝ, ‖paperFT f r‖ ^ 2
  set N := ∫ u, ‖f u‖ ^ 2
  have hπ : (Real.pi : ℂ) ≠ 0 := by exact_mod_cast Real.pi_ne_zero
  have h : (Φ : ℂ) = ((2 * Real.pi * N : ℝ) : ℂ) := by
    push_cast
    rw [hinv]
    field_simp
  exact_mod_cast h

/-- `(∫|f|)² ≤ 2a ∫|f|²` for `f` supported in `[−a, a]` (Cauchy–Schwarz via AM–GM). -/
lemma sq_integral_norm_le {a : ℝ} (ha : 0 < a) {f : ℝ → ℂ} (hf : TestFn a f) :
    (∫ u, ‖f u‖) ^ 2 ≤ 2 * a * ∫ u, ‖f u‖ ^ 2 := by
  have hc := hf.continuous
  have hcs := hf.hasCompactSupport
  have hzero : ∀ u, u ∉ Icc (-a) a → f u = 0 := fun u hu => by
    by_contra h
    exact hu (hf.2 (Function.mem_support.mpr h))
  have hI1 : Integrable (fun u => ‖f u‖) := hc.norm.integrable_of_hasCompactSupport hcs.norm
  have hcs2 : HasCompactSupport (fun u => ‖f u‖ ^ 2) :=
    HasCompactSupport.intro (isCompact_Icc (a := -a) (b := a))
      (fun u (hu : u ∉ Icc (-a) a) => by simp [hzero u hu])
  have hI2 : Integrable (fun u => ‖f u‖ ^ 2) :=
    (hc.norm.pow 2).integrable_of_hasCompactSupport hcs2
  have hA' : ∫ u, ‖f u‖ = ∫ u in Icc (-a) a, ‖f u‖ :=
    (setIntegral_eq_integral_of_forall_compl_eq_zero (fun u hu => by simp [hzero u hu])).symm
  have hN' : ∫ u, ‖f u‖ ^ 2 = ∫ u in Icc (-a) a, ‖f u‖ ^ 2 :=
    (setIntegral_eq_integral_of_forall_compl_eq_zero (fun u hu => by simp [hzero u hu])).symm
  have key : ∀ t : ℝ, 0 < t →
      ∫ u, ‖f u‖ ≤ (t * (∫ u, ‖f u‖ ^ 2) + t⁻¹ * (2 * a)) / 2 := by
    intro t ht
    have ht' := ht.ne'
    rw [hA', hN']
    have hint : IntegrableOn (fun u => (t * ‖f u‖ ^ 2 + t⁻¹) / 2) (Icc (-a) a) :=
      ((hI2.integrableOn.const_mul t).add continuous_const.integrableOn_Icc).div_const 2
    calc ∫ u in Icc (-a) a, ‖f u‖ ≤ ∫ u in Icc (-a) a, (t * ‖f u‖ ^ 2 + t⁻¹) / 2 := by
          refine setIntegral_mono_on hI1.integrableOn hint measurableSet_Icc (fun u _ => ?_)
          have e : t * ‖f u‖ ^ 2 + t⁻¹ - 2 * ‖f u‖ = t⁻¹ * (t * ‖f u‖ - 1) ^ 2 := by
            field_simp
            ring
          have : 0 ≤ t⁻¹ * (t * ‖f u‖ - 1) ^ 2 := by positivity
          linarith
      _ = (t * (∫ u in Icc (-a) a, ‖f u‖ ^ 2) + t⁻¹ * (2 * a)) / 2 := by
          rw [integral_div, integral_add (hI2.integrableOn.const_mul t)
            continuous_const.integrableOn_Icc, integral_const_mul, setIntegral_const,
            Real.volume_real_Icc_of_le (by linarith), smul_eq_mul]
          ring
  set A := ∫ u, ‖f u‖
  set N := ∫ u, ‖f u‖ ^ 2
  have hA0 : 0 ≤ A := integral_nonneg fun u => norm_nonneg _
  have hN0 : 0 ≤ N := integral_nonneg fun u => by positivity
  rcases eq_or_lt_of_le hA0 with h | h
  · rw [← h, zero_pow two_ne_zero]
    exact mul_nonneg (by linarith) hN0
  rcases eq_or_lt_of_le hN0 with hN | hN
  · exfalso
    have k := key (4 * a / A) (by positivity)
    rw [← hN, inv_div] at k
    have : A / (4 * a) * (2 * a) = A / 2 := by field_simp; ring
    rw [mul_zero, zero_add, this] at k
    linarith
  · have k := key (A / N) (by positivity)
    have h1 : A / N * N = A := div_mul_cancel₀ A hN.ne'
    have h2 : (A / N)⁻¹ * (2 * a) = 2 * a * N / A := by rw [inv_div]; ring
    rw [h1, h2] at k
    have k2 : A ≤ 2 * a * N / A := by linarith
    rw [le_div_iff₀ h] at k2
    nlinarith

/-! ### The pole term -/

lemma norm_exp_sub_le {a u : ℝ} (ha1 : a ≤ 1) (hu : u ∈ Icc (-a) a) :
    ‖cexp (I * (I / 2) * (u : ℂ)) - cexp (I * (-(I / 2)) * (u : ℂ))‖ ≤ 2 * a := by
  have h1 : I * (I / 2) * (u : ℂ) = ((-u / 2 : ℝ) : ℂ) := by
    push_cast
    linear_combination ((u : ℂ) / 2) * I_sq
  have h2 : I * (-(I / 2)) * (u : ℂ) = ((u / 2 : ℝ) : ℂ) := by
    push_cast
    linear_combination (-(u : ℂ) / 2) * I_sq
  rw [h1, h2, ← Complex.ofReal_exp, ← Complex.ofReal_exp, ← Complex.ofReal_sub, Complex.norm_real,
    Real.norm_eq_abs]
  have hu' : |u| ≤ a := abs_le.2 hu
  have e1 : |Real.exp (-u / 2) - 1| ≤ 2 * |-u / 2| :=
    Real.abs_exp_sub_one_le (by rw [abs_div, abs_neg]; norm_num; linarith)
  have e2 : |Real.exp (u / 2) - 1| ≤ 2 * |u / 2| :=
    Real.abs_exp_sub_one_le (by rw [abs_div]; norm_num; linarith)
  have e3 : |Real.exp (-u / 2) - Real.exp (u / 2)|
      ≤ |Real.exp (-u / 2) - 1| + |Real.exp (u / 2) - 1| := by
    have := abs_sub_le (Real.exp (-u / 2)) 1 (Real.exp (u / 2))
    rwa [abs_sub_comm 1] at this
  rw [abs_div, abs_neg, abs_two] at e1
  rw [abs_div, abs_two] at e2
  linarith

/-- `Re(pole term) = 2 Re(F(i/2) · conj F(−i/2))`. -/
lemma pole_re_eq {a : ℝ} {f : ℝ → ℂ} (hf : TestFn a f) :
    (poleTerm (autocorr f)).re = 2 * (paperFT f (I / 2) * conj (paperFT f (-(I / 2)))).re := by
  have hc1 : (starRingEnd ℂ) (I / 2) = -(I / 2) := by
    rw [map_div₀, Complex.conj_I, map_ofNat, neg_div]
  have hc2 : (starRingEnd ℂ) (-(I / 2)) = I / 2 := by rw [map_neg, hc1, neg_neg]
  have hpole : poleTerm (autocorr f) = paperFT f (I / 2) * conj (paperFT f (-(I / 2)))
      + paperFT f (-(I / 2)) * conj (paperFT f (I / 2)) := by
    unfold poleTerm
    rw [ft_eq_paperFT, ft_eq_paperFT, paperFT_autocorr hf, paperFT_autocorr hf, hc1, hc2]
  have hsym : (paperFT f (-(I / 2)) * conj (paperFT f (I / 2))).re
      = (paperFT f (I / 2) * conj (paperFT f (-(I / 2)))).re := by
    rw [show paperFT f (-(I / 2)) * conj (paperFT f (I / 2))
        = conj (paperFT f (I / 2) * conj (paperFT f (-(I / 2)))) by simp [mul_comm]]
    exact Complex.conj_re _
  rw [hpole, Complex.add_re, hsym]
  ring

/-- `Re(pole term) ≥ −(2a ∫|f|)²`. -/
lemma pole_re_ge {a : ℝ} (ha1 : a ≤ 1) {f : ℝ → ℂ} (hf : TestFn a f) :
    -(2 * a * ∫ u, ‖f u‖) ^ 2 ≤ (poleTerm (autocorr f)).re := by
  have hc := hf.continuous
  have hcs := hf.hasCompactSupport
  have hzero : ∀ u, u ∉ Icc (-a) a → f u = 0 := fun u hu => by
    by_contra h
    exact hu (hf.2 (Function.mem_support.mpr h))
  have hre : (poleTerm (autocorr f)).re = ‖paperFT f (I / 2)‖ ^ 2 + ‖paperFT f (-(I / 2))‖ ^ 2
      - ‖paperFT f (I / 2) - paperFT f (-(I / 2))‖ ^ 2 := by
    rw [pole_re_eq hf, ← Complex.normSq_eq_norm_sq, ← Complex.normSq_eq_norm_sq,
      ← Complex.normSq_eq_norm_sq, Complex.normSq_sub]
    ring
  have hI1 : Integrable (fun u => ‖f u‖) := hc.norm.integrable_of_hasCompactSupport hcs.norm
  have hPM : ‖paperFT f (I / 2) - paperFT f (-(I / 2))‖ ≤ 2 * a * ∫ u, ‖f u‖ := by
    unfold paperFT
    rw [← integral_sub (Fayed.WeilSuff.integrable_mul_cexp hc hcs _)
      (Fayed.WeilSuff.integrable_mul_cexp hc hcs _), ← integral_const_mul]
    refine norm_integral_le_of_norm_le (hI1.const_mul (2 * a)) (Eventually.of_forall fun u => ?_)
    by_cases hu : u ∈ Icc (-a) a
    · rw [← mul_sub, norm_mul, mul_comm (2 * a)]
      exact mul_le_mul_of_nonneg_left (norm_exp_sub_le ha1 hu) (norm_nonneg _)
    · simp [hzero u hu]
  have h2 := pow_le_pow_left₀ (norm_nonneg _) hPM 2
  rw [hre]
  nlinarith [sq_nonneg ‖paperFT f (I / 2)‖, sq_nonneg ‖paperFT f (-(I / 2))‖]

/-! ### The archimedean term -/

lemma gammaBracket_eq (r : ℝ) : gammaBracket r = 2 * Real.pi * mu r := by
  have hπ := Real.pi_pos
  unfold gammaBracket mu
  field_simp

lemma arch_re {a : ℝ} {f : ℝ → ℂ} (hf : TestFn a f) :
    Integrable (fun r : ℝ => ‖paperFT f r‖ ^ 2 * mu r) ∧
    (archTermFix (autocorr f)).re = ∫ r : ℝ, ‖paperFT f r‖ ^ 2 * mu r := by
  have hk := hf.autocorr_contDiff
  have hkc := hf.autocorr_hasCompactSupport
  have hmu := integrable_paperFT_mul_mu hk hkc gammaFacts
  have hpt : ∀ r : ℝ, paperFT (autocorr f) r * (mu r : ℂ)
      = ((‖paperFT f r‖ ^ 2 * mu r : ℝ) : ℂ) := by
    intro r
    rw [paperFT_autocorr_real hf]
    push_cast
    ring
  refine ⟨?_, ?_⟩
  · refine hmu.re.congr (Eventually.of_forall fun r => ?_)
    simp only [hpt, RCLike.re_to_complex, Complex.ofReal_re]
  · have h := arch_eq hk hkc
    have e : archTermFix (autocorr f)
        = (1 / (2 * Real.pi) : ℂ) * ∫ r : ℝ, paperFT (autocorr f) r * (gammaBracket r : ℂ) := by
      unfold archTermFix
      exact h.symm
    rw [e]
    have hpt2 : ∀ r : ℝ, paperFT (autocorr f) r * (gammaBracket r : ℂ)
        = ((2 * Real.pi * (‖paperFT f r‖ ^ 2 * mu r) : ℝ) : ℂ) := by
      intro r
      rw [paperFT_autocorr_real hf, gammaBracket_eq]
      push_cast
      ring
    simp only [hpt2, integral_complex_ofReal]
    rw [integral_const_mul,
      show (1 / (2 * Real.pi) : ℂ) = ((1 / (2 * Real.pi) : ℝ) : ℂ) by norm_cast,
      ← Complex.ofReal_mul, Complex.ofReal_re, one_div, ← mul_assoc,
      inv_mul_cancel₀ (by positivity), one_mul]

/-! ### The bathtub (layer-cake) bound -/

lemma bathtub {φ w : ℝ → ℝ} {T c v0 v1 v2 v3 v4 r1 r2 r3 r4 : ℝ}
    (hφ0 : ∀ r, 0 ≤ φ r) (hφc : ∀ r, φ r ≤ c) (hφi : Integrable φ)
    (hφw : Integrable (fun r => φ r * w r)) (hT : ∫ r, φ r = T)
    (h0 : ∀ r, v0 ≤ w r) (h1 : ∀ r, r1 ≤ |r| → v1 ≤ w r) (h2 : ∀ r, r2 ≤ |r| → v2 ≤ w r)
    (h3 : ∀ r, r3 ≤ |r| → v3 ≤ w r) (h4 : ∀ r, r4 ≤ |r| → v4 ≤ w r)
    (h01 : v0 ≤ v1) (h12 : v1 ≤ v2) (h23 : v2 ≤ v3) (h34 : v3 ≤ v4)
    (hr1 : 0 ≤ r1) (hr12 : r1 ≤ r2) (hr23 : r2 ≤ r3) (hr34 : r3 ≤ r4) :
    v4 * T - c * (2 * ((v1 - v0) * r1 + (v2 - v1) * r2 + (v3 - v2) * r3 + (v4 - v3) * r4))
      ≤ ∫ r, φ r * w r := by
  set ind : ℝ → ℝ → ℝ := fun R => (Ioo (-R) R).indicator φ with hind
  have hI : ∀ R, Integrable (ind R) := fun R => hφi.indicator measurableSet_Ioo
  have hb : ∀ R, 0 ≤ R → ∫ r, ind R r ≤ c * (2 * R) := by
    intro R hR
    simp only [hind]
    rw [integral_indicator measurableSet_Ioo]
    calc ∫ r in Ioo (-R) R, φ r ≤ ∫ r in Ioo (-R) R, c :=
          setIntegral_mono_on hφi.integrableOn
            ((continuous_const.integrableOn_Icc).mono_set Ioo_subset_Icc_self)
            measurableSet_Ioo (fun r _ => hφc r)
      _ = c * (2 * R) := by
          rw [setIntegral_const, Real.volume_real_Ioo_of_le (by linarith), smul_eq_mul]
          ring
  have hpt : ∀ r, v4 * φ r - ((v1 - v0) * ind r1 r + (v2 - v1) * ind r2 r
      + (v3 - v2) * ind r3 r + (v4 - v3) * ind r4 r) ≤ φ r * w r := by
    intro r
    have hp := hφ0 r
    simp only [hind, Set.indicator_apply, Set.mem_Ioo, ← abs_lt]
    rcases lt_or_ge |r| r1 with q1 | q1
    · rw [ite_eq_left q1, ite_eq_left (by linarith), ite_eq_left (by linarith), ite_eq_left (by linarith)]
      nlinarith [mul_le_mul_of_nonneg_left (h0 r) hp]
    rcases lt_or_ge |r| r2 with q2 | q2
    · rw [ite_eq_right (not_lt.2 q1), ite_eq_left q2, ite_eq_left (by linarith), ite_eq_left (by linarith)]
      nlinarith [mul_le_mul_of_nonneg_left (h1 r q1) hp]
    rcases lt_or_ge |r| r3 with q3 | q3
    · rw [ite_eq_right (not_lt.2 q1), ite_eq_right (not_lt.2 q2), ite_eq_left q3, ite_eq_left (by linarith)]
      nlinarith [mul_le_mul_of_nonneg_left (h2 r q2) hp]
    rcases lt_or_ge |r| r4 with q4 | q4
    · rw [ite_eq_right (not_lt.2 q1), ite_eq_right (not_lt.2 q2), ite_eq_right (not_lt.2 q3), ite_eq_left q4]
      nlinarith [mul_le_mul_of_nonneg_left (h3 r q3) hp]
    · rw [ite_eq_right (not_lt.2 q1), ite_eq_right (not_lt.2 q2), ite_eq_right (not_lt.2 q3), ite_eq_right (not_lt.2 q4)]
      nlinarith [mul_le_mul_of_nonneg_left (h4 r q4) hp]
  have hS1 : Integrable (fun r => (v1 - v0) * ind r1 r + (v2 - v1) * ind r2 r) :=
    ((hI r1).const_mul (v1 - v0)).add ((hI r2).const_mul (v2 - v1))
  have hS2 : Integrable (fun r => (v1 - v0) * ind r1 r + (v2 - v1) * ind r2 r
      + (v3 - v2) * ind r3 r) := hS1.add ((hI r3).const_mul (v3 - v2))
  have hS3 : Integrable (fun r => (v1 - v0) * ind r1 r + (v2 - v1) * ind r2 r
      + (v3 - v2) * ind r3 r + (v4 - v3) * ind r4 r) := hS2.add ((hI r4).const_mul (v4 - v3))
  have hint : Integrable (fun r => v4 * φ r - ((v1 - v0) * ind r1 r + (v2 - v1) * ind r2 r
      + (v3 - v2) * ind r3 r + (v4 - v3) * ind r4 r)) := (hφi.const_mul v4).sub hS3
  have e : ∫ r, (v4 * φ r - ((v1 - v0) * ind r1 r + (v2 - v1) * ind r2 r
      + (v3 - v2) * ind r3 r + (v4 - v3) * ind r4 r))
      = v4 * T - ((v1 - v0) * (∫ r, ind r1 r) + (v2 - v1) * (∫ r, ind r2 r)
        + (v3 - v2) * (∫ r, ind r3 r) + (v4 - v3) * (∫ r, ind r4 r)) := by
    rw [integral_sub (hφi.const_mul v4) hS3, integral_add hS2 ((hI r4).const_mul _),
      integral_add hS1 ((hI r3).const_mul _), integral_add ((hI r1).const_mul _)
      ((hI r2).const_mul _)]
    simp only [integral_const_mul, hT]
  have hmono : ∫ r, (v4 * φ r - ((v1 - v0) * ind r1 r + (v2 - v1) * ind r2 r
      + (v3 - v2) * ind r3 r + (v4 - v3) * ind r4 r)) ≤ ∫ r, φ r * w r :=
    integral_mono hint hφw hpt
  rw [e] at hmono
  have b1 := mul_le_mul_of_nonneg_left (hb r1 hr1) (sub_nonneg.2 h01)
  have b2 := mul_le_mul_of_nonneg_left (hb r2 (by linarith)) (sub_nonneg.2 h12)
  have b3 := mul_le_mul_of_nonneg_left (hb r3 (by linarith)) (sub_nonneg.2 h23)
  have b4 := mul_le_mul_of_nonneg_left (hb r4 (by linarith)) (sub_nonneg.2 h34)
  nlinarith

/-! ### Core: positivity from a 4-step lower bound for `μ` -/

theorem core {a : ℝ} (ha : 0 < a) (ha1 : a ≤ 1) (hlog : 2 * a < Real.log 2)
    {v0 v1 v2 v3 v4 r1 r2 r3 r4 : ℝ}
    (h0 : ∀ r, v0 ≤ mu r) (h1 : ∀ r, r1 ≤ |r| → v1 ≤ mu r) (h2 : ∀ r, r2 ≤ |r| → v2 ≤ mu r)
    (h3 : ∀ r, r3 ≤ |r| → v3 ≤ mu r) (h4 : ∀ r, r4 ≤ |r| → v4 ≤ mu r)
    (h01 : v0 ≤ v1) (h12 : v1 ≤ v2) (h23 : v2 ≤ v3) (h34 : v3 ≤ v4)
    (hr1 : 0 ≤ r1) (hr12 : r1 ≤ r2) (hr23 : r2 ≤ r3) (hr34 : r3 ≤ r4)
    (hnum : 8 * a ^ 3 + 2 * a * (2 * ((v1 - v0) * r1 + (v2 - v1) * r2 + (v3 - v2) * r3
      + (v4 - v3) * r4)) ≤ v4 * (2 * Real.pi))
    {f : ℝ → ℂ} (hf : TestFn a f) : 0 ≤ (weilQ f).re := by
  have hk := hf.autocorr_contDiff
  have hkc := hf.autocorr_hasCompactSupport
  have hFk := integrable_fourier_of_contDiff_two hk hkc
  set N := ∫ u, ‖f u‖ ^ 2 with hN
  set A := ∫ u, ‖f u‖ with hA
  have hN0 : 0 ≤ N := integral_nonneg fun u => by positivity
  have hA2 : A ^ 2 ≤ 2 * a * N := sq_integral_norm_le ha hf
  set φ : ℝ → ℝ := fun r => ‖paperFT f r‖ ^ 2 with hφ
  have hφ0 : ∀ r, 0 ≤ φ r := fun r => by positivity
  have hφc : ∀ r, φ r ≤ 2 * a * N := fun r =>
    (pow_le_pow_left₀ (norm_nonneg _) (norm_paperFT_le f r) 2).trans hA2
  have hφi : Integrable φ := by
    refine (integrable_paperFT_ofReal hFk).re.congr (Eventually.of_forall fun r => ?_)
    simp only [hφ, paperFT_autocorr_real hf, RCLike.re_to_complex, Complex.ofReal_re]
  have hT : ∫ r, φ r = 2 * Real.pi * N := plancherel hf
  obtain ⟨hφw, harch⟩ := arch_re hf
  have hbath := bathtub (w := mu) hφ0 hφc hφi hφw hT h0 h1 h2 h3 h4 h01 h12 h23 h34
    hr1 hr12 hr23 hr34
  have hpole := pole_re_ge ha1 hf
  have hprime := primeTerm_eq_zero hf.2 hlog
  have hQ : (weilQ f).re = (poleTerm (autocorr f)).re + (archTermFix (autocorr f)).re
      - (primeTerm (autocorr f)).re := by
    rw [weilQ_eq_fix]
    unfold weilQFix
    simp only [Complex.add_re, Complex.sub_re]
  rw [hQ, hprime, harch, Complex.zero_re, sub_zero]
  have hp2 : (2 * a * A) ^ 2 ≤ 8 * a ^ 3 * N := by
    have : (2 * a * A) ^ 2 = 4 * a ^ 2 * A ^ 2 := by ring
    rw [this]
    nlinarith [mul_le_mul_of_nonneg_left hA2 (by positivity : (0 : ℝ) ≤ 4 * a ^ 2)]
  have hfin := mul_le_mul_of_nonneg_left hnum hN0
  nlinarith

/-! ### Lower bounds for `μ` -/

lemma mu_abs_mono {r0 r : ℝ} (h0 : 0 ≤ r0) (h : r0 ≤ |r|) : mu r0 ≤ mu r := by
  rcases le_total 0 r with hr | hr
  · rw [abs_of_nonneg hr] at h
    exact MuFields.mu_monotoneOn (mem_Ici.2 h0) (mem_Ici.2 hr) h
  · rw [abs_of_nonpos hr] at h
    rw [← mu_even r]
    exact MuFields.mu_monotoneOn (mem_Ici.2 h0) (mem_Ici.2 (by linarith)) h

lemma mu_ge_stirling {r0 q : ℝ} (h1 : 1 ≤ r0) (hq : Real.exp q * (2 * Real.pi) ≤ r0) :
    (q - 20 / r0 ^ 2) / (2 * Real.pi) ≤ mu r0 := by
  have hπ := Real.pi_pos
  have ht : (1 : ℝ) / 2 ≤ |r0 / 2| := by rw [abs_of_pos (by linarith)]; linarith
  have h := StirlingVert.re_digamma_stirling' (a := 1 / 4) (by norm_num) (by norm_num) ht
  rw [MuFields.mu_eq]
  have hD := (abs_le.1 h).1
  have hlog : q + Real.log Real.pi ≤ Real.log |r0 / 2| := by
    rw [abs_of_pos (by linarith), ← Real.log_exp q, ← Real.log_mul (Real.exp_pos q).ne' hπ.ne']
    exact Real.log_le_log (by positivity) (by linarith)
  have h20 : 5 / (r0 / 2) ^ 2 = 20 / r0 ^ 2 := by field_simp; ring
  rw [h20] at hD
  rw [show (1 / (2 * Real.pi)) * (Complex.digamma ((((1 : ℝ) / 4 : ℝ) : ℂ)
      + Complex.I * ((r0 / 2 : ℝ) : ℂ))).re - Real.log Real.pi / (2 * Real.pi)
      = ((Complex.digamma ((((1 : ℝ) / 4 : ℝ) : ℂ) + Complex.I * ((r0 / 2 : ℝ) : ℂ))).re
        - Real.log Real.pi) / (2 * Real.pi) by ring]
  apply div_le_div_of_nonneg_right _ (by positivity)
  linarith

lemma exp_half_lt : Real.exp (1 / 2) < 1.65 := by
  have h : Real.exp (1 / 2) ^ 2 = Real.exp 1 := by rw [← Real.exp_nat_mul]; norm_num
  nlinarith [Real.exp_one_lt_d9, Real.exp_pos (1 / 2)]

lemma exp_half_gt : 1.6 < Real.exp (1 / 2) := by
  have h : Real.exp (1 / 2) ^ 2 = Real.exp 1 := by rw [← Real.exp_nat_mul]; norm_num
  nlinarith [Real.exp_one_gt_d9, Real.exp_pos (1 / 2)]

/-! ### Main theorem -/

/-- **Weil positivity on the window `1/16`** (test functions supported in `[−1/32, 1/32]`). -/
theorem weilPos_sixteenth : WeilPos (1 / 16) := by
  intro f hf
  have hf' : TestFn (1 / 32) f := by
    have e : (1 / 16 : ℝ) / 2 = 1 / 32 := by norm_num
    rwa [e] at hf
  have hπ := Real.pi_pos
  have hπ1 := Real.pi_gt_d4
  have hπ2 := Real.pi_lt_d4
  have hlog : 2 * (1 / 32 : ℝ) < Real.log 2 := by
    have := Real.log_two_gt_d9
    linarith
  have hlo := exp_half_gt
  have hhi := exp_half_lt
  have he := Real.exp_one_lt_d9
  have hep := Real.exp_pos (1 / 2)
  -- step radii and levels
  have k1 : Real.exp (-(1 / 2)) * (2 * Real.pi) ≤ 4 := by
    rw [Real.exp_neg, inv_mul_le_iff₀ hep]
    nlinarith
  have k2 : Real.exp (1 / 2) * (2 * Real.pi) ≤ 21 / 2 := by nlinarith
  have k3 : Real.exp 1 * (2 * Real.pi) ≤ 69 / 4 := by nlinarith [Real.exp_pos 1]
  have hE : Real.exp (3 / 2) < 4.486 := by
    rw [show (3 / 2 : ℝ) = 1 + 1 / 2 by norm_num, Real.exp_add]
    nlinarith [Real.exp_pos 1]
  have k4 : Real.exp (3 / 2) * (2 * Real.pi) ≤ 113 / 4 := by
    nlinarith [Real.exp_pos (3 / 2)]
  have l1 := mu_ge_stirling (by norm_num) k1
  have l2 := mu_ge_stirling (by norm_num) k2
  have l3 := mu_ge_stirling (by norm_num) k3
  have l4 := mu_ge_stirling (by norm_num) k4
  have h2π : (0 : ℝ) < 2 * Real.pi := by positivity
  refine core (a := 1 / 32) (by norm_num) (by norm_num) hlog
    (v0 := -1) (v1 := (-(1 / 2) - 20 / 4 ^ 2) / (2 * Real.pi))
    (v2 := (1 / 2 - 20 / (21 / 2) ^ 2) / (2 * Real.pi))
    (v3 := (1 - 20 / (69 / 4) ^ 2) / (2 * Real.pi))
    (v4 := (3 / 2 - 20 / (113 / 4) ^ 2) / (2 * Real.pi))
    (r1 := 4) (r2 := 21 / 2) (r3 := 69 / 4) (r4 := 113 / 4)
    (fun r => (MuFields.neg_one_lt_mu_zero).le.trans (MuFields.mu_zero_le r))
    (fun r hr => l1.trans (mu_abs_mono (by norm_num) hr))
    (fun r hr => l2.trans (mu_abs_mono (by norm_num) hr))
    (fun r hr => l3.trans (mu_abs_mono (by norm_num) hr))
    (fun r hr => l4.trans (mu_abs_mono (by norm_num) hr))
    ?_ ?_ ?_ ?_ (by norm_num) (by norm_num) (by norm_num) (by norm_num) ?_ hf'
  · rw [le_div_iff₀ h2π]; nlinarith
  · exact div_le_div_of_nonneg_right (by norm_num) h2π.le
  · exact div_le_div_of_nonneg_right (by norm_num) h2π.le
  · exact div_le_div_of_nonneg_right (by norm_num) h2π.le
  · -- numeric inequality, linear in `1/(2π)`
    have hv4 : (3 / 2 - 20 / (113 / 4) ^ 2) / (2 * Real.pi) * (2 * Real.pi)
        = 3 / 2 - 20 / (113 / 4) ^ 2 := div_mul_cancel₀ _ h2π.ne'
    rw [hv4]
    set s := 1 / (2 * Real.pi) with hs
    have hs1 : s ≤ 1 / (2 * 3.1415) := by
      rw [hs]; exact one_div_le_one_div_of_le (by norm_num) (by linarith)
    have hs0 : 0 ≤ s := by positivity
    have ex : ∀ x : ℝ, x / (2 * Real.pi) = x * s := fun x => by rw [hs]; ring
    simp only [ex]
    nlinarith

/-! ### Prime-free form of the functional (interface for windows up to `log 2`) -/

/-- Below `log 2`: `Re Q(f) = 2 Re(F(i/2) · conj F(−i/2)) + ∫ |F(r)|² μ(r) dr`, `F = paperFT f`. -/
theorem weilQ_re_eq {a : ℝ} (hlog : 2 * a < Real.log 2) {f : ℝ → ℂ} (hf : TestFn a f) :
    (weilQ f).re = 2 * (paperFT f (I / 2) * conj (paperFT f (-(I / 2)))).re
      + ∫ r : ℝ, ‖paperFT f r‖ ^ 2 * mu r := by
  rw [weilQ_eq_fix]
  unfold weilQFix
  simp only [Complex.add_re, Complex.sub_re]
  rw [primeTerm_eq_zero hf.2 hlog, (arch_re hf).2, pole_re_eq hf, Complex.zero_re, sub_zero]

/-- Yoshida's theorem (window `log 2`, strict) is exactly positivity of this prime-free form. -/
theorem weilPos_iff_primeFree {L : ℝ} (hL : L < Real.log 2) :
    WeilPos L ↔ ∀ f : ℝ → ℂ, TestFn (L / 2) f →
      0 ≤ 2 * (paperFT f (I / 2) * conj (paperFT f (-(I / 2)))).re
        + ∫ r : ℝ, ‖paperFT f r‖ ^ 2 * mu r := by
  have h2 : 2 * (L / 2) < Real.log 2 := by linarith
  exact forall_congr' fun f => imp_congr_right fun hf => by rw [weilQ_re_eq h2 hf]

/-- **Yoshida-type small window**: Weil positivity for every window `0 < L ≤ 1/16`. -/
theorem weilPos_small (L : ℝ) (hL : L ≤ 1 / 16) : WeilPos L :=
  WeilPos.mono hL weilPos_sixteenth

end Fayed.Yoshida
