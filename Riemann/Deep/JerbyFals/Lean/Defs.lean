import Mathlib.Analysis.SpecialFunctions.Gamma.Digamma
import Mathlib.NumberTheory.LSeries.RiemannZeta
import Mathlib.MeasureTheory.Integral.IntervalIntegral.FundThmCalculus
import Mathlib.Analysis.Calculus.Deriv.MeanValue
import Mathlib.Analysis.Complex.RealDeriv
import Mathlib.Analysis.SpecialFunctions.ExpDeriv
import Mathlib.Analysis.Complex.CauchyIntegral
import Mathlib.Analysis.Calculus.FDeriv.Analytic

/-!
# Jerby's Gram-point repulsion conjecture (arXiv 2406.06548, conj:rep): definitions

* `thetaDeriv t = ½ Re ψ(1/4 + it/2) − ½ log π`, ψ = Γ′/Γ (`Complex.digamma`).
* `theta t = ∫₀ᵗ thetaDeriv`: the Riemann–Siegel θ, the continuous branch of
  arg Γ(1/4 + it/2) − (t/2) log π with θ(0) = 0 (Phase.lean proves the phase identity).
* `hardyZ t = Re (e^{iθ(t)} ζ(1/2 + it))`; `hardyZ'` is Arb's formula for Z′ (`hasDerivAt_hardyZ`).
* Gram points, good/bad, `JerbyConj` (Gram points of index n ± 1 exist and all are good) and
  `JerbyConjEx` (given Gram points of index n ± 1 good): Ex ⇒ JerbyConj; equal under `ThetaDerivPos`.
-/

open Complex Set
open scoped Real

namespace JerbyLean

noncomputable section

/-- θ′(t) = ½ Re ψ(1/4 + it/2) − ½ log π (certify_cx.py: `digamma().real/2 - log(pi)/2`). -/
def thetaDeriv (t : ℝ) : ℝ := (digamma (1 / 4 + t / 2 * I)).re / 2 - Real.log π / 2

/-- The Riemann–Siegel theta function θ(t) = ∫₀ᵗ θ′(u) du. -/
def theta (t : ℝ) : ℝ := ∫ u in (0 : ℝ)..t, thetaDeriv u

/-- Hardy's Z-function Z(t) = e^{iθ(t)} ζ(1/2 + it) (real part; the imaginary part is 0,
`hardyZ_eq` in Phase.lean). -/
def hardyZ (t : ℝ) : ℝ := (exp (I * theta t) * riemannZeta (1 / 2 + t * I)).re

/-- Z′(t) as certify_cx.py computes it: Re (e^{iθ} (iθ′ ζ(s) + i ζ′(s))), s = 1/2 + it. -/
def hardyZ' (t : ℝ) : ℝ :=
  (exp (I * theta t) * (I * thetaDeriv t * riemannZeta (1 / 2 + t * I) +
    I * deriv riemannZeta (1 / 2 + t * I))).re

/-! ## Calculus -/

/-- The half-plane Re s > 0. -/
def rhp : Set ℂ := {s | 0 < s.re}

lemma isOpen_rhp : IsOpen rhp := isOpen_lt continuous_const continuous_re

lemma ne_neg_nat_of_re_pos {s : ℂ} (hs : 0 < s.re) : ∀ m : ℕ, s ≠ -m := by
  intro m hm
  rw [hm, neg_re, natCast_re] at hs
  linarith [(Nat.cast_nonneg m : (0 : ℝ) ≤ m)]

lemma differentiableOn_Gamma_rhp : DifferentiableOn ℂ Gamma rhp := fun s hs =>
  (differentiableAt_Gamma s (ne_neg_nat_of_re_pos hs)).differentiableWithinAt

lemma continuousOn_digamma_rhp : ContinuousOn digamma rhp := by
  have hA : AnalyticOnNhd ℂ Gamma rhp := differentiableOn_Gamma_rhp.analyticOnNhd isOpen_rhp
  have h : ContinuousOn (fun s => deriv Gamma s / Gamma s) rhp :=
    hA.deriv.continuousOn.div differentiableOn_Gamma_rhp.continuousOn
      (fun s hs => Gamma_ne_zero_of_re_pos hs)
  have e : digamma = fun s => deriv Gamma s / Gamma s := by
    funext s; rw [digamma_def, logDeriv_apply]
  rw [e]; exact h

lemma re_quarter_line (t : ℝ) : (1 / 4 + t / 2 * I : ℂ).re = 1 / 4 := by simp

lemma continuous_thetaDeriv : Continuous thetaDeriv := by
  have hp : Continuous (fun t : ℝ => (1 / 4 + t / 2 * I : ℂ)) := by fun_prop
  have h := continuousOn_digamma_rhp.comp_continuous hp
    (fun t => by show 0 < (1 / 4 + (t : ℂ) / 2 * I).re; rw [re_quarter_line]; norm_num)
  exact ((continuous_re.comp h).div_const 2).sub continuous_const

lemma hasDerivAt_theta (t : ℝ) : HasDerivAt theta (thetaDeriv t) t :=
  intervalIntegral.integral_hasDerivAt_right (continuous_thetaDeriv.intervalIntegrable _ _)
    continuous_thetaDeriv.aestronglyMeasurable.stronglyMeasurableAtFilter
    continuous_thetaDeriv.continuousAt

lemma continuous_theta : Continuous theta :=
  continuous_iff_continuousAt.2 fun t => (hasDerivAt_theta t).continuousAt

lemma theta_zero : theta 0 = 0 := by simp [theta]

lemma hasDerivAt_zeta_line (t : ℝ) :
    HasDerivAt (fun u : ℝ => riemannZeta (1 / 2 + u * I))
      (deriv riemannZeta (1 / 2 + t * I) * I) t := by
  have hs : (1 / 2 + t * I : ℂ) ≠ 1 := by
    intro h
    have := congrArg re h
    norm_num at this
  have hl : HasDerivAt (fun u : ℝ => (1 / 2 + u * I : ℂ)) I t := by
    have := ((hasDerivAt_id t).ofReal_comp.mul_const I).const_add (1 / 2 : ℂ)
    simpa using this
  exact (differentiableAt_riemannZeta hs).hasDerivAt.comp t hl

/-- Z′ = `hardyZ'` (the derivative Arb computed is the derivative of Z). -/
lemma hasDerivAt_hardyZ (t : ℝ) : HasDerivAt hardyZ (hardyZ' t) t := by
  have h1 : HasDerivAt (fun u : ℝ => exp (I * theta u))
      (exp (I * theta t) * (I * thetaDeriv t)) t :=
    ((hasDerivAt_theta t).ofReal_comp.const_mul I).cexp
  have h2 := h1.mul (hasDerivAt_zeta_line t)
  have h3 := Complex.reCLM.hasFDerivAt.comp_hasDerivAt t h2
  simp only [Function.comp_def, Complex.reCLM_apply] at h3
  have e : hardyZ' t = (exp (I * theta t) * (I * thetaDeriv t) * riemannZeta (1 / 2 + t * I) +
      exp (I * theta t) * (deriv riemannZeta (1 / 2 + t * I) * I)).re := by
    unfold hardyZ'; congr 1; ring
  rw [e]; exact h3

lemma deriv_hardyZ (t : ℝ) : deriv hardyZ t = hardyZ' t := (hasDerivAt_hardyZ t).deriv

/-! ## Gram points and the conjecture -/

/-- Classical fact (named hypothesis): θ′(t) > 0 for t ≥ 7, i.e. Re ψ(1/4 + it/2) > log π
(θ decreases to its minimum at t ≈ 6.29, then increases). -/
def ThetaDerivPos : Prop := ∀ t : ℝ, 7 ≤ t → 0 < thetaDeriv t

lemma strictMonoOn_theta (h : ThetaDerivPos) : StrictMonoOn theta (Ici 7) := by
  refine strictMonoOn_of_deriv_pos (convex_Ici 7) continuous_theta.continuousOn fun x hx => ?_
  rw [interior_Ici] at hx
  rw [(hasDerivAt_theta x).deriv]
  exact h x (le_of_lt hx)

/-- t is the Gram point g_n: θ(g_n) = πn (2406.06548 l. 164, 207), t ≥ 7. -/
def IsGram (n : ℤ) (t : ℝ) : Prop := 7 ≤ t ∧ theta t = n * π

/-- Good: Gram's law (−1)^n Z(g_n) > 0 holds (l. 212, l. 523). -/
def Good (n : ℤ) (t : ℝ) : Prop := 0 < (-1 : ℝ) ^ n * hardyZ t

/-- Bad: not good (l. 523, "Otherwise, it is considered bad"). -/
def Bad (n : ℤ) (t : ℝ) : Prop := ¬ Good n t

/-- conj:rep (l. 313–315): if g_n is a bad Gram point with good consecutive neighbours g_{n±1},
then |Z′(g_n)| > 4|Z(g_n)|. Neighbours: Gram points of index n ± 1 exist, and all of them are good.
The weakest reading (`jerbyConjEx_imp_jerbyConj`). -/
def JerbyConj : Prop :=
  ∀ (n : ℤ) (g : ℝ), IsGram n g → Bad n g →
    ((∃ g', IsGram (n - 1) g') ∧ ∀ g', IsGram (n - 1) g' → Good (n - 1) g') →
    ((∃ g', IsGram (n + 1) g') ∧ ∀ g', IsGram (n + 1) g' → Good (n + 1) g') →
    4 * |hardyZ g| < |deriv hardyZ g|

/-- conj:rep with the neighbours given: Gram points gl, gr of index n − 1, n + 1. -/
def JerbyConjEx : Prop :=
  ∀ (n : ℤ) (g gl gr : ℝ), IsGram n g → IsGram (n - 1) gl → IsGram (n + 1) gr →
    Bad n g → Good (n - 1) gl → Good (n + 1) gr → 4 * |hardyZ g| < |deriv hardyZ g|

lemma isGram_unique (h : ThetaDerivPos) {n : ℤ} {a b : ℝ} (ha : IsGram n a)
    (hb : IsGram n b) : a = b :=
  (strictMonoOn_theta h).injOn ha.1 hb.1 (ha.2.trans hb.2.symm)

/-- The given-neighbours reading implies `JerbyConj`, with no hypothesis: `JerbyConj` is weaker. -/
lemma jerbyConjEx_imp_jerbyConj : JerbyConjEx → JerbyConj := by
  rintro H n g hg hb ⟨⟨gl, hl⟩, hgl⟩ ⟨⟨gr, hr⟩, hgr⟩
  exact H n g gl gr hg hl hr hb (hgl gl hl) (hgr gr hr)

/-- With θ increasing (Gram points unique) the two readings agree. -/
lemma jerbyConj_imp_ex (h : ThetaDerivPos) : JerbyConj → JerbyConjEx := by
  intro H n g gl gr hg hl hr hb hgl hgr
  refine H n g hg hb ⟨⟨gl, hl⟩, fun g' h' => ?_⟩ ⟨⟨gr, hr⟩, fun g' h' => ?_⟩
  · rwa [isGram_unique h h' hl]
  · rwa [isGram_unique h h' hr]

/-- With θ increasing, Gram points of index n − 1, n, n + 1 are consecutive: g_{n−1} < g_n <
g_{n+1}, and every Gram point in (g_{n−1}, g_{n+1}) is g_n. -/
lemma consecutive (h : ThetaDerivPos) {n : ℤ} {a g b : ℝ} (ha : IsGram (n - 1) a)
    (hg : IsGram n g) (hb : IsGram (n + 1) b) :
    a < g ∧ g < b ∧ ∀ (m : ℤ) (t : ℝ), IsGram m t → a < t → t < b → m = n ∧ t = g := by
  have hm := strictMonoOn_theta h
  have hπ := Real.pi_pos
  refine ⟨?_, ?_, ?_⟩
  · refine (hm.lt_iff_lt ha.1 hg.1).1 ?_
    rw [ha.2, hg.2]; push_cast; linarith
  · refine (hm.lt_iff_lt hg.1 hb.1).1 ?_
    rw [hg.2, hb.2]; push_cast; linarith
  · intro m t ht hat htb
    have h1 : theta a < theta t := hm ha.1 ht.1 hat
    have h2 : theta t < theta b := hm ht.1 hb.1 htb
    rw [ha.2, ht.2] at h1
    rw [ht.2, hb.2] at h2
    push_cast at h1 h2
    have h1' : (n : ℝ) - 1 < m := by
      by_contra hc; have hc' := not_lt.1 hc; nlinarith
    have h2' : (m : ℝ) < n + 1 := by
      by_contra hc; have hc' := not_lt.1 hc; nlinarith
    have e1 : n - 1 < m := by exact_mod_cast h1'
    have e2 : m < n + 1 := by exact_mod_cast h2'
    have hmn : m = n := by omega
    subst hmn
    exact ⟨rfl, isGram_unique h ht hg⟩

lemma good_of_even_pos {n : ℤ} {t : ℝ} (hn : Even n) (h : 0 < hardyZ t) : Good n t := by
  unfold Good; rw [hn.neg_one_zpow, one_mul]; exact h

lemma good_of_odd_neg {n : ℤ} {t : ℝ} (hn : Odd n) (h : hardyZ t < 0) : Good n t := by
  unfold Good; rw [hn.neg_one_zpow]; linarith

lemma gramLaw_neg_of_odd_pos {n : ℤ} {t : ℝ} (hn : Odd n) (h : 0 < hardyZ t) :
    (-1 : ℝ) ^ n * hardyZ t < 0 := by
  rw [hn.neg_one_zpow]; linarith

lemma gramLaw_neg_of_even_neg {n : ℤ} {t : ℝ} (hn : Even n) (h : hardyZ t < 0) :
    (-1 : ℝ) ^ n * hardyZ t < 0 := by
  rw [hn.neg_one_zpow, one_mul]; exact h

lemma bad_of_gramLaw_neg {n : ℤ} {t : ℝ} (h : (-1 : ℝ) ^ n * hardyZ t < 0) : Bad n t :=
  fun hg => lt_asymm h hg

/-- A good–bad–good Gram triple with |Z′(g_n)| ≤ 4|Z(g_n)| refutes `JerbyConjEx`. -/
theorem not_jerbyConjEx_of {n : ℤ} {gl g gr : ℝ} (hl : IsGram (n - 1) gl) (hg : IsGram n g)
    (hr : IsGram (n + 1) gr) (hgl : Good (n - 1) gl) (hbad : Bad n g) (hgr : Good (n + 1) gr)
    (hμ : |deriv hardyZ g| ≤ 4 * |hardyZ g|) : ¬ JerbyConjEx :=
  fun H => absurd (H n g gl gr hg hl hr hbad hgl hgr) (not_lt.2 hμ)

end

end JerbyLean
