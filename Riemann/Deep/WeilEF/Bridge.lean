/-
Riemann/Deep/WeilEF/Bridge.lean — our Weil functional (Riemann/Deep/Weil/Defs.lean) equals the
literature right-hand side of the explicit formula proved in the ported zeta23 development
(adapted from anthropics/formal-math zeta23, Apache-2.0), hence equals the zero sum.

* `weilQ_eq_literatureRHS`: `weilQ f = literatureRHS (autocorr f)` for every test function
  (the `− g(0) log π` in our archimedean term is `(1/2π)∫ h · log π` by Fourier inversion).
* `weilQ_eq_zeroSum`: `weilQ f = Σ_ρ m_ρ F(γ_ρ) · conj F(conj γ_ρ)`, `F = paperFT f`,
  `γ_ρ = (ρ − 1/2)/i`, the sum absolutely convergent (Mathlib's ζ zeros, `analyticOrderAt`).
-/
import Riemann.Deep.Weil.Defs
import Riemann.Deep.WeilEF.Z23.WeilEF.Main
import Riemann.Deep.WeilEF.Z23.GammaFacts.Complete

noncomputable section

open Complex MeasureTheory Set
open scoped ComplexConjugate

namespace Riemann.Deep.WeilEF

open Riemann.Deep.Weil Zeta23 Zeta23.EF

lemma ft_eq_paperFT (g : ℝ → ℂ) (z : ℂ) : ft g z = paperFT g z := rfl

/-! ### The (fixed) archimedean term

`archTerm` in `Riemann/Deep/Weil/Defs.lean` used to read `(1/2π) * ∫ r, h r * Re ψ(..) - g 0 * log π`,
whose binder swallowed `- g 0 * log π` (see `Bug.lean`, which keeps a frozen copy `archTermOld`).
It is now parenthesised and agrees with `archTermFix` definitionally. -/

/-- The intended archimedean term `(1/2π) ∫ h(r) Re ψ(1/4 + ir/2) dr − g(0) log π`. -/
def archTermFix (g : ℝ → ℂ) : ℂ :=
  (1 / (2 * Real.pi) : ℂ) *
      (∫ r : ℝ, ft g r * ((logDeriv Complex.Gamma (1 / 4 + I * r / 2)).re : ℂ)) -
    g 0 * (Real.log Real.pi : ℂ)

/-- The intended Weil functional. -/
def weilQFix (f : ℝ → ℂ) : ℂ :=
  poleTerm (autocorr f) + archTermFix (autocorr f) - primeTerm (autocorr f)

/-- The intended window positivity. -/
def WeilPosFix (L : ℝ) : Prop :=
  ∀ f : ℝ → ℂ, TestFn (L / 2) f → 0 ≤ (weilQFix f).re

/-- After the parenthesisation fix, `archTerm` is `archTermFix` by `rfl`. -/
lemma archTerm_eq_fix (g : ℝ → ℂ) : archTerm g = archTermFix g := rfl

lemma weilQ_eq_fix (f : ℝ → ℂ) : weilQ f = weilQFix f := rfl

/-- `WeilPos` (Defs.lean) is the honest Weil positivity `WeilPosFix`. -/
lemma weilPos_iff_fix (L : ℝ) : WeilPos L ↔ WeilPosFix L := Iff.rfl

lemma autocorr_eq_weilTest (f : ℝ → ℂ) : autocorr f = weilTest f f := by
  funext x
  simp only [autocorr, weilTest, convolution_def, ContinuousLinearMap.mul_apply', tilde, neg_sub]

/-- Regularity facts of a test function. -/
lemma _root_.Riemann.Deep.Weil.TestFn.contDiff_two {a : ℝ} {f : ℝ → ℂ} (hf : TestFn a f) : ContDiff ℝ 2 f :=
  hf.1.of_le (WithTop.coe_le_coe.mpr (le_top : (2 : ℕ∞) ≤ ⊤))

lemma _root_.Riemann.Deep.Weil.TestFn.continuous {a : ℝ} {f : ℝ → ℂ} (hf : TestFn a f) : Continuous f :=
  hf.contDiff_two.continuous

lemma _root_.Riemann.Deep.Weil.TestFn.hasCompactSupport {a : ℝ} {f : ℝ → ℂ} (hf : TestFn a f) : HasCompactSupport f :=
  hasCompactSupport_of_tsupport_subset (L := a)
    (closure_minimal (hf.2.trans (Icc_subset_Icc (le_refl _) (le_refl _))) isClosed_Icc)

lemma _root_.Riemann.Deep.Weil.TestFn.autocorr_contDiff {a : ℝ} {f : ℝ → ℂ} (hf : TestFn a f) :
    ContDiff ℝ 2 (autocorr f) := by
  rw [autocorr_eq_weilTest]
  exact weilTest_contDiff hf.contDiff_two hf.continuous hf.hasCompactSupport

lemma _root_.Riemann.Deep.Weil.TestFn.autocorr_hasCompactSupport {a : ℝ} {f : ℝ → ℂ} (hf : TestFn a f) :
    HasCompactSupport (autocorr f) := by
  rw [autocorr_eq_weilTest]
  exact weilTest_hasCompactSupport hf.hasCompactSupport hf.hasCompactSupport

/-- The archimedean term: `(1/2π)∫ h·[Re ψ − log π] = (1/2π)∫ h·Re ψ − k(0) log π`. -/
lemma arch_eq {k : ℝ → ℂ} (hk : ContDiff ℝ 2 k) (hkc : HasCompactSupport k) :
    (1 / (2 * Real.pi) : ℂ) * ∫ r : ℝ, paperFT k r * (gammaBracket r : ℂ)
      = (1 / (2 * Real.pi) : ℂ) * (∫ r : ℝ, ft k r *
          ((logDeriv Complex.Gamma (1 / 4 + I * r / 2)).re : ℂ)) - k 0 * (Real.log Real.pi : ℂ) := by
  have hFk := integrable_fourier_of_contDiff_two hk hkc
  have hI : Integrable (fun r : ℝ => paperFT k r) := integrable_paperFT_ofReal hFk
  have hmu := integrable_paperFT_mul_mu hk hkc gammaFacts
  have hIg : Integrable (fun r : ℝ => paperFT k r * (gammaBracket r : ℂ)) := by
    have := hmu.const_mul ((2 * Real.pi : ℝ) : ℂ)
    refine this.congr (Filter.Eventually.of_forall fun r => ?_)
    simp only [mu, gammaBracket]
    push_cast
    field_simp
  have hIl : Integrable (fun r : ℝ => (Real.log Real.pi : ℂ) * paperFT k r) := hI.const_mul _
  have hpt : (fun r : ℝ => ft k r * ((logDeriv Complex.Gamma (1 / 4 + I * r / 2)).re : ℂ))
      = fun r : ℝ => paperFT k r * (gammaBracket r : ℂ) + (Real.log Real.pi : ℂ) * paperFT k r := by
    funext r
    simp only [ft_eq_paperFT, gammaBracket, Complex.digamma_def]
    push_cast
    ring
  have hinv := paper_inversion hk.continuous (hk.continuous.integrable_of_hasCompactSupport hkc) hFk 0
  simp only [ofReal_zero, mul_zero, Complex.exp_zero, mul_one] at hinv
  rw [hpt, integral_add hIg hIl, integral_const_mul, hinv]
  ring

/-- **Bridge.** Our (intended) `weilQ` is the literature right-hand side at `k = autocorr f`. -/
theorem weilQFix_eq_literatureRHS {a : ℝ} {f : ℝ → ℂ} (hf : TestFn a f) :
    weilQFix f = literatureRHS (autocorr f) := by
  have hk := hf.autocorr_contDiff
  have hkc := hf.autocorr_hasCompactSupport
  unfold weilQFix poleTerm archTermFix primeTerm literatureRHS
  rw [arch_eq hk hkc, ft_eq_paperFT, ft_eq_paperFT, neg_div]
  ring

/-- `paperFT (autocorr f) z = F z · conj (F (conj z))`. -/
lemma paperFT_autocorr {a : ℝ} {f : ℝ → ℂ} (hf : TestFn a f) (z : ℂ) :
    paperFT (autocorr f) z = paperFT f z * conj (paperFT f (conj z)) := by
  rw [autocorr_eq_weilTest]
  exact paperFT_weilTest hf.continuous hf.continuous hf.hasCompactSupport hf.hasCompactSupport z

/-- **Explicit formula for `weilQ`** (sorry-free; Mathlib's ζ zeros with `analyticOrderAt`
multiplicities): the zero sum is summable and equals `weilQ f`. -/
theorem weilQFix_eq_zeroSum {a : ℝ} {f : ℝ → ℂ} (hf : TestFn a f) :
    Summable (fun ρ : zetaZeroConfig.carrier => (zetaZeroConfig.mult ρ : ℂ) *
        (paperFT f (gammaOf ρ) * conj (paperFT f (conj (gammaOf ρ))))) ∧
    weilQFix f = ∑' ρ : zetaZeroConfig.carrier, (zetaZeroConfig.mult ρ : ℂ) *
        (paperFT f (gammaOf ρ) * conj (paperFT f (conj (gammaOf ρ)))) := by
  obtain ⟨hs, he⟩ := WeilEF.EF_lit_zetaZeroConfig (autocorr f) hf.autocorr_contDiff
    hf.autocorr_hasCompactSupport
  simp only [paperFT_autocorr hf] at hs he
  exact ⟨hs, (weilQFix_eq_literatureRHS hf).trans he.symm⟩

end Riemann.Deep.WeilEF
