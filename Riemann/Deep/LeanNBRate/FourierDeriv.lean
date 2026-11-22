import Riemann.Deep.LeanNBRate.Kernel
import Mathlib.Analysis.Calculus.IteratedDeriv.Lemmas

/-!
# Derivatives of `θ ↦ ∫ f(u) e^{iθu} du` for test functions `f`

`fint f θ = ∫ f(u) e^{iθu} du = 2π · kern f (−θ)`, and
`iteratedDeriv n (fint f) θ = fint ((iu)^n f) θ`.
-/

open MeasureTheory Filter Topology Set Complex FourierTransform

noncomputable section

namespace Riemann.Deep.LeanNBRate

/-- `fint f θ = ∫ f(u) e^{iθu} du`. -/
def fint (f : ℝ → ℂ) (θ : ℝ) : ℂ := ∫ u : ℝ, f u * cexp (I * θ * u)

theorem fint_eq_kern (f : ℝ → ℂ) (θ : ℝ) : fint f θ = (2 * Real.pi : ℂ) * kern f (-θ) := by
  rw [kern_eq]
  have hπ : (2 * Real.pi : ℂ) ≠ 0 := by exact_mod_cast (show (0:ℝ) < 2 * Real.pi by positivity).ne'
  rw [← mul_assoc, mul_inv_cancel₀ hπ, one_mul]
  unfold fint
  refine integral_congr_ae (Eventually.of_forall fun u => ?_)
  simp only; push_cast; ring_nf

/-- The test function `u ↦ (iu)^n f(u)`. -/
def pmul (n : ℕ) (f : ℝ → ℂ) (u : ℝ) : ℂ := (I * u) ^ n * f u

theorem IsTest.pmul_test {f : ℝ → ℂ} (hf : IsTest f) (n : ℕ) : IsTest (pmul n f) := by
  refine ⟨?_, ?_⟩
  · unfold Riemann.Deep.LeanNBRate.pmul
    exact ((contDiff_const.mul ofRealCLM.contDiff).pow n).mul hf.smooth
  · unfold Riemann.Deep.LeanNBRate.pmul
    exact hf.supp.mul_left

/-- `iteratedDeriv n (fint f) = fint ((iu)^n f)`. -/
theorem IsTest.iteratedDeriv_fint {f : ℝ → ℂ} (hf : IsTest f) (n : ℕ) :
    iteratedDeriv n (fint f) = fint (pmul n f) := by
  have hπ : (0 : ℝ) < 2 * Real.pi := by positivity
  have hπc : (2 * Real.pi : ℂ) ≠ 0 := by exact_mod_cast hπ.ne'
  -- fint f θ = 2π · (2π)⁻¹ 𝓕 f (−θ/(2π)) = 𝓕 f ((−(2π)⁻¹) θ)
  have hrep : ∀ g : ℝ → ℂ, fint g = fun θ => 𝓕 g ((-(2 * Real.pi)⁻¹) * θ) := by
    intro g; ext θ
    rw [fint_eq_kern]; unfold kern
    rw [← mul_assoc, mul_inv_cancel₀ hπc, one_mul]
    congr 1; field_simp
  have hint : ∀ k : ℕ, k ≤ (⊤ : ℕ∞) → Integrable (fun x : ℝ => x ^ k • f x) := by
    intro k _
    have := (hf.pmul_test k).integrable
    have h2 : Integrable (fun x : ℝ => (I ^ k)⁻¹ * pmul k f x) := this.const_mul _
    refine h2.congr (Eventually.of_forall fun x => ?_)
    simp only [Riemann.Deep.LeanNBRate.pmul]
    rw [Complex.real_smul, mul_pow, ← mul_assoc, ← mul_assoc,
      inv_mul_cancel₀ (pow_ne_zero _ I_ne_zero), one_mul]
    push_cast; ring
  have hF := Real.iteratedDeriv_fourier (f := f) (N := (⊤ : ℕ∞)) (n := n) hint le_top
  have hsmooth : ContDiff ℝ n (𝓕 f) := by
    have : 𝓕 f = ((𝓕 hf.toS : SchwartzMap ℝ ℂ) : ℝ → ℂ) := by rw [SchwartzMap.fourier_coe]; rfl
    rw [this]; exact (𝓕 hf.toS).smooth n
  rw [hrep f, hrep (pmul n f)]
  ext θ
  have hcomp := congrFun (iteratedDeriv_comp_const_smul (𝕜 := ℝ) (f := 𝓕 f) hsmooth
    (-(2 * Real.pi)⁻¹)) θ
  rw [hcomp, hF]
  rw [Real.fourier_real_eq, Real.fourier_real_eq, ← integral_smul]
  refine integral_congr_ae (Eventually.of_forall fun x => ?_)
  simp only [smul_eq_mul, Riemann.Deep.LeanNBRate.pmul, Complex.real_smul, Circle.smul_def]
  have : ((-(2 * Real.pi)⁻¹ : ℝ) : ℂ) ^ n * (-2 * Real.pi * I * x) ^ n = (I * x) ^ n := by
    rw [← mul_pow]; congr 1; push_cast; field_simp
  rw [← this]; push_cast; ring

end Riemann.Deep.LeanNBRate
