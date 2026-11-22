import Riemann.Deep.LeanNBRate.ExpSum
import Mathlib.Analysis.Fourier.Inversion
import Mathlib.Analysis.Fourier.FourierTransformDeriv
import Mathlib.Analysis.Distribution.SchwartzSpace.Fourier

/-!
# Fourier kernels of test functions

For a smooth compactly supported `G : ℝ → ℂ`, `kern G v = (2π)⁻¹ ∫ G(θ) e^{−iθv} dθ`; then
`∫ kern G v · (iv)^l e^{iθv} dv = G^{(l)}(θ)` (Fourier inversion) and `|v|^N |kern G v| ≤ (2π)⁻¹ ∫ |G^{(N)}|`.
-/

open MeasureTheory Filter Topology Set Complex FourierTransform

noncomputable section

namespace Riemann.Deep.LeanNBRate

/-- Smooth compactly supported complex functions on `ℝ`. -/
structure IsTest (G : ℝ → ℂ) : Prop where
  smooth : ContDiff ℝ (⊤ : ℕ∞) G
  supp : HasCompactSupport G

namespace IsTest

variable {G : ℝ → ℂ} (hG : IsTest G)
include hG

theorem continuous : Continuous G := hG.smooth.continuous

theorem integrable : Integrable G := hG.continuous.integrable_of_hasCompactSupport hG.supp

theorem iterDeriv (n : ℕ) : IsTest (_root_.iteratedDeriv n G) := by
  rw [iteratedDeriv_eq_iterate]
  refine ⟨hG.smooth.iterate_deriv n, ?_⟩
  induction n with
  | zero => simpa using hG.supp
  | succ n ih => rw [Function.iterate_succ_apply']; exact ih.deriv

/-- The Schwartz function attached to a test function. -/
def toS : SchwartzMap ℝ ℂ := hG.supp.toSchwartzMap (by simpa using hG.smooth)

theorem toS_apply (x : ℝ) : hG.toS x = G x := rfl

theorem integrable_fourier : Integrable (𝓕 G) := by
  have : 𝓕 G = ((𝓕 hG.toS : SchwartzMap ℝ ℂ) : ℝ → ℂ) := by
    rw [SchwartzMap.fourier_coe]; rfl
  rw [this]; exact (𝓕 hG.toS).integrable

theorem fourierInv_fourier : 𝓕⁻ (𝓕 G) = G :=
  hG.continuous.fourierInv_fourier_eq hG.integrable hG.integrable_fourier

end IsTest

/-- `kern G v = (2π)⁻¹ ∫ G(θ) e^{−iθv} dθ`. -/
def kern (G : ℝ → ℂ) (v : ℝ) : ℂ := (2 * Real.pi : ℂ)⁻¹ * 𝓕 G (v / (2 * Real.pi))

theorem kern_eq (G : ℝ → ℂ) (v : ℝ) :
    kern G v = (2 * Real.pi : ℂ)⁻¹ * ∫ θ : ℝ, cexp (-(I * θ * v)) * G θ := by
  unfold kern
  rw [Real.fourier_real_eq_integral_exp_smul]
  congr 1
  refine integral_congr_ae (Eventually.of_forall fun θ => ?_)
  simp only [smul_eq_mul]
  congr 1
  congr 1
  have hπ : (Real.pi : ℂ) ≠ 0 := by exact_mod_cast Real.pi_ne_zero
  push_cast
  field_simp

theorem norm_kern_le (G : ℝ → ℂ) (v : ℝ) :
    ‖kern G v‖ ≤ (2 * Real.pi)⁻¹ * ∫ θ : ℝ, ‖G θ‖ := by
  rw [kern_eq, norm_mul]
  have h2 : ‖(2 * Real.pi : ℂ)⁻¹‖ = (2 * Real.pi)⁻¹ := by
    rw [norm_inv]; congr 1
    rw [show (2 * Real.pi : ℂ) = ((2 * Real.pi : ℝ) : ℂ) by push_cast; ring, Complex.norm_real,
      Real.norm_eq_abs, abs_of_pos (by positivity)]
  rw [h2]
  gcongr
  refine (norm_integral_le_integral_norm _).trans (le_of_eq ?_)
  congr 1; ext θ
  rw [norm_mul, show (-(I * θ * v) : ℂ) = ((-(θ * v) : ℝ) : ℂ) * I by push_cast; ring,
    Complex.norm_exp_ofReal_mul_I, one_mul]

namespace IsTest

variable {G : ℝ → ℂ} (hG : IsTest G)
include hG

/-- Fourier inversion: `∫ kern G v · e^{iθv} dv = G θ`. -/
theorem integral_kern_mul_exp (θ : ℝ) : ∫ v : ℝ, kern G v * cexp (I * θ * v) = G θ := by
  have hπ : (0 : ℝ) < 2 * Real.pi := by positivity
  set g : ℝ → ℂ := fun ξ => (2 * Real.pi : ℂ)⁻¹ * (𝓕 G ξ * cexp (I * θ * (2 * Real.pi * ξ)))
  have hfun : (fun v : ℝ => kern G v * cexp (I * θ * v)) = fun v => g (v / (2 * Real.pi)) := by
    ext v; simp only [kern, g]
    have : (2 * Real.pi * (v / (2 * Real.pi)) : ℂ) = v := by
      have : (2 * Real.pi : ℂ) ≠ 0 := by exact_mod_cast hπ.ne'
      push_cast; field_simp
    rw [show ((v / (2 * Real.pi) : ℝ) : ℂ) = (v : ℂ) / (2 * Real.pi) by push_cast; ring]
    rw [this]; ring
  rw [hfun, Measure.integral_comp_div, abs_of_pos hπ]
  simp only [g, integral_const_mul]
  have hinv : ∫ ξ : ℝ, 𝓕 G ξ * cexp (I * θ * (2 * Real.pi * ξ)) = 𝓕⁻ (𝓕 G) θ := by
    rw [Real.fourierInv_eq']
    refine integral_congr_ae (Eventually.of_forall fun ξ => ?_)
    have hin : (inner ℝ ξ θ : ℝ) = ξ * θ := by
      simp [RCLike.inner_apply, mul_comm]
    simp only [smul_eq_mul, hin]
    rw [mul_comm]
    congr 1; congr 1
    push_cast; ring
  rw [hinv, hG.fourierInv_fourier]
  rw [Complex.real_smul]
  have : (2 * Real.pi : ℂ) ≠ 0 := by exact_mod_cast hπ.ne'
  push_cast
  field_simp

/-- `kern (G^{(n)}) v = (iv)^n · kern G v`. -/
theorem kern_iteratedDeriv (n : ℕ) (v : ℝ) :
    kern (_root_.iteratedDeriv n G) v = (I * v) ^ n * kern G v := by
  unfold kern
  have h := Real.fourier_iteratedDeriv (f := G) (N := (⊤ : ℕ∞)) (n := n) hG.smooth
    (fun k _ => (hG.iterDeriv k).integrable) le_top
  rw [h]
  simp only [smul_eq_mul]
  have hπ : (2 * Real.pi : ℂ) ≠ 0 := by exact_mod_cast (show (0:ℝ) < 2 * Real.pi by positivity).ne'
  rw [show ((v / (2 * Real.pi) : ℝ) : ℂ) = (v : ℂ) / (2 * Real.pi) by push_cast; ring]
  have : (2 * Real.pi * I * ((v : ℂ) / (2 * Real.pi))) = I * v := by field_simp
  rw [this]; ring

/-- `∫ kern G v · (iv)^l e^{iθv} dv = G^{(l)}(θ)`. -/
theorem integral_kern_pow_mul_exp (l : ℕ) (θ : ℝ) :
    ∫ v : ℝ, kern G v * (I * v) ^ l * cexp (I * θ * v) = _root_.iteratedDeriv l G θ := by
  have := (hG.iterDeriv l).integral_kern_mul_exp θ
  rw [← this]
  refine integral_congr_ae (Eventually.of_forall fun v => ?_)
  simp only
  rw [hG.kern_iteratedDeriv l v]; ring

/-- `|v|^N |kern G v| ≤ (2π)⁻¹ ∫ |G^{(N)}|`. -/
theorem norm_pow_mul_kern_le (N : ℕ) (v : ℝ) :
    |v| ^ N * ‖kern G v‖ ≤ (2 * Real.pi)⁻¹ * ∫ θ : ℝ, ‖_root_.iteratedDeriv N G θ‖ := by
  have h := norm_kern_le (_root_.iteratedDeriv N G) v
  rw [hG.kern_iteratedDeriv N v, norm_mul, norm_pow, norm_mul, Complex.norm_I, one_mul,
    Complex.norm_real, Real.norm_eq_abs] at h
  exact h

end IsTest

end Riemann.Deep.LeanNBRate
