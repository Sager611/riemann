import Riemann.Deep.LeanNBRate.Dual3
import Mathlib.Analysis.Calculus.Deriv.Star

/-!
# Dual problem: the Green form of N(·, p), interpolation at the low zeros, growth of E_low
-/

open MeasureTheory Filter Topology Set Complex ComplexConjugate
open scoped InnerProductSpace

noncomputable section

namespace Riemann.Deep.LeanNBRate

/-- The Green data of p: density −p″ + ¼p and endpoint masses. -/
def gDens (p p' p'' : ℝ → ℂ) (u : ℝ) : ℂ := -p'' u + (1 / 4 : ℂ) * p u
def gA (p p' : ℝ → ℂ) : ℂ := (1 / 2 : ℂ) * p 0 - p' 0
def gB (ℓ : ℝ) (p p' : ℝ → ℂ) : ℂ := (1 / 2 : ℂ) * p ℓ + p' ℓ

theorem hasDerivAt_conj {g g' : ℝ → ℂ} {x : ℝ} (h : HasDerivAt g (g' x) x) :
    HasDerivAt (fun y => conj (g y)) (conj (g' x)) x := by
  exact h.star

/-- **Green form**: N(f, p) = ∫ f·conj(−p″+¼p) + f(0)·conj(½p(0) − p′(0)) + f(ℓ)·conj(½p(ℓ) + p′(ℓ)). -/
theorem NForm_green {ℓ : ℝ} (f f' p p' p'' : ℝ → ℂ) (hf : ∀ x, HasDerivAt f (f' x) x)
    (hp' : ∀ x, HasDerivAt p' (p'' x) x) (hcf : Continuous f) (hcf' : Continuous f') (hcp : Continuous p)
    (hcp' : Continuous p') (hcp'' : Continuous p'') :
    NForm ℓ f f' p p' = (∫ u in (0 : ℝ)..ℓ, f u * conj (gDens p p' p'' u)) +
      f 0 * conj (gA p p') + f ℓ * conj (gB ℓ p p') := by
  have hibp := intervalIntegral.integral_mul_deriv_eq_deriv_mul (a := 0) (b := ℓ) (u := f) (u' := f')
    (v := fun y => conj (p' y)) (v' := fun y => conj (p'' y)) (fun x _ => hf x)
    (fun x _ => hasDerivAt_conj (hp' x)) (hcf'.intervalIntegrable _ _)
    ((Complex.continuous_conj.comp hcp'').intervalIntegrable _ _)
  unfold NForm oubil gDens gA gB
  have i1 : IntervalIntegrable (fun u => f' u * conj (p' u)) volume 0 ℓ :=
    (hcf'.mul (Complex.continuous_conj.comp hcp')).intervalIntegrable _ _
  have i2 : IntervalIntegrable (fun u => (1 / 4 : ℂ) * (f u * conj (p u))) volume 0 ℓ :=
    ((hcf.mul (Complex.continuous_conj.comp hcp)).intervalIntegrable _ _).const_mul _
  have i3 : IntervalIntegrable (fun u => f u * conj (p'' u)) volume 0 ℓ :=
    (hcf.mul (Complex.continuous_conj.comp hcp'')).intervalIntegrable _ _
  have i4 : IntervalIntegrable (fun u => f u * conj (p u)) volume 0 ℓ :=
    (hcf.mul (Complex.continuous_conj.comp hcp)).intervalIntegrable _ _
  rw [intervalIntegral.integral_add i1 i2]
  have e1 : ∫ u in (0 : ℝ)..ℓ, f' u * conj (p' u) =
      f ℓ * conj (p' ℓ) - f 0 * conj (p' 0) - ∫ u in (0 : ℝ)..ℓ, f u * conj (p'' u) := by
    rw [hibp]; ring
  have e2 : ∫ u in (0 : ℝ)..ℓ, f u * conj (-p'' u + (1 / 4 : ℂ) * p u) =
      -(∫ u in (0 : ℝ)..ℓ, f u * conj (p'' u)) + ∫ u in (0 : ℝ)..ℓ, (1 / 4 : ℂ) * (f u * conj (p u)) := by
    have e : ∀ u, f u * conj (-p'' u + (1 / 4 : ℂ) * p u) =
        -(f u * conj (p'' u)) + (1 / 4 : ℂ) * (f u * conj (p u)) := by
      intro u; simp only [map_add, map_neg, map_mul, map_div₀, map_one, map_ofNat]; ring
    simp_rw [e]
    rw [intervalIntegral.integral_add (f := fun u => -(f u * conj (p'' u)))
      (g := fun u => (1 / 4 : ℂ) * (f u * conj (p u))) i3.neg i2, intervalIntegral.integral_neg]
  rw [e1, e2]
  simp only [map_add, map_sub, map_mul, map_div₀, map_one, map_ofNat]
  ring

/-- The total mass of the Green data. -/
def gMass (ℓ : ℝ) (p p' p'' : ℝ → ℂ) : ℝ :=
  (∫ u in (0 : ℝ)..ℓ, ‖gDens p p' p'' u‖) + ‖gA p p'‖ + ‖gB ℓ p p'‖

theorem norm_NForm_le_sup {ℓ : ℝ} (hℓ : 0 ≤ ℓ) (f f' p p' p'' : ℝ → ℂ) (hf : ∀ x, HasDerivAt f (f' x) x)
    (hp' : ∀ x, HasDerivAt p' (p'' x) x) (hcf : Continuous f) (hcf' : Continuous f') (hcp : Continuous p)
    (hcp' : Continuous p') (hcp'' : Continuous p'') {K : ℝ} (hK : ∀ u ∈ Icc 0 ℓ, ‖f u‖ ≤ K) :
    ‖NForm ℓ f f' p p'‖ ≤ K * gMass ℓ p p' p'' := by
  rw [NForm_green f f' p p' p'' hf hp' hcf hcf' hcp hcp' hcp'']
  have hK0 : 0 ≤ K := (norm_nonneg _).trans (hK 0 ⟨le_rfl, hℓ⟩)
  have hgc : Continuous (gDens p p' p'') := by unfold gDens; fun_prop
  have h1 : ‖∫ u in (0 : ℝ)..ℓ, f u * conj (gDens p p' p'' u)‖ ≤ K * ∫ u in (0 : ℝ)..ℓ, ‖gDens p p' p'' u‖ := by
    refine (intervalIntegral.norm_integral_le_integral_norm hℓ).trans ?_
    rw [← intervalIntegral.integral_const_mul]
    refine intervalIntegral.integral_mono_on hℓ ?_ ?_ ?_
    · exact ((hcf.mul (Complex.continuous_conj.comp hgc)).norm).intervalIntegrable _ _
    · exact (continuous_const.mul hgc.norm).intervalIntegrable _ _
    · intro u hu
      rw [norm_mul, Complex.norm_conj]
      exact mul_le_mul_of_nonneg_right (hK u hu) (norm_nonneg _)
  have h2 : ‖f 0 * conj (gA p p')‖ ≤ K * ‖gA p p'‖ := by
    rw [norm_mul, Complex.norm_conj]; exact mul_le_mul_of_nonneg_right (hK 0 ⟨le_rfl, hℓ⟩) (norm_nonneg _)
  have h3 : ‖f ℓ * conj (gB ℓ p p')‖ ≤ K * ‖gB ℓ p p'‖ := by
    rw [norm_mul, Complex.norm_conj]; exact mul_le_mul_of_nonneg_right (hK ℓ ⟨hℓ, le_rfl⟩) (norm_nonneg _)
  calc _ ≤ ‖∫ u in (0 : ℝ)..ℓ, f u * conj (gDens p p' p'' u)‖ + ‖f 0 * conj (gA p p')‖ +
        ‖f ℓ * conj (gB ℓ p p')‖ := norm_add₃_le
    _ ≤ K * (∫ u in (0 : ℝ)..ℓ, ‖gDens p p' p'' u‖) + K * ‖gA p p'‖ + K * ‖gB ℓ p p'‖ := by gcongr
    _ = K * gMass ℓ p p' p'' := by unfold gMass; ring

end Riemann.Deep.LeanNBRate

