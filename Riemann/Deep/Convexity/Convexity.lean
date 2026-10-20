import Riemann.Deep.Convexity.MaxMod
import Riemann.Fayed.Paper

/-!
# Discharging Fayed.Inputs.convexity from Mathlib alone

convexity_holds is literally the field Inputs.convexity (checked by building an Inputs below).
Proof: partial summation + identity theorem (Partial), functional equation with exact
|Γ(1-it)|² (FE), maximum modulus on rectangles of ζ(s)(-is)^{(s-1-δ)/2}e^{εs²} with a
bootstrap on the left edge (Bootstrap, MaxMod), δ = min(1, 2ε), ε' = 1/t².
-/

open Complex
open scoped Real ComplexConjugate

namespace Fayed.Convexity

theorem convexity_holds : ∀ η : ℝ, 0 < η → η < 1 / 2 → ∀ ε : ℝ, 0 < ε → ∃ C : ℝ, ∀ T : ℝ,
    1 ≤ |T| → ‖riemannZeta (((1 - η : ℝ) : ℂ) + T * I)‖ ≤ C * |T| ^ (η / 2 + ε) := by
  intro η hη0 hη1 ε hε
  refine ⟨10 * Real.exp 1 * Kc (min 1 (2 * ε)), fun T hT => ?_⟩
  rcases le_or_gt 0 T with h | h
  · rw [abs_of_nonneg h] at hT ⊢
    exact convexity_pos hη0 hη1 hε hT
  · rw [abs_of_neg h] at hT ⊢
    have e : ((1 - η : ℝ) : ℂ) + T * I = conj (((1 - η : ℝ) : ℂ) + ((-T : ℝ) : ℂ) * I) := by
      apply Complex.ext <;> simp
    rw [e, norm_zeta_conj]
    exact convexity_pos hη0 hη1 hε hT

/-- Fayed's literature inputs minus the convexity bound (now a theorem). -/
structure InputsNoConvexity : Prop where
  em17 : ∀ s : ℂ, 0 < s.re → s.re < 1 → ∃ C : ℝ, ∀ u : ℝ, 1 ≤ u →
    ‖S u s - (riemannZeta s + (u : ℂ) ^ (1 - s) / (1 - s) - EMcorr u s)‖ ≤ C * u ^ (-(s.re + 3))
  perron33 : ∀ s : ℂ, 0 < s.re → s.re < 1 → ∀ η : ℝ, 0 < η → η < s.re → η < 1 - s.re →
    ∀ γ : ℝ, 0 < γ → γ < 1 / 4 → ∀ ε : ℝ, 0 < ε → ∃ C U : ℝ, ∀ M₀ : ℕ, ∀ u ∈ Dγ γ M₀, U ≤ u →
      ‖S u s - (J η u s + (u : ℂ) ^ (1 - s) / (1 - s))‖ ≤ C * u ^ (-4 + ε)
  leibniz44 : ∀ η σ₀ t₀ : ℝ, 0 < η → η < σ₀ → η < 1 - σ₀ → ∀ u : ℝ, 0 < u →
    HasDerivAt (Ffun η σ₀ t₀) (Fder η σ₀ t₀ u) u
  strip : ∀ s : ℂ, riemannZeta s = 0 → (¬∃ n : ℕ, s = -2 * (n + 1)) → s ≠ 1 →
    0 < s.re ∧ s.re < 1
  reflect : ∀ s : ℂ, riemannZeta s = 0 → 0 < s.re → s.re < 1 → riemannZeta (1 - conj s) = 0

/-- Rebuild the full Inputs: the convexity field is supplied by convexity_holds. -/
theorem inputs_of (h : InputsNoConvexity) : Inputs :=
  ⟨h.em17, h.perron33, h.leibniz44, convexity_holds, h.strip, h.reflect⟩

/-- Fayed's argument without the convexity hypothesis. -/
theorem fayed_RH_noConvexity (h : InputsNoConvexity) (thm1 : Theorem1C) : RiemannHypothesis :=
  fayed_RH (inputs_of h) thm1

end Fayed.Convexity

#print axioms Fayed.Convexity.convexity_holds
#print axioms Fayed.Convexity.fayed_RH_noConvexity
