import Riemann.Deep.Inputs.Leibniz
import Riemann.Deep.Inputs.EM
import Riemann.Deep.Convexity.Convexity
import Riemann.Fayed.StripReflect

/-!
# Fayed.Inputs: five of six fields are theorems

em17 (EM), leibniz44 (Leibniz), convexity (Convexity), strip, reflect (StripReflect) are proved.
perron33 (truncated Perron + contour shift) is the one remaining named hypothesis.
-/

open Complex
open scoped Real ComplexConjugate

namespace Fayed.InputsProof

/-- The perron33 field of Inputs, as a named hypothesis. -/
def Perron33 : Prop :=
  ∀ s : ℂ, 0 < s.re → s.re < 1 → ∀ η : ℝ, 0 < η → η < s.re → η < 1 - s.re →
    ∀ γ : ℝ, 0 < γ → γ < 1 / 4 → ∀ ε : ℝ, 0 < ε → ∃ C U : ℝ, ∀ M₀ : ℕ, ∀ u ∈ Dγ γ M₀, U ≤ u →
      ‖S u s - (J η u s + (u : ℂ) ^ (1 - s) / (1 - s))‖ ≤ C * u ^ (-4 + ε)

theorem inputs_of_perron33 (h : Perron33) : Inputs :=
  ⟨em17_holds, h, leibniz44_holds, Fayed.Convexity.convexity_holds, strip_holds, reflect_holds⟩

/-- Inputs is equivalent to its perron33 field alone. -/
theorem inputs_iff_perron33 : Inputs ↔ Perron33 :=
  ⟨fun h => h.perron33, inputs_of_perron33⟩

/-- Fayed's reduction, conditional only on the truncated-Perron step. -/
theorem RH_of_theorem1C_perron33 (h : Perron33) : Theorem1C → RiemannHypothesis :=
  fayed_RH (inputs_of_perron33 h)

end Fayed.InputsProof

#print axioms Fayed.InputsProof.em17_holds
#print axioms Fayed.InputsProof.leibniz44_holds
#print axioms Fayed.InputsProof.RH_of_theorem1C_perron33
