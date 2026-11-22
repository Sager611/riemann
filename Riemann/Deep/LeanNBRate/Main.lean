import Riemann.Deep.LeanNBRate.Final5

/-!
# DISCOVERIES #1 — the continuous Nyman–Beurling rate under RH (main theorem)

Assume RH. With D(λ) the L²(0,∞) distance from 1_{(0,1)} to span{ {θ/x} : λ ≤ θ ≤ 1 } and
L = log(1/λ): L·D(λ)² → Σ_ρ m_ρ²/|ρ|² as λ → 0⁺ (distinct nontrivial zeros, m_ρ = multiplicity),
the sum equals 2 + γ - log 4π when all zeros are simple, and Nyman's constrained distance has the
same limit. No hypothesis on ζ'(ρ) is used. The literature inputs are the Props bundled in NBInputs.
-/

open Filter Topology

namespace Riemann.Deep.LeanNBRate

/-- The literature inputs (each a Prop stating a published theorem, see LitInputs and LitInputsU1). -/
structure NBInputs : Prop where
  burnol : Burnol2002_Thm
  turan : Nazarov1993_Turan
  gg : GoldstonGonek2007_Thm1
  count : Titchmarsh_Thm9_2
  heights : Titchmarsh_Thm9_7
  davenport : Davenport_Ch12
  mellinFract : Titchmarsh_2_1_5
  plancherel : Rudin1987_Thm9_13

/-- **The continuous Nyman–Beurling rate under RH.** -/
theorem nb_rate (hI : NBInputs) (hRH : RiemannHypothesis) :
    Summable (fun ρ : NTZ => (mult (ρ : ℂ) : ℝ) ^ 2 / ‖(ρ : ℂ)‖ ^ 2) ∧
    Tendsto (fun lam : ℝ => Real.log (1 / lam) * nbDist lam ^ 2) (𝓝[>] 0) (𝓝 zeroSum) ∧
    ((∀ ρ : NTZ, mult (ρ : ℂ) = 1) →
      zeroSum = 2 + Real.eulerMascheroniConstant - Real.log (4 * Real.pi)) ∧
    Tendsto (fun lam : ℝ => Real.log (1 / lam) * nymanDist lam ^ 2) (𝓝[>] 0) (𝓝 zeroSum) :=
  mainClaim_of_upper hI.burnol hI.count hI.davenport
    (upper_bound hI.plancherel hI.mellinFract hI.heights hI.turan hI.gg hI.count) hRH

#print axioms nb_rate

end Riemann.Deep.LeanNBRate

