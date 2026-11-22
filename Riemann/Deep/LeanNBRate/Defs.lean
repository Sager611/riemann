import Mathlib.NumberTheory.LSeries.RiemannZeta
import Mathlib.NumberTheory.Harmonic.EulerMascheroni
import Mathlib.MeasureTheory.Function.LpSpace.Indicator
import Mathlib.Topology.MetricSpace.HausdorffDistance
import Mathlib.Analysis.Analytic.Order

/-!
# Continuous Nyman–Beurling distance: definitions and the main statement

Page: `discoveries/nyman-beurling-continuous-rate.md`, **Claim**:
Assume RH. Let D(λ) = dist_{L²(0,∞)}(1_{(0,1)}, span{ {θ/x} : λ ≤ θ ≤ 1 }) and L = log(1/λ). Then
L·D(λ)² → Σ_ρ m_ρ²/|ρ|² as λ → 0 (sum over the distinct nontrivial zeros, m_ρ the multiplicity;
= 2 + γ − log 4π = 0.0461… if all zeros are simple). Nyman's constrained variant has the same limit.
No hypothesis on ζ′(ρ) is used.
-/

open MeasureTheory Filter Topology Set

noncomputable section

namespace Riemann.Deep.LeanNBRate

/-- Lebesgue measure on `(0, ∞)`. -/
abbrev volPos : Measure ℝ := volume.restrict (Ioi (0 : ℝ))

/-- The real Hilbert space `L²(0, ∞)`. -/
abbrev L2pos : Type := Lp ℝ 2 volPos

/-- The function `x ↦ {θ / x}` (`Int.fract` is the fractional part). -/
def nbFun (θ : ℝ) (x : ℝ) : ℝ := Int.fract (θ / x)

theorem volPos_Ioo_ne_top : volPos (Ioo (0 : ℝ) 1) ≠ ⊤ := by
  refine ne_of_lt (lt_of_le_of_lt (Measure.restrict_apply_le (Ioi (0 : ℝ)) (Ioo (0 : ℝ) 1)) ?_)
  rw [Real.volume_Ioo]
  exact ENNReal.ofReal_lt_top

/-- The class of the indicator function `1_{(0,1)}` in `L²(0, ∞)`. -/
def chiVec : L2pos :=
  indicatorConstLp 2 measurableSet_Ioo volPos_Ioo_ne_top (1 : ℝ)

/-- The classes in `L²(0,∞)` of the functions `x ↦ {θ/x}`, `λ ≤ θ ≤ 1`. -/
def nbSet (lam : ℝ) : Set L2pos :=
  {f | ∃ θ ∈ Icc lam 1, (f : ℝ → ℝ) =ᵐ[volPos] nbFun θ}

/-- `D(λ)`: the distance in `L²(0,∞)` from `1_{(0,1)}` to `span { {θ/x} : λ ≤ θ ≤ 1 }`. -/
def nbDist (lam : ℝ) : ℝ :=
  Metric.infDist chiVec (Submodule.span ℝ (nbSet lam) : Set L2pos)

/-- Nyman's constrained variant: the classes of `Σ_k c_k {θ_k/x}` with `λ ≤ θ_k ≤ 1` and `Σ_k c_k θ_k = 0`. -/
def nymanSet (lam : ℝ) : Set L2pos :=
  {f | ∃ (n : ℕ) (c θ : Fin n → ℝ), (∀ k, θ k ∈ Icc lam 1) ∧ ∑ k, c k * θ k = 0 ∧
    (f : ℝ → ℝ) =ᵐ[volPos] fun x => ∑ k, c k * nbFun (θ k) x}

/-- The constrained distance `D_c(λ)` from `1_{(0,1)}` to Nyman's constrained set. -/
def nymanDist (lam : ℝ) : ℝ := Metric.infDist chiVec (nymanSet lam)

/-- The nontrivial zeros: `ζ(ρ) = 0`, `0 < Re ρ < 1`. -/
def NTZ : Set ℂ := {ρ | riemannZeta ρ = 0 ∧ 0 < ρ.re ∧ ρ.re < 1}

/-- The multiplicity `m_ρ`: Mathlib's analytic order of `riemannZeta` at `ρ`. -/
def mult (ρ : ℂ) : ℕ := analyticOrderNatAt riemannZeta ρ

/-- `Σ_ρ m_ρ² / |ρ|²` over the distinct nontrivial zeros. -/
def zeroSum : ℝ := ∑' ρ : NTZ, (mult (ρ : ℂ) : ℝ) ^ 2 / ‖(ρ : ℂ)‖ ^ 2

/-- The page's Claim, parts (a)–(c), with the summability of the zero sum. -/
def MainClaim : Prop :=
  RiemannHypothesis →
    Summable (fun ρ : NTZ => (mult (ρ : ℂ) : ℝ) ^ 2 / ‖(ρ : ℂ)‖ ^ 2) ∧
    Tendsto (fun lam : ℝ => Real.log (1 / lam) * nbDist lam ^ 2) (𝓝[>] 0) (𝓝 zeroSum) ∧
    ((∀ ρ : NTZ, mult (ρ : ℂ) = 1) →
      zeroSum = 2 + Real.eulerMascheroniConstant - Real.log (4 * Real.pi)) ∧
    Tendsto (fun lam : ℝ => Real.log (1 / lam) * nymanDist lam ^ 2) (𝓝[>] 0) (𝓝 zeroSum)

end Riemann.Deep.LeanNBRate
