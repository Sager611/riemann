module

public import Mathlib.NumberTheory.LSeries.RiemannZeta
public import Mathlib.NumberTheory.LSeries.ZetaZeros
public import Mathlib.Analysis.SpecialFunctions.Gamma.Digamma
public import Mathlib.NumberTheory.Harmonic.EulerMascheroni
public import Mathlib.MeasureTheory.Integral.IntervalIntegral.Basic

/-!
# LeanLi: the objects of the Claim of `discoveries/li-remainder.md`

All objects are built from Mathlib's `riemannZeta`, `analyticOrderNatAt`, `Complex.digamma`,
`Real.eulerMascheroniConstant`; the theorem statements use Mathlib's `RiemannHypothesis`.
-/

@[expose] public section

open Complex Filter Topology

namespace LeanLi

/-- The nontrivial zeros of `ζ`, in the sense of Mathlib's `RiemannHypothesis`:
`ζ s = 0`, `s` is not a trivial zero `-2(n+1)`, and `s ≠ 1`. -/
def NontrivialZeros : Set ℂ :=
  {s | riemannZeta s = 0 ∧ (¬∃ n : ℕ, s = -2 * (n + 1)) ∧ s ≠ 1}

/-- The multiplicity of `s` as a zero of `ζ`: Mathlib's order of vanishing. -/
noncomputable def mult (s : ℂ) : ℕ := analyticOrderNatAt riemannZeta s

/-- Partial Li sum over the nontrivial zeros `ρ` with `|ρ| ≤ T`, counted with multiplicity:
`Σ_{|ρ| ≤ T} m_ρ [1 - (1 - 1/ρ)ⁿ]`. Lagarias (arXiv math/0404394, TeX l.212–218): the sum runs
over the nontrivial zeros counted with multiplicity and "the (conditionally convergent) sum is to be
interpreted as `lim_{T → ∞} Σ_{ρ : |ρ| ≤ T}`". The index set is finite. -/
noncomputable def liPartial (n : ℕ) (T : ℝ) : ℂ :=
  ∑ᶠ ρ ∈ {ρ : ℂ | ρ ∈ NontrivialZeros ∧ ‖ρ‖ ≤ T}, (mult ρ : ℂ) * (1 - (1 - 1 / ρ) ^ n)

/-- The Li coefficient `λ_n = Σ_ρ [1 - (1 - 1/ρ)ⁿ] := lim_{T → ∞} liPartial n T`
(the limit exists and is real: see `Zeros.lean`). -/
noncomputable def liCoeff (n : ℕ) : ℂ := limUnder atTop (liPartial n)

/-- `λ_n` as a real number (`liCoeff n` is real). -/
noncomputable def lam (n : ℕ) : ℝ := (liCoeff n).re

/-- The trend `(n/2)(log n - 1 - log 2π + γ)`. -/
noncomputable def liTrend (n : ℕ) : ℝ :=
  (n / 2 : ℝ) * (Real.log n - 1 - Real.log (2 * Real.pi) + Real.eulerMascheroniConstant)

/-- The remainder `R(n) = λ_n - (n/2)(log n - 1 - log 2π + γ)`. -/
noncomputable def R (n : ℕ) : ℝ := lam n - liTrend n

/-- `c₁ = (3/5)^{1/2} (4/5)² / (2π√2)`. -/
noncomputable def c₁ : ℝ := √(3 / 5) * (4 / 5) ^ 2 / (2 * Real.pi * √2)

/-- Keiper's `y_m`, defined by `λ_m/m = ½(log m + γ - log 2π - 1) + y_m`. -/
noncomputable def keiperY (m : ℕ) : ℝ :=
  lam m / m - (1 / 2) * (Real.log m + Real.eulerMascheroniConstant - Real.log (2 * Real.pi) - 1)

/-- Voros's `S̄_n` (arXiv math/0506326, eq. (TN), TeX l.462–477): the contribution of the explicit
factor of `Ξ(s) = Γ(s/2) π^{-s/2} s (s-1) ζ(s)` to `λ_n` (Maślanka's trend),
`S̄_n = 1 - ((log 4π + γ)/2) n + Σ_{j=2}^n C(n,j) (-1)^j (1 - 2^{-j}) ζ(j)`. -/
noncomputable def vorosTrend (n : ℕ) : ℝ :=
  1 - (Real.log (4 * Real.pi) + Real.eulerMascheroniConstant) / 2 * n +
    ∑ j ∈ Finset.Icc 2 n, (n.choose j : ℝ) * (-1) ^ j * (1 - (2 : ℝ)⁻¹ ^ j) * (riemannZeta j).re

/-- The nontrivial zeros with ordinate `0 < γ ≤ T`. -/
def zerosUpTo (T : ℝ) : Set ℂ := {ρ | ρ ∈ NontrivialZeros ∧ 0 < ρ.im ∧ ρ.im ≤ T}

/-- The nontrivial zeros with ordinate `0 < γ < T`. -/
def zerosBelow (T : ℝ) : Set ℂ := {ρ | ρ ∈ NontrivialZeros ∧ 0 < ρ.im ∧ ρ.im < T}

/-- `N(T) = #{ρ : 0 < γ ≤ T}`, counted with multiplicity (right-continuous). -/
noncomputable def zeroCount (T : ℝ) : ℕ := ∑ᶠ ρ ∈ zerosUpTo T, mult ρ

/-- Titchmarsh's `N(T)` (Brent–Platt–Trudgian, arXiv 2009.13791, TeX l.263–275): at an ordinate
`T` the average `(N(T-0) + N(T+0))/2`; elsewhere `N(T)`. -/
noncomputable def zeroCountAvg (T : ℝ) : ℝ :=
  (((∑ᶠ ρ ∈ zerosBelow T, mult ρ : ℕ) : ℝ) + zeroCount T) / 2

/-- The Riemann–Siegel theta function `θ(t) = arg Γ(1/4 + it/2) - (t/2) log π` (continuous,
`θ(0) = 0`; e.g. arXiv 2004.00926, TeX l.65–69), written as the integral of its derivative
`d/du arg Γ(1/4 + iu/2) = ½ Re ψ(1/4 + iu/2)`, `ψ = Γ'/Γ` (Mathlib's `Complex.digamma`). -/
noncomputable def rsTheta (t : ℝ) : ℝ :=
  ∫ u in (0 : ℝ)..t, ((digamma (1 / 4 + (u / 2 : ℝ) * I)).re / 2 - Real.log Real.pi / 2)

/-- `S(T) = N(T) - θ(T)/π - 1`: by the Riemann–von Mangoldt formula `N(T) = θ(T)/π + 1 + S(T)`
(Titchmarsh §9.3) this is `π⁻¹ arg ζ(½ + iT)` (with Titchmarsh's conventions). -/
noncomputable def zetaS (T : ℝ) : ℝ := zeroCountAvg T - rsTheta T / Real.pi - 1

/-- Montgomery's form factor with multiplicities,
`F(α, T) = (2π/(T log T)) Σ_{0<γ,γ'≤T} T^{iα(γ-γ')} w(γ - γ')`, `w(u) = 4/(4 + u²)`. -/
noncomputable def montgomeryF (α T : ℝ) : ℂ :=
  (2 * Real.pi / (T * Real.log T) : ℝ) * ∑ᶠ ρ ∈ zerosUpTo T, ∑ᶠ ρ' ∈ zerosUpTo T,
    ((mult ρ * mult ρ' : ℕ) : ℂ) * Complex.exp (I * (α * (ρ.im - ρ'.im) * Real.log T)) *
      ((4 / (4 + (ρ.im - ρ'.im) ^ 2) : ℝ) : ℂ)

end LeanLi

