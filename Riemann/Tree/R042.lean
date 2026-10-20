import Mathlib

/-!
# Round 42 — beyond two thirds (Alpöge–Furman 2026, > 2/3 of zeros simple and on the line)

| question | answer | verdict |
|----------|--------|---------|
| new idea | a new mean-value estimate for a longer mollifier | (inspect zeta23 for assumed inputs) |
| limit of the method | not known to reach 100%; Farmer's θ = ∞ would | CIRCULAR |
| density one + zero-density ⇒ no off-line zeros? | rare ≠ none | DEAD |
| Levinson toolkit off the line | Selberg N(1/2+δ,T) ≪ T^{1-δ/4} log T | density only |
| consequences for Möbius / Weil form | none: off-line zeros dominate | DEAD |
| other L-functions; Euler-product blindness | DH-type combinations: 100% on the line, off-line zeros | DEAD |
| proportion → 1 ⇒ quasi-RH or Lindelöf? | no known implication | DEAD |
| reuse zeta23 for density / zero-free regions | partly (counting, bookkeeping) | open engineering |
-/

namespace Fayed.R42

/-- Zeros of ζ in the critical strip with 0 < Im s ≤ T (counted without multiplicity). -/
noncomputable def Nstrip (T : ℝ) : ℕ :=
  Set.ncard {s : ℂ | riemannZeta s = 0 ∧ 0 < s.re ∧ s.re < 1 ∧ 0 < s.im ∧ s.im ≤ T}

/-- Zeros on the critical line with 0 < t ≤ T. -/
noncomputable def Nline (T : ℝ) : ℕ :=
  Set.ncard {t : ℝ | riemannZeta ((1 / 2 : ℂ) + (t : ℂ) * Complex.I) = 0 ∧ 0 < t ∧ t ≤ T}

/-- A proportion c of the zeros up to height T lies on the line, for all large T. -/
def ProportionOnLine (c : ℝ) : Prop :=
  ∀ᶠ T in Filter.atTop, c * (Nstrip T : ℝ) ≤ (Nline T : ℝ)

/-- Density one: every proportion c < 1. Known not to imply RH for Davenport–Heilbronn-type
functions (Bombieri–Hejhal), so no Euler-product-blind argument can get RH from it. -/
def DensityOne : Prop := ∀ c : ℝ, c < 1 → ProportionOnLine c

end Fayed.R42
