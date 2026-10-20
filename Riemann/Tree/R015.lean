import Mathlib

/-!
# Round 15 — resurrection of quickly killed branches

| id | branch revisited | status after review |
|----|------------------|---------------------|
| R15.1 | Turán power sums on Σ x^ρ | DEAD: tails and Ω-oscillation absorb a dominant zero |
| R15.2 | Rankin–Selberg positivity (Goldfeld–Hoffstein–Lieman style) | DEAD: real-parameter positivity, no complex-zero leverage |
| R15.3 | AFE with |χ(ρ)| ≠ 1 across many shifts | DEAD: large sieve gives density, not exclusion (cancellation_exists) |
| R15.4 | self-repulsion vs Selberg's short-interval zeros | DEAD: abundance ≠ repulsion |
| R15.5 | amoebas / truncated ζ_N zeros as N → ∞ | DEAD: no stable limit near σ = 1/2 |
| R15.6 | Fayed with a B̄₂-shaped G, second computation | CIRCULAR: same identity twice |
| R15.7 | zeros vs subconvexity (Bourgain 13/84) | DEAD: upper bounds compatible with zeros |
| R15.8 | Connes–Consani–Moscovici prolate operators | CIRCULAR: spectral bridge is RH |
-/

namespace Fayed.R15

/-- R15.3: at a point, any χ and T admit a cancelling S — the relation S + χT + E = 0
constrains nothing without a lower bound on one of the polynomials. -/
theorem cancellation_exists (χ T E : ℂ) : ∃ S : ℂ, S + χ * T + E = 0 :=
  ⟨-(χ * T + E), by ring⟩

end Fayed.R15
