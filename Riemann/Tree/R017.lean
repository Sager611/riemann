import Mathlib

/-!
# Round 17 — the axiom ladder: which properties allow off-line zeros?

| rung | class | known | status |
|------|-------|-------|--------|
| R17.1 | FE + Dirichlet series (no Euler product) | Davenport–Heilbronn 1936, Epstein: off-line zeros, even Re > 1 | DEAD rung |
| R17.2 | Euler product + continuation (no FE) | Diamond–Montgomery–Vorhauer 2006: zeros → Re s = 1 | DEAD rung |
| R17.3 | Euler product + FE (Selberg class) | no off-line zero known; Kaczorowski–Kulas 2006 only in the extended class | OPEN rung (≡ GRH for the class) |
| R17.4 | degree classification (Kaczorowski–Perelli) | degree 1 = Dirichlet L | CIRCULAR |
| R17.5 | Selberg orthonormality conjectures | control primes, not zeros | DEAD |
| R17.6 | random Euler products + FE | no RH theorem | DEAD |
| R17.7 | Hamburger + Euler product rigidity | identifies, does not locate | DEAD |
| R17.8 | finite Euler-factor changes | break FE; regions need all primes | DEAD |

Conclusion: the only rung with no counter-model is "Euler product + functional equation";
every proof must use both, and the known methods use the Euler product only near σ = 1.
-/

namespace Fayed.R17

/-- An off-critical zero of a function F in the strip. -/
def HasOffCriticalZero (F : ℂ → ℂ) : Prop :=
  ∃ s, F s = 0 ∧ 0 < s.re ∧ s.re < 1 ∧ s.re ≠ 1 / 2

end Fayed.R17
