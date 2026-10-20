import Mathlib

/-!
# Round 40 — Fourier optimization and bootstrap: positivity of the spectral side is RH

| approach | gives | obstruction | verdict |
|----------|-------|-------------|---------|
| Cohn–Elkies LP bound in the explicit formula | ordinate bounds | off-line quadruples contribute with either sign | CIRCULAR |
| Viazovska-style magic functions (ξ = Mellin of θ) | — | all-window certificate ≡ RH | CIRCULAR |
| Carneiro–Milinovich–Soundararajan | prime-gap bounds under RH | unconditional version needs zero control | CIRCULAR |
| Beurling–Selberg extremal functions | density, large sieve | density ≠ exclusion | DEAD |
| explicit-formula bootstrap (Bonifacio–Mazáč–Pal style) | finite-window exclusion | spectral side not positive without RH | CIRCULAR |
| Lasserre / SOS hierarchy | finite windows | no uniform certificate | CIRCULAR |
| conformal-bootstrap analogy (FE = crossing) | — | no unitarity: zero terms not squares | CIRCULAR |
| de Branges extremal functions (Holt–Vaaler) | bounds on S(t), gaps | conditional on RH | CIRCULAR |

In a bootstrap the spectral side must be a sum of squares. For ζ, an on-line zero contributes |f̂(γ)|² ≥ 0,
but an off-line quadruple contributes a cross term f̂(γ+iδ)·conj f̂(γ-iδ), which can be negative:
-/

namespace Fayed.R40

/-- Cross terms are sign-indefinite (the missing "unitarity" of the ζ bootstrap). -/
theorem cross_term_indefinite : ∃ a b : ℂ, (a * (starRingEnd ℂ) b).re < 0 :=
  ⟨1, -1, by norm_num⟩

end Fayed.R40
