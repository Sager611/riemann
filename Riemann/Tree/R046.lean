import Riemann.Fayed.StripReflect

/-!
# Round 46 — deformation paths from a provable-RH object to ζ: none controllable

| object | RH status | path to ζ | what breaks | verdict |
|--------|-----------|-----------|-------------|---------|
| Selberg zeta, compact surface | Selberg 1956 (up to exceptional eigenvalues) | none canonical | degeneration, zeros from infinity | DEAD |
| Ihara zeta, Ramanujan graph | Hashimoto 1989 | none | degree / FE change | DEAD |
| q-zeta (KKW, Cherednik), q → 1 | unproved | q → 1 | no RH endpoint, no FE | DEAD |
| CUE characteristic polynomials | tautological (unitary) | statistics only | no zero-preserving limit | DEAD |
| curves over F_q, q → 1 / g → ∞ | Weil 1948 per curve | statistics only | no zero tracking | DEAD |
| de Bruijn–Newman H_t, t: 0.2 → 0 | RH for t ≥ 0.2 | exact path | collisions before 0 ⇔ ¬RH | CIRCULAR |
| partial sums / AFE truncations | no | N → ∞ | no exact FE, zeros unpaired | DEAD |
| Epstein zeta | h ≥ 2: off-line zeros (Lamzouri 2019) | lattice degeneration | wrong endpoint | DEAD |

Mechanism behind every deformation argument (machine-checked below): an off-line zero of ζ in the strip
has a distinct partner 1 - conj ρ, so a simple zero on the line cannot leave it alone; zeros leave only
by colliding — exactly what de Bruijn–Newman flow does for t < Λ.
-/

namespace Fayed.R46

/-- Off the critical line, the reflected point is a different point. -/
theorem reflect_ne (ρ : ℂ) (h : ρ.re ≠ 1 / 2) : 1 - (starRingEnd ℂ) ρ ≠ ρ := by
  intro heq
  apply h
  have := congrArg Complex.re heq
  simp at this
  linarith

/-- An off-line zero of ζ in the critical strip has a distinct partner zero. -/
theorem zeta_offline_pair (ρ : ℂ) (hz : riemannZeta ρ = 0) (h0 : 0 < ρ.re) (h1 : ρ.re < 1)
    (h : ρ.re ≠ 1 / 2) :
    1 - (starRingEnd ℂ) ρ ≠ ρ ∧ riemannZeta (1 - (starRingEnd ℂ) ρ) = 0 :=
  ⟨reflect_ne ρ h, Fayed.reflect_holds ρ hz h0 h1⟩

end Fayed.R46

#print axioms Fayed.R46.zeta_offline_pair
