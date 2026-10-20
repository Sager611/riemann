import Mathlib

/-!
# Round 13 — exotic RH-equivalents and their native tools

| id | equivalent | native tool | status |
|----|-----------|-------------|--------|
| R13.1 | Balazard–Saias–Yor ∫ log|ζ(½+it)| dt/(¼+t²) = 0 | Poisson–Jensen | CIRCULAR (integral = Σ_{Re ρ>½} log|ρ/(1-ρ)| ≥ 0) |
| R13.2 | Salem integral equation | Wiener–Tauberian (gave PNT on Re s = 1) | CIRCULAR: whole strip |
| R13.3 | Hardy-space cyclicity / outerness | Beurling's invariant subspaces | CIRCULAR |
| R13.4 | Riesz series O(x^{1/4+ε}) | Mellin / Tauberian | CIRCULAR |
| R13.5 | Volchkov double integral | Poisson–Jensen | CIRCULAR |
| R13.6 | Hilbert transform on the line (S(t)) | singular integrals | CIRCULAR |
| R13.7 | Bohr–Jessen value distribution | random Euler products | DEAD: not RH-sensitive |
| R13.8 | Jensen/Turán inequalities for Ξ | Pólya–Schur | DEAD: proven finite families do not force hyperbolicity |
| R13.9 | complexity / Möbius randomness | ergodic theory | DEAD: no RH-equivalence |
| R13.10 | Franel–Landau Farey discrepancy | SL₂(ℤ) spectral theory | CIRCULAR: spectral gap too weak |
-/

open Real

namespace Fayed.R13

/-- R13.1: Balazard–Saias–Yor (1999), cited: RH ⟺ this integral vanishes. -/
def BSYIntegralZero : Prop :=
  ∫ t : ℝ, Real.log ‖riemannZeta ((1 / 2 : ℂ) + t * Complex.I)‖ / (1 / 4 + t ^ 2) = 0

/-- R13.4: Riesz (1916), cited: RH ⟺ R(x) = O(x^{1/4+ε}) for every ε > 0. -/
noncomputable def rieszSeries (x : ℝ) : ℂ :=
  ∑' k : ℕ, ((-x) ^ (k + 1) : ℝ) / ((Nat.factorial k : ℂ) * riemannZeta (2 * ((k : ℂ) + 1)))

def RieszBound : Prop :=
  ∀ ε : ℝ, 0 < ε → ∃ C : ℝ, ∀ x : ℝ, 1 ≤ x → ‖rieszSeries x‖ ≤ C * x ^ (1 / 4 + ε)

end Fayed.R13
