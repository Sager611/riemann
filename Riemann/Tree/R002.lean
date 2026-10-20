import Riemann.Fayed.Theorem1

/-!
# Round 2 — frontier R2: a global σ-sensitive input passing the Davenport–Heilbronn filter

| id | direction | status |
|----|-----------|--------|
| R2.1 | explicit formula + Λ ≥ 0 + Landau | DEAD: positivity ⇒ singularity at abscissa 1 only; Ω-results |
| R2.2 | de la Vallée Poussin trig positivity | DEAD: pole at 1 dominates only near σ = 1 (dlvp_nonneg) |
| R2.3 | Perron on 1/ζ, Mertens M(x) | CIRCULAR: M(x) = O(x^{1/2+ε}) ≡ RH (MertensZeroFree target) |
| R2.4 | Hamburger conductor-1 rigidity | CIRCULAR: identifies ζ, says nothing on zeros |
| R2.5 | Jensen / Laguerre–Pólya / Lee–Yang | DEAD: fixed-degree hyperbolicity is local; needs uniformity |
| R2.6 | Fourier quasicrystals | DEAD: log-primes not uniformly discrete |
| R2.7 | prime sums Σ_p p^{-s} | DEAD: diverge for real σ ≤ 1 (prime_rpow_not_summable) |
| R2.8 | Connes trace formula / Weil positivity | CIRCULAR: positivity ≡ RH |
| R2.9 | double filter DH + Beurling, minimal axioms | DEAD: Selberg-class GRH no easier |
| R2.10 | xp-type Hamiltonians | DEAD: no self-adjoint realization; zeros put in by hand |

Meta-lesson: every known global input is either too weak (positivity near σ = 1)
or RH-equivalent. R3 looks at partial theorems that could be pushed.
-/

open Real

namespace Fayed.R2

/-- R2.2: de la Vallée Poussin, 3 + 4 cos θ + cos 2θ = 2(1 + cos θ)² ≥ 0. -/
theorem dlvp_nonneg (θ : ℝ) : 0 ≤ 3 + 4 * cos θ + cos (2 * θ) := by
  rw [cos_two_mul]; nlinarith [sq_nonneg (1 + cos θ)]

/-- R2.1: Λ(n) ≥ 0 — the free positivity of the explicit formula. -/
theorem vonMangoldt_nonneg' (n : ℕ) : 0 ≤ ArithmeticFunction.vonMangoldt n :=
  ArithmeticFunction.vonMangoldt_nonneg

/-- R2.7 no-go: Σ_p p^{-σ} diverges for every σ ≤ 1, so "prime sums converge on σ > 1/2"
is false as stated (only conditional or regularized versions are RH-equivalent). -/
theorem prime_rpow_not_summable {σ : ℝ} (hσ : σ ≤ 1) :
    ¬ Summable (fun p : Nat.Primes => (p : ℝ) ^ (-σ)) := by
  rw [Nat.Primes.summable_rpow]; linarith

/-- R2.3: Mertens bound |M(N)| ≤ C N^θ. -/
def MertensBound (θ : ℝ) : Prop :=
  ∃ C : ℝ, ∀ N : ℕ, |∑ n ∈ Finset.Icc 1 N, (ArithmeticFunction.moebius n : ℝ)| ≤ C * (N : ℝ) ^ θ

/-- R2.3, Littlewood's direction (true; formalization target): a Mertens bound gives a
zero-free half-plane. -/
def MertensZeroFree : Prop :=
  ∀ θ : ℝ, MertensBound θ → ∀ s : ℂ, θ < s.re → riemannZeta s ≠ 0

/-- R2.3, cited equivalence (stated, not asserted). -/
def RHIffMertens : Prop := RiemannHypothesis ↔ ∀ ε > 0, MertensBound (1 / 2 + ε)

end Fayed.R2
