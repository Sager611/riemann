import Mathlib

/-!
# Round 11 — attacking the recurring obstructions

| id | obstruction | status |
|----|-------------|--------|
| R11.1 | σ-blindness of identities | HEURISTIC; FE alone allows zeros anywhere (fe_allows_any_zero) |
| R11.2 | Davenport–Heilbronn filter | passable (Euler product + Λ ≥ 0 give classical regions) but only near σ = 1 |
| R11.3 | positivity ≡ RH | DEAD: restricted positivity gives only classical regions |
| R11.4 | density cannot exclude one zero | THEOREM (trivial counting) |
| R11.5 | finite verification vs all heights | HEURISTIC: no explicit search bound known (bounded_search_decides) |
| R11.6 | Euler product only for Re s > 1 | converges on Re s = 1 (no zero info); interior ⇒ beyond RH |
| R11.7 | Mertens exponent gap | M(x) = Ω±(x^{1/2}) rules out δ > 1/2 only |
| R11.8 | no self-adjoint operator | CIRCULAR: spectrum by fiat is trivial (spectrum_by_fiat) |
-/

open scoped ComplexConjugate

namespace Fayed.R11

/-- R11.1: the functional-equation symmetry s ↦ 1 - s alone allows zeros anywhere:
(s - a)(s - (1 - a)) is symmetric with zeros a, 1 - a for every a. -/
theorem fe_allows_any_zero (a : ℂ) :
    (∀ s : ℂ, (1 - s - a) * (1 - s - (1 - a)) = (s - a) * (s - (1 - a))) ∧
    (a - a) * (a - (1 - a)) = 0 ∧ ((1 - a) - a) * ((1 - a) - (1 - a)) = 0 := by
  refine ⟨fun s => by ring, by ring, by ring⟩

/-- R11.5: a bounded counterexample search would decide RH — the missing piece is the bound. -/
theorem bounded_search_decides {P : ℕ → Prop} {RH : Prop} (B : ℕ)
    (hb : ¬ RH → ∃ n, n ≤ B ∧ ¬ P n) (hok : ∀ n, n ≤ B → P n) : RH := by
  by_contra h
  obtain ⟨n, hn, hP⟩ := hb h
  exact hP (hok n hn)

/-- R11.8: any finite real spectrum is realized by a Hermitian (diagonal) matrix with the
right power traces, so "an operator with the zeros as spectrum" has no content by itself. -/
theorem spectrum_by_fiat {n : ℕ} (γ : Fin n → ℝ) :
    ∃ H : Matrix (Fin n) (Fin n) ℂ, H.IsHermitian ∧
      ∀ k : ℕ, (H ^ k).trace = ∑ i, ((γ i : ℂ)) ^ k := by
  refine ⟨Matrix.diagonal (fun i => (γ i : ℂ)), ?_, fun k => ?_⟩
  · rw [Matrix.IsHermitian, Matrix.diagonal_conjTranspose]
    congr 1
    funext i
    simp
  · rw [Matrix.diagonal_pow, Matrix.trace_diagonal]
    rfl

end Fayed.R11
