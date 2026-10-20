import Mathlib

/-!
# Round 5 — breadth: other frontiers, disproof side, meta

| id | direction | status |
|----|-----------|--------|
| R5.1 | Levinson–Conrey mollifiers (PRZZ 2019: > 41.67% on the line) | OPEN but insufficient: 100% ≠ RH; long mollifiers need off-diagonal control |
| R5.2 | Li coefficients, finite positivity | DEAD: finite λ_n mix all zeros; no zero exclusion |
| R5.3 | disproof by computation (> 3·10¹²) | DEAD: no targeted height; scans cannot certify off-line zeros |
| R5.4 | logic: Π⁰₁, independence | DEAD: independence unprovable from inside; gives no proof |
| R5.5 | Lax–Phillips scattering on the modular surface | CIRCULAR: causality ≡ RH |
| R5.6 | de Branges spaces / Suzuki screw functions | DEAD: Conrey–Li refute de Branges positivity |
| R5.7 | F₁ absolute geometry | DEAD: no surface with Frobenius + Hodge index |
| R5.8 | optimal transport / Fourier interpolation | CIRCULAR: reduces to Weil positivity |
-/

open Real

namespace Fayed.R5

/-- An off-line zero in the critical strip. -/
def OffLineZero : Prop := ∃ s : ℂ, 0 < s.re ∧ s.re < 1 ∧ s.re ≠ 1 / 2 ∧ riemannZeta s = 0

/-- RH ⟺ no off-line zero, given that nontrivial zeros lie in the strip. -/
theorem RH_iff_not_offLineZero
    (hstrip : ∀ s : ℂ, riemannZeta s = 0 → (¬∃ n : ℕ, s = -2 * (n + 1)) → s ≠ 1 →
      0 < s.re ∧ s.re < 1) :
    RiemannHypothesis ↔ ¬ OffLineZero := by
  constructor
  · rintro h ⟨s, h0, h1, hne, hz⟩
    refine hne (h s hz ?_ ?_)
    · rintro ⟨n, rfl⟩
      have : (-2 * ((n : ℂ) + 1)).re = -2 * ((n : ℝ) + 1) := by simp
      rw [this] at h0
      have : (0 : ℝ) ≤ n := Nat.cast_nonneg n
      linarith
    · rintro rfl
      simp at h1
  · intro h s hz hnt hs1
    obtain ⟨h0, h1⟩ := hstrip s hz hnt hs1
    by_contra hne
    exact h ⟨s, h0, h1, hne, hz⟩

/-- R5.4: Lagarias (2002), cited: RH ⟺ σ(n) ≤ H_n + exp(H_n) log H_n for all n ≥ 1
(stated, not asserted). A Π⁰₁ statement. -/
def LagariasCriterion : Prop :=
  RiemannHypothesis ↔ ∀ n : ℕ, 1 ≤ n →
    (ArithmeticFunction.sigma 1 n : ℝ) ≤ (harmonic n : ℝ) + exp (harmonic n) * log (harmonic n)

end Fayed.R5
