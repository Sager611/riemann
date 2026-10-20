import Riemann.Fayed.Theorem1

/-!
# Round 3 — frontier R3: push existing partial theorems; finite faces of RH

| id | direction | status |
|----|-----------|--------|
| R3.1 | Weil positivity for small support (Yoshida 1992; Connes–Consani 2021) | PROMISING → R4: window beyond log 2, prime 2 enters |
| R3.2 | zeros of partial sums ζ_N (Turán; Montgomery 1983; Gonek–Ledoan 2009) | DEAD as route: ζ_N zeros ≠ ζ zeros |
| R3.3 | Robin σ(n) < e^γ n log log n via colossally abundant n | CIRCULAR: CA control needs RH-quality prime counts |
| R3.4 | Bagchi strong recurrence / Voronin universality | CIRCULAR |
| R3.5 | Guth–Maynard density → zero-free by amplification | DEAD: density bounds never exclude one zero |
| R3.6 | Deligne tensor-power trick | DEAD: GL₁ tensor powers trivial |
| R3.7 | Lean-certified zero verification (Turing) | DEAD as route (finite T); RHUpTo stated |
| R3.8 | o-minimality / Pila–Wilkie | DEAD: zeros are not algebraic points |
| R3.9 | Connes–Consani Riemann–Roch for Spec ℤ | CIRCULAR: no Hodge index on Spec ℤ × Spec ℤ |
| R3.10 | adelic renormalized heat flow | DEAD: no functorial finite-prime flow |
-/

open Real

namespace Fayed.R3

/-- R3.7: no off-line zero up to height T. Platt (2021): holds for T = 3·10¹² (cited). -/
def RHUpTo (T : ℝ) : Prop :=
  ∀ s : ℂ, 0 < s.re → s.re < 1 → |s.im| ≤ T → riemannZeta s = 0 → s.re = 1 / 2

/-- RH ⟺ RHUpTo T for all T, given that nontrivial zeros lie in the strip. -/
theorem RH_iff_forall_RHUpTo
    (hstrip : ∀ s : ℂ, riemannZeta s = 0 → (¬∃ n : ℕ, s = -2 * (n + 1)) → s ≠ 1 →
      0 < s.re ∧ s.re < 1) :
    RiemannHypothesis ↔ ∀ T, RHUpTo T := by
  constructor
  · intro h T s h0 h1 _ hz
    refine h s hz ?_ ?_
    · rintro ⟨n, rfl⟩
      have : (-2 * ((n : ℂ) + 1)).re = -2 * ((n : ℝ) + 1) := by simp
      rw [this] at h0
      have : (0 : ℝ) ≤ n := Nat.cast_nonneg n
      linarith
    · rintro rfl
      simp at h1
  · intro h s hz hnt hs1
    obtain ⟨h0, h1⟩ := hstrip s hz hnt hs1
    exact h |s.im| s h0 h1 le_rfl hz

/-- R3.3: Robin's inequality at n. -/
def RobinIneq (n : ℕ) : Prop :=
  (ArithmeticFunction.sigma 1 n : ℝ) < exp eulerMascheroniConstant * n * log (log n)

/-- R3.3, Robin (1984), cited: RH ⟺ Robin's inequality for all n > 5040 (stated, not asserted). -/
def RobinCriterion : Prop := RiemannHypothesis ↔ ∀ n : ℕ, 5040 < n → RobinIneq n

end Fayed.R3
