import Mathlib

/-!
# Round 36 — random models: none transfers to ζ

| model | proven | verdict |
|-------|--------|---------|
| Wintner 1944 random multiplicative ±1 | O(x^{1/2+ε}) almost surely | DEAD (μ is one sample) |
| Harper 2020 | better than square-root cancellation in L¹ | DEAD |
| Selberg CLT, Bohr–Jessen | value distribution | DEAD (distribution ≠ exclusion) |
| Voronin universality, Bagchi 1981 | RH ⇔ strong recurrence | CIRCULAR (definition below) |
| Saksman–Webb GMC | chaos limit on the line | DEAD |
| Fyodorov–Hiary–Keating maxima | upper bounds on the line | DEAD |
| Montgomery / Rudnick–Sarnak statistics | mostly under RH | CIRCULAR |
| Odlyzko–te Riele, Ng 2004 | Mertens conjecture false; M(x)/√x law under RH | CIRCULAR (Landau: off-line zero ⇒ M large, known) |
-/

namespace Fayed.R36

/-- Bagchi's strong recurrence: ζ approximates its own vertical shifts on compact subsets of the
right half of the critical strip, for a set of shifts of positive lower density. -/
def StronglyRecurrent : Prop :=
  ∀ K : Set ℂ, IsCompact K → (∀ s ∈ K, 1 / 2 < s.re ∧ s.re < 1) → ∀ ε : ℝ, 0 < ε →
    0 < Filter.liminf (fun T : ℝ => (MeasureTheory.volume {τ : ℝ | τ ∈ Set.Icc 0 T ∧
      ∀ s ∈ K, ‖riemannZeta (s + (τ : ℂ) * Complex.I) - riemannZeta s‖ < ε}).toReal / T)
      Filter.atTop

/-- Bagchi (1981): strong recurrence is equivalent to RH (cited statement, not proved here). -/
def BagchiCriterion : Prop := StronglyRecurrent ↔ RiemannHypothesis

end Fayed.R36
