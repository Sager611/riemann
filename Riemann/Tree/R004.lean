import Mathlib

/-!
# Round 4 — frontier R3.1: Weil positivity beyond the log 2 window

RH ⟺ W(f * f̃) ≥ 0 ∀ f. For supp f ⊆ (-L/2, L/2) only prime powers n ≤ e^L enter.
Proven: L ≤ log 2 (Yoshida 1992; Bombieri 2000; Connes–Consani 2021). Deep dive: Riemann/Deep/Weil/.

| id | sub-question | status |
|----|--------------|--------|
| R4.1 | literature: largest proven full window | frontier still log 2; prolate works are finite-prime → PROMISING (certify > log 2) |
| R4.2 | perturbation: p = 2 term vs archimedean gap | DEAD as pure norm bound: term norm ≈ 0.49, no known gap that large; needs numerics |
| R4.3 | sharp high-frequency cutoff | DEAD: compact support excludes exact band limits; use finite-rank projections |
| R4.4 | window ⇒ zero-free region / Li cutoff | DEAD: modulation reaches every height at fixed L |
| R4.5 | prolate interlacing ⇒ monotone λ_min(L) | DEAD: window change alters space and prime terms |
| R4.6 | function-field window induction | DEAD: Hodge index is global, not window by window |
| R4.7 | induction on primes | DEAD: prime insertions are indefinite rank-two forms (antitone only one way) |
| R4.8 | Lean feasibility | Mathlib has Schwartz/Fourier, digamma; no explicit formula: multi-step project |
-/

namespace Fayed.R4

/-- A windowed quadratic form: supp L = test functions supported in the window of length L. -/
structure WindowedForm (V : Type*) [AddCommGroup V] [Module ℝ V] where
  supp : ℝ → Submodule ℝ V
  mono : Monotone supp
  Q : V → ℝ

variable {V : Type*} [AddCommGroup V] [Module ℝ V]

/-- Positivity on window L. -/
def WindowedForm.Pos (W : WindowedForm V) (L : ℝ) : Prop := ∀ v ∈ W.supp L, 0 ≤ W.Q v

/-- Positivity passes to smaller windows (only this direction: R4.7). -/
theorem WindowedForm.pos_antitone (W : WindowedForm V) : Antitone W.Pos :=
  fun _ _ h hpos v hv => hpos v (W.mono h hv)

/-- R4.7 no-go, abstract: a nonnegative form minus a positive multiple of an indefinite
term need not stay nonnegative, so positivity cannot be carried across a prime insertion
without a quantitative gap. -/
theorem indefinite_update_breaks {a : ℝ} (ha : 0 < a) :
    ∃ Q B : ℝ → ℝ, (∀ v, 0 ≤ Q v) ∧ (∃ v, B v < 0) ∧ ¬ ∀ v, 0 ≤ Q v - a * B v := by
  refine ⟨fun _ => 0, fun v => -v, fun _ => le_rfl, ⟨1, by norm_num⟩, fun h => ?_⟩
  have := h (-1)
  nlinarith

end Fayed.R4
