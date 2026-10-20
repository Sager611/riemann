import Riemann.Deep.Mertens.ZeroFree
import Riemann.Deep.Redheffer.Det
import Riemann.Tree.R007

/-!
# Formal backbone: machine-checked reductions of RH

* Littlewood: Mertens bounds M(x) = O(x^{1/2+ε}) for all ε imply RH.
* Redheffer: det R_n = M(n), so determinant bounds imply RH.

Pieces: Fayed.R2.mertensZeroFree (Deep/Mertens), Fayed.Redheffer.det_redheffer (Deep/Redheffer),
Fayed.R7.RH_iff_QuasiRH_half (Tree/R007, via Mathlib's functional equation).
The converses (RH ⇒ bounds) are the hard directions and are not formalized.
-/

namespace Fayed.Backbone

open Fayed.R2 Fayed.R7

/-- Littlewood's criterion, sufficiency: Mertens bounds for every exponent > 1/2 imply RH. -/
theorem RH_of_mertens (h : ∀ ε > 0, MertensBound (1 / 2 + ε)) : RiemannHypothesis := by
  rw [RH_iff_QuasiRH_half]
  intro s hs _
  have hε : 0 < (s.re - 1 / 2) / 2 := by linarith
  exact Fayed.R2.mertensZeroFree (1 / 2 + (s.re - 1 / 2) / 2) (h _ hε) s (by linarith)

/-- Redheffer's criterion, sufficiency: |det R_n| = O(n^{1/2+ε}) for every ε > 0 implies RH. -/
theorem RH_of_redheffer
    (h : ∀ ε : ℝ, 0 < ε → ∃ C : ℝ, ∀ n : ℕ, 0 < n →
      |((Fayed.Redheffer.redheffer n).det : ℝ)| ≤ C * (n : ℝ) ^ (1 / 2 + ε)) :
    RiemannHypothesis := by
  apply RH_of_mertens
  intro ε hε
  obtain ⟨C, hC⟩ := h ε hε
  refine ⟨C, fun N => ?_⟩
  rcases Nat.eq_zero_or_pos N with rfl | hN
  · have h0 : ((0 : ℕ) : ℝ) ^ (1 / 2 + ε) = 0 := by
      rw [Nat.cast_zero]; exact Real.zero_rpow (by positivity)
    rw [h0, mul_zero]
    simp
  · have := hC N hN
    rw [Fayed.Redheffer.det_redheffer N hN] at this
    push_cast at this
    exact this

end Fayed.Backbone

#print axioms Fayed.Backbone.RH_of_mertens
#print axioms Fayed.Backbone.RH_of_redheffer
