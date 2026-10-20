import Riemann.Tree.R009
import Riemann.Tree.R007

/-!
# Lagarias positivity ⇒ RH (sufficiency direction)

Lagarias (1999): RH ⟺ Re ξ'(s)/ξ(s) > 0 on Re s > 1/2. Here: the direction ⇐.

* `xi_eq`: ξ(s) = s(s-1)/2 · Λ(s) for s ≠ 0, 1 (Mathlib's Λ = Λ₀ - 1/s - 1/(1-s)).
* `xi_zero_iff`: on Re s > 0, s ≠ 1, ξ(s) = 0 ↔ ζ(s) = 0 (ζ = Λ / Γℝ, Γℝ ≠ 0).
* `RH_of_LagariasPositivity`: a zero ρ with Re ρ > 1/2 gives ξ(ρ) = 0, so in Lean
  ξ'(ρ)/ξ(ρ) = ξ'(ρ)/0 = 0 has real part 0, contradicting positivity; conclude via
  `Fayed.R7.RH_iff_QuasiRH_half`.
-/

open Complex

namespace Fayed.Lagarias

open Fayed.R9

/-- ξ(s) = s(s-1)/2 · Λ(s) away from the poles of Λ. -/
theorem xi_eq (s : ℂ) (h0 : s ≠ 0) (h1 : s ≠ 1) :
    xi s = s * (s - 1) / 2 * completedRiemannZeta s := by
  have h1' : (1 : ℂ) - s ≠ 0 := sub_ne_zero.mpr (Ne.symm h1)
  rw [xi, completedRiemannZeta_eq]
  field_simp
  ring

/-- On Re s > 0, s ≠ 1: zeros of ξ are exactly zeros of ζ. -/
theorem xi_zero_iff {s : ℂ} (hs : 0 < s.re) (h1 : s ≠ 1) :
    xi s = 0 ↔ riemannZeta s = 0 := by
  have h0 : s ≠ 0 := by rintro rfl; simp at hs
  have h1' : s - 1 ≠ 0 := sub_ne_zero.mpr h1
  have hG : Gammaℝ s ≠ 0 := Gammaℝ_ne_zero_of_re_pos hs
  rw [xi_eq s h0 h1, riemannZeta_def_of_ne_zero h0, div_eq_zero_iff, or_iff_left hG]
  constructor
  · intro h
    rcases mul_eq_zero.mp h with h | h
    · exfalso
      rcases div_eq_zero_iff.mp h with h | h
      · exact (mul_ne_zero h0 h1') h
      · exact two_ne_zero h
    · exact h
  · intro h; rw [h, mul_zero]

/-- Lagarias positivity implies RH.

The Lean statement `LagariasPositivity` is used only at zeros of ξ, where Lean's
convention x / 0 = 0 makes the real part 0, already contradicting `0 < _`. The
mathematical content is stronger: at a zero ρ of order m, ξ'/ξ has a simple pole with
residue m > 0, so Re ξ'/ξ takes negative values arbitrarily close to ρ (on the side
Re(s - ρ) < 0); positivity on the open half-plane therefore fails near ρ with genuine
values, not just the junk value at ρ itself. -/
theorem RH_of_LagariasPositivity (h : LagariasPositivity) : RiemannHypothesis := by
  refine Fayed.R7.RH_iff_QuasiRH_half.mpr ?_
  intro s hs hs1 hz
  have hs' : 1 / 2 < s.re := by linarith
  have hxi : xi s = 0 := (xi_zero_iff (by linarith) hs1).mpr hz
  have := h s hs'
  rw [hxi, div_zero, zero_re] at this
  exact lt_irrefl _ this

end Fayed.Lagarias

#print axioms Fayed.Lagarias.xi_eq
#print axioms Fayed.Lagarias.xi_zero_iff
#print axioms Fayed.Lagarias.RH_of_LagariasPositivity
