import Riemann.Deep.PNTError.General
import Riemann.Deep.PNTError.MertensMain
import Riemann.Deep.QRH.Bridge

/-!
# ψ(x) = x + O(x^{7/8+ε}) from OpenAI's zero-free half-plane Re s > 7/8

Main module of Riemann/Deep/PNTError/. The analytic work is zero-free-region-agnostic
(Basic.lean: Perron + contour shift; LogDeriv.lean: ζ'/ζ ≪ log t via Borel–Carathéodory;
General.lean: assembly). This file only feeds in
`Fayed.QRH.zeta_ne_zero_of_seven_eighths_lt_re` (the ported OpenAI proof, Riemann/Deep/QRH/).
-/

open Complex Filter Asymptotics
open scoped Real

namespace Fayed.PNTError

/-- OpenAI's theorem in this directory's vocabulary: ζ has no zero with Re s > 7/8. -/
theorem zeroFree_seven_eighths : ZeroFree (7 / 8) :=
  fun s hs => Fayed.QRH.zeta_ne_zero_of_seven_eighths_lt_re s hs

/-- **Rung 1.** For every ε > 0, ζ'/ζ(s) ≪_ε log|t| on Re s ≥ 7/8 + ε, |t| ≥ 2. -/
theorem logDeriv_bound_seven_eighths {ε : ℝ} (hε : 0 < ε) :
    ∃ C : ℝ, 0 ≤ C ∧ ∀ s : ℂ, 7 / 8 + ε ≤ s.re → 2 ≤ |s.im| →
      ‖deriv riemannZeta s / riemannZeta s‖ ≤ C * Real.log |s.im| :=
  logDeriv_bound_two zeroFree_seven_eighths (by norm_num) (by norm_num) hε

/-- **Rung 2: power-saving prime number theorem.** For every ε > 0 there is C with
|ψ(x) - x| ≤ C x^{7/8+ε} for all x ≥ 2. -/
theorem psi_sub_self_le {ε : ℝ} (hε : 0 < ε) :
    ∃ C : ℝ, ∀ x : ℝ, 2 ≤ x → |Chebyshev.psi x - x| ≤ C * x ^ (7 / 8 + ε) :=
  psi_error_of_zeroFree zeroFree_seven_eighths (by norm_num) (by norm_num) hε

/-- Rung 2 as a Landau statement: ψ(x) - x = O(x^{7/8+ε}). -/
theorem psi_sub_self_isBigO {ε : ℝ} (hε : 0 < ε) :
    (fun x : ℝ => Chebyshev.psi x - x) =O[atTop] (fun x : ℝ => x ^ (7 / 8 + ε)) := by
  obtain ⟨C, hC⟩ := psi_sub_self_le hε
  refine IsBigO.of_bound C ?_
  filter_upwards [eventually_ge_atTop (2 : ℝ)] with x hx
  rw [Real.norm_eq_abs, Real.norm_eq_abs, abs_of_nonneg (Real.rpow_nonneg (by linarith) _)]
  exact hC x hx

/-- **Rung 3: Mertens function.** For every ε > 0 there is C with |M(x)| ≤ C x^{7/8+ε}
for all x ≥ 2, where M(x) = Σ_{n ≤ x} μ(n). -/
theorem mertens_le {ε : ℝ} (hε : 0 < ε) :
    ∃ C : ℝ, ∀ x : ℝ, 2 ≤ x → |mertens x| ≤ C * x ^ (7 / 8 + ε) :=
  mertens_error_of_zeroFree zeroFree_seven_eighths (by norm_num) (by norm_num) hε

/-- 1/ζ(s) ≪_η |t|^η on Re s ≥ 7/8 + ε, |t| ≥ 4. -/
theorem inv_zeta_bound_seven_eighths {ε η : ℝ} (hε : 0 < ε) (hη : 0 < η) :
    ∃ C : ℝ, 0 ≤ C ∧ ∀ s : ℂ, 7 / 8 + ε ≤ s.re → 4 ≤ |s.im| →
      ‖(riemannZeta s)⁻¹‖ ≤ C * |s.im| ^ η :=
  inv_zeta_bound zeroFree_seven_eighths (by norm_num) (by norm_num) hε hη

end Fayed.PNTError
