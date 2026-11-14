/-
E1Lean.Inst — UNCOMPILED (its imports exceed the 300 MB compile cap; do not import from Riemann.lean).
The intended instance of `Fns`: Mathlib's real functions; every field is the named Mathlib theorem
(names and binders checked against the Mathlib source of this toolchain, tag v4.35.0-rc2).
-/
import Mathlib.Analysis.SpecialFunctions.Trigonometric.Bounds
import Mathlib.Analysis.Complex.ExponentialBounds
import Mathlib.Analysis.SpecialFunctions.Log.Basic
import Mathlib.Analysis.SpecialFunctions.Pow.Real
import Riemann.Deep.SantiAudit.Lean.Arch

namespace E1Lean

/-- Mathlib's √, sin, cos, exp, log, rpow. -/
noncomputable def Fns.real : Fns where
  sqrt := Real.sqrt
  sin := Real.sin
  cos := Real.cos
  exp := Real.exp
  log := Real.log
  rpow := fun T θ => T ^ θ
  sqrt_nonneg := Real.sqrt_nonneg
  mul_self_sqrt := fun _ h => Real.mul_self_sqrt h
  sin_bound := fun _ h => Real.sin_bound h
  cos_bound := fun _ h => Real.cos_bound h
  sin_two_mul := Real.sin_two_mul
  cos_two_mul := Real.cos_two_mul
  sin_sq_add_cos_sq := Real.sin_sq_add_cos_sq
  sin_sub := Real.sin_sub
  sin_lt := fun _ h => Real.sin_lt h
  sin_gt_sub_cube := fun _ h => Real.sin_gt_sub_cube h
  exp_pos := Real.exp_pos
  exp_zero := Real.exp_zero
  exp_add := Real.exp_add
  exp_one_gt_d9 := Real.exp_one_gt_d9
  quadratic_le_exp_of_nonneg := fun _ h => Real.quadratic_le_exp_of_nonneg h
  log_le_iff_le_exp := fun _ _ h => Real.log_le_iff_le_exp h
  exp_one_lt_d9 := Real.exp_one_lt_d9
  exp_bound_div_one_sub_of_interval := fun _ h1 h2 => Real.exp_bound_div_one_sub_of_interval h1 h2

/-- Our cot is Mathlib's. -/
theorem cot_real (x : ℝ) : cot Fns.real x = Real.cot x := (Real.cot_eq_cos_div_sin x).symm

/-- Our c(θ) is Wang's c(θ) = 2 − θ/2 − cot(θ/√2)/√2 with Mathlib's functions. -/
theorem cW_real (θ : ℝ) : cW Fns.real θ = 2 - θ / 2 - Real.cot (θ / Real.sqrt 2) / Real.sqrt 2 := by
  simp only [cW, cot_real]; rfl

/-- Theorem E1 over Mathlib's real functions. -/
theorem theoremE1_real (I : Inputs) (h1 : H1 I) (h2 : H2 Fns.real I) (h3 : H3 Fns.real I) (h4 : H4 I)
    (θ : ℝ) (hθ1 : 1 / 2 < θ) (hθ2 : θ < 1) :
    LiminfGe (ratio Fns.real I.S I.N θ) (h Fns.real θ (κ₁ Fns.real I)) ∧
      h Fns.real (1 / 2) (κ₁ Fns.real I) ≤ h Fns.real θ (κ₁ Fns.real I) ∧
      (0.0055 : ℝ) < h Fns.real (1 / 2) (κ₁ Fns.real I) :=
  theoremE1 Fns.real I h1 h2 h3 h4 θ hθ1 hθ2

end E1Lean
