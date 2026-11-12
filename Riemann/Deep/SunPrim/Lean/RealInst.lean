import Riemann.Deep.SunPrim.Lean.Conj310
import Riemann.Deep.SunPrim.Lean.A228728
import Riemann.Deep.SunPrim.Lean.Roots
import Mathlib.Basic.Real.Basic

/-! The refutations stated for an arbitrary field, specialised to `ℝ` (the sources' setting). -/

namespace SunLean

theorem not_conj_3_10_ii_real : ¬ Conj_3_10_ii ℝ := not_conj_3_10_ii ℝ

theorem not_conj_A228728_real : ¬ Conj_A228728 ℝ := not_conj_A228728 ℝ

/-- For every positive-root function on `ℝ`, e.g. `x ↦ x ^ (1 / n : ℝ)`. -/
theorem not_conj_4_30_iii_real (rt : ℕ → ℝ → ℝ) (hrt : IsRootFn ℝ rt) : ¬ Conj_4_30_iii ℝ rt :=
  not_conj_4_30_iii ℝ rt hrt

/-- For every positive-root function on `ℝ`, e.g. `x ↦ x ^ (1 / n : ℝ)`. -/
theorem not_conj_3_1_real (rt : ℕ → ℝ → ℝ) (hrt : IsRootFn ℝ rt) : ¬ Conj_3_1 ℝ rt :=
  not_conj_3_1 ℝ rt hrt

end SunLean

#print axioms SunLean.not_conj_3_10_ii_real
#print axioms SunLean.not_conj_A228728_real
#print axioms SunLean.not_conj_4_30_iii_real
#print axioms SunLean.not_conj_3_1_real
