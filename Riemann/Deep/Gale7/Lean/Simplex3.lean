import Riemann.Deep.Gale7.Lean.Fam
import Mathlib.Analysis.Real.Sqrt

/-!
# Gale7: the regular tetrahedron of unit vectors in ℝ³

`u_i(k) = c_ik √ρ_k` with `c_ik = 1` (i < k), `−k` (i = k), `0` (i > k), `ρ_k = (m+1)/(m k (k+1))` (k = 1…m);
here m = 3, ρ = (2/3, 2/9, 1/9) = w/9: the coordinates of points7_sym69.txt before rounding. The Gram identity is checked on the integers `c_ik`, `w_k` by
kernel evaluation (`decide +kernel`).
-/

namespace Gale7

/-- integer coefficients c_ik for m = 3 -/
def c3 : Fin 4 → Fin 3 → ℤ := ![![1, 1, 1], ![-1, 1, 1], ![0, -2, 1], ![0, 0, -3]]
/-- w_k = 9 ρ_k -/
def w3 : Fin 3 → ℤ := ![6, 2, 1]
/-- unit vertex vectors of a regular tetrahedron in ℝ³ -/
noncomputable def U3 : Fin 4 → Fin 3 → ℝ := fun i t => (c3 i t : ℝ) * Real.sqrt ((w3 t : ℝ) / 9)

theorem gram3Z : ∀ i k : Fin 4, ∑ t, c3 i t * c3 k t * w3 t = if i = k then 9 else -3 := by
  decide +kernel

theorem w3_nonneg : ∀ t : Fin 3, 0 ≤ w3 t := by
  decide +kernel

theorem U3_simplex : IsUnitSimplex U3 := by
  intro i k
  have e : ∀ t, U3 i t * U3 k t = ((c3 i t * c3 k t * w3 t : ℤ) : ℝ) * (1 / 9) := by
    intro t
    have hw : (0 : ℝ) ≤ (w3 t : ℝ) / 9 := div_nonneg (Int.cast_nonneg (w3_nonneg t)) (by norm_num)
    have h := Real.mul_self_sqrt hw
    simp only [U3]
    rw [Int.cast_mul, Int.cast_mul]
    calc (c3 i t : ℝ) * Real.sqrt ((w3 t : ℝ) / 9) * ((c3 k t : ℝ) * Real.sqrt ((w3 t : ℝ) / 9))
        = (c3 i t : ℝ) * (c3 k t) * (Real.sqrt ((w3 t : ℝ) / 9) * Real.sqrt ((w3 t : ℝ) / 9)) := by
          ring
      _ = (c3 i t : ℝ) * (c3 k t) * ((w3 t : ℝ) / 9) := by rw [h]
      _ = (c3 i t : ℝ) * (c3 k t) * (w3 t) * (1 / 9) := by ring
  rw [Finset.sum_congr rfl (fun t _ => e t), ← Finset.sum_mul, ← Int.cast_sum, gram3Z i k]
  split_ifs <;> norm_num

end Gale7
