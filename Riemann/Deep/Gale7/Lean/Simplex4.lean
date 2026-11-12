import Riemann.Deep.Gale7.Lean.Fam
import Mathlib.Analysis.Real.Sqrt

/-!
# Gale7: the regular 4-simplex of unit vectors in ℝ⁴

`u_i(k) = c_ik √ρ_k` with `c_ik = 1` (i < k), `−k` (i = k), `0` (i > k) and `ρ_k = 5/(4k(k+1)) = w_k/48`
(k = 1…4), as in Simplex3 and points8.txt. The Gram identity is checked on the integers `c_ik`, `w_k` by
kernel evaluation (`decide +kernel`).
-/

namespace Gale7

/-- integer coefficients c_ik for m = 4 -/
def c4 : Fin 5 → Fin 4 → ℤ :=
  ![![1, 1, 1, 1], ![-1, 1, 1, 1], ![0, -2, 1, 1], ![0, 0, -3, 1], ![0, 0, 0, -4]]
/-- w_k = 48 ρ_k -/
def w4 : Fin 4 → ℤ := ![30, 10, 5, 3]
/-- unit vertex vectors of a regular 4-simplex in ℝ⁴ -/
noncomputable def U4 : Fin 5 → Fin 4 → ℝ := fun i t => (c4 i t : ℝ) * Real.sqrt ((w4 t : ℝ) / 48)

theorem gram4Z : ∀ i k : Fin 5, ∑ t, c4 i t * c4 k t * w4 t = if i = k then 48 else -12 := by
  decide +kernel

theorem w4_nonneg : ∀ t : Fin 4, 0 ≤ w4 t := by
  decide +kernel

theorem U4_simplex : IsUnitSimplex U4 := by
  intro i k
  have e : ∀ t, U4 i t * U4 k t = ((c4 i t * c4 k t * w4 t : ℤ) : ℝ) * (1 / 48) := by
    intro t
    have hw : (0 : ℝ) ≤ (w4 t : ℝ) / 48 := div_nonneg (Int.cast_nonneg (w4_nonneg t)) (by norm_num)
    have h := Real.mul_self_sqrt hw
    simp only [U4]
    rw [Int.cast_mul, Int.cast_mul]
    calc (c4 i t : ℝ) * Real.sqrt ((w4 t : ℝ) / 48) * ((c4 k t : ℝ) * Real.sqrt ((w4 t : ℝ) / 48))
        = (c4 i t : ℝ) * (c4 k t) * (Real.sqrt ((w4 t : ℝ) / 48) * Real.sqrt ((w4 t : ℝ) / 48)) := by
          ring
      _ = (c4 i t : ℝ) * (c4 k t) * ((w4 t : ℝ) / 48) := by rw [h]
      _ = (c4 i t : ℝ) * (c4 k t) * (w4 t) * (1 / 48) := by ring
  rw [Finset.sum_congr rfl (fun t _ => e t), ← Finset.sum_mul, ← Int.cast_sum, gram4Z i k]
  split_ifs <;> norm_num

end Gale7
