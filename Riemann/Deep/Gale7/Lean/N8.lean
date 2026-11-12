import Riemann.Deep.Gale7.Lean.Crit
import Riemann.Deep.Gale7.Lean.Q8
import Riemann.Deep.Gale7.Lean.Simplex4

/-!
# Gale7: Gale's conjecture fails in ℝ⁸

The member V8 of the family (NOTE §2), ℝ⁴ ⊕ ℝ⁴, parameters of family_params.txt row n = 8
(decimals as exact rationals), 85 points. Squared distances: rational class values `clsQ`, checked
by kernel evaluation in Q8 (`decide +kernel`): all ≤ D², all edge classes of `H(5, 5)` ≥ t, and σ8² D² < t.
-/

namespace Gale7

open Vx

/-- V8 ⊂ ℝ⁸ -/
noncomputable def pt8 : Vx 5 5 → Fin 8 → ℝ := famPt U4 U4 rA8 rB8 a8 b8 c8 d8 e8 f8
/-- the point set V8 -/
noncomputable def X8 : Finset (Fin 8 → ℝ) := Finset.univ.image pt8

theorem sqd8 (v w : Vx 5 5) : sqd (pt8 v) (pt8 w) =
    (clsQ 4 4 rA8 rB8 a8 b8 c8 d8 e8 f8 (orb v) (orb w) (decide (ix1 v = ix1 w))
      (decide (ix2 v = ix2 w)) : ℝ) :=
  sqd_famPt_cls (m₁ := 4) (m₂ := 4) U4_simplex U4_simplex _ _ _ _ _ _ _ _ v w

theorem hD8 (v w : Vx 5 5) : sqd (pt8 v) (pt8 w) ≤ (D8 : ℝ) := by
  rw [sqd8]
  exact Rat.cast_le.mpr (cls8_le _ _ _ _)

theorem hH8 (v w : Vx 5 5) (h : Adj v w) : (t8 : ℝ) ≤ sqd (pt8 v) (pt8 w) := by
  rw [sqd8]
  exact Rat.cast_le.mpr (cls8_ge _ _ _ _ (adj_adjD h))

theorem hσ8 : sigma 8 ^ 2 * (D8 : ℝ) < (t8 : ℝ) := by
  rw [sigma_even_sq (by norm_num)]
  norm_num [D8, t8]

/-- σ₈ = √(9/10) = 0.9486832981… < 0.9486834 -/
theorem sigma8_lt : sigma 8 < 0.9486834 := by
  rw [sigma, if_pos (by norm_num), Real.sqrt_lt' (by norm_num)]
  norm_num

/-- √(t₈/D₈) = 0.9499470859… > 0.949947 -/
theorem ratio8_gt : (0.949947 : ℝ) < Real.sqrt ((t8 : ℝ) / (D8 : ℝ)) := by
  rw [Real.lt_sqrt (by norm_num)]
  norm_num [t8, D8]

/-- **Gale's conjecture fails in ℝ⁸**: V8 has no partition into 9 parts of diameter ≤ σ8 · diam V8. -/
theorem not_galeSplit_X8 : ¬ GaleSplit 8 X8 :=
  not_galeSplit_of_H (by norm_num) pt8 (by norm_num [D8]) hD8 hH8 hσ8

/-- every partition of V8 into 9 parts has a part with two points at squared distance ≥ t8 -/
theorem X8_long (c : (Fin 8 → ℝ) → Fin 9) :
    ∃ x ∈ X8, ∃ y ∈ X8, c x = c y ∧ (t8 : ℝ) ≤ sqd x y :=
  long_pair (by norm_num) pt8 hH8 c

/-- diam V8 = √D8 (attained by Q₀, Q₁) -/
theorem diam_X8 : diam X8 = Real.sqrt (D8 : ℝ) := by
  apply le_antisymm (diam_image_le pt8 hD8)
  have hQ : sqd (pt8 (Q 0)) (pt8 (Q 1)) = (D8 : ℝ) := by
    rw [sqd8, Rat.cast_inj]
    decide +kernel
  have h := le_diam (X := X8) (Finset.mem_image_of_mem pt8 (Finset.mem_univ (Q 0)))
    (Finset.mem_image_of_mem pt8 (Finset.mem_univ (Q 1)))
  rw [eucl, hQ] at h
  exact h

end Gale7
