import Riemann.Deep.Gale7.Lean.Crit
import Riemann.Deep.Gale7.Lean.Q7
import Riemann.Deep.Gale7.Lean.Simplex3
import Riemann.Deep.Gale7.Lean.Simplex4

/-!
# Gale7: Gale's conjecture fails in ℝ⁷

The member V7 of the family (NOTE §2), ℝ³ ⊕ ℝ⁴, parameters of family_params.txt row n = 7
(decimals as exact rationals), 69 points. Squared distances: rational class values `clsQ`, checked
by kernel evaluation in Q7 (`decide +kernel`): all ≤ D², all edge classes of `H(4, 5)` ≥ t, and σ7² D² < t.
-/

namespace Gale7

open Vx

/-- V7 ⊂ ℝ⁷ -/
noncomputable def pt7 : Vx 4 5 → Fin 7 → ℝ := famPt U3 U4 rA7 rB7 a7 b7 c7 d7 e7 f7
/-- the point set V7 -/
noncomputable def X7 : Finset (Fin 7 → ℝ) := Finset.univ.image pt7

theorem sqd7 (v w : Vx 4 5) : sqd (pt7 v) (pt7 w) =
    (clsQ 3 4 rA7 rB7 a7 b7 c7 d7 e7 f7 (orb v) (orb w) (decide (ix1 v = ix1 w))
      (decide (ix2 v = ix2 w)) : ℝ) :=
  sqd_famPt_cls (m₁ := 3) (m₂ := 4) U3_simplex U4_simplex _ _ _ _ _ _ _ _ v w

theorem hD7 (v w : Vx 4 5) : sqd (pt7 v) (pt7 w) ≤ (D7 : ℝ) := by
  rw [sqd7]
  exact Rat.cast_le.mpr (cls7_le _ _ _ _)

theorem hH7 (v w : Vx 4 5) (h : Adj v w) : (t7 : ℝ) ≤ sqd (pt7 v) (pt7 w) := by
  rw [sqd7]
  exact Rat.cast_le.mpr (cls7_ge _ _ _ _ (adj_adjD h))

theorem hσ7 : sigma 7 ^ 2 * (D7 : ℝ) < (t7 : ℝ) := by
  rw [sigma_odd_sq (by norm_num) (by norm_num)]
  have h1 : Real.sqrt ((((7 : ℕ) : ℝ) - 1) / (((7 : ℕ) : ℝ) + 3)) < 31 / 40 := by
    rw [Real.sqrt_lt' (by norm_num)]
    norm_num
  have hD : (0 : ℝ) < (D7 : ℝ) := by norm_num [D7]
  calc (1 / 2 + Real.sqrt ((((7 : ℕ) : ℝ) - 1) / (((7 : ℕ) : ℝ) + 3)) / 2) * (D7 : ℝ)
      < (1 / 2 + 31 / 40 / 2) * (D7 : ℝ) := by
        apply mul_lt_mul_of_pos_right _ hD
        linarith
    _ < (t7 : ℝ) := by norm_num [D7, t7]

/-- σ₇ = 0.9419651451… < 0.9419652 -/
theorem sigma7_lt : sigma 7 < 0.9419652 := by
  rw [sigma, if_neg (by norm_num), Real.sqrt_lt' (by norm_num)]
  have h : Real.sqrt ((((7 : ℕ) : ℝ) - 1) / (((7 : ℕ) : ℝ) + 3)) < 77459686 / 100000000 := by
    rw [Real.sqrt_lt' (by norm_num)]
    norm_num
  have e : (0.9419652 : ℝ) ^ 2 = 88729843801104 / 100000000000000 := by norm_num
  rw [e]
  linarith

/-- √(t₇/D₇) = 0.9433983553… > 0.9433983 -/
theorem ratio7_gt : (0.9433983 : ℝ) < Real.sqrt ((t7 : ℝ) / (D7 : ℝ)) := by
  rw [Real.lt_sqrt (by norm_num)]
  norm_num [t7, D7]

/-- **Gale's conjecture fails in ℝ⁷**: V7 has no partition into 8 parts of diameter ≤ σ7 · diam V7. -/
theorem not_galeSplit_X7 : ¬ GaleSplit 7 X7 :=
  not_galeSplit_of_H (by norm_num) pt7 (by norm_num [D7]) hD7 hH7 hσ7

/-- every partition of V7 into 8 parts has a part with two points at squared distance ≥ t7 -/
theorem X7_long (c : (Fin 7 → ℝ) → Fin 8) :
    ∃ x ∈ X7, ∃ y ∈ X7, c x = c y ∧ (t7 : ℝ) ≤ sqd x y :=
  long_pair (by norm_num) pt7 hH7 c

/-- diam V7 = √D7 (attained by Q₀, Q₁) -/
theorem diam_X7 : diam X7 = Real.sqrt (D7 : ℝ) := by
  apply le_antisymm (diam_image_le pt7 hD7)
  have hQ : sqd (pt7 (Q 0)) (pt7 (Q 1)) = (D7 : ℝ) := by
    rw [sqd7, Rat.cast_inj]
    decide +kernel
  have h := le_diam (X := X7) (Finset.mem_image_of_mem pt7 (Finset.mem_univ (Q 0)))
    (Finset.mem_image_of_mem pt7 (Finset.mem_univ (Q 1)))
  rw [eucl, hQ] at h
  exact h

end Gale7
