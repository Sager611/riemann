import Riemann.Deep.Gale7.Lean.Desc
import Mathlib.Basic.Real.Basic
import Mathlib.Algebra.BigOperators.Fin
import Mathlib.Data.Rat.Cast.Order
import Mathlib.Data.Rat.Cast.Lemmas
import Mathlib.Algebra.Module.Pi
import Mathlib.Tactic.Ring

/-!
# Gale7: the family Vₙ (NOTE §2) and its squared distances

`sqd`: squared Euclidean distance in ℝⁿ (coordinates `Fin n → ℝ`). `famPt U W rA rB a b c d e f`: the points
`P_i = r_A u_i`, `Q_j = r_B v_j`, `A_ij = a u_i − b v_j`, `B_ij = −c u_i + d v_j`, `C_ij = e u_i + f v_j` of
ℝ^{m₁} ⊕ ℝ^{m₂} (`Fin.append`), `U`, `W` the unit vertex vectors of regular simplices (`IsUnitSimplex`).
`sqd_famPt`: the class formula; `sqd_famPt_cls`: for rational parameters, the rational class value `clsQ`.
-/

namespace Gale7

open Vx

/-- Squared Euclidean distance in ℝⁿ. -/
noncomputable def sqd {n : ℕ} (x y : Fin n → ℝ) : ℝ := ∑ k, (x k - y k) ^ 2

theorem sqd_nonneg {n : ℕ} (x y : Fin n → ℝ) : 0 ≤ sqd x y :=
  Finset.sum_nonneg fun k _ => sq_nonneg (x k - y k)

theorem sqd_self {n : ℕ} (x : Fin n → ℝ) : sqd x x = 0 := by
  simp [sqd]

/-- `U i` (i = 0…m) are unit vectors with ⟨u_i, u_k⟩ = −1/m for i ≠ k: the vertices of a regular simplex
centred at 0 with circumradius 1 (NOTE §2). -/
def IsUnitSimplex {m : ℕ} (U : Fin (m + 1) → Fin m → ℝ) : Prop :=
  ∀ i k, ∑ t, U i t * U k t = if i = k then 1 else -1 / (m : ℝ)

variable {m₁ m₂ : ℕ}

section
variable (U : Fin (m₁ + 1) → Fin m₁ → ℝ) (W : Fin (m₂ + 1) → Fin m₂ → ℝ) (rA rB a b c d e f : ℝ)

/-- The family Vₙ (NOTE §2), ℝⁿ = ℝ^{m₁} ⊕ ℝ^{m₂}. -/
noncomputable def famPt : Vx (m₁ + 1) (m₂ + 1) → Fin (m₁ + m₂) → ℝ
  | P i => Fin.append (rA • U i) 0
  | Q j => Fin.append 0 (rB • W j)
  | A i j => Fin.append (a • U i) ((-b) • W j)
  | B i j => Fin.append ((-c) • U i) (d • W j)
  | C i j => Fin.append (e • U i) (f • W j)

end

/-- coefficient of the first-block vector -/
def cf1 (rA a c e : ℝ) : Vx (m₁ + 1) (m₂ + 1) → ℝ
  | P _ => rA
  | Q _ => 0
  | A _ _ => a
  | B _ _ => -c
  | C _ _ => e

/-- coefficient of the second-block vector -/
def cf2 (rB b d f : ℝ) : Vx (m₁ + 1) (m₂ + 1) → ℝ
  | P _ => 0
  | Q _ => rB
  | A _ _ => -b
  | B _ _ => d
  | C _ _ => f

theorem famPt_eq (U : Fin (m₁ + 1) → Fin m₁ → ℝ) (W : Fin (m₂ + 1) → Fin m₂ → ℝ)
    (rA rB a b c d e f : ℝ) (v : Vx (m₁ + 1) (m₂ + 1)) :
    famPt U W rA rB a b c d e f v =
      Fin.append (cf1 rA a c e v • U (ix1 v)) (cf2 rB b d f v • W (ix2 v)) := by
  cases v <;> simp only [famPt, cf1, cf2, ix1, ix2, zero_smul]

theorem sqd_append {k l : ℕ} (x x' : Fin k → ℝ) (y y' : Fin l → ℝ) :
    sqd (Fin.append x y) (Fin.append x' y') = sqd x x' + sqd y y' := by
  simp only [sqd, Fin.sum_univ_add, Fin.append_left, Fin.append_right]

theorem sqd_smul {m : ℕ} (α α' : ℝ) (u u' : Fin m → ℝ) :
    sqd (α • u) (α' • u') = α ^ 2 * (∑ t, u t * u t) + α' ^ 2 * (∑ t, u' t * u' t)
      - 2 * α * α' * (∑ t, u t * u' t) := by
  have h : ∀ t, (α * u t - α' * u' t) ^ 2 =
      α ^ 2 * (u t * u t) + α' ^ 2 * (u' t * u' t) - 2 * α * α' * (u t * u' t) := fun t => by ring
  simp only [sqd, Pi.smul_apply, smul_eq_mul, h, Finset.sum_add_distrib, Finset.sum_sub_distrib,
    Finset.mul_sum]

/-- Squared distances of Vₙ: the class formula (NOTE §2, family.py `class_d2`). -/
theorem sqd_famPt {U : Fin (m₁ + 1) → Fin m₁ → ℝ} {W : Fin (m₂ + 1) → Fin m₂ → ℝ}
    (hU : IsUnitSimplex U) (hW : IsUnitSimplex W) (rA rB a b c d e f : ℝ)
    (v w : Vx (m₁ + 1) (m₂ + 1)) :
    sqd (famPt U W rA rB a b c d e f v) (famPt U W rA rB a b c d e f w) =
      cf1 rA a c e v ^ 2 + cf1 rA a c e w ^ 2
        - 2 * cf1 rA a c e v * cf1 rA a c e w * (if ix1 v = ix1 w then 1 else -1 / (m₁ : ℝ))
      + (cf2 rB b d f v ^ 2 + cf2 rB b d f w ^ 2
        - 2 * cf2 rB b d f v * cf2 rB b d f w * (if ix2 v = ix2 w then 1 else -1 / (m₂ : ℝ))) := by
  rw [famPt_eq, famPt_eq, sqd_append, sqd_smul, sqd_smul, hU, hU, hU, hW, hW, hW,
    if_pos rfl, if_pos rfl, if_pos rfl, if_pos rfl]
  ring

theorem cf1_cast (rA a c e : ℚ) (v : Vx (m₁ + 1) (m₂ + 1)) :
    cf1 (rA : ℝ) (a : ℝ) (c : ℝ) (e : ℝ) v = (k1Q rA a c e (orb v) : ℝ) := by
  cases v <;> simp only [cf1, k1Q, orb, Rat.cast_neg, Rat.cast_zero]

theorem cf2_cast (rB b d f : ℚ) (v : Vx (m₁ + 1) (m₂ + 1)) :
    cf2 (rB : ℝ) (b : ℝ) (d : ℝ) (f : ℝ) v = (k2Q rB b d f (orb v) : ℝ) := by
  cases v <;> simp only [cf2, k2Q, orb, Rat.cast_neg, Rat.cast_zero]

theorem ite_cast (P : Prop) [Decidable P] (m : ℕ) :
    (if P then (1 : ℝ) else -1 / (m : ℝ)) = ((if decide P then (1 : ℚ) else -1 / (m : ℚ) : ℚ) : ℝ) := by
  by_cases h : P <;> simp [h]

/-- Squared distances of Vₙ with rational parameters are the rational class values. -/
theorem sqd_famPt_cls {U : Fin (m₁ + 1) → Fin m₁ → ℝ} {W : Fin (m₂ + 1) → Fin m₂ → ℝ}
    (hU : IsUnitSimplex U) (hW : IsUnitSimplex W) (rA rB a b c d e f : ℚ)
    (v w : Vx (m₁ + 1) (m₂ + 1)) :
    sqd (famPt U W rA rB a b c d e f v) (famPt U W rA rB a b c d e f w) =
      (clsQ m₁ m₂ rA rB a b c d e f (orb v) (orb w) (decide (ix1 v = ix1 w))
        (decide (ix2 v = ix2 w)) : ℝ) := by
  rw [sqd_famPt hU hW, cf1_cast, cf1_cast, cf2_cast, cf2_cast, ite_cast, ite_cast]
  simp only [clsQ, Rat.cast_add, Rat.cast_sub, Rat.cast_mul, Rat.cast_pow, Rat.cast_ofNat]

end Gale7
