import Mathlib.Algebra.Order.Field.Rat
import Mathlib.Data.Nat.Choose.Basic
import Mathlib.Algebra.BigOperators.Group.Finset.Basic
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.Ring
import Mathlib.Tactic.NormNum
import Mathlib.Tactic.Positivity

/-!
# Fasi–Higham conjecture, m = 3: definitions

Fasi–Higham (SIMAX 40 (2019), Conj. 3.2): q_{k,m}(z) = Σ_{j=0}^{m} (−1)^j C(m,j) (k+m−j)!/(k+m)! z^j,
the denominator of the [k/m] Padé approximant of e^z; conjecture: for odd m every Taylor coefficient
of 1/q_{k,m} at 0 is positive. Here m = 3. qcoef k m j is the coefficient of z^j of q_{k,m};
c k n is defined by the recurrence q · Σ c_n z^n = 1, and conv (Conv.lean) proves that identity.
-/

namespace FasiLean

open Finset

/-- Coefficient of z^j in q_{k,m}: (−1)^j C(m,j) (k+m−j)!/(k+m)! for j ≤ m, and 0 for j > m. -/
def qcoef (k m j : ℕ) : ℚ :=
  if j ≤ m then (-1) ^ j * (m.choose j : ℚ) * ((k + m - j).factorial : ℚ) / ((k + m).factorial : ℚ)
  else 0

/-- Shifted coefficients: cs k 0 = c_{−2} = 0, cs k 1 = c_{−1} = 0, cs k 2 = c_0 = 1 and, for n ≥ 0,
c_{n+1} = 3/(k+3) c_n − 3/((k+2)(k+3)) c_{n−1} + 1/((k+1)(k+2)(k+3)) c_{n−2}. -/
def cs (k : ℕ) : ℕ → ℚ
  | 0 => 0
  | 1 => 0
  | 2 => 1
  | n + 3 => 3 / (k + 3) * cs k (n + 2) - 3 / ((k + 2) * (k + 3)) * cs k (n + 1)
      + 1 / ((k + 1) * (k + 2) * (k + 3)) * cs k n

/-- c_n(k), the n-th Taylor coefficient of 1/q_{k,3} at 0 (see conv and PS.lean). -/
def c (k n : ℕ) : ℚ := cs k (n + 2)

theorem qcoef_zero (k : ℕ) : qcoef k 3 0 = 1 := by
  have h : ((k + 3).factorial : ℚ) ≠ 0 := by positivity
  simp only [qcoef, ite_eq_left (show 0 ≤ 3 by norm_num), Nat.sub_zero, pow_zero, Nat.choose_zero_right]
  field_simp
  norm_num

theorem qcoef_one (k : ℕ) : qcoef k 3 1 = -3 / (k + 3) := by
  have h : ((k + 2).factorial : ℚ) ≠ 0 := by positivity
  simp only [qcoef, ite_eq_left (show 1 ≤ 3 by norm_num), show k + 3 - 1 = k + 2 by omega,
    show k + 3 = (k + 2) + 1 by omega, Nat.factorial_succ]
  push_cast
  rw [show (Nat.choose 3 1 : ℚ) = 3 by norm_num [Nat.choose]]
  field_simp
  ring

theorem qcoef_two (k : ℕ) : qcoef k 3 2 = 3 / ((k + 2) * (k + 3)) := by
  have h : ((k + 1).factorial : ℚ) ≠ 0 := by positivity
  simp only [qcoef, ite_eq_left (show 2 ≤ 3 by norm_num), show k + 3 - 2 = k + 1 by omega,
    show k + 3 = (k + 1) + 1 + 1 by omega, Nat.factorial_succ]
  push_cast
  rw [show (Nat.choose 3 2 : ℚ) = 3 by norm_num [Nat.choose]]
  field_simp
  ring

theorem qcoef_three (k : ℕ) : qcoef k 3 3 = -1 / ((k + 1) * (k + 2) * (k + 3)) := by
  have h : ((k).factorial : ℚ) ≠ 0 := by positivity
  simp only [qcoef, ite_eq_left (show 3 ≤ 3 by norm_num), show k + 3 - 3 = k by omega,
    show k + 3 = k + 1 + 1 + 1 by omega, Nat.factorial_succ]
  push_cast
  rw [show (Nat.choose 3 3 : ℚ) = 1 by norm_num [Nat.choose]]
  field_simp
  ring

theorem qcoef_big (k j : ℕ) : qcoef k 3 (j + 4) = 0 := by
  simp only [qcoef, ite_eq_right (show ¬ (j + 4 ≤ 3) by omega)]

theorem cs_rec (k n : ℕ) : cs k (n + 3) = 3 / (k + 3) * cs k (n + 2)
    - 3 / ((k + 2) * (k + 3)) * cs k (n + 1) + 1 / ((k + 1) * (k + 2) * (k + 3)) * cs k n := rfl

theorem conv_succ3 (k n : ℕ) :
    ∑ j ∈ range (n + 4), qcoef k 3 j * c k (n + 3 - j) = 0 := by
  have hz : ∑ i ∈ range n, qcoef k 3 (i + 1 + 1 + 1 + 1) * c k (n + 3 - (i + 1 + 1 + 1 + 1)) = 0 := by
    apply Finset.sum_eq_zero
    intro i _
    rw [show i + 1 + 1 + 1 + 1 = i + 4 by omega, qcoef_big, zero_mul]
  rw [Finset.sum_range_succ', Finset.sum_range_succ', Finset.sum_range_succ', Finset.sum_range_succ', hz]
  rw [show 0 + 1 + 1 + 1 = 3 by rfl, show 0 + 1 + 1 = 2 by rfl, show (0 : ℕ) + 1 = 1 by rfl,
    show n + 3 - 3 = n by omega, show n + 3 - 2 = n + 1 by omega, show n + 3 - 1 = n + 2 by omega,
    Nat.sub_zero, qcoef_zero, qcoef_one, qcoef_two, qcoef_three]
  simp only [c]
  rw [show n + 3 + 2 = (n + 2) + 3 by omega, cs_rec k (n + 2)]
  rw [show n + 2 + 2 = n + 4 by omega, show n + 1 + 2 = n + 3 by omega]
  ring

/-- q_{k,3} · Σ_n c_n z^n = 1, coefficientwise: Σ_{j ≤ n} q_j c_{n−j} = 1 if n = 0, else 0. -/
theorem conv (k n : ℕ) :
    ∑ j ∈ range (n + 1), qcoef k 3 j * c k (n - j) = if n = 0 then 1 else 0 := by
  match n with
  | 0 => simp [qcoef_zero, c, cs]
  | 1 =>
    simp only [Finset.sum_range_succ, Finset.sum_range_zero, qcoef_zero, qcoef_one, c, cs]
    norm_num [cs]
    ring
  | 2 =>
    simp only [Finset.sum_range_succ, Finset.sum_range_zero, qcoef_zero, qcoef_one, qcoef_two, c, cs]
    norm_num [cs]
    ring
  | n + 3 =>
    rw [ite_eq_right (by omega)]
    exact conv_succ3 k n

end FasiLean
