import Mathlib.Data.Int.Interval
import Mathlib.Data.Int.GCD
import Mathlib.Data.Set.Card
import Mathlib.Analysis.Real.Sqrt

/-!
# Combos NOTE §1, Lemma 1 (primitive lattice points), NOTE l.67–75

For Q(x,y) = a x² + b x y + c y² with 0 < a, |b| ≤ a ≤ c (every reduced form satisfies this;
only |D| = 4ac − b² ≥ 3a² is used) and N ≥ 1:
P_Q(N) := #{(x,y) ∈ ℤ² : gcd(x,y) = 1, Q(x,y) ≤ N} ≤ 2 + 8N|D|^{−1/2} + 4N^{1/2}|D|^{−1/4},
stated with square roots (|D|^{−1/2} = 1/√|D|, N^{1/2} = √N, |D|^{−1/4} = 1/√√|D|).
Proof as in the NOTE: rows |y| ≤ V = (4aN/|D|)^{1/2}, row y = 0 holds only (±1, 0),
each other row at most 2(N/a)^{1/2} + 1 points.
-/

open Finset

namespace RingClass

/-- Q(x, y) = a x² + b x y + c y². -/
def Qf (a b c x y : ℤ) : ℤ := a * x ^ 2 + b * x * y + c * y ^ 2

/-- The set counted by P_Q(N) (NOTE l.69). -/
def PQset (a b c : ℤ) (N : ℝ) : Set (ℤ × ℤ) :=
  {p | Int.gcd p.1 p.2 = 1 ∧ (Qf a b c p.1 p.2 : ℝ) ≤ N}

/-- P_Q(N) = #{(x,y) ∈ ℤ² : gcd(x,y) = 1, Q(x,y) ≤ N} (NOTE l.69); the set is finite
(PQset_finite), so ncard is its cardinality. -/
noncomputable def PQ (a b c : ℤ) (N : ℝ) : ℕ := (PQset a b c N).ncard

/-- At most L + 1 integers lie in [l, l + L]. -/
lemma card_Icc_le (l L : ℝ) (hL : 0 ≤ L) :
    ((Icc ⌈l⌉ ⌊l + L⌋).card : ℝ) ≤ L + 1 := by
  rw [Int.card_Icc]
  have h1 : ((⌊l + L⌋ : ℤ) : ℝ) ≤ l + L := Int.floor_le _
  have h2 : l ≤ ((⌈l⌉ : ℤ) : ℝ) := Int.le_ceil _
  rcases le_or_gt 0 (⌊l + L⌋ + 1 - ⌈l⌉) with h | h
  · have e : (((⌊l + L⌋ + 1 - ⌈l⌉).toNat : ℕ) : ℝ) = ((⌊l + L⌋ + 1 - ⌈l⌉ : ℤ) : ℝ) := by
      rw [← Int.cast_natCast, Int.toNat_of_nonneg h]
    rw [e]; push_cast; linarith
  · rw [Int.toNat_eq_zero.mpr h.le]; simp only [Nat.cast_zero]; linarith

/-- the candidate rows y ≠ 0, |y| ≤ V -/
noncomputable def rows (V : ℝ) : Finset ℤ := (Icc ⌈-V⌉ ⌊-V + 2 * V⌋).erase 0

/-- left end of the x-interval of row y (p plays (aN)^{1/2}) -/
noncomputable def lft (a b : ℤ) (p : ℝ) (y : ℤ) : ℝ := (-(b : ℝ) * y - 2 * p) / (2 * a)

/-- the points of row y -/
noncomputable def row (a b : ℤ) (p L : ℝ) (y : ℤ) : Finset (ℤ × ℤ) :=
  (Icc ⌈lft a b p y⌉ ⌊lft a b p y + L⌋).image (fun x => (x, y))

lemma card_rows_le (V : ℝ) (hV : 0 ≤ V) : ((rows V).card : ℝ) ≤ 2 * V := by
  have h0 : (0 : ℤ) ∈ Icc ⌈-V⌉ ⌊-V + 2 * V⌋ := by
    rw [Finset.mem_Icc, Int.ceil_le, Int.le_floor]; constructor <;> push_cast <;> linarith
  have hc := card_Icc_le (-V) (2 * V) (by linarith)
  have hpos : 1 ≤ (Icc ⌈-V⌉ ⌊-V + 2 * V⌋).card := Finset.card_pos.mpr ⟨0, h0⟩
  rw [rows, Finset.card_erase_of_mem h0, Nat.cast_sub hpos]
  push_cast; linarith

lemma disc_facts (a b c : ℝ) (ha : 0 < a) (hb : |b| ≤ a) (hac : a ≤ c) :
    3 * a ^ 2 ≤ 4 * a * c - b ^ 2 := by
  have hb2 : b ^ 2 ≤ a ^ 2 := sq_le_sq' (abs_le.mp hb).1 (abs_le.mp hb).2
  have h4 : a ^ 2 ≤ a * c := by rw [sq]; exact mul_le_mul_of_nonneg_left hac ha.le
  linarith

end RingClass
