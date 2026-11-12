import Mathlib.Basic.Real.Basic
import Mathlib.Algebra.BigOperators.Group.Finset.Basic
import Mathlib.Algebra.Order.BigOperators.Group.Finset
import Mathlib.Algebra.Order.BigOperators.GroupWithZero.Finset
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.Positivity
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.NormNum
import Mathlib.Tactic.Ring

/-!
# Abstract logarithm (light file: no Mathlib analysis imported)

Under the 290 MB memory cap no file importing Mathlib's Real.log (Mathlib.Analysis.SpecialFunctions.Log.Basic,
peak footprint ≥ 326 MB alone) can be compiled. The light files therefore work with a function L : ℝ → ℝ
satisfying LogLike L: exactly the properties of Real.log used, each a Mathlib theorem
(Real.log_mul, Real.log_le_log, Real.one_sub_inv_le_log_of_pos, and the two series bounds proved from
Real.hasSum_log_one_add in LogSeries.lean). RealLog.lean proves LogLike Real.log.
-/

namespace OddRobin

open Finset

/-- Partial sum of the series log (1 + a) = sum 2 y^(2k+1)/(2k+1), y = a/(a+2). -/
noncomputable def lsum (a : ℝ) (m : ℕ) : ℝ :=
  ∑ k ∈ range m, (2 : ℝ) * (1 / (2 * k + 1)) * (a / (a + 2)) ^ (2 * k + 1)

/-- Geometric tail bound of that series after m terms. -/
noncomputable def ltail (a : ℝ) (m : ℕ) : ℝ :=
  2 * (a / (a + 2)) ^ (2 * m + 1) / ((2 * m + 1) * (1 - (a / (a + 2)) ^ 2))

/-- The properties of Real.log used by the proof (all are Mathlib theorems for Real.log). -/
structure LogLike (L : ℝ → ℝ) : Prop where
  map_mul : ∀ x y : ℝ, 0 < x → 0 < y → L (x * y) = L x + L y
  mono : ∀ x y : ℝ, 0 < x → x ≤ y → L x ≤ L y
  one_sub_inv_le : ∀ x : ℝ, 0 < x → 1 - x⁻¹ ≤ L x
  series_lo : ∀ a : ℝ, 0 ≤ a → ∀ m : ℕ, lsum a m ≤ L (1 + a)
  series_hi : ∀ a : ℝ, 0 ≤ a → ∀ m : ℕ, L (1 + a) ≤ lsum a m + ltail a m

variable {L : ℝ → ℝ}

namespace LogLike

lemma map_one (hL : LogLike L) : L 1 = 0 := by
  have h := hL.map_mul 1 1 one_pos one_pos
  rw [one_mul] at h
  linarith

lemma map_div (hL : LogLike L) {x y : ℝ} (hx : 0 < x) (hy : 0 < y) : L (x / y) = L x - L y := by
  have h := hL.map_mul (x / y) y (div_pos hx hy) hy
  rw [div_mul_cancel₀ x hy.ne'] at h
  linarith

lemma map_pow (hL : LogLike L) {x : ℝ} (hx : 0 < x) (k : ℕ) : L (x ^ k) = k * L x := by
  induction k with
  | zero => simp [hL.map_one]
  | succ k ih =>
    rw [pow_succ, hL.map_mul _ _ (pow_pos hx k) hx, ih]
    push_cast
    ring

lemma map_prod (hL : LogLike L) {ι : Type*} (s : Finset ι) (f : ι → ℝ) (hf : ∀ i ∈ s, 0 < f i) :
    L (∏ i ∈ s, f i) = ∑ i ∈ s, L (f i) := by
  classical
  induction s using Finset.induction_on with
  | empty => simp [hL.map_one]
  | insert a s ha ih =>
    rw [Finset.prod_insert ha, Finset.sum_insert ha,
      hL.map_mul _ _ (hf a (Finset.mem_insert_self a s))
        (Finset.prod_pos (fun i hi => hf i (Finset.mem_insert_of_mem hi))),
      ih (fun i hi => hf i (Finset.mem_insert_of_mem hi))]

lemma le_sub_one (hL : LogLike L) {x : ℝ} (hx : 0 < x) : L x ≤ x - 1 := by
  have h := hL.one_sub_inv_le x⁻¹ (inv_pos.2 hx)
  have h2 := hL.map_mul x x⁻¹ hx (inv_pos.2 hx)
  rw [mul_inv_cancel₀ hx.ne', hL.map_one] at h2
  rw [inv_inv] at h
  linarith

lemma pos_of_one_lt (hL : LogLike L) {x : ℝ} (hx : 1 < x) : 0 < L x := by
  have h := hL.one_sub_inv_le x (by linarith)
  have : x⁻¹ < 1 := inv_lt_one_of_one_lt₀ hx
  linarith

lemma strictMono (hL : LogLike L) {x y : ℝ} (hx : 0 < x) (hxy : x < y) : L x < L y := by
  have hy : 0 < y := by linarith
  have h := hL.map_mul x (y / x) hx (div_pos hy hx)
  rw [mul_div_cancel₀ y hx.ne'] at h
  have : 0 < L (y / x) := hL.pos_of_one_lt (by rw [lt_div_iff₀ hx]; linarith)
  linarith

lemma lt_iff (hL : LogLike L) {x y : ℝ} (hx : 0 < x) (hy : 0 < y) : L x < L y ↔ x < y := by
  constructor
  · intro h
    by_contra hc
    have := hL.mono y x hy (not_lt.1 hc)
    linarith
  · exact hL.strictMono hx

lemma two_pow_mul (hL : LogLike L) (k : ℕ) {r : ℝ} (hr : 0 < r) :
    L (2 ^ k * r) = k * L 2 + L r := by
  rw [hL.map_mul _ _ (by positivity) hr, hL.map_pow (by norm_num) k]

/-- The case m = 1 of the series bounds: 2y ≤ L (1 + a). -/
lemma two_y_le (hL : LogLike L) {a : ℝ} (ha : 0 ≤ a) : 2 * (a / (a + 2)) ≤ L (1 + a) := by
  have h := hL.series_lo a ha 1
  simpa [lsum] using h

/-- The case m = 1 of the series bounds: L (1 + a) ≤ 2y + 2y^3/(3(1 − y^2)). -/
lemma le_two_y_add (hL : LogLike L) {a : ℝ} (ha : 0 ≤ a) :
    L (1 + a) ≤ 2 * (a / (a + 2)) + 2 * (a / (a + 2)) ^ 3 / (3 * (1 - (a / (a + 2)) ^ 2)) := by
  have h := hL.series_hi a ha 1
  norm_num [lsum, ltail] at h ⊢
  linarith

/-- Upper bound for L q, q = 1 + a (used with explicit rationals, closed by norm_num). -/
lemma le_of (hL : LogLike L) {q a U : ℝ} (hq : q = 1 + a) (ha : 0 ≤ a) (m : ℕ)
    (h : lsum a m + ltail a m ≤ U) : L q ≤ U :=
  hq ▸ (hL.series_hi a ha m).trans h

/-- Lower bound for L q, q = 1 + a. -/
lemma ge_of (hL : LogLike L) {q a V : ℝ} (hq : q = 1 + a) (ha : 0 ≤ a) (m : ℕ) (h : V ≤ lsum a m) :
    V ≤ L q :=
  hq ▸ h.trans (hL.series_lo a ha m)

end LogLike

end OddRobin

