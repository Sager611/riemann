import Mathlib.Algebra.Order.Field.Basic
import Mathlib.Data.Nat.Prime.Nth
import Mathlib.Data.Nat.Fib.Basic
import Mathlib.Order.Monotone.Basic
import Mathlib.Tactic.NormNum
import Mathlib.Tactic.Ring

/-!
# Sun's conjectures on `ⁿ√·` of prime and Lucas sequences: 1211.1588 Conj. 4.30(iii); 1208.2683 Conj. 3.1

The sources are about `ℝ` and `ⁿ√x = x ^ (1/n)`. `Mathlib.Analysis.SpecialFunctions.Pow.Real` does not
load under the 290 MB cap (killed at 318 MB), so the root is a parameter: any `rt : ℕ → K → K` on a
linearly ordered field `K` with `rt n x > 0` and `(rt n x) ^ n = x` for `n ≥ 1`, `x > 0` (`IsRootFn`).
Positive `n`-th roots are unique in an ordered field, so for `K = ℝ` every such `rt` agrees with
`x ↦ x ^ (1/n : ℝ)` on `x > 0`, and the refutations below hold for the real conjectures.
-/

namespace SunLean

/-- `rt n x` is a positive `n`-th root of every `x > 0`, for every `n ≥ 1`. -/
def IsRootFn (K : Type) [Field K] [LinearOrder K] (rt : ℕ → K → K) : Prop :=
  ∀ n : ℕ, 0 < n → ∀ x : K, 0 < x → 0 < rt n x ∧ rt n x ^ n = x

/-- `p_n`, the `n`-th prime (`p_1 = 2`). -/
noncomputable def pr (n : ℕ) : ℕ := Nat.nth Nat.Prime (n - 1)

/-- `p^{(m)}_n` of arXiv 1211.1588 (TeX l. 1339-1340): `p^{(1)}_n = p_n`, `p^{(m+1)}_n = p^{(m)}_{p_n}`
(here with the auxiliary `p^{(0)}_n = n`, which gives exactly these two equations). -/
noncomputable def iterPrime : ℕ → ℕ → ℕ
  | 0, n => n
  | m + 1, n => iterPrime m (pr n)

theorem iterPrime_one (n : ℕ) : iterPrime 1 n = pr n := rfl
theorem iterPrime_succ (m n : ℕ) : iterPrime (m + 1) n = iterPrime m (pr n) := rfl

/-- Z.-W. Sun, arXiv 1211.1588 v29, Conjecture 4.30(iii) (TeX l. 1349-1352): for all `m, n ≥ 1`,
`ⁿ⁺¹√(p^{(m+1)}_{n+1}) / ⁿ√(p^{(m+1)}_n) < ⁿ⁺¹√(p^{(m)}_{n+1}) / ⁿ√(p^{(m)}_n) < 1`. -/
def Conj_4_30_iii (K : Type) [Field K] [LinearOrder K] (rt : ℕ → K → K) : Prop :=
  ∀ m n : ℕ, 0 < m → 0 < n →
    rt (n + 1) (iterPrime (m + 1) (n + 1) : K) / rt n (iterPrime (m + 1) n : K) <
        rt (n + 1) (iterPrime m (n + 1) : K) / rt n (iterPrime m n : K) ∧
      rt (n + 1) (iterPrime m (n + 1) : K) / rt n (iterPrime m n : K) < 1

/-- `m = 1`, `n = 2`: `p^{(2)}_3 = 11`, `p^{(2)}_2 = 5`, `p^{(1)}_3 = 5`, `p^{(1)}_2 = 3`, and
`∛11 / √5 > ∛5 / √3` since `11² · 3³ = 3267 > 5⁵ = 3125`. -/
theorem not_conj_4_30_iii (K : Type) [Field K] [LinearOrder K] [IsStrictOrderedRing K]
    (rt : ℕ → K → K) (hrt : IsRootFn K rt) : ¬ Conj_4_30_iii K rt := by
  intro h
  have h1 := (h 1 2 one_pos two_pos).1
  have e1 : iterPrime (1 + 1) (2 + 1) = 11 := by
    simp [iterPrime, pr, Nat.nth_prime_two_eq_five, Nat.nth_prime_four_eq_eleven]
  have e2 : iterPrime (1 + 1) 2 = 5 := by
    simp [iterPrime, pr, Nat.nth_prime_one_eq_three, Nat.nth_prime_two_eq_five]
  have e3 : iterPrime 1 (2 + 1) = 5 := by simp [iterPrime, pr, Nat.nth_prime_two_eq_five]
  have e4 : iterPrime 1 2 = 3 := by simp [iterPrime, pr, Nat.nth_prime_one_eq_three]
  rw [e1, e2, e3, e4] at h1
  norm_num at h1
  obtain ⟨ha, ha3⟩ := hrt 3 (by norm_num) 11 (by norm_num)
  obtain ⟨hb, hb2⟩ := hrt 2 (by norm_num) 5 (by norm_num)
  obtain ⟨hc, hc3⟩ := hrt 3 (by norm_num) 5 (by norm_num)
  obtain ⟨hd, hd2⟩ := hrt 2 (by norm_num) 3 (by norm_num)
  have key : rt 3 11 * rt 2 3 < rt 3 5 * rt 2 5 := (div_lt_div_iff₀ hb hd).1 h1
  have k6 := pow_lt_pow_left₀ key (mul_pos ha hd).le (by norm_num : (6 : ℕ) ≠ 0)
  have L : (rt 3 11 * rt 2 3) ^ 6 = (rt 3 11 ^ 3) ^ 2 * (rt 2 3 ^ 2) ^ 3 := by ring
  have R : (rt 3 5 * rt 2 5) ^ 6 = (rt 3 5 ^ 3) ^ 2 * (rt 2 5 ^ 2) ^ 3 := by ring
  rw [L, R, ha3, hd2, hc3, hb2] at k6
  norm_num at k6

/-- `u_0 = 0`, `u_1 = 1`, `u_(n+1) = A u_n - B u_(n-1)`. -/
def lucasU (A B : ℤ) : ℕ → ℤ
  | 0 => 0
  | 1 => 1
  | n + 2 => A * lucasU A B (n + 1) - B * lucasU A B n

/-- Z.-W. Sun, arXiv 1208.2683 v9, Conjecture 3.1 (TeX l. 462-470), all three parts:
`(ⁿ√F_n)_{n ≥ 2}` strictly increasing; `(ⁿ⁺¹√F_(n+1) / ⁿ√F_n)_{n ≥ 4}` strictly decreasing; and for all
integers `A > 1`, `B ≠ 0` with `A² > 4B` and (`A > 2` or `B ≥ -9`), `(ⁿ⁺¹√u_(n+1) / ⁿ√u_n)_{n ≥ 4}` is
strictly decreasing with limit 1. -/
def Conj_3_1 (K : Type) [Field K] [LinearOrder K] (rt : ℕ → K → K) : Prop :=
  StrictMonoOn (fun n : ℕ => rt n (Nat.fib n : K)) {n | 2 ≤ n} ∧
  StrictAntiOn (fun n : ℕ => rt (n + 1) (Nat.fib (n + 1) : K) / rt n (Nat.fib n : K)) {n | 4 ≤ n} ∧
  ∀ A B : ℤ, 1 < A → B ≠ 0 → 4 * B < A ^ 2 → (2 < A ∨ -9 ≤ B) →
    StrictAntiOn (fun n : ℕ => rt (n + 1) (lucasU A B (n + 1) : K) / rt n (lucasU A B n : K))
        {n | 4 ≤ n} ∧
      ∀ ε : K, 0 < ε → ∃ N : ℕ, ∀ n ≥ N,
        |rt (n + 1) (lucasU A B (n + 1) : K) / rt n (lucasU A B n : K) - 1| < ε

/-- `(A, B) = (3, -26)` (`A > 2`, `A² = 9 > -104`): `u_5, u_6, u_7 = 1459, 9135, 65339` and
`u_6^70 < u_5^42 u_7^30`, so `⁶√u_6 / ⁵√u_5 < ⁷√u_7 / ⁶√u_6`: not decreasing at `n = 5`. -/
theorem not_conj_3_1 (K : Type) [Field K] [LinearOrder K] [IsStrictOrderedRing K]
    (rt : ℕ → K → K) (hrt : IsRootFn K rt) : ¬ Conj_3_1 K rt := by
  intro h
  have hd := (h.2.2 3 (-26) (by norm_num) (by norm_num) (by norm_num) (Or.inl (by norm_num))).1
  have h56 := hd (show 5 ∈ {n : ℕ | 4 ≤ n} by simp) (show 6 ∈ {n : ℕ | 4 ≤ n} by simp)
    (by norm_num : (5 : ℕ) < 6)
  have u5 : lucasU 3 (-26) 5 = 1459 := by decide
  have u6 : lucasU 3 (-26) 6 = 9135 := by decide
  have u7 : lucasU 3 (-26) 7 = 65339 := by decide
  simp only [Nat.reduceAdd, u5, u6, u7] at h56
  norm_num at h56
  obtain ⟨ha, ha5⟩ := hrt 5 (by norm_num) 1459 (by norm_num)
  obtain ⟨hb, hb6⟩ := hrt 6 (by norm_num) 9135 (by norm_num)
  obtain ⟨hc, hc7⟩ := hrt 7 (by norm_num) 65339 (by norm_num)
  have key : rt 7 65339 * rt 5 1459 < rt 6 9135 * rt 6 9135 := (div_lt_div_iff₀ hb ha).1 h56
  have k := pow_lt_pow_left₀ key (mul_pos hc ha).le (by norm_num : (210 : ℕ) ≠ 0)
  have L : (rt 7 65339 * rt 5 1459) ^ 210 = (rt 7 65339 ^ 7) ^ 30 * (rt 5 1459 ^ 5) ^ 42 := by ring
  have R : (rt 6 9135 * rt 6 9135) ^ 210 = (rt 6 9135 ^ 6) ^ 70 := by ring
  rw [L, R, hc7, ha5, hb6] at k
  norm_num at k

end SunLean

#print axioms SunLean.not_conj_4_30_iii
#print axioms SunLean.not_conj_3_1
