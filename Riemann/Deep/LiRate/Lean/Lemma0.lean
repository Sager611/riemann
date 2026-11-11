import Riemann.Deep.LiRate.Lean.Cor
import Mathlib.NumberTheory.Harmonic.EulerMascheroni

/-!
# LiRate: the trend `t(n)` and Lemma 0, elementary part (proofs.md lines 7, 9–11)

`liTrend n = (n/2)(log n - κ)`, `κ = 1 + log 2π - γ_E` (proofs.md line 7). Lemma 0: if `Σ λ_n z^n`
has radius of convergence `< 1` (when RH fails, the singularities `ρ/(ρ-1)`, `Re ρ < 1/2`, of
`log ξ(1/(1-z))` lie in `|z| < 1`), then `λ_n r^n` is unbounded for some `r < 1`
(Cauchy–Hadamard); since `t(n)` is polynomially bounded, `|R(n)| ≥ e^{cn}` infinitely often.
-/

open Finset Real Filter

namespace LiRate

/-- `t(n) = (n/2)(log n - 1 - log 2π + γ_E)` (proofs.md line 7). -/
noncomputable def liTrend (n : ℕ) : ℝ :=
  n / 2 * (Real.log n - (1 + Real.log (2 * π) - eulerMascheroniConstant))

lemma liTrend_zero : liTrend 0 = 0 := by simp [liTrend]

lemma liTrend_bound (n : ℕ) : |liTrend n| ≤ 4 * (n : ℝ) ^ 2 := by
  rcases Nat.eq_zero_or_pos n with h | h
  · subst h; simp [liTrend]
  have hn1 : (1 : ℝ) ≤ n := by exact_mod_cast h
  have hlog0 : 0 ≤ Real.log n := Real.log_nonneg hn1
  have hlog1 : Real.log n ≤ n := by linarith [Real.log_le_sub_one_of_pos (by linarith : (0 : ℝ) < n)]
  have h2π0 : 0 < Real.log (2 * π) := Real.log_pos (by linarith [Real.pi_gt_three])
  have h2π1 : Real.log (2 * π) ≤ 2 * π - 1 :=
    Real.log_le_sub_one_of_pos (by positivity)
  have hγ1 := one_half_lt_eulerMascheroniConstant
  have hγ2 := eulerMascheroniConstant_lt_two_thirds
  have hπ := Real.pi_lt_d2
  have hκ : |Real.log n - (1 + Real.log (2 * π) - eulerMascheroniConstant)| ≤ n + 6 := by
    rw [abs_le]; constructor <;> linarith
  unfold liTrend
  rw [abs_mul, abs_of_nonneg (by positivity : (0 : ℝ) ≤ n / 2)]
  calc (n : ℝ) / 2 * |Real.log n - (1 + Real.log (2 * π) - eulerMascheroniConstant)|
      ≤ n / 2 * (n + 6) := mul_le_mul_of_nonneg_left hκ (by positivity)
    _ ≤ 4 * (n : ℝ) ^ 2 := by nlinarith

/-- **Lemma 0**, elementary part: `λ_n r^n` unbounded for some `0 < r < 1` and `|t(n)| ≤ C n²`
give `|λ_n - t(n)| ≥ e^{cn}` for infinitely many `n`. -/
theorem lemma0_of_unbounded (lam t : ℕ → ℝ) {r : ℝ} (hr0 : 0 < r) (hr1 : r < 1)
    (hunb : ¬ ∃ B : ℝ, ∀ n, |lam n| * r ^ n ≤ B) (C : ℝ) (ht : ∀ n, |t n| ≤ C * (n : ℝ) ^ 2) :
    Lemma0 (fun n => lam n - t n) := by
  have hfreq : ∃ᶠ n : ℕ in atTop, 1 < |lam n| * r ^ n := by
    rw [frequently_atTop]
    intro n₀
    by_contra h
    push Not at h
    apply hunb
    refine ⟨max 1 (∑ m ∈ range n₀, |lam m| * r ^ m), fun n => ?_⟩
    by_cases hn : n₀ ≤ n
    · exact (h n hn).trans (le_max_left _ _)
    · push Not at hn
      exact le_trans (single_le_sum (f := fun m => |lam m| * r ^ m)
        (fun m _ => mul_nonneg (abs_nonneg _) (pow_nonneg hr0.le _)) (mem_range.2 hn))
        (le_max_right _ _)
  have hsr : 0 < √r := Real.sqrt_pos.2 hr0
  have hsr1 : √r < 1 := by rw [Real.sqrt_lt' one_pos]; simpa using hr1
  set q : ℝ := 1 / √r with hq
  have hq1 : 1 < q := by rw [hq, lt_div_iff₀ hsr]; simpa using hsr1
  have hq0 : 0 < q := by linarith
  have hrne : r ≠ 0 := hr0.ne'
  have hrq : r * (q * q) = 1 := by
    rw [hq, div_mul_div_comm, one_mul, Real.mul_self_sqrt hr0.le]; field_simp
  refine ⟨Real.log q, Real.log_pos hq1, ?_⟩
  have hpoly : ∀ᶠ n : ℕ in atTop, (|C| + 1) * (n : ℝ) ^ 2 < q ^ n := by
    have := (tendsto_pow_const_div_const_pow_of_one_lt 2 hq1).eventually
      (gt_mem_nhds (show (0 : ℝ) < 1 / (|C| + 1) by positivity))
    filter_upwards [this] with n hn
    rw [div_lt_div_iff₀ (pow_pos hq0 n) (by positivity)] at hn
    linarith
  have h2 : ∀ᶠ n : ℕ in atTop, (2 : ℝ) ≤ q ^ n :=
    (tendsto_pow_atTop_atTop_of_one_lt hq1).eventually_ge_atTop 2
  refine (hfreq.and_eventually (hpoly.and h2)).mono fun n hn => ?_
  obtain ⟨h1, h3, h4⟩ := hn
  show Real.exp (Real.log q * n) ≤ |lam n - t n|
  rw [mul_comm, Real.exp_nat_mul, Real.exp_log hq0]
  have hqn : 0 < q ^ n * q ^ n := mul_pos (pow_pos hq0 n) (pow_pos hq0 n)
  have hrn : r ^ n * (q ^ n * q ^ n) = 1 := by
    rw [← mul_pow, ← mul_pow, hrq, one_pow]
  have hlam : q ^ n * q ^ n < |lam n| := by
    have := mul_lt_mul_of_pos_right h1 hqn
    calc q ^ n * q ^ n = 1 * (q ^ n * q ^ n) := (one_mul _).symm
      _ < |lam n| * r ^ n * (q ^ n * q ^ n) := this
      _ = |lam n| * (r ^ n * (q ^ n * q ^ n)) := by ring
      _ = |lam n| := by rw [hrn, mul_one]
  have ht' : |t n| < q ^ n := by
    have h5 := ht n
    have h6 : C * (n : ℝ) ^ 2 ≤ (|C| + 1) * (n : ℝ) ^ 2 :=
      mul_le_mul_of_nonneg_right (by linarith [le_abs_self C]) (sq_nonneg _)
    linarith
  have habs : |lam n| - |t n| ≤ |lam n - t n| := abs_sub_abs_le_abs_sub _ _
  have hqq : 2 * q ^ n ≤ q ^ n * q ^ n := mul_le_mul_of_nonneg_right h4 (pow_pos hq0 n).le
  linarith

end LiRate
