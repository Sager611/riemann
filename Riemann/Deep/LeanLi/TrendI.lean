module

public import Riemann.Deep.LeanLi.A0
public import Riemann.Deep.LeanLi.PhiFacts
public import Riemann.Deep.LeanLi.IFine
public import Mathlib.NumberTheory.Harmonic.Bounds

/-!
# LeanLi: the trend transform of the model trend `I(n)`: `τ_I(a) = ∫₀^π B(a,v) a₀(v) dv`
-/

@[expose] public section

open Real Finset MeasureTheory Set

namespace LeanLi

/-- The integrand of `Imod n` is integrable on `[0, π]`. -/
lemma intervalIntegrable_Imod_integrand (n : ℕ) :
    IntervalIntegrable (fun v => 2 * (1 - cos (n * v)) * a0 v) volume 0 π := by
  have hF : IntervalIntegrable (fun v => (1 / (2 * π)) *
      ((n + 2 * ∑ j ∈ Finset.Ico 1 n, ((n : ℝ) - j) * cos (j * v)) * ell v)) volume 0 π := by
    have := (intervalIntegrable_ell 0 π).continuousOn_mul
      (g := fun v => (n + 2 * ∑ j ∈ Finset.Ico 1 n, ((n : ℝ) - j) * cos (j * v)))
      (Continuous.continuousOn (by fun_prop))
    exact this.const_mul _
  refine hF.congr_uIoo fun v hv => ?_
  rw [uIoo_of_le pi_pos.le] at hv
  exact (integrand_eq n (sin_half_pos hv.1 hv.2.le).ne').symm

lemma B_eq_sum (N : ℕ) (w : ℕ → ℝ) (a v : ℝ) :
    LR.B N w a v = 2 * ∑ n ∈ Icc 1 N, w n * cos (n * a) * (2 * (1 - cos (n * v))) := by
  rw [← LR.li_pair]
  unfold LR.ip
  simp

lemma B_mul_a0_eq (N : ℕ) (w : ℕ → ℝ) (a v : ℝ) :
    LR.B N w a v * a0 v = 2 * ∑ n ∈ Icc 1 N, w n * cos (n * a) * (2 * (1 - cos (n * v)) * a0 v) := by
  rw [B_eq_sum, mul_assoc, Finset.sum_mul]
  congr 1
  refine Finset.sum_congr rfl fun n _ => ?_
  ring

lemma intervalIntegrable_B_a0 (N : ℕ) (w : ℕ → ℝ) (a : ℝ) :
    IntervalIntegrable (fun v => LR.B N w a v * a0 v) volume 0 π := by
  have h : IntervalIntegrable (fun v => 2 * ∑ n ∈ Icc 1 N, w n * cos (n * a) *
      (2 * (1 - cos (n * v)) * a0 v)) volume 0 π := by
    have := IntervalIntegrable.sum (Icc 1 N) fun n _ =>
      (intervalIntegrable_Imod_integrand n).const_mul (w n * cos (n * a))
    refine (this.const_mul 2).congr_uIoo fun v _ => ?_
    simp [Finset.sum_apply]
  exact h.congr_uIoo fun v _ => (B_mul_a0_eq N w a v).symm

/-- `τ_I(a) = ∫₀^π B(a, v) a₀(v) dv`. -/
theorem trend_Iclosed_eq (N : ℕ) (w : ℕ → ℝ) (a : ℝ) :
    LR.trend N w Iclosed a = ∫ v in (0 : ℝ)..π, LR.B N w a v * a0 v := by
  have hI0 : Iclosed 0 = 0 := by simp [Iclosed]
  unfold LR.trend LR.ip
  rw [hI0, mul_zero, zero_add]
  simp_rw [B_mul_a0_eq]
  rw [intervalIntegral.integral_const_mul, intervalIntegral.integral_finsetSum
    (fun n _ => (intervalIntegrable_Imod_integrand n).const_mul _)]
  congr 1
  refine Finset.sum_congr rfl fun n _ => ?_
  rw [intervalIntegral.integral_const_mul, ← Imod_eq_closed]
  rfl

/-- `|τ_t(a) - τ_I(a)| ≤ 8 (1 + log N)` (from `|t(n) - I(n)| ≤ 4/n`). -/
theorem trend_sub_le (N : ℕ) (w : ℕ → ℝ) (hw0 : ∀ n, 0 ≤ w n) (hw1 : ∀ n, w n ≤ 1) (a : ℝ) :
    |LR.trend N w trendT a - LR.trend N w Iclosed a| ≤ 8 * (1 + Real.log N) := by
  have hI0 : Iclosed 0 = 0 := by simp [Iclosed]
  have ht0 : trendT 0 = 0 := by simp [trendT]
  unfold LR.trend LR.ip
  rw [hI0, ht0]
  simp only [mul_zero, zero_add]
  rw [← mul_sub, ← Finset.sum_sub_distrib]
  have hterm : ∀ n ∈ Finset.Icc 1 N, |w n * cos (n * a) * trendT n - w n * cos (n * a) * Iclosed n| ≤ 4 / n := by
    intro n hn
    rw [← mul_sub, abs_mul, abs_mul, abs_sub_comm]
    have h1 := Iclosed_sub_trend_le_div n (Finset.mem_Icc.1 hn).1
    have h2 : |w n| ≤ 1 := by rw [abs_of_nonneg (hw0 n)]; exact hw1 n
    have h3 := abs_cos_le_one (n * a)
    calc |w n| * |cos (n * a)| * |Iclosed n - trendT n| ≤ 1 * 1 * (4 / n) := by
          gcongr
      _ = 4 / n := by ring
  have hsum : ∑ n ∈ Finset.Icc 1 N, (4 : ℝ) / n ≤ 4 * (1 + Real.log N) := by
    have := harmonic_le_one_add_log N
    have e : ∑ n ∈ Finset.Icc 1 N, (1 : ℝ) / n = harmonic N := by
      rw [harmonic_eq_sum_Icc]; push_cast; simp [one_div]
    simp_rw [div_eq_mul_one_div (4 : ℝ)]
    rw [← Finset.mul_sum, e]
    linarith
  rw [abs_mul, abs_of_pos (by norm_num : (0 : ℝ) < 2)]
  calc 2 * |∑ n ∈ Icc 1 N, (w n * cos (n * a) * trendT n - w n * cos (n * a) * Iclosed n)|
      ≤ 2 * ∑ n ∈ Icc 1 N, |w n * cos (n * a) * trendT n - w n * cos (n * a) * Iclosed n| := by
        gcongr; exact Finset.abs_sum_le_sum_abs _ _
    _ ≤ 2 * ∑ n ∈ Icc 1 N, (4 : ℝ) / n := by gcongr with n hn; exact hterm n hn
    _ ≤ 2 * (4 * (1 + Real.log N)) := by gcongr
    _ = 8 * (1 + Real.log N) := by ring

end LeanLi
