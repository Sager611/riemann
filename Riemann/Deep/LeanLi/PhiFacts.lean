module

public import Riemann.Deep.LeanLi.Weights
public import Mathlib.Analysis.SpecialFunctions.Integrals.Basic

/-!
# LeanLi: facts on `Φ(x) = w₀ + 2 Σ_{n≤N} w_n cos(nx)` used by the trend estimate
-/

@[expose] public section

open Real Finset MeasureTheory

namespace LeanLi

lemma integral_cos_nat_mul (n : ℕ) (hn : 1 ≤ n) : ∫ x in (-π)..π, cos (n * x) = 0 := by
  have hn0 : (n : ℝ) ≠ 0 := by exact_mod_cast (by omega : n ≠ 0)
  rw [intervalIntegral.integral_comp_mul_left (fun x => cos x) hn0, integral_cos]
  rw [show (n : ℝ) * -π = -(n * π) by ring, sin_neg, sin_nat_mul_pi]
  simp

lemma intervalIntegrable_phi (N : ℕ) (w : ℕ → ℝ) (a b : ℝ) :
    IntervalIntegrable (LR.Phi N w) volume a b := by
  apply Continuous.intervalIntegrable
  unfold LR.Phi
  fun_prop

lemma integral_phi (N : ℕ) (w : ℕ → ℝ) : ∫ x in (-π)..π, LR.Phi N w x = 2 * π * w 0 := by
  unfold LR.Phi
  have hint : ∀ n ∈ Icc 1 N, IntervalIntegrable (fun x : ℝ => w n * cos (n * x)) volume (-π) π :=
    fun n _ => (Continuous.intervalIntegrable (by fun_prop) _ _)
  have hsum : IntervalIntegrable (fun x : ℝ => 2 * ∑ n ∈ Icc 1 N, w n * cos (n * x)) volume (-π) π := by
    apply Continuous.intervalIntegrable; fun_prop
  rw [intervalIntegral.integral_add (by simp) hsum,
    intervalIntegral.integral_const, intervalIntegral.integral_const_mul,
    intervalIntegral.integral_finsetSum hint]
  have : ∀ n ∈ Icc 1 N, ∫ x in (-π)..π, w n * cos (n * x) = 0 := by
    intro n hn
    rw [intervalIntegral.integral_const_mul, integral_cos_nat_mul n (Finset.mem_Icc.1 hn).1, mul_zero]
  rw [Finset.sum_eq_zero this]
  simp only [sub_neg_eq_add, smul_eq_mul, mul_zero, add_zero]
  ring

lemma phi_neg (N : ℕ) (w : ℕ → ℝ) (x : ℝ) : LR.Phi N w (-x) = LR.Phi N w x := by
  unfold LR.Phi; simp [mul_neg, cos_neg]

lemma phi_two_pi_sub (N : ℕ) (w : ℕ → ℝ) (x : ℝ) : LR.Phi N w (2 * π - x) = LR.Phi N w x := by
  unfold LR.Phi
  congr 2
  refine Finset.sum_congr rfl fun n _ => ?_
  rw [show (n : ℝ) * (2 * π - x) = -(n * x) + (n : ℕ) * (2 * π) by ring, cos_add_nat_mul_two_pi, cos_neg]

/-- Off a `1/(2√N)`-neighbourhood of `0`, `Φ ≤ C 2³⁰/N¹⁴` (from the decay with `2K = 30`). -/
lemma phi_far {Φ : ℝ → ℝ} {C : ℝ} {N : ℕ} (hN : 1 ≤ N)
    (hdec : ∀ x, |x| ≤ π → Φ x * (1 + N * |x|) ^ 30 ≤ C * N) (hΦ : ∀ x, 0 ≤ Φ x)
    {x : ℝ} (hx : 1 / (2 * √(N : ℝ)) ≤ |x|) (hx' : |x| ≤ π) :
    Φ x ≤ C * 2 ^ 30 / (N : ℝ) ^ 14 := by
  have hNpos : (0 : ℝ) < N := by exact_mod_cast hN
  have hs : 0 < √(N : ℝ) := Real.sqrt_pos.2 hNpos
  have hss : √(N : ℝ) * √(N : ℝ) = N := Real.mul_self_sqrt hNpos.le
  have h1 : √(N : ℝ) / 2 ≤ N * |x| := by
    have := mul_le_mul_of_nonneg_left hx hNpos.le
    have e : (N : ℝ) * (1 / (2 * √(N : ℝ))) = √(N : ℝ) / 2 := by
      field_simp; rw [Real.sq_sqrt hNpos.le]
    linarith
  have h2 : (√(N : ℝ) / 2) ^ 30 ≤ (1 + N * |x|) ^ 30 :=
    pow_le_pow_left₀ (by positivity) (by linarith) 30
  have h3 : (√(N : ℝ) / 2) ^ 30 = (N : ℝ) ^ 15 / 2 ^ 30 := by
    rw [div_pow, show (30 : ℕ) = 2 * 15 by norm_num, pow_mul, Real.sq_sqrt hNpos.le]
  have hpos : 0 < (N : ℝ) ^ 15 / 2 ^ 30 := by positivity
  have h4 : Φ x * ((N : ℝ) ^ 15 / 2 ^ 30) ≤ C * N := by
    calc Φ x * ((N : ℝ) ^ 15 / 2 ^ 30) ≤ Φ x * (1 + N * |x|) ^ 30 := by
          rw [← h3]; exact mul_le_mul_of_nonneg_left h2 (hΦ x)
      _ ≤ C * N := hdec x hx'
  rw [le_div_iff₀ (by positivity)]
  have e : Φ x * (N : ℝ) ^ 14 * N = Φ x * ((N : ℝ) ^ 15 / 2 ^ 30) * 2 ^ 30 := by
    field_simp
  have : Φ x * (N : ℝ) ^ 14 * N ≤ C * 2 ^ 30 * N := by
    rw [e]; nlinarith
  exact le_of_mul_le_mul_right this hNpos

end LeanLi
