import Riemann.Deep.LeanNBRate.GammaBd
import Riemann.Deep.LeanNBRate.LitInputs
import Mathlib.NumberTheory.LSeries.Dirichlet
import Mathlib.NumberTheory.ZetaValues
import Mathlib.Analysis.Real.Pi.Bounds
import Mathlib.Analysis.Complex.Liouville
import Mathlib.Analysis.SpecialFunctions.Integrals.Basic

/-!
# Bounds and facts for the Riemann zeta function used in the construction (U1)
-/

open Complex Filter Topology Set MeasureTheory
open scoped Real

noncomputable section

namespace Riemann.Deep.LeanNBRate

theorem Gamma_add_nat_eq_prod (z : ℂ) (n : ℕ) (hz : ∀ j < n, z + j ≠ 0) :
    Complex.Gamma (z + n) = Complex.Gamma z * ∏ j ∈ Finset.range n, (z + j) := by
  induction n with
  | zero => simp
  | succ n ih =>
    rw [Finset.prod_range_succ, ← mul_assoc, ← ih (fun j hj => hz j (by omega))]
    rw [show z + ((n + 1 : ℕ) : ℂ) = (z + n) + 1 by push_cast; ring]
    rw [Complex.Gamma_add_one _ (hz n (by omega))]
    ring

theorem norm_cos_le_exp (z : ℂ) : ‖Complex.cos z‖ ≤ Real.exp |z.im| := by
  rw [Complex.cos]
  have h1 : ‖cexp (z * I)‖ = Real.exp (-z.im) := by rw [Complex.norm_exp]; simp
  have h2 : ‖cexp (-z * I)‖ = Real.exp z.im := by rw [Complex.norm_exp]; simp
  calc ‖(cexp (z * I) + cexp (-z * I)) / 2‖ ≤ (‖cexp (z * I)‖ + ‖cexp (-z * I)‖) / 2 := by
        rw [norm_div]; simp only [Complex.norm_ofNat]
        exact div_le_div_of_nonneg_right (norm_add_le _ _) (by norm_num)
    _ = (Real.exp (-z.im) + Real.exp z.im) / 2 := by rw [h1, h2]
    _ ≤ Real.exp |z.im| := by
        have a1 : Real.exp (-z.im) ≤ Real.exp |z.im| := Real.exp_le_exp.mpr (neg_le_abs _)
        have a2 : Real.exp z.im ≤ Real.exp |z.im| := Real.exp_le_exp.mpr (le_abs_self _)
        linarith

/-- |Gamma(s) cos(pi s/2)| prod |s+j| <= sqrt(2 pi)(|s|+m+1) when Re(s+m) in [1/2,3/2]. -/
theorem norm_Gamma_mul_cos_le {s : ℂ} (m : ℕ) (hm1 : 1 / 2 ≤ s.re + m) (hm2 : s.re + m ≤ 3 / 2)
    (hs : ∀ j < m, s + j ≠ 0) :
    ‖Complex.Gamma s * Complex.cos (π * s / 2)‖ * ∏ j ∈ Finset.range m, ‖s + j‖ ≤
      Real.sqrt (2 * π) * (‖s‖ + m + 1) := by
  have hG := Gamma_add_nat_eq_prod s m hs
  have hre : (s + m).re = s.re + m := by simp
  have hb := norm_Gamma_strip (z := s + m) (by rw [hre]; exact hm1) (by rw [hre]; exact hm2)
  rw [hG, norm_mul, norm_prod] at hb
  have him : (s + (m : ℂ)).im = s.im := by simp
  rw [him] at hb
  have hc := norm_cos_le_exp (π * s / 2)
  have hcim : (π * s / 2 : ℂ).im = π / 2 * s.im := by simp; ring
  rw [hcim, abs_mul, abs_of_pos (by positivity : (0 : ℝ) < π / 2)] at hc
  have hn : ‖s + (m : ℂ) + 1‖ ≤ ‖s‖ + m + 1 := by
    calc ‖s + (m : ℂ) + 1‖ ≤ ‖s‖ + ‖(m : ℂ)‖ + ‖(1 : ℂ)‖ := norm_add₃_le
      _ = ‖s‖ + m + 1 := by simp
  rw [norm_mul]
  calc ‖Complex.Gamma s‖ * ‖Complex.cos (π * s / 2)‖ * ∏ j ∈ Finset.range m, ‖s + j‖
      = (‖Complex.Gamma s‖ * ∏ j ∈ Finset.range m, ‖s + j‖) * ‖Complex.cos (π * s / 2)‖ := by ring
    _ ≤ (Real.sqrt (2 * π) * ‖s + (m : ℂ) + 1‖ * Real.exp (-(π / 2 * |s.im|))) *
        Real.exp (π / 2 * |s.im|) := mul_le_mul hb hc (norm_nonneg _) (by positivity)
    _ = Real.sqrt (2 * π) * ‖s + (m : ℂ) + 1‖ := by
        rw [mul_assoc, ← Real.exp_add, neg_add_cancel, Real.exp_zero, mul_one]
    _ ≤ Real.sqrt (2 * π) * (‖s‖ + m + 1) := mul_le_mul_of_nonneg_left hn (by positivity)

/-- |zeta(s) - 1| <= pi^2/6 - 1 for Re s >= 2. -/
theorem norm_zeta_sub_one_le {s : ℂ} (hs : 2 ≤ s.re) : ‖riemannZeta s - 1‖ ≤ π ^ 2 / 6 - 1 := by
  have h1 : 1 < s.re := by linarith
  have hsum : Summable fun n : ℕ => 1 / ((n : ℂ) + 1) ^ s := by
    have := Complex.summable_one_div_nat_cpow.mpr h1
    exact (summable_nat_add_iff 1).mpr this |>.congr fun n => by push_cast; ring_nf
  rw [zeta_eq_tsum_one_div_nat_add_one_cpow h1, hsum.tsum_eq_zero_add]
  simp only [Nat.cast_zero, zero_add, Complex.one_cpow, div_one, add_sub_cancel_left]
  have hn : ∀ n : ℕ, ‖1 / (((n + 1 : ℕ) : ℂ) + 1) ^ s‖ = 1 / ((n : ℝ) + 2) ^ s.re := by
    intro n
    rw [norm_div, norm_one, show (((n + 1 : ℕ) : ℂ) + 1) = (((n : ℝ) + 2 : ℝ) : ℂ) by push_cast; ring,
      Complex.norm_cpow_eq_rpow_re_of_pos (by positivity)]
  have hb : ∀ n : ℕ, 1 / ((n : ℝ) + 2) ^ s.re ≤ 1 / ((n : ℝ) + 2) ^ 2 := by
    intro n
    rw [← Real.rpow_natCast]
    apply one_div_le_one_div_of_le (by positivity)
    exact Real.rpow_le_rpow_of_exponent_le (by linarith [(Nat.cast_nonneg n : (0 : ℝ) ≤ n)])
      (by norm_num; exact hs)
  have hz2 : HasSum (fun n : ℕ => 1 / ((n : ℝ) + 2) ^ 2) (π ^ 2 / 6 - 1) := by
    have := (hasSum_nat_add_iff' 2).mpr hasSum_zeta_two
    simp only [Finset.sum_range_succ, Finset.sum_range_zero, Nat.cast_zero, Nat.cast_one] at this
    norm_num at this
    convert this using 1
    funext n; push_cast; ring_nf
  have hs2 : Summable fun n : ℕ => ‖1 / (((n + 1 : ℕ) : ℂ) + 1) ^ s‖ := by
    refine Summable.of_nonneg_of_le (fun n => norm_nonneg _) (fun n => ?_) hz2.summable
    rw [hn]; exact hb n
  calc ‖∑' n : ℕ, 1 / (((n + 1 : ℕ) : ℂ) + 1) ^ s‖ ≤ ∑' n : ℕ, ‖1 / (((n + 1 : ℕ) : ℂ) + 1) ^ s‖ :=
        norm_tsum_le_tsum_norm hs2
    _ ≤ ∑' n : ℕ, 1 / ((n : ℝ) + 2) ^ 2 :=
        hs2.tsum_le_tsum (fun n => by rw [hn]; exact hb n) hz2.summable
    _ = π ^ 2 / 6 - 1 := hz2.tsum_eq

theorem pi_sq_div_six_lt : π ^ 2 / 6 < 5 / 3 := by
  have := Real.pi_lt_d2
  have h0 := Real.pi_pos
  nlinarith

theorem norm_zeta_ge_third {s : ℂ} (hs : 2 ≤ s.re) : 1 / 3 ≤ ‖riemannZeta s‖ := by
  have h := norm_zeta_sub_one_le hs
  have h2 := pi_sq_div_six_lt
  have h3 : ‖(1 : ℂ)‖ ≤ ‖riemannZeta s‖ + ‖riemannZeta s - 1‖ := by
    have := norm_sub_le (riemannZeta s) (riemannZeta s - 1)
    simpa using this
  simp at h3
  linarith

theorem norm_inv_zeta_le_three {s : ℂ} (hs : 2 ≤ s.re) : ‖(riemannZeta s)⁻¹‖ ≤ 3 := by
  have h := norm_zeta_ge_third hs
  rw [norm_inv]
  rw [inv_le_comm₀ (by linarith) (by norm_num)]
  linarith

theorem norm_zeta_le_two {s : ℂ} (hs : 2 ≤ s.re) : ‖riemannZeta s‖ ≤ 2 := by
  have h := norm_zeta_sub_one_le hs
  have h2 := pi_sq_div_six_lt
  have h3 := norm_le_insert' (riemannZeta s) 1
  have : ‖riemannZeta s‖ ≤ ‖(1 : ℂ)‖ + ‖riemannZeta s - 1‖ := by
    have := norm_add_le (1 : ℂ) (riemannZeta s - 1); simpa using this
  simp at this
  linarith

/-- In the left half-plane off the non-positive integers: zeta(s) != 0 and
1/zeta(s) = 2 (2 pi)^{-s} Gamma(s) cos(pi s/2)/zeta(1-s). -/
theorem inv_zeta_left {s : ℂ} (hs : s.re < 0) (hsn : ∀ n : ℕ, s ≠ -n) :
    riemannZeta s ≠ 0 ∧ (riemannZeta s)⁻¹ =
      2 * (2 * π : ℂ) ^ (-s) * Complex.Gamma s * Complex.cos (π * s / 2) * (riemannZeta (1 - s))⁻¹ := by
  have hs1 : s ≠ 1 := fun h => by rw [h] at hs; simp at hs; linarith
  have h := riemannZeta_one_sub hsn hs1
  have h1 : riemannZeta (1 - s) ≠ 0 := riemannZeta_ne_zero_of_one_lt_re (by simp; linarith)
  have h0 : riemannZeta s ≠ 0 := by
    intro h0; rw [h0, mul_zero] at h; exact h1 h
  refine ⟨h0, ?_⟩
  have hX : (2 * (2 * π : ℂ) ^ (-s) * Complex.Gamma s * Complex.cos (π * s / 2)) ≠ 0 := by
    intro hX; apply h1; rw [h, hX, zero_mul]
  rw [h, mul_inv, ← mul_assoc, mul_inv_cancel₀ hX, one_mul]

theorem cos_pi_mul_nat (n : ℕ) : Complex.cos (π * (n : ℂ)) = (-1) ^ n := by
  have := Real.cos_nat_mul_pi n
  rw [show (π : ℂ) * (n : ℂ) = ((n * π : ℝ) : ℂ) by push_cast; ring, ← Complex.ofReal_cos, this]
  push_cast; ring

theorem zeta_neg_odd_ne_zero (k : ℕ) : riemannZeta (-(2 * k + 1 : ℂ)) ≠ 0 := by
  have hsn : ∀ n : ℕ, ((2 * k + 2 : ℕ) : ℂ) ≠ -n := by
    intro n h; have := congrArg Complex.re h; simp at this; linarith [(Nat.cast_nonneg n : (0:ℝ) ≤ n)]
  have hs1 : ((2 * k + 2 : ℕ) : ℂ) ≠ 1 := by
    intro h; have := congrArg Complex.re h; simp at this; linarith [(Nat.cast_nonneg k : (0:ℝ) ≤ k)]
  have h := riemannZeta_one_sub hsn hs1
  have e : (1 : ℂ) - ((2 * k + 2 : ℕ) : ℂ) = -(2 * k + 1 : ℂ) := by push_cast; ring
  rw [e] at h
  rw [h]
  have hc : Complex.cos (π * ((2 * k + 2 : ℕ) : ℂ) / 2) = (-1) ^ (k + 1) := by
    rw [show (π : ℂ) * ((2 * k + 2 : ℕ) : ℂ) / 2 = π * ((k + 1 : ℕ) : ℂ) by push_cast; ring, cos_pi_mul_nat]
  rw [hc]
  have hG : Complex.Gamma ((2 * k + 2 : ℕ) : ℂ) ≠ 0 := by
    rw [show ((2 * k + 2 : ℕ) : ℂ) = ((2 * k + 1 : ℕ) : ℂ) + 1 by push_cast; ring,
      Complex.Gamma_nat_eq_factorial]
    exact_mod_cast Nat.factorial_ne_zero _
  have hz : riemannZeta ((2 * k + 2 : ℕ) : ℂ) ≠ 0 :=
    riemannZeta_ne_zero_of_one_lt_re (by simp; linarith [(Nat.cast_nonneg k : (0:ℝ) ≤ k)])
  have hb : (2 * π : ℂ) ≠ 0 := by exact_mod_cast (by positivity : (2 * π : ℝ) ≠ 0)
  have hp : (2 * π : ℂ) ^ (-((2 * k + 2 : ℕ) : ℂ)) ≠ 0 := by
    rw [Complex.cpow_def_of_ne_zero hb]; exact Complex.exp_ne_zero _
  exact mul_ne_zero (mul_ne_zero (mul_ne_zero (mul_ne_zero two_ne_zero hp) hG)
    (pow_ne_zero _ (by norm_num))) hz

/-- Zeros of zeta with negative real part are the trivial zeros. -/
theorem zeta_eq_zero_left {s : ℂ} (h : riemannZeta s = 0) (hs : s.re < 0) :
    ∃ k : ℕ, s = -2 * (k + 1) := by
  by_cases hsn : ∀ n : ℕ, s ≠ -n
  · exact absurd h (inv_zeta_left hs hsn).1
  push Not at hsn
  obtain ⟨n, rfl⟩ := hsn
  rcases Nat.even_or_odd n with ⟨k, hk⟩ | ⟨k, hk⟩
  · subst hk
    rcases k with _ | k
    · simp at hs
    · exact ⟨k, by push_cast; ring⟩
  · subst hk
    refine absurd ?_ (zeta_neg_odd_ne_zero k)
    rw [show (-(2 * (k : ℂ) + 1)) = -((2 * k + 1 : ℕ) : ℂ) by push_cast; ring]
    exact h

end Riemann.Deep.LeanNBRate

