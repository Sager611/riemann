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


theorem norm_two_pi_cpow (s : ℂ) : ‖(2 * π : ℂ) ^ s‖ = (2 * π) ^ s.re := by
  rw [show (2 * π : ℂ) = ((2 * π : ℝ) : ℂ) by push_cast; ring]
  exact Complex.norm_cpow_eq_rpow_re_of_pos (by positivity) s

/-- Left half-plane: |1/zeta(s)| prod_{j<m}|s+j| <= 6 sqrt(2 pi) (2 pi)^{-Re s} (|s|+m+1). -/
theorem norm_inv_zeta_left {s : ℂ} (hs : s.re ≤ -1) (hsn : ∀ n : ℕ, s ≠ -n) (m : ℕ)
    (hm1 : 1 / 2 ≤ s.re + m) (hm2 : s.re + m ≤ 3 / 2) :
    ‖(riemannZeta s)⁻¹‖ * ∏ j ∈ Finset.range m, ‖s + j‖ ≤
      6 * Real.sqrt (2 * π) * (2 * π) ^ (-s.re) * (‖s‖ + m + 1) := by
  have hj : ∀ j < m, s + j ≠ 0 := fun j _ h => hsn j (by rw [← sub_eq_zero, sub_neg_eq_add, h])
  rw [(inv_zeta_left (by linarith) hsn).2]
  have h3 : ‖(riemannZeta (1 - s))⁻¹‖ ≤ 3 := norm_inv_zeta_le_three (by simp; linarith)
  have hg := norm_Gamma_mul_cos_le m hm1 hm2 hj
  rw [norm_mul, norm_mul, norm_mul, norm_mul, norm_two_pi_cpow]
  simp only [Complex.norm_ofNat, neg_re]
  calc 2 * (2 * π) ^ (-s.re) * ‖Complex.Gamma s‖ * ‖Complex.cos (π * s / 2)‖ *
        ‖(riemannZeta (1 - s))⁻¹‖ * ∏ j ∈ Finset.range m, ‖s + j‖
      = 2 * (2 * π) ^ (-s.re) * ‖(riemannZeta (1 - s))⁻¹‖ *
        ((‖Complex.Gamma s‖ * ‖Complex.cos (π * s / 2)‖) * ∏ j ∈ Finset.range m, ‖s + j‖) := by ring
    _ ≤ 2 * (2 * π) ^ (-s.re) * 3 * (Real.sqrt (2 * π) * (‖s‖ + m + 1)) := by
        rw [← norm_mul]
        gcongr
    _ = 6 * Real.sqrt (2 * π) * (2 * π) ^ (-s.re) * (‖s‖ + m + 1) := by ring

theorem sin_pi_odd_half (k : ℕ) : Complex.sin (π * ((2 * k + 3 : ℕ) : ℂ) / 2) = (-1) ^ (k + 1) := by
  rw [show (π : ℂ) * ((2 * k + 3 : ℕ) : ℂ) / 2 = π * ((k + 1 : ℕ) : ℂ) + π / 2 by push_cast; ring,
    Complex.sin_add_pi_div_two, cos_pi_mul_nat]

theorem sin_pi_mul_nat (n : ℕ) : Complex.sin (π * (n : ℂ)) = 0 := by
  have := Real.sin_nat_mul_pi n
  rw [show (π : ℂ) * (n : ℂ) = ((n * π : ℝ) : ℂ) by push_cast; ring, ← Complex.ofReal_sin, this]
  simp

theorem cos_pi_odd_half (k : ℕ) : Complex.cos (π * ((2 * k + 3 : ℕ) : ℂ) / 2) = 0 := by
  rw [show (π : ℂ) * ((2 * k + 3 : ℕ) : ℂ) / 2 = π * ((k + 1 : ℕ) : ℂ) + π / 2 by push_cast; ring,
    Complex.cos_add_pi_div_two, sin_pi_mul_nat, neg_zero]

/-- The derivative of zeta at the trivial zero -2(k+1). -/
theorem deriv_zeta_neg_two (k : ℕ) :
    deriv riemannZeta (-2 * ((k : ℂ) + 1)) =
      2 * ((2 * π : ℂ) ^ (2 * k + 3))⁻¹ * ((2 * k + 2).factorial : ℂ) * riemannZeta ((2 * k + 3 : ℕ) : ℂ) *
        ((-1) ^ (k + 1) * (π / 2)) := by
  set s₀ : ℂ := ((2 * k + 3 : ℕ) : ℂ) with hs₀
  have h1s : (1 : ℂ) - s₀ = -2 * ((k : ℂ) + 1) := by rw [hs₀]; push_cast; ring
  have hre : 1 < s₀.re := by rw [hs₀]; simp; linarith [(Nat.cast_nonneg k : (0 : ℝ) ≤ k)]
  set B : ℂ → ℂ := fun s => 2 * (2 * π : ℂ) ^ (-s) * Complex.Gamma s * riemannZeta s with hB
  set c : ℂ → ℂ := fun s => Complex.cos (π * s / 2) with hc
  have hev : (fun s => riemannZeta (1 - s)) =ᶠ[𝓝 s₀] fun s => B s * c s := by
    filter_upwards [(isOpen_lt continuous_const Complex.continuous_re).mem_nhds hre] with s hs
    have hsn : ∀ n : ℕ, s ≠ -n := fun n h => by
      rw [h] at hs; simp at hs; linarith [(Nat.cast_nonneg n : (0 : ℝ) ≤ n)]
    have hs1 : s ≠ 1 := fun h => by rw [h] at hs; simp at hs
    rw [riemannZeta_one_sub hsn hs1, hB, hc]; ring
  have hb : (2 * π : ℂ) ≠ 0 := by exact_mod_cast (by positivity : (2 * π : ℝ) ≠ 0)
  have hBd : DifferentiableAt ℂ B s₀ := by
    have hG : DifferentiableAt ℂ Complex.Gamma s₀ := Complex.differentiableAt_Gamma _ fun m h => by
      have := congrArg Complex.re h; rw [hs₀] at this; simp at this
      linarith [(Nat.cast_nonneg m : (0 : ℝ) ≤ m), (Nat.cast_nonneg k : (0 : ℝ) ≤ k)]
    have hZ : DifferentiableAt ℂ riemannZeta s₀ := differentiableAt_riemannZeta fun h => by
      rw [h] at hre; simp at hre
    have hP : DifferentiableAt ℂ (fun s : ℂ => (2 * π : ℂ) ^ (-s)) s₀ :=
      (differentiableAt_id.neg).const_cpow (Or.inl hb)
    exact ((hP.const_mul 2).mul hG).mul hZ
  have hcd : HasDerivAt c (-Complex.sin (π * s₀ / 2) * (π / 2)) s₀ := by
    have := (((hasDerivAt_id s₀).const_mul (π : ℂ)).div_const 2).ccos
    simpa only [id, mul_one] using this
  have hv : HasDerivAt (fun s => B s * c s) (deriv B s₀ * c s₀ + B s₀ * (-Complex.sin (π * s₀ / 2) * (π / 2))) s₀ :=
    hBd.hasDerivAt.mul hcd
  have hu : HasDerivAt (fun s => riemannZeta (1 - s)) (deriv riemannZeta (1 - s₀) * (-1)) s₀ := by
    have hz : DifferentiableAt ℂ riemannZeta (1 - s₀) := differentiableAt_riemannZeta fun h => by
      rw [h1s] at h; have := congrArg Complex.re h; simp at this; linarith [(Nat.cast_nonneg k : (0 : ℝ) ≤ k)]
    exact hz.hasDerivAt.comp s₀ ((hasDerivAt_id s₀).const_sub 1)
  have heq := hu.unique (hv.congr_of_eventuallyEq hev)
  have hc0 : c s₀ = 0 := by rw [hc, hs₀]; exact cos_pi_odd_half k
  rw [hc0, mul_zero, zero_add] at heq
  rw [← h1s]
  have e2 : deriv riemannZeta (1 - s₀) = B s₀ * (Complex.sin (π * s₀ / 2) * (π / 2)) := by
    linear_combination -heq
  rw [e2, hB]
  simp only
  rw [hs₀, sin_pi_odd_half, Complex.cpow_neg, Complex.cpow_natCast,
    show ((2 * k + 3 : ℕ) : ℂ) = ((2 * k + 2 : ℕ) : ℂ) + 1 by push_cast; ring, Complex.Gamma_nat_eq_factorial]
    <;> ring

theorem norm_deriv_zeta_neg_two_ge (k : ℕ) :
    ((2 * k + 2).factorial : ℝ) / (6 * (2 * π) ^ (2 * k + 2)) ≤ ‖deriv riemannZeta (-2 * ((k : ℂ) + 1))‖ := by
  rw [deriv_zeta_neg_two]
  have hz := norm_zeta_ge_third (s := ((2 * k + 3 : ℕ) : ℂ)) (by simp; linarith [(Nat.cast_nonneg k : (0 : ℝ) ≤ k)])
  simp only [norm_mul, norm_inv, norm_pow, norm_neg, norm_one, one_pow, Complex.norm_ofNat, norm_div,
    Complex.norm_natCast, one_mul]
  have hpi' : ‖(π : ℂ)‖ = π := by rw [Complex.norm_real, Real.norm_eq_abs, abs_of_pos Real.pi_pos]
  rw [hpi']
  have hp : (0 : ℝ) < (2 * π) ^ (2 * k + 2) := by positivity
  have hf : (0 : ℝ) < ((2 * k + 2).factorial : ℝ) := by exact_mod_cast Nat.factorial_pos _
  rw [div_le_iff₀ (by positivity), pow_succ]
  have h2pi : (0 : ℝ) < 2 * π := by positivity
  field_simp
  nlinarith [mul_pos hf hp, Real.pi_pos]

/-- From Titchmarsh (2.1.5): |zeta(s)| <= |s|(1/sigma + 1/(1-sigma)) for 0 < sigma < 1. -/
theorem norm_zeta_le_of_strip (h215 : Titchmarsh_2_1_5) {s : ℂ} (h0 : 0 < s.re) (h1 : s.re < 1) :
    ‖riemannZeta s‖ ≤ ‖s‖ * (1 / s.re + 1 / (1 - s.re)) := by
  rw [h215 s h0 h1, norm_mul]
  refine mul_le_mul_of_nonneg_left ?_ (norm_nonneg _)
  set g : ℝ → ℝ := (Ioc (0 : ℝ) 1).indicator (fun x => x ^ (-s.re)) + (Ioi (1 : ℝ)).indicator
    (fun x => x ^ (-s.re - 1)) with hg
  have hgi1 : IntegrableOn (fun x : ℝ => x ^ (-s.re)) (Ioc 0 1) := by
    have := intervalIntegral.intervalIntegrable_rpow' (a := 0) (b := 1) (r := -s.re) (by linarith)
    exact (intervalIntegrable_iff_integrableOn_Ioc_of_le zero_le_one).mp this
  have hgi2 : IntegrableOn (fun x : ℝ => x ^ (-s.re - 1)) (Ioi 1) :=
    integrableOn_Ioi_rpow_of_lt (by linarith) one_pos
  have hgint : Integrable g (volume.restrict (Ioi 0)) := by
    refine Integrable.add ?_ ?_
    · exact (integrable_indicator_iff measurableSet_Ioc).mpr hgi1 |>.restrict
    · exact (integrable_indicator_iff measurableSet_Ioi).mpr hgi2 |>.restrict
  have hle : ∀ᵐ x : ℝ ∂(volume.restrict (Ioi (0 : ℝ))),
      ‖((⌊x⌋ : ℂ) - x) * (x : ℂ) ^ (-s - 1)‖ ≤ g x := by
    filter_upwards [ae_restrict_mem measurableSet_Ioi] with x hx
    have hx0 : (0 : ℝ) < x := hx
    rw [norm_mul, Complex.norm_cpow_eq_rpow_re_of_pos hx0]
    have hfl : ‖((⌊x⌋ : ℂ) - x)‖ = Int.fract x := by
      rw [show ((⌊x⌋ : ℂ) - x) = (((⌊x⌋ : ℝ) - x : ℝ) : ℂ) by push_cast; ring, Complex.norm_real,
        Real.norm_eq_abs, abs_sub_comm, Int.self_sub_floor, abs_of_nonneg (Int.fract_nonneg x)]
    rw [hfl]
    have hre : (-s - 1).re = -s.re - 1 := by simp
    rw [hre, hg]
    simp only [Pi.add_apply]
    by_cases hx1 : x ≤ 1
    · rw [indicator_of_mem (show x ∈ Ioc (0 : ℝ) 1 from ⟨hx0, hx1⟩),
        indicator_of_notMem (show x ∉ Ioi (1 : ℝ) from fun h => by simp at h; linarith), add_zero]
      have hf : Int.fract x ≤ x := by
        have h0' := Int.floor_nonneg.mpr hx0.le
        have h' : (0 : ℝ) ≤ (⌊x⌋ : ℝ) := by exact_mod_cast h0'
        rw [← Int.self_sub_floor]; linarith
      calc Int.fract x * x ^ (-s.re - 1) ≤ x * x ^ (-s.re - 1) :=
            mul_le_mul_of_nonneg_right hf (by positivity)
        _ = x ^ (-s.re) := by
            rw [Real.rpow_sub hx0, Real.rpow_one]; field_simp
    · push Not at hx1
      rw [indicator_of_notMem (show x ∉ Ioc (0 : ℝ) 1 from fun h => by simp at h; linarith [h.2]),
        indicator_of_mem (show x ∈ Ioi (1 : ℝ) from hx1), zero_add]
      calc Int.fract x * x ^ (-s.re - 1) ≤ 1 * x ^ (-s.re - 1) :=
            mul_le_mul_of_nonneg_right (Int.fract_lt_one x).le (by positivity)
        _ = x ^ (-s.re - 1) := one_mul _
  refine (norm_integral_le_of_norm_le hgint hle).trans (le_of_eq ?_)
  rw [show (∫ x in Ioi (0 : ℝ), g x) = ∫ x in Ioi (0 : ℝ), ((Ioc (0 : ℝ) 1).indicator (fun x => x ^ (-s.re)) x +
    (Ioi (1 : ℝ)).indicator (fun x => x ^ (-s.re - 1)) x) from rfl]
  rw [integral_add ((integrable_indicator_iff measurableSet_Ioc).mpr hgi1 |>.restrict)
    ((integrable_indicator_iff measurableSet_Ioi).mpr hgi2 |>.restrict)]
  rw [setIntegral_indicator measurableSet_Ioc, setIntegral_indicator measurableSet_Ioi]
  rw [show Ioi (0 : ℝ) ∩ Ioc 0 1 = Ioc 0 1 from inter_eq_right.mpr Ioc_subset_Ioi_self,
    show Ioi (0 : ℝ) ∩ Ioi 1 = Ioi 1 from inter_eq_right.mpr (Ioi_subset_Ioi zero_le_one)]
  rw [integral_Ioi_rpow_of_lt (by linarith) one_pos, ← intervalIntegral.integral_of_le zero_le_one,
    integral_rpow (Or.inl (by linarith))]
  simp only [Real.one_rpow, Real.zero_rpow (show -s.re + 1 ≠ 0 by linarith)]
  have h1' : (1 : ℝ) - s.re ≠ 0 := by linarith
  have h0' : s.re ≠ 0 := h0.ne'
  rw [show -s.re + 1 = 1 - s.re by ring, show -s.re - 1 + 1 = -s.re by ring]
  first | (field_simp; done) | (field_simp; ring)

/-- Cauchy estimate: |zeta'(z)| <= 64|z| + 8 for 3/8 <= Re z <= 5/8. -/
theorem norm_deriv_zeta_le (h215 : Titchmarsh_2_1_5) {z : ℂ} (h1 : 3 / 8 ≤ z.re) (h2 : z.re ≤ 5 / 8) :
    ‖deriv riemannZeta z‖ ≤ 64 * ‖z‖ + 8 := by
  have hd : DiffContOnCl ℂ riemannZeta (Metric.ball z (1 / 8)) := by
    refine DifferentiableOn.diffContOnCl ?_
    rw [closure_ball z (by norm_num)]
    intro w hw
    refine (differentiableAt_riemannZeta fun h => ?_).differentiableWithinAt
    rw [Metric.mem_closedBall, h, dist_comm] at hw
    have := (Complex.abs_re_le_norm (z - 1)).trans (show ‖z - 1‖ ≤ 1 / 8 by rw [← dist_eq_norm]; exact hw)
    simp at this
    rw [abs_le] at this
    linarith [this.1]
  have hb : ∀ w ∈ Metric.sphere z (1 / 8), ‖riemannZeta w‖ ≤ 8 * ‖z‖ + 1 := by
    intro w hw
    rw [mem_sphere_iff_norm] at hw
    have hre : |w.re - z.re| ≤ 1 / 8 := by
      have := Complex.abs_re_le_norm (w - z); simp at this; linarith
    rw [abs_le] at hre
    have hw0 : 0 < w.re := by linarith
    have hw1 : w.re < 1 := by linarith
    refine (norm_zeta_le_of_strip h215 hw0 hw1).trans ?_
    have hn : ‖w‖ ≤ ‖z‖ + 1 / 8 := by
      calc ‖w‖ = ‖z + (w - z)‖ := by ring_nf
        _ ≤ ‖z‖ + ‖w - z‖ := norm_add_le _ _
        _ = ‖z‖ + 1 / 8 := by rw [hw]
    have hs : 1 / w.re + 1 / (1 - w.re) ≤ 8 := by
      have a1 : 1 / w.re ≤ 4 := by rw [div_le_iff₀ hw0]; linarith
      have a2 : 1 / (1 - w.re) ≤ 4 := by rw [div_le_iff₀ (by linarith)]; linarith
      linarith
    calc ‖w‖ * (1 / w.re + 1 / (1 - w.re)) ≤ (‖z‖ + 1 / 8) * 8 :=
          mul_le_mul hn hs (by positivity) (by positivity)
      _ = 8 * ‖z‖ + 1 := by ring
  have := Complex.norm_deriv_le_of_forall_mem_sphere_norm_le (by norm_num) hd hb
  calc ‖deriv riemannZeta z‖ ≤ (8 * ‖z‖ + 1) / (1 / 8) := this
    _ = 64 * ‖z‖ + 8 := by ring


end Riemann.Deep.LeanNBRate
