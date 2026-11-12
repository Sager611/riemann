import Riemann.Deep.FasiHigham.Lean.Defs
import Mathlib.RingTheory.MvPowerSeries.Inverse
import Mathlib.RingTheory.PowerSeries.Basic

/-!
# Fasi–Higham, m = 3: c_n(k) is the n-th coefficient of the power-series inverse of q_{k,3}

qPS k m is q_{k,m} as an element of Mathlib's PowerSeries ℚ (coefficient function qcoef k m,
which vanishes beyond z^m); (qPS k 3)⁻¹ is Mathlib's inverse (constant coefficient 1).
-/

namespace FasiLean

open Finset

/-- q_{k,m} as a formal power series over ℚ (a polynomial: coefficients vanish beyond z^m). -/
noncomputable def qPS (k m : ℕ) : PowerSeries ℚ := PowerSeries.mk (qcoef k m)

theorem qPS_const (k : ℕ) : MvPowerSeries.constantCoeff (qPS k 3) ≠ 0 := by
  show PowerSeries.constantCoeff (qPS k 3) ≠ 0
  rw [qPS, PowerSeries.constantCoeff_mk, qcoef_zero]
  exact one_ne_zero

theorem qPS_mul_iff (k : ℕ) (b : ℕ → ℚ) : qPS k 3 * PowerSeries.mk b = 1 ↔
    ∀ n, ∑ j ∈ range (n + 1), qcoef k 3 j * b (n - j) = if n = 0 then 1 else 0 := by
  constructor
  · intro h n
    have := congrArg (PowerSeries.coeff n) h
    rw [PowerSeries.coeff_mul, Finset.Nat.sum_antidiagonal_eq_sum_range_succ_mk,
      PowerSeries.coeff_one] at this
    simpa only [qPS, PowerSeries.coeff_mk] using this
  · intro h
    ext n
    rw [PowerSeries.coeff_mul, Finset.Nat.sum_antidiagonal_eq_sum_range_succ_mk, PowerSeries.coeff_one]
    simpa only [qPS, PowerSeries.coeff_mk] using h n

theorem inv_qPS (k : ℕ) : (qPS k 3)⁻¹ = PowerSeries.mk (c k) := by
  rw [MvPowerSeries.inv_eq_iff_mul_eq_one (qPS_const k), mul_comm]
  exact (qPS_mul_iff k (c k)).2 (conv k)

/-- c_n(k) = [z^n] (q_{k,3})⁻¹ in PowerSeries ℚ. -/
theorem coeff_inv (k n : ℕ) : PowerSeries.coeff n (qPS k 3)⁻¹ = c k n := by
  rw [inv_qPS, PowerSeries.coeff_mk]

/-- Uniqueness: every b with Σ_{j ≤ n} q_j b_{n−j} = [n = 0] for all n is c. -/
theorem conv_unique (k : ℕ) (b : ℕ → ℚ)
    (hb : ∀ n, ∑ j ∈ range (n + 1), qcoef k 3 j * b (n - j) = if n = 0 then 1 else 0) :
    b = c k := by
  have h1 : (qPS k 3)⁻¹ = PowerSeries.mk b := by
    rw [MvPowerSeries.inv_eq_iff_mul_eq_one (qPS_const k), mul_comm]
    exact (qPS_mul_iff k b).2 hb
  funext n
  have := congrArg (PowerSeries.coeff n) (h1.symm.trans (inv_qPS k))
  simpa only [PowerSeries.coeff_mk] using this

end FasiLean
