import Riemann.Deep.LeanNBRate.Contour6
import Mathlib.NumberTheory.LSeries.ZetaZeros

/-!
# U1 contour argument, part 7: far-left limit, summability of the residues, the middle rectangle
-/

open Complex Filter Topology Set MeasureTheory intervalIntegral
open scoped Real Interval Nat

noncomputable section

namespace Riemann.Deep.LeanNBRate

variable {H : ℂ → ℂ} {M β D : ℝ} {N : ℕ}

theorem tendsto_lin_exp_zero {κ : ℝ} (hκ : 0 < κ) :
    Tendsto (fun K : ℕ => ((4 * K + 5 : ℝ)) * Real.exp (-(κ * (2 * K + 3 / 2)))) atTop (𝓝 0) := by
  set r := Real.exp (-(2 * κ)) with hr_def
  have hr0 : 0 ≤ r := (Real.exp_pos _).le
  have hr1 : r < 1 := Real.exp_lt_one_iff.mpr (by linarith)
  have hrn : ‖r‖ < 1 := by rw [Real.norm_eq_abs, abs_of_nonneg hr0]; exact hr1
  have h1 := summable_pow_mul_geometric_of_norm_lt_one (R := ℝ) 1 hrn
  have h0 := summable_geometric_of_lt_one hr0 hr1
  have hs : Summable fun K : ℕ => (4 * ((K : ℝ) ^ 1 * r ^ K) + 5 * r ^ K) * Real.exp (-(3 / 2 * κ)) :=
    ((h1.mul_left 4).add (h0.mul_left 5)).mul_right _
  refine hs.tendsto_atTop_zero.congr fun K => ?_
  have e : Real.exp (-(κ * (2 * K + 3 / 2))) = r ^ K * Real.exp (-(3 / 2 * κ)) := by
    rw [hr_def, ← Real.exp_nat_mul, ← Real.exp_add]; congr 1; ring
  rw [e]; ring

/-- Far left: VI(Psi, hl K) -> 0 when L > beta + log(2 pi). -/
theorem tendsto_VI_Psi_hl (hH : HBound H M β D N) {L : ℝ} (hL : 0 ≤ L) {s : ℂ} (hs : 3 < s.re)
    (hLβ : β + Real.log (2 * π) < L) : Tendsto (fun K : ℕ => VI (Psi H L s) (hl K)) atTop (𝓝 0) := by
  set κ := L - β - Real.log (2 * π) with hκ
  have hκ0 : 0 < κ := by rw [hκ]; linarith
  set c₁ := M * (12 * Real.sqrt (2 * π)) * (D / min D 1) ^ N * Real.exp (β / 2 - L * s.re) with hc₁
  set cI := ∫ y : ℝ, ((1 + |y|) ^ 2)⁻¹ with hcI
  have hCK : ∀ K : ℕ, M * Real.exp (β * (1 / 2 - hl K)) * (12 * Real.sqrt (2 * π) * (2 * π) ^ (-hl K) *
      (4 * K + 5)) * Real.exp (-(L * (s.re - hl K))) * (D / min D 1) ^ N =
      c₁ * (((4 * K + 5 : ℝ)) * Real.exp (-(κ * (2 * K + 3 / 2)))) := by
    intro K
    have h2pi : (0 : ℝ) < 2 * π := by positivity
    rw [Real.rpow_def_of_pos h2pi, hc₁, hκ]
    have e : Real.exp (β * (1 / 2 - hl K)) * Real.exp (Real.log (2 * π) * -hl K) *
        Real.exp (-(L * (s.re - hl K))) = Real.exp (β / 2 - L * s.re) *
        Real.exp (-((L - β - Real.log (2 * π)) * (2 * K + 3 / 2))) := by
      rw [← Real.exp_add, ← Real.exp_add, ← Real.exp_add]; congr 1; unfold hl; ring
    calc M * Real.exp (β * (1 / 2 - hl K)) * (12 * Real.sqrt (2 * π) * Real.exp (Real.log (2 * π) * -hl K) *
        (4 * K + 5)) * Real.exp (-(L * (s.re - hl K))) * (D / min D 1) ^ N
        = M * (12 * Real.sqrt (2 * π)) * (D / min D 1) ^ N * (4 * K + 5) *
          (Real.exp (β * (1 / 2 - hl K)) * Real.exp (Real.log (2 * π) * -hl K) *
          Real.exp (-(L * (s.re - hl K)))) := by ring
      _ = _ := by rw [e]; ring
  have hbound : ∀ K : ℕ, ‖VI (Psi H L s) (hl K)‖ ≤ cI * (c₁ * (((4 * K + 5 : ℝ)) * Real.exp (-(κ * (2 * K + 3 / 2))))) := by
    intro K
    rw [VI]
    refine (MeasureTheory.norm_integral_le_integral_norm _).trans ?_
    rw [← hCK K, mul_comm, ← MeasureTheory.integral_const_mul]
    refine MeasureTheory.integral_mono (integrable_Psi_hl hH hL hs K).norm
      (integrable_inv_one_add_abs_sq.const_mul _) fun y => ?_
    exact norm_Psi_hl hH hL hs K y
  have ht := (tendsto_lin_exp_zero hκ0).const_mul c₁ |>.const_mul cI
  simp only [mul_zero] at ht
  exact squeeze_zero_norm hbound ht

theorem summable_Ares (hH : HBound H M β D N) {L : ℝ} (hL : 0 ≤ L) {s : ℂ} (hs : 3 < s.re) :
    Summable fun j : ℕ => Ares H L s j := by
  set r := Real.exp (2 * β + 2 * Real.log (2 * π) - 2 * L) with hr
  set c₂ := 6 * M * Real.exp (5 / 2 * β + 2 * Real.log (2 * π) - L * (s.re + 2)) with hc₂
  have hM := hH.M_nonneg
  have hD := hH.D_pos
  refine Summable.of_norm_bounded (g := fun j : ℕ => c₂ * (r ^ j / j !)) ((Real.summable_pow_div_factorial r).mul_left c₂) ?_
  intro j
  have htz : (tz j).re = -2 * j - 2 := tz_re j
  have hj0 : (0 : ℝ) ≤ j := Nat.cast_nonneg j
  have hsz : 1 ≤ ‖s - tz j‖ := by
    have h1 := Complex.abs_re_le_norm (s - tz j)
    rw [sub_re, htz] at h1
    have h2 : (0 : ℝ) < s.re - (-2 * j - 2) := by linarith
    rw [abs_of_pos h2] at h1
    linarith
  have hk : ‖kap L (s - tz j)‖ = Real.exp (-(L * (s.re + 2 * j + 2))) := by
    rw [norm_kap, sub_re, htz]; congr 2; ring
  have hHb := hH.bound (tz j)
  have hw : (D / (D + ‖tz j - 1 / 2‖)) ^ N ≤ 1 :=
    pow_le_one₀ (by positivity) (div_le_one_of_le₀ (by linarith [norm_nonneg (tz j - 1 / 2)]) (by positivity))
  have hmax : max 0 (1 / 2 - (tz j).re) = 2 * j + 5 / 2 := by
    rw [htz, max_eq_right (by linarith [(Nat.cast_nonneg j : (0:ℝ) ≤ j)])]; ring
  rw [hmax] at hHb
  have hd := norm_deriv_zeta_neg_two_ge j
  have hfpos : (0 : ℝ) < ((2 * j + 2).factorial : ℝ) := by exact_mod_cast Nat.factorial_pos _
  have hdpos : (0 : ℝ) < ((2 * j + 2).factorial : ℝ) / (6 * (2 * π) ^ (2 * j + 2)) := by positivity
  have hrk : ‖rk H j‖ ≤ M * Real.exp (β * (2 * j + 5 / 2)) * (6 * (2 * π) ^ (2 * j + 2) / (2 * j + 2).factorial) := by
    rw [rk, norm_div]
    have h1 : ‖H (-2 * ((j : ℂ) + 1))‖ ≤ M * Real.exp (β * (2 * j + 5 / 2)) := by
      have := hHb; rw [tz] at this
      calc ‖H (-2 * ((j : ℂ) + 1))‖ ≤ M * Real.exp (β * (2 * j + 5 / 2)) * (D / (D + ‖tz j - 1 / 2‖)) ^ N := this
        _ ≤ M * Real.exp (β * (2 * j + 5 / 2)) * 1 := by gcongr
        _ = _ := mul_one _
    rw [div_le_iff₀ (lt_of_lt_of_le hdpos hd)]
    calc ‖H (-2 * ((j : ℂ) + 1))‖ ≤ M * Real.exp (β * (2 * j + 5 / 2)) := h1
      _ = M * Real.exp (β * (2 * j + 5 / 2)) * (6 * (2 * π) ^ (2 * j + 2) / (2 * j + 2).factorial) *
          (((2 * j + 2).factorial : ℝ) / (6 * (2 * π) ^ (2 * j + 2))) := by
          field_simp
      _ ≤ _ := by gcongr
  have hfac : ((j ! : ℕ) : ℝ) ≤ ((2 * j + 2).factorial : ℝ) := by
    exact_mod_cast Nat.factorial_le (by omega)
  have hjpos : (0 : ℝ) < (j ! : ℝ) := by exact_mod_cast Nat.factorial_pos _
  rw [Ares, norm_div, norm_mul, hk]
  calc ‖rk H j‖ * Real.exp (-(L * (s.re + 2 * j + 2))) / ‖s - tz j‖
      ≤ ‖rk H j‖ * Real.exp (-(L * (s.re + 2 * j + 2))) := div_le_self (by positivity) hsz
    _ ≤ M * Real.exp (β * (2 * j + 5 / 2)) * (6 * (2 * π) ^ (2 * j + 2) / (2 * j + 2).factorial) *
        Real.exp (-(L * (s.re + 2 * j + 2))) := by gcongr
    _ ≤ M * Real.exp (β * (2 * j + 5 / 2)) * (6 * (2 * π) ^ (2 * j + 2) / j !) *
        Real.exp (-(L * (s.re + 2 * j + 2))) := by gcongr
    _ = c₂ * (r ^ j / j !) := by
        have h2pi : (0 : ℝ) < 2 * π := by positivity
        have ep : (2 * π) ^ (2 * j + 2) = Real.exp (((2 * j + 2 : ℕ) : ℝ) * Real.log (2 * π)) := by
          rw [Real.exp_nat_mul, Real.exp_log h2pi]
        have er : r ^ j = Real.exp (j * (2 * β + 2 * Real.log (2 * π) - 2 * L)) := by
          rw [hr, Real.exp_nat_mul]
        rw [ep, er, hc₂]
        have e : Real.exp (β * (2 * j + 5 / 2)) * Real.exp (((2 * j + 2 : ℕ) : ℝ) * Real.log (2 * π)) *
            Real.exp (-(L * (s.re + 2 * j + 2))) = Real.exp (5 / 2 * β + 2 * Real.log (2 * π) - L * (s.re + 2)) *
            Real.exp (j * (2 * β + 2 * Real.log (2 * π) - 2 * L)) := by
          rw [← Real.exp_add, ← Real.exp_add, ← Real.exp_add]; congr 1; push_cast; ring
        calc M * Real.exp (β * (2 * j + 5 / 2)) * (6 * Real.exp (((2 * j + 2 : ℕ) : ℝ) * Real.log (2 * π)) / j !) *
            Real.exp (-(L * (s.re + 2 * j + 2)))
            = 6 * M * (Real.exp (β * (2 * j + 5 / 2)) * Real.exp (((2 * j + 2 : ℕ) : ℝ) * Real.log (2 * π)) *
              Real.exp (-(L * (s.re + 2 * j + 2)))) / j ! := by ring
          _ = _ := by rw [e]; ring

/-- The left half-plane: I VI(Psi, hl 0) = 2 pi i R_tr(s). -/
theorem VI_Psi_hl0 (hH : HBound H M β D N) {L : ℝ} (hL : 0 ≤ L) {s : ℂ} (hs : 3 < s.re)
    (hLβ : β + Real.log (2 * π) < L) : I * VI (Psi H L s) (hl 0) = 2 * π * I * Rtr H L s := by
  have hsum := (summable_Ares hH hL hs).hasSum
  rw [Rtr_eq]
  have h1 : Tendsto (fun K : ℕ => 2 * π * I * ∑ j ∈ Finset.range K, Ares H L s j) atTop
      (𝓝 (2 * π * I * ∑' j, Ares H L s j)) := hsum.tendsto_sum_nat.const_mul _
  have h2 : Tendsto (fun K : ℕ => I * VI (Psi H L s) (hl 0) - I * VI (Psi H L s) (hl K)) atTop
      (𝓝 (I * VI (Psi H L s) (hl 0) - I * 0)) :=
    tendsto_const_nhds.sub ((tendsto_VI_Psi_hl hH hL hs hLβ).const_mul I)
  rw [mul_zero, sub_zero] at h2
  have h3 := h2.congr fun K => VI_Psi_telescope hH hL hs K
  exact tendsto_nhds_unique h3 h1

end Riemann.Deep.LeanNBRate

