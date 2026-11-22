import Riemann.Deep.LeanNBRate.Final2

/-!
# Numeric bounds: the contour mass, Q and R_tr by the explicit Qb and Rb
-/

open MeasureTheory Filter Topology Set Complex
open scoped Real

noncomputable section

namespace Riemann.Deep.LeanNBRate

theorem Z_le {H : ℂ → ℂ} {M β D T₂ L s₁ : ℝ} {N n : ℕ} (hH : HBound H M β D N) (hT₂ : 0 < T₂)
    (hL : 1 ≤ L) (hM : M ≤ 2 ^ n * (1 + 2 * s₁) * L ^ 4) (hβ : β ≤ L) (hq : D / (D + T₂) ≤ 2 / L ^ 3)
    (hD3 : 3 + 2 * D ≤ 5 * L ^ 3) (hN : N - 2 = n - 3) (hn3 : 3 ≤ n) :
    ZU H T₂ + ZL H T₂ ≤ 80 * (1 + 2 * s₁) * L ^ 7 * Real.exp (L * (3 / 2)) * (4 / L ^ 3) ^ (n - 3) := by
  have hD := hH.D_pos
  have hM0 := hH.M_nonneg
  have hq0 : 0 ≤ D / (D + T₂) := div_nonneg hD.le (by linarith)
  have hE : Real.exp (β * (3 / 2)) ≤ Real.exp (L * (3 / 2)) := Real.exp_le_exp.mpr (by linarith)
  have hpow : (D / (D + T₂)) ^ (N - 2) ≤ (2 / L ^ 3) ^ (n - 3) := by
    rw [hN]; exact pow_le_pow_left₀ hq0 hq _
  have hone : 0 ≤ 3 + 2 * D := by linarith
  have hZU := ZU_le hH hT₂.le
  have hZL := ZL_le hH hT₂.le
  have hbd : M * Real.exp (β * (3 / 2)) * (D / (D + T₂)) ^ (N - 2) * (3 + 2 * D) ≤
      2 ^ n * (1 + 2 * s₁) * L ^ 4 * Real.exp (L * (3 / 2)) * (2 / L ^ 3) ^ (n - 3) * (5 * L ^ 3) := by
    have h1 : 0 ≤ (D / (D + T₂)) ^ (N - 2) := pow_nonneg hq0 _
    have h2 : 0 ≤ (2 / L ^ 3) ^ (n - 3) := by positivity
    have h3 : 0 ≤ 2 ^ n * (1 + 2 * s₁) * L ^ 4 := hM0.trans hM
    exact mul_le_mul (mul_le_mul (mul_le_mul hM hE (Real.exp_pos _).le h3) hpow h1 (by positivity)) hD3 hone
      (by positivity)
  have e : (2 : ℝ) ^ n * (2 / L ^ 3) ^ (n - 3) = 8 * (4 / L ^ 3) ^ (n - 3) := by
    obtain ⟨k, hk⟩ : ∃ k, n = k + 3 := ⟨n - 3, by omega⟩
    rw [hk, Nat.add_sub_cancel, pow_add]
    have h4 : (4 / L ^ 3) = 2 * (2 / L ^ 3) := by ring
    rw [h4, mul_pow]
    ring
  have hfin : 2 ^ n * (1 + 2 * s₁) * L ^ 4 * Real.exp (L * (3 / 2)) * (2 / L ^ 3) ^ (n - 3) * (5 * L ^ 3) =
      40 * (1 + 2 * s₁) * L ^ 7 * Real.exp (L * (3 / 2)) * (4 / L ^ 3) ^ (n - 3) := by
    have : 2 ^ n * (1 + 2 * s₁) * L ^ 4 * Real.exp (L * (3 / 2)) * (2 / L ^ 3) ^ (n - 3) * (5 * L ^ 3) =
        ((2 : ℝ) ^ n * (2 / L ^ 3) ^ (n - 3)) * (5 * (1 + 2 * s₁) * L ^ 7 * Real.exp (L * (3 / 2))) := by ring
    rw [this, e]; ring
  rw [hfin] at hbd
  linarith

theorem Bc_le {B C A T₂ L : ℝ} (hBdef : B = C * T₂ ^ A) (hT₂1 : 1 ≤ T₂) (hT₂ : T₂ ≤ 2 * L ^ 6) (hL : 1 ≤ L) :
    Bc B ≤ KB C A * L ^ (6 * kA A) := by
  have hA : A ≤ (kA A : ℝ) := (le_max_left A 0).trans (Nat.le_ceil _)
  have h1 : T₂ ^ A ≤ T₂ ^ ((kA A : ℕ) : ℝ) := Real.rpow_le_rpow_of_exponent_le hT₂1 hA
  rw [Real.rpow_natCast] at h1
  have h2 : T₂ ^ kA A ≤ (2 * L ^ 6) ^ kA A := pow_le_pow_left₀ (by linarith) hT₂ _
  have h3 : (2 * L ^ 6) ^ kA A = 2 ^ kA A * L ^ (6 * kA A) := by rw [mul_pow, ← pow_mul]
  have hTA : 0 ≤ T₂ ^ A := Real.rpow_nonneg (by linarith) A
  have hB : B ≤ |C| * (2 ^ kA A * L ^ (6 * kA A)) := by
    rw [hBdef]
    calc C * T₂ ^ A ≤ |C| * T₂ ^ A := mul_le_mul_of_nonneg_right (le_abs_self C) hTA
      _ ≤ |C| * T₂ ^ kA A := mul_le_mul_of_nonneg_left h1 (abs_nonneg C)
      _ ≤ |C| * (2 * L ^ 6) ^ kA A := mul_le_mul_of_nonneg_left h2 (abs_nonneg C)
      _ = _ := by rw [h3]
  have hpi : π ≤ 4 := Real.pi_le_four
  have hs : Real.sqrt (2 * π) ≤ 3 := by
    have h4 : Real.sqrt (2 * π) ^ 2 = 2 * π := Real.sq_sqrt (by positivity)
    nlinarith [Real.sqrt_nonneg (2 * π)]
  have h60 : 60 * π * Real.sqrt (2 * π) ≤ 720 := by
    have := mul_le_mul hpi hs (Real.sqrt_nonneg _) (by norm_num)
    nlinarith
  have hLp : 1 ≤ L ^ (6 * kA A) := one_le_pow₀ hL
  unfold Bc KB
  nlinarith [abs_nonneg C, pow_pos (by norm_num : (0 : ℝ) < 2) (kA A)]

theorem QK_le {H : ℂ → ℂ} {L T₂ B C A s₁ : ℝ} {n : ℕ} (hL : 2 ≤ L) (h5 : 4 * Real.exp 5 ≤ L ^ 3)
    (hLp : (pE A : ℝ) ≤ L) (hn : L + 3 ≤ ((n - 3 : ℕ) : ℝ)) (hB0 : 0 ≤ B)
    (hBc : Bc B ≤ KB C A * L ^ (6 * kA A)) (h64 : 64 * (T₂ + 3) + 24 ≤ 200 * L ^ 6)
    (hZ : ZU H T₂ + ZL H T₂ ≤ 80 * (1 + 2 * s₁) * L ^ 7 * Real.exp (L * (3 / 2)) * (4 / L ^ 3) ^ (n - 3))
    (hs₁ : s₁ = Real.sqrt (zeroSum + 1)) (hT₂ : 0 ≤ T₂) :
    QK H L T₂ B ≤ Qb C A L := by
  have hL1 : 1 ≤ L := by linarith
  have hpow := pow_four_div_le (m := n - 3) (p := pE A) hL h5 (by linarith) hLp
  have hBc0 := Bc_nonneg hB0
  have hKB : 0 ≤ KB C A := by unfold KB; positivity
  have hX : 0 ≤ 64 * (T₂ + 3) + 24 := by linarith
  have hs0 : 0 ≤ s₁ := by rw [hs₁]; exact Real.sqrt_nonneg _
  have h8 : 8 + Bc B * (64 * (T₂ + 3) + 24) ≤ (8 + 200 * KB C A) * L ^ (6 * kA A + 6) := by
    have h1 : Bc B * (64 * (T₂ + 3) + 24) ≤ KB C A * L ^ (6 * kA A) * (200 * L ^ 6) :=
      mul_le_mul hBc h64 hX (by positivity)
    have h2 : L ^ (6 * kA A + 6) = L ^ (6 * kA A) * L ^ 6 := pow_add _ _ _
    have h3 : 1 ≤ L ^ (6 * kA A + 6) := one_le_pow₀ hL1
    rw [h2] at h3 ⊢
    nlinarith
  have hZ' : ZU H T₂ + ZL H T₂ ≤ 80 * (1 + 2 * s₁) * L ^ 7 * Real.exp (L * (3 / 2)) *
      (Real.exp (-(5 * (L - pE A))) / L ^ pE A) := by
    refine hZ.trans (mul_le_mul_of_nonneg_left hpow (by positivity))
  have hE := Real.exp_pos (L * (3 / 2))
  have hZ0 : 0 ≤ ZU H T₂ + ZL H T₂ := add_nonneg (ZU_nonneg H T₂) (ZL_nonneg H T₂)
  unfold QK
  calc (2 * π)⁻¹ * (Real.exp (L * (3 / 2)) * (8 + Bc B * (64 * (T₂ + 3) + 24))) * (ZU H T₂ + ZL H T₂)
      ≤ (2 * π)⁻¹ * (Real.exp (L * (3 / 2)) * ((8 + 200 * KB C A) * L ^ (6 * kA A + 6))) *
          (80 * (1 + 2 * s₁) * L ^ 7 * Real.exp (L * (3 / 2)) * (Real.exp (-(5 * (L - pE A))) / L ^ pE A)) := by
        gcongr
    _ = Qb C A L := by
        have hp : L ^ (6 * kA A + 6) * L ^ 7 = L ^ pE A := by rw [← pow_add]; rfl
        have hLp0 : L ^ pE A ≠ 0 := pow_ne_zero _ (by linarith)
        have he : Real.exp (L * (3 / 2)) * Real.exp (L * (3 / 2)) * Real.exp (-(5 * (L - pE A))) =
            Real.exp (5 * pE A) * Real.exp (-(2 * L)) := by
          rw [← Real.exp_add, ← Real.exp_add, ← Real.exp_add]; congr 1; ring
        unfold Qb KQ
        rw [← hs₁]
        calc (2 * π)⁻¹ * (Real.exp (L * (3 / 2)) * ((8 + 200 * KB C A) * L ^ (6 * kA A + 6))) *
              (80 * (1 + 2 * s₁) * L ^ 7 * Real.exp (L * (3 / 2)) * (Real.exp (-(5 * (L - pE A))) / L ^ pE A))
            = (2 * π)⁻¹ * (8 + 200 * KB C A) * 80 * (1 + 2 * s₁) * ((L ^ (6 * kA A + 6) * L ^ 7) / L ^ pE A) *
                (Real.exp (L * (3 / 2)) * Real.exp (L * (3 / 2)) * Real.exp (-(5 * (L - pE A)))) := by ring
          _ = _ := by rw [hp, div_self hLp0, he]; ring

end Riemann.Deep.LeanNBRate

