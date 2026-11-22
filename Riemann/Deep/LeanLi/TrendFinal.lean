module

public import Riemann.Deep.LeanLi.TrendLow

/-!
# LeanLi: the trend estimate (H4') for `t(n) = (n/2)(log n - 1 - log 2π + γ)`
-/

@[expose] public section

open Real Finset MeasureTheory Set

namespace LeanLi

/-- The constant of the trend estimate. -/
noncomputable def Kc (C δ₀ : ℝ) : ℝ :=
  (4 + 2 * π ^ 2) + ((π / 4) * (5 + 5 * π) + (π / 4) * (C * 2 ^ 30) * L1) +
    (3 * π / 2) * (C * 2 ^ 30) * L1 + (3 * π ^ 2 / 2) * (C * 2 ^ 30) * (5 * π + log (π / δ₀) + 1) +
    (π / 2) * (C * 2 ^ 30) * L1 + 2 * C * π ^ 2 * (π / δ₀ + 1 / 2 + 6 * π + log (2 * π / δ₀) + 1) + 8

lemma log_le_self' {x : ℝ} (hx : 0 < x) : log x ≤ x := by
  have := log_le_sub_one_of_pos hx; linarith

lemma tinyB_mul_pow_le {C : ℝ} (hC : 0 ≤ C) {N : ℕ} (hN : 1 ≤ N) {k : ℕ} (hk : k ≤ 14) :
    tinyB C N * (N : ℝ) ^ k ≤ C * 2 ^ 30 := by
  have hN1 : (1 : ℝ) ≤ N := by exact_mod_cast hN
  unfold tinyB
  rw [div_mul_eq_mul_div, div_le_iff₀ (by positivity)]
  have : (N : ℝ) ^ k ≤ (N : ℝ) ^ 14 := pow_le_pow_right₀ hN1 hk
  nlinarith [mul_le_mul_of_nonneg_left this (by positivity : (0 : ℝ) ≤ C * 2 ^ 30)]

set_option maxHeartbeats 4000000 in
theorem err_le {N : ℕ} {w : ℕ → ℝ} {C : ℝ} (hg : PhiGood N w C) (hN : 2 ≤ N) {a δ₀ : ℝ}
    (ha : 1 / √(N : ℝ) ≤ a) (ha' : a ≤ π - δ₀) (hδ₀ : 1 / (2 * √(N : ℝ)) ≤ δ₀ / 2) :
    (3 + |log (2 * π * a)| +
      ((π / 4) * (2 + log π + log (4 * π) + 3 * log N) + π * tinyB C N * (N : ℝ) ^ 6 / 4 * L1) +
      3 * π / 2 * tinyB C N * N * L1 +
      (2 * π * tinyB C N * |a0 a| + M1 a δ₀ * (2 * C / N) +
        tinyB C N * (π * N / 2 * L1 + π * |a0 a|))) + 8 * (1 + log N) ≤
      Kc C δ₀ * (√(N : ℝ) * (1 + log N)) := by
  have hN1 : 1 ≤ N := by omega
  obtain ⟨he0, he1, he2, he3⟩ := eps_facts hN1
  have hπ := pi_gt_three
  have hπ' := pi_pos
  have hNpos : (0 : ℝ) < N := by exact_mod_cast hN1
  have hN1' : (1 : ℝ) ≤ N := by exact_mod_cast hN1
  have hs1 : 1 ≤ √(N : ℝ) := by rw [Real.one_le_sqrt]; exact hN1'
  have hs0 : 0 < √(N : ℝ) := by linarith
  have hlogN : 0 ≤ log N := log_nonneg hN1'
  have ha0 : 0 < a := lt_of_lt_of_le (by positivity) ha
  have hδ0 : 0 < δ₀ := by linarith
  have hC := hg.C_nonneg hN1
  have hL1 := L1_nonneg
  set X := √(N : ℝ) * (1 + log N) with hX
  have hX1 : 1 ≤ X := by nlinarith
  have hXl : 1 + log N ≤ X := by nlinarith
  have hXs : √(N : ℝ) ≤ X := by nlinarith
  have hXsl : √(N : ℝ) * log N ≤ X := by nlinarith
  -- |log(2πa)| ≤ 2π² + log N
  have hlog2 : |log (2 * π * a)| ≤ 2 * π ^ 2 + log N := by
    rcases le_or_gt 1 (2 * π * a) with h | h
    · rw [abs_of_nonneg (log_nonneg h)]
      have := log_le_self' (by positivity : 0 < 2 * π * a)
      nlinarith
    · rw [abs_of_neg (log_neg (by positivity) h)]
      have h1 : 1 / (2 * π * a) ≤ √(N : ℝ) := by
        rw [div_le_iff₀ (by positivity)]
        rw [div_le_iff₀ hs0] at ha
        nlinarith
      have h2 : log (1 / (2 * π * a)) ≤ log √(N : ℝ) := log_le_log (by positivity) h1
      rw [one_div, log_inv, Real.log_sqrt hNpos.le] at h2
      nlinarith
  -- tiny bounds
  have t6 := tinyB_mul_pow_le hC hN1 (k := 6) (by norm_num)
  have t1 : tinyB C N * N ≤ C * 2 ^ 30 := by simpa using tinyB_mul_pow_le hC hN1 (k := 1) (by norm_num)
  have hτ := tinyB_nonneg hg hN1
  -- |ℓ(a)| and |a₀(a)|
  have hell : |ell a| ≤ 5 * π + log (π / δ₀) + log N := by
    have h := abs_ell_le ha0 (by linarith)
    have l1 : log (π / a) ≤ π + log N := by
      have : π / a ≤ π * √(N : ℝ) := by
        rw [div_le_iff₀ ha0]; rw [div_le_iff₀ hs0] at ha; nlinarith
      have h2 : log (π / a) ≤ log (π * √(N : ℝ)) := log_le_log (by positivity) this
      rw [log_mul hπ'.ne' hs0.ne', Real.log_sqrt hNpos.le] at h2
      have := log_le_self' hπ'
      linarith
    have l2 : log (π / (π - a)) ≤ log (π / δ₀) :=
      log_le_log (div_pos hπ' (by linarith)) (div_le_div_of_nonneg_left hπ'.le hδ0 (by linarith))
    have l3 : log (4 * π) ≤ 4 * π := log_le_self' (by positivity)
    linarith
  have ha0N := abs_a0_le_N hN1 (v := a) (by linarith) (by linarith)
  have hlogd : 0 ≤ log (π / δ₀) := log_nonneg (by rw [le_div_iff₀ hδ0]; linarith)
  have hlogd2 : 0 ≤ log (2 * π / δ₀) := log_nonneg (by rw [le_div_iff₀ hδ0]; linarith)
  have hA : tinyB C N * |a0 a| ≤ C * 2 ^ 30 * (π / 2) * (5 * π + log (π / δ₀) + log N) := by
    calc tinyB C N * |a0 a| ≤ tinyB C N * (π * N / 2 * |ell a|) := mul_le_mul_of_nonneg_left ha0N hτ
      _ = (tinyB C N * N) * (π / 2) * |ell a| := by ring
      _ ≤ (C * 2 ^ 30) * (π / 2) * (5 * π + log (π / δ₀) + log N) := by
          apply mul_le_mul (mul_le_mul_of_nonneg_right t1 (by positivity)) hell (abs_nonneg _)
          positivity
  -- M₁ term
  have hM : M1 a δ₀ * (2 * C / N) ≤ 2 * C * π ^ 2 * (π / δ₀ + 1 / 2 + 6 * π + log (2 * π / δ₀) + log N) * √(N : ℝ) := by
    have l1 : log (2 * π / a) ≤ 2 * π + log N := by
      have : 2 * π / a ≤ 2 * π * √(N : ℝ) := by
        rw [div_le_iff₀ ha0]; rw [div_le_iff₀ hs0] at ha; nlinarith
      have h2 : log (2 * π / a) ≤ log (2 * π * √(N : ℝ)) := log_le_log (by positivity) this
      rw [log_mul (by positivity) hs0.ne', Real.log_sqrt hNpos.le] at h2
      have := log_le_self' (by positivity : (0:ℝ) < 2 * π)
      linarith
    have l3 : log (4 * π) ≤ 4 * π := log_le_self' (by positivity)
    have hnum : π / δ₀ + 1 / 2 + (log (2 * π / a) + log (2 * π / δ₀) + log (4 * π)) ≤
        π / δ₀ + 1 / 2 + 6 * π + log (2 * π / δ₀) + log N := by linarith
    have hnum0 : 0 ≤ π / δ₀ + 1 / 2 + (log (2 * π / a) + log (2 * π / δ₀) + log (4 * π)) := by
      have : 0 ≤ log (2 * π / a) := log_nonneg (by rw [le_div_iff₀ ha0]; linarith)
      have : 0 ≤ log (4 * π) := log_nonneg (by linarith)
      positivity
    have hden : 8 * π * (a / (2 * π)) ^ 3 = a ^ 3 / π ^ 2 := by field_simp; ring
    have ha3 : 1 / a ^ 3 ≤ (N : ℝ) * √(N : ℝ) := by
      have hsa : 1 ≤ a * √(N : ℝ) := by rw [div_le_iff₀ hs0] at ha; linarith
      have : 1 ≤ (a * √(N : ℝ)) ^ 3 := one_le_pow₀ hsa
      rw [div_le_iff₀ (by positivity)]
      have e : (a * √(N : ℝ)) ^ 3 = (N : ℝ) * √(N : ℝ) * a ^ 3 := by
        rw [mul_pow, show (3 : ℕ) = 2 + 1 by norm_num, pow_succ (√(N : ℝ)), Real.sq_sqrt hNpos.le]; ring
      linarith
    unfold M1
    rw [hden]
    calc (π / δ₀ + 1 / 2 + (log (2 * π / a) + log (2 * π / δ₀) + log (4 * π))) / (a ^ 3 / π ^ 2) *
          (2 * C / N)
        = (π / δ₀ + 1 / 2 + (log (2 * π / a) + log (2 * π / δ₀) + log (4 * π))) * π ^ 2 * (1 / a ^ 3) *
          (2 * C / N) := by field_simp
      _ ≤ (π / δ₀ + 1 / 2 + 6 * π + log (2 * π / δ₀) + log N) * π ^ 2 * ((N : ℝ) * √(N : ℝ)) *
          (2 * C / N) := by
          apply mul_le_mul_of_nonneg_right _ (by positivity)
          apply mul_le_mul (mul_le_mul_of_nonneg_right hnum (by positivity)) ha3 (by positivity)
          positivity
      _ = 2 * C * π ^ 2 * (π / δ₀ + 1 / 2 + 6 * π + log (2 * π / δ₀) + log N) * √(N : ℝ) := by
          field_simp
  -- assemble
  have hlp : log π ≤ π := log_le_self' hπ'
  have hl4 : log (4 * π) ≤ 4 * π := log_le_self' (by positivity)
  have hlX : log N ≤ X := by linarith
  have k1 : 3 + |log (2 * π * a)| ≤ (4 + 2 * π ^ 2) * X := by
    have : 2 * π ^ 2 ≤ 2 * π ^ 2 * X := le_mul_of_one_le_right (by positivity) hX1
    have e : (4 + 2 * π ^ 2) * X = 4 * X + 2 * π ^ 2 * X := by ring
    rw [e]; linarith
  have k2 : (π / 4) * (2 + log π + log (4 * π) + 3 * log N) + π * tinyB C N * (N : ℝ) ^ 6 / 4 * L1 ≤
      ((π / 4) * (5 + 5 * π) + (π / 4) * (C * 2 ^ 30) * L1) * X := by
    have e1 : (π / 4) * (2 + log π + log (4 * π) + 3 * log N) ≤ (π / 4) * (5 + 5 * π) * X := by
      have : 2 + log π + log (4 * π) + 3 * log N ≤ (5 + 5 * π) * X := by
        have h5 : 5 * π ≤ 5 * π * X := le_mul_of_one_le_right (by positivity) hX1
        have e : (5 + 5 * π) * X = 5 * X + 5 * π * X := by ring
        rw [e]; linarith
      have := mul_le_mul_of_nonneg_left this (by positivity : (0 : ℝ) ≤ π / 4)
      linarith
    have e2 : π * tinyB C N * (N : ℝ) ^ 6 / 4 * L1 ≤ (π / 4) * (C * 2 ^ 30) * L1 * X := by
      have : π * tinyB C N * (N : ℝ) ^ 6 / 4 * L1 = (π / 4) * (tinyB C N * (N : ℝ) ^ 6) * L1 := by ring
      rw [this]
      have h0 : 0 ≤ (π / 4) * (C * 2 ^ 30) * L1 := by positivity
      calc (π / 4) * (tinyB C N * (N : ℝ) ^ 6) * L1 ≤ (π / 4) * (C * 2 ^ 30) * L1 := by
            gcongr
        _ ≤ (π / 4) * (C * 2 ^ 30) * L1 * X := le_mul_of_one_le_right h0 hX1
    linarith
  have k3 : 3 * π / 2 * tinyB C N * N * L1 ≤ (3 * π / 2) * (C * 2 ^ 30) * L1 * X := by
    have h0 : 0 ≤ (3 * π / 2) * (C * 2 ^ 30) * L1 := by positivity
    calc 3 * π / 2 * tinyB C N * N * L1 = (3 * π / 2) * (tinyB C N * N) * L1 := by ring
      _ ≤ (3 * π / 2) * (C * 2 ^ 30) * L1 := by gcongr
      _ ≤ (3 * π / 2) * (C * 2 ^ 30) * L1 * X := le_mul_of_one_le_right h0 hX1
  have k4 : 2 * π * tinyB C N * |a0 a| + tinyB C N * (π * |a0 a|) ≤
      (3 * π ^ 2 / 2) * (C * 2 ^ 30) * (5 * π + log (π / δ₀) + 1) * X := by
    have e : 2 * π * tinyB C N * |a0 a| + tinyB C N * (π * |a0 a|) = 3 * π * (tinyB C N * |a0 a|) := by ring
    rw [e]
    have : 5 * π + log (π / δ₀) + log N ≤ (5 * π + log (π / δ₀) + 1) * X := by
      have h5 : 5 * π + log (π / δ₀) ≤ (5 * π + log (π / δ₀)) * X :=
        le_mul_of_one_le_right (by positivity) hX1
      have e : (5 * π + log (π / δ₀) + 1) * X = (5 * π + log (π / δ₀)) * X + X := by ring
      rw [e]; linarith
    calc 3 * π * (tinyB C N * |a0 a|) ≤ 3 * π * (C * 2 ^ 30 * (π / 2) * (5 * π + log (π / δ₀) + log N)) :=
          mul_le_mul_of_nonneg_left hA (by positivity)
      _ = (3 * π ^ 2 / 2) * (C * 2 ^ 30) * (5 * π + log (π / δ₀) + log N) := by ring
      _ ≤ (3 * π ^ 2 / 2) * (C * 2 ^ 30) * ((5 * π + log (π / δ₀) + 1) * X) :=
          mul_le_mul_of_nonneg_left this (by positivity)
      _ = _ := by ring
  have k5 : tinyB C N * (π * N / 2 * L1) ≤ (π / 2) * (C * 2 ^ 30) * L1 * X := by
    have h0 : 0 ≤ (π / 2) * (C * 2 ^ 30) * L1 := by positivity
    calc tinyB C N * (π * N / 2 * L1) = (π / 2) * (tinyB C N * N) * L1 := by ring
      _ ≤ (π / 2) * (C * 2 ^ 30) * L1 := by gcongr
      _ ≤ (π / 2) * (C * 2 ^ 30) * L1 * X := le_mul_of_one_le_right h0 hX1
  have k6 : M1 a δ₀ * (2 * C / N) ≤ 2 * C * π ^ 2 * (π / δ₀ + 1 / 2 + 6 * π + log (2 * π / δ₀) + 1) * X := by
    refine hM.trans ?_
    have hc0 : 0 ≤ π / δ₀ + 1 / 2 + 6 * π + log (2 * π / δ₀) := by positivity
    have : (π / δ₀ + 1 / 2 + 6 * π + log (2 * π / δ₀) + log N) * √(N : ℝ) ≤
        (π / δ₀ + 1 / 2 + 6 * π + log (2 * π / δ₀) + 1) * X := by
      set c := π / δ₀ + 1 / 2 + 6 * π + log (2 * π / δ₀)
      have e : (c + 1) * X = c * √(N : ℝ) + c * (√(N : ℝ) * log N) + √(N : ℝ) + √(N : ℝ) * log N := by
        rw [hX]; ring
      have e2 : (c + log N) * √(N : ℝ) = c * √(N : ℝ) + √(N : ℝ) * log N := by ring
      rw [e, e2]
      have : 0 ≤ c * (√(N : ℝ) * log N) := by positivity
      linarith
    calc 2 * C * π ^ 2 * (π / δ₀ + 1 / 2 + 6 * π + log (2 * π / δ₀) + log N) * √(N : ℝ) =
          2 * C * π ^ 2 * ((π / δ₀ + 1 / 2 + 6 * π + log (2 * π / δ₀) + log N) * √(N : ℝ)) := by ring
      _ ≤ 2 * C * π ^ 2 * ((π / δ₀ + 1 / 2 + 6 * π + log (2 * π / δ₀) + 1) * X) :=
          mul_le_mul_of_nonneg_left this (by positivity)
      _ = _ := by ring
  have k7 : 8 * (1 + log N) ≤ 8 * X := by linarith
  unfold Kc
  linarith [k1, k2, k3, k4, k5, k6, k7]

/-- **(H4') for `t`**: `τ_t(a) ≥ log(2πa)/a² - K √N (1 + log N)`. -/
theorem trendT_lower {N : ℕ} {w : ℕ → ℝ} {C : ℝ} (hg : PhiGood N w C) (hN : 2 ≤ N) {a δ₀ : ℝ}
    (ha : 1 / √(N : ℝ) ≤ a) (ha' : a ≤ π - δ₀) (hδ₀ : 1 / (2 * √(N : ℝ)) ≤ δ₀ / 2) :
    log (2 * π * a) / a ^ 2 - Kc C δ₀ * (√(N : ℝ) * (1 + log N)) ≤ LR.trend N w trendT a := by
  have h1 := trendI_lower hg hN ha ha' hδ₀
  have h2 := trend_sub_le N w hg.nonneg hg.le_one a
  have h3 := err_le hg hN ha ha' hδ₀
  have h2' := (abs_le.1 h2).1
  linarith

end LeanLi
