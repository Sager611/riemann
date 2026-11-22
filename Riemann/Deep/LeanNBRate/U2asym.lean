import Riemann.Deep.LeanNBRate.U2aux

/-!
# U2 asymptotics: the parameter conditions for large ℓ = q⁴
-/

open MeasureTheory Filter Topology Set

noncomputable section

namespace Riemann.Deep.LeanNBRate

theorem tendsto_logpow_div (b : ℕ) (C : ℝ) :
    Tendsto (fun q : ℝ => (Real.log q + C) ^ b / q) atTop (𝓝 0) := by
  have h1 : Tendsto (fun q : ℝ => Real.log q ^ b / (1 * q + 0)) atTop (𝓝 0) :=
    Real.tendsto_pow_log_div_mul_add_atTop 1 0 b one_ne_zero
  simp only [one_mul, add_zero] at h1
  have h2 : Tendsto (fun q : ℝ => 2 ^ b * (Real.log q ^ b / q)) atTop (𝓝 0) := by
    simpa using h1.const_mul (2 ^ b)
  refine tendsto_of_tendsto_of_tendsto_of_le_of_le' tendsto_const_nhds h2 ?_ ?_
  · filter_upwards [eventually_ge_atTop (Real.exp |C|), eventually_gt_atTop 0] with q hq hq0
    have : |C| ≤ Real.log q := by rwa [Real.le_log_iff_exp_le hq0]
    have : 0 ≤ Real.log q + C := by have := neg_abs_le C; linarith
    positivity
  · filter_upwards [eventually_ge_atTop (Real.exp |C|), eventually_gt_atTop 0] with q hq hq0
    have hl : |C| ≤ Real.log q := by rwa [Real.le_log_iff_exp_le hq0]
    have h0 : 0 ≤ Real.log q + C := by have := neg_abs_le C; linarith
    have h3 : Real.log q + C ≤ 2 * Real.log q := by have := le_abs_self C; linarith
    have := pow_le_pow_left₀ h0 h3 b
    rw [mul_pow] at this
    rw [← mul_div_assoc]
    exact div_le_div_of_nonneg_right this hq0.le

/-- The multiplicity cutoff for high groups at ℓ = q⁴. -/
def KuOf (ε₁ q : ℝ) : ℕ := ⌊30 * ε₁ * Real.log q⌋₊ + 1

theorem KuOf_le {ε₁ q : ℝ} (hε : 0 ≤ ε₁) (hε30 : 30 * ε₁ ≤ 1) (hq : 1 ≤ q) :
    (KuOf ε₁ q : ℝ) ≤ Real.log q + 1 := by
  unfold KuOf
  push_cast
  have hl := Real.log_nonneg hq
  have h1 : (⌊30 * ε₁ * Real.log q⌋₊ : ℝ) ≤ 30 * ε₁ * Real.log q := Nat.floor_le (by positivity)
  have h2 : 30 * ε₁ * Real.log q ≤ Real.log q := by nlinarith
  linarith

theorem KuOf_ge {ε₁ q T : ℝ} (hε : 0 ≤ ε₁) (hq : 2 ≤ q) (hT0 : 0 ≤ T) (hT : T ≤ q ^ 28) :
    ε₁ * Real.log (T + 2) ≤ KuOf ε₁ q := by
  unfold KuOf
  have hq0 : 0 < q := by linarith
  have h1 : T + 2 ≤ q ^ 30 := by
    have h4 : q ^ 28 * 4 ≤ q ^ 30 := by
      have e : q ^ 30 = q ^ 28 * q ^ 2 := by ring
      rw [e]; gcongr; nlinarith
    have : 1 ≤ q ^ 28 := one_le_pow₀ (by linarith)
    linarith
  have h2 : Real.log (T + 2) ≤ 30 * Real.log q := by
    calc Real.log (T + 2) ≤ Real.log (q ^ 30) := Real.log_le_log (by linarith) h1
      _ = 30 * Real.log q := by rw [Real.log_pow]; norm_num
  have h3 := Nat.lt_floor_add_one (30 * ε₁ * Real.log q)
  push_cast
  nlinarith

theorem pow_KuOf_le {κ ε₁ q : ℝ} (hκ : 1 ≤ κ) (hε : 0 ≤ ε₁) (hsmall : 180 * ε₁ * Real.log (2 * κ) ≤ 1)
    (hq : 1 ≤ q) : (2 * κ) ^ (6 * KuOf ε₁ q) ≤ (2 * κ) ^ 6 * q := by
  unfold KuOf
  have h2κ : 0 < 2 * κ := by linarith
  have hl := Real.log_nonneg (by linarith : (1 : ℝ) ≤ 2 * κ)
  have hlq := Real.log_nonneg hq
  set n := ⌊30 * ε₁ * Real.log q⌋₊
  have e : (2 * κ) ^ (6 * (n + 1)) = (2 * κ) ^ 6 * (2 * κ) ^ (6 * n) := by ring
  rw [e]
  have h6 : 0 < (2 * κ) ^ 6 := by positivity
  apply mul_le_mul_of_nonneg_left _ h6.le
  have hn : (n : ℝ) ≤ 30 * ε₁ * Real.log q := Nat.floor_le (by positivity)
  rw [← Real.exp_log (pow_pos h2κ _), ← Real.exp_log (by linarith : (0 : ℝ) < q)]
  apply Real.exp_le_exp.mpr
  rw [Real.log_pow]
  push_cast
  calc (6 * (n : ℝ)) * Real.log (2 * κ) ≤ 6 * (30 * ε₁ * Real.log q) * Real.log (2 * κ) := by gcongr
    _ = (180 * ε₁ * Real.log (2 * κ)) * Real.log q := by ring
    _ ≤ 1 * Real.log q := by gcongr
    _ = Real.log q := one_mul _

section conds

variable {q Y δ₀ : ℝ} (hq : 1 ≤ q) (hY : 1 ≤ Y) (hX : Y ^ 90 / q ≤ δ₀)
include hq hY hX

theorem Ypow_le {k : ℕ} (hk : k ≤ 90) : Y ^ k ≤ δ₀ * q := by
  have hq0 : 0 < q := by linarith
  have h1 : Y ^ k ≤ Y ^ 90 := pow_le_pow_right₀ hY hk
  have h2 : Y ^ 90 ≤ δ₀ * q := by rwa [div_le_iff₀ hq0] at hX
  linarith

theorem condA {M : ℕ} (hδ : δ₀ ≤ 1 / 16) (hMY : (M : ℝ) + 1 ≤ Y) (hAY : A₀ ≤ Y) :
    4 * A₀ / (1 / q ^ 2) * ((M : ℝ) + 1) ^ 2 ≤ q ^ 4 := by
  have hq0 : 0 < q := by linarith
  have hA0 := A₀_pos
  have e : 4 * A₀ / (1 / q ^ 2) = 4 * A₀ * q ^ 2 := by field_simp
  rw [e]
  have h3 := Ypow_le hq hY hX (k := 3) (by norm_num)
  have hM0 : (0 : ℝ) ≤ (M : ℝ) + 1 := by positivity
  have h1 : ((M : ℝ) + 1) ^ 2 ≤ Y ^ 2 := pow_le_pow_left₀ hM0 hMY 2
  have h2 : 4 * A₀ * q ^ 2 * ((M : ℝ) + 1) ^ 2 ≤ 4 * Y ^ 3 * q ^ 2 := by
    have : 4 * A₀ * q ^ 2 * ((M : ℝ) + 1) ^ 2 ≤ 4 * Y * q ^ 2 * Y ^ 2 := by gcongr
    calc _ ≤ 4 * Y * q ^ 2 * Y ^ 2 := this
      _ = 4 * Y ^ 3 * q ^ 2 := by ring
  have hδ0 : 0 ≤ δ₀ := by
    have : 0 ≤ Y ^ 90 / q := by positivity
    linarith
  calc 4 * A₀ * q ^ 2 * ((M : ℝ) + 1) ^ 2 ≤ 4 * Y ^ 3 * q ^ 2 := h2
    _ ≤ 4 * (δ₀ * q) * q ^ 2 := by gcongr
    _ ≤ 4 * (1 / 16 * q) * q ^ 2 := by gcongr
    _ ≤ q ^ 4 := by nlinarith [pow_le_pow_right₀ hq (by norm_num : 3 ≤ 4)]

theorem condB {M : ℕ} (hδ : δ₀ ≤ 1 / 16) (hMY : (M : ℝ) + 1 ≤ Y) (hAY : A₀ ≤ Y) :
    (M : ℝ) * (4 * A₀ / (1 / q ^ 3 / 4)) * ((M : ℝ) + 1) ^ 2 ≤ q ^ 4 := by
  have hq0 : 0 < q := by linarith
  have hA0 := A₀_pos
  have e : 4 * A₀ / (1 / q ^ 3 / 4) = 16 * A₀ * q ^ 3 := by field_simp; ring
  rw [e]
  have h4 := Ypow_le hq hY hX (k := 4) (by norm_num)
  have hM0 : (0 : ℝ) ≤ (M : ℝ) + 1 := by positivity
  have hMle : (M : ℝ) ≤ Y := by linarith
  have h1 : ((M : ℝ) + 1) ^ 2 ≤ Y ^ 2 := pow_le_pow_left₀ hM0 hMY 2
  have hδ0 : 0 ≤ δ₀ := by
    have : 0 ≤ Y ^ 90 / q := by positivity
    linarith
  calc (M : ℝ) * (16 * A₀ * q ^ 3) * ((M : ℝ) + 1) ^ 2 ≤ Y * (16 * Y * q ^ 3) * Y ^ 2 := by gcongr
    _ = 16 * Y ^ 4 * q ^ 3 := by ring
    _ ≤ 16 * (δ₀ * q) * q ^ 3 := by gcongr
    _ ≤ 16 * (1 / 16 * q) * q ^ 3 := by gcongr
    _ = q ^ 4 := by ring

theorem condD {M : ℕ} (hδ : δ₀ ≤ 1) (hMY : (M : ℝ) + 1 ≤ Y) : (M : ℝ) ≤ q := by
  have h1 := Ypow_le hq hY hX (k := 1) (by norm_num)
  have hq0 : 0 < q := by linarith
  have : δ₀ * q ≤ q := by nlinarith
  rw [pow_one] at h1
  linarith

theorem condE {M : ℕ} {ε₁ : ℝ} (hδ : δ₀ ≤ ε₁ / 4) (hMY : (M : ℝ) + 1 ≤ Y) :
    ((M : ℝ) + 2) * (1 / q ^ 3) ≤ ε₁ / 2 := by
  have h1 := Ypow_le hq hY hX (k := 1) (by norm_num)
  rw [pow_one] at h1
  have hq0 : 0 < q := by linarith
  have hδ0 : 0 ≤ δ₀ := by
    have : 0 ≤ Y ^ 90 / q := by positivity
    linarith
  have h2 : (M : ℝ) + 2 ≤ 2 * (δ₀ * q) := by linarith
  have h3 : 1 / q ^ 3 ≤ 1 / q := by
    apply div_le_div_of_nonneg_left (by norm_num) hq0
    calc q = q ^ 1 := (pow_one q).symm
      _ ≤ q ^ 3 := pow_le_pow_right₀ hq (by norm_num)
  have h4 : 0 ≤ (M : ℝ) + 2 := by positivity
  calc ((M : ℝ) + 2) * (1 / q ^ 3) ≤ (2 * (δ₀ * q)) * (1 / q) := by gcongr
    _ = 2 * δ₀ := by field_simp
    _ ≤ ε₁ / 2 := by linarith

end conds

theorem smallB_le {κ q Y w : ℝ} {M : ℕ} (hκ : 1 ≤ κ) (hq : 1 ≤ q) (hw0 : 0 ≤ w) (hw : w ≤ 1 / q ^ 2)
    (hpow : (2 * κ) ^ (6 * M) ≤ (2 * κ) ^ 6 * q) (hMY : (M : ℝ) + 1 ≤ Y) :
    smallB κ M w (1 / q ^ 2) ≤ 64 * CB ^ 2 * (2 * κ) ^ 6 * Y ^ 6 / q ^ 3 := by
  unfold smallB
  have hq0 : 0 < q := by linarith
  have h2κ : 1 ≤ 2 * κ := by linarith
  have hCB := CB_nonneg
  set η : ℝ := 1 / q ^ 2
  have hη : 0 < η := by positivity
  have hwη : (w + η) ^ 2 / η ≤ 4 * η := by
    rw [div_le_iff₀ hη]
    have : w + η ≤ 2 * η := by linarith
    have h0 : 0 ≤ w + η := by linarith
    nlinarith
  have e1 : (2 * κ) ^ (2 * (M - 1)) ≤ (2 * κ) ^ (2 * M) := pow_le_pow_right₀ h2κ (by omega)
  have e2 : (2 * κ) ^ (M - 1) ≤ (2 * κ) ^ M := pow_le_pow_right₀ h2κ (by omega)
  have hM0 : (0 : ℝ) ≤ (M : ℝ) + 1 := by positivity
  have e3 : ((M : ℝ) + 1) ^ 3 ≤ Y ^ 3 := pow_le_pow_left₀ hM0 hMY 3
  have hY0 : 0 ≤ Y := by linarith
  calc 4 * (2 * κ) ^ (2 * (M - 1)) * ((2 * κ) ^ (M - 1) * CB * ((M : ℝ) + 1) ^ 3 * ((w + η) ^ 2 / η)) ^ 2 *
        ((2 * κ) ^ (M - 1)) ^ 2
      ≤ 4 * (2 * κ) ^ (2 * M) * ((2 * κ) ^ M * CB * Y ^ 3 * (4 * η)) ^ 2 * ((2 * κ) ^ M) ^ 2 := by
        gcongr
    _ = 64 * CB ^ 2 * (2 * κ) ^ (6 * M) * Y ^ 6 * η ^ 2 := by ring
    _ ≤ 64 * CB ^ 2 * ((2 * κ) ^ 6 * q) * Y ^ 6 * η ^ 2 := by gcongr
    _ = 64 * CB ^ 2 * (2 * κ) ^ 6 * Y ^ 6 / q ^ 3 := by simp only [η]; field_simp

theorem bernC_le {κ q Y w : ℝ} {M : ℕ} (hκ : 1 ≤ κ) (hq : 1 ≤ q) (hw0 : 0 ≤ w) (hw : w ≤ 1 / q ^ 2)
    (hpow : (2 * κ) ^ (6 * M) ≤ (2 * κ) ^ 6 * q) (hMY : (M : ℝ) + 1 ≤ Y) :
    q ^ 4 * bernC κ M w (1 / q ^ 2) ≤ 4 * CB * (2 * κ) ^ 6 * q ^ 3 * Y ^ 3 := by
  unfold bernC
  have hq0 : 0 < q := by linarith
  have h2κ : 1 ≤ 2 * κ := by linarith
  have hCB := CB_nonneg
  set η : ℝ := 1 / q ^ 2
  have hη : 0 < η := by positivity
  have hwη : (w + η) ^ 2 / η ≤ 4 * η := by
    rw [div_le_iff₀ hη]
    have : w + η ≤ 2 * η := by linarith
    have h0 : 0 ≤ w + η := by linarith
    nlinarith
  have e2 : (2 * κ) ^ (M - 1) ≤ (2 * κ) ^ (6 * M) := pow_le_pow_right₀ h2κ (by omega)
  have hM0 : (0 : ℝ) ≤ (M : ℝ) + 1 := by positivity
  have e3 : ((M : ℝ) + 1) ^ 3 ≤ Y ^ 3 := pow_le_pow_left₀ hM0 hMY 3
  have hY0 : 0 ≤ Y := by linarith
  calc q ^ 4 * ((2 * κ) ^ (M - 1) * CB * ((M : ℝ) + 1) ^ 3 * ((w + η) ^ 2 / η))
      ≤ q ^ 4 * (((2 * κ) ^ 6 * q) * CB * Y ^ 3 * (4 * η)) := by gcongr; exact e2.trans hpow
    _ = 4 * CB * (2 * κ) ^ 6 * q ^ 3 * Y ^ 3 := by simp only [η]; field_simp

/-- The tail constant at ℓ = q⁴, s₀ = q⁻³, N = 40. -/
def KT : ℝ := (2 * Real.pi)⁻¹ ^ 2 * (32 * ((B₀ * B₀) * (B₀ * B₀)) * 3200 ^ 40) ^ 2 * 4 * 2 ^ 78 *
  (8 * Real.pi)

theorem KT_nonneg : 0 ≤ KT := by unfold KT; have := B₀_pos; positivity

theorem tailC_bound {κ q Y w : ℝ} {M : ℕ} (hκ : 1 ≤ κ) (hq : 1 ≤ q) (hw0 : 0 ≤ w)
    (hw : w + 1 / q ^ 3 / 4 ≤ Y / q ^ 3) (hMY : (M : ℝ) + 1 + 8 * A₀ ≤ Y) (hA : 1 + 16 * A₀ ≤ 2 * Y) :
    tailC κ M w (q ^ 4) (1 / q ^ 3) 40 ≤ KT * (2 * κ) ^ (6 * M) * Y ^ 84 / q ^ 77 := by
  have hq0 : 0 < q := by linarith
  have h2κ : 1 ≤ 2 * κ := by linarith
  have hB := B₀_pos
  have hA0 := A₀_pos
  have hY0 : 0 ≤ Y := by linarith
  set B4 := (B₀ * B₀) * (B₀ * B₀)
  have hB4 : 0 < B4 := by positivity
  have e1 : (2 * κ) ^ (M - 1) ≤ (2 * κ) ^ M := pow_le_pow_right₀ h2κ (Nat.sub_le M 1)
  have e2 : (2 * κ) ^ (4 * (M - 1)) ≤ (2 * κ) ^ (4 * M) := pow_le_pow_right₀ h2κ (by omega)
  have p2 : q ^ 4 + 4 * A₀ / (1 / q ^ 3 / 4) ≤ 2 * Y * q ^ 4 := by
    have e : 4 * A₀ / (1 / q ^ 3 / 4) = 16 * A₀ * q ^ 3 := by field_simp; ring
    rw [e]
    have h3 : q ^ 3 ≤ q ^ 4 := pow_le_pow_right₀ hq (by norm_num)
    have h4 : 16 * A₀ * q ^ 3 ≤ 16 * A₀ * q ^ 4 := mul_le_mul_of_nonneg_left h3 (by positivity)
    have h5 : (1 + 16 * A₀) * q ^ 4 ≤ 2 * Y * q ^ 4 := mul_le_mul_of_nonneg_right hA (by positivity)
    linarith
  have p3 : 2 * ((M : ℝ) + 1 + 8 * A₀) * ((40 : ℕ) : ℝ) ^ 2 / q ^ 4 ≤ 3200 * Y / q ^ 4 := by
    apply div_le_div_of_nonneg_right _ (by positivity)
    push_cast
    linarith
  have p30 : 0 ≤ 2 * ((M : ℝ) + 1 + 8 * A₀) * ((40 : ℕ) : ℝ) ^ 2 / q ^ 4 := by positivity
  have hX : (2 * κ) ^ (M - 1) * (16 * (w + 1 / q ^ 3 / 4) * B4 * q ^ 4 * (q ^ 4 + 4 * A₀ / (1 / q ^ 3 / 4)) *
      (2 * ((M : ℝ) + 1 + 8 * A₀) * ((40 : ℕ) : ℝ) ^ 2 / q ^ 4) ^ 40) ≤
      (2 * κ) ^ M * (16 * (Y / q ^ 3) * B4 * q ^ 4 * (2 * Y * q ^ 4) * (3200 * Y / q ^ 4) ^ 40) := by
    have hp : (2 * ((M : ℝ) + 1 + 8 * A₀) * ((40 : ℕ) : ℝ) ^ 2 / q ^ 4) ^ 40 ≤ (3200 * Y / q ^ 4) ^ 40 :=
      pow_le_pow_left₀ p30 p3 40
    have h0 : 0 ≤ 16 * (w + 1 / q ^ 3 / 4) := by positivity
    have h00 : 0 ≤ q ^ 4 + 4 * A₀ / (1 / q ^ 3 / 4) := by positivity
    gcongr
  have hX0 : 0 ≤ (2 * κ) ^ (M - 1) * (16 * (w + 1 / q ^ 3 / 4) * B4 * q ^ 4 *
      (q ^ 4 + 4 * A₀ / (1 / q ^ 3 / 4)) * (2 * ((M : ℝ) + 1 + 8 * A₀) * ((40 : ℕ) : ℝ) ^ 2 / q ^ 4) ^ 40) := by
    positivity
  unfold tailC
  have e78 : 2 * 40 - 2 = 78 := by norm_num
  rw [e78]
  calc (2 * Real.pi)⁻¹ ^ 2 * ((2 * κ) ^ (M - 1) * (16 * (w + 1 / q ^ 3 / 4) * B4 * q ^ 4 *
        (q ^ 4 + 4 * A₀ / (1 / q ^ 3 / 4)) * (2 * ((M : ℝ) + 1 + 8 * A₀) * ((40 : ℕ) : ℝ) ^ 2 / q ^ 4) ^ 40)) ^ 2 *
        (4 * (2 * κ) ^ (4 * (M - 1)) / q ^ 4) * ((1 / q ^ 3 / 2) ^ 78)⁻¹ * (8 * Real.pi / (1 / q ^ 3))
      ≤ (2 * Real.pi)⁻¹ ^ 2 * ((2 * κ) ^ M * (16 * (Y / q ^ 3) * B4 * q ^ 4 * (2 * Y * q ^ 4) *
        (3200 * Y / q ^ 4) ^ 40)) ^ 2 * (4 * (2 * κ) ^ (4 * M) / q ^ 4) * ((1 / q ^ 3 / 2) ^ 78)⁻¹ *
        (8 * Real.pi / (1 / q ^ 3)) := by
        gcongr
    _ = KT * (2 * κ) ^ (6 * M) * Y ^ 84 / q ^ 77 := by
        unfold KT
        field_simp
        ring

theorem phiE_sq_le {κ q Y w H₀ mk : ℝ} {M : ℕ} (hκ : 1 ≤ κ) (hq : 1 ≤ q) (hw0 : 0 ≤ w) (hw : w ≤ 1 / q ^ 2)
    (hpow : (2 * κ) ^ (6 * M) ≤ (2 * κ) ^ 6 * q) (hMY : (M : ℝ) + 1 ≤ Y) (hH₀ : 0 ≤ H₀) (hmk : 0 ≤ mk) :
    phiE κ H₀ mk (q ^ 4) M w (1 / q ^ 2) ^ 2 ≤
      (6 * (H₀ + mk + 1) * (2 + mk + 4 * CB * (2 * κ) ^ 6 * q ^ 3 * Y ^ 3) / q ^ 2) ^ 2 *
        (2 * κ) ^ (4 * (M - 1)) := by
  unfold phiE
  have hq0 : 0 < q := by linarith
  have h2κ : 0 < 2 * κ := by linarith
  have hsq : Real.sqrt (4 * (2 * κ) ^ (4 * (M - 1)) / q ^ 4) = 2 * (2 * κ) ^ (2 * (M - 1)) / q ^ 2 := by
    rw [show 4 * (2 * κ) ^ (4 * (M - 1)) / q ^ 4 = (2 * (2 * κ) ^ (2 * (M - 1)) / q ^ 2) ^ 2 by
      rw [div_pow, mul_pow, ← pow_mul]; ring_nf]
    exact Real.sqrt_sq (by positivity)
  rw [hsq]
  have hb := bernC_le hκ hq hw0 hw hpow hMY
  have hbern0 : 0 ≤ bernC κ M w (1 / q ^ 2) := by
    unfold bernC; have := CB_nonneg; have : (0 : ℝ) < 1 / q ^ 2 := by positivity
    positivity
  have h1 : 3 * (H₀ + mk + 1) * (2 + mk + q ^ 4 * bernC κ M w (1 / q ^ 2)) * (2 * (2 * κ) ^ (2 * (M - 1)) / q ^ 2) ≤
      3 * (H₀ + mk + 1) * (2 + mk + 4 * CB * (2 * κ) ^ 6 * q ^ 3 * Y ^ 3) * (2 * (2 * κ) ^ (2 * (M - 1)) / q ^ 2) := by
    gcongr
  have h0 : 0 ≤ 3 * (H₀ + mk + 1) * (2 + mk + q ^ 4 * bernC κ M w (1 / q ^ 2)) * (2 * (2 * κ) ^ (2 * (M - 1)) / q ^ 2) := by
    positivity
  calc (3 * (H₀ + mk + 1) * (2 + mk + q ^ 4 * bernC κ M w (1 / q ^ 2)) * (2 * (2 * κ) ^ (2 * (M - 1)) / q ^ 2)) ^ 2
      ≤ (3 * (H₀ + mk + 1) * (2 + mk + 4 * CB * (2 * κ) ^ 6 * q ^ 3 * Y ^ 3) * (2 * (2 * κ) ^ (2 * (M - 1)) / q ^ 2)) ^ 2 :=
        pow_le_pow_left₀ h0 h1 2
    _ = (6 * (H₀ + mk + 1) * (2 + mk + 4 * CB * (2 * κ) ^ 6 * q ^ 3 * Y ^ 3) / q ^ 2) ^ 2 *
        (2 * κ) ^ (4 * (M - 1)) := by
        have hx : (2 * κ) ^ (4 * (M - 1)) = ((2 * κ) ^ (2 * (M - 1))) ^ 2 := by
          rw [← pow_mul]; congr 1; ring
        rw [hx]
        ring


end Riemann.Deep.LeanNBRate

