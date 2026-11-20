/-
E1PLean.Det — the detector threshold (REPORT step 5; FULLPROOF §8: θ_det = 0.494636567599178 ± 7·10⁻¹⁷, k*(θ_det) = κ):
k*(0.4946366) < κ₂ < k*(0.4946365), so h(θ; κ₂) > 0 for θ ∈ [0.4946366, 1] and h(θ; κ₂) ≤ 0 for θ ∈ (0, 0.4946365];
both ends lie below 371/750 = 0.494666…, the moment edge (the moment edge binds).
-/
import Riemann.Deep.B3Cert.Lean.DetTrig
import Riemann.Deep.B3Cert.Lean.Arch

namespace E1PLean
open E1Lean

variable (F : Fns) (I : Inputs)

/-- c(0.4946366) ≥ −0.185884824698. -/
theorem cW_detp_ge : (-185884824698 / 1000000000000 : ℝ) ≤ cW F (4946366 / 10000000) := by
  obtain ⟨hs, hss, hs1⟩ := sqrt2_facts F
  obtain ⟨hsl, hsh⟩ := sqrt2_bounds F
  have ht := trig_det_lo F
  have hcot : cot F (349760894081 / 1000000000000) ≤ (1370773535391953308165304452301 / 500000000000000000000000000000 : ℝ) := by
    unfold cot
    calc F.cos (349760894081 / 1000000000000) / F.sin (349760894081 / 1000000000000) ≤ (469727337451832017944062175461 / 500000000000000000000000000000 : ℝ) / (171336594019346940584443320563 / 500000000000000000000000000000 : ℝ) :=
          div_le_div₀ (by norm_num) ht.2 (by norm_num) ht.1
      _ ≤ (1370773535391953308165304452301 / 500000000000000000000000000000 : ℝ) := by norm_num
  set x₀ := (4946366 / 10000000 : ℝ) / F.sqrt 2 with hx₀
  have hxlo : (349760894081 / 1000000000000 : ℝ) ≤ x₀ := by rw [hx₀, le_div_iff₀ hs]; linarith
  have hx1 : x₀ ≤ 1 := by rw [hx₀, div_le_iff₀ hs]; linarith
  have hcotx : cot F x₀ ≤ (1370773535391953308165304452301 / 500000000000000000000000000000 : ℝ) := by
    rcases hxlo.lt_or_eq with hlt | heq
    · have := cot_sub_gt F (by norm_num) hlt hx1; linarith
    · rw [← heq]; exact hcot
  unfold cW
  rw [← hx₀]
  have h1 : cot F x₀ / F.sqrt 2 ≤ (1370773535391953308165304452301 / 500000000000000000000000000000 : ℝ) / F.sqrt 2 := div_le_div_of_nonneg_right hcotx hs.le
  have h2 : (1370773535391953308165304452301 / 500000000000000000000000000000 : ℝ) / F.sqrt 2 ≤ (1370773535391953308165304452301 / 500000000000000000000000000000 : ℝ) / (141421356237 / 100000000000) :=
    div_le_div_of_nonneg_left (by norm_num) (by norm_num) hsl
  have h3 : (-185884824698 / 1000000000000 : ℝ) ≤ 2 - 4946366 / 10000000 / 2 - (1370773535391953308165304452301 / 500000000000000000000000000000 : ℝ) / (141421356237 / 100000000000) := by norm_num
  linarith

/-- c(0.4946365) ≤ −0.185885200141. -/
theorem cW_detm_le : cW F (4946365 / 10000000) ≤ (-185885200141 / 1000000000000 : ℝ) := by
  obtain ⟨hs, hss, hs1⟩ := sqrt2_facts F
  obtain ⟨hsl, hsh⟩ := sqrt2_bounds F
  obtain ⟨hc, hsn, hsp⟩ := trig_det_hi F
  have hcot : (685386918118050252079830259917 / 250000000000000000000000000000 : ℝ) ≤ cot F (174880411687 / 500000000000) := by
    unfold cot
    calc (685386918118050252079830259917 / 250000000000000000000000000000 : ℝ) ≤ (234863674778899578598089897037 / 250000000000000000000000000000 : ℝ) / (171336560831794961546611225969 / 500000000000000000000000000000 : ℝ) := by norm_num
      _ ≤ F.cos (174880411687 / 500000000000) / F.sin (174880411687 / 500000000000) := div_le_div₀ (by linarith) hc hsp hsn
  set x₀ := (4946365 / 10000000 : ℝ) / F.sqrt 2 with hx₀
  have hx0 : 0 < x₀ := div_pos (by norm_num) hs
  have hxhi : x₀ ≤ (174880411687 / 500000000000 : ℝ) := by rw [hx₀, div_le_iff₀ hs]; linarith
  have hcotx : (685386918118050252079830259917 / 250000000000000000000000000000 : ℝ) ≤ cot F x₀ := by
    rcases hxhi.lt_or_eq with hlt | heq
    · have := cot_sub_gt F hx0 hlt (by norm_num); linarith
    · rw [heq]; exact hcot
  unfold cW
  rw [← hx₀]
  have h1 : (685386918118050252079830259917 / 250000000000000000000000000000 : ℝ) / F.sqrt 2 ≤ cot F x₀ / F.sqrt 2 := div_le_div_of_nonneg_right hcotx hs.le
  have h2 : (685386918118050252079830259917 / 250000000000000000000000000000 : ℝ) / (141421356238 / 100000000000) ≤ (685386918118050252079830259917 / 250000000000000000000000000000 : ℝ) / F.sqrt 2 :=
    div_le_div_of_nonneg_left (by norm_num) hs hsh
  have h3 : 2 - 4946365 / 10000000 / 2 - (685386918118050252079830259917 / 250000000000000000000000000000 : ℝ) / (141421356238 / 100000000000) ≤ (-185885200141 / 1000000000000 : ℝ) := by norm_num
  linarith

/-- k* is antitone on (0, 1] (c increasing: cW_lt; c < 1: cW_lt_one). -/
theorem kstar_anti {θ₁ θ₂ : ℝ} (h1 : 0 < θ₁) (h12 : θ₁ ≤ θ₂) (h2 : θ₂ ≤ 1) : kstar F θ₂ ≤ kstar F θ₁ := by
  have hc2 := cW_lt_one F (lt_of_lt_of_le h1 h12) h2
  have hcc : cW F θ₁ ≤ cW F θ₂ := by
    rcases h12.lt_or_eq with hlt | heq
    · exact (cW_lt F h1 hlt h2).le
    · rw [heq]
  unfold kstar
  have := div_le_div_of_nonneg_left (by norm_num : (0 : ℝ) ≤ 2) (by linarith) (by linarith : 2 - cW F θ₂ ≤ 2 - cW F θ₁)
  linarith

/-- k*(0.4946366) ≤ 0.0850387096 and k*(0.4946365) ≥ 0.0850388666 (Arb at θ_det: k* = κ = 0.08503876039). -/
theorem kstar_det_bounds : kstar F (4946366 / 10000000) ≤ 850387096 / 10000000000 ∧
    (850388666 / 10000000000 : ℝ) ≤ kstar F (4946365 / 10000000) := by
  have hp := cW_detp_ge F
  have hm := cW_detm_le F
  have hc1 := cW_lt_one F (by norm_num : (0 : ℝ) < 4946366 / 10000000) (by norm_num)
  have ks1 : kstar F (4946366 / 10000000) ≤ 1 - 2 / (2 - (-185884824698 / 1000000000000 : ℝ)) := by
    unfold kstar
    have := div_le_div_of_nonneg_left (by norm_num : (0 : ℝ) ≤ 2) (by linarith)
      (by linarith : 2 - cW F (4946366 / 10000000) ≤ 2 - (-185884824698 / 1000000000000 : ℝ))
    linarith
  have ks2 : 1 - 2 / (2 - (-185885200141 / 1000000000000 : ℝ)) ≤ kstar F (4946365 / 10000000) := by
    unfold kstar
    have := div_le_div_of_nonneg_left (by norm_num : (0 : ℝ) ≤ 2) (by norm_num)
      (by linarith : 2 - (-185885200141 / 1000000000000 : ℝ) ≤ 2 - cW F (4946365 / 10000000))
    linarith
  have g1 : 1 - 2 / (2 - (-185884824698 / 1000000000000 : ℝ)) ≤ (850387096 / 10000000000 : ℝ) := by norm_num
  have g2 : (850388666 / 10000000000 : ℝ) ≤ 1 - 2 / (2 - (-185885200141 / 1000000000000 : ℝ)) := by norm_num
  exact ⟨ks1.trans g1, g2.trans ks2⟩

/-- **Detector threshold** (from H4′): k*(0.4946366) < κ₂ < k*(0.4946365), i.e. θ_det ∈ (0.4946365, 0.4946366):
  h(θ; κ₂) > 0 on [0.4946366, 1], h(θ; κ₂) ≤ 0 on (0, 0.4946365]; 0.4946366 < 371/750. -/
theorem detector_threshold (h4 : H4P I) :
    kstar F (4946366 / 10000000) < κ₂ F I ∧ κ₂ F I < kstar F (4946365 / 10000000) ∧
    (∀ θ : ℝ, 4946366 / 10000000 ≤ θ → θ ≤ 1 → 0 < h F θ (κ₂ F I)) ∧
    (∀ θ : ℝ, 0 < θ → θ ≤ 4946365 / 10000000 → h F θ (κ₂ F I) ≤ 0) ∧ (4946366 / 10000000 : ℝ) < 371 / 750 := by
  obtain ⟨hk, hk1⟩ := kappa2_bounds F h4.1 h4.2
  have hku := kappa2_upper F h4.1
  obtain ⟨kp, km⟩ := kstar_det_bounds F
  have e1 : kstar F (4946366 / 10000000) < κ₂ F I := by
    unfold κ₂; linarith
  have e2 : κ₂ F I < kstar F (4946365 / 10000000) := by
    unfold κ₂; linarith
  refine ⟨e1, e2, ?_, ?_, by norm_num⟩
  · intro θ hθ1 hθ2
    have ha := kstar_anti F (by norm_num) hθ1 hθ2
    exact (h_pos_iff F (by linarith) hθ2 hk1).2 (lt_of_le_of_lt ha e1)
  · intro θ hθ0 hθ1
    have ha := kstar_anti F hθ0 hθ1 (by norm_num)
    by_contra hpos
    rw [not_le] at hpos
    have := (h_pos_iff F hθ0 (by linarith) hk1).1 hpos
    linarith

end E1PLean
