/-
E1PLean.Num — E1′ constants (B3Cert REPORT 'Theorem', FULLPROOF §8; Arb values b3c_arb.log):
c(371/750) ∈ [−0.1857718405, −0.18577184] (Arb −0.18577184013932); κ = 1 − (44/5)⁻¹ log c ∈ [0.0850387602, 0.0850387611]
for c in the H4′ interval (Arb 0.08503876039177); k*(371/750) ∈ [0.0849914142, 0.0849914146] (Arb 0.08499141434976);
h(371/750; κ) ∈ (0.000033545, 0.000033547] (Arb 3.354573606·10⁻⁵).
-/
import Riemann.Deep.B3Cert.Lean.Trig
import Riemann.Deep.B3Cert.Lean.Exp
import Riemann.Deep.SantiAudit.Lean.Bounds2

namespace E1PLean
open E1Lean

variable (F : Fns)

/-- c(371/750) ≥ −0.1857718405. -/
theorem cW_edge_ge : (-1857718405 / 10000000000 : ℝ) ≤ cW F (371 / 750) := by
  obtain ⟨hs, hss, hs1⟩ := sqrt2_facts F
  obtain ⟨hsl, hsh⟩ := sqrt2_bounds F
  have ht := trig_xlo2 F
  have hcot : cot F (13991286177 / 40000000000) ≤ (685341506634698256223664973907 / 250000000000000000000000000000 : ℝ) := by
    unfold cot
    calc F.cos (13991286177 / 40000000000) / F.sin (13991286177 / 40000000000) ≤ (939447389341493878758616448191 / 1000000000000000000000000000000 : ℝ) / (2677290321021843573842093297 / 7812500000000000000000000000 : ℝ) :=
          div_le_div₀ (by norm_num) ht.2 (by norm_num) ht.1
      _ ≤ (685341506634698256223664973907 / 250000000000000000000000000000 : ℝ) := by norm_num
  set x₀ := (371 / 750 : ℝ) / F.sqrt 2 with hx₀
  have hxlo : (13991286177 / 40000000000 : ℝ) ≤ x₀ := by rw [hx₀, le_div_iff₀ hs]; linarith
  have hx1 : x₀ ≤ 1 := by rw [hx₀, div_le_iff₀ hs]; linarith
  have hcotx : cot F x₀ ≤ (685341506634698256223664973907 / 250000000000000000000000000000 : ℝ) := by
    rcases hxlo.lt_or_eq with hlt | heq
    · have := cot_sub_gt F (by norm_num) hlt hx1; linarith
    · rw [← heq]; exact hcot
  unfold cW
  rw [← hx₀]
  have h1 : cot F x₀ / F.sqrt 2 ≤ (685341506634698256223664973907 / 250000000000000000000000000000 : ℝ) / F.sqrt 2 := div_le_div_of_nonneg_right hcotx hs.le
  have h2 : (685341506634698256223664973907 / 250000000000000000000000000000 : ℝ) / F.sqrt 2 ≤ (685341506634698256223664973907 / 250000000000000000000000000000 : ℝ) / (141421356237 / 100000000000) :=
    div_le_div_of_nonneg_left (by norm_num) (by norm_num) hsl
  have h3 : (-1857718405 / 10000000000 : ℝ) ≤ 2 - 371 / 750 / 2 - (685341506634698256223664973907 / 250000000000000000000000000000 : ℝ) / (141421356237 / 100000000000) := by norm_num
  linarith

/-- c(371/750) ≤ −0.18577184. -/
theorem cW_edge_le : cW F (371 / 750) ≤ (-18577184 / 100000000 : ℝ) := by
  obtain ⟨hs, hss, hs1⟩ := sqrt2_facts F
  obtain ⟨hsl, hsh⟩ := sqrt2_bounds F
  obtain ⟨hc, hsn, hsp⟩ := trig_xhi2 F
  have hcot : (685341506513710769096219558427 / 250000000000000000000000000000 : ℝ) ≤ cot F (87445538607 / 250000000000) := by
    unfold cot
    calc (685341506513710769096219558427 / 250000000000000000000000000000 : ℝ) ≤ (18788947786460058914970169321 / 20000000000000000000000000000 : ℝ) / (13707726445781940956107589929 / 40000000000000000000000000000 : ℝ) := by norm_num
      _ ≤ F.cos (87445538607 / 250000000000) / F.sin (87445538607 / 250000000000) := div_le_div₀ (by linarith) hc hsp hsn
  set x₀ := (371 / 750 : ℝ) / F.sqrt 2 with hx₀
  have hx0 : 0 < x₀ := div_pos (by norm_num) hs
  have hxhi : x₀ ≤ (87445538607 / 250000000000 : ℝ) := by rw [hx₀, div_le_iff₀ hs]; linarith
  have hcotx : (685341506513710769096219558427 / 250000000000000000000000000000 : ℝ) ≤ cot F x₀ := by
    rcases hxhi.lt_or_eq with hlt | heq
    · have := cot_sub_gt F hx0 hlt (by norm_num); linarith
    · rw [heq]; exact hcot
  unfold cW
  rw [← hx₀]
  have h1 : (685341506513710769096219558427 / 250000000000000000000000000000 : ℝ) / F.sqrt 2 ≤ cot F x₀ / F.sqrt 2 := div_le_div_of_nonneg_right hcotx hs.le
  have h2 : (685341506513710769096219558427 / 250000000000000000000000000000 : ℝ) / (141421356238 / 100000000000) ≤ (685341506513710769096219558427 / 250000000000000000000000000000 : ℝ) / F.sqrt 2 :=
    div_le_div_of_nonneg_left (by norm_num) hs hsh
  have h3 : 2 - 371 / 750 / 2 - (685341506513710769096219558427 / 250000000000000000000000000000 : ℝ) / (141421356238 / 100000000000) ≤ (-18577184 / 100000000 : ℝ) := by norm_num
  linarith

/-- For c in the H4′ interval: 0.0850387602 ≤ κ = 1 − (44/5)⁻¹ log c < 1. -/
theorem kappa2_bounds {c : ℝ} (hlo : (3138997965019943478277 / 1000000000000000000 : ℝ) ≤ c) (hhi : c ≤ (3138997965019943478278 / 1000000000000000000 : ℝ)) :
    (850387602 / 10000000000 : ℝ) ≤ kappa F c (44 / 5) ∧ kappa F c (44 / 5) < 1 := by
  have hc0 : 0 < c := lt_of_lt_of_le (by norm_num) hlo
  have hlog : F.log c ≤ 8 + 258294547 / 5000000000 := (F.log_le_iff_le_exp c _ hc0).2 (hhi.trans (exp_Lhi2 F))
  have hlog0 : 0 < F.log c := by
    by_contra hneg
    rw [not_lt] at hneg
    have h1 := (F.log_le_iff_le_exp c 0 hc0).1 hneg
    rw [F.exp_zero] at h1
    linarith
  unfold kappa
  constructor
  · have h1 : (44 / 5 : ℝ)⁻¹ * F.log c ≤ (44 / 5 : ℝ)⁻¹ * (8 + 258294547 / 5000000000) :=
      mul_le_mul_of_nonneg_left hlog (by norm_num)
    have h2 : (850387602 / 10000000000 : ℝ) ≤ 1 - (44 / 5 : ℝ)⁻¹ * (8 + 258294547 / 5000000000) := by norm_num
    linarith
  · have h1 : 0 < (44 / 5 : ℝ)⁻¹ * F.log c := mul_pos (by norm_num) hlog0
    linarith

/-- For c ≥ c_lo: κ = 1 − (44/5)⁻¹ log c ≤ 0.0850387611. -/
theorem kappa2_upper {c : ℝ} (hlo : (3138997965019943478277 / 1000000000000000000 : ℝ) ≤ c) : kappa F c (44 / 5) ≤ 850387611 / 10000000000 := by
  have hc0 : 0 < c := lt_of_lt_of_le (by norm_num) hlo
  have hl := le_log_of_exp_le F hc0 ((exp_Llo2 F).trans hlo)
  unfold kappa
  have h1 : (44 / 5 : ℝ)⁻¹ * (8 + 51658903 / 1000000000) ≤ (44 / 5 : ℝ)⁻¹ * F.log c :=
    mul_le_mul_of_nonneg_left hl (by norm_num)
  have h2 : 1 - (44 / 5 : ℝ)⁻¹ * (8 + 51658903 / 1000000000) ≤ (850387611 / 10000000000 : ℝ) := by norm_num
  linarith

/-- h(371/750; κ) > 0.000033545 for every κ ∈ [0.0850387602, 1). -/
theorem h_edge_gt {κ : ℝ} (hk : (850387602 / 10000000000 : ℝ) ≤ κ) (hk1 : κ < 1) :
    (0.000033545 : ℝ) < h F (371 / 750) κ := by
  have hcw := cW_edge_ge F
  have hc1 := cW_lt_one F (by norm_num : (0 : ℝ) < 371 / 750) (by norm_num)
  have m1 := hc_mono_k F hc1 hk hk1
  have m2 := hc_mono_c F (by norm_num : (850387602 / 10000000000 : ℝ) < 1) hcw hc1
  have base : (0.000033545 : ℝ) < hc F (-1857718405 / 10000000000 : ℝ) (850387602 / 10000000000) := by
    unfold hc
    have hA : (0 : ℝ) ≤ (1 - 850387602 / 10000000000) * (9 - 850387602 / 10000000000 - 8 * (-1857718405 / 10000000000 : ℝ)) := by norm_num
    have := sqrt_lt_of_lt_sq F hA (B := 3 + 850387602 / 10000000000 - 4 * 0.000033545) (by norm_num) (by norm_num)
    linarith
  unfold h
  linarith

/-- h(371/750; κ) ≤ 0.000033547 for κ ≤ 0.0850387611. -/
theorem h_edge_le {κ : ℝ} (hk : κ ≤ 850387611 / 10000000000) : h F (371 / 750) κ ≤ 0.000033547 := by
  have hc1 := cW_lt_one F (by norm_num : (0 : ℝ) < 371 / 750) (by norm_num)
  have m1 := hc_mono_k F hc1 hk (by norm_num)
  have m2 := hc_mono_c F (by norm_num : (850387611 / 10000000000 : ℝ) < 1) (cW_edge_le F) (by norm_num)
  have top : hc F (-18577184 / 100000000 : ℝ) (850387611 / 10000000000) ≤ 0.000033547 := by
    unfold hc
    have := le_sqrt_of_sq_le F (a := (1 - 850387611 / 10000000000) * (9 - 850387611 / 10000000000 - 8 * (-18577184 / 100000000 : ℝ)))
      (B := 3 + 850387611 / 10000000000 - 4 * 0.000033547) (by norm_num)
    linarith
  unfold h
  linarith

/-- k*(371/750) = 1 − 2/(2 − c(371/750)) ∈ [0.0849914142, 0.0849914146]. -/
theorem kstar_edge_bounds : (849914142 / 10000000000 : ℝ) ≤ kstar F (371 / 750) ∧
    kstar F (371 / 750) ≤ 849914146 / 10000000000 := by
  have hlo := cW_edge_ge F
  have hhi := cW_edge_le F
  have ks1 : kstar F (371 / 750) ≤ 1 - 2 / (2 - (-1857718405 / 10000000000 : ℝ)) := by
    unfold kstar
    have := div_le_div_of_nonneg_left (by norm_num : (0 : ℝ) ≤ 2) (by linarith)
      (by linarith : 2 - cW F (371 / 750) ≤ 2 - (-1857718405 / 10000000000 : ℝ))
    linarith
  have ks2 : 1 - 2 / (2 - (-18577184 / 100000000 : ℝ)) ≤ kstar F (371 / 750) := by
    unfold kstar
    have := div_le_div_of_nonneg_left (by norm_num : (0 : ℝ) ≤ 2) (by norm_num)
      (by linarith : 2 - (-18577184 / 100000000 : ℝ) ≤ 2 - cW F (371 / 750))
    linarith
  have g1 : (849914142 / 10000000000 : ℝ) ≤ 1 - 2 / (2 - (-18577184 / 100000000 : ℝ)) := by norm_num
  have g2 : 1 - 2 / (2 - (-1857718405 / 10000000000 : ℝ)) ≤ (849914146 / 10000000000 : ℝ) := by norm_num
  exact ⟨g1.trans ks2, ks1.trans g2⟩

end E1PLean
