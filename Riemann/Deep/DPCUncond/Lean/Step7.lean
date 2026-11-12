import Mathlib.Basic.Real.Basic
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.Positivity
import Mathlib.Tactic.Ring

/-!
# dpcuncond Step 7 (new step (e)): the contradiction from abstract inequalities

All analytic objects at the fixed Step-0 data are real functions of T (T → ∞ as explicit ∀ η > 0 ∃ T₀ ∀ T ≥ T₀
statements): I T = Re ∫_ℝ G_T|R_T|²P(V), M T = Re M_T, Hp T = Re 𝓗_T, S T = Σ|c_n|², R4 T = ∫|R_T|⁴,
P4 T = ∫|P(V)|⁴, Bad T = ∫_{E_T} G_T|R_T|²|P(V)|, mE T = meas E_T, L T = log T, h T = 2πφ/log T (h·L ≤ H₀ = 2πφ),
Nc T = 𝒩(T), Off T = #{off-line zeros, T ≤ γ ≤ 2T}, Φ ε m = Φ(ε) (a-sum truncated at m). Hypotheses, in
PROOF.md's order: h3 (Step 3 identity), hM (DPC Prop 1 + Lemma 6 main term), h6 (Step 6 hairpin bound),
hΦ (Φ(ε) → 0, Step 6f), h2 (Step 2), h1 (Step 1: Re P(V) ≥ 0 outside E_T), hHol (Hölder 2-4-4), h5 (Step 5),
hcov (meas E_T ≤ h𝒩(T) + h). Conclusion: 𝒩(T) > ε₀ T L for all large T.
-/

namespace DPCLean

/-- Step 7 at one height T: Steps 3, 6, 1 and the main term give Bad ≥ (δ₀S₀/4)T. -/
lemma step7_bad {I M Hp S Bad T δ₀ S₀ Φm : ℝ} (hT0 : 0 ≤ T) (hδ₀ : 0 < δ₀)
    (hHp : |Hp| ≤ Φm * T + δ₀ * S₀ / 16 * T) (hΦm : Φm ≤ δ₀ * S₀ / 2)
    (hI : |I - (M + Hp)| ≤ δ₀ * S₀ / 16 * T) (hSS : S₀ ≤ S)
    (hMM : M ≤ -3 * δ₀ * (T / 2) * S * (1 - 1 / 3)) (hB : -Bad - δ₀ * S₀ / 8 * T ≤ I) :
    δ₀ * S₀ / 4 * T ≤ Bad := by
  have hab := le_abs_self Hp
  have hΦT : Φm * T ≤ δ₀ * S₀ / 2 * T := mul_le_mul_of_nonneg_right hΦm hT0
  have hI' := (abs_le.1 hI).2
  have e : -3 * δ₀ * (T / 2) * S * (1 - 1 / 3) = -(δ₀ * T * S) := by ring
  have h4 : δ₀ * T * S₀ ≤ δ₀ * T * S := mul_le_mul_of_nonneg_left hSS (by positivity)
  linarith

/-- Step 7 at one height T: Hölder (Bad⁴ ≤ KG·R4²·P4·mE), Step 5 and the covering bound turn Bad ≥ bT into
b⁴ T L ≤ Kc H₀ (𝒩 + 1), Kc = |KG| K₁² |K₂| + 1. -/
lemma step7_holder {Bad R4 P4 mE LT hT NcT KG K₁ K₂ H₀ T b : ℝ} (hT0 : 0 < T) (hb : 0 ≤ b)
    (hBad : b * T ≤ Bad) (hB4 : Bad ^ 4 ≤ KG * R4 ^ 2 * P4 * mE) (hR0 : 0 ≤ R4) (hR1 : R4 ≤ K₁ * T)
    (hP0 : 0 ≤ P4) (hP1 : P4 ≤ K₂ * T) (hm0 : 0 ≤ mE) (hL0 : 0 ≤ LT) (hhL : hT * LT ≤ H₀)
    (hmE : mE ≤ hT * NcT + hT) (hN : 0 ≤ NcT) :
    b ^ 4 * T * LT ≤ (|KG| * K₁ ^ 2 * |K₂| + 1) * H₀ * (NcT + 1) := by
  have hX0 : 0 ≤ b * T := mul_nonneg hb hT0.le
  have hpow : (b * T) ^ 4 ≤ Bad ^ 4 := pow_le_pow_left₀ hX0 hBad 4
  have hR2 : R4 ^ 2 ≤ K₁ ^ 2 * T ^ 2 := by
    calc R4 ^ 2 ≤ (K₁ * T) ^ 2 := pow_le_pow_left₀ hR0 hR1 2
      _ = K₁ ^ 2 * T ^ 2 := by ring
  have hP' : P4 ≤ |K₂| * T := le_trans hP1 (mul_le_mul_of_nonneg_right (le_abs_self K₂) hT0.le)
  have hB4' : Bad ^ 4 ≤ (|KG| * K₁ ^ 2 * |K₂| + 1) * T ^ 3 * mE := by
    have hY : 0 ≤ R4 ^ 2 * P4 * mE := by positivity
    calc Bad ^ 4 ≤ KG * R4 ^ 2 * P4 * mE := hB4
      _ = KG * (R4 ^ 2 * P4 * mE) := by ring
      _ ≤ |KG| * (R4 ^ 2 * P4 * mE) := mul_le_mul_of_nonneg_right (le_abs_self KG) hY
      _ = |KG| * R4 ^ 2 * P4 * mE := by ring
      _ ≤ |KG| * (K₁ ^ 2 * T ^ 2) * (|K₂| * T) * mE := by
          apply mul_le_mul_of_nonneg_right _ hm0
          exact mul_le_mul (mul_le_mul_of_nonneg_left hR2 (abs_nonneg KG)) hP' hP0 (by positivity)
      _ = (|KG| * K₁ ^ 2 * |K₂|) * (T ^ 3 * mE) := by ring
      _ ≤ (|KG| * K₁ ^ 2 * |K₂| + 1) * (T ^ 3 * mE) := by
          apply mul_le_mul_of_nonneg_right _ (by positivity)
          linarith
      _ = (|KG| * K₁ ^ 2 * |K₂| + 1) * T ^ 3 * mE := by ring
  have hmL : mE * LT ≤ H₀ * (NcT + 1) := by
    calc mE * LT ≤ (hT * NcT + hT) * LT := mul_le_mul_of_nonneg_right hmE hL0
      _ = (hT * LT) * (NcT + 1) := by ring
      _ ≤ H₀ * (NcT + 1) := mul_le_mul_of_nonneg_right hhL (by linarith)
  have hKT : 0 ≤ (|KG| * K₁ ^ 2 * |K₂| + 1) * T ^ 3 := by positivity
  have hchain : T ^ 3 * (b ^ 4 * T * LT) ≤ T ^ 3 * ((|KG| * K₁ ^ 2 * |K₂| + 1) * H₀ * (NcT + 1)) := by
    calc T ^ 3 * (b ^ 4 * T * LT) = (b * T) ^ 4 * LT := by ring
      _ ≤ Bad ^ 4 * LT := mul_le_mul_of_nonneg_right hpow hL0
      _ ≤ ((|KG| * K₁ ^ 2 * |K₂| + 1) * T ^ 3 * mE) * LT := mul_le_mul_of_nonneg_right hB4' hL0
      _ = (|KG| * K₁ ^ 2 * |K₂| + 1) * T ^ 3 * (mE * LT) := by ring
      _ ≤ (|KG| * K₁ ^ 2 * |K₂| + 1) * T ^ 3 * (H₀ * (NcT + 1)) := mul_le_mul_of_nonneg_left hmL hKT
      _ = T ^ 3 * ((|KG| * K₁ ^ 2 * |K₂| + 1) * H₀ * (NcT + 1)) := by ring
  exact le_of_mul_le_mul_left hchain (by positivity)

/-- **Step 7** (PROOF.md, DPC l.766–787 with Steps 1–6), from abstract inequalities. -/
theorem step7 {I M Hp S R4 P4 Bad mE L h Nc Off : ℝ → ℝ} {Φ : ℝ → ℕ → ℝ} {H₀ : ℝ}
    (hH₀ : 0 < H₀)
    (hL : ∃ T₀ : ℝ, ∀ T : ℝ, T₀ ≤ T → 1 ≤ L T)
    (hNc : ∀ T : ℝ, 0 ≤ Nc T)
    (h3 : ∀ η : ℝ, 0 < η → ∃ T₀ : ℝ, ∀ T : ℝ, T₀ ≤ T → |I T - (M T + Hp T)| ≤ η * T)
    (hM : ∃ δ₀ S₀ : ℝ, 0 < δ₀ ∧ 0 < S₀ ∧ ∀ η : ℝ, 0 < η → ∃ T₀ : ℝ, ∀ T : ℝ, T₀ ≤ T →
      S₀ ≤ S T ∧ M T ≤ -3 * δ₀ * (T / 2) * S T * (1 - η))
    (h6 : ∀ ε : ℝ, 0 < ε → ε < 1 → ∀ η : ℝ, 0 < η → ∃ T₀ : ℝ, ∀ T : ℝ, T₀ ≤ T →
      Off T ≤ 2 * ε * T * L T → ∃ m : ℕ, |Hp T| ≤ Φ ε m * T + η * T)
    (hΦ : ∀ τ : ℝ, 0 < τ → ∃ ε₁ : ℝ, 0 < ε₁ ∧ ∀ ε : ℝ, 0 < ε → ε ≤ ε₁ → ∀ m : ℕ, Φ ε m ≤ τ)
    (h2 : ∀ T : ℝ, Off T ≤ 2 * Nc T)
    (h1 : ∀ η : ℝ, 0 < η → ∃ T₀ : ℝ, ∀ T : ℝ, T₀ ≤ T → -Bad T - η * T ≤ I T)
    (hHol : ∃ KG T₀ : ℝ, ∀ T : ℝ, T₀ ≤ T → 0 ≤ Bad T ∧ Bad T ^ 4 ≤ KG * R4 T ^ 2 * P4 T * mE T)
    (h5 : ∃ K₁ K₂ T₀ : ℝ, ∀ T : ℝ, T₀ ≤ T →
      0 ≤ R4 T ∧ R4 T ≤ K₁ * T ∧ 0 ≤ P4 T ∧ P4 T ≤ K₂ * T)
    (hcov : ∃ T₀ : ℝ, ∀ T : ℝ, T₀ ≤ T →
      0 ≤ mE T ∧ 0 ≤ h T ∧ h T * L T ≤ H₀ ∧ mE T ≤ h T * Nc T + h T) :
    ∃ ε₀ : ℝ, 0 < ε₀ ∧ ∃ T₀ : ℝ, ∀ T : ℝ, T₀ ≤ T → ε₀ * T * L T < Nc T := by
  obtain ⟨TL, hTL⟩ := hL
  obtain ⟨δ₀, S₀, hδ₀, hS₀, hMT⟩ := hM
  obtain ⟨KG, TH, hTH⟩ := hHol
  obtain ⟨K₁, K₂, T5, hT5⟩ := h5
  obtain ⟨TC, hTC⟩ := hcov
  have hd : 0 < δ₀ * S₀ := mul_pos hδ₀ hS₀
  obtain ⟨ε₁, hε₁, hΦε⟩ := hΦ (δ₀ * S₀ / 2) (by positivity)
  obtain ⟨Kc, hKc⟩ : ∃ Kc : ℝ, Kc = |KG| * K₁ ^ 2 * |K₂| + 1 := ⟨_, rfl⟩
  have hKc0 : 0 < Kc := by rw [hKc]; positivity
  obtain ⟨D4, hD4⟩ : ∃ D4 : ℝ, D4 = (δ₀ * S₀ / 4) ^ 4 := ⟨_, rfl⟩
  have hD40 : 0 < D4 := by rw [hD4]; positivity
  have hKH : 0 < Kc * H₀ := mul_pos hKc0 hH₀
  obtain ⟨ε₀, hε₀⟩ : ∃ ε₀ : ℝ, ε₀ = min ε₁ (min (1 / 2) (D4 / (2 * (Kc * H₀)))) := ⟨_, rfl⟩
  have hε₀0 : 0 < ε₀ := by
    rw [hε₀]
    exact lt_min hε₁ (lt_min (by norm_num) (by positivity))
  have hε₀1 : ε₀ ≤ ε₁ := by rw [hε₀]; exact min_le_left _ _
  have hε₀2 : ε₀ ≤ 1 / 2 := by rw [hε₀]; exact le_trans (min_le_right _ _) (min_le_left _ _)
  have hε₀3 : ε₀ ≤ D4 / (2 * (Kc * H₀)) := by
    rw [hε₀]; exact le_trans (min_le_right _ _) (min_le_right _ _)
  have hKHε : Kc * H₀ * ε₀ ≤ D4 / 2 := by
    have h4 := mul_le_mul_of_nonneg_left hε₀3 hKH.le
    have e : Kc * H₀ * (D4 / (2 * (Kc * H₀))) = D4 / 2 := by
      rw [mul_div_assoc', div_eq_div_iff (by positivity) (by norm_num)]
      ring
    linarith
  obtain ⟨T3, hT3⟩ := h3 (δ₀ * S₀ / 16) (by positivity)
  obtain ⟨TM, hTM⟩ := hMT (1 / 3) (by norm_num)
  obtain ⟨T6, hT6⟩ := h6 ε₀ hε₀0 (by linarith) (δ₀ * S₀ / 16) (by positivity)
  obtain ⟨T1, hT1⟩ := h1 (δ₀ * S₀ / 8) (by positivity)
  refine ⟨ε₀, hε₀0, max (max (max TL TH) (max T5 TC))
    (max (max T3 TM) (max (max T6 T1) (2 * Kc * H₀ / D4 + 1))), fun T hT => ?_⟩
  simp only [max_le_iff] at hT
  obtain ⟨⟨⟨hTL', hTH'⟩, hT5', hTC'⟩, ⟨hT3', hTM'⟩, ⟨hT6', hT1'⟩, hTbig⟩ := hT
  have hT0 : 0 < T := by
    have h4 : 0 < 2 * Kc * H₀ / D4 := by positivity
    linarith
  have hLT := hTL T hTL'
  have hL0 : 0 ≤ L T := by linarith
  rcases lt_or_ge (ε₀ * T * L T) (Nc T) with hlt | hge
  · exact hlt
  exfalso
  have hOff : Off T ≤ 2 * ε₀ * T * L T := by
    have h4 := h2 T
    linarith
  obtain ⟨m, hm⟩ := hT6 T hT6' hOff
  obtain ⟨hSS, hMM⟩ := hTM T hTM'
  have hBad := step7_bad hT0.le hδ₀ hm (hΦε ε₀ hε₀0 hε₀1 m) (hT3 T hT3') hSS hMM (hT1 T hT1')
  obtain ⟨hB0, hB4⟩ := hTH T hTH'
  obtain ⟨hR0, hR1, hP0, hP1⟩ := hT5 T hT5'
  obtain ⟨hm0, hh0, hhL, hmE⟩ := hTC T hTC'
  have hH := step7_holder hT0 (by positivity) hBad hB4 hR0 hR1 hP0 hP1 hm0 hL0 hhL hmE (hNc T)
  rw [← hKc, ← hD4] at hH
  have h5' : Kc * H₀ * (Nc T + 1) ≤ Kc * H₀ * (ε₀ * T * L T + 1) :=
    mul_le_mul_of_nonneg_left (by linarith) hKH.le
  have hTL0 : 0 ≤ T * L T := mul_nonneg hT0.le hL0
  have h8 : Kc * H₀ * ε₀ * (T * L T) ≤ D4 / 2 * (T * L T) := mul_le_mul_of_nonneg_right hKHε hTL0
  have h9 : Kc * H₀ < D4 / 2 * T := by
    have h4 : 2 * Kc * H₀ / D4 < T := by linarith
    rw [div_lt_iff₀ hD40] at h4
    linarith
  have h10 : D4 / 2 * T ≤ D4 / 2 * (T * L T) :=
    mul_le_mul_of_nonneg_left (le_mul_of_one_le_right hT0.le hLT) (by positivity)
  have e : Kc * H₀ * (ε₀ * T * L T + 1) = Kc * H₀ * ε₀ * (T * L T) + Kc * H₀ := by ring
  have e2 : Kc * H₀ * (Nc T + 1) = Kc * (H₀ * (Nc T + 1)) := by ring
  linarith

end DPCLean

