import Riemann.Deep.DPCBeyond.Lean.Rpow
import Riemann.Deep.DPCUncond.Lean.Step7

/-!
# dpcbeyond Step 7′ (new step (e)): Hölder with ‖R_T‖_q, q > 2, in place of ‖R_T‖₄

PROOF.md Step 7′: with 2 < q < κ/θ (M3) and 2/q + 1/s < 1, Re P(V) ≥ 0 outside E_T (Step 1) and Step 6e′,
δ₁T ≤ ∫_{E_T} G_T|R_T|²|P(V)| ≤ ‖R_T‖_q²‖P(V)‖_s meas(E_T)^{1−2/q−1/s} ≪ T^{2/q+1/s} meas(E_T)^{1−2/q−1/s}, so
meas(E_T) ≥ c₂T with c₂ independent of ε; the rest of Step 7 is unchanged. step7B_meas: the measure bound at one
height (c₂ = (δ₁/K_B)^{1/α}); step7B: dpclean's step7 (Riemann/Deep/DPCUncond/Lean/Step7.lean) with its fourth-moment
Hölder hypothesis replaced by the q-th moment one; dpclean's step7_bad is reused verbatim.
-/

namespace DPCBLean

open DPCLean

variable {ex lg : ℝ → ℝ}

/-- Step 7′ at one height T: bT ≤ Bad ≤ K_B T^{1−α} meas(E_T)^α, α > 0, gives meas(E_T) ≥ (b/K_B)^{1/α} T. -/
lemma step7B_meas (h : ExpLog ex lg) {Bad mE T b KB α : ℝ} (hT0 : 0 < T) (hb : 0 < b) (hKB : 0 < KB)
    (hα : 0 < α) (hBad : b * T ≤ Bad) (hHol : Bad ≤ KB * rp ex lg T (1 - α) * rp ex lg mE α) :
    rp ex lg (b / KB) (1 / α) * T ≤ mE := by
  have hbT : 0 < b * T := mul_pos hb hT0
  have hm0 : 0 < mE := by
    by_contra hn
    rw [rp_of_nonpos (not_lt.1 hn)] at hHol
    linarith only [hbT, hBad, hHol]
  have e1 : rp ex lg T (1 - α) = T * rp ex lg T (-α) := by
    rw [sub_eq_add_neg, rp_add h hT0, rp_one h hT0]
  have hTa := rp_neg h hT0 α
  have hc0 : 0 < b / KB := div_pos hb hKB
  obtain ⟨c₂, hc₂⟩ : ∃ c : ℝ, c = rp ex lg (b / KB) (1 / α) := ⟨_, rfl⟩
  rw [← hc₂]
  have hc₂0 : 0 < c₂ := by rw [hc₂]; exact rp_pos h hc0 _
  have e2 : rp ex lg (c₂ * T) α = b / KB * rp ex lg T α := by
    rw [rp_mul h hc₂0 hT0, hc₂, rp_rp h hc0, one_div_mul_cancel hα.ne', rp_one h hc0]
  have h3 : b / KB * rp ex lg T α ≤ rp ex lg mE α := by
    have hTa0 : 0 < rp ex lg T α := rp_pos h hT0 α
    have h4 : T * b ≤ T * (KB * rp ex lg T (-α) * rp ex lg mE α) := by
      calc T * b = b * T := by ring
        _ ≤ Bad := hBad
        _ ≤ KB * rp ex lg T (1 - α) * rp ex lg mE α := hHol
        _ = T * (KB * rp ex lg T (-α) * rp ex lg mE α) := by rw [e1]; ring
    have h5 : b ≤ KB * rp ex lg T (-α) * rp ex lg mE α := le_of_mul_le_mul_left h4 hT0
    rw [div_mul_eq_mul_div, div_le_iff₀ hKB]
    calc b * rp ex lg T α ≤ KB * rp ex lg T (-α) * rp ex lg mE α * rp ex lg T α :=
          mul_le_mul_of_nonneg_right h5 hTa0.le
      _ = rp ex lg mE α * KB * (rp ex lg T (-α) * rp ex lg T α) := by ring
      _ = rp ex lg mE α * KB := by rw [hTa, mul_one]
  rw [← e2] at h3
  by_contra hn
  have h7 := rp_lt_rp h hm0 (not_le.1 hn) hα
  linarith only [h3, h7]

/-- **Step 7′** (PROOF.md; DPC l.766–787 with Steps 1–6′): dpclean's step7 with Hölder at exponents
(q/2, s, 1/α), α = 1 − 2/q − 1/s > 0, and ∫|R_T|^q ≤ K₁T (Lemma M3), ∫|P(V)|^s ≤ K₂T (5c). -/
theorem step7B (hEL : ExpLog ex lg) {I M Hp S Rq Ps Bad mE L h Nc Off : ℝ → ℝ} {Φ : ℝ → ℕ → ℝ}
    {H₀ q s : ℝ} (hH₀ : 0 < H₀) (hq : 2 < q) (hs : 0 < s) (hqs : 2 / q + 1 / s < 1)
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
    (hHol : ∃ KG T₀ : ℝ, ∀ T : ℝ, T₀ ≤ T → Bad T ≤ KG * rp ex lg (Rq T) (2 / q) *
      rp ex lg (Ps T) (1 / s) * rp ex lg (mE T) (1 - 2 / q - 1 / s))
    (h5 : ∃ K₁ K₂ T₀ : ℝ, ∀ T : ℝ, T₀ ≤ T → Rq T ≤ K₁ * T ∧ Ps T ≤ K₂ * T)
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
  have hα : 0 < 1 - 2 / q - 1 / s := by linarith only [hqs]
  have hK₁ : 0 < |K₁| + 1 := by positivity
  have hK₂ : 0 < |K₂| + 1 := by positivity
  obtain ⟨KB, hKB⟩ : ∃ KB : ℝ, KB = (|KG| + 1) * rp ex lg (|K₁| + 1) (2 / q) *
      rp ex lg (|K₂| + 1) (1 / s) := ⟨_, rfl⟩
  have hKB0 : 0 < KB := by
    rw [hKB]
    exact mul_pos (mul_pos (by positivity) (rp_pos hEL hK₁ _)) (rp_pos hEL hK₂ _)
  have hb : 0 < δ₀ * S₀ / 4 := by positivity
  obtain ⟨c₂, hc₂⟩ : ∃ c : ℝ, c = rp ex lg (δ₀ * S₀ / 4 / KB) (1 / (1 - 2 / q - 1 / s)) := ⟨_, rfl⟩
  have hc₂0 : 0 < c₂ := by rw [hc₂]; exact rp_pos hEL (div_pos hb hKB0) _
  obtain ⟨ε₀, hε₀⟩ : ∃ ε₀ : ℝ, ε₀ = min ε₁ (min (1 / 2) (c₂ / (2 * H₀))) := ⟨_, rfl⟩
  have hε₀0 : 0 < ε₀ := by
    rw [hε₀]
    exact lt_min hε₁ (lt_min (by norm_num) (by positivity))
  have hε₀1 : ε₀ ≤ ε₁ := by rw [hε₀]; exact min_le_left _ _
  have hε₀2 : ε₀ ≤ 1 / 2 := by rw [hε₀]; exact le_trans (min_le_right _ _) (min_le_left _ _)
  have hε₀3 : ε₀ ≤ c₂ / (2 * H₀) := by
    rw [hε₀]; exact le_trans (min_le_right _ _) (min_le_right _ _)
  have hHε : H₀ * ε₀ ≤ c₂ / 2 := by
    have h4 := mul_le_mul_of_nonneg_left hε₀3 hH₀.le
    have e : H₀ * (c₂ / (2 * H₀)) = c₂ / 2 := by field_simp
    linarith only [h4, e]
  obtain ⟨T3, hT3⟩ := h3 (δ₀ * S₀ / 16) (by positivity)
  obtain ⟨TM, hTM⟩ := hMT (1 / 3) (by norm_num)
  obtain ⟨T6, hT6⟩ := h6 ε₀ hε₀0 (by linarith only [hε₀2]) (δ₀ * S₀ / 16) (by positivity)
  obtain ⟨T1, hT1⟩ := h1 (δ₀ * S₀ / 8) (by positivity)
  refine ⟨ε₀, hε₀0, max (max (max TL TH) (max T5 TC))
    (max (max T3 TM) (max (max T6 T1) (2 * H₀ / c₂ + 1))), fun T hT => ?_⟩
  simp only [max_le_iff] at hT
  obtain ⟨⟨⟨hTL', hTH'⟩, hT5', hTC'⟩, ⟨hT3', hTM'⟩, ⟨hT6', hT1'⟩, hTbig⟩ := hT
  have hT0 : 0 < T := by
    have h4 : 0 < 2 * H₀ / c₂ := by positivity
    linarith only [h4, hTbig]
  have hLT := hTL T hTL'
  have hL0 : 0 ≤ L T := by linarith only [hLT]
  rcases lt_or_ge (ε₀ * T * L T) (Nc T) with hlt | hge
  · exact hlt
  exfalso
  have hOff : Off T ≤ 2 * ε₀ * T * L T := by
    have h4 := h2 T
    linarith only [h4, hge]
  obtain ⟨m, hm⟩ := hT6 T hT6' hOff
  obtain ⟨hSS, hMM⟩ := hTM T hTM'
  have hBad := step7_bad hT0.le hδ₀ hm (hΦε ε₀ hε₀0 hε₀1 m) (hT3 T hT3') hSS hMM (hT1 T hT1')
  have hB := hTH T hTH'
  obtain ⟨hR1, hP1⟩ := hT5 T hT5'
  obtain ⟨hm0, hh0, hhL, hmE⟩ := hTC T hTC'
  have hR' : rp ex lg (Rq T) (2 / q) ≤ rp ex lg ((|K₁| + 1) * T) (2 / q) :=
    rp_le_rp hEL (le_trans hR1 (mul_le_mul_of_nonneg_right (by linarith only [le_abs_self K₁]) hT0.le))
      (by positivity)
  have hP' : rp ex lg (Ps T) (1 / s) ≤ rp ex lg ((|K₂| + 1) * T) (1 / s) :=
    rp_le_rp hEL (le_trans hP1 (mul_le_mul_of_nonneg_right (by linarith only [le_abs_self K₂]) hT0.le))
      (by positivity)
  have hX0 : 0 ≤ rp ex lg (Rq T) (2 / q) * rp ex lg (Ps T) (1 / s) *
      rp ex lg (mE T) (1 - 2 / q - 1 / s) :=
    mul_nonneg (mul_nonneg (rp_nonneg hEL _ _) (rp_nonneg hEL _ _)) (rp_nonneg hEL _ _)
  have hHol' : Bad T ≤ KB * rp ex lg T (1 - (1 - 2 / q - 1 / s)) *
      rp ex lg (mE T) (1 - 2 / q - 1 / s) := by
    calc Bad T ≤ KG * rp ex lg (Rq T) (2 / q) * rp ex lg (Ps T) (1 / s) *
          rp ex lg (mE T) (1 - 2 / q - 1 / s) := hB
      _ = KG * (rp ex lg (Rq T) (2 / q) * rp ex lg (Ps T) (1 / s) *
          rp ex lg (mE T) (1 - 2 / q - 1 / s)) := by ring
      _ ≤ (|KG| + 1) * (rp ex lg (Rq T) (2 / q) * rp ex lg (Ps T) (1 / s) *
          rp ex lg (mE T) (1 - 2 / q - 1 / s)) :=
          mul_le_mul_of_nonneg_right (by linarith only [le_abs_self KG]) hX0
      _ ≤ (|KG| + 1) * (rp ex lg ((|K₁| + 1) * T) (2 / q) * rp ex lg ((|K₂| + 1) * T) (1 / s) *
          rp ex lg (mE T) (1 - 2 / q - 1 / s)) := by
          apply mul_le_mul_of_nonneg_left _ (by positivity)
          apply mul_le_mul_of_nonneg_right _ (rp_nonneg hEL _ _)
          exact mul_le_mul hR' hP' (rp_nonneg hEL _ _) (rp_nonneg hEL _ _)
      _ = KB * rp ex lg T (1 - (1 - 2 / q - 1 / s)) * rp ex lg (mE T) (1 - 2 / q - 1 / s) := by
          rw [hKB, rp_mul hEL hK₁ hT0, rp_mul hEL hK₂ hT0,
            show 1 - (1 - 2 / q - 1 / s) = 2 / q + 1 / s by ring, rp_add hEL hT0]
          ring
  have hmeas := step7B_meas hEL hT0 hb hKB0 hα hBad hHol'
  rw [← hc₂] at hmeas
  have h8 : c₂ * T * L T ≤ H₀ * (Nc T + 1) := by
    calc c₂ * T * L T ≤ mE T * L T := mul_le_mul_of_nonneg_right hmeas hL0
      _ ≤ (h T * Nc T + h T) * L T := mul_le_mul_of_nonneg_right hmE hL0
      _ = (h T * L T) * (Nc T + 1) := by ring
      _ ≤ H₀ * (Nc T + 1) := mul_le_mul_of_nonneg_right hhL (by linarith only [hNc T])
  have h9 : H₀ * Nc T ≤ H₀ * (ε₀ * T * L T) := mul_le_mul_of_nonneg_left hge hH₀.le
  have hTL0 : 0 ≤ T * L T := mul_nonneg hT0.le hL0
  have h10 : H₀ * (ε₀ * T * L T) ≤ c₂ / 2 * (T * L T) := by
    calc H₀ * (ε₀ * T * L T) = H₀ * ε₀ * (T * L T) := by ring
      _ ≤ c₂ / 2 * (T * L T) := mul_le_mul_of_nonneg_right hHε hTL0
  have h11 : H₀ < c₂ / 2 * T := by
    have h4 : 2 * H₀ / c₂ < T := by linarith only [hTbig]
    rw [div_lt_iff₀ hc₂0] at h4
    linarith only [h4]
  have h12 : c₂ / 2 * T ≤ c₂ / 2 * (T * L T) :=
    mul_le_mul_of_nonneg_left (le_mul_of_one_le_right hT0.le hLT) (by positivity)
  have e3 : c₂ * T * L T = 2 * (c₂ / 2 * (T * L T)) := by ring
  have e4 : H₀ * (Nc T + 1) = H₀ * Nc T + H₀ := by ring
  linarith only [h8, h9, h10, h11, h12, e3, e4]

end DPCBLean

