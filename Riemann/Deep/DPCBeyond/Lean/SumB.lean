import Riemann.Deep.DPCBeyond.Lean.HairpinB

/-!
# dpcbeyond Step 6e′ (d)–(e), the class sum (new step (d), part 2)

step6B_of_classes: from DensityK and Step6ClassesB (HairpinB.lean), |𝓗_T| ≤ (K₁Phi(2θ, 1 − 1/w) + K₂Phi(0, 1/v))T + o(T)
whenever #Z_off ≤ 2εTL (Step 2's output); dpclean's hairpin_sum (Riemann/Deep/DPCUncond/Lean/Hairpin.lean) applied
twice. PhiB_small: Φ(ε) → 0 by dpclean's Phi_small twice (PROOF.md 6e′ (e)).
-/

namespace DPCBLean

open Finset PPLean DPCLean

variable {ex lg : ℝ → ℝ}

/-- **Step 6e′ from the density theorem and the class bounds** (6e′ (d) bookkeeping and the class sum proved):
|𝓗_T| ≤ (K₁ Phi(2θ, 1 − 1/w) + K₂ Phi(0, 1/v)) T + o(T) under Step 2's #Z_off ≤ 2εTL. -/
theorem step6B_of_classes (hEL : ExpLog ex lg) {Z : Zeros} {β : ℕ → ℝ} {κ θ : ℝ} {Hp : ℝ → ℝ}
    (hθ : 0 ≤ θ) (hden : DensityK Z β ex lg κ) (h6c : Step6ClassesB Z β ex lg κ θ Hp)
    {u w : ℝ} (hu : 1 < u) (hw : 1 < w) (hv : 0 < 1 - 1 / u - 1 / w) (huθ : 2 * u * θ < κ) :
    ∃ (C : ℕ) (K₁ K₂ : ℝ), 0 ≤ K₁ ∧ 0 ≤ K₂ ∧ ∀ c ε : ℝ, 0 < ε → ε < 1 → ∀ η : ℝ, 0 < η →
      ∃ T₀ : ℝ, ∀ T : ℝ, T₀ ≤ T → ((offSet Z β T).card : ℝ) ≤ 2 * ε * T * lg T →
        ∃ m : ℕ, |Hp T| ≤ (K₁ * Phi ex lg C (2 * θ) κ (1 - 1 / w) c ε m +
          K₂ * Phi ex lg C 0 κ (1 - 1 / u - 1 / w) c ε m) * T + η * T := by
  obtain ⟨C, K, hK, Ha, Y₁, Y₂, MY, MO, Mt, hcl⟩ := h6c u w hu hw hv huθ
  obtain ⟨Cd, Td, hCd, hdT⟩ := hden
  have hw0 : 0 < w := by linarith
  have hE1 : 0 ≤ 1 - 1 / w := by
    have h1 : 1 / w < 1 := (div_lt_one hw0).2 hw
    linarith
  have hK0 : 0 ≤ K * (K + 1) ^ 2 := by positivity
  refine ⟨2 * C + 1, K * (K + 1) ^ 2 * max 1 (ex (lg Cd * (1 - 1 / w))),
    K * (K + 1) ^ 2 * max 1 (ex (lg Cd * (1 - 1 / u - 1 / w))),
    mul_nonneg hK0 (le_trans zero_le_one (le_max_left _ _)),
    mul_nonneg hK0 (le_trans zero_le_one (le_max_left _ _)), fun c ε hε _ η hη => ?_⟩
  obtain ⟨T₀, hT₀⟩ := hcl η hη
  refine ⟨max (max T₀ Td) (ex 1), fun T hT hoff => ?_⟩
  have hT1 : T₀ ≤ T := le_trans (le_trans (le_max_left _ _) (le_max_left _ _)) hT
  have hT2 : Td ≤ T := le_trans (le_trans (le_max_right _ _) (le_max_left _ _)) hT
  have hT3 : ex 1 ≤ T := le_trans (le_max_right _ _) hT
  have hT0 : 0 < T := lt_of_lt_of_le (hEL.pos 1) hT3
  have hlg : 1 ≤ lg T := hEL.le_lg hT3
  have hTL : 0 ≤ T * lg T := mul_nonneg hT0.le (by linarith)
  obtain ⟨hdec, hMO, hcls⟩ := hT₀ T hT1
  refine ⟨Mt T, ?_⟩
  have hs0 : ∀ a : ℕ, 0 < min (3 * ε) (Cd * ex (-(κ * a))) :=
    fun a => lt_min (by positivity) (mul_pos hCd (hEL.pos _))
  have hb : ∀ a : ℕ, Ha T a ≤
      K * (K + 1) ^ 2 * T * ((a : ℝ) + 2) ^ (2 * C) * ex (2 * θ) ^ (a + 1) *
        ex (lg (min (3 * ε) (Cd * ex (-(κ * a)))) * (1 - 1 / w)) +
      K * (K + 1) ^ 2 * T * ((a : ℝ) + 2) ^ (2 * C) * ex 0 ^ (a + 1) *
        ex (lg (min (3 * ε) (Cd * ex (-(κ * a)))) * (1 - 1 / u - 1 / w)) := by
    intro a
    obtain ⟨hHa, hMY, hY₁, hY₂⟩ := hcls a
    have ha0 : (0 : ℝ) ≤ a := Nat.cast_nonneg a
    have hB : 1 ≤ ((a : ℝ) + 2) ^ C := one_le_pow₀ (by linarith)
    have h := class_bound hEL hT0 hlg hK hB (pow_nonneg (hEL.pos _).le _) (hs0 a) hw hu hv
      (class_count hε hTL hoff (hdT T hT2) a) hMO hMY hY₁ hY₂ hHa
    refine le_trans h (le_of_eq ?_)
    rw [rp_of_pos (hs0 a), rp_of_pos (hs0 a), hEL.zero, one_pow, pow_mul',
      mul_comm (1 - 1 / w) (lg _), mul_comm (1 - 1 / u - 1 / w) (lg _)]
    ring
  have hsum1 := hairpin_sum hEL (C := 2 * C) (θ := 2 * θ) (κ := κ) (E := 1 - 1 / w) (c := c)
    (ε := ε) (K := K * (K + 1) ^ 2) (T := T) (Cd := Cd) (Mt T)
    (fun a => K * (K + 1) ^ 2 * T * ((a : ℝ) + 2) ^ (2 * C) * ex (2 * θ) ^ (a + 1) *
      ex (lg (min (3 * ε) (Cd * ex (-(κ * a)))) * (1 - 1 / w)))
    (fun a => min (3 * ε) (Cd * ex (-(κ * a)))) (by linarith) hE1 hK0 hT0.le hCd hs0
    (fun a => min_le_left _ _) (fun a => min_le_right _ _) (fun a => le_rfl)
  have hsum2 := hairpin_sum hEL (C := 2 * C) (θ := 0) (κ := κ) (E := 1 - 1 / u - 1 / w) (c := c)
    (ε := ε) (K := K * (K + 1) ^ 2) (T := T) (Cd := Cd) (Mt T)
    (fun a => K * (K + 1) ^ 2 * T * ((a : ℝ) + 2) ^ (2 * C) * ex 0 ^ (a + 1) *
      ex (lg (min (3 * ε) (Cd * ex (-(κ * a)))) * (1 - 1 / u - 1 / w)))
    (fun a => min (3 * ε) (Cd * ex (-(κ * a)))) le_rfl hv.le hK0 hT0.le hCd hs0
    (fun a => min_le_left _ _) (fun a => min_le_right _ _) (fun a => le_rfl)
  have hsplit : ∑ a ∈ range (Mt T + 1), Ha T a ≤
      ∑ a ∈ range (Mt T + 1), K * (K + 1) ^ 2 * T * ((a : ℝ) + 2) ^ (2 * C) * ex (2 * θ) ^ (a + 1) *
        ex (lg (min (3 * ε) (Cd * ex (-(κ * a)))) * (1 - 1 / w)) +
      ∑ a ∈ range (Mt T + 1), K * (K + 1) ^ 2 * T * ((a : ℝ) + 2) ^ (2 * C) * ex 0 ^ (a + 1) *
        ex (lg (min (3 * ε) (Cd * ex (-(κ * a)))) * (1 - 1 / u - 1 / w)) := by
    rw [← sum_add_distrib]
    exact sum_le_sum (fun a _ => hb a)
  rw [add_mul]
  linarith only [hdec, hsplit, hsum1, hsum2]

/-- **6e′ (e)**: Φ(ε) = K₁ Phi(2θ, 1 − 1/w) + K₂ Phi(0, 1/v) → 0 as ε → 0, uniformly in the truncation,
when 2θ < κ(1 − 1/w), c > 0 and 2θc < 1 − 1/w (PROOF.md: A = c log(1/ε), 0 < c < (1 − 1/w)/(2θ)). -/
theorem PhiB_small (hEL : ExpLog ex lg) (C : ℕ) {κ θ w u c K₁ K₂ : ℝ} (hK₁ : 0 ≤ K₁) (hK₂ : 0 ≤ K₂)
    (hκ : 0 < κ) (hw : 1 < w) (hv : 0 < 1 - 1 / u - 1 / w) (hθw : 2 * θ < κ * (1 - 1 / w))
    (hc : 0 < c) (hcθ : c * (2 * θ) < 1 - 1 / w) :
    ∀ τ : ℝ, 0 < τ → ∃ ε₁ : ℝ, 0 < ε₁ ∧ ∀ ε : ℝ, 0 < ε → ε ≤ ε₁ → ∀ m : ℕ,
      K₁ * Phi ex lg C (2 * θ) κ (1 - 1 / w) c ε m + K₂ * Phi ex lg C 0 κ (1 - 1 / u - 1 / w) c ε m
        ≤ τ := by
  intro τ hτ
  have hw0 : 0 < w := by linarith
  have hE1 : 0 ≤ 1 - 1 / w := by
    have h1 : 1 / w < 1 := (div_lt_one hw0).2 hw
    linarith
  obtain ⟨ε₁, h₁, hb₁⟩ := Phi_small hEL C hE1 hθw hc hcθ (τ / (2 * (K₁ + 1))) (by positivity)
  obtain ⟨ε₂, h₂, hb₂⟩ := Phi_small hEL C hv.le (mul_pos hκ hv) hc (by rw [mul_zero]; exact hv)
    (τ / (2 * (K₂ + 1))) (by positivity)
  refine ⟨min ε₁ ε₂, lt_min h₁ h₂, fun ε hε hεε m => ?_⟩
  have b1 := hb₁ ε hε (le_trans hεε (min_le_left _ _)) m
  have b2 := hb₂ ε hε (le_trans hεε (min_le_right _ _)) m
  have c1 := mul_le_mul_of_nonneg_left b1 hK₁
  have c2 := mul_le_mul_of_nonneg_left b2 hK₂
  have d1 : K₁ * (τ / (2 * (K₁ + 1))) ≤ τ / 2 := by
    rw [mul_div_assoc', div_le_iff₀ (by positivity)]
    have e : τ / 2 * (2 * (K₁ + 1)) = K₁ * τ + τ := by ring
    rw [e]
    linarith
  have d2 : K₂ * (τ / (2 * (K₂ + 1))) ≤ τ / 2 := by
    rw [mul_div_assoc', div_le_iff₀ (by positivity)]
    have e : τ / 2 * (2 * (K₂ + 1)) = K₂ * τ + τ := by ring
    rw [e]
    linarith
  linarith

end DPCBLean

