import Mathlib.Order.LiminfLimsup
import Mathlib.Algebra.Order.Archimedean.Real.Basic
import Riemann.Deep.ExplorePointProc.Lean.Model

/-!
# DPos is D(λ) > 0 for Lamzouri's D(λ) = liminf_{T→∞} #(small gaps up to T)/N(T) (TeX l.183)

Dlam Z λ := liminf_{T→∞} #smallGaps(λ, T)/N(T) (Mathlib's Filter.liminf on ℝ, filter atTop).
dpos_iff: DPos Z λ ↔ 0 < Dlam Z λ, given N(T) → ∞ (which ScaleToInf and RvM imply; nToInf).
-/

namespace PPLean

open Filter Finset

/-- Lamzouri's D(λ), TeX l.183. -/
noncomputable def Dlam (Z : Zeros) (lam : ℝ) : ℝ :=
  liminf (fun T : ℝ => ((smallGaps Z lam T).card : ℝ) / (Z.N T : ℝ)) atTop

lemma card_smallGaps_le (Z : Zeros) (lam T : ℝ) : (smallGaps Z lam T).card ≤ Z.N T := by
  have h := card_le_card (filter_subset (fun k => (Z.γ (k + 1) - Z.γ k) * Z.sc T ≤ lam)
    (Icc 1 (Z.N T)))
  rw [Nat.card_Icc] at h
  unfold smallGaps
  omega

lemma ratio_nonneg (Z : Zeros) (lam T : ℝ) :
    0 ≤ ((smallGaps Z lam T).card : ℝ) / (Z.N T : ℝ) := by positivity

lemma ratio_le_one (Z : Zeros) (lam T : ℝ) :
    ((smallGaps Z lam T).card : ℝ) / (Z.N T : ℝ) ≤ 1 := by
  rcases eq_or_lt_of_le (Nat.cast_nonneg (Z.N T) : (0 : ℝ) ≤ (Z.N T : ℝ)) with h0 | h0
  · rw [← h0, div_zero]
    norm_num
  · rw [div_le_one h0]
    exact_mod_cast card_smallGaps_le Z lam T

/-- N(T) → ∞ follows from log T/(2π) → ∞ and Riemann–von Mangoldt. -/
lemma nToInf {Z : Zeros} (hsc : ScaleToInf Z) (hRvM : RvM Z) :
    ∀ K : ℝ, ∃ T₀ : ℝ, ∀ T : ℝ, T₀ ≤ T → K ≤ (Z.N T : ℝ) := by
  intro K
  obtain ⟨T₁, hT₁⟩ := hsc 1
  obtain ⟨T₂, hT₂⟩ := hRvM 1 one_pos
  refine ⟨max (max T₁ T₂) (max 1 (2 * K)), fun T hT => ?_⟩
  obtain ⟨⟨h1, h2⟩, h3, h4⟩ : (T₁ ≤ T ∧ T₂ ≤ T) ∧ 1 ≤ T ∧ 2 * K ≤ T := by
    simpa only [max_le_iff] using hT
  have hs := hT₁ T h1
  have hr := hT₂ T h2
  have h' := mul_le_mul_of_nonneg_left hs (by linarith : (0 : ℝ) ≤ T)
  linarith

/-- **DPos Z λ ⟺ 0 < D(λ)** with Lamzouri's liminf definition. -/
theorem dpos_iff {Z : Zeros} (hN : ∀ K : ℝ, ∃ T₀ : ℝ, ∀ T : ℝ, T₀ ≤ T → K ≤ (Z.N T : ℝ))
    (lam : ℝ) : DPos Z lam ↔ 0 < Dlam Z lam := by
  have hb : IsBoundedUnder (· ≥ ·) atTop
      (fun T : ℝ => ((smallGaps Z lam T).card : ℝ) / (Z.N T : ℝ)) :=
    isBoundedUnder_of ⟨0, fun T => ratio_nonneg Z lam T⟩
  have hcb : IsCoboundedUnder (· ≥ ·) atTop
      (fun T : ℝ => ((smallGaps Z lam T).card : ℝ) / (Z.N T : ℝ)) :=
    isCoboundedUnder_ge_of_le atTop (fun T => ratio_le_one Z lam T)
  constructor
  · rintro ⟨δ, hδ, T₀, hT⟩
    obtain ⟨T₁, hT₁⟩ := hN 1
    have hev : ∀ᶠ T in atTop, δ ≤ ((smallGaps Z lam T).card : ℝ) / (Z.N T : ℝ) := by
      filter_upwards [eventually_ge_atTop (max T₀ T₁)] with T hT'
      have hNpos : 0 < (Z.N T : ℝ) := by
        have := hT₁ T (le_trans (le_max_right _ _) hT')
        linarith
      rw [le_div_iff₀ hNpos]
      exact hT T (le_trans (le_max_left _ _) hT')
    have := le_liminf_of_le hcb hev
    unfold Dlam
    linarith
  · intro hD
    have hev := eventually_lt_of_lt_liminf (half_lt_self hD) hb
    obtain ⟨T₀, hT₀⟩ := eventually_atTop.1 hev
    refine ⟨Dlam Z lam / 2, by positivity, T₀, fun T hT => ?_⟩
    have h := hT₀ T hT
    rcases eq_or_lt_of_le (Nat.cast_nonneg (Z.N T) : (0 : ℝ) ≤ (Z.N T : ℝ)) with h0 | h0
    · rw [← h0, mul_zero]
      exact Nat.cast_nonneg _
    · rw [lt_div_iff₀ h0] at h
      exact h.le

end PPLean
