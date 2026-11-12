import Riemann.Deep.DPCUncond.Lean.DLam

/-!
# dpcbeyond: a Theorem-U-type bound with threshold φ₀ gives D(λ) > 0 for every λ > φ₀ (new step (f), part 1)

UFrom Z φ₀: for every φ ∈ (φ₀, 1), 𝒩(T) ≥ c(φ) T sc T for large T (sc T = log T/(2π)); Theorem U is UFrom Z (2/3),
Theorem U′ is UFrom Z (max{2/(2+κ), 1/2}). dpos_of_from, nToInf_of_from, dlam_pos_of_from: dpclean's dpos_of_U,
nToInf_of_U, dlam_pos_of_U (Riemann/Deep/DPCUncond/Lean/DLam.lean) with 2/3 replaced by φ₀ ∈ [0, 1) (proofs copied
with that change; ScaleHalf, RvMUp, gapSet_half_sub, card_gapSet_le and pplean's DPos, Dlam, dpos_iff reused).
-/

namespace DPCBLean

open Finset PPLean DPCLean

/-- 𝒩_φ(T) ≥ c(φ) T sc T for every φ ∈ (φ₀, 1), for large T. -/
def UFrom (Z : Zeros) (φ₀ : ℝ) : Prop :=
  ∀ φ : ℝ, φ₀ < φ → φ < 1 → ∃ c : ℝ, 0 < c ∧ ∃ T₀ : ℝ, ∀ T : ℝ, T₀ ≤ T →
    c * T * Z.sc T ≤ ((gapSet Z φ T).card : ℝ)

/-- UFrom Z φ₀ ⇒ D(λ) > 0 (DPos form) for every λ > φ₀. -/
theorem dpos_of_from {Z : Zeros} (hZ : Z.Ordered) {φ₀ : ℝ} (hφ₀ : 0 ≤ φ₀) (hφ₀1 : φ₀ < 1)
    (hU : UFrom Z φ₀) (hs1 : ScaleToInf Z) (hs2 : ScaleHalf Z) (hR : RvMUp Z) {lam : ℝ}
    (hlam : φ₀ < lam) : DPos Z lam := by
  obtain ⟨φ, hφ⟩ : ∃ φ : ℝ, φ = (φ₀ + min lam 1) / 2 := ⟨_, rfl⟩
  have hm1 : min lam 1 ≤ lam := min_le_left _ _
  have hm2 : min lam 1 ≤ 1 := min_le_right _ _
  have hm3 : φ₀ < min lam 1 := lt_min hlam hφ₀1
  have hφ1 : φ₀ < φ := by rw [hφ]; linarith
  have hφ2 : φ < 1 := by rw [hφ]; linarith
  have hφl : φ < lam := by rw [hφ]; linarith
  have hφ0 : 0 < φ := by linarith
  have hlam0 : 0 < lam := by linarith
  obtain ⟨c, hc, T₁, hT₁⟩ := hU φ hφ1 hφ2
  obtain ⟨k, hk⟩ := hs2
  obtain ⟨K, T₂, hT₂⟩ := hR
  obtain ⟨T₃, hT₃⟩ := hs1 (max 1 (max (2 * k) (φ * k / (lam - φ))))
  refine ⟨c / (4 * (|K| + 1)), by positivity,
    max (max (2 * T₁) (2 * T₃)) (max T₂ (max T₃ 1)), fun T hT => ?_⟩
  simp only [max_le_iff] at hT
  obtain ⟨⟨hT1, hT3a⟩, hT2, hT3b, hT1'⟩ := hT
  have hT0 : 0 < T := by linarith
  have hsT := hT₃ T hT3b
  have hsT2 := hT₃ (T / 2) (by linarith)
  simp only [max_le_iff] at hsT hsT2
  obtain ⟨hsT1, hsTk, -⟩ := hsT
  obtain ⟨-, -, hsT2φ⟩ := hsT2
  have hkT := hk T hT0
  have hscale : φ * Z.sc T ≤ lam * Z.sc (T / 2) := by
    have h1 : φ * k ≤ Z.sc (T / 2) * (lam - φ) := by
      rwa [div_le_iff₀ (by linarith)] at hsT2φ
    have h2 : φ * Z.sc T ≤ φ * (Z.sc (T / 2) + k) := mul_le_mul_of_nonneg_left hkT hφ0.le
    linarith
  have hcard : ((gapSet Z φ (T / 2)).card : ℝ) ≤ ((smallGaps Z lam T).card : ℝ) := by
    exact_mod_cast card_le_card (gapSet_half_sub hZ hφ0 hlam0 hscale)
  have hU' := hT₁ (T / 2) (by linarith)
  have hhalf : Z.sc T / 2 ≤ Z.sc (T / 2) := by linarith
  have hTs : 0 ≤ T * Z.sc T := mul_nonneg hT0.le (by linarith)
  have hN' : (Z.N T : ℝ) ≤ (|K| + 1) * (T * Z.sc T) := by
    have h1 := hT₂ T hT2
    have h2 : K * T * Z.sc T ≤ (|K| + 1) * (T * Z.sc T) := by
      calc K * T * Z.sc T = K * (T * Z.sc T) := by ring
        _ ≤ (|K| + 1) * (T * Z.sc T) :=
          mul_le_mul_of_nonneg_right (by linarith [le_abs_self K]) hTs
    linarith
  have hfin : c / (4 * (|K| + 1)) * (Z.N T : ℝ) ≤ c / 4 * (T * Z.sc T) := by
    calc c / (4 * (|K| + 1)) * (Z.N T : ℝ) ≤ c / (4 * (|K| + 1)) * ((|K| + 1) * (T * Z.sc T)) :=
          mul_le_mul_of_nonneg_left hN' (by positivity)
      _ = c / 4 * (T * Z.sc T) := by
          rw [div_mul_eq_mul_div, div_eq_iff (by positivity)]
          ring
  have hmid : c * (T / 2) * (Z.sc T / 2) ≤ c * (T / 2) * Z.sc (T / 2) :=
    mul_le_mul_of_nonneg_left hhalf (by positivity)
  linarith

/-- N(T) → ∞, from UFrom itself (N(T) ≥ 𝒩_φ(T/2), φ = (φ₀ + 1)/2). -/
lemma nToInf_of_from {Z : Zeros} {φ₀ : ℝ} (hφ₀1 : φ₀ < 1) (hU : UFrom Z φ₀) (hs1 : ScaleToInf Z) :
    ∀ K : ℝ, ∃ T₀ : ℝ, ∀ T : ℝ, T₀ ≤ T → K ≤ (Z.N T : ℝ) := by
  intro K
  obtain ⟨c, hc, T₁, hT₁⟩ := hU ((φ₀ + 1) / 2) (by linarith) (by linarith)
  obtain ⟨T₃, hT₃⟩ := hs1 1
  refine ⟨max (max (2 * T₁) (2 * T₃)) (2 * |K| / c), fun T hT => ?_⟩
  simp only [max_le_iff] at hT
  obtain ⟨⟨h1, h3⟩, h4⟩ := hT
  have hs := hT₃ (T / 2) (by linarith)
  have hU' := hT₁ (T / 2) (by linarith)
  have hcard : ((gapSet Z ((φ₀ + 1) / 2) (T / 2)).card : ℝ) ≤ (Z.N T : ℝ) := by
    have h := card_gapSet_le Z ((φ₀ + 1) / 2) (T / 2)
    rw [show 2 * (T / 2) = T by ring] at h
    exact_mod_cast h
  have hK0 : 0 ≤ 2 * |K| / c := by positivity
  have h5 : c * (T / 2) ≤ c * (T / 2) * Z.sc (T / 2) :=
    le_mul_of_one_le_right (by nlinarith) hs
  have h6 : |K| ≤ c * (T / 2) := by
    rw [div_le_iff₀ hc] at h4
    linarith
  linarith [le_abs_self K]

/-- UFrom Z φ₀ ⇒ 0 < D(λ) (Lamzouri's liminf) for every λ > φ₀. -/
theorem dlam_pos_of_from {Z : Zeros} (hZ : Z.Ordered) {φ₀ : ℝ} (hφ₀ : 0 ≤ φ₀) (hφ₀1 : φ₀ < 1)
    (hU : UFrom Z φ₀) (hs1 : ScaleToInf Z) (hs2 : ScaleHalf Z) (hR : RvMUp Z) {lam : ℝ}
    (hlam : φ₀ < lam) : 0 < Dlam Z lam :=
  (dpos_iff (nToInf_of_from hφ₀1 hU hs1) lam).1 (dpos_of_from hZ hφ₀ hφ₀1 hU hs1 hs2 hR hlam)

end DPCBLean

