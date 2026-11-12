import Riemann.Deep.DPCUncond.Lean.TheoremU
import Riemann.Deep.ExplorePointProc.Lean.Liminf

/-!
# dpcuncond (new step (f)): Theorem U ⇒ D(λ) > 0 for every λ > 2/3

pplean's zero model and Lamzouri's D(λ) = liminf #smallGaps(λ, T)/N(T) (TeX l.183; PPLean.Dlam, PPLean.dpos_iff,
imported from Riemann/Deep/ExplorePointProc/Lean). dpos_of_U: choose 2/3 < φ < min(λ, 1) and apply Theorem U on
[T/2, T]: for large T, 2πφ/log(T/2) < 2πλ/log T (ScaleHalf), so 𝒩_φ(T/2) ≤ #smallGaps(λ, T), and
𝒩_φ(T/2) ≥ c(T/2) log(T/2)/(2π) ≫ N(T) (RvMUp). D_pos_two_thirds: the whole chain from the named inputs.
-/

namespace DPCLean

open Finset PPLean

/-- sc T − sc(T/2) is bounded (= log 2/(2π)). -/
def ScaleHalf (Z : Zeros) : Prop := ∃ k : ℝ, ∀ T : ℝ, 0 < T → Z.sc T ≤ Z.sc (T / 2) + k

/-- Riemann–von Mangoldt, upper half: N(T) ≤ K T log T/(2π) for large T (Lz l.177). -/
def RvMUp (Z : Zeros) : Prop := ∃ K T₀ : ℝ, ∀ T : ℝ, T₀ ≤ T → (Z.N T : ℝ) ≤ K * T * Z.sc T

lemma scaleToInf_of {Z : Zeros} {ex lg : ℝ → ℝ} (hEL : ExpLog ex lg) {tp : ℝ} (htp : 0 < tp)
    (hsc : ∀ T : ℝ, 0 < T → Z.sc T * tp = lg T) : ScaleToInf Z := by
  intro K
  refine ⟨ex (K * tp), fun T hT => ?_⟩
  have hT0 : 0 < T := lt_of_lt_of_le (hEL.pos _) hT
  have h1 := hEL.le_lg hT
  rw [← hsc T hT0] at h1
  exact le_of_mul_le_mul_right h1 htp

lemma scaleHalf_of {Z : Zeros} {ex lg : ℝ → ℝ} (hEL : ExpLog ex lg) {tp : ℝ} (htp : 0 < tp)
    (hsc : ∀ T : ℝ, 0 < T → Z.sc T * tp = lg T) : ScaleHalf Z := by
  refine ⟨lg 2 / tp, fun T hT => ?_⟩
  have hT2 : 0 < T / 2 := by linarith
  have h1 : lg T = lg 2 + lg (T / 2) := by
    rw [← hEL.lg_mul (show (0 : ℝ) < 2 by norm_num) hT2]
    congr 1
    ring
  have h2 := hsc T hT
  have h3 := hsc (T / 2) hT2
  have h4 : Z.sc T * tp = (Z.sc (T / 2) + lg 2 / tp) * tp := by
    rw [h2, h1, ← h3, add_mul, div_mul_cancel₀ _ htp.ne']
    ring
  exact le_of_eq (mul_right_cancel₀ htp.ne' h4)

/-- Gaps of 𝒩_φ(T/2) are small gaps at height T when φ sc T ≤ λ sc(T/2). -/
lemma gapSet_half_sub {Z : Zeros} (hZ : Z.Ordered) {φ lam T : ℝ} (hφ : 0 < φ) (hlam : 0 < lam)
    (hscale : φ * Z.sc T ≤ lam * Z.sc (T / 2)) : gapSet Z φ (T / 2) ⊆ smallGaps Z lam T := by
  intro n hn
  simp only [gapSet, mem_filter, mem_Icc] at hn
  obtain ⟨⟨hn1, hnN⟩, -, -, hgap⟩ := hn
  rw [show 2 * (T / 2) = T by ring] at hnN
  simp only [smallGaps, mem_filter, mem_Icc]
  refine ⟨⟨hn1, hnN⟩, ?_⟩
  have hg : 0 ≤ Z.γ (n + 1) - Z.γ n := sub_nonneg.2 (hZ.mono n hn1)
  have h1 : φ * ((Z.γ (n + 1) - Z.γ n) * Z.sc T) ≤ φ * lam := by
    calc φ * ((Z.γ (n + 1) - Z.γ n) * Z.sc T) = (Z.γ (n + 1) - Z.γ n) * (φ * Z.sc T) := by ring
      _ ≤ (Z.γ (n + 1) - Z.γ n) * (lam * Z.sc (T / 2)) := mul_le_mul_of_nonneg_left hscale hg
      _ = lam * ((Z.γ (n + 1) - Z.γ n) * Z.sc (T / 2)) := by ring
      _ ≤ lam * φ := mul_le_mul_of_nonneg_left hgap.le hlam.le
      _ = φ * lam := by ring
  exact le_of_mul_le_mul_left h1 hφ

/-- **Theorem U ⇒ D(λ) > 0** (DPos form) for every λ > 2/3. -/
theorem dpos_of_U {Z : Zeros} (hZ : Z.Ordered) (hU : TheoremU Z) (hs1 : ScaleToInf Z)
    (hs2 : ScaleHalf Z) (hR : RvMUp Z) {lam : ℝ} (hlam : 2 / 3 < lam) : DPos Z lam := by
  obtain ⟨φ, hφ⟩ : ∃ φ : ℝ, φ = (2 / 3 + min lam 1) / 2 := ⟨_, rfl⟩
  have hm1 : min lam 1 ≤ lam := min_le_left _ _
  have hm2 : min lam 1 ≤ 1 := min_le_right _ _
  have hm3 : 2 / 3 < min lam 1 := lt_min hlam (by norm_num)
  have hφ1 : 2 / 3 < φ := by rw [hφ]; linarith
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

/-- N(T) → ∞, from Theorem U itself (N(T) ≥ 𝒩_{5/6}(T/2)). -/
lemma nToInf_of_U {Z : Zeros} (hU : TheoremU Z) (hs1 : ScaleToInf Z) :
    ∀ K : ℝ, ∃ T₀ : ℝ, ∀ T : ℝ, T₀ ≤ T → K ≤ (Z.N T : ℝ) := by
  intro K
  obtain ⟨c, hc, T₁, hT₁⟩ := hU (5 / 6) (by norm_num) (by norm_num)
  obtain ⟨T₃, hT₃⟩ := hs1 1
  refine ⟨max (max (2 * T₁) (2 * T₃)) (2 * |K| / c), fun T hT => ?_⟩
  simp only [max_le_iff] at hT
  obtain ⟨⟨h1, h3⟩, h4⟩ := hT
  have hs := hT₃ (T / 2) (by linarith)
  have hU' := hT₁ (T / 2) (by linarith)
  have hcard : ((gapSet Z (5 / 6) (T / 2)).card : ℝ) ≤ (Z.N T : ℝ) := by
    have h := card_gapSet_le Z (5 / 6) (T / 2)
    rw [show 2 * (T / 2) = T by ring] at h
    exact_mod_cast h
  have hK0 : 0 ≤ 2 * |K| / c := by positivity
  have h5 : c * (T / 2) ≤ c * (T / 2) * Z.sc (T / 2) :=
    le_mul_of_one_le_right (by nlinarith) hs
  have h6 : |K| ≤ c * (T / 2) := by
    rw [div_le_iff₀ hc] at h4
    linarith
  linarith [le_abs_self K]

/-- Theorem U ⇒ 0 < D(λ) (Lamzouri's liminf) for every λ > 2/3. -/
theorem dlam_pos_of_U {Z : Zeros} (hZ : Z.Ordered) (hU : TheoremU Z) (hs1 : ScaleToInf Z)
    (hs2 : ScaleHalf Z) (hR : RvMUp Z) {lam : ℝ} (hlam : 2 / 3 < lam) : 0 < Dlam Z lam :=
  (dpos_iff (nToInf_of_U hU hs1) lam).1 (dpos_of_U hZ hU hs1 hs2 hR hlam)

/-- **D(λ) > 0 for every λ > 2/3** (Lamzouri's liminf D), from the zero model, ExpLog, Riemann–von Mangoldt's
upper bound and the named analytic inputs of DPC's method (DPCInputs). -/
theorem D_pos_two_thirds {Z : Zeros} (hZ : Z.Ordered) {β : ℕ → ℝ} (hβ : Reflect Z β) {ex lg : ℝ → ℝ}
    (hEL : ExpLog ex lg) {tp : ℝ} (htp : 0 < tp) (hsc : ∀ T : ℝ, 0 < T → Z.sc T * tp = lg T)
    (hR : RvMUp Z) (hin : DPCInputs Z β ex lg) : ∀ lam : ℝ, 2 / 3 < lam → 0 < Dlam Z lam :=
  fun _ hlam => dlam_pos_of_U hZ (theoremU hZ hβ hEL htp hsc hin) (scaleToInf_of hEL htp hsc)
    (scaleHalf_of hEL htp hsc) hR hlam

end DPCLean

