import Riemann.Deep.DPCBeyond.Lean.TheoremUB

/-!
# dpcbeyond: Theorem U′ ⇒ D(λ) > 0 for every λ > 7/11 (new step (f), part 2)

D_pos_of_kappa: Theorem U′ with slope κ gives 0 < D(λ) (Lamzouri's liminf, pplean's Dlam) for every
λ > max{2/(2+κ), 1/2}. theoremUB_conrey, D_pos_seven_elevenths: with Conrey's slope κ = 8/7 − δ for every
δ ∈ (0, 8/7) (DensityK and DPCInputsB at each such κ), 𝒩(T) ≫ T log T for every φ ∈ (7/11, 1) and 0 < D(λ)
for every λ > 7/11 (ParamsB.conrey_threshold picks δ).
-/

namespace DPCBLean

open PPLean DPCLean

variable {ex lg : ℝ → ℝ}

/-- Theorem U′ ⇒ D(λ) > 0 for every λ > max{2/(2+κ), 1/2}. -/
theorem D_pos_of_kappa {Z : Zeros} (hZ : Z.Ordered) {β : ℕ → ℝ} (hβ : Reflect Z β)
    (hEL : ExpLog ex lg) {tp : ℝ} (htp : 0 < tp) (hsc : ∀ T : ℝ, 0 < T → Z.sc T * tp = lg T)
    (hR : RvMUp Z) {κ : ℝ} (hκ : 0 < κ) (hden : DensityK Z β ex lg κ)
    (hin : DPCInputsB Z β ex lg κ) :
    ∀ lam : ℝ, max (2 / (2 + κ)) (1 / 2) < lam → 0 < Dlam Z lam := by
  intro lam hlam
  have h0 : 0 ≤ max (2 / (2 + κ)) (1 / 2) := le_trans (by norm_num) (le_max_right _ _)
  have h1 : max (2 / (2 + κ)) (1 / 2) < 1 :=
    max_lt ((div_lt_one (by linarith)).2 (by linarith)) (by norm_num)
  exact dlam_pos_of_from hZ h0 h1 (theoremUB hZ hβ hEL htp hsc hκ hden hin)
    (scaleToInf_of hEL htp hsc) (scaleHalf_of hEL htp hsc) hR hlam

/-- **Theorem U′ with Conrey's slope**: 𝒩(T) ≥ c(φ) T log T for every φ ∈ (7/11, 1), from the density and the
inputs at every κ = 8/7 − δ, δ ∈ (0, 8/7). -/
theorem theoremUB_conrey {Z : Zeros} (hZ : Z.Ordered) {β : ℕ → ℝ} (hβ : Reflect Z β)
    (hEL : ExpLog ex lg) {tp : ℝ} (htp : 0 < tp) (hsc : ∀ T : ℝ, 0 < T → Z.sc T * tp = lg T)
    (hden : ∀ δ : ℝ, 0 < δ → δ < 8 / 7 → DensityK Z β ex lg (8 / 7 - δ))
    (hin : ∀ δ : ℝ, 0 < δ → δ < 8 / 7 → DPCInputsB Z β ex lg (8 / 7 - δ)) : UFrom Z (7 / 11) := by
  intro φ hφ1 hφ2
  obtain ⟨δ, hδ0, hδ1, hth⟩ := conrey_threshold hφ1
  exact theoremUB hZ hβ hEL htp hsc (by linarith) (hden δ hδ0 hδ1) (hin δ hδ0 hδ1) φ hth hφ2

/-- **D(λ) > 0 for every λ > 7/11** (Lamzouri's liminf D), from the zero model, ExpLog, Riemann–von Mangoldt's
upper bound, Conrey's density slope 8/7 − δ and the named inputs of DPC's method beyond T^{1/2}. -/
theorem D_pos_seven_elevenths {Z : Zeros} (hZ : Z.Ordered) {β : ℕ → ℝ} (hβ : Reflect Z β)
    (hEL : ExpLog ex lg) {tp : ℝ} (htp : 0 < tp) (hsc : ∀ T : ℝ, 0 < T → Z.sc T * tp = lg T)
    (hR : RvMUp Z) (hden : ∀ δ : ℝ, 0 < δ → δ < 8 / 7 → DensityK Z β ex lg (8 / 7 - δ))
    (hin : ∀ δ : ℝ, 0 < δ → δ < 8 / 7 → DPCInputsB Z β ex lg (8 / 7 - δ)) :
    ∀ lam : ℝ, 7 / 11 < lam → 0 < Dlam Z lam := by
  intro lam hlam
  obtain ⟨δ, hδ0, hδ1, hth⟩ := conrey_threshold hlam
  exact D_pos_of_kappa hZ hβ hEL htp hsc hR (by linarith) (hden δ hδ0 hδ1) (hin δ hδ0 hδ1) lam hth

end DPCBLean

