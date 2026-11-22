import Riemann.Deep.LeanNBRate.KBd
import Riemann.Deep.LeanNBRate.RtrBd

/-!
# The formula on Re s > 1 and the bound for G(1) (the Nyman constant c₀)
-/

open Complex Filter Topology Set MeasureTheory
open scoped Real

noncomputable section

namespace Riemann.Deep.LeanNBRate

variable {H : ℂ → ℂ} {M β D : ℝ} {N : ℕ}

/-- The formula ζG = H - K - ζR_tr on Re s > 1 (identity theorem). -/
theorem formula_right (hH : HBound H M β D N) {T₂ : ℝ} (hZ : HZeros H T₂) (hT1 : 1 ≤ T₂) {L : ℝ}
    (hL : 0 ≤ L) (hLβ : β + Real.log (2 * π) < L) {s : ℂ} (hs : 1 < s.re) :
    riemannZeta s * GF H L s = H s - Kf H L T₂ s - riemannZeta s * Rtr H L s := by
  set Φ : ℂ → ℂ := fun s => riemannZeta s * GF H L s + riemannZeta s * Rtr H L s + Kf H L T₂ s - H s with hΦ
  have hU₁o : IsOpen {s : ℂ | 1 < s.re} := isOpen_lt continuous_const Complex.continuous_re
  have hΦd : DifferentiableOn ℂ Φ {s : ℂ | 1 < s.re} := by
    intro z hz
    have hz' : 1 < z.re := hz
    have hz1 : z ≠ 1 := fun h => by rw [h] at hz'; simp at hz'
    refine DifferentiableAt.differentiableWithinAt ?_
    have hR : DifferentiableAt ℂ (Rtr H L) z :=
      (differentiableOn_Rtr hH hL).differentiableAt
        ((isOpen_lt continuous_const Complex.continuous_re).mem_nhds (show -1 < z.re by linarith))
    exact ((((differentiableAt_riemannZeta hz1).mul (differentiable_GF hH hL z)).add
      ((differentiableAt_riemannZeta hz1).mul hR)).add (differentiableAt_Kf hH hZ hT1 hL hz1)).sub (hH.diff z)
  have hΦ3 : ∀ z : ℂ, 3 < z.re → Φ z = 0 := by
    intro z hz
    have hG := GF_eq_main hH hZ hL hLβ hz
    have hK := Kf_eq_zeta_Jf hH hZ hT1 hL hz
    have hz0 : riemannZeta z ≠ 0 := riemannZeta_ne_zero_of_one_lt_re (by linarith)
    have hF : riemannZeta z * Fq H z = H z := by rw [Fq]; field_simp
    simp only [hΦ, hG, hK]
    linear_combination hF
  have h1 : EqOn Φ (fun _ => 0) {s : ℂ | 1 < s.re} := by
    refine (hΦd.analyticOnNhd hU₁o).eqOn_of_preconnected_of_eventuallyEq analyticOnNhd_const
      (convex_halfSpace_re_gt 1).isPreconnected (z₀ := 4) (by show (1 : ℝ) < (4 : ℂ).re; norm_num) ?_
    filter_upwards [(isOpen_lt continuous_const Complex.continuous_re).mem_nhds
      (show (3 : ℝ) < (4 : ℂ).re by norm_num)] with z hz
    exact hΦ3 z hz
  have := h1 hs
  simp only [hΦ] at this
  linear_combination this

theorem tendsto_inv_zeta_right :
    Tendsto (fun σ : ℝ => (riemannZeta σ)⁻¹) (𝓝[>] 1) (𝓝 0) := by
  have h1 : Tendsto (fun s : ℂ => (riemannZeta s)⁻¹) (𝓝[≠] 1) (𝓝 0) := by
    have hA : Tendsto (fun s : ℂ => s - 1) (𝓝[≠] 1) (𝓝 0) := by
      have : Tendsto (fun s : ℂ => s - 1) (𝓝 1) (𝓝 (1 - 1)) :=
        ((continuous_id.sub continuous_const).tendsto 1)
      rw [sub_self] at this
      exact this.mono_left nhdsWithin_le_nhds
    have hB := riemannZeta_residue_one.inv₀ one_ne_zero
    have := hA.mul hB
    rw [zero_mul] at this
    refine this.congr' ?_
    filter_upwards [self_mem_nhdsWithin] with s hs
    have hs' : s - 1 ≠ 0 := sub_ne_zero.mpr hs
    rw [mul_inv, ← mul_assoc, mul_inv_cancel₀ hs', one_mul]
  have hc : Tendsto (fun σ : ℝ => (σ : ℂ)) (𝓝[>] (1 : ℝ)) (𝓝 (1 : ℂ)) := by
    have := (Complex.continuous_ofReal.tendsto (1 : ℝ)).mono_left (nhdsWithin_le_nhds (s := Ioi (1 : ℝ)))
    simpa using this
  have hne : ∀ᶠ σ : ℝ in 𝓝[>] 1, (σ : ℂ) ∈ ({1}ᶜ : Set ℂ) := by
    filter_upwards [self_mem_nhdsWithin] with σ hσ
    simp only [mem_compl_iff, mem_singleton_iff]
    intro h
    have h' := congrArg Complex.re h
    simp at h'
    have : (1 : ℝ) < σ := hσ
    linarith
  exact h1.comp (tendsto_nhdsWithin_iff.mpr ⟨hc, hne⟩)

/-- **The Nyman constant**: |G(1)| is bounded by the contour mass of H and the R_tr bound. -/
theorem norm_GF_one_le (hH : HBound H M β D N) {T₂ : ℝ} (hZ : HZeros H T₂) (hT1 : 1 ≤ T₂) {L : ℝ}
    (hL : 0 ≤ L) (hLβ : β + Real.log (2 * π) < L) {B : ℝ} (hB0 : 0 ≤ B)
    (hBu : ∀ x ∈ Icc (-1 : ℝ) 2, ‖(riemannZeta (x + T₂ * I))⁻¹‖ ≤ B)
    (hBl : ∀ x ∈ Icc (-1 : ℝ) 2, ‖(riemannZeta (x + ((-T₂ : ℝ) : ℂ) * I))⁻¹‖ ≤ B) {R₁ : ℝ}
    (hR : ∀ σ : ℝ, 1 < σ → ‖Rtr H L σ‖ ≤ R₁) :
    ‖GF H L 1‖ ≤ (2 * π)⁻¹ * (Real.exp (L * (3 / 2)) * Bc B) * (ZU H T₂ + ZL H T₂) + R₁ := by
  have hid : ∀ σ : ℝ, 1 < σ → GF H L σ =
      (riemannZeta σ)⁻¹ * H σ - (riemannZeta σ)⁻¹ * Kf H L T₂ σ - Rtr H L σ := by
    intro σ hσ
    have hz : riemannZeta σ ≠ 0 := riemannZeta_ne_zero_of_one_lt_re (by simpa using hσ)
    have h := formula_right hH hZ hT1 hL hLβ (s := σ) (by simpa using hσ)
    have e : (riemannZeta σ)⁻¹ * riemannZeta σ = 1 := inv_mul_cancel₀ hz
    linear_combination (riemannZeta σ)⁻¹ * h - (GF H L σ + Rtr H L σ) * e
  obtain ⟨Z, hZdef⟩ : ∃ Z : ℝ, Z = ZU H T₂ + ZL H T₂ := ⟨_, rfl⟩
  obtain ⟨E, hEdef⟩ : ∃ E : ℝ, E = Real.exp (L * (3 / 2)) := ⟨_, rfl⟩
  rw [← hZdef, ← hEdef]
  have hbd : ∀ σ : ℝ, 1 < σ → ‖GF H L σ‖ ≤
      ‖H σ‖ * ‖(riemannZeta σ)⁻¹‖ + (2 * π)⁻¹ * (E * (Bc B + ‖(riemannZeta σ)⁻¹‖)) * Z + R₁ := by
    intro σ hσ
    rw [hid σ hσ]
    have hK := norm_Kf_div_zeta hH hT1 hL hB0 hBu hBl (σ := σ) hσ.le
    rw [← hZdef, ← hEdef] at hK
    have h1 := norm_sub_le ((riemannZeta σ)⁻¹ * H σ - (riemannZeta σ)⁻¹ * Kf H L T₂ σ) (Rtr H L σ)
    have h2 := norm_sub_le ((riemannZeta σ)⁻¹ * H σ) ((riemannZeta σ)⁻¹ * Kf H L T₂ σ)
    rw [norm_mul] at h2
    have := hR σ hσ
    nlinarith [norm_nonneg (H σ), norm_nonneg (riemannZeta σ)⁻¹]
  have hGc : Tendsto (fun σ : ℝ => ‖GF H L σ‖) (𝓝[>] 1) (𝓝 ‖GF H L 1‖) := by
    have h := (((differentiable_GF hH hL).continuous.comp Complex.continuous_ofReal).norm.tendsto 1).mono_left
      (nhdsWithin_le_nhds (s := Ioi (1 : ℝ)))
    simpa using h
  have hHc : Tendsto (fun σ : ℝ => ‖H σ‖) (𝓝[>] 1) (𝓝 ‖H 1‖) := by
    have h := ((hH.diff.continuous.comp Complex.continuous_ofReal).norm.tendsto 1).mono_left
      (nhdsWithin_le_nhds (s := Ioi (1 : ℝ)))
    simpa using h
  have hinv : Tendsto (fun σ : ℝ => ‖(riemannZeta σ)⁻¹‖) (𝓝[>] 1) (𝓝 0) := by
    have := tendsto_inv_zeta_right.norm
    simpa using this
  have hlim : Tendsto (fun σ : ℝ => ‖H σ‖ * ‖(riemannZeta σ)⁻¹‖ +
      (2 * π)⁻¹ * (E * (Bc B + ‖(riemannZeta σ)⁻¹‖)) * Z + R₁) (𝓝[>] 1)
      (𝓝 (‖H 1‖ * 0 + (2 * π)⁻¹ * (E * (Bc B + 0)) * Z + R₁)) :=
    ((hHc.mul hinv).add (((tendsto_const_nhds.add hinv).const_mul E).const_mul (2 * π)⁻¹ |>.mul_const Z)).add
      tendsto_const_nhds
  have hev : ∀ᶠ σ : ℝ in 𝓝[>] 1, ‖GF H L σ‖ ≤
      ‖H σ‖ * ‖(riemannZeta σ)⁻¹‖ + (2 * π)⁻¹ * (E * (Bc B + ‖(riemannZeta σ)⁻¹‖)) * Z + R₁ := by
    filter_upwards [self_mem_nhdsWithin] with σ hσ using hbd σ hσ
  have := le_of_tendsto_of_tendsto hGc hlim hev
  simpa using this

end Riemann.Deep.LeanNBRate

