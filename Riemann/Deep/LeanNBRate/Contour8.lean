import Riemann.Deep.LeanNBRate.Contour7

/-!
# U1 contour argument, part 8: the middle rectangle and the main identity G = F - R_tr - J (Re s > 3)
-/

open Complex Filter Topology Set MeasureTheory intervalIntegral
open scoped Real Interval

noncomputable section

namespace Riemann.Deep.LeanNBRate

variable {H : ℂ → ℂ} {M β D : ℝ} {N : ℕ}

/-- Hypotheses on the zeros: a good height T2 and removability of H/zeta at the zeros below it. -/
structure HZeros (H : ℂ → ℂ) (T₂ : ℝ) : Prop where
  T₂_pos : 0 < T₂
  good : ∀ σ ∈ Icc (-1 : ℝ) 2, riemannZeta (σ + T₂ * I) ≠ 0
  rem : ∀ ρ : ℂ, riemannZeta ρ = 0 → 0 ≤ ρ.re → ρ.re < 1 → |ρ.im| < T₂ →
    ∃ g : ℂ → ℂ, ContinuousAt g ρ ∧ ∀ᶠ z in 𝓝[≠] ρ, Fq H z = g z

theorem HZeros.good_neg {T₂ : ℝ} (hZ : HZeros H T₂) :
    ∀ σ ∈ Icc (-1 : ℝ) 2, riemannZeta (σ + ((-T₂ : ℝ) : ℂ) * I) ≠ 0 := by
  intro σ hσ h
  apply hZ.good σ hσ
  have e : (σ : ℂ) + T₂ * I = (starRingEnd ℂ) ((σ : ℂ) + ((-T₂ : ℝ) : ℂ) * I) := by
    apply Complex.ext <;> simp
  rw [e, riemannZeta_conj, h, map_zero]

theorem tendsto_inv_zeta_one : Tendsto (fun z => (riemannZeta z)⁻¹) (𝓝[≠] 1) (𝓝 0) := by
  have h1 := riemannZeta_residue_one
  have h2 : Tendsto (fun z : ℂ => z - 1) (𝓝[≠] (1 : ℂ)) (𝓝 0) :=
    tendsto_sub_nhds_zero_iff.mpr (tendsto_id.mono_left nhdsWithin_le_nhds)
  have h3 := h2.mul (h1.inv₀ one_ne_zero)
  rw [zero_mul] at h3
  refine h3.congr' ?_
  filter_upwards [self_mem_nhdsWithin] with z hz
  have hz' : z - 1 ≠ 0 := sub_ne_zero.mpr hz
  rw [mul_inv, ← mul_assoc, mul_inv_cancel₀ hz', one_mul]

theorem tendsto_Psi_one (hH : HBound H M β D N) (L : ℝ) {s : ℂ} (hs : 3 < s.re) :
    Tendsto (Psi H L s) (𝓝[≠] 1) (𝓝 0) := by
  have hs1 : s - 1 ≠ 0 := fun h => by have := congrArg Complex.re h; simp at this; linarith
  have hQ : ContinuousAt (fun z => H z * kap L (s - z) / (s - z)) 1 := by
    unfold kap
    exact ((hH.diff 1).continuousAt.mul (by fun_prop)).div (continuousAt_const.sub continuousAt_id) hs1
  have h := (hQ.tendsto.mono_left nhdsWithin_le_nhds).mul tendsto_inv_zeta_one
  rw [mul_zero] at h
  refine h.congr fun z => ?_
  simp only [Psi, Fq]
  ring

/-- The middle rectangle integral vanishes (removable singularities only). -/
theorem RectI_mid_zero (hH : HBound H M β D N) {T₂ : ℝ} (hZ : HZeros H T₂) (L : ℝ) {s : ℂ}
    (hs : 3 < s.re) : RectI (Psi H L s) ⟨-1, -T₂⟩ ⟨2, T₂⟩ = 0 := by
  classical
  have hT := hZ.T₂_pos
  set R : Set ℂ := [[(-1 : ℝ), 2]] ×ℂ [[-T₂, T₂]] with hR
  have hRc : IsCompact R := IsCompact.reProdIm isCompact_uIcc isCompact_uIcc
  have hfin := hRc.inter_riemannZetaZeros_finite
  set P : Finset ℂ := hfin.toFinset ∪ {1} with hP
  have hmemP : ∀ x, x ∈ P ↔ (x ∈ R ∧ riemannZeta x = 0) ∨ x = 1 := by
    intro x
    simp only [hP, Finset.mem_union, Set.Finite.mem_toFinset, Set.mem_inter_iff, mem_riemannZetaZeros,
      Finset.mem_singleton]
  have hzeroR : ∀ x ∈ R, riemannZeta x = 0 → 0 ≤ x.re ∧ x.re < 1 ∧ |x.im| < T₂ := by
    intro x hx hz
    rw [hR, Complex.mem_reProdIm, Set.uIcc_of_le (by norm_num), Set.uIcc_of_le (by linarith)] at hx
    have hre0 : 0 ≤ x.re := by
      by_contra hneg; push Not at hneg
      obtain ⟨k, hk⟩ := zeta_eq_zero_left hz hneg
      have := congrArg Complex.re hk; simp at this
      linarith [hx.1.1, (Nat.cast_nonneg k : (0:ℝ) ≤ k)]
    have hre1 : x.re < 1 := by
      by_contra hge; push Not at hge
      exact riemannZeta_ne_zero_of_one_le_re hge hz
    refine ⟨hre0, hre1, ?_⟩
    rw [abs_lt]
    have hxe : x = (x.re : ℂ) + x.im * I := (Complex.re_add_im x).symm
    have hσ : x.re ∈ Icc (-1 : ℝ) 2 := ⟨by linarith, by linarith⟩
    constructor
    · rcases lt_or_eq_of_le hx.2.1 with h | h
      · exact h
      · exfalso; apply hZ.good_neg x.re hσ; rw [h]; rw [hxe] at hz; simpa [h] using hz
    · rcases lt_or_eq_of_le hx.2.2 with h | h
      · exact h
      · exfalso; apply hZ.good x.re hσ; rw [hxe] at hz; simpa [h] using hz
  have hRmem : ∀ x : ℂ, x ∈ R ↔ -1 ≤ x.re ∧ x.re ≤ 2 ∧ -T₂ ≤ x.im ∧ x.im ≤ T₂ := by
    intro x
    rw [hR, Complex.mem_reProdIm, Set.uIcc_of_le (by norm_num), Set.uIcc_of_le (by linarith)]
    simp only [mem_Icc]; tauto
  have h := RectI_eq_sum_residues (f := Psi H L s) (z := ⟨-1, -T₂⟩) (w := ⟨2, T₂⟩) (by norm_num)
    (by simp; linarith) P (fun _ => 0) ?_ ?_ ?_ ?_
  · simpa using h
  · intro p hp
    rcases (hmemP p).mp hp with ⟨hpR, hpz⟩ | rfl
    · obtain ⟨a1, a2, a3⟩ := hzeroR p hpR hpz
      rw [abs_lt] at a3
      exact ⟨by simp; linarith, by simp; linarith, by simp; linarith, by simp; linarith⟩
    · exact ⟨by norm_num, by norm_num, by simp; linarith, by simp; linarith⟩
  · intro x hx hxP
    have hxR : x ∈ R := hx
    have hx1 : x ≠ 1 := fun h => hxP ((hmemP x).mpr (Or.inr h))
    have hx0 : riemannZeta x ≠ 0 := fun h => hxP ((hmemP x).mpr (Or.inl ⟨hxR, h⟩))
    have hxs : x ≠ s := fun h => by
      rw [hRmem] at hxR; rw [h] at hxR; linarith [hxR.2.1]
    exact continuousAt_Psi hH L hx1 hx0 hxs
  · intro x hx hxP
    have hxR : x ∈ R := by
      rw [Complex.mem_reProdIm] at hx; rw [hRmem]; simp at hx
      exact ⟨hx.1.1.le, hx.1.2.le, hx.2.1.le, hx.2.2.le⟩
    have hx1 : x ≠ 1 := fun h => hxP ((hmemP x).mpr (Or.inr h))
    have hx0 : riemannZeta x ≠ 0 := fun h => hxP ((hmemP x).mpr (Or.inl ⟨hxR, h⟩))
    have hxs : x ≠ s := fun h => by
      rw [hRmem] at hxR; rw [h] at hxR; linarith [hxR.2.1]
    exact differentiableAt_Psi hH L hx1 hx0 hxs
  · intro p hp
    rcases (hmemP p).mp hp with ⟨hpR, hpz⟩ | rfl
    · obtain ⟨a1, a2, a3⟩ := hzeroR p hpR hpz
      obtain ⟨g₀, hg₀c, hg₀e⟩ := hZ.rem p hpz a1 a2 a3
      have hps : s - p ≠ 0 := fun h => by
        have := congrArg Complex.re h; simp at this; linarith
      refine ⟨fun z => g₀ z * kap L (s - z) / (s - z), ?_, ?_⟩
      · unfold kap
        exact (hg₀c.mul (by fun_prop)).div (continuousAt_const.sub continuousAt_id) hps
      · filter_upwards [hg₀e] with z hz
        simp only [Psi, hz, zero_div, zero_add]
    · refine ⟨Function.update (Psi H L s) 1 0, continuousAt_update_same.mpr (tendsto_Psi_one hH L hs), ?_⟩
      filter_upwards [self_mem_nhdsWithin] with z hz
      rw [Function.update_of_ne hz, zero_div, zero_add]

/-- **The main identity** (U1): for Re s > 3, G(s) = F(s) - R_tr(s) - J(s). -/
theorem GF_eq_main (hH : HBound H M β D N) {T₂ : ℝ} (hZ : HZeros H T₂) {L : ℝ} (hL : 0 ≤ L)
    (hLβ : β + Real.log (2 * π) < L) {s : ℂ} (hs : 3 < s.re) :
    GF H L s = Fq H s - Rtr H L s - Jf H L T₂ s := by
  have h1 := GF_eq hH hL hs
  have h2 := integral_Psi0_line hH hs
  have h3 := VI_Psi_three_two hH hL hs
  have h4 := VI_sub_VI_eq (f := Psi H L s) (a := -1) (b := 2) T₂ (integrable_Psi_neg_one hH hL hs)
    (by simpa using integrable_Psi_line_right hH hL hs (σ := 2) (by norm_num) (by norm_num))
  rw [RectI_mid_zero hH hZ L hs, zero_add] at h4
  have h5 := VI_Psi_neg_one hH hL hs
  have h6 := VI_Psi_hl0 hH hL hs hLβ
  have e3 : VI (Psi H L s) 3 = ∫ y : ℝ, Psi H L s (3 + y * I) := by simp [VI]
  have hI : (I : ℂ) ≠ 0 := Complex.I_ne_zero
  have key : I * VI (Psi H L s) 3 = 2 * π * I * Rtr H L s + (UpperU (Psi H L s) (-1) 2 T₂ -
      LowerU (Psi H L s) (-1) 2 T₂) := by
    linear_combination h3 + h4 + h5 + h6
  rw [h1, h2, ← e3]
  have hpi : (π : ℂ) ≠ 0 := by exact_mod_cast Real.pi_ne_zero
  have hV : VI (Psi H L s) 3 = 2 * π * Rtr H L s + (UpperU (Psi H L s) (-1) 2 T₂ -
      LowerU (Psi H L s) (-1) 2 T₂) / I := by
    field_simp
    linear_combination key
  rw [hV, Jf]
  field_simp
  ring

end Riemann.Deep.LeanNBRate

