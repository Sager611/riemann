import Riemann.Deep.LeanNBRate.ContourBd

/-!
# K on the critical line and K(σ)/ζ(σ) for real σ ≥ 1: bounds by the contour mass of H
-/

open Complex Filter Topology Set MeasureTheory Metric
open scoped Real

noncomputable section

namespace Riemann.Deep.LeanNBRate

variable {H : ℂ → ℂ} {M β D : ℝ} {N : ℕ}

theorem re_mem_of_near {s z : ℂ} (hs : s.re = 1 / 2) (hz : z ∈ closedBall s (1 / 8)) :
    3 / 8 ≤ z.re ∧ z.re ≤ 5 / 8 := by
  have h := Complex.abs_re_le_norm (z - s)
  rw [mem_closedBall, dist_eq_norm] at hz
  rw [sub_re, hs] at h
  have h2 := abs_le.mp (h.trans hz)
  constructor <;> linarith [h2.1, h2.2]

theorem norm_mul_norm_inv_le (z : ℂ) : ‖z‖ * ‖z⁻¹‖ ≤ 1 := by
  rw [← norm_mul]
  by_cases h0 : z = 0
  · simp [h0]
  · rw [mul_inv_cancel₀ h0, norm_one]

/-- dslope bound near the critical line: far away by the triangle inequality, near by Cauchy. -/
theorem norm_dslope_div_le (h215 : Titchmarsh_2_1_5) {s w : ℂ} (hs : s.re = 1 / 2) {Bw R : ℝ}
    (hBw0 : 0 ≤ Bw) (hR : 0 ≤ R) (hBw : ‖(riemannZeta w)⁻¹‖ ≤ Bw) (hnear : ‖s - w‖ < 1 / 8 → ‖s‖ ≤ R) :
    ‖dslope riemannZeta w s / riemannZeta w‖ ≤ (8 + Bw * (64 * R + 24)) * (1 + ‖riemannZeta s‖) := by
  have hzs := norm_nonneg (riemannZeta s)
  have hiw := norm_nonneg (riemannZeta w)⁻¹
  rw [div_eq_mul_inv, norm_mul]
  by_cases hfar : 1 / 8 ≤ ‖s - w‖
  · have hsw : s ≠ w := fun h => by rw [h, sub_self, norm_zero] at hfar; norm_num at hfar
    rw [dslope_of_ne _ hsw, slope_def_field]
    have h1 : ‖(riemannZeta s - riemannZeta w) / (s - w)‖ ≤ 8 * (‖riemannZeta s‖ + ‖riemannZeta w‖) := by
      rw [norm_div, div_le_iff₀ (by linarith)]
      have := norm_sub_le (riemannZeta s) (riemannZeta w)
      nlinarith [norm_nonneg (riemannZeta s - riemannZeta w), norm_nonneg (riemannZeta w)]
    have h2 := norm_mul_norm_inv_le (riemannZeta w)
    calc ‖(riemannZeta s - riemannZeta w) / (s - w)‖ * ‖(riemannZeta w)⁻¹‖
        ≤ 8 * (‖riemannZeta s‖ + ‖riemannZeta w‖) * ‖(riemannZeta w)⁻¹‖ :=
          mul_le_mul_of_nonneg_right h1 hiw
      _ = 8 * (‖riemannZeta s‖ * ‖(riemannZeta w)⁻¹‖ + ‖riemannZeta w‖ * ‖(riemannZeta w)⁻¹‖) := by ring
      _ ≤ 8 * (‖riemannZeta s‖ * Bw + 1) := by
          have := mul_le_mul_of_nonneg_left hBw hzs
          linarith
      _ ≤ _ := by
          nlinarith [mul_nonneg (mul_nonneg hBw0 hR) hzs, mul_nonneg hBw0 hzs, mul_nonneg hBw0 hR]
  · push Not at hfar
    have hsR := hnear hfar
    have hsS : s ∈ closedBall s (1 / 8) := mem_closedBall_self (by norm_num)
    have hwS : w ∈ closedBall s (1 / 8) := by
      rw [mem_closedBall, dist_eq_norm, ← norm_neg, neg_sub]; linarith
    have hdiff : ∀ z ∈ closedBall s (1 / 8), DifferentiableAt ℂ riemannZeta z := fun z hz =>
      differentiableAt_riemannZeta (fun h => by
        have := (re_mem_of_near hs hz).2; rw [h] at this; norm_num at this)
    have hbd : ∀ z ∈ closedBall s (1 / 8), ‖deriv riemannZeta z‖ ≤ 64 * R + 16 := by
      intro z hz
      have hm := re_mem_of_near hs hz
      have h1 := norm_deriv_zeta_le h215 hm.1 hm.2
      have h2 : ‖z‖ ≤ ‖s‖ + 1 / 8 := by
        have h3 := hz; rw [mem_closedBall, dist_eq_norm] at h3
        have := norm_sub_norm_le z s
        linarith
      linarith
    have hdsl : ‖dslope riemannZeta w s‖ ≤ 64 * R + 16 := by
      by_cases hsw : s = w
      · rw [hsw, dslope_same]; exact hbd w hwS
      · rw [dslope_of_ne _ hsw, slope_def_field, norm_div]
        have := Convex.norm_image_sub_le_of_norm_deriv_le hdiff hbd (convex_closedBall s (1 / 8)) hwS hsS
        have hpos : 0 < ‖s - w‖ := norm_pos_iff.mpr (sub_ne_zero.mpr hsw)
        rw [div_le_iff₀ hpos]; exact this
    have h0 : 0 ≤ 64 * R + 16 := by linarith
    calc ‖dslope riemannZeta w s‖ * ‖(riemannZeta w)⁻¹‖ ≤ (64 * R + 16) * Bw :=
          mul_le_mul hdsl hBw hiw h0
      _ ≤ _ := by
          nlinarith [mul_nonneg (mul_nonneg hBw0 hR) hzs, mul_nonneg hBw0 hzs, mul_nonneg hBw0 hR]

theorem norm_kap_le_three_halves {L : ℝ} (hL : 0 ≤ L) {s w : ℂ} (hs : 1 / 2 ≤ s.re) (hw : w.re ≤ 2) :
    ‖kap L (s - w)‖ ≤ Real.exp (L * (3 / 2)) := by
  rw [norm_kap, Real.exp_le_exp, sub_re]; nlinarith

theorem norm_kint_le (h215 : Titchmarsh_2_1_5) {L : ℝ} (hL : 0 ≤ L) {s w : ℂ} (hs : s.re = 1 / 2)
    (hw2 : w.re ≤ 2) {Bw R : ℝ} (hBw0 : 0 ≤ Bw) (hR : 0 ≤ R) (hBw : ‖(riemannZeta w)⁻¹‖ ≤ Bw)
    (hnear : ‖s - w‖ < 1 / 8 → ‖s‖ ≤ R) :
    ‖kint H L s w‖ ≤ (Real.exp (L * (3 / 2)) * (8 + Bw * (64 * R + 24)) * (1 + ‖riemannZeta s‖)) * ‖H w‖ := by
  have hk := norm_kap_le_three_halves hL (s := s) (w := w) (by rw [hs]) hw2
  have hd := norm_dslope_div_le h215 hs hBw0 hR hBw hnear
  rw [kint, mul_div_assoc, norm_mul, norm_mul]
  have h0 : 0 ≤ (8 + Bw * (64 * R + 24)) * (1 + ‖riemannZeta s‖) := by positivity
  calc ‖H w‖ * ‖kap L (s - w)‖ * ‖dslope riemannZeta w s / riemannZeta w‖
      ≤ ‖H w‖ * Real.exp (L * (3 / 2)) * ((8 + Bw * (64 * R + 24)) * (1 + ‖riemannZeta s‖)) := by
        gcongr
    _ = _ := by ring

theorem norm_two_pi_I_inv : ‖(2 * (π : ℂ) * I)⁻¹‖ = (2 * π)⁻¹ := by
  rw [norm_inv, norm_mul, norm_mul, Complex.norm_I, mul_one, Complex.norm_real,
    Real.norm_of_nonneg Real.pi_pos.le]
  norm_num

/-- The 1/ζ bound on the whole contour. -/
def Bc (B : ℝ) : ℝ := B + 3 + 60 * π * Real.sqrt (2 * π)

theorem Bc_nonneg {B : ℝ} (hB : 0 ≤ B) : 0 ≤ Bc B := by unfold Bc; positivity

theorem le_Bc {B x : ℝ} (h : x ≤ B) : x ≤ Bc B := by
  unfold Bc; have : 0 ≤ 60 * π * Real.sqrt (2 * π) := by positivity
  linarith

theorem inv_zeta_right_le {B : ℝ} (hB : 0 ≤ B) (y : ℝ) : ‖(riemannZeta (((2 : ℝ) : ℂ) + y * I))⁻¹‖ ≤ Bc B := by
  have h := norm_inv_zeta_le_three (s := ((2 : ℝ) : ℂ) + y * I) (by simp)
  unfold Bc; have : 0 ≤ 60 * π * Real.sqrt (2 * π) := by positivity
  linarith

theorem inv_zeta_left_le {B : ℝ} (hB : 0 ≤ B) {y : ℝ} (hy : 1 ≤ |y|) :
    ‖(riemannZeta (((-1 : ℝ) : ℂ) + y * I))⁻¹‖ ≤ Bc B := by
  have h := norm_inv_zeta_neg_one_line hy
  unfold Bc; linarith

/-- **K on the critical line** is bounded by the contour mass of H times (1 + |ζ|). -/
theorem norm_Kf_line (h215 : Titchmarsh_2_1_5) (hH : HBound H M β D N) {T₂ : ℝ} (hT1 : 1 ≤ T₂) {L : ℝ}
    (hL : 0 ≤ L) {B : ℝ} (hB0 : 0 ≤ B)
    (hBu : ∀ x ∈ Icc (-1 : ℝ) 2, ‖(riemannZeta (x + T₂ * I))⁻¹‖ ≤ B)
    (hBl : ∀ x ∈ Icc (-1 : ℝ) 2, ‖(riemannZeta (x + ((-T₂ : ℝ) : ℂ) * I))⁻¹‖ ≤ B) (t : ℝ) :
    ‖Kf H L T₂ (1 / 2 + t * I)‖ ≤ (2 * π)⁻¹ * (Real.exp (L * (3 / 2)) * (8 + Bc B * (64 * (T₂ + 3) + 24)) *
      (1 + ‖riemannZeta (1 / 2 + t * I)‖)) * (ZU H T₂ + ZL H T₂) := by
  obtain ⟨s, hsdef⟩ : ∃ s : ℂ, s = 1 / 2 + t * I := ⟨_, rfl⟩
  rw [← hsdef]
  have hs : s.re = 1 / 2 := by rw [hsdef]; simp
  have hBc0 := Bc_nonneg hB0
  have hR0 : (0 : ℝ) ≤ T₂ + 3 := by linarith
  obtain ⟨C, hCdef⟩ : ∃ C : ℝ, C = Real.exp (L * (3 / 2)) * (8 + Bc B * (64 * (T₂ + 3) + 24)) *
      (1 + ‖riemannZeta s‖) := ⟨_, rfl⟩
  rw [← hCdef]
  have hC : 0 ≤ C := by rw [hCdef]; positivity
  have hgen : ∀ w : ℂ, w.re ≤ 2 → ‖(riemannZeta w)⁻¹‖ ≤ Bc B → (‖s - w‖ < 1 / 8 → ‖s‖ ≤ T₂ + 3) →
      ‖kint H L s w‖ ≤ C * ‖H w‖ := fun w hw2 hBw hn => by
    rw [hCdef]; exact norm_kint_le h215 hL hs hw2 hBc0 hR0 hBw hn
  have hvert : ∀ (σ y : ℝ), (σ = 2 ∨ σ = -1) → ‖s - ((σ : ℂ) + y * I)‖ < 1 / 8 → ‖s‖ ≤ T₂ + 3 := by
    intro σ y hσ h
    exfalso
    have h1 := Complex.abs_re_le_norm (s - ((σ : ℂ) + y * I))
    rw [sub_re, hs] at h1
    simp only [add_re, ofReal_re, mul_re, ofReal_im, I_re, I_im, mul_zero, mul_one, sub_zero] at h1
    rcases hσ with rfl | rfl <;> norm_num at h1 h <;> linarith
  have hhor : ∀ (x y : ℝ), x ∈ Icc (-1 : ℝ) 2 → |y| = T₂ → ‖s - ((x : ℂ) + y * I)‖ < 1 / 8 →
      ‖s‖ ≤ T₂ + 3 := by
    intro x y hx hy h
    have h1 : ‖(x : ℂ) + y * I‖ ≤ |x| + |y| := by
      have := Complex.norm_le_abs_re_add_abs_im ((x : ℂ) + y * I); simpa using this
    have h2 : |x| ≤ 2 := abs_le.mpr ⟨by linarith [hx.1], hx.2⟩
    have h3 := norm_le_norm_add_norm_sub' s ((x : ℂ) + y * I)
    linarith
  have hT0 : 0 < T₂ := by linarith
  have hU := norm_UpperU_le hH (f := kint H L s) hC (T := T₂)
    (fun x hx => hgen _ (by simp; exact hx.2) (le_Bc (hBu x hx)) (hhor x T₂ hx (abs_of_pos hT0)))
    (fun y _ => hgen _ (by norm_num) (inv_zeta_right_le hB0 y) (hvert 2 y (Or.inl rfl)))
    (fun y hy => hgen _ (by norm_num)
      (inv_zeta_left_le hB0 (by have : T₂ < y := hy; rw [abs_of_pos (by linarith)]; linarith))
      (hvert (-1) y (Or.inr rfl)))
  have hLo := norm_LowerU_le hH (f := kint H L s) hC (T := T₂)
    (fun x hx => hgen _ (by simp; exact hx.2) (le_Bc (hBl x hx))
      (hhor x (-T₂) hx (by rw [abs_neg, abs_of_pos hT0])))
    (fun y _ => hgen _ (by norm_num) (inv_zeta_right_le hB0 y) (hvert 2 y (Or.inl rfl)))
    (fun y hy => hgen _ (by norm_num)
      (inv_zeta_left_le hB0 (by have : y ≤ -T₂ := hy; rw [abs_of_neg (by linarith)]; linarith))
      (hvert (-1) y (Or.inr rfl)))
  unfold Kf
  rw [norm_mul, norm_two_pi_I_inv]
  have := norm_sub_le (UpperU (kint H L s) (-1) 2 T₂) (LowerU (kint H L s) (-1) 2 T₂)
  have hpi : 0 ≤ (2 * π)⁻¹ := by positivity
  calc (2 * π)⁻¹ * ‖UpperU (kint H L s) (-1) 2 T₂ - LowerU (kint H L s) (-1) 2 T₂‖
      ≤ (2 * π)⁻¹ * (C * ZU H T₂ + C * ZL H T₂) := by gcongr; linarith
    _ = (2 * π)⁻¹ * C * (ZU H T₂ + ZL H T₂) := by ring

theorem norm_inv_mul_kint_le {L : ℝ} (hL : 0 ≤ L) {σ : ℝ} (hσ : 1 ≤ σ) {w : ℂ} (hw2 : w.re ≤ 2)
    (him : 1 ≤ |w.im|) {Bw : ℝ} (hBw : ‖(riemannZeta w)⁻¹‖ ≤ Bw) :
    ‖(riemannZeta σ)⁻¹ * kint H L σ w‖ ≤ (Real.exp (L * (3 / 2)) * (Bw + ‖(riemannZeta σ)⁻¹‖)) * ‖H w‖ := by
  have hsw1 : 1 ≤ ‖(σ : ℂ) - w‖ := by
    have h := Complex.abs_im_le_norm ((σ : ℂ) - w)
    rw [sub_im, ofReal_im, zero_sub, abs_neg] at h
    linarith
  have hsw : (σ : ℂ) ≠ w := fun h => by rw [h, sub_self, norm_zero] at hsw1; norm_num at hsw1
  have hk : ‖kap L ((σ : ℂ) - w)‖ ≤ Real.exp (L * (3 / 2)) :=
    norm_kap_le_three_halves hL (by simp; linarith) hw2
  have hdsl : ‖dslope riemannZeta w σ‖ ≤ ‖riemannZeta σ‖ + ‖riemannZeta w‖ := by
    rw [dslope_of_ne _ hsw, slope_def_field, norm_div]
    exact (div_le_self (norm_nonneg _) hsw1).trans (norm_sub_le _ _)
  have h1 := norm_mul_norm_inv_le (riemannZeta σ)
  have h2 := norm_mul_norm_inv_le (riemannZeta w)
  have hA := norm_nonneg (riemannZeta (σ : ℂ))⁻¹
  have hWi := norm_nonneg (riemannZeta w)⁻¹
  have hHw := norm_nonneg (H w)
  have hE := Real.exp_pos (L * (3 / 2))
  rw [kint, div_eq_mul_inv, norm_mul, norm_mul, norm_mul, norm_mul]
  calc ‖(riemannZeta σ)⁻¹‖ * (‖H w‖ * ‖kap L (σ - w)‖ * ‖dslope riemannZeta w σ‖ * ‖(riemannZeta w)⁻¹‖)
      ≤ ‖(riemannZeta σ)⁻¹‖ * (‖H w‖ * Real.exp (L * (3 / 2)) * (‖riemannZeta σ‖ + ‖riemannZeta w‖) *
          ‖(riemannZeta w)⁻¹‖) := by gcongr
    _ = ‖H w‖ * Real.exp (L * (3 / 2)) * ((‖riemannZeta σ‖ * ‖(riemannZeta σ)⁻¹‖) * ‖(riemannZeta w)⁻¹‖ +
          ‖(riemannZeta σ)⁻¹‖ * (‖riemannZeta w‖ * ‖(riemannZeta w)⁻¹‖)) := by ring
    _ ≤ ‖H w‖ * Real.exp (L * (3 / 2)) * (1 * Bw + ‖(riemannZeta σ)⁻¹‖ * 1) := by
        gcongr
    _ = _ := by ring

/-- **K(σ)/ζ(σ)** for real σ ≥ 1, bounded by the contour mass of H. -/
theorem norm_Kf_div_zeta (hH : HBound H M β D N) {T₂ : ℝ} (hT1 : 1 ≤ T₂) {L : ℝ} (hL : 0 ≤ L) {B : ℝ}
    (hB0 : 0 ≤ B) (hBu : ∀ x ∈ Icc (-1 : ℝ) 2, ‖(riemannZeta (x + T₂ * I))⁻¹‖ ≤ B)
    (hBl : ∀ x ∈ Icc (-1 : ℝ) 2, ‖(riemannZeta (x + ((-T₂ : ℝ) : ℂ) * I))⁻¹‖ ≤ B) {σ : ℝ} (hσ : 1 ≤ σ) :
    ‖(riemannZeta σ)⁻¹ * Kf H L T₂ σ‖ ≤
      (2 * π)⁻¹ * (Real.exp (L * (3 / 2)) * (Bc B + ‖(riemannZeta σ)⁻¹‖)) * (ZU H T₂ + ZL H T₂) := by
  obtain ⟨C, hCdef⟩ : ∃ C : ℝ, C = Real.exp (L * (3 / 2)) * (Bc B + ‖(riemannZeta σ)⁻¹‖) := ⟨_, rfl⟩
  rw [← hCdef]
  have hBc0 := Bc_nonneg hB0
  have hC : 0 ≤ C := by rw [hCdef]; positivity
  have hT0 : 0 < T₂ := by linarith
  have hgen : ∀ w : ℂ, w.re ≤ 2 → 1 ≤ |w.im| → ‖(riemannZeta w)⁻¹‖ ≤ Bc B →
      ‖(riemannZeta σ)⁻¹ * kint H L σ w‖ ≤ C * ‖H w‖ := fun w hw2 him hBw => by
    rw [hCdef]; exact norm_inv_mul_kint_le hL hσ hw2 him hBw
  have hU := norm_UpperU_le hH (f := fun w => (riemannZeta σ)⁻¹ * kint H L σ w) hC (T := T₂)
    (fun x hx => hgen _ (by simp; exact hx.2) (by simp [abs_of_pos hT0]; exact hT1) (le_Bc (hBu x hx)))
    (fun y hy => hgen _ (by norm_num) (by have : T₂ < y := hy; simp [abs_of_pos (by linarith : (0:ℝ) < y)]; linarith)
      (inv_zeta_right_le hB0 y))
    (fun y hy => hgen _ (by norm_num) (by have : T₂ < y := hy; simp [abs_of_pos (by linarith : (0:ℝ) < y)]; linarith)
      (inv_zeta_left_le hB0 (by have : T₂ < y := hy; rw [abs_of_pos (by linarith)]; linarith)))
  have hLo := norm_LowerU_le hH (f := fun w => (riemannZeta σ)⁻¹ * kint H L σ w) hC (T := T₂)
    (fun x hx => hgen _ (by simp; exact hx.2) (by simp [abs_of_pos hT0]; exact hT1) (le_Bc (hBl x hx)))
    (fun y hy => hgen _ (by norm_num) (by have : y ≤ -T₂ := hy; simp [abs_of_neg (by linarith : y < 0)]; linarith)
      (inv_zeta_right_le hB0 y))
    (fun y hy => hgen _ (by norm_num) (by have : y ≤ -T₂ := hy; simp [abs_of_neg (by linarith : y < 0)]; linarith)
      (inv_zeta_left_le hB0 (by have : y ≤ -T₂ := hy; rw [abs_of_neg (by linarith)]; linarith)))
  rw [UpperU_const_mul] at hU
  rw [LowerU_const_mul] at hLo
  unfold Kf
  have e : (riemannZeta σ)⁻¹ * ((2 * (π : ℂ) * I)⁻¹ * (UpperU (kint H L σ) (-1) 2 T₂ - LowerU (kint H L σ) (-1) 2 T₂)) =
      (2 * (π : ℂ) * I)⁻¹ * ((riemannZeta σ)⁻¹ * UpperU (kint H L σ) (-1) 2 T₂ -
        (riemannZeta σ)⁻¹ * LowerU (kint H L σ) (-1) 2 T₂) := by ring
  rw [e, norm_mul, norm_two_pi_I_inv]
  have := norm_sub_le ((riemannZeta σ)⁻¹ * UpperU (kint H L σ) (-1) 2 T₂)
    ((riemannZeta σ)⁻¹ * LowerU (kint H L σ) (-1) 2 T₂)
  have hpi : 0 ≤ (2 * π)⁻¹ := by positivity
  calc (2 * π)⁻¹ * ‖(riemannZeta σ)⁻¹ * UpperU (kint H L σ) (-1) 2 T₂ -
        (riemannZeta σ)⁻¹ * LowerU (kint H L σ) (-1) 2 T₂‖
      ≤ (2 * π)⁻¹ * (C * ZU H T₂ + C * ZL H T₂) := by gcongr; linarith
    _ = (2 * π)⁻¹ * C * (ZU H T₂ + ZL H T₂) := by ring

end Riemann.Deep.LeanNBRate

