import Riemann.Deep.LeanNBRate.HConc3
import Riemann.Deep.LeanNBRate.GF1
import Riemann.Deep.LeanNBRate.LineE
import Riemann.Deep.LeanNBRate.LaplaceClosure

/-!
# The Nyman closure bound for the concrete construction (all terms on the critical line)

D_N(e^{-L})² ≤ (1+δ) N(p) + 3(1+1/δ)(2π)⁻¹ (8π/T₁ + 4πQ² + (2Q² + (R+G₁)²) Z_ζ).
-/

open MeasureTheory Filter Topology Set Complex ComplexConjugate
open scoped Real

noncomputable section

namespace Riemann.Deep.LeanNBRate

theorem integral_inv_sq_add {T : ℝ} (hT : 0 < T) :
    Integrable (fun t : ℝ => (T ^ 2 + t ^ 2)⁻¹) ∧ ∫ t : ℝ, (T ^ 2 + t ^ 2)⁻¹ = π / T := by
  have e : (fun t : ℝ => (T ^ 2 + t ^ 2)⁻¹) = fun t => (T ^ 2)⁻¹ * (1 + (T⁻¹ * t) ^ 2)⁻¹ := by
    funext t
    rw [← mul_inv]
    congr 1
    field_simp
  rw [e]
  refine ⟨(integrable_inv_one_add_mul_sq (inv_ne_zero hT.ne')).const_mul _, ?_⟩
  rw [MeasureTheory.integral_const_mul, integral_univ_inv_one_add_mul_sq, abs_inv, abs_of_pos hT]
  field_simp

/-- The K-line constant. -/
def QK (H : ℂ → ℂ) (L T₂ B : ℝ) : ℝ :=
  (2 * π)⁻¹ * (Real.exp (L * (3 / 2)) * (8 + Bc B * (64 * (T₂ + 3) + 24))) * (ZU H T₂ + ZL H T₂)

/-- The bound for G(1). -/
def G1b (H : ℂ → ℂ) (L T₂ B Rt : ℝ) : ℝ :=
  (2 * π)⁻¹ * (Real.exp (L * (3 / 2)) * Bc B) * (ZU H T₂ + ZL H T₂) + Rt

/-- The integral of |ζ(s)/s|² over the critical line. -/
def Zz : ℝ := ∫ t : ℝ, ‖riemannZeta (1 / 2 + t * I) / (1 / 2 + t * I)‖ ^ 2

theorem ZU_nonneg (H : ℂ → ℂ) (T : ℝ) : 0 ≤ ZU H T := by
  unfold ZU
  have h1 : 0 ≤ ∫ x in (-1 : ℝ)..2, ‖H (x + T * I)‖ :=
    intervalIntegral.integral_nonneg (by norm_num) fun _ _ => norm_nonneg _
  have h2 : 0 ≤ ∫ y in Ioi T, ‖H (((2 : ℝ) : ℂ) + y * I)‖ :=
    setIntegral_nonneg measurableSet_Ioi fun _ _ => norm_nonneg _
  have h3 : 0 ≤ ∫ y in Ioi T, ‖H (((-1 : ℝ) : ℂ) + y * I)‖ :=
    setIntegral_nonneg measurableSet_Ioi fun _ _ => norm_nonneg _
  linarith

theorem ZL_nonneg (H : ℂ → ℂ) (T : ℝ) : 0 ≤ ZL H T := by
  unfold ZL
  have h1 : 0 ≤ ∫ x in (-1 : ℝ)..2, ‖H (x + ((-T : ℝ) : ℂ) * I)‖ :=
    intervalIntegral.integral_nonneg (by norm_num) fun _ _ => norm_nonneg _
  have h2 : 0 ≤ ∫ y in Iic (-T), ‖H (((2 : ℝ) : ℂ) + y * I)‖ :=
    setIntegral_nonneg measurableSet_Iic fun _ _ => norm_nonneg _
  have h3 : 0 ≤ ∫ y in Iic (-T), ‖H (((-1 : ℝ) : ℂ) + y * I)‖ :=
    setIntegral_nonneg measurableSet_Iic fun _ _ => norm_nonneg _
  linarith

theorem sq_add_le_eps {x y δ : ℝ} (hδ : 0 < δ) : (x + y) ^ 2 ≤ (1 + δ) * x ^ 2 + (1 + δ⁻¹) * y ^ 2 := by
  have key : δ * (δ * x ^ 2 + δ⁻¹ * y ^ 2 - 2 * x * y) = (δ * x - y) ^ 2 := by
    field_simp
    ring
  have h1 : 0 ≤ δ * x ^ 2 + δ⁻¹ * y ^ 2 - 2 * x * y := by
    have h2 : 0 ≤ (δ * x - y) ^ 2 := sq_nonneg _
    rw [← key] at h2
    exact (mul_nonneg_iff_of_pos_left hδ).mp h2
  nlinarith

theorem sq_add3_le (x y z : ℝ) : (x + y + z) ^ 2 ≤ 3 * (x ^ 2 + y ^ 2 + z ^ 2) := by
  nlinarith [sq_nonneg (x - y), sq_nonneg (y - z), sq_nonneg (x - z)]

theorem norm_line_sq (t : ℝ) : ‖(1 / 2 : ℂ) + t * I‖ ^ 2 = (1 / 2) ^ 2 + t ^ 2 := by
  rw [Complex.sq_norm, Complex.normSq_apply]; simp; ring

set_option maxHeartbeats 2000000 in
/-- **The closure bound.** -/
theorem nyman_le (hP : Rudin1987_Thm9_13) (h215 : Titchmarsh_2_1_5) {ℓ a T₁ L T₂ M β D B Rt δ : ℝ}
    {n N : ℕ} (r : ExpSum) (hℓ : 0 ≤ ℓ) (ha : 0 < a) (hT₁ : 0 < T₁) (hna : (n : ℝ) * a = 2 / T₁)
    (hL : 0 < L) (hH : HBound (Hc ℓ a n r) M β D N) (hZ : HZeros (Hc ℓ a n r) T₂) (hT1 : 1 ≤ T₂)
    (hLβ : β + Real.log (2 * π) < L) (hB0 : 0 ≤ B)
    (hBu : ∀ x ∈ Icc (-1 : ℝ) 2, ‖(riemannZeta (x + T₂ * I))⁻¹‖ ≤ B)
    (hBl : ∀ x ∈ Icc (-1 : ℝ) 2, ‖(riemannZeta (x + ((-T₂ : ℝ) : ℂ) * I))⁻¹‖ ≤ B)
    (hRt : ∀ s : ℂ, 1 / 2 ≤ s.re → ‖Rtr (Hc ℓ a n r) L s‖ ≤ Rt) (hδ : 0 < δ) :
    nymanDist (Real.exp (-L)) ^ 2 ≤ (1 + δ) * Nself ℓ r.eval (r.dmod 0).eval +
      3 * (1 + δ⁻¹) * (2 * π)⁻¹ * (8 * π / T₁ + 4 * π * QK (Hc ℓ a n r) L T₂ B ^ 2 +
        (2 * QK (Hc ℓ a n r) L T₂ B ^ 2 + (Rt + G1b (Hc ℓ a n r) L T₂ B Rt) ^ 2) * Zz) := by
  obtain ⟨H, hHdef⟩ : ∃ H, H = Hc ℓ a n r := ⟨_, rfl⟩
  rw [← hHdef] at hH hZ hRt ⊢
  obtain ⟨Q, hQ⟩ : ∃ Q, Q = QK H L T₂ B := ⟨_, rfl⟩
  obtain ⟨G₁, hG₁⟩ : ∃ G₁, G₁ = G1b H L T₂ B Rt := ⟨_, rfl⟩
  rw [← hQ, ← hG₁]
  obtain ⟨Np, hNp⟩ : ∃ Np, Np = Nself ℓ r.eval (r.dmod 0).eval := ⟨_, rfl⟩
  rw [← hNp]
  have hLE := line_Elow hP hℓ r.eval (r.dmod 0).eval ((r.dmod 0).dmod 0).eval r.hasDerivAt_eval
    (r.dmod 0).hasDerivAt_eval r.continuous_eval (r.dmod 0).continuous_eval ((r.dmod 0).dmod 0).continuous_eval
  rw [← hNp] at hLE
  have hZm := memLp_zeta_line hP h215
  have hEi : Integrable fun t : ℝ => ‖Elow ℓ r.eval (r.dmod 0).eval (1 / 2 + t * I) / (1 / 2 + t * I)‖ ^ 2 :=
    (memLp_two_iff_integrable_sq_norm hLE.1.aestronglyMeasurable).mp hLE.1
  have hZi : Integrable fun t : ℝ => ‖riemannZeta (1 / 2 + t * I) / (1 / 2 + t * I)‖ ^ 2 :=
    (memLp_two_iff_integrable_sq_norm hZm.aestronglyMeasurable).mp hZm
  obtain ⟨hI1, hV1⟩ := integral_inv_sq_add hT₁
  obtain ⟨hI2, hV2⟩ := integral_inv_sq_add (by norm_num : (0 : ℝ) < 1 / 2)
  have hG1 : ‖GF H L 1‖ ≤ G₁ := by
    rw [hG₁]
    exact norm_GF_one_le hH hZ hT1 hL.le hLβ hB0 hBu hBl (fun σ hσ => hRt σ (by simp; linarith))
  have hGn : 0 ≤ G₁ := (norm_nonneg _).trans hG1
  have hRt0 : 0 ≤ Rt := (norm_nonneg _).trans (hRt (1 / 2 : ℂ) (by norm_num))
  have hQ0 : 0 ≤ Q := by
    rw [hQ]; unfold QK
    have h1 := Bc_nonneg hB0
    have h2 := ZU_nonneg H T₂
    have h3 := ZL_nonneg H T₂
    have h4 : 0 ≤ Real.exp (L * (3 / 2)) * (8 + Bc B * (64 * (T₂ + 3) + 24)) := by
      have : 0 ≤ Bc B * (64 * (T₂ + 3) + 24) := mul_nonneg h1 (by linarith)
      have := Real.exp_pos (L * (3 / 2))
      positivity
    have h5 : (0 : ℝ) ≤ (2 * π)⁻¹ := by positivity
    exact mul_nonneg (mul_nonneg h5 h4) (by linarith)
  set Fb : ℝ → ℝ := fun t => (1 + δ) *
      ‖Elow ℓ r.eval (r.dmod 0).eval (1 / 2 + t * I) / (1 / 2 + t * I)‖ ^ 2 +
      3 * (1 + δ⁻¹) * (8 * (T₁ ^ 2 + t ^ 2)⁻¹ + 2 * Q ^ 2 * ((1 / 2) ^ 2 + t ^ 2)⁻¹ +
        (2 * Q ^ 2 + (Rt + G₁) ^ 2) * ‖riemannZeta (1 / 2 + t * I) / (1 / 2 + t * I)‖ ^ 2) with hFb
  have hFi : Integrable Fb := by
    rw [hFb]
    exact (hEi.const_mul _).add ((((hI1.const_mul 8).add (hI2.const_mul _)).add (hZi.const_mul _)).const_mul _)
  have hpt : ∀ t : ℝ, ‖(1 + riemannZeta (1 / 2 + t * I) * (GF H L 1 + ∫ u in (0 : ℝ)..L,
      -gF H u * cexp (-(u * (1 / 2 + t * I))))) / (1 / 2 + t * I)‖ ^ 2 ≤ Fb t := by
    intro t
    obtain ⟨s, hsdef⟩ : ∃ s : ℂ, s = 1 / 2 + t * I := ⟨_, rfl⟩
    rw [← hsdef]
    have hs0 : s ≠ 0 := by rw [hsdef]; exact line_ne_zero t
    have hsre : s.re = 1 / 2 := by rw [hsdef]; simp
    have hsn : ‖s‖ ^ 2 = (1 / 2) ^ 2 + t ^ 2 := by rw [hsdef]; exact norm_line_sq t
    have hW : 0 < (1 / 2 : ℝ) ^ 2 + t ^ 2 := by positivity
    have hGs : GF H L 1 + ∫ u in (0 : ℝ)..L, -gF H u * cexp (-(u * s)) = GF H L 1 - GF H L s := by
      have : ∫ u in (0 : ℝ)..L, -gF H u * cexp (-(u * s)) = -GF H L s := by
        unfold GF; rw [← intervalIntegral.integral_neg]
        refine intervalIntegral.integral_congr fun u _ => ?_
        simp only [neg_mul]
      rw [this]; ring
    rw [hGs]
    have hlf := line_formula hH hZ hT1 hL.le hLβ (s := s) (by rw [hsre]; norm_num) (by rw [hsre]; norm_num)
    obtain ⟨E, hE⟩ : ∃ E, E = Elow ℓ r.eval (r.dmod 0).eval s := ⟨_, rfl⟩
    obtain ⟨Φ, hΦ⟩ : ∃ Φ, Φ = Phi a n s := ⟨_, rfl⟩
    have hHs : H s = (1 - E) * Φ := by rw [hHdef, hE, hΦ]; rfl
    rw [hHs] at hlf
    obtain ⟨ζs, hζs⟩ : ∃ ζs, ζs = riemannZeta s := ⟨_, rfl⟩
    rw [← hζs] at hlf ⊢
    have hdec : (1 + ζs * (GF H L 1 - GF H L s)) / s =
        Φ * (E / s) + ((1 - Φ) / s + Kf H L T₂ s / s + ζs * (Rtr H L s + GF H L 1) / s) := by
      field_simp
      linear_combination (-1 : ℂ) * hlf
    rw [hdec]
    have hΦ1 : ‖Φ‖ ≤ 1 := by rw [hΦ, hsdef]; exact norm_Phi_line_le ha.le n t
    have h1Φ : ‖1 - Φ‖ ≤ |t| / T₁ := by
      rw [hΦ, hsdef]
      have := norm_one_sub_Phi_line_le ha.le n t
      calc ‖1 - Phi a n (1 / 2 + t * I)‖ ≤ n * (a * |t| / 2) := this
        _ = (n * a) * |t| / 2 := by ring
        _ = |t| / T₁ := by rw [hna]; field_simp
    have h1Φ2 : ‖1 - Φ‖ ≤ 2 := by
      have := norm_sub_le (1 : ℂ) Φ; rw [norm_one] at this; linarith
    have hK := norm_Kf_line h215 hH hT1 hL.le hB0 hBu hBl t
    rw [← hsdef, ← hζs] at hK
    have hK' : ‖Kf H L T₂ s‖ ≤ Q * (1 + ‖ζs‖) := by
      rw [hQ]; unfold QK
      calc ‖Kf H L T₂ s‖ ≤ _ := hK
        _ = _ := by ring
    have hR := hRt s (by rw [hsre])
    have hA : ‖Φ * (E / s)‖ ≤ ‖E / s‖ := by
      rw [norm_mul]; exact mul_le_of_le_one_left (norm_nonneg _) hΦ1
    have htri : ‖Φ * (E / s) + ((1 - Φ) / s + Kf H L T₂ s / s + ζs * (Rtr H L s + GF H L 1) / s)‖ ≤
        ‖E / s‖ + (‖(1 - Φ) / s‖ + ‖Kf H L T₂ s / s‖ + ‖ζs * (Rtr H L s + GF H L 1) / s‖) := by
      have h1 := norm_add_le (Φ * (E / s)) ((1 - Φ) / s + Kf H L T₂ s / s + ζs * (Rtr H L s + GF H L 1) / s)
      have h2 : ‖(1 - Φ) / s + Kf H L T₂ s / s + ζs * (Rtr H L s + GF H L 1) / s‖ ≤
          ‖(1 - Φ) / s‖ + ‖Kf H L T₂ s / s‖ + ‖ζs * (Rtr H L s + GF H L 1) / s‖ := norm_add₃_le
      linarith
    have hsq := sq_add_le_eps (x := ‖E / s‖)
      (y := ‖(1 - Φ) / s‖ + ‖Kf H L T₂ s / s‖ + ‖ζs * (Rtr H L s + GF H L 1) / s‖) hδ
    have h3 := sq_add3_le ‖(1 - Φ) / s‖ ‖Kf H L T₂ s / s‖ ‖ζs * (Rtr H L s + GF H L 1) / s‖
    have hb1 : ‖(1 - Φ) / s‖ ^ 2 ≤ 8 * (T₁ ^ 2 + t ^ 2)⁻¹ := by
      rw [norm_div, div_pow, hsn]
      have hX1 : ‖1 - Φ‖ ^ 2 * T₁ ^ 2 ≤ t ^ 2 := by
        have h4 : ‖1 - Φ‖ * T₁ ≤ |t| := by
          have := (le_div_iff₀ hT₁).mp h1Φ; linarith
        have h5 := pow_le_pow_left₀ (by positivity) h4 2
        rw [mul_pow, sq_abs] at h5
        exact h5
      have hX2 : ‖1 - Φ‖ ^ 2 ≤ 4 := by
        have := pow_le_pow_left₀ (norm_nonneg _) h1Φ2 2; norm_num at this; linarith
      rw [← div_eq_mul_inv, div_le_div_iff₀ hW (by positivity)]
      have h6 : ‖1 - Φ‖ ^ 2 * t ^ 2 ≤ 4 * t ^ 2 := mul_le_mul_of_nonneg_right hX2 (sq_nonneg t)
      calc ‖1 - Φ‖ ^ 2 * (T₁ ^ 2 + t ^ 2) = ‖1 - Φ‖ ^ 2 * T₁ ^ 2 + ‖1 - Φ‖ ^ 2 * t ^ 2 := by ring
        _ ≤ t ^ 2 + 4 * t ^ 2 := add_le_add hX1 h6
        _ ≤ 8 * ((1 / 2) ^ 2 + t ^ 2) := by nlinarith [sq_nonneg t]
    have hb2 : ‖Kf H L T₂ s / s‖ ^ 2 ≤ 2 * Q ^ 2 * ((1 / 2) ^ 2 + t ^ 2)⁻¹ + 2 * Q ^ 2 * ‖ζs / s‖ ^ 2 := by
      have hKs : ‖Kf H L T₂ s‖ ^ 2 ≤ 2 * Q ^ 2 * (1 + ‖ζs‖ ^ 2) := by
        have h4 := pow_le_pow_left₀ (norm_nonneg _) hK' 2
        have h5 : (Q * (1 + ‖ζs‖)) ^ 2 ≤ 2 * Q ^ 2 * (1 + ‖ζs‖ ^ 2) := by
          have : (1 + ‖ζs‖) ^ 2 ≤ 2 * (1 + ‖ζs‖ ^ 2) := by nlinarith [sq_nonneg (1 - ‖ζs‖)]
          rw [mul_pow]
          have hq2 : 0 ≤ Q ^ 2 := sq_nonneg Q
          exact (mul_le_mul_of_nonneg_left this hq2).trans (le_of_eq (by ring))
        linarith
      have hz : ‖ζs / s‖ ^ 2 = ‖ζs‖ ^ 2 / ((1 / 2) ^ 2 + t ^ 2) := by
        rw [norm_div, div_pow ‖ζs‖ ‖s‖ 2, hsn]
      have hk : ‖Kf H L T₂ s / s‖ ^ 2 = ‖Kf H L T₂ s‖ ^ 2 / ((1 / 2) ^ 2 + t ^ 2) := by
        rw [norm_div, div_pow ‖Kf H L T₂ s‖ ‖s‖ 2, hsn]
      have e : 2 * Q ^ 2 * ((1 / 2) ^ 2 + t ^ 2)⁻¹ + 2 * Q ^ 2 * (‖ζs‖ ^ 2 / ((1 / 2) ^ 2 + t ^ 2)) =
          2 * Q ^ 2 * (1 + ‖ζs‖ ^ 2) / ((1 / 2) ^ 2 + t ^ 2) := by
        field_simp
      rw [hz, hk, e]
      exact div_le_div_of_nonneg_right hKs hW.le
    have hb3 : ‖ζs * (Rtr H L s + GF H L 1) / s‖ ^ 2 ≤ (Rt + G₁) ^ 2 * ‖ζs / s‖ ^ 2 := by
      have e : ζs * (Rtr H L s + GF H L 1) / s = (ζs / s) * (Rtr H L s + GF H L 1) := by ring
      rw [e, norm_mul, mul_pow, mul_comm]
      have h4 : ‖Rtr H L s + GF H L 1‖ ≤ Rt + G₁ := (norm_add_le _ _).trans (add_le_add hR hG1)
      have h5 := pow_le_pow_left₀ (norm_nonneg _) h4 2
      exact mul_le_mul_of_nonneg_right h5 (by positivity)
    have hsum : ‖(1 - Φ) / s‖ ^ 2 + ‖Kf H L T₂ s / s‖ ^ 2 + ‖ζs * (Rtr H L s + GF H L 1) / s‖ ^ 2 ≤
        8 * (T₁ ^ 2 + t ^ 2)⁻¹ + 2 * Q ^ 2 * ((1 / 2) ^ 2 + t ^ 2)⁻¹ +
          (2 * Q ^ 2 + (Rt + G₁) ^ 2) * ‖ζs / s‖ ^ 2 := by linarith
    have hc3 : 0 ≤ 3 * (1 + δ⁻¹) := by positivity
    have hFt : Fb t = (1 + δ) * ‖E / s‖ ^ 2 + 3 * (1 + δ⁻¹) * (8 * (T₁ ^ 2 + t ^ 2)⁻¹ +
        2 * Q ^ 2 * ((1 / 2) ^ 2 + t ^ 2)⁻¹ + (2 * Q ^ 2 + (Rt + G₁) ^ 2) * ‖ζs / s‖ ^ 2) := by
      rw [hFb, hE, hζs, hsdef]
    rw [hFt]
    have h6 := pow_le_pow_left₀ (norm_nonneg _) htri 2
    have h7 := mul_le_mul_of_nonneg_left hsum hc3
    have h8 : (1 + δ⁻¹) * (‖(1 - Φ) / s‖ + ‖Kf H L T₂ s / s‖ + ‖ζs * (Rtr H L s + GF H L 1) / s‖) ^ 2 ≤
        (1 + δ⁻¹) * (3 * (‖(1 - Φ) / s‖ ^ 2 + ‖Kf H L T₂ s / s‖ ^ 2 +
          ‖ζs * (Rtr H L s + GF H L 1) / s‖ ^ 2)) := mul_le_mul_of_nonneg_left h3 (by positivity)
    nlinarith
  have hgc : Continuous fun u : ℝ => -gF H u := (continuous_gF hH).neg
  have hc0 : GF H L 1 + ∫ u in (0 : ℝ)..L, (fun u : ℝ => -gF H u) u * cexp (-(u : ℂ)) = 0 := by
    have : ∫ u in (0 : ℝ)..L, (fun u : ℝ => -gF H u) u * cexp (-(u : ℂ)) = -GF H L 1 := by
      unfold GF
      rw [← intervalIntegral.integral_neg]
      refine intervalIntegral.integral_congr fun u _ => ?_
      simp only [mul_one, neg_mul]
    rw [this]; ring
  have hcl := nymanDist_sq_le_laplace hP h215 hL hgc (GF H L 1) hc0
  refine hcl.trans ?_
  refine (mul_le_mul_of_nonneg_left (integral_mono_of_nonneg
    (Filter.Eventually.of_forall fun t => by positivity) hFi (Filter.Eventually.of_forall hpt))
    (by positivity)).trans ?_
  have hval : ∫ t, Fb t = (1 + δ) * (2 * π * Np) + 3 * (1 + δ⁻¹) * (8 * (π / T₁) + 2 * Q ^ 2 * (π / (1 / 2)) +
      (2 * Q ^ 2 + (Rt + G₁) ^ 2) * Zz) := by
    have i1 : Integrable (fun t : ℝ => (1 + δ) *
        ‖Elow ℓ r.eval (r.dmod 0).eval (1 / 2 + t * I) / (1 / 2 + t * I)‖ ^ 2) := hEi.const_mul _
    have j1 : Integrable (fun t : ℝ => 8 * (T₁ ^ 2 + t ^ 2)⁻¹ + 2 * Q ^ 2 * ((1 / 2) ^ 2 + t ^ 2)⁻¹) :=
      (hI1.const_mul 8).add (hI2.const_mul _)
    have j2 : Integrable (fun t : ℝ => (2 * Q ^ 2 + (Rt + G₁) ^ 2) *
        ‖riemannZeta (1 / 2 + t * I) / (1 / 2 + t * I)‖ ^ 2) := hZi.const_mul _
    have i2 : Integrable (fun t : ℝ => 3 * (1 + δ⁻¹) * (8 * (T₁ ^ 2 + t ^ 2)⁻¹ +
        2 * Q ^ 2 * ((1 / 2) ^ 2 + t ^ 2)⁻¹ + (2 * Q ^ 2 + (Rt + G₁) ^ 2) *
          ‖riemannZeta (1 / 2 + t * I) / (1 / 2 + t * I)‖ ^ 2)) := (j1.add j2).const_mul _
    rw [hFb]
    rw [integral_add i1 i2, MeasureTheory.integral_const_mul, hLE.2, MeasureTheory.integral_const_mul,
      integral_add j1 j2, integral_add (hI1.const_mul 8) (hI2.const_mul _), MeasureTheory.integral_const_mul,
      MeasureTheory.integral_const_mul, MeasureTheory.integral_const_mul, hV1, hV2]
    rfl
  rw [hval]
  have hpi : 0 < π := Real.pi_pos
  apply le_of_eq
  field_simp
  ring

end Riemann.Deep.LeanNBRate

