import Riemann.Deep.FirstZero.Cert
import Riemann.Deep.FirstZero.GammaRatio
import Riemann.Deep.Hardy.Repr

/-!
# The first zero of ζ on the critical line, certified in Lean

**Theorem.** There is t ∈ (14.13, 14.14) with ζ(1/2 + it) = 0.

Proof. Λ(1/2 + it) = π^{-s/2}Γ(s/2)ζ(s) is real (Fayed.Hardy.Lam_conj) and continuous in t.
At t₁ = 14.1301 and t₂ = 14.1399 the order-8 Euler–Maclaurin enclosures (Cert.lean) give
ζ(s₂)·conj ζ(s₁) ≈ -1.50·10⁻⁵ with error ≤ 4.3·10⁻⁶, while Γℝ(s₂)/Γℝ(s₁) is within 0.3 of 1
(GammaRatio.lean). Hence Λ(s₁)Λ(s₂) = |Γℝ(s₁)|²·Re(ζ(s₂) conj ζ(s₁) Γℝ(s₂)/Γℝ(s₁)) < 0, and the
intermediate value theorem gives a zero of Λ, hence of ζ, in between.
-/

open Complex Set

namespace Fayed.FirstZero

lemma re_mul_conj_neg {z1 z2 ρ : ℂ} {a1 b1 r1 a2 b2 r2 R : ℝ}
    (h1 : ‖z1 - (a1 + b1 * I)‖ ≤ r1) (h2 : ‖z2 - (a2 + b2 * I)‖ ≤ r2) (hρ : ‖ρ - 1‖ ≤ 3 / 10)
    (hR : (|a2| + |b2|) * r1 + (|a1| + |-b1|) * r2 + r2 * r1 + |a2 * a1 - b2 * -b1 - (a2 * a1 + b2 * b1)|
      + |a2 * -b1 + b2 * a1 - (b2 * a1 - a2 * b1)| ≤ R)
    (hnum : (a2 * a1 + b2 * b1) + R + 3 / 10 * (|a2 * a1 + b2 * b1| + |b2 * a1 - a2 * b1| + R) < 0) :
    (z2 * (starRingEnd ℂ) z1 * ρ).re < 0 := by
  have hc : ‖(starRingEnd ℂ) z1 - (a1 + ((-b1 : ℝ) : ℂ) * I)‖ ≤ r1 := by
    rw [← norm_conj]
    convert h1 using 2
    simp [map_sub, map_add, map_mul, Complex.conj_ofReal, Complex.conj_I]
  have hw := ball_mul (e := a2 * a1 + b2 * b1) (f := b2 * a1 - a2 * b1) (R := R) h2 hc hR
  set w := z2 * (starRingEnd ℂ) z1
  set w0 : ℂ := ((a2 * a1 + b2 * b1 : ℝ) : ℂ) + ((b2 * a1 - a2 * b1 : ℝ) : ℂ) * I
  have hre : w.re ≤ (a2 * a1 + b2 * b1) + R := by
    have := (re_le_norm (w - w0)).trans hw
    simp [w0] at this
    linarith
  have hnorm : ‖w‖ ≤ |a2 * a1 + b2 * b1| + |b2 * a1 - a2 * b1| + R := by
    have h0 : ‖w0‖ ≤ |a2 * a1 + b2 * b1| + |b2 * a1 - a2 * b1| :=
      (norm_le_abs_re_add_abs_im _).trans (le_of_eq (by simp [w0]))
    have := norm_le_norm_add_norm_sub' w w0
    linarith
  have e : (w * ρ).re = w.re + (w * (ρ - 1)).re := by
    rw [mul_sub, mul_one, sub_re]; ring
  have h3 : (w * (ρ - 1)).re ≤ ‖w‖ * (3 / 10) := by
    refine (re_le_norm _).trans ?_
    rw [norm_mul]
    exact mul_le_mul_of_nonneg_left hρ (norm_nonneg _)
  rw [e]
  nlinarith [norm_nonneg w]

/-- Γℝ(1/2 + i t₂) / Γℝ(1/2 + i t₁) is within 3/10 of 1. -/
lemma Gammaℝ_ratio :
    ‖Gammaℝ ((1 / 2 : ℂ) + ((141399 / 10000 : ℝ) : ℂ) * I) /
      Gammaℝ ((1 / 2 : ℂ) + ((141301 / 10000 : ℝ) : ℂ) * I) - 1‖ ≤ 3 / 10 := by
  set z : ℂ := ((1 / 2 : ℂ) + ((141301 / 10000 : ℝ) : ℂ) * I) / 2 with hz
  set δ : ℝ := 49 / 10000 with hδ
  have hs2 : ((1 / 2 : ℂ) + ((141399 / 10000 : ℝ) : ℂ) * I) / 2 = z + δ * I := by
    rw [hz, hδ]; push_cast; ring
  have hzre : z.re = 1 / 4 := by rw [hz]; simp; norm_num
  have hzim : z.im = 141301 / 20000 := by rw [hz]; simp; norm_num
  have hG := Gamma_ratio_le z hzre (by rw [hzim]; norm_num) (by rw [hzim]; norm_num)
    (δ := δ) (by norm_num [hδ]) (by norm_num [hδ])
  have hpi0 : ((Real.pi : ℝ) : ℂ) ≠ 0 := by exact_mod_cast Real.pi_ne_zero
  have hP : (Real.pi : ℂ) ^ (-((1 / 2 : ℂ) + ((141399 / 10000 : ℝ) : ℂ) * I) / 2) =
      (Real.pi : ℂ) ^ (-((1 / 2 : ℂ) + ((141301 / 10000 : ℝ) : ℂ) * I) / 2) *
        cexp (I * ((-(δ * Real.log Real.pi)) : ℝ)) := by
    have hPI : (Real.pi : ℂ) ^ (-((δ : ℂ) * I)) = cexp (I * ((-(δ * Real.log Real.pi)) : ℝ)) := by
      rw [cpow_def_of_ne_zero hpi0, ← ofReal_log Real.pi_pos.le]
      congr 1; push_cast; ring
    rw [show -((1 / 2 : ℂ) + ((141399 / 10000 : ℝ) : ℂ) * I) / 2 =
      -((1 / 2 : ℂ) + ((141301 / 10000 : ℝ) : ℂ) * I) / 2 + (-((δ : ℂ) * I)) by
        rw [hδ]; push_cast; ring,
      cpow_add _ _ hpi0, hPI]
  have hlogpi : Real.log Real.pi ≤ 22 / 10 := by
    have := Real.log_le_sub_one_of_pos Real.pi_pos
    have := Real.pi_lt_d2
    norm_num at *; linarith
  have hlogpi0 : 0 ≤ Real.log Real.pi := Real.log_nonneg (by linarith [Real.pi_gt_three])
  have hPb : ‖cexp (I * ((-(δ * Real.log Real.pi)) : ℝ)) - 1‖ ≤ 11 / 1000 := by
    refine (Real.norm_exp_I_mul_ofReal_sub_one_le).trans ?_
    rw [Real.norm_eq_abs, abs_neg, abs_of_nonneg (by rw [hδ]; positivity), hδ]
    nlinarith
  have hGne : Gamma z ≠ 0 := Gamma_ne_zero_of_re_pos (by rw [hzre]; norm_num)
  have hPne : (Real.pi : ℂ) ^ (-((1 / 2 : ℂ) + ((141301 / 10000 : ℝ) : ℂ) * I) / 2) ≠ 0 := by
    rw [Ne, cpow_eq_zero_iff]; tauto
  have hratio : Gammaℝ ((1 / 2 : ℂ) + ((141399 / 10000 : ℝ) : ℂ) * I) /
      Gammaℝ ((1 / 2 : ℂ) + ((141301 / 10000 : ℝ) : ℂ) * I) =
      cexp (I * ((-(δ * Real.log Real.pi)) : ℝ)) * (Gamma (z + δ * I) / Gamma z) := by
    simp only [Gammaℝ]
    rw [hP, hs2, ← hz]
    field_simp
  rw [hratio]
  set P := cexp (I * ((-(δ * Real.log Real.pi)) : ℝ))
  set G := Gamma (z + δ * I) / Gamma z
  have e : P * G - 1 = (P - 1) * G + (G - 1) := by ring
  have hGn : ‖G‖ ≤ 5 / 4 := by
    have := norm_le_norm_add_norm_sub' G 1
    simp at this; linarith
  rw [e]
  calc ‖(P - 1) * G + (G - 1)‖ ≤ ‖P - 1‖ * ‖G‖ + ‖G - 1‖ :=
        (norm_add_le _ _).trans (by rw [norm_mul])
    _ ≤ 11 / 1000 * (5 / 4) + 1 / 4 := by gcongr
    _ ≤ 3 / 10 := by norm_num

lemma completed_eq_mul {s : ℂ} (hs : 0 < s.re) :
    completedRiemannZeta s = riemannZeta s * Gammaℝ s := by
  have hs0 : s ≠ 0 := by rintro rfl; simp at hs
  rw [riemannZeta_def_of_ne_zero hs0, div_mul_cancel₀ _ (Gammaℝ_ne_zero_of_re_pos hs)]

lemma Lam_im (t : ℝ) : (completedRiemannZeta ((1 / 2 : ℂ) + (t : ℂ) * I)).im = 0 :=
  (conj_eq_iff_im.mp (Fayed.Hardy.Lam_conj t))

/-- **The first zero of ζ on the critical line lies in (14.13, 14.14).** -/
theorem first_zero : ∃ t : ℝ, 14.13 < t ∧ t < 14.14 ∧
    riemannZeta ((1 / 2 : ℂ) + t * Complex.I) = 0 := by
  set t1 : ℝ := 141301 / 10000 with ht1
  set t2 : ℝ := 141399 / 10000 with ht2
  set f : ℝ → ℝ := fun t => (completedRiemannZeta ((1 / 2 : ℂ) + (t : ℂ) * I)).re with hf
  have hre : ∀ t : ℝ, 0 < ((1 / 2 : ℂ) + (t : ℂ) * I).re := fun t => by simp
  -- the sign change
  have hsign : f t1 * f t2 < 0 := by
    have h1 := zeta_t1
    have h2 := zeta_t2
    have hρ := Gammaℝ_ratio
    have key := re_mul_conj_neg h1 h2 hρ (R := 43 / 10000000) (by norm_num) (by norm_num)
    set s1 : ℂ := (1 / 2 : ℂ) + ((141301 / 10000 : ℝ) : ℂ) * I
    set s2 : ℂ := (1 / 2 : ℂ) + ((141399 / 10000 : ℝ) : ℂ) * I
    have hG1 : Gammaℝ s1 ≠ 0 := Gammaℝ_ne_zero_of_re_pos (hre _)
    have hL1 : completedRiemannZeta s1 = riemannZeta s1 * Gammaℝ s1 := completed_eq_mul (hre _)
    have hL2 : completedRiemannZeta s2 = riemannZeta s2 * Gammaℝ s2 := completed_eq_mul (hre _)
    have hi1 := Lam_im (141301 / 10000)
    have hi2 := Lam_im (141399 / 10000)
    -- Λ₂ · conj Λ₁ = (ζ₂ conj ζ₁ ρ) · normSq Γℝ(s₁)
    have hprod : completedRiemannZeta s2 * (starRingEnd ℂ) (completedRiemannZeta s1) =
        (riemannZeta s2 * (starRingEnd ℂ) (riemannZeta s1) * (Gammaℝ s2 / Gammaℝ s1)) *
          ((normSq (Gammaℝ s1) : ℝ) : ℂ) := by
      rw [hL1, hL2, ← mul_conj, map_mul]
      field_simp
    have hns : 0 < normSq (Gammaℝ s1) := normSq_pos.mpr hG1
    have hre_prod : (completedRiemannZeta s2 * (starRingEnd ℂ) (completedRiemannZeta s1)).re =
        f t2 * f t1 := by
      simp only [hf, ht1, ht2, mul_re, conj_re, conj_im]
      rw [hi1, hi2]; ring
    rw [hprod, re_mul_ofReal] at hre_prod
    have : f t2 * f t1 < 0 := by rw [← hre_prod]; exact mul_neg_of_neg_of_pos key hns
    linarith [mul_comm (f t1) (f t2)]
  -- continuity
  have hcont : ContinuousOn f (Icc t1 t2) := by
    intro t _
    have hs0 : ((1 / 2 : ℂ) + (t : ℂ) * I) ≠ 0 := by
      intro h; have := congrArg Complex.re h; simp at this
    have hs1 : ((1 / 2 : ℂ) + (t : ℂ) * I) ≠ 1 := by
      intro h; have := congrArg Complex.re h; norm_num at this
    have hc1 : ContinuousAt (fun t : ℝ => (1 / 2 : ℂ) + (t : ℂ) * I) t := by fun_prop
    have hc2 := (differentiableAt_completedZeta hs0 hs1).continuousAt
    have hc3 := ContinuousAt.comp (g := completedRiemannZeta)
      (f := fun t : ℝ => (1 / 2 : ℂ) + (t : ℂ) * I) (x := t) hc2 hc1
    exact (continuous_re.continuousAt.comp hc3).continuousWithinAt
  have hle : t1 ≤ t2 := by norm_num [ht1, ht2]
  -- intermediate value theorem
  obtain ⟨t, ht, hft⟩ : ∃ t ∈ Ioo t1 t2, f t = 0 := by
    rcases lt_or_gt_of_ne (show f t1 ≠ 0 by intro h; rw [h, zero_mul] at hsign; exact lt_irrefl _ hsign)
      with hneg | hpos
    · have h2 : 0 < f t2 := by
        by_contra h; push Not at h; nlinarith
      exact intermediate_value_Ioo hle hcont ⟨hneg, h2⟩
    · have h2 : f t2 < 0 := by
        by_contra h; push Not at h; nlinarith
      exact intermediate_value_Ioo' hle hcont ⟨h2, hpos⟩
  refine ⟨t, ?_, ?_, ?_⟩
  · have := ht.1; rw [ht1] at this; norm_num at this ⊢; linarith
  · have := ht.2; rw [ht2] at this; norm_num at this ⊢; linarith
  · have hΛ : completedRiemannZeta ((1 / 2 : ℂ) + (t : ℂ) * I) = 0 :=
      Complex.ext hft (by simpa using Lam_im t)
    have hs0 : ((1 / 2 : ℂ) + (t : ℂ) * I) ≠ 0 := by
      intro h; have := congrArg Complex.re h; simp at this
    rw [riemannZeta_def_of_ne_zero hs0, hΛ, zero_div]

end Fayed.FirstZero
