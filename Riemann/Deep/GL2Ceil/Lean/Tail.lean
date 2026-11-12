import Riemann.Deep.GL2Ceil.Lean.Cells
import Riemann.Deep.GL2Ceil.Lean.Closed
import Riemann.Deep.GL2Ceil.Lean.Gd
import Mathlib.Analysis.Complex.ExponentialBounds

/-!
# GL2Ceil (gl2lean): the small-R tails proved (DERIVATION §4.1-4.2, derivative-free variant)

(A) For 0 < t ≤ 1, γ(1-t)R²/t ≤ tα: J_R(Q) ≥ (1-t) J_0(Q) (e^{2Ry} ≥ 1 and (a+b)² ≥ (1-t)a² - ((1-t)/t)b²;
    replaces the slope bound of §4.2).
(B) §4.1: if 4αγ ≥ 1 then J_0(Q) ≥ 2 - 1/(3/4 + α/4) (Q = β/2 + h, h odd about 1/2,
    αh² + γh'² ≥ ±hh', FTC on [0,1/2] and [1/2,1]); 4 C_λ B_λ ≥ 1 from the inner minimum (ClBl_ge).
(C) Numbers: R ≤ 1/2, all λ ≥ 0 (t = 0.0289, e < 2.7182818286) and R ∈ [1/2, 19/32], λ ≥ 3
    (t = 1/45, C_3 ≥ 3/2, B_3 ≤ 1/5). So TailSmallR and the first half of TailLargeLam are theorems.
-/

open Real MeasureTheory Set

noncomputable section

namespace GL2Ceil

lemma sq_shift (a b t : ℝ) (ht : 0 < t) : (1 - t) * a ^ 2 - (1 - t) / t * b ^ 2 ≤ (a + b) ^ 2 := by
  have h0 : 0 ≤ (t * a + b) ^ 2 / t := div_nonneg (sq_nonneg _) ht.le
  have e : (t * a + b) ^ 2 / t = (a + b) ^ 2 - ((1 - t) * a ^ 2 - (1 - t) / t * b ^ 2) := by
    field_simp
    ring
  linarith

/-- (A): J_R ≥ (1 - t) J_0. -/
theorem J_R_ge {R α γ t : ℝ} (hR : 0 ≤ R) (hα : 0 ≤ α) (hγ : 0 ≤ γ) (ht : 0 < t)
    (hcond : γ * ((1 - t) / t) * R ^ 2 ≤ t * α) {Q : ℝ → ℝ} (hQ : QAdm Q) :
    (1 - t) * Jag 0 α γ Q ≤ Jag R α γ Q := by
  obtain ⟨hQd, hQc⟩ := contDiff_one_iff_deriv.mp hQ.smooth
  have hQ0 := hQd.continuous
  have hw0 : ∀ y, w 0 Q y = Q y := fun y => by simp [w]
  have hwd0 : ∀ y, wd 0 Q y = deriv Q y := fun y => by simp [wd]
  have hE : ∀ y ∈ Icc (0:ℝ) 1, 1 ≤ Real.exp (R * y) :=
    fun y hy => Real.one_le_exp (mul_nonneg hR hy.1)
  have hW0 : W0 0 Q ≤ W0 R Q := by
    unfold W0
    apply intervalIntegral.integral_mono_on zero_le_one
      (((w_cont hQ0).fun_pow 2).intervalIntegrable _ _) (((w_cont hQ0).fun_pow 2).intervalIntegrable _ _)
    intro y hy
    rw [hw0]
    simp only [w]
    have h1 := hE y hy
    nlinarith [mul_nonneg (mul_nonneg (sub_nonneg.2 h1) (by linarith : (0:ℝ) ≤ Real.exp (R * y) + 1))
      (sq_nonneg (Q y))]
  have hW1 : (1 - t) * W1 0 Q - (1 - t) / t * R ^ 2 * W0 0 Q ≤ W1 R Q := by
    have hlhs : (1 - t) * W1 0 Q - (1 - t) / t * R ^ 2 * W0 0 Q
        = ∫ y in (0:ℝ)..1, ((1 - t) * deriv Q y ^ 2 - (1 - t) / t * R ^ 2 * Q y ^ 2) := by
      unfold W0 W1
      simp_rw [hw0, hwd0]
      rw [intervalIntegral.integral_sub (((hQc.fun_pow 2).intervalIntegrable 0 1).const_mul _)
        (((hQ0.fun_pow 2).intervalIntegrable 0 1).const_mul _), intervalIntegral.integral_const_mul,
        intervalIntegral.integral_const_mul]
    rw [hlhs]
    unfold W1
    apply intervalIntegral.integral_mono_on zero_le_one
      ((((hQc.fun_pow 2).intervalIntegrable 0 1).const_mul _).sub
        (((hQ0.fun_pow 2).intervalIntegrable 0 1).const_mul _))
      (((wd_cont hQ0 hQc).fun_pow 2).intervalIntegrable _ _)
    intro y hy
    have h1 := hE y hy
    have hs := sq_shift (deriv Q y) (R * Q y) t ht
    simp only [wd]
    have h2 : (deriv Q y + R * Q y) ^ 2 ≤ (Real.exp (R * y) * (deriv Q y + R * Q y)) ^ 2 := by
      nlinarith [mul_nonneg (mul_nonneg (sub_nonneg.2 h1)
        (by linarith : (0:ℝ) ≤ Real.exp (R * y) + 1)) (sq_nonneg (deriv Q y + R * Q y))]
    nlinarith [hs, h2]
  have hw1 : w 0 Q 1 ^ 2 ≤ w R Q 1 ^ 2 := by
    rw [hw0]
    simp only [w]
    have h1 := hE 1 ⟨zero_le_one, le_refl 1⟩
    nlinarith [mul_nonneg (mul_nonneg (sub_nonneg.2 h1)
      (by linarith : (0:ℝ) ≤ Real.exp (R * 1) + 1)) (sq_nonneg (Q 1))]
  have hW00 := W0_nonneg 0 Q
  have hW10 := W1_nonneg 0 Q
  have hprod : 0 ≤ (t * α - γ * ((1 - t) / t) * R ^ 2) * W0 0 Q :=
    mul_nonneg (sub_nonneg.2 hcond) hW00
  have hγW : γ * ((1 - t) * W1 0 Q - (1 - t) / t * R ^ 2 * W0 0 Q) ≤ γ * W1 R Q :=
    mul_le_mul_of_nonneg_left hW1 hγ
  have hαW : α * W0 0 Q ≤ α * W0 R Q := mul_le_mul_of_nonneg_left hW0 hα
  have htw : 0 ≤ t * (1 + w 0 Q 1 ^ 2) := mul_nonneg ht.le (by positivity)
  unfold Jag
  nlinarith [hw1, hprod, hγW, hαW, htw]

/-- (B), §4.1: J_0(Q) ≥ 2 - 1/(3/4 + α/4) when 4αγ ≥ 1. -/
theorem J0_ge {α γ : ℝ} (hα : 0 < α) (hαγ : 1 ≤ 4 * α * γ) {Q : ℝ → ℝ}
    (hQ : QAdm Q) : 2 - 1 / (3 / 4 + α / 4) ≤ Jag 0 α γ Q := by
  obtain ⟨hQd, hQc⟩ := contDiff_one_iff_deriv.mp hQ.smooth
  have hQ0 := hQd.continuous
  obtain ⟨β, hβ⟩ := hQ.symm
  have hQ1 : Q 1 = β - 1 := by
    have := hβ 0 ⟨le_refl 0, zero_le_one⟩
    rw [sub_zero, hQ.zero] at this
    linarith
  have hQh : Q (1 / 2) = β / 2 := by
    have := hβ (1 / 2) ⟨by norm_num, by norm_num⟩
    rw [show (1:ℝ) - 1 / 2 = 1 / 2 by norm_num] at this
    linarith
  have hw0 : ∀ y, w 0 Q y = Q y := fun y => by simp [w]
  have hwd0 : ∀ y, wd 0 Q y = deriv Q y := fun y => by simp [wd]
  have hch : Continuous fun y => Q y - β / 2 := by fun_prop
  have hB1 : ∫ y in (0:ℝ)..1, (Q y - β / 2) = 0 := by
    have h := intervalIntegral.integral_comp_sub_left (fun y => Q y - β / 2) (a := 0) (b := 1) 1
    simp only [sub_self, sub_zero] at h
    have h2 : ∫ x in (0:ℝ)..1, (Q (1 - x) - β / 2) = ∫ x in (0:ℝ)..1, -(Q x - β / 2) := by
      apply intervalIntegral.integral_congr
      intro x hx
      rw [uIcc_of_le zero_le_one] at hx
      have := hβ x hx
      simp only
      linarith
    rw [h2, intervalIntegral.integral_neg] at h
    linarith
  have hB2 : W0 0 Q = β ^ 2 / 4 + ∫ y in (0:ℝ)..1, (Q y - β / 2) ^ 2 := by
    unfold W0
    simp_rw [hw0]
    have h := integral_lin3 (β ^ 2 / 4) β 1 (continuous_const : Continuous fun _ : ℝ => (1:ℝ))
      hch (hch.fun_pow 2)
    rw [show (∫ y in (0:ℝ)..1, Q y ^ 2)
        = ∫ y in (0:ℝ)..1, (β ^ 2 / 4 * 1 + β * (Q y - β / 2) + 1 * (Q y - β / 2) ^ 2) from
      intervalIntegral.integral_congr (fun y _ => by ring), h, hB1]
    simp
  have hB3 : (1 - β / 2) ^ 2 ≤ α * (∫ y in (0:ℝ)..1, (Q y - β / 2) ^ 2)
      + γ * ∫ y in (0:ℝ)..1, deriv Q y ^ 2 := by
    have hf : Continuous fun y => α * (Q y - β / 2) ^ 2 + γ * deriv Q y ^ 2 := by fun_prop
    have hsum : α * (∫ y in (0:ℝ)..1, (Q y - β / 2) ^ 2) + γ * ∫ y in (0:ℝ)..1, deriv Q y ^ 2
        = ∫ y in (0:ℝ)..1, (α * (Q y - β / 2) ^ 2 + γ * deriv Q y ^ 2) := by
      rw [intervalIntegral.integral_add (((hch.fun_pow 2).intervalIntegrable 0 1).const_mul α)
        (((hQc.fun_pow 2).intervalIntegrable 0 1).const_mul γ), intervalIntegral.integral_const_mul,
        intervalIntegral.integral_const_mul]
    rw [hsum, ← intervalIntegral.integral_add_adjacent_intervals (b := 1 / 2)
      (hf.intervalIntegrable _ _) (hf.intervalIntegrable _ _)]
    have hprim : ∀ x, HasDerivAt (fun y => (Q y - β / 2) * (Q y - β / 2) / 2)
        ((Q x - β / 2) * deriv Q x) x := by
      intro x
      have h1 := (hQd x).hasDerivAt.sub_const (β / 2)
      have := (h1.fun_mul h1).div_const 2
      convert this using 1
      ring
    have hgi : IntervalIntegrable (fun y => (Q y - β / 2) * deriv Q y) volume 0 (1 / 2) :=
      (hch.fun_mul hQc).intervalIntegrable _ _
    have hgi' : IntervalIntegrable (fun y => (Q y - β / 2) * deriv Q y) volume (1 / 2) 1 :=
      (hch.fun_mul hQc).intervalIntegrable _ _
    have hint1 : ∫ y in (0:ℝ)..(1 / 2), (Q y - β / 2) * deriv Q y
        = (Q (1 / 2) - β / 2) * (Q (1 / 2) - β / 2) / 2 - (Q 0 - β / 2) * (Q 0 - β / 2) / 2 :=
      intervalIntegral.integral_eq_sub_of_hasDerivAt (fun x _ => hprim x) hgi
    have hint2 : ∫ y in (1 / 2 : ℝ)..1, (Q y - β / 2) * deriv Q y
        = (Q 1 - β / 2) * (Q 1 - β / 2) / 2 - (Q (1 / 2) - β / 2) * (Q (1 / 2) - β / 2) / 2 :=
      intervalIntegral.integral_eq_sub_of_hasDerivAt (fun x _ => hprim x) hgi'
    have hp1 : ∫ y in (0:ℝ)..(1 / 2), -((Q y - β / 2) * deriv Q y)
        ≤ ∫ y in (0:ℝ)..(1 / 2), (α * (Q y - β / 2) ^ 2 + γ * deriv Q y ^ 2) := by
      apply intervalIntegral.integral_mono_on (by norm_num) hgi.neg (hf.intervalIntegrable _ _)
      intro y _
      show -((Q y - β / 2) * deriv Q y) ≤ α * (Q y - β / 2) ^ 2 + γ * deriv Q y ^ 2
      nlinarith [sq_nonneg (α * (Q y - β / 2) + deriv Q y / 2),
        mul_nonneg (sub_nonneg.2 hαγ) (sq_nonneg (deriv Q y)), hα]
    have hp2 : ∫ y in (1 / 2 : ℝ)..1, (Q y - β / 2) * deriv Q y
        ≤ ∫ y in (1 / 2 : ℝ)..1, (α * (Q y - β / 2) ^ 2 + γ * deriv Q y ^ 2) := by
      apply intervalIntegral.integral_mono_on (by norm_num) hgi' (hf.intervalIntegrable _ _)
      intro y _
      nlinarith [sq_nonneg (α * (Q y - β / 2) - deriv Q y / 2),
        mul_nonneg (sub_nonneg.2 hαγ) (sq_nonneg (deriv Q y)), hα]
    rw [intervalIntegral.integral_neg, hint1, hQh, hQ.zero] at hp1
    rw [hint2, hQh, hQ1] at hp2
    nlinarith [hp1, hp2]
  have hW1eq : W1 0 Q = ∫ y in (0:ℝ)..1, deriv Q y ^ 2 := by
    unfold W1
    simp_rw [hwd0]
  unfold Jag
  rw [hw0, hB2, hW1eq, hQ1]
  have ha : 0 < 3 / 4 + α / 4 := by linarith
  have key : (1 + (β - 1) ^ 2) / 2 + α * (β ^ 2 / 4) + (1 - β / 2) ^ 2 - (2 - 1 / (3 / 4 + α / 4))
      = (3 / 4 + α / 4) * (β - 1 / (3 / 4 + α / 4)) ^ 2 := by
    field_simp
    ring
  nlinarith [hB3, key, mul_nonneg ha.le (sq_nonneg (β - 1 / (3 / 4 + α / 4)))]

/-- 4 C_λ B_λ ≥ 1 (DERIVATION §1.5: C_P B_P > 1/4). -/
lemma ClBl_ge {lam : ℝ} (hlam : 0 ≤ lam) : 1 ≤ 4 * Cl lam * Bl lam := by
  have hC := Cl_ge_one hlam
  have hB := Bl_pos lam
  have hB' := hB.ne'
  have hk2 : Real.sqrt (Cl lam / Bl lam) ^ 2 = Cl lam / Bl lam :=
    Real.sq_sqrt (div_nonneg (by linarith) hB.le)
  have hkpos : 0 < Real.sqrt (Cl lam / Bl lam) := Real.sqrt_pos.2 (div_pos (by linarith) hB)
  have hm := inner_min (Pl_adm lam) (a := Bl lam) (b := Cl lam) (lam := Real.sqrt (Cl lam / Bl lam))
    hB.le (by linarith) (by rw [hk2]; field_simp)
  rw [CP_Pl, BP_Pl] at hm
  have hadd := Cl_add (Real.sqrt (Cl lam / Bl lam))
  have hPld : Pld (Real.sqrt (Cl lam / Bl lam)) 1 = Real.sqrt (Cl lam / Bl lam)
      * Real.cosh (Real.sqrt (Cl lam / Bl lam)) / Real.sinh (Real.sqrt (Cl lam / Bl lam)) := by
    simp [Pld, hkpos.ne']
  have hs := Real.sinh_pos_iff.2 hkpos
  have hcoth : 1 ≤ Real.cosh (Real.sqrt (Cl lam / Bl lam)) / Real.sinh (Real.sqrt (Cl lam / Bl lam)) := by
    rw [le_div_iff₀ hs]
    linarith [Real.sinh_lt_cosh (Real.sqrt (Cl lam / Bl lam))]
  have e1 : Bl lam * Cl (Real.sqrt (Cl lam / Bl lam)) + Cl lam * Bl (Real.sqrt (Cl lam / Bl lam))
      = Bl lam * Pld (Real.sqrt (Cl lam / Bl lam)) 1 := by
    rw [← hadd, hk2]
    field_simp
  have h1 : Bl lam * Real.sqrt (Cl lam / Bl lam) ≤ 2 * (Cl lam * Bl lam) := by
    have h0 : Real.sqrt (Cl lam / Bl lam) ≤ Pld (Real.sqrt (Cl lam / Bl lam)) 1 := by
      rw [hPld, mul_div_assoc]
      nlinarith [hcoth, hkpos]
    have := mul_le_mul_of_nonneg_left h0 hB.le
    nlinarith [hm, e1]
  have h2 : (Bl lam * Real.sqrt (Cl lam / Bl lam)) ^ 2 = Cl lam * Bl lam := by
    rw [mul_pow, hk2]
    field_simp
  have hmpos := mul_pos hB hkpos
  have hm1 : 1 / 2 ≤ Bl lam * Real.sqrt (Cl lam / Bl lam) := by nlinarith [h1, h2, hmpos]
  nlinarith [hm1, h2]

lemma Cl_three : 3 / 2 ≤ Cl 3 := by
  rw [Cl_closed (by norm_num : (3:ℝ) ≠ 0)]
  have hs := Real.sinh_pos_iff.2 (by norm_num : (0:ℝ) < 3)
  have hc : 1 ≤ Real.cosh 3 / Real.sinh 3 := by
    rw [le_div_iff₀ hs]
    linarith [Real.sinh_lt_cosh 3]
  have : 0 ≤ (3:ℝ) ^ 2 / (2 * Real.sinh 3 ^ 2) := by positivity
  nlinarith

lemma Bl_three : Bl 3 ≤ 1 / 5 := by
  rw [Bl_closed (by norm_num : (3:ℝ) ≠ 0), Real.sinh_two_mul]
  have hs := Real.sinh_pos_iff.2 (by norm_num : (0:ℝ) < 3)
  have hE : 4 ≤ Real.exp 3 := by linarith [Real.add_one_le_exp 3]
  have hE0 := Real.exp_pos 3
  have hsc : 5 * Real.cosh 3 ≤ 6 * Real.sinh 3 := by
    rw [Real.cosh_eq, Real.sinh_eq, Real.exp_neg]
    have hinv : Real.exp 3 * (Real.exp 3)⁻¹ = 1 := mul_inv_cancel₀ hE0.ne'
    have hipos := inv_pos.2 hE0
    nlinarith [hinv, hE, mul_le_mul hE hE (by norm_num) hE0.le, hipos]
  rw [div_le_iff₀ (by positivity)]
  nlinarith [mul_le_mul_of_nonneg_left hsc hs.le]

lemma exp_half_le : Real.exp (1 / 2) ≤ 1648722 / 1000000 := by
  have he := Real.exp_one_lt_d9
  have hsq : Real.exp (1 / 2) ^ 2 = Real.exp 1 := by
    rw [← Real.exp_nat_mul]
    norm_num
  have hpos := Real.exp_pos (1 / 2)
  nlinarith [hsq, he, hpos]

/-- §4.1-4.2: TailSmallR holds (R ≤ 1/2, every λ ≥ 0). -/
theorem tailSmallR : TailSmallR := by
  intro R lam hR0 hR hlam Q hQ
  have hC := Cl_ge_one hlam
  have hB := Bl_pos lam
  have hB3 : Bl lam ≤ 1 / 3 := Bl_zero ▸ Bl_anti le_rfl hlam
  have hx : (0:ℝ) < xi0 := by norm_num [xi0]
  have hα : 0 < Cl lam / xi0 := div_pos (by linarith) hx
  have hγ : 0 ≤ xi0 * Bl lam := mul_nonneg hx.le hB.le
  have hαγ : 1 ≤ 4 * (Cl lam / xi0) * (xi0 * Bl lam) := by
    have h := ClBl_ge hlam
    have e : 4 * (Cl lam / xi0) * (xi0 * Bl lam) = 4 * Cl lam * Bl lam := by
      field_simp
    rw [e]
    exact h
  have hcond : xi0 * Bl lam * ((1 - 289 / 10000) / (289 / 10000)) * R ^ 2
      ≤ 289 / 10000 * (Cl lam / xi0) := by
    have hR2 : R ^ 2 ≤ 1 / 4 := by nlinarith
    have hK : (0:ℝ) ≤ xi0 * ((1 - 289 / 10000) / (289 / 10000)) := by norm_num [xi0]
    calc xi0 * Bl lam * ((1 - 289 / 10000) / (289 / 10000)) * R ^ 2
        = (xi0 * ((1 - 289 / 10000) / (289 / 10000))) * (Bl lam * R ^ 2) := by ring
      _ ≤ (xi0 * ((1 - 289 / 10000) / (289 / 10000))) * (1 / 3 * (1 / 4)) :=
          mul_le_mul_of_nonneg_left (mul_le_mul hB3 hR2 (sq_nonneg R) (by norm_num)) hK
      _ ≤ 289 / 10000 * (1 / xi0) := by norm_num [xi0]
      _ ≤ 289 / 10000 * (Cl lam / xi0) := by
          unfold xi0
          linarith
  have h1 := J_R_ge (R := R) hR0.le hα.le hγ (by norm_num) hcond hQ
  have h2 := J0_ge hα hαγ hQ
  have h3 : (562:ℝ) / 331 ≤ 2 - 1 / (3 / 4 + Cl lam / xi0 / 4) := by
    have h31 : (331:ℝ) / 100 ≤ 3 / 4 + Cl lam / xi0 / 4 := by
      unfold xi0
      linarith
    have h32 : 1 / (3 / 4 + Cl lam / xi0 / 4) ≤ 100 / 331 := by
      rw [div_le_div_iff₀ (by linarith) (by norm_num)]
      linarith
    linarith
  have h4 : target R ≤ Real.exp (1 / 2) := by
    unfold target kup
    apply Real.exp_le_exp.2
    nlinarith
  have h5 := exp_half_le
  show target R ≤ Jag R (Cl lam / xi0) (xi0 * Bl lam) Q
  nlinarith [h1, h2, h3, h4, h5]

lemma exp_nineteenth_le : Real.exp (1 / 19) ≤ 19 / 18 := by
  have h := Real.add_one_le_exp (-(1 / 19))
  rw [Real.exp_neg] at h
  have hpos := Real.exp_pos (1 / 19)
  rw [le_div_iff₀ (by norm_num)]
  have hinv : Real.exp (1 / 19) * (Real.exp (1 / 19))⁻¹ = 1 := mul_inv_cancel₀ hpos.ne'
  nlinarith [h, hinv, hpos]

/-- §4.2 strip part: λ ≥ 3, R ∈ [1/2, 19/32]. -/
theorem tailStripSmall :
    ∀ R lam : ℝ, 1 / 2 ≤ R → R ≤ 19 / 32 → 3 ≤ lam → MGe xi0 R lam (target R) := by
  intro R lam hR1 hR2 hlam Q hQ
  have hlam0 : 0 ≤ lam := by linarith
  have hC3 : 3 / 2 ≤ Cl lam := Cl_three.trans (Cl_mono (by norm_num) hlam)
  have hB3 : Bl lam ≤ 1 / 5 := (Bl_anti (by norm_num) hlam).trans Bl_three
  have hB := Bl_pos lam
  have hx : (0:ℝ) < xi0 := by norm_num [xi0]
  have hα : 0 < Cl lam / xi0 := div_pos (by linarith) hx
  have hγ : 0 ≤ xi0 * Bl lam := mul_nonneg hx.le hB.le
  have hαγ : 1 ≤ 4 * (Cl lam / xi0) * (xi0 * Bl lam) := by
    have h := ClBl_ge hlam0
    have e : 4 * (Cl lam / xi0) * (xi0 * Bl lam) = 4 * Cl lam * Bl lam := by
      field_simp
    rw [e]
    exact h
  have hcond : xi0 * Bl lam * ((1 - 1 / 45) / (1 / 45)) * R ^ 2 ≤ 1 / 45 * (Cl lam / xi0) := by
    have hR2' : R ^ 2 ≤ (19 / 32) ^ 2 := by nlinarith
    have hK : (0:ℝ) ≤ xi0 * ((1 - 1 / 45) / (1 / 45)) := by norm_num [xi0]
    calc xi0 * Bl lam * ((1 - 1 / 45) / (1 / 45)) * R ^ 2
        = (xi0 * ((1 - 1 / 45) / (1 / 45))) * (Bl lam * R ^ 2) := by ring
      _ ≤ (xi0 * ((1 - 1 / 45) / (1 / 45))) * (1 / 5 * (19 / 32) ^ 2) :=
          mul_le_mul_of_nonneg_left (mul_le_mul hB3 hR2' (sq_nonneg R) (by norm_num)) hK
      _ ≤ 1 / 45 * (3 / 2 / xi0) := by norm_num [xi0]
      _ ≤ 1 / 45 * (Cl lam / xi0) := by
          unfold xi0
          linarith
  have h1 := J_R_ge (R := R) (by linarith) hα.le hγ (by norm_num) hcond hQ
  have h2 := J0_ge hα hαγ hQ
  have h3 : (818:ℝ) / 459 ≤ 2 - 1 / (3 / 4 + Cl lam / xi0 / 4) := by
    have h31 : (459:ℝ) / 100 ≤ 3 / 4 + Cl lam / xi0 / 4 := by
      unfold xi0
      linarith
    have h32 : 1 / (3 / 4 + Cl lam / xi0 / 4) ≤ 100 / 459 := by
      rw [div_le_div_iff₀ (by linarith) (by norm_num)]
      linarith
    linarith
  have h4 : target R ≤ Real.exp (1 / 2) * Real.exp (1 / 19) := by
    rw [← Real.exp_add]
    unfold target kup
    apply Real.exp_le_exp.2
    nlinarith
  have h5 := exp_half_le
  have h6 := exp_nineteenth_le
  have h7 : Real.exp (1 / 2) * Real.exp (1 / 19) ≤ 1648722 / 1000000 * (19 / 18) :=
    mul_le_mul h5 h6 (Real.exp_pos _).le (by norm_num)
  show target R ≤ Jag R (Cl lam / xi0) (xi0 * Bl lam) Q
  nlinarith [h1, h2, h3, h4, h7]

/-- Main theorem with the analytic tails proved: from Cover, CertCells (1276 JSON cells), the crude
strip inequality on [19/32, 28] and the crude inequality for R ≥ 28 (all checked by ceil.py, Arb). -/
theorem ceiling_bernard_final
    (hcov : Cover jsonCells) (hcert : CertCells jsonCells)
    (hstrip : ∀ R : ℝ, 19 / 32 ≤ R → R ≤ 28 → target R ≤ crude R (Cl 3 / xi0 * G R))
    (hR : TailLargeR)
    (P : Polynomial ℝ) (hP0 : P.eval 0 = 0) (hP1 : P.eval 1 = 1) (M : ℕ) (lam : ℕ → ℝ)
    {R ξ : ℝ} (hR0 : 0 < R) (hξ : 0 < ξ) (hξ0 : ξ ≤ 25 / 256) :
    kappaC (fun x => P.eval x) (QB M lam) R ξ ≤ 0.07005633 :=
  ceiling_bernard_json hcov hcert tailSmallR ⟨tailStripSmall, hstrip⟩ hR P hP0 hP1 M lam hR0 hξ hξ0

/-- Same, real C¹ (in particular real polynomial) P, Q of the class. -/
theorem ceiling_final (cells : List Cell) (hcov : Cover cells) (hcert : CertCells cells)
    (hstrip : ∀ R : ℝ, 19 / 32 ≤ R → R ≤ 28 → target R ≤ crude R (Cl 3 / xi0 * G R))
    (hR : TailLargeR) {P Q : ℝ → ℝ} (hP : PAdm P) (hQ : QAdm Q) {R ξ : ℝ} (hR0 : 0 < R)
    (hξ : 0 < ξ) (hξ0 : ξ ≤ 25 / 256) : kappa P Q R ξ ≤ 0.07005633 :=
  ceiling cells hcov hcert tailSmallR ⟨tailStripSmall, hstrip⟩ hR hP hQ hR0 hξ hξ0

end GL2Ceil




