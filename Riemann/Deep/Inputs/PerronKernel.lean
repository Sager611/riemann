import Riemann.Deep.WeilEF.Z23.FromPNTPlus.ResidueCalcOnRectangles

/-!
# Truncated Perron kernel (Apostol, Lemma 11.12, explicit form)

For L ≠ 0, c > 0, T > 0:
‖(1/2π) ∫_{-T}^{T} e^{L(c+iy)}/(c+iy) dy - [L > 0]‖ ≤ 2 e^{cL} / (π |L| T).

Proof: rectangle contour. L > 0: rectangle [-T, c] × [-T, T], simple pole at 0 with residue 1
(`ResidueTheoremOnRectangleWithSimplePole'`). L < 0: rectangle [c, c+T] × [-T, T], no pole
(`HolomorphicOn.vanishesOnRectangle`). The three other sides are bounded directly.
-/

open Complex Set MeasureTheory intervalIntegral Filter Topology
open scoped Real Interval

namespace Fayed.InputsProof.PerronAux

noncomputable section

/-- The Perron kernel e^{Lz}/z. -/
def kern (L : ℝ) (z : ℂ) : ℂ := cexp (L * z) / z

lemma norm_kern (L x y : ℝ) :
    ‖kern L (x + y * I)‖ = Real.exp (L * x) / ‖(x : ℂ) + y * I‖ := by
  rw [kern, norm_div, Complex.norm_exp]
  congr 2
  simp

lemma norm_kern_le_im {L x y T : ℝ} (hT : 0 < T) (hy : T ≤ |y|) :
    ‖kern L (x + y * I)‖ ≤ Real.exp (L * x) / T := by
  rw [norm_kern]
  apply div_le_div_of_nonneg_left (Real.exp_pos _).le hT
  refine hy.trans ?_
  simpa using abs_im_le_norm ((x : ℂ) + y * I)

lemma norm_kern_le_re {L x y : ℝ} (hx : 0 < |x|) :
    ‖kern L (x + y * I)‖ ≤ Real.exp (L * x) / |x| := by
  rw [norm_kern]
  apply div_le_div_of_nonneg_left (Real.exp_pos _).le hx
  simpa using abs_re_le_norm ((x : ℂ) + y * I)

lemma integral_exp_mul {L a b : ℝ} (hL : L ≠ 0) :
    ∫ x in a..b, Real.exp (L * x) = (Real.exp (L * b) - Real.exp (L * a)) / L := by
  rw [intervalIntegral.integral_comp_mul_left (fun x => Real.exp x) hL, integral_exp,
    smul_eq_mul]
  field_simp

lemma norm_hint_le {L a b y T : ℝ} (hL : L ≠ 0) (hab : a ≤ b) (hT : 0 < T) (hy : T ≤ |y|) :
    ‖∫ x in a..b, kern L (x + y * I)‖ ≤ (Real.exp (L * b) - Real.exp (L * a)) / (L * T) := by
  have h := intervalIntegral.norm_integral_le_of_norm_le (μ := volume) hab
    (g := fun x => Real.exp (L * x) / T)
    (ae_of_all _ (fun x _ => norm_kern_le_im (x := x) (L := L) hT hy))
    (by
      apply Continuous.intervalIntegrable
      fun_prop)
  refine h.trans (le_of_eq ?_)
  rw [intervalIntegral.integral_div, integral_exp_mul hL]
  field_simp

lemma norm_vint_le {L X T : ℝ} (hX : 0 < |X|) (hT : 0 < T) :
    ‖∫ y in (-T)..T, kern L (X + y * I)‖ ≤ 2 * T * Real.exp (L * X) / |X| := by
  have h := intervalIntegral.norm_integral_le_of_norm_le_const (a := -T) (b := T)
    (C := Real.exp (L * X) / |X|) (f := fun y => kern L (X + y * I))
    (fun y _ => norm_kern_le_re hX)
  refine h.trans (le_of_eq ?_)
  rw [abs_of_pos (show (0:ℝ) < T - -T by linarith)]
  ring

lemma exp_neg_le_inv {a : ℝ} (ha : 0 < a) : Real.exp (-a) ≤ 1 / a := by
  rw [Real.exp_neg, one_div]
  exact inv_anti₀ ha (by linarith [Real.add_one_le_exp a])

lemma kern_differentiableAt (L : ℝ) {s : ℂ} (hs : s ≠ 0) : DifferentiableAt ℂ (kern L) s := by
  unfold kern
  exact DifferentiableAt.div (by fun_prop) differentiableAt_id hs

theorem perron_kernel {L c T : ℝ} (hL : L ≠ 0) (hc : 0 < c) (hT : 0 < T) :
    ‖((1 / (2 * π) : ℝ) : ℂ) * (∫ y in (-T)..T, kern L ((c : ℂ) + y * I)) -
        (if 0 < L then 1 else 0)‖ ≤ 2 * Real.exp (c * L) / (π * |L| * T) := by
  have hπ : 0 < π := Real.pi_pos
  set K := Real.exp (c * L) / (|L| * T) with hK
  have hLabs : 0 < |L| := abs_pos.mpr hL
  have hgoal : 2 * Real.exp (c * L) / (π * |L| * T) = 1 / (2 * π) * (4 * K) := by
    rw [hK]; field_simp; ring
  rw [hgoal]
  have hnorm : ∀ V : ℂ, ‖((1 / (2 * π) : ℝ) : ℂ) * V‖ = 1 / (2 * π) * ‖V‖ := by
    intro V
    rw [norm_mul, Complex.norm_real, Real.norm_eq_abs, abs_of_pos (by positivity)]
  rcases lt_or_gt_of_ne hL with hneg | hpos
  · -- L < 0: no pole in [c, c+T] × [-T, T]
    rw [ite_eq_right (not_lt.mpr hneg.le), sub_zero, hnorm]
    have hhol : HolomorphicOn (kern L) {s : ℂ | 0 < s.re} := by
      intro s hs
      have : s ≠ 0 := by intro h; rw [h] at hs; simp at hs
      exact (kern_differentiableAt L this).differentiableWithinAt
    have hsub : Rectangle (⟨c, -T⟩ : ℂ) ⟨c + T, T⟩ ⊆ {s : ℂ | 0 < s.re} := by
      intro s hs
      have h1 : s.re ∈ [[c, c + T]] := hs.1
      rw [uIcc_of_le (by linarith)] at h1
      exact lt_of_lt_of_le hc h1.1
    have h0 := hhol.vanishesOnRectangle hsub
    simp only [RectangleIntegral, HIntegral, VIntegral, smul_eq_mul] at h0
    set Hb := ∫ x in c..c + T, kern L (x + (-T : ℝ) * I)
    set Ht := ∫ x in c..c + T, kern L (x + (T : ℝ) * I)
    set VR := ∫ y in (-T)..T, kern L (((c + T : ℝ) : ℂ) + y * I)
    set Vc := ∫ y in (-T)..T, kern L ((c : ℂ) + y * I)
    have hVc : Vc = VR - I * Hb + I * Ht := by linear_combination I * h0 + (Vc - VR) * I_sq
    have hb1 : ‖Hb‖ ≤ K := by
      refine (norm_hint_le hL (by linarith) hT (by rw [abs_neg, abs_of_pos hT])).trans ?_
      have e : (Real.exp (L * (c + T)) - Real.exp (L * c)) / (L * T) =
          (Real.exp (L * c) - Real.exp (L * (c + T))) / (|L| * T) := by
        rw [abs_of_neg hneg]; field_simp; ring
      rw [e, hK, mul_comm L c]
      exact div_le_div_of_nonneg_right (by linarith [Real.exp_pos (L * (c + T))]) (by positivity)
    have hb2 : ‖Ht‖ ≤ K := by
      refine (norm_hint_le hL (by linarith) hT (by rw [abs_of_pos hT])).trans ?_
      have e : (Real.exp (L * (c + T)) - Real.exp (L * c)) / (L * T) =
          (Real.exp (L * c) - Real.exp (L * (c + T))) / (|L| * T) := by
        rw [abs_of_neg hneg]; field_simp; ring
      rw [e, hK, mul_comm L c]
      exact div_le_div_of_nonneg_right (by linarith [Real.exp_pos (L * (c + T))]) (by positivity)
    have hb3 : ‖VR‖ ≤ 2 * K := by
      have hcT : 0 < c + T := by linarith
      refine (norm_vint_le (by rw [abs_of_pos hcT]; exact hcT) hT).trans ?_
      rw [abs_of_pos hcT]
      have hE : Real.exp (L * (c + T)) = Real.exp (c * L) * Real.exp (-(|L| * T)) := by
        rw [← Real.exp_add, abs_of_neg hneg]; ring_nf
      have hE2 : Real.exp (-(|L| * T)) ≤ 1 / (|L| * T) := exp_neg_le_inv (by positivity)
      have step1 : 2 * T * Real.exp (L * (c + T)) / (c + T) ≤ 2 * Real.exp (L * (c + T)) := by
        rw [div_le_iff₀ hcT]
        nlinarith [Real.exp_pos (L * (c + T))]
      refine step1.trans ?_
      rw [hE, hK, ← mul_one_div (Real.exp (c * L))]
      nlinarith [Real.exp_pos (c * L)]
    rw [hVc]
    have htri : ‖VR - I * Hb + I * Ht‖ ≤ ‖VR‖ + ‖Hb‖ + ‖Ht‖ := by
      refine (norm_add_le _ _).trans ?_
      refine (add_le_add (norm_sub_le _ _) le_rfl).trans ?_
      simp
    have : 0 ≤ 1 / (2 * π) := by positivity
    calc 1 / (2 * π) * ‖VR - I * Hb + I * Ht‖ ≤ 1 / (2 * π) * (‖VR‖ + ‖Hb‖ + ‖Ht‖) :=
          mul_le_mul_of_nonneg_left htri this
      _ ≤ 1 / (2 * π) * (4 * K) := by
          apply mul_le_mul_of_nonneg_left _ this
          linarith
  · -- L > 0: simple pole at 0 with residue 1 inside [-T, c] × [-T, T]
    simp only [hpos, ↓reduceIte]
    have hres : RectangleIntegral' (kern L) (⟨-T, -T⟩ : ℂ) ⟨c, T⟩ = 1 := by
      refine ResidueTheoremOnRectangleWithSimplePole' (p := 0) (A := 1)
        (by simp only; linarith) (by simp only; linarith) ?_ ?_ ?_
      · rw [rectangle_mem_nhds_iff]
        simp only [mem_reProdIm, zero_re, zero_im]
        rw [uIoo_of_le (by linarith), uIoo_of_le (by linarith)]
        exact ⟨⟨by linarith, hc⟩, ⟨by linarith, hT⟩⟩
      · intro s hs
        exact (kern_differentiableAt L hs.2).differentiableWithinAt
      · have hd : HasDerivAt (fun s : ℂ => cexp (L * s)) (cexp (L * 0) * (L * 1)) 0 :=
          ((hasDerivAt_id (0 : ℂ)).const_mul (L : ℂ)).cexp
        have hT' : Tendsto (kern L - fun s => 1 / (s - 0)) (𝓝[≠] 0)
            (𝓝 (cexp (L * 0) * (L * 1))) := by
          refine hd.tendsto_slope_zero.congr (fun t => ?_)
          simp only [kern, zero_add, mul_zero, Complex.exp_zero, smul_eq_mul, Pi.sub_apply,
            sub_zero]
          ring
        exact hT'.isBigO_one ℂ
    simp only [RectangleIntegral', RectangleIntegral, HIntegral, VIntegral, smul_eq_mul] at hres
    set Hb := ∫ x in (-T)..c, kern L (x + (-T : ℝ) * I)
    set Ht := ∫ x in (-T)..c, kern L (x + (T : ℝ) * I)
    set VL := ∫ y in (-T)..T, kern L (((-T : ℝ) : ℂ) + y * I)
    set Vc := ∫ y in (-T)..T, kern L ((c : ℂ) + y * I)
    have hπ' : (π : ℂ) ≠ 0 := by exact_mod_cast hπ.ne'
    have h' : Hb - Ht + I * Vc - I * VL = 2 * π * I := by
      field_simp at hres
      linear_combination hres
    have hVc : ((1 / (2 * π) : ℝ) : ℂ) * Vc - 1 =
        ((1 / (2 * π) : ℝ) : ℂ) * (VL + I * Hb - I * Ht) := by
      push_cast
      field_simp
      linear_combination (-I) * h' + (Vc - VL) * I_sq + (-2 * (π : ℂ)) * I_sq
    rw [hVc, hnorm]
    have hb1 : ‖Hb‖ ≤ K := by
      refine (norm_hint_le hL (by linarith) hT (by rw [abs_neg, abs_of_pos hT])).trans ?_
      rw [hK, abs_of_pos hpos, mul_comm L c]
      exact div_le_div_of_nonneg_right (by linarith [Real.exp_pos (L * -T)]) (by positivity)
    have hb2 : ‖Ht‖ ≤ K := by
      refine (norm_hint_le hL (by linarith) hT (by rw [abs_of_pos hT])).trans ?_
      rw [hK, abs_of_pos hpos, mul_comm L c]
      exact div_le_div_of_nonneg_right (by linarith [Real.exp_pos (L * -T)]) (by positivity)
    have hb3 : ‖VL‖ ≤ 2 * K := by
      refine (norm_vint_le (by rw [abs_neg, abs_of_pos hT]; exact hT) hT).trans ?_
      rw [abs_neg, abs_of_pos hT]
      have hE2 : Real.exp (-(L * T)) ≤ 1 / (L * T) := exp_neg_le_inv (by positivity)
      have hc1 : 1 ≤ Real.exp (c * L) := Real.one_le_exp (by positivity)
      have e : 2 * T * Real.exp (L * -T) / T = 2 * Real.exp (-(L * T)) := by
        field_simp
      rw [e, hK, abs_of_pos hpos]
      have : 1 / (L * T) ≤ Real.exp (c * L) / (L * T) :=
        div_le_div_of_nonneg_right hc1 (by positivity)
      linarith
    have htri : ‖VL + I * Hb - I * Ht‖ ≤ ‖VL‖ + ‖Hb‖ + ‖Ht‖ := by
      refine (norm_sub_le _ _).trans ?_
      refine (add_le_add (norm_add_le _ _) le_rfl).trans ?_
      simp
    have : 0 ≤ 1 / (2 * π) := by positivity
    calc 1 / (2 * π) * ‖VL + I * Hb - I * Ht‖ ≤ 1 / (2 * π) * (‖VL‖ + ‖Hb‖ + ‖Ht‖) :=
          mul_le_mul_of_nonneg_left htri this
      _ ≤ 1 / (2 * π) * (4 * K) := by
          apply mul_le_mul_of_nonneg_left _ this
          linarith

end

end Fayed.InputsProof.PerronAux
