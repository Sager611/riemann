import Riemann.Deep.FirstZero.GammaRatio
import Mathlib.Analysis.SpecialFunctions.Gamma.Deligne
import Mathlib.Analysis.Real.Pi.Bounds

/-!
# Phase stability of Γ on a short vertical segment, any real part > 0

Generalises Riemann/Deep/FirstZero/GammaRatio.lean (Re z = 1/4, Im z ≥ 7) to any z with Re z > 0,
given lower bounds ‖z + iδ + j‖ ≥ (j+1)/2 and ‖z + iδ + j‖ ≥ m0 > 0. Used at z = 3/4 + i t/2
for the Γℝ(s + 1) factor of the odd character χ₋₄.
-/

open Complex Filter Topology

namespace Fayed.Beta

open Fayed.FirstZero

/-- One paired factor of the Gamma-ratio product. -/
lemma pair_bound' (z : ℂ) (hwz : ∀ j : ℕ, ((j : ℝ) + 1) / 2 ≤ ‖z + j‖) {δ : ℝ} (hδ0 : 0 ≤ δ)
    (hδ1 : δ ≤ 1) (j : ℕ) :
    ‖(1 - δ * I / (z + j)) *
        cexp (δ * I * ((Real.log ((j : ℝ) + 2) - Real.log ((j : ℝ) + 1) : ℝ) : ℂ)) - 1‖
      ≤ (3 * δ ^ 2 + δ + 2 * δ * ‖z - 1‖) / ((j : ℝ) + 1) ^ 2 := by
  set w : ℂ := z + j with hw
  set L : ℝ := Real.log ((j : ℝ) + 2) - Real.log ((j : ℝ) + 1) with hL
  obtain ⟨hL0, hL1⟩ := log_succ_ratio_bounds j
  rw [← hL] at hL0 hL1
  have hj1 : (0 : ℝ) < (j : ℝ) + 1 := by positivity
  have hwn : ((j : ℝ) + 1) / 2 ≤ ‖w‖ := hwz j
  have hwpos : 0 < ‖w‖ := lt_of_lt_of_le (by positivity) hwn
  have hw0 : w ≠ 0 := norm_pos_iff.mp hwpos
  have hLpos : 0 ≤ L := le_trans (by positivity) hL0
  set x : ℂ := δ * I * (L : ℂ) with hx
  have hxn : ‖x‖ = δ * L := by
    simp only [hx, norm_mul, Complex.norm_real, Real.norm_eq_abs, Complex.norm_I, mul_one,
      abs_of_nonneg hδ0, abs_of_nonneg hLpos]
  have hdL : δ * L ≤ δ / ((j : ℝ) + 1) := by
    rw [div_eq_mul_one_div]; exact mul_le_mul_of_nonneg_left hL1 hδ0
  have hdL1 : δ * L ≤ 1 := by
    refine hdL.trans ?_
    rw [div_le_one hj1]; linarith [(Nat.cast_nonneg j : (0 : ℝ) ≤ j)]
  -- the pieces
  have hA : ‖δ * I / w‖ ≤ 2 * δ / ((j : ℝ) + 1) := by
    rw [norm_div, norm_mul, Complex.norm_real, Real.norm_eq_abs, Complex.norm_I, mul_one,
      abs_of_nonneg hδ0, div_le_div_iff₀ hwpos hj1]
    nlinarith
  have hB : ‖cexp x - 1‖ ≤ δ / ((j : ℝ) + 1) := by
    have hx' : x = I * ((δ * L : ℝ) : ℂ) := by rw [hx]; push_cast; ring
    rw [hx']
    refine (Real.norm_exp_I_mul_ofReal_sub_one_le).trans ?_
    rw [Real.norm_eq_abs, abs_of_nonneg (by positivity)]
    exact hdL
  have hBx : ‖cexp x - 1 - x‖ ≤ (δ / ((j : ℝ) + 1)) ^ 2 := by
    refine (norm_exp_sub_one_sub_id_le (by rw [hxn]; exact hdL1)).trans ?_
    rw [hxn]
    exact pow_le_pow_left₀ (by positivity) hdL 2
  have hLd : |L - 1 / ((j : ℝ) + 1)| ≤ 1 / ((j : ℝ) + 1) ^ 2 := by
    have e1 : 1 / ((j : ℝ) + 1) - 1 / ((j : ℝ) + 2) = 1 / (((j : ℝ) + 1) * ((j : ℝ) + 2)) := by
      field_simp; ring
    have e2 : 1 / (((j : ℝ) + 1) * ((j : ℝ) + 2)) ≤ 1 / ((j : ℝ) + 1) ^ 2 :=
      one_div_le_one_div_of_le (by positivity) (by nlinarith)
    rw [abs_le]; constructor <;> linarith
  have hC : ‖δ * I * ((L : ℂ) - 1 / ((j : ℂ) + 1))‖ ≤ δ / ((j : ℝ) + 1) ^ 2 := by
    have hcast : (L : ℂ) - 1 / ((j : ℂ) + 1) = ((L - 1 / ((j : ℝ) + 1) : ℝ) : ℂ) := by push_cast; ring
    rw [hcast, norm_mul, norm_mul, Complex.norm_real, Complex.norm_real, Real.norm_eq_abs,
      Real.norm_eq_abs, Complex.norm_I, mul_one, abs_of_nonneg hδ0]
    calc δ * |L - 1 / ((j : ℝ) + 1)| ≤ δ * (1 / ((j : ℝ) + 1) ^ 2) :=
          mul_le_mul_of_nonneg_left hLd hδ0
      _ = _ := by ring
  have hD : ‖δ * I * (1 / ((j : ℂ) + 1) - 1 / w)‖ ≤ 2 * δ * ‖z - 1‖ / ((j : ℝ) + 1) ^ 2 := by
    have hj0 : ((j : ℂ) + 1) ≠ 0 := by
      have : ((j : ℂ) + 1) = (((j : ℝ) + 1 : ℝ) : ℂ) := by push_cast; ring
      rw [this]; exact_mod_cast hj1.ne'
    have e : 1 / ((j : ℂ) + 1) - 1 / w = (z - 1) / (((j : ℂ) + 1) * w) := by
      rw [div_sub_div _ _ hj0 hw0, hw]; ring
    have hn : ‖(j : ℂ) + 1‖ = (j : ℝ) + 1 := by
      have : ((j : ℂ) + 1) = (((j : ℝ) + 1 : ℝ) : ℂ) := by push_cast; ring
      rw [this, Complex.norm_real, Real.norm_eq_abs, abs_of_pos hj1]
    rw [e, norm_mul, norm_mul, norm_div, norm_mul, hn, Complex.norm_real, Real.norm_eq_abs,
      Complex.norm_I, mul_one, abs_of_nonneg hδ0]
    have hz := norm_nonneg (z - 1)
    have h2 : ‖z - 1‖ / (((j : ℝ) + 1) * ‖w‖) ≤ 2 * ‖z - 1‖ / ((j : ℝ) + 1) ^ 2 := by
      rw [div_le_div_iff₀ (by positivity) (by positivity)]
      have h3 : ‖z - 1‖ * ((j : ℝ) + 1) * ((j : ℝ) + 1) ≤ ‖z - 1‖ * ((j : ℝ) + 1) * (2 * ‖w‖) :=
        mul_le_mul_of_nonneg_left (by linarith) (by positivity)
      nlinarith
    calc δ * (‖z - 1‖ / (((j : ℝ) + 1) * ‖w‖)) ≤ δ * (2 * ‖z - 1‖ / ((j : ℝ) + 1) ^ 2) :=
          mul_le_mul_of_nonneg_left h2 hδ0
      _ = _ := by ring
  have hdec : (1 - δ * I / w) * cexp x - 1 =
      (-(δ * I / w)) * (cexp x - 1) + (cexp x - 1 - x) + δ * I * ((L : ℂ) - 1 / ((j : ℂ) + 1))
        + δ * I * (1 / ((j : ℂ) + 1) - 1 / w) := by
    rw [hx]; ring
  rw [hdec]
  have hAB : ‖(-(δ * I / w)) * (cexp x - 1)‖ ≤ 2 * δ ^ 2 / ((j : ℝ) + 1) ^ 2 := by
    rw [norm_mul, norm_neg]
    calc ‖δ * I / w‖ * ‖cexp x - 1‖ ≤ (2 * δ / ((j : ℝ) + 1)) * (δ / ((j : ℝ) + 1)) :=
          mul_le_mul hA hB (norm_nonneg _) (by positivity)
      _ = _ := by field_simp
  have hBx' : (δ / ((j : ℝ) + 1)) ^ 2 = δ ^ 2 / ((j : ℝ) + 1) ^ 2 := by rw [div_pow]
  calc _ ≤ ‖(-(δ * I / w)) * (cexp x - 1)‖ + ‖cexp x - 1 - x‖
        + ‖δ * I * ((L : ℂ) - 1 / ((j : ℂ) + 1))‖ + ‖δ * I * (1 / ((j : ℂ) + 1) - 1 / w)‖ := by
        refine (norm_add_le _ _).trans ?_
        gcongr
        refine (norm_add_le _ _).trans ?_
        gcongr
        exact norm_add_le _ _
    _ ≤ 2 * δ ^ 2 / ((j : ℝ) + 1) ^ 2 + δ ^ 2 / ((j : ℝ) + 1) ^ 2 + δ / ((j : ℝ) + 1) ^ 2
        + 2 * δ * ‖z - 1‖ / ((j : ℝ) + 1) ^ 2 := by
        rw [← hBx']; gcongr
    _ = _ := by ring


/-- **Γ moves little along a short vertical segment.** -/
theorem Gamma_ratio_gen' (z : ℂ) (hzre : 0 < z.re) {δ : ℝ} (hδ0 : 0 ≤ δ) (hδ : δ ≤ 1)
    (hw : ∀ j : ℕ, ((j : ℝ) + 1) / 2 ≤ ‖z + δ * I + j‖) {m0 : ℝ} (hm0 : 0 < m0)
    (hlow : ∀ j : ℕ, m0 ≤ ‖z + δ * I + j‖) {B : ℝ} (hB : ‖z + δ * I - 1‖ ≤ B) {S : ℝ}
    (hS : 2 * (3 * δ ^ 2 + δ + 2 * δ * B) + 2 * (δ / m0) ≤ S) :
    ‖Gamma (z + δ * I) / Gamma z - 1‖ ≤ Real.exp S - 1 := by
  obtain ⟨z2, hz2⟩ : ∃ z2, z2 = z + δ * I := ⟨_, rfl⟩
  rw [← hz2]
  rw [← hz2] at hw hlow hB
  have hG : Gamma z ≠ 0 := Gamma_ne_zero_of_re_pos hzre
  have hlim := (((GammaSeq_tendsto_Gamma z2).div (GammaSeq_tendsto_Gamma z) hG).sub_const 1).norm
  refine le_of_tendsto hlim (Filter.eventually_atTop.2 ⟨1, fun n hn => ?_⟩)
  obtain ⟨m, rfl⟩ : ∃ m, n = m + 1 := ⟨n - 1, by omega⟩
  rw [Pi.div_apply, GammaSeq_ratio z z2 hzre δ hz2 m]
  set K := 3 * δ ^ 2 + δ + 2 * δ * ‖z2 - 1‖ with hK
  set e : ℕ → ℂ := fun j => (1 - δ * I / (z2 + j)) *
    cexp (δ * I * ((Real.log ((j : ℝ) + 2) - Real.log ((j : ℝ) + 1) : ℝ) : ℂ)) - 1 with he_def
  have he : ∀ j, ‖e j‖ ≤ K / ((j : ℝ) + 1) ^ 2 :=
    fun j => pair_bound' z2 hw hδ0 hδ j
  have hK0 : 0 ≤ K := by positivity
  have hS2 : ∑ j ∈ Finset.range m, ‖e j‖ ≤ 2 * K := by
    calc ∑ j ∈ Finset.range m, ‖e j‖ ≤ ∑ j ∈ Finset.range m, K / ((j : ℝ) + 1) ^ 2 :=
          Finset.sum_le_sum fun j _ => he j
      _ = K * ∑ j ∈ Finset.range m, 1 / ((j : ℝ) + 1) ^ 2 := by
          rw [Finset.mul_sum]; refine Finset.sum_congr rfl fun j _ => ?_; ring
      _ ≤ K * 2 := mul_le_mul_of_nonneg_left (sum_range_inv_sq_le_two m) hK0
      _ = 2 * K := by ring
  have hX := norm_prod_one_add_sub_one_le e m
  have hprod_eq : ∏ j ∈ Finset.range m, ((1 - δ * I / (z2 + j)) *
      cexp (δ * I * ((Real.log ((j : ℝ) + 2) - Real.log ((j : ℝ) + 1) : ℝ) : ℂ))) =
      ∏ j ∈ Finset.range m, (1 + e j) := by
    refine Finset.prod_congr rfl fun j _ => ?_
    simp only [he_def]; ring
  rw [hprod_eq]
  have ha : ∀ k : ℕ, ‖-(δ * I / (z2 + k))‖ ≤ δ / m0 := by
    intro k
    rw [norm_neg, norm_div, norm_mul, Complex.norm_real, Real.norm_eq_abs, Complex.norm_I, mul_one,
      abs_of_nonneg hδ0]
    exact div_le_div_of_nonneg_left hδ0 hm0 (hlow k)
  have step1 := norm_mul_one_add_sub_one_le (e := -(δ * I / (z2 + m))) hX
  have step2 := norm_mul_one_add_sub_one_le (e := -(δ * I / (z2 + ((m + 1 : ℕ) : ℂ)))) step1
  have hfin : (∏ j ∈ Finset.range m, (1 + e j)) * (1 - δ * I / (z2 + m)) *
      (1 - δ * I / (z2 + ((m + 1 : ℕ) : ℂ))) =
      ((∏ j ∈ Finset.range m, (1 + e j)) * (1 + -(δ * I / (z2 + m)))) *
        (1 + -(δ * I / (z2 + ((m + 1 : ℕ) : ℂ)))) := by ring
  rw [hfin]
  refine step2.trans ?_
  have h1 := ha m
  have h2 := ha (m + 1)
  have hKb : K ≤ 3 * δ ^ 2 + δ + 2 * δ * B := by
    rw [hK]; nlinarith [mul_le_mul_of_nonneg_left hB hδ0]
  have hsum : ∑ j ∈ Finset.range m, ‖e j‖ + ‖-(δ * I / (z2 + m))‖
      + ‖-(δ * I / (z2 + ((m + 1 : ℕ) : ℂ)))‖ ≤ S := by
    push_cast at h2 ⊢
    linarith
  have hmono := Real.exp_le_exp.mpr hsum
  linarith

/-- e^{iδ log 4}·Γℝ(s₂+1)/Γℝ(s₁+1) is close to 1 for s_k = 1/2 + i t_k, t₂ = t₁ + 2δ. -/
lemma GammaR1_ratio {t1 t2 δ : ℝ} (ht1 : 6 ≤ t1) (ht2 : t2 = t1 + 2 * δ) (hδ0 : 0 ≤ δ)
    (hδ : δ ≤ 1) {S : ℝ} (hS : 2 * (3 * δ ^ 2 + δ + 2 * δ * (1 / 4 + t1 / 2 + δ)) + 2 * (δ / 3) ≤ S)
    (hS0 : 0 < S) (hS1 : S < 1) {c : ℝ}
    (hc : 3 / 10 * δ * (1 / (1 - S)) + (1 / (1 - S) - 1) ≤ c) :
    ‖cexp (I * ((δ * Real.log 4 : ℝ) : ℂ)) *
      (Gammaℝ ((1 / 2 : ℂ) + (t2 : ℂ) * I + 1) / Gammaℝ ((1 / 2 : ℂ) + (t1 : ℂ) * I + 1)) - 1‖
      ≤ c := by
  set z : ℂ := ((1 / 2 : ℂ) + (t1 : ℂ) * I + 1) / 2 with hz
  have hs2 : ((1 / 2 : ℂ) + (t2 : ℂ) * I + 1) / 2 = z + δ * I := by
    rw [hz, ht2]; push_cast; ring
  have hzre : z.re = 3 / 4 := by rw [hz]; simp; norm_num
  have hzim : z.im = t1 / 2 := by rw [hz]; simp
  have hz2re : ∀ j : ℕ, (z + δ * I + j).re = 3 / 4 + j := fun j => by simp [hzre]
  have hz2im : ∀ j : ℕ, (z + δ * I + j).im = t1 / 2 + δ := fun j => by simp [hzim]
  have hw : ∀ j : ℕ, ((j : ℝ) + 1) / 2 ≤ ‖z + δ * I + j‖ := by
    intro j
    have h := abs_re_le_norm (z + δ * I + j)
    rw [hz2re, abs_of_nonneg (by positivity)] at h
    have : (0 : ℝ) ≤ j := Nat.cast_nonneg j
    linarith
  have hlow : ∀ j : ℕ, 3 ≤ ‖z + δ * I + j‖ := by
    intro j
    have h := abs_im_le_norm (z + δ * I + j)
    rw [hz2im, abs_of_nonneg (by linarith)] at h
    linarith
  have hB : ‖z + δ * I - 1‖ ≤ 1 / 4 + t1 / 2 + δ := by
    refine (norm_le_abs_re_add_abs_im _).trans ?_
    have h1 : (z + δ * I - 1).re = -1 / 4 := by simp [hzre]; norm_num
    have h2 : (z + δ * I - 1).im = t1 / 2 + δ := by simp [hzim]
    rw [h1, h2, abs_of_nonneg (show (0 : ℝ) ≤ t1 / 2 + δ by linarith)]
    norm_num; linarith
  have hG0 := Gamma_ratio_gen' z (by rw [hzre]; norm_num) hδ0 hδ hw (by norm_num) hlow hB hS
  have hexp : Real.exp S < 1 / (1 - S) := Real.exp_bound_div_one_sub_of_interval' hS0 hS1
  have hG : ‖Gamma (z + δ * I) / Gamma z - 1‖ ≤ 1 / (1 - S) - 1 := by linarith
  have hpi0 : ((Real.pi : ℝ) : ℂ) ≠ 0 := by exact_mod_cast Real.pi_ne_zero
  have hP : (Real.pi : ℂ) ^ (-((1 / 2 : ℂ) + (t2 : ℂ) * I + 1) / 2) =
      (Real.pi : ℂ) ^ (-((1 / 2 : ℂ) + (t1 : ℂ) * I + 1) / 2) *
        cexp (I * ((-(δ * Real.log Real.pi)) : ℝ)) := by
    have hPI : (Real.pi : ℂ) ^ (-((δ : ℂ) * I)) = cexp (I * ((-(δ * Real.log Real.pi)) : ℝ)) := by
      rw [cpow_def_of_ne_zero hpi0, ← ofReal_log Real.pi_pos.le]
      congr 1; push_cast; ring
    rw [show -((1 / 2 : ℂ) + (t2 : ℂ) * I + 1) / 2 =
      -((1 / 2 : ℂ) + (t1 : ℂ) * I + 1) / 2 + (-((δ : ℂ) * I)) by
        rw [ht2]; push_cast; ring,
      cpow_add _ _ hpi0, hPI]
  have hGne : Gamma z ≠ 0 := Gamma_ne_zero_of_re_pos (by rw [hzre]; norm_num)
  have hPne : (Real.pi : ℂ) ^ (-((1 / 2 : ℂ) + (t1 : ℂ) * I + 1) / 2) ≠ 0 := by
    rw [Ne, cpow_eq_zero_iff]; tauto
  have hratio : cexp (I * ((δ * Real.log 4 : ℝ) : ℂ)) *
      (Gammaℝ ((1 / 2 : ℂ) + (t2 : ℂ) * I + 1) / Gammaℝ ((1 / 2 : ℂ) + (t1 : ℂ) * I + 1)) =
      cexp (I * ((δ * (Real.log 4 - Real.log Real.pi) : ℝ) : ℂ)) * (Gamma (z + δ * I) / Gamma z) := by
    simp only [Gammaℝ]
    rw [hP, hs2, ← hz]
    have e2 : cexp (I * ((δ * (Real.log 4 - Real.log Real.pi) : ℝ) : ℂ)) =
        cexp (I * ((δ * Real.log 4 : ℝ) : ℂ)) * cexp (I * ((-(δ * Real.log Real.pi)) : ℝ)) := by
      rw [← Complex.exp_add]; congr 1; push_cast; ring
    rw [e2]
    field_simp
  rw [hratio]
  have hl4 : Real.log 4 - Real.log Real.pi ≤ 3 / 10 := by
    rw [← Real.log_div (by norm_num) Real.pi_ne_zero]
    have := Real.log_le_sub_one_of_pos (x := 4 / Real.pi) (by positivity)
    have hp := Real.pi_gt_d2
    norm_num at hp
    have h54 : 4 / Real.pi ≤ 13 / 10 := by rw [div_le_iff₀ Real.pi_pos]; linarith
    linarith
  have hl4' : 0 ≤ Real.log 4 - Real.log Real.pi := by
    rw [sub_nonneg]; exact Real.log_le_log Real.pi_pos (by linarith [Real.pi_lt_d2])
  set P := cexp (I * ((δ * (Real.log 4 - Real.log Real.pi) : ℝ) : ℂ))
  have hPb : ‖P - 1‖ ≤ 3 / 10 * δ := by
    refine (Real.norm_exp_I_mul_ofReal_sub_one_le).trans ?_
    rw [Real.norm_eq_abs, abs_of_nonneg (by positivity)]
    nlinarith
  set G := Gamma (z + δ * I) / Gamma z
  have e : P * G - 1 = (P - 1) * G + (G - 1) := by ring
  have hGn : ‖G‖ ≤ 1 / (1 - S) := by
    have := norm_le_norm_add_norm_sub' G 1
    simp at this; linarith
  rw [e]
  calc ‖(P - 1) * G + (G - 1)‖ ≤ ‖P - 1‖ * ‖G‖ + ‖G - 1‖ :=
        (norm_add_le _ _).trans (by rw [norm_mul])
    _ ≤ 3 / 10 * δ * (1 / (1 - S)) + (1 / (1 - S) - 1) := by gcongr
    _ ≤ c := hc

end Fayed.Beta
