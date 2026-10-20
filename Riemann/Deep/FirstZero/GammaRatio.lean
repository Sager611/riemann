import Mathlib.Analysis.SpecialFunctions.Gamma.Beta
import Mathlib.Analysis.SpecialFunctions.Trigonometric.Bounds

/-!
# Phase stability of Γ on a short vertical segment

‖Γ(z + iδ)/Γ(z) - 1‖ ≤ 1/4 for Re z = 1/4, 7 ≤ Im z ≤ 71/10, 0 ≤ δ ≤ 1/200.
Proof: Euler's limit Γ = lim GammaSeq (Mathlib). The ratio of the two GammaSeq terms is
n^{iδ} Π_j (z+j)/(z+iδ+j); telescoping n^{iδ} = Π_j exp(iδ log((j+2)/(j+1))) and pairing
factors gives terms 1 + O(δ |z| / (j+1)²), whose product stays within exp(0.18) - 1 of 1.
-/

open Complex Filter Topology

namespace Fayed.FirstZero

lemma norm_mul_one_add_sub_one_le {x e : ℂ} {S : ℝ} (h : ‖x - 1‖ ≤ Real.exp S - 1) :
    ‖x * (1 + e) - 1‖ ≤ Real.exp (S + ‖e‖) - 1 := by
  have h1 : x * (1 + e) - 1 = (x - 1) * (1 + e) + e := by ring
  rw [h1]
  have hE := Real.exp_pos S
  have he := norm_nonneg e
  have hexp := Real.add_one_le_exp ‖e‖
  have h1e : ‖(1 : ℂ) + e‖ ≤ 1 + ‖e‖ := (norm_add_le _ _).trans (by simp)
  calc ‖(x - 1) * (1 + e) + e‖ ≤ ‖x - 1‖ * (1 + ‖e‖) + ‖e‖ := by
        refine (norm_add_le _ _).trans ?_
        rw [norm_mul]
        gcongr
    _ ≤ (Real.exp S - 1) * (1 + ‖e‖) + ‖e‖ := by gcongr
    _ = Real.exp S * (1 + ‖e‖) - 1 := by ring
    _ ≤ Real.exp S * Real.exp ‖e‖ - 1 := by gcongr; linarith
    _ = Real.exp (S + ‖e‖) - 1 := by rw [Real.exp_add]

lemma norm_prod_one_add_sub_one_le (e : ℕ → ℂ) (m : ℕ) :
    ‖∏ j ∈ Finset.range m, (1 + e j) - 1‖ ≤ Real.exp (∑ j ∈ Finset.range m, ‖e j‖) - 1 := by
  induction m with
  | zero => simp
  | succ m ih =>
    rw [Finset.prod_range_succ, Finset.sum_range_succ]
    exact norm_mul_one_add_sub_one_le ih

lemma sum_range_inv_sq_le_two (m : ℕ) :
    ∑ j ∈ Finset.range m, 1 / ((j : ℝ) + 1) ^ 2 ≤ 2 := by
  have key : ∀ m : ℕ, ∑ j ∈ Finset.range (m + 1), 1 / ((j : ℝ) + 1) ^ 2 ≤ 2 - 1 / ((m : ℝ) + 1) := by
    intro m
    induction m with
    | zero => norm_num
    | succ m ih =>
      rw [Finset.sum_range_succ]
      have hm : (0 : ℝ) ≤ m := Nat.cast_nonneg m
      have e1 : 1 / ((m : ℝ) + 1) - 1 / ((m : ℝ) + 2) = 1 / (((m : ℝ) + 1) * ((m : ℝ) + 2)) := by
        field_simp; ring
      have e2 : 1 / ((m : ℝ) + 2) ^ 2 ≤ 1 / (((m : ℝ) + 1) * ((m : ℝ) + 2)) :=
        one_div_le_one_div_of_le (by positivity) (by nlinarith)
      push_cast
      have e3 : ((m : ℝ) + 1 + 1) = (m : ℝ) + 2 := by ring
      rw [e3]
      linarith
  cases m with
  | zero => simp
  | succ m =>
    have : (0 : ℝ) < 1 / ((m : ℝ) + 1) := by positivity
    linarith [key m]

lemma log_succ_ratio_bounds (j : ℕ) :
    1 / ((j : ℝ) + 2) ≤ Real.log ((j : ℝ) + 2) - Real.log ((j : ℝ) + 1) ∧
      Real.log ((j : ℝ) + 2) - Real.log ((j : ℝ) + 1) ≤ 1 / ((j : ℝ) + 1) := by
  have h1 : (0 : ℝ) < (j : ℝ) + 1 := by positivity
  have h2 : (0 : ℝ) < (j : ℝ) + 2 := by positivity
  rw [← Real.log_div h2.ne' h1.ne']
  constructor
  · have := Real.one_sub_inv_le_log_of_pos (x := ((j : ℝ) + 2) / ((j : ℝ) + 1)) (by positivity)
    rw [inv_div] at this
    have e : 1 - ((j : ℝ) + 1) / ((j : ℝ) + 2) = 1 / ((j : ℝ) + 2) := by field_simp; ring
    linarith
  · have := Real.log_le_sub_one_of_pos (x := ((j : ℝ) + 2) / ((j : ℝ) + 1)) (by positivity)
    have e : ((j : ℝ) + 2) / ((j : ℝ) + 1) - 1 = 1 / ((j : ℝ) + 1) := by field_simp; ring
    linarith

lemma norm_shift_ge (z : ℂ) (hre : z.re = 1 / 4) (him : 7 ≤ z.im) (j : ℕ) :
    ((j : ℝ) + 1) / 2 ≤ ‖z + j‖ := by
  have hsq : ‖z + j‖ ^ 2 = (1 / 4 + (j : ℝ)) ^ 2 + z.im ^ 2 := by
    rw [Complex.sq_norm, Complex.normSq_apply]
    simp only [add_re, natCast_re, add_im, natCast_im, add_zero, hre]
    ring
  have hj : (0 : ℝ) ≤ j := Nat.cast_nonneg j
  nlinarith [norm_nonneg (z + j)]

lemma norm_shift_ge_im (z : ℂ) (him : 7 ≤ z.im) (j : ℕ) : 7 ≤ ‖z + j‖ := by
  have h := abs_im_le_norm (z + j)
  simp only [add_im, natCast_im, add_zero] at h
  rw [abs_of_nonneg (by linarith)] at h
  linarith

/-- One paired factor of the Gamma-ratio product. -/
lemma pair_bound (z : ℂ) (hre : z.re = 1 / 4) (him : 7 ≤ z.im) {δ : ℝ} (hδ0 : 0 ≤ δ)
    (hδ1 : δ ≤ 1) (j : ℕ) :
    ‖(1 - δ * I / (z + j)) *
        cexp (δ * I * ((Real.log ((j : ℝ) + 2) - Real.log ((j : ℝ) + 1) : ℝ) : ℂ)) - 1‖
      ≤ (3 * δ ^ 2 + δ + 2 * δ * ‖z - 1‖) / ((j : ℝ) + 1) ^ 2 := by
  set w : ℂ := z + j with hw
  set L : ℝ := Real.log ((j : ℝ) + 2) - Real.log ((j : ℝ) + 1) with hL
  obtain ⟨hL0, hL1⟩ := log_succ_ratio_bounds j
  rw [← hL] at hL0 hL1
  have hj1 : (0 : ℝ) < (j : ℝ) + 1 := by positivity
  have hwn : ((j : ℝ) + 1) / 2 ≤ ‖w‖ := norm_shift_ge z hre him j
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

lemma GammaSeq_ratio (z z2 : ℂ) (hz : 0 < z.re) (δ : ℝ) (hz2 : z2 = z + δ * I) (m : ℕ) :
    GammaSeq z2 (m + 1) / GammaSeq z (m + 1) =
      (∏ j ∈ Finset.range m, ((1 - δ * I / (z2 + j)) *
        cexp (δ * I * ((Real.log ((j : ℝ) + 2) - Real.log ((j : ℝ) + 1) : ℝ) : ℂ)))) *
      (1 - δ * I / (z2 + m)) * (1 - δ * I / (z2 + ((m + 1 : ℕ) : ℂ))) := by
  have hne : ∀ (w : ℂ), 0 < w.re → ∀ j : ℕ, w + j ≠ 0 := fun w hw j h => by
    have := congrArg Complex.re h
    simp only [add_re, natCast_re, zero_re] at this
    linarith [(Nat.cast_nonneg j : (0 : ℝ) ≤ j)]
  have hz2re : 0 < z2.re := by rw [hz2]; simpa using hz
  have hP1 : ∏ j ∈ Finset.range (m + 1 + 1), (z + j) ≠ 0 :=
    Finset.prod_ne_zero_iff.mpr fun j _ => hne z hz j
  have hP2 : ∏ j ∈ Finset.range (m + 1 + 1), (z2 + j) ≠ 0 :=
    Finset.prod_ne_zero_iff.mpr fun j _ => hne z2 hz2re j
  have hn0 : ((m + 1 : ℕ) : ℂ) ≠ 0 := by exact_mod_cast Nat.succ_ne_zero m
  have hfac : (((m + 1).factorial : ℕ) : ℂ) ≠ 0 := by exact_mod_cast (Nat.factorial_pos _).ne'
  have hnz : ((m + 1 : ℕ) : ℂ) ^ z ≠ 0 := by
    rw [Ne, cpow_eq_zero_iff]; tauto
  have hpow : ((m + 1 : ℕ) : ℂ) ^ z2 = ((m + 1 : ℕ) : ℂ) ^ z *
      ∏ j ∈ Finset.range m,
        cexp (δ * I * ((Real.log ((j : ℝ) + 2) - Real.log ((j : ℝ) + 1) : ℝ) : ℂ)) := by
    have hlog : Complex.log ((m + 1 : ℕ) : ℂ) = ((Real.log ((m : ℝ) + 1) : ℝ) : ℂ) := by
      rw [← Complex.natCast_log]; push_cast; ring_nf
    have hsum : Real.log ((m : ℝ) + 1) =
        ∑ j ∈ Finset.range m, (Real.log ((j : ℝ) + 2) - Real.log ((j : ℝ) + 1)) := by
      have := Finset.sum_range_sub (fun j : ℕ => Real.log ((j : ℝ) + 1)) m
      simp only [Nat.cast_add, Nat.cast_one, Nat.cast_zero, zero_add, Real.log_one,
        sub_zero] at this
      rw [← this]
      refine Finset.sum_congr rfl fun j _ => ?_
      ring_nf
    have hpowI : ((m + 1 : ℕ) : ℂ) ^ ((δ : ℂ) * I) = ∏ j ∈ Finset.range m,
        cexp (δ * I * ((Real.log ((j : ℝ) + 2) - Real.log ((j : ℝ) + 1) : ℝ) : ℂ)) := by
      rw [cpow_def_of_ne_zero hn0, hlog, hsum, ← Complex.exp_sum]
      congr 1
      push_cast
      rw [Finset.sum_mul]
      refine Finset.sum_congr rfl fun j _ => ?_
      ring
    rw [hz2, cpow_add _ _ hn0, hpowI]
  have hfactor : ∀ j : ℕ, (z + j) / (z2 + j) = 1 - δ * I / (z2 + j) := by
    intro j
    rw [one_sub_div (hne z2 hz2re j), hz2]
    congr 1; ring
  have hprod : (∏ j ∈ Finset.range (m + 1 + 1), (z + j)) / ∏ j ∈ Finset.range (m + 1 + 1), (z2 + j)
      = ∏ j ∈ Finset.range (m + 1 + 1), (1 - δ * I / (z2 + j)) := by
    rw [← Finset.prod_div_distrib]
    exact Finset.prod_congr rfl fun j _ => hfactor j
  unfold GammaSeq
  rw [hpow]
  have e1 : ((m + 1 : ℕ) : ℂ) ^ z *
        (∏ j ∈ Finset.range m,
          cexp (δ * I * ((Real.log ((j : ℝ) + 2) - Real.log ((j : ℝ) + 1) : ℝ) : ℂ))) *
        (((m + 1).factorial : ℕ) : ℂ) / (∏ j ∈ Finset.range (m + 1 + 1), (z2 + j)) /
      (((m + 1 : ℕ) : ℂ) ^ z * (((m + 1).factorial : ℕ) : ℂ) /
        ∏ j ∈ Finset.range (m + 1 + 1), (z + j)) =
      ((∏ j ∈ Finset.range (m + 1 + 1), (z + j)) / ∏ j ∈ Finset.range (m + 1 + 1), (z2 + j)) *
        ∏ j ∈ Finset.range m,
          cexp (δ * I * ((Real.log ((j : ℝ) + 2) - Real.log ((j : ℝ) + 1) : ℝ) : ℂ)) := by
    field_simp
  rw [e1, hprod, Finset.prod_range_succ, Finset.prod_range_succ, mul_comm,
    ← mul_assoc, ← mul_assoc, ← Finset.prod_mul_distrib]
  congr 2
  exact Finset.prod_congr rfl fun j _ => by ring

/-- **Γ moves little along a short vertical segment at height ≈ 7.** -/
theorem Gamma_ratio_le (z : ℂ) (hre : z.re = 1 / 4) (him : 7 ≤ z.im) (him' : z.im ≤ 71 / 10)
    {δ : ℝ} (hδ0 : 0 ≤ δ) (hδ : δ ≤ 1 / 200) :
    ‖Gamma (z + δ * I) / Gamma z - 1‖ ≤ 1 / 4 := by
  obtain ⟨z2, hz2⟩ : ∃ z2, z2 = z + δ * I := ⟨_, rfl⟩
  rw [← hz2]
  have hre2 : z2.re = 1 / 4 := by rw [hz2]; simp [hre]
  have him2 : 7 ≤ z2.im := by rw [hz2]; simp; linarith
  have hz1 : ‖z2 - 1‖ ≤ 8 := by
    refine (norm_le_abs_re_add_abs_im _).trans ?_
    have h1 : (z2 - 1).re = -3 / 4 := by rw [sub_re, hre2]; norm_num
    have h2 : (z2 - 1).im = z.im + δ := by rw [hz2]; simp
    rw [h1, h2, abs_of_nonneg (show 0 ≤ z.im + δ by linarith)]
    norm_num; linarith
  have hzre : 0 < z.re := by rw [hre]; norm_num
  have hG : Gamma z ≠ 0 := Gamma_ne_zero_of_re_pos hzre
  have hlim := (((GammaSeq_tendsto_Gamma z2).div (GammaSeq_tendsto_Gamma z) hG).sub_const 1).norm
  refine le_of_tendsto hlim (eventually_atTop.2 ⟨1, fun n hn => ?_⟩)
  obtain ⟨m, rfl⟩ : ∃ m, n = m + 1 := ⟨n - 1, by omega⟩
  rw [Pi.div_apply, GammaSeq_ratio z z2 hzre δ hz2 m]
  set K := 3 * δ ^ 2 + δ + 2 * δ * ‖z2 - 1‖ with hK
  set e : ℕ → ℂ := fun j => (1 - δ * I / (z2 + j)) *
    cexp (δ * I * ((Real.log ((j : ℝ) + 2) - Real.log ((j : ℝ) + 1) : ℝ) : ℂ)) - 1 with he_def
  have he : ∀ j, ‖e j‖ ≤ K / ((j : ℝ) + 1) ^ 2 :=
    fun j => pair_bound z2 hre2 him2 hδ0 (by linarith) j
  have hK0 : 0 ≤ K := by positivity
  have hS : ∑ j ∈ Finset.range m, ‖e j‖ ≤ 2 * K := by
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
  have ha : ∀ k : ℕ, ‖-(δ * I / (z2 + k))‖ ≤ δ / 7 := by
    intro k
    have h7 := norm_shift_ge_im z2 him2 k
    rw [norm_neg, norm_div, norm_mul, Complex.norm_real, Real.norm_eq_abs, Complex.norm_I, mul_one,
      abs_of_nonneg hδ0]
    exact div_le_div_of_nonneg_left hδ0 (by norm_num) h7
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
  have hKb : K ≤ 3 * δ ^ 2 + δ + 2 * δ * 8 := by
    rw [hK]; nlinarith [mul_le_mul_of_nonneg_left hz1 hδ0]
  have hsum : ∑ j ∈ Finset.range m, ‖e j‖ + ‖-(δ * I / (z2 + m))‖
      + ‖-(δ * I / (z2 + ((m + 1 : ℕ) : ℂ)))‖ ≤ 9 / 50 := by
    nlinarith
  have hexp : Real.exp (9 / 50) < 1 / (1 - 9 / 50) :=
    Real.exp_bound_div_one_sub_of_interval' (by norm_num) (by norm_num)
  have hmono := Real.exp_le_exp.mpr hsum
  norm_num at hexp
  linarith

end Fayed.FirstZero
