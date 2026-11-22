import Riemann.Deep.LeanNBRate.Contour3

/-!
# U1 contour argument, part 4: bounds for Psi, the strips [2,3] and [-3/2,-1]
-/

open Complex Filter Topology Set MeasureTheory intervalIntegral
open scoped Real Interval

noncomputable section

namespace Riemann.Deep.LeanNBRate

variable {H : ℂ → ℂ} {M β D : ℝ} {N : ℕ}

theorem integrable_of_cont_le_inv_sq {g : ℝ → ℂ} (hg : Continuous g) {C Y₀ : ℝ}
    (h : ∀ y, Y₀ ≤ |y| → ‖g y‖ ≤ C * ((1 + |y|) ^ 2)⁻¹) : Integrable g := by
  rw [← integrableOn_univ, show (univ : Set ℝ) = Icc (-Y₀) Y₀ ∪ (Icc (-Y₀) Y₀)ᶜ from (union_compl_self _).symm]
  refine IntegrableOn.union (hg.continuousOn.integrableOn_Icc) ?_
  refine ((integrable_inv_one_add_abs_sq.const_mul C).integrableOn).mono' hg.aestronglyMeasurable.restrict ?_
  filter_upwards [ae_restrict_mem (measurableSet_Icc.compl)] with y hy
  apply h y
  simp only [mem_compl_iff, mem_Icc, not_and_or, not_le] at hy
  rcases hy with hy | hy
  · linarith [neg_le_abs y]
  · linarith [le_abs_self y]

theorem continuousAt_Psi (hH : HBound H M β D N) (L : ℝ) {s z : ℂ} (hz1 : z ≠ 1)
    (hz0 : riemannZeta z ≠ 0) (hzs : z ≠ s) : ContinuousAt (Psi H L s) z := by
  unfold Psi kap
  exact ((continuousAt_Fq hH hz1 hz0).mul (by fun_prop)).div (continuousAt_const.sub continuousAt_id)
    (sub_ne_zero.mpr (Ne.symm hzs))

theorem differentiableAt_Psi (hH : HBound H M β D N) (L : ℝ) {s z : ℂ} (hz1 : z ≠ 1)
    (hz0 : riemannZeta z ≠ 0) (hzs : z ≠ s) : DifferentiableAt ℂ (Psi H L s) z := by
  unfold Psi kap
  exact ((differentiableAt_Fq hH hz1 hz0).mul (by fun_prop)).div
    ((differentiableAt_const _).sub differentiableAt_id) (sub_ne_zero.mpr (Ne.symm hzs))

theorem norm_Psi_eq (H : ℂ → ℂ) (L : ℝ) (s z : ℂ) :
    ‖Psi H L s z‖ = ‖H z‖ * ‖(riemannZeta z)⁻¹‖ * Real.exp (-(L * (s - z).re)) / ‖s - z‖ := by
  rw [Psi, Fq, norm_div, norm_mul, norm_div, norm_kap, norm_inv]
  ring

theorem norm_Psi_right (hH : HBound H M β D N) {L : ℝ} (hL : 0 ≤ L) {s z : ℂ} (hs : 3 < s.re)
    (hz1 : 2 ≤ z.re) (hz2 : z.re ≤ 3) :
    ‖Psi H L s z‖ ≤ 3 * M * (D / min D 1) ^ 2 / (s.re - 3) * ((1 + |z.im|) ^ 2)⁻¹ := by
  have hd : s.re - 3 ≤ ‖s - z‖ := by
    have := Complex.abs_re_le_norm (s - z); simp at this
    have h' : s.re - z.re ≤ |s.re - z.re| := le_abs_self _
    linarith
  have hk : ‖kap L (s - z)‖ ≤ 1 := norm_kap_le_one hL (by simp; linarith)
  rw [Psi, norm_div, norm_mul]
  have hF := norm_Fq_right_le hH hz1
  have hM := hH.M_nonneg
  calc ‖Fq H z‖ * ‖kap L (s - z)‖ / ‖s - z‖ ≤ (3 * M * (D / min D 1) ^ 2 * ((1 + |z.im|) ^ 2)⁻¹) * 1 /
        (s.re - 3) := by gcongr
    _ = _ := by ring

theorem continuousAt_Psi_right (hH : HBound H M β D N) (L : ℝ) {s z : ℂ} (hs : 3 < s.re)
    (hz1 : 2 ≤ z.re) (hz2 : z.re ≤ 3) : ContinuousAt (Psi H L s) z :=
  continuousAt_Psi hH L (ne_one_of_two_le hz1) (riemannZeta_ne_zero_of_one_lt_re (by linarith))
    (fun h => by rw [h] at hz2; linarith)

theorem integrable_Psi_line_right (hH : HBound H M β D N) {L : ℝ} (hL : 0 ≤ L) {s : ℂ} (hs : 3 < s.re)
    {σ : ℝ} (h1 : 2 ≤ σ) (h2 : σ ≤ 3) : Integrable fun y : ℝ => Psi H L s (σ + y * I) := by
  have hc : Continuous fun y : ℝ => Psi H L s (σ + y * I) := by
    rw [continuous_iff_continuousAt]
    intro y
    exact ContinuousAt.comp (f := fun y : ℝ => (σ : ℂ) + y * I) (x := y)
      (continuousAt_Psi_right hH L hs (by simp; exact h1) (by simp; exact h2))
      (by fun_prop : Continuous fun y : ℝ => (σ : ℂ) + y * I).continuousAt
  refine integrable_of_le_inv_sq (C := 3 * M * (D / min D 1) ^ 2 / (s.re - 3))
    hc.aestronglyMeasurable (fun y => ?_)
  have := norm_Psi_right hH hL hs (z := (σ : ℂ) + y * I) (by simp; exact h1) (by simp; exact h2)
  simpa using this

/-- Step C: the strip [2,3] (no poles). -/
theorem VI_Psi_three_two (hH : HBound H M β D N) {L : ℝ} (hL : 0 ≤ L) {s : ℂ} (hs : 3 < s.re) :
    I * VI (Psi H L s) 3 - I * VI (Psi H L s) 2 = 0 := by
  have h := strip_residues (f := Psi H L s) (a := 2) (b := 3) (by norm_num) ∅ (fun _ => 0)
    (by simp) (fun x h1 h2 _ => continuousAt_Psi_right hH L hs h1 h2)
    (fun x h1 h2 _ => differentiableAt_Psi hH L (ne_one_of_two_le h1.le)
      (riemannZeta_ne_zero_of_one_lt_re (by linarith)) (fun h => by rw [h] at h2; linarith))
    (by simp) (by simpa using integrable_Psi_line_right hH hL hs (σ := 2) (by norm_num) (by norm_num))
    (by simpa using integrable_Psi_line_right hH hL hs (σ := 3) (by norm_num) (by norm_num))
    (decay_of_le_inv_sq (C := 3 * M * (D / min D 1) ^ 2 / (s.re - 3)) fun z h1 h2 _ =>
      norm_Psi_right hH hL hs h1 h2)
  simpa using h

/-- Generic left-region bound. -/
theorem norm_Psi_left (hH : HBound H M β D N) {L : ℝ} (hL : 0 ≤ L) {s z : ℂ} (hs : 3 < s.re)
    (hz : z.re ≤ -1) (hzn : ∀ n : ℕ, z ≠ -n) (m : ℕ) (hm1 : 1 / 2 ≤ z.re + m) (hm2 : z.re + m ≤ 3 / 2)
    {P₀ : ℝ} (hP₀ : 0 < P₀) (hP : P₀ ≤ ∏ j ∈ Finset.range m, ‖z + j‖) :
    ‖Psi H L s z‖ ≤ M * Real.exp (β * (1 / 2 - z.re)) * (D / (D + |z.im|)) ^ N *
      (6 * Real.sqrt (2 * π) * (2 * π) ^ (-z.re) * (‖z‖ + m + 1) / P₀) * Real.exp (-(L * (s.re - z.re))) := by
  rw [norm_Psi_eq]
  have hHb := hH.bound z
  have hmax : max 0 (1 / 2 - z.re) = 1 / 2 - z.re := max_eq_right (by linarith)
  rw [hmax] at hHb
  have hw := weight_mono hH.D_pos N (z := z)
  have hzi := norm_inv_zeta_left hz hzn m hm1 hm2
  have hPpos : 0 < ∏ j ∈ Finset.range m, ‖z + j‖ := lt_of_lt_of_le hP₀ hP
  have hzi' : ‖(riemannZeta z)⁻¹‖ ≤ 6 * Real.sqrt (2 * π) * (2 * π) ^ (-z.re) * (‖z‖ + m + 1) / P₀ := by
    rw [le_div_iff₀ hP₀]
    calc ‖(riemannZeta z)⁻¹‖ * P₀ ≤ ‖(riemannZeta z)⁻¹‖ * ∏ j ∈ Finset.range m, ‖z + j‖ :=
          mul_le_mul_of_nonneg_left hP (norm_nonneg _)
      _ ≤ _ := hzi
  have hsz : 1 ≤ ‖s - z‖ := by
    have := Complex.abs_re_le_norm (s - z); simp at this
    have h' : s.re - z.re ≤ |s.re - z.re| := le_abs_self _
    linarith
  have hre : (s - z).re = s.re - z.re := by simp
  rw [hre]
  have hM := hH.M_nonneg
  have hD := hH.D_pos
  calc ‖H z‖ * ‖(riemannZeta z)⁻¹‖ * Real.exp (-(L * (s.re - z.re))) / ‖s - z‖
      ≤ ‖H z‖ * ‖(riemannZeta z)⁻¹‖ * Real.exp (-(L * (s.re - z.re))) :=
        div_le_self (by positivity) hsz
    _ ≤ (M * Real.exp (β * (1 / 2 - z.re)) * (D / (D + |z.im|)) ^ N) *
        (6 * Real.sqrt (2 * π) * (2 * π) ^ (-z.re) * (‖z‖ + m + 1) / P₀) * Real.exp (-(L * (s.re - z.re))) := by
        gcongr
        exact hHb.trans (mul_le_mul_of_nonneg_left hw (by positivity))

theorem floor_bounds {x : ℝ} (hx : x ≤ -1) : 1 / 2 ≤ x + ⌊3 / 2 - x⌋₊ ∧ x + ⌊3 / 2 - x⌋₊ ≤ 3 / 2 := by
  have h0 : 0 ≤ 3 / 2 - x := by linarith
  have h1 := Nat.floor_le h0
  have h2 := Nat.lt_floor_add_one (3 / 2 - x)
  constructor <;> linarith

theorem one_le_prod_of_im {z : ℂ} (hz : 1 ≤ |z.im|) (m : ℕ) : 1 ≤ ∏ j ∈ Finset.range m, ‖z + j‖ := by
  induction m with
  | zero => simp
  | succ m ih =>
    rw [Finset.prod_range_succ]
    have : 1 ≤ ‖z + (m : ℂ)‖ := by
      have := Complex.abs_im_le_norm (z + m); simp at this; linarith
    nlinarith

theorem not_nat_of_im {z : ℂ} (hz : 1 ≤ |z.im|) : ∀ n : ℕ, z ≠ -n := fun n h => by
  rw [h] at hz; simp at hz; linarith

/-- Left strips, |Im z| >= 1: |Psi z| <= C/(1+|Im z|)^2 on a < Re z < b <= -1. -/
theorem norm_Psi_left_strip (hH : HBound H M β D N) {L : ℝ} (hL : 0 ≤ L) {s : ℂ} (hs : 3 < s.re)
    {a b : ℝ} (hab : a ≤ b) (hb : b ≤ -1) {z : ℂ} (hz1 : a ≤ z.re) (hz2 : z.re ≤ b) (hzi : 1 ≤ |z.im|) :
    ‖Psi H L s z‖ ≤ M * Real.exp (β * (1 / 2 - a)) * (6 * Real.sqrt (2 * π) * (2 * π) ^ (-a) *
      (|a| + (3 / 2 - a) + 1 + 1)) * (D / min D 1) ^ N * ((1 + |z.im|) ^ 2)⁻¹ := by
  obtain ⟨hm1, hm2⟩ := floor_bounds (hz2.trans hb)
  set m := ⌊3 / 2 - z.re⌋₊ with hm
  have h := norm_Psi_left hH hL hs (hz2.trans hb) (not_nat_of_im hzi) m hm1 hm2 one_pos
    (one_le_prod_of_im hzi m)
  have hM := hH.M_nonneg
  have hD := hH.D_pos
  have hmle : (m : ℝ) ≤ 3 / 2 - a := by
    have := Nat.floor_le (show 0 ≤ 3 / 2 - z.re by linarith); linarith
  have hnz : ‖z‖ ≤ |a| + |z.im| := by
    calc ‖z‖ ≤ |z.re| + |z.im| := Complex.norm_le_abs_re_add_abs_im z
      _ ≤ |a| + |z.im| := by
          have : |z.re| ≤ |a| := by rw [abs_of_neg (by linarith), abs_of_neg (by linarith)]; linarith
          linarith
  have e1 : Real.exp (β * (1 / 2 - z.re)) ≤ Real.exp (β * (1 / 2 - a)) :=
    Real.exp_le_exp.mpr (mul_le_mul_of_nonneg_left (by linarith) hH.β_nonneg)
  have e2 : (2 * π) ^ (-z.re) ≤ (2 * π) ^ (-a) :=
    Real.rpow_le_rpow_of_exponent_le (by linarith [Real.two_le_pi]) (by linarith)
  have e3 : Real.exp (-(L * (s.re - z.re))) ≤ 1 := by
    rw [Real.exp_le_one_iff]; nlinarith
  have e4 : ‖z‖ + m + 1 ≤ (|a| + (3 / 2 - a) + 1 + 1) * (1 + |z.im|) := by
    have : 0 ≤ |a| + (3 / 2 - a) + 1 := by have := abs_nonneg a; linarith
    nlinarith [abs_nonneg z.im]
  have hwl := weight_le hD (by linarith [hH.N_ge] : 4 ≤ N) z.im
  have hy1 : 1 + |z.im| ≤ (1 + |z.im|) ^ 2 := by nlinarith [abs_nonneg z.im]
  have hK : 0 ≤ |a| + (3 / 2 - a) + 1 + 1 := by have := abs_nonneg a; linarith
  calc ‖Psi H L s z‖ ≤ M * Real.exp (β * (1 / 2 - z.re)) * (D / (D + |z.im|)) ^ N *
        (6 * Real.sqrt (2 * π) * (2 * π) ^ (-z.re) * (‖z‖ + m + 1) / 1) * Real.exp (-(L * (s.re - z.re))) := h
    _ ≤ M * Real.exp (β * (1 / 2 - a)) * (D / (D + |z.im|)) ^ N *
        (6 * Real.sqrt (2 * π) * (2 * π) ^ (-a) * ((|a| + (3 / 2 - a) + 1 + 1) * (1 + |z.im|) ^ 2) / 1) * 1 := by
        gcongr
        exact e4.trans (mul_le_mul_of_nonneg_left hy1 hK)
    _ = M * Real.exp (β * (1 / 2 - a)) * (6 * Real.sqrt (2 * π) * (2 * π) ^ (-a) *
        (|a| + (3 / 2 - a) + 1 + 1)) * ((1 + |z.im|) ^ 2 * (D / (D + |z.im|)) ^ N) := by ring
    _ ≤ M * Real.exp (β * (1 / 2 - a)) * (6 * Real.sqrt (2 * π) * (2 * π) ^ (-a) *
        (|a| + (3 / 2 - a) + 1 + 1)) * ((D / min D 1) ^ N * ((1 + |z.im|) ^ 2)⁻¹) := by gcongr
    _ = _ := by ring

end Riemann.Deep.LeanNBRate

