import Riemann.Deep.Hardy.Moments

/-!
# Hardy's theorem, part 5: behaviour of the kernel derivatives at the cusp

* iteratedDeriv_Hk_eq: Hk^{(m)} = Rk^{(m)} - Ck^{(m)} on the strip;
* iteratedDeriv_Ck_even: Ck^{(2n)} = 16^{-n} Ck;
* Rk_cusp_bound: near u = -iπ/2 the theta part is exp-small: ‖Rk u‖ ≤ K/δ · e^{-1/(72δ)}
  on the circle of radius δ/2 about (-π/2 + δ) i (theta transformation at the cusp τ = -1);
* tendsto_iteratedDeriv_Rk: Rk^{(m)}(-2αi) → 0 as α → π/4⁻ (Cauchy estimates).
-/

open Complex Real MeasureTheory Set Filter Topology
open scoped ComplexConjugate Nat

namespace Fayed.Hardy

lemma contDiffAt_Rk {z : ℂ} (hz : z ∈ strip) (n : ℕ) : ContDiffAt ℂ n Rk z :=
  (differentiableOn_Rk.analyticOnNhd isOpen_strip z hz).contDiffAt

lemma contDiffAt_Ck (z : ℂ) (n : ℕ) : ContDiffAt ℂ n Ck z :=
  ((analyticOnNhd_univ_iff_differentiable.mpr differentiable_Ck) z trivial).contDiffAt

lemma iteratedDeriv_Hk_eq (m : ℕ) {z : ℂ} (hz : z ∈ strip) :
    iteratedDeriv m Hk z = iteratedDeriv m Rk z - iteratedDeriv m Ck z := by
  have : Hk = Rk - Ck := rfl
  rw [this, iteratedDeriv_sub (contDiffAt_Rk hz m) (contDiffAt_Ck z m)]

lemma hasDerivAt_Ck (u : ℂ) : HasDerivAt Ck ((1 / 2) * Complex.sinh (u / 4)) u := by
  have h1 : HasDerivAt (fun u : ℂ => u / 4) (1 / 4) u := (hasDerivAt_id u).div_const 4
  have h2 := (h1.ccosh).const_mul 2
  show HasDerivAt (fun u : ℂ => 2 * Complex.cosh (u / 4)) _ u
  convert h2 using 1; ring

lemma hasDerivAt_sinh4 (u : ℂ) :
    HasDerivAt (fun u : ℂ => Complex.sinh (u / 4)) ((1 / 4) * Complex.cosh (u / 4)) u := by
  have h1 : HasDerivAt (fun u : ℂ => u / 4) (1 / 4) u := (hasDerivAt_id u).div_const 4
  have h2 := h1.csinh
  convert h2 using 1; ring

lemma iteratedDeriv_Ck_even (n : ℕ) :
    iteratedDeriv (2 * n) Ck = fun u => (1 / 16 : ℂ) ^ n * Ck u := by
  induction n with
  | zero => ext u; simp
  | succ n ih =>
    have h1 : iteratedDeriv (2 * n + 1) Ck =
        fun u => (1 / 16 : ℂ) ^ n * ((1 / 2) * Complex.sinh (u / 4)) := by
      rw [iteratedDeriv_succ, ih]; ext u
      exact ((hasDerivAt_Ck u).const_mul _).deriv
    rw [show 2 * (n + 1) = 2 * n + 1 + 1 by ring, iteratedDeriv_succ, h1]
    ext u
    have := (((hasDerivAt_sinh4 u).const_mul (1 / 2 : ℂ)).const_mul ((1 / 16 : ℂ) ^ n)).deriv
    rw [this, Ck, pow_succ]; ring

lemma Ck_neg_two_mul_I (α : ℝ) : Ck (-2 * α * I) = ((2 * Real.cos (α / 2) : ℝ) : ℂ) := by
  rw [Ck, show -2 * (α : ℂ) * I / 4 = ((-(α / 2) : ℝ) : ℂ) * I by push_cast; ring,
    Complex.cosh_mul_I, ← ofReal_cos, Real.cos_neg]
  push_cast; ring

/-- The shift τ ↦ τ + 1 moves theta to the characteristic 1/2. -/
lemma jacobiTheta₂_zero_eq_half (τ : ℂ) : jacobiTheta₂ 0 τ = jacobiTheta₂ (1 / 2) (τ + 1) := by
  unfold jacobiTheta₂
  congr 1; ext n
  unfold jacobiTheta₂_term
  obtain ⟨k, hk⟩ := Int.even_mul_succ_self n
  have : 2 * π * I * n * (1 / 2) + π * I * n ^ 2 * (τ + 1) =
      (2 * π * I * n * 0 + π * I * n ^ 2 * τ) + k * (2 * π * I) := by
    have hk' : ((n : ℂ) * (n + 1)) = (k : ℂ) + k := by exact_mod_cast hk
    linear_combination (π * I) * hk'
  rw [this, Complex.exp_add (_ + _) (k * (2 * π * I)), Complex.exp_int_mul_two_pi_mul_I, mul_one]

/-- A summable majorant for the theta terms. -/
noncomputable def Sθ : ℝ := ∑' n : ℤ, rexp (-π * (1 * (n : ℝ) ^ 2 - 2 * (1 / 2) * |(n : ℝ)|))

lemma summable_Sθ : Summable (fun n : ℤ => rexp (-π * (1 * (n : ℝ) ^ 2 - 2 * (1 / 2) * |(n : ℝ)|))) := by
  simpa using summable_pow_mul_jacobiTheta₂_term_bound (1 / 2) one_pos 0

lemma int_abs_le_sq (n : ℤ) : |(n : ℝ)| ≤ (n : ℝ) ^ 2 := by
  have : |n| ≤ n ^ 2 := by
    rcases eq_or_ne n 0 with h | h
    · simp [h]
    · rw [← sq_abs]; nlinarith [Int.one_le_abs h]
  have h := (Int.cast_le (R := ℝ)).mpr this
  push_cast at h; exact h

/-- Bound for theta with characteristic 1/(2v) at -1/v when Im(-1/v) ≥ 1. -/
lemma norm_theta_aux {Y : ℝ} {z τ : ℂ} (hτ : τ.im = Y) (hz : z.im = -Y / 2) (hY : 1 ≤ Y) :
    ‖jacobiTheta₂ z τ‖ ≤ Sθ := by
  have hb : ∀ n : ℤ, ‖jacobiTheta₂_term n z τ‖ ≤
      rexp (-π * (1 * (n : ℝ) ^ 2 - 2 * (1 / 2) * |(n : ℝ)|)) := by
    intro n
    rw [norm_jacobiTheta₂_term, hτ, hz]
    apply Real.exp_le_exp.mpr
    have h1 := int_abs_le_sq n
    have h2 : (n : ℝ) ≤ |(n : ℝ)| := le_abs_self _
    have h3 : 0 ≤ (n : ℝ) ^ 2 - |(n : ℝ)| := by linarith
    have h4 : (n : ℝ) ^ 2 - |(n : ℝ)| ≤ Y * ((n : ℝ) ^ 2 - n) := by nlinarith
    have h5 := mul_le_mul_of_nonneg_left h4 Real.pi_pos.le
    nlinarith [h5]
  have hs : Summable (fun n : ℤ => ‖jacobiTheta₂_term n z τ‖) :=
    Summable.of_nonneg_of_le (fun _ => norm_nonneg _) hb summable_Sθ
  calc ‖jacobiTheta₂ z τ‖ ≤ ∑' n : ℤ, ‖jacobiTheta₂_term n z τ‖ := norm_tsum_le_tsum_norm hs
    _ ≤ Sθ := hs.tsum_le_tsum hb summable_Sθ

lemma Sθ_nonneg : 0 ≤ Sθ := tsum_nonneg (fun _ => (Real.exp_pos _).le)

/-- The theta part is exponentially small near the cusp. -/
lemma Rk_cusp_bound {δ : ℝ} (hδ0 : 0 < δ) (hδ : δ ≤ 1 / 100) {u : ℂ}
    (hu : ‖u - ((-π / 2 + δ : ℝ) : ℂ) * I‖ ≤ δ / 2) :
    ‖Rk u‖ ≤ (2 * π * rexp 1 * Sθ / δ) * rexp (-(1 / (72 * δ))) := by
  have hπ3 := Real.pi_gt_three
  have hπ4 := Real.pi_le_four
  set ζ := u + ((π / 2 : ℝ) : ℂ) * I with hζ
  have hζ1 : ‖ζ - (δ : ℂ) * I‖ ≤ δ / 2 := by
    convert hu using 2; rw [hζ]; push_cast; ring
  have hre : |ζ.re| ≤ δ / 2 := by
    have := (abs_re_le_norm (ζ - (δ : ℂ) * I)).trans hζ1
    simpa using this
  have him : |ζ.im - δ| ≤ δ / 2 := by
    have := (abs_im_le_norm (ζ - (δ : ℂ) * I)).trans hζ1
    simpa using this
  have hb1 : δ / 2 ≤ ζ.im := by linarith [abs_le.mp him]
  have hb2 : ζ.im ≤ 3 * δ / 2 := by linarith [abs_le.mp him]
  have ha2 : ζ.re ≤ δ / 2 := (abs_le.mp hre).2
  have hζn : ‖ζ‖ ≤ 3 * δ / 2 := by
    have hn : ‖(δ : ℂ) * I‖ = δ := by simp [abs_of_pos hδ0]
    calc ‖ζ‖ = ‖(ζ - (δ : ℂ) * I) + (δ : ℂ) * I‖ := by ring_nf
      _ ≤ ‖ζ - (δ : ℂ) * I‖ + ‖(δ : ℂ) * I‖ := norm_add_le _ _
      _ ≤ 3 * δ / 2 := by rw [hn]; linarith
  -- the theta argument
  have hI : cexp (((π / 2 : ℝ) : ℂ) * I) = I := by
    rw [Complex.exp_mul_I, ← ofReal_cos, ← ofReal_sin, Real.cos_pi_div_two, Real.sin_pi_div_two]
    simp
  have hexp : cexp (-u) = I * cexp (-ζ) := by
    rw [show -u = ((π / 2 : ℝ) : ℂ) * I + -ζ by rw [hζ]; ring, Complex.exp_add, hI]
  set v := 1 - cexp (-ζ) with hv
  have hθ : theta (cexp (-u)) = jacobiTheta₂ (1 / 2) v := by
    rw [theta, hexp, ← mul_assoc, I_mul_I, jacobiTheta₂_zero_eq_half, hv]; congr 1; ring
  have hvim : v.im = rexp (-ζ.re) * Real.sin ζ.im := by
    rw [hv, sub_im, one_im, exp_im, neg_re, neg_im, Real.sin_neg]; ring
  have hsin : δ / π ≤ Real.sin ζ.im := by
    have := Real.mul_le_sin (by linarith : 0 ≤ ζ.im) (by linarith : ζ.im ≤ π / 2)
    have h2 : δ / π ≤ 2 / π * ζ.im := by
      rw [div_le_iff₀ (by positivity)]; field_simp; linarith
    linarith
  have hea : 1 / 2 ≤ rexp (-ζ.re) := by linarith [Real.add_one_le_exp (-ζ.re)]
  have hvim_lb : δ / (2 * π) ≤ v.im := by
    rw [hvim]
    calc δ / (2 * π) = 1 / 2 * (δ / π) := by field_simp
      _ ≤ rexp (-ζ.re) * Real.sin ζ.im :=
        mul_le_mul hea hsin (by positivity) (Real.exp_pos _).le
  have hvpos : 0 < v.im := lt_of_lt_of_le (by positivity) hvim_lb
  have hvn_ub : ‖v‖ ≤ 3 * δ := by
    rw [hv, norm_sub_rev]
    have h1 : ‖-ζ‖ ≤ 1 := by rw [norm_neg]; linarith
    have := Complex.norm_exp_sub_one_le h1
    rw [norm_neg] at this; linarith
  have hvn_lb : δ / (2 * π) ≤ ‖v‖ := hvim_lb.trans ((le_abs_self _).trans (abs_im_le_norm v))
  have hvn_pos : 0 < ‖v‖ := lt_of_lt_of_le (by positivity) hvn_lb
  have hv0 : v ≠ 0 := norm_pos_iff.mp hvn_pos
  have hns : Complex.normSq v = ‖v‖ ^ 2 := Complex.normSq_eq_norm_sq v
  set Y := v.im / Complex.normSq v with hY
  have hY_lb : 1 / (18 * π * δ) ≤ Y := by
    rw [hY, hns, div_le_div_iff₀ (by positivity) (by positivity)]
    have : ‖v‖ ^ 2 ≤ 9 * δ ^ 2 := by nlinarith [norm_nonneg v]
    have h2 : δ / (2 * π) * (18 * π * δ) = 9 * δ ^ 2 := by field_simp; ring
    nlinarith [mul_le_mul_of_nonneg_right hvim_lb (by positivity : (0 : ℝ) ≤ 18 * π * δ)]
  have hY1 : 1 ≤ Y := by
    refine le_trans ?_ hY_lb
    rw [le_div_iff₀ (by positivity)]; nlinarith
  -- functional equation
  have hFE := jacobiTheta₂_functional_equation (1 / 2) v
  have hn1 : ‖1 / (-I * v) ^ (1 / 2 : ℂ)‖ ≤ 2 * π / δ := by
    rw [norm_div, norm_one, show (1 / 2 : ℂ) = ((1 / 2 : ℝ) : ℂ) by push_cast; ring,
      norm_cpow_real, norm_mul, norm_neg, norm_I, one_mul]
    have h1 : ‖v‖ ≤ ‖v‖ ^ (1 / 2 : ℝ) := by
      have := Real.rpow_le_rpow_of_exponent_ge hvn_pos (by linarith) (by norm_num : (1 / 2 : ℝ) ≤ 1)
      rwa [Real.rpow_one] at this
    rw [div_le_div_iff₀ (Real.rpow_pos_of_pos hvn_pos _) hδ0]
    have h2 : δ ≤ 2 * π * ‖v‖ := by
      rw [div_le_iff₀ (by positivity)] at hvn_lb; linarith
    nlinarith
  have hre2 : (-(π : ℂ) * I * (1 / 2) ^ 2 / v).re = -(π / 4) * Y := by
    rw [show -(π : ℂ) * I * (1 / 2) ^ 2 / v = ((-(π / 4) : ℝ) : ℂ) * (I * v⁻¹) by push_cast; ring,
      re_ofReal_mul, I_mul_re, inv_im, hY]
    ring
  have hn2 : ‖cexp (-(π : ℂ) * I * (1 / 2) ^ 2 / v)‖ ≤ rexp (-(1 / (72 * δ))) := by
    rw [Complex.norm_exp, hre2]
    apply Real.exp_le_exp.mpr
    have : π / 4 * (1 / (18 * π * δ)) = 1 / (72 * δ) := by field_simp; ring
    nlinarith [Real.pi_pos]
  have hτim : (-1 / v).im = Y := by
    rw [hY, div_eq_mul_inv]; simp [inv_im]; ring
  have hzim : (1 / 2 / v).im = -Y / 2 := by
    rw [hY, div_eq_mul_inv]; simp [inv_im]; ring
  have hn3 := norm_theta_aux hτim hzim hY1
  have hθn : ‖theta (cexp (-u))‖ ≤ (2 * π / δ) * rexp (-(1 / (72 * δ))) * Sθ := by
    rw [hθ, hFE, norm_mul, norm_mul]
    gcongr
  have hu4 : ‖cexp (-u / 4)‖ ≤ rexp 1 := by
    rw [Complex.norm_exp]; apply Real.exp_le_exp.mpr
    have : u.re = ζ.re := by rw [hζ]; simp
    have : (-u / 4).re = -u.re / 4 := by simp [div_ofNat_re]
    rw [this]; linarith [(abs_le.mp hre).1]
  rw [Rk, norm_mul]
  calc ‖cexp (-u / 4)‖ * ‖theta (cexp (-u))‖
      ≤ rexp 1 * ((2 * π / δ) * rexp (-(1 / (72 * δ))) * Sθ) :=
        mul_le_mul hu4 hθn (norm_nonneg _) (Real.exp_pos _).le
    _ = (2 * π * rexp 1 * Sθ / δ) * rexp (-(1 / (72 * δ))) := by ring


/-- Cauchy estimate at the cusp: all derivatives of the theta part vanish in the limit. -/
theorem tendsto_iteratedDeriv_Rk (m : ℕ) :
    Tendsto (fun α : ℝ => iteratedDeriv m Rk (-2 * α * I)) (𝓝[<] (π / 4)) (𝓝 0) := by
  have hπ3 := Real.pi_gt_three
  set K := 2 * π * rexp 1 * Sθ with hK
  have hbound : ∀ α ∈ Ioo (π / 4 - 1 / 200) (π / 4), ‖iteratedDeriv m Rk (-2 * α * I)‖ ≤
      m ! * (K * 2 ^ m) * ((π / 2 - 2 * α)⁻¹ ^ (m + 1) * rexp (-((π / 2 - 2 * α)⁻¹) / 72)) := by
    intro α hα
    set δ := π / 2 - 2 * α with hδdef
    have hδ0 : 0 < δ := by rw [hδdef]; linarith [hα.2]
    have hδ : δ ≤ 1 / 100 := by rw [hδdef]; linarith [hα.1]
    have hc : (-2 * (α : ℂ) * I) = ((-π / 2 + δ : ℝ) : ℂ) * I := by
      rw [hδdef]; push_cast; ring
    have hball : Metric.closedBall (-2 * (α : ℂ) * I) (δ / 2) ⊆ strip := by
      intro z hz
      rw [Metric.mem_closedBall, dist_eq_norm, hc] at hz
      have := (abs_im_le_norm _).trans hz
      simp only [sub_im, mul_im, ofReal_re, I_im, mul_one, ofReal_im, I_re, mul_zero,
        add_zero] at this
      show |z.im| < π / 2
      rw [abs_lt]; constructor <;> linarith [(abs_le.mp this).1, (abs_le.mp this).2]
    have hd := differentiableOn_Rk.diffContOnCl_ball hball
    have hest := norm_iteratedDeriv_le_of_forall_mem_sphere_norm_le m (half_pos hδ0) hd
      (C := (K / δ) * rexp (-(1 / (72 * δ)))) (fun z hz => by
        rw [mem_sphere_iff_norm, hc] at hz
        exact Rk_cusp_bound hδ0 hδ hz.le)
    refine hest.trans (le_of_eq ?_)
    rw [show -(1 / (72 * δ)) = -δ⁻¹ / 72 by field_simp, inv_pow, div_pow]
    field_simp
    ring
  have h1 : Tendsto (fun α : ℝ => π / 2 - 2 * α) (𝓝[<] (π / 4)) (𝓝[>] 0) := by
    refine tendsto_nhdsWithin_iff.mpr ⟨?_, ?_⟩
    · have : Tendsto (fun α : ℝ => π / 2 - 2 * α) (𝓝 (π / 4)) (𝓝 (π / 2 - 2 * (π / 4))) :=
        ((continuous_const.sub (continuous_const.mul continuous_id)).tendsto _)
      rw [show π / 2 - 2 * (π / 4) = 0 by ring] at this
      exact this.mono_left nhdsWithin_le_nhds
    · filter_upwards [self_mem_nhdsWithin] with α hα
      simp only [mem_Ioi]; linarith [mem_Iio.mp hα]
  have hδlim := tendsto_inv_nhdsGT_zero.comp h1
  have hg : Tendsto (fun x : ℝ => x ^ (m + 1) * rexp (-x / 72)) atTop (𝓝 0) := by
    have h2 := ((tendsto_pow_mul_exp_neg_atTop_nhds_zero (m + 1)).comp
      (tendsto_id.atTop_div_const (by norm_num : (0 : ℝ) < 72))).const_mul ((72 : ℝ) ^ (m + 1))
    rw [mul_zero] at h2
    refine h2.congr (fun x => ?_)
    simp only [Function.comp_apply, id]
    rw [div_pow, neg_div]; field_simp
  have hB := (hg.comp hδlim).const_mul (m ! * (K * 2 ^ m))
  rw [mul_zero] at hB
  refine squeeze_zero_norm' ?_ hB
  filter_upwards [Ioo_mem_nhdsLT (by linarith [Real.pi_pos] : π / 4 - 1 / 200 < π / 4)] with α hα
  exact hbound α hα

end Fayed.Hardy
