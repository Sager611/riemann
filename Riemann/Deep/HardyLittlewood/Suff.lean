import Riemann.Tree.R007

/-!
# Hardy–Littlewood (1916) criterion, sufficiency:
# F(x) = O(x^{-1/4+ε}) for all ε > 0 implies RH

Classical form: `F(x) = Σ_{k≥1} (-x)^k / (k! ζ(2k+1)) = Σ_n (μ(n)/n) e^{-x/n²}`.

Lean form. The literal Möbius series `∑' n, (μ n / n) * exp (-(x/n²))` is NOT absolutely
summable (its terms are ≈ μ(n)/n), so Lean's `tsum` makes it identically `0`; with that
definition the statement would just be RH itself. We therefore use the absolutely convergent
form
  `hlF x = Σ_n (μ(n)/n) (e^{-x/n²} - 1)`,
which equals `Σ_{k≥1} (-x)^k/(k! ζ(2k+1))` unconditionally (expand `e^{-x/n²} - 1`, swap the
absolutely convergent double sum, `Σ_n μ(n) n^{-2k-1} = 1/ζ(2k+1)`), and equals the classical
conditionally convergent `Σ_n (μ(n)/n) e^{-x/n²}` because `Σ_n μ(n)/n = 0` (Landau). Neither
identity is formalized here; the criterion is stated for `hlF`.

Proof:
* `hlF' = -g`, `g x = Σ μ(n) n^{-3} e^{-x/n²}`; `mellin_G`: `M[g](z) = Γ(z) L(μ, 3-2z)`
  for `0 < Re z < 1` (Mathlib's `hasSum_mellin`).
* `mellin_F_eq`: integration by parts, `M[hlF](s) = M[g](s+1)/s = Γ(s) L(μ, 1-2s)` for
  `-1 < Re s < -1/2` (unconditional bounds `|hlF x| ≤ S_b x^b`, `0 < b ≤ 1`).
* `differentiableAt_M`: the hypothesis makes `M[hlF]` holomorphic on `-1 < Re s < 1/4`.
* identity theorem on two convex domains avoiding `s = 0`: `ζ(1-2s) M(s) = Γ(s) ≠ 0` on
  `0 < Re s < 1/4`, i.e. `ζ(w) ≠ 0` for `1/2 < Re w < 1`; then `RH_iff_QuasiRH_half`.
-/

open Complex Asymptotics Filter Topology MeasureTheory Set

namespace Fayed.HL

/-- n-th term of `hlF`. -/
noncomputable def term (n : ℕ) (x : ℝ) : ℝ :=
  (ArithmeticFunction.moebius n : ℝ) / (n : ℝ) * (Real.exp (-(x / (n : ℝ) ^ 2)) - 1)

/-- Hardy–Littlewood function, absolutely convergent Möbius form:
`Σ (μ(n)/n)(e^{-x/n²} - 1) = Σ_{k≥1} (-x)^k/(k! ζ(2k+1))`. -/
noncomputable def hlF (x : ℝ) : ℝ :=
  ∑' n : ℕ, (ArithmeticFunction.moebius n : ℝ) / (n : ℝ) * (Real.exp (-(x / (n : ℝ) ^ 2)) - 1)

lemma hlF_eq (x : ℝ) : hlF x = ∑' n, term n x := rfl

/-- n-th term of `g = -hlF'`. -/
noncomputable def gterm (n : ℕ) (x : ℝ) : ℝ :=
  (ArithmeticFunction.moebius n : ℝ) / (n : ℝ) ^ 3 * Real.exp (-(x / (n : ℝ) ^ 2))

/-- `g = -hlF' = Σ μ(n) n^{-3} e^{-x/n²}`. -/
noncomputable def g (x : ℝ) : ℝ := ∑' n, gterm n x

/-- Complexified `hlF`. -/
noncomputable def F (x : ℝ) : ℂ := (hlF x : ℂ)

/-- Complexified `g`. -/
noncomputable def G (x : ℝ) : ℂ := (g x : ℂ)

/-- M(s) = ∫_0^∞ hlF(x) x^{s-1} dx. -/
noncomputable def M (s : ℂ) : ℂ := mellin F s

/-- Dirichlet coefficients μ(n)/n³. -/
noncomputable def a (n : ℕ) : ℂ := (ArithmeticFunction.moebius n : ℂ) / (n : ℂ) ^ 3

/-- Exponential rates 1/n². -/
noncomputable def p (n : ℕ) : ℝ := ((n : ℝ) ^ 2)⁻¹

/-- `S b = Σ 1/n^{1+2b}`. -/
noncomputable def S (b : ℝ) : ℝ := ∑' n : ℕ, 1 / (n : ℝ) ^ (1 + 2 * b)

lemma abs_mu_le (n : ℕ) : |(ArithmeticFunction.moebius n : ℝ)| ≤ 1 := by
  exact_mod_cast ArithmeticFunction.abs_moebius_le_one

lemma summable_inv_sq : Summable (fun n : ℕ => 1 / (n : ℝ) ^ 2) :=
  Real.summable_one_div_nat_pow.mpr one_lt_two

lemma summable_inv_cube : Summable (fun n : ℕ => 1 / (n : ℝ) ^ 3) :=
  Real.summable_one_div_nat_pow.mpr (by norm_num)

lemma summable_S {b : ℝ} (hb : 0 < b) : Summable (fun n : ℕ => 1 / (n : ℝ) ^ (1 + 2 * b)) :=
  Real.summable_one_div_nat_rpow.mpr (by linarith)

lemma exp_neg_le_rpow_neg {y c : ℝ} (hy : 0 < y) (hc0 : 0 ≤ c) (hc1 : c ≤ 1) :
    Real.exp (-y) ≤ y ^ (-c) := by
  rcases le_or_gt 1 y with h | h
  · have h1 : Real.exp (-y) ≤ y⁻¹ := by
      rw [Real.exp_neg]; exact inv_anti₀ hy (by linarith [Real.add_one_le_exp y])
    have h2 : y⁻¹ ≤ y ^ (-c) := by
      rw [← Real.rpow_neg_one]; exact Real.rpow_le_rpow_of_exponent_le h (by linarith)
    linarith
  · have h1 : 1 ≤ y ^ (-c) := Real.one_le_rpow_of_pos_of_le_one_of_nonpos hy h.le (by linarith)
    have h2 : Real.exp (-y) ≤ 1 := Real.exp_le_one_iff.mpr (by linarith)
    linarith

lemma abs_exp_neg_sub_one_le {y b : ℝ} (hy : 0 ≤ y) (hb0 : 0 < b) (hb1 : b ≤ 1) :
    |Real.exp (-y) - 1| ≤ y ^ b := by
  have he : Real.exp (-y) ≤ 1 := Real.exp_le_one_iff.mpr (by linarith)
  rw [abs_of_nonpos (by linarith)]
  rcases hy.eq_or_lt with h0 | hy'
  · subst h0; simp [Real.zero_rpow hb0.ne']
  rcases le_or_gt y 1 with h | h
  · have h1 : 1 - y ≤ Real.exp (-y) := by linarith [Real.add_one_le_exp (-y)]
    have h2 : y ≤ y ^ b := by
      simpa using Real.rpow_le_rpow_of_exponent_ge hy' h hb1
    linarith
  · have : 1 ≤ y ^ b := Real.one_le_rpow h.le hb0.le
    linarith [Real.exp_pos (-y)]

lemma abs_term_le (n : ℕ) {b x : ℝ} (hb0 : 0 < b) (hb1 : b ≤ 1) (hx : 0 ≤ x) :
    |term n x| ≤ x ^ b * (1 / (n : ℝ) ^ (1 + 2 * b)) := by
  rcases eq_or_ne n 0 with rfl | hn
  · simp only [term, ArithmeticFunction.map_zero, Int.cast_zero, zero_div, zero_mul, abs_zero]
    positivity
  have hn' : (0 : ℝ) < n := by exact_mod_cast Nat.pos_of_ne_zero hn
  have hy : 0 ≤ x / (n : ℝ) ^ 2 := by positivity
  have h1 := abs_exp_neg_sub_one_le hy hb0 hb1
  have h2 : (x / (n : ℝ) ^ 2) ^ b = x ^ b / (n : ℝ) ^ (2 * b) := by
    rw [Real.div_rpow hx (by positivity), ← Real.rpow_natCast, ← Real.rpow_mul hn'.le]
    norm_num
  have h3 : (n : ℝ) ^ (1 + 2 * b) = n * (n : ℝ) ^ (2 * b) := by
    rw [Real.rpow_add hn', Real.rpow_one]
  have h4 : (0 : ℝ) < (n : ℝ) ^ (2 * b) := Real.rpow_pos_of_pos hn' _
  rw [term, abs_mul, abs_div, Nat.abs_cast]
  calc |(ArithmeticFunction.moebius n : ℝ)| / n * |Real.exp (-(x / (n : ℝ) ^ 2)) - 1|
      ≤ 1 / n * (x / (n : ℝ) ^ 2) ^ b :=
        mul_le_mul (by gcongr; exact abs_mu_le n) h1 (abs_nonneg _) (by positivity)
    _ = x ^ b * (1 / (n : ℝ) ^ (1 + 2 * b)) := by
        rw [h2, h3]; field_simp

lemma summable_term {x : ℝ} (hx : 0 ≤ x) : Summable (fun n => term n x) :=
  Summable.of_norm_bounded ((summable_S one_pos).mul_left (x ^ (1 : ℝ))) fun n => by
    rw [Real.norm_eq_abs]; exact abs_term_le n one_pos le_rfl hx

lemma abs_hlF_le {b x : ℝ} (hb0 : 0 < b) (hb1 : b ≤ 1) (hx : 0 ≤ x) :
    |hlF x| ≤ x ^ b * S b := by
  rw [← Real.norm_eq_abs, hlF_eq]
  exact tsum_of_norm_bounded ((summable_S hb0).hasSum.mul_left (x ^ b))
    (fun n => by rw [Real.norm_eq_abs]; exact abs_term_le n hb0 hb1 hx)

lemma abs_gterm_le_inv_cube (n : ℕ) {x : ℝ} (hx : 0 ≤ x) : |gterm n x| ≤ 1 / (n : ℝ) ^ 3 := by
  have he : Real.exp (-(x / (n : ℝ) ^ 2)) ≤ 1 :=
    Real.exp_le_one_iff.mpr (by have : 0 ≤ x / (n : ℝ) ^ 2 := by positivity
                                linarith)
  rw [gterm, abs_mul, abs_div, abs_pow, Nat.abs_cast, abs_of_pos (Real.exp_pos _)]
  calc |(ArithmeticFunction.moebius n : ℝ)| / (n : ℝ) ^ 3 * Real.exp (-(x / (n : ℝ) ^ 2))
      ≤ 1 / (n : ℝ) ^ 3 * 1 :=
        mul_le_mul (by gcongr; exact abs_mu_le n) he (Real.exp_pos _).le (by positivity)
    _ = 1 / (n : ℝ) ^ 3 := mul_one _

lemma abs_gterm_le (n : ℕ) {x : ℝ} (hx : 0 < x) :
    |gterm n x| ≤ x ^ (-(1 / 2) : ℝ) * (1 / (n : ℝ) ^ 2) := by
  rcases eq_or_ne n 0 with rfl | hn
  · simp only [gterm, ArithmeticFunction.map_zero, Int.cast_zero, zero_div, zero_mul, abs_zero]
    positivity
  have hn' : (0 : ℝ) < n := by exact_mod_cast Nat.pos_of_ne_zero hn
  have he := exp_neg_le_rpow_neg (y := x / (n : ℝ) ^ 2) (by positivity) (c := 1 / 2)
    (by norm_num) (by norm_num)
  have hsq : ((n : ℝ) ^ 2) ^ (1 / 2 : ℝ) = n := by
    rw [show (1 / 2 : ℝ) = ((2 : ℕ) : ℝ)⁻¹ by norm_num,
      Real.pow_rpow_inv_natCast hn'.le (by norm_num)]
  have h2 : (x / (n : ℝ) ^ 2) ^ (-(1 / 2) : ℝ) = x ^ (-(1 / 2) : ℝ) * n := by
    rw [Real.div_rpow hx.le (by positivity), Real.rpow_neg (by positivity : (0 : ℝ) ≤ (n : ℝ) ^ 2),
      hsq]
    field_simp
  rw [gterm, abs_mul, abs_div, abs_pow, Nat.abs_cast, abs_of_pos (Real.exp_pos _)]
  calc |(ArithmeticFunction.moebius n : ℝ)| / (n : ℝ) ^ 3 * Real.exp (-(x / (n : ℝ) ^ 2))
      ≤ 1 / (n : ℝ) ^ 3 * (x / (n : ℝ) ^ 2) ^ (-(1 / 2) : ℝ) :=
        mul_le_mul (by gcongr; exact abs_mu_le n) he (Real.exp_pos _).le (by positivity)
    _ = x ^ (-(1 / 2) : ℝ) * (1 / (n : ℝ) ^ 2) := by
        rw [h2]; field_simp

lemma hasDerivAt_term (n : ℕ) (y : ℝ) : HasDerivAt (term n) (-gterm n y) y := by
  have := ((((hasDerivAt_id y).div_const ((n : ℝ) ^ 2)).neg.exp).sub_const 1).const_mul
    ((ArithmeticFunction.moebius n : ℝ) / (n : ℝ))
  convert this using 1
  · funext z; simp only [term, Pi.neg_apply, id]
  · simp only [gterm, Pi.neg_apply, id]; ring

lemma hasDerivAt_hlF {x : ℝ} (hx : 0 < x) : HasDerivAt hlF (-g x) x := by
  have := hasDerivAt_tsum_of_isPreconnected (u := fun n : ℕ => 1 / (n : ℝ) ^ 3) (g := term)
    (g' := fun n y => -gterm n y) (t := Ioi 0) (y₀ := 1) summable_inv_cube isOpen_Ioi
    isPreconnected_Ioi (fun n y _ => hasDerivAt_term n y)
    (fun n y hy => by rw [norm_neg, Real.norm_eq_abs]; exact abs_gterm_le_inv_cube n (le_of_lt hy))
    (by norm_num : (1 : ℝ) ∈ Ioi 0) (summable_term zero_le_one) hx
  rw [tsum_neg] at this
  exact this

lemma locallyIntegrableOn_F : LocallyIntegrableOn F (Ioi 0) := by
  refine ContinuousOn.locallyIntegrableOn (fun x hx => ?_) measurableSet_Ioi
  exact (continuous_ofReal.continuousAt.comp (hasDerivAt_hlF hx).continuousAt).continuousWithinAt

lemma continuousOn_g : ContinuousOn g (Ioi 0) := by
  refine continuousOn_tsum (u := fun n : ℕ => 1 / (n : ℝ) ^ 3) (fun n => ?_) summable_inv_cube
    (fun n x hx => ?_)
  · have : Continuous (gterm n) := by unfold gterm; fun_prop
    exact this.continuousOn
  · rw [Real.norm_eq_abs]; exact abs_gterm_le_inv_cube n (le_of_lt hx)

lemma locallyIntegrableOn_G : LocallyIntegrableOn G (Ioi 0) :=
  (continuous_ofReal.comp_continuousOn continuousOn_g).locallyIntegrableOn measurableSet_Ioi

lemma isBigO_F_atTop {θ : ℝ} (h : ∃ C : ℝ, ∀ x : ℝ, 1 ≤ x → |hlF x| ≤ C * x ^ θ) :
    F =O[atTop] (fun t : ℝ => t ^ (-(-θ))) := by
  obtain ⟨C, hC⟩ := h
  rw [neg_neg]
  refine IsBigO.of_bound C ?_
  filter_upwards [eventually_ge_atTop 1] with x hx
  rw [F, Complex.norm_real, Real.norm_eq_abs,
    Real.norm_of_nonneg (Real.rpow_nonneg (by linarith) _)]
  exact hC x hx

lemma isBigO_F_zero : F =O[𝓝[>] 0] (fun t : ℝ => t ^ (-(-1 : ℝ))) := by
  refine IsBigO.of_bound (S 1) ?_
  filter_upwards [self_mem_nhdsWithin] with x (hx : 0 < x)
  rw [neg_neg, F, Complex.norm_real, Real.norm_eq_abs,
    Real.norm_of_nonneg (Real.rpow_nonneg hx.le _), mul_comm]
  exact abs_hlF_le one_pos le_rfl hx.le

lemma isBigO_G_atTop : G =O[atTop] (fun t : ℝ => t ^ (-(1 / 2 : ℝ))) := by
  refine IsBigO.of_bound (∑' n : ℕ, 1 / (n : ℝ) ^ 2) ?_
  filter_upwards [eventually_ge_atTop 1] with x hx
  have hx0 : 0 < x := by linarith
  rw [G, Complex.norm_real, Real.norm_of_nonneg (Real.rpow_nonneg hx0.le _), mul_comm, g]
  exact tsum_of_norm_bounded (summable_inv_sq.hasSum.mul_left _)
    (fun n => by rw [Real.norm_eq_abs]; exact abs_gterm_le n hx0)

lemma isBigO_G_zero : G =O[𝓝[>] 0] (fun t : ℝ => t ^ (-(0 : ℝ))) := by
  refine IsBigO.of_bound (∑' n : ℕ, 1 / (n : ℝ) ^ 3) ?_
  filter_upwards [self_mem_nhdsWithin] with x (hx : 0 < x)
  rw [neg_zero, Real.rpow_zero, norm_one, mul_one, G, Complex.norm_real, g]
  exact tsum_of_norm_bounded summable_inv_cube.hasSum
    (fun n => by rw [Real.norm_eq_abs]; exact abs_gterm_le_inv_cube n hx.le)

lemma differentiableAt_M
    (h : ∀ ε : ℝ, 0 < ε → ∃ C : ℝ, ∀ x : ℝ, 1 ≤ x → |hlF x| ≤ C * x ^ (-(1 / 4 : ℝ) + ε))
    {s : ℂ} (h1 : -1 < s.re) (h2 : s.re < 1 / 4) : DifferentiableAt ℂ M s := by
  have hε : 0 < (1 / 4 - s.re) / 2 := by linarith
  exact mellin_differentiableAt_of_isBigO_rpow (b := -1) locallyIntegrableOn_F
    (isBigO_F_atTop (h _ hε)) (by linarith) isBigO_F_zero h1

lemma norm_a_le (n : ℕ) : ‖a n‖ ≤ 1 / (n : ℝ) ^ 3 := by
  rw [a, norm_div, norm_pow, Complex.norm_natCast, Complex.norm_intCast]
  exact div_le_div_of_nonneg_right (abs_mu_le n) (by positivity)

lemma hasSum_G {t : ℝ} (ht : 0 ≤ t) :
    HasSum (fun n => a n * (Real.exp (-p n * t) : ℂ)) (G t) := by
  have hs : Summable (fun n => gterm n t) :=
    Summable.of_norm_bounded summable_inv_cube
      (fun n => by rw [Real.norm_eq_abs]; exact abs_gterm_le_inv_cube n ht)
  have := Complex.hasSum_ofReal.mpr hs.hasSum
  rw [show G t = ((∑' n, gterm n t : ℝ) : ℂ) from rfl]
  convert this using 1
  funext n
  rw [show -p n * t = -(t / (n : ℝ) ^ 2) by rw [p]; ring]
  simp only [a, gterm]
  push_cast
  ring

lemma summable_coeff {z : ℂ} (h1 : z.re < 1) :
    Summable fun n => ‖a n‖ / p n ^ z.re := by
  refine Summable.of_nonneg_of_le
    (fun n => div_nonneg (norm_nonneg _) (Real.rpow_nonneg (by rw [p]; positivity) _)) (fun n => ?_)
    (Real.summable_one_div_nat_rpow.mpr (show (1 : ℝ) < 3 - 2 * z.re by linarith))
  rcases eq_or_ne n 0 with rfl | hn
  · simp [a]; positivity
  have hn' : (0 : ℝ) < n := by exact_mod_cast Nat.pos_of_ne_zero hn
  rw [p, Real.inv_rpow (by positivity), div_inv_eq_mul]
  have h2 : ((n : ℝ) ^ 2) ^ z.re = (n : ℝ) ^ (2 * z.re) := by
    rw [← Real.rpow_natCast, ← Real.rpow_mul hn'.le]
    norm_num
  have h3 : (n : ℝ) ^ (3 - 2 * z.re) = (n : ℝ) ^ 3 / (n : ℝ) ^ (2 * z.re) := by
    rw [Real.rpow_sub hn']
    norm_num
  have h4 : (0 : ℝ) < (n : ℝ) ^ (2 * z.re) := Real.rpow_pos_of_pos hn' _
  rw [h2, h3]
  calc ‖a n‖ * (n : ℝ) ^ (2 * z.re)
      ≤ 1 / (n : ℝ) ^ 3 * (n : ℝ) ^ (2 * z.re) := by gcongr; exact norm_a_le n
    _ = 1 / ((n : ℝ) ^ 3 / (n : ℝ) ^ (2 * z.re)) := by field_simp

/-- Mellin transform of `g`: Γ(z) L(μ, 3 - 2z). -/
lemma mellin_G {z : ℂ} (h0 : 0 < z.re) (h1 : z.re < 1) :
    mellin G z = Gamma z * LSeries (fun n => (ArithmeticFunction.moebius n : ℂ)) (3 - 2 * z) := by
  have hp : ∀ n, a n = 0 ∨ 0 < p n := by
    intro n
    rcases eq_or_ne n 0 with rfl | hn
    · left; simp [a]
    · right; rw [p]; have : (0 : ℝ) < n := by exact_mod_cast Nat.pos_of_ne_zero hn
      positivity
  have key := hasSum_mellin hp h0 (fun t ht => hasSum_G (le_of_lt ht)) (summable_coeff h1)
  rw [← key.tsum_eq, LSeries, ← tsum_mul_left]
  congr 1
  funext n
  rcases eq_or_ne n 0 with rfl | hn
  · simp [a, LSeries.term_zero]
  rw [LSeries.term_of_ne_zero hn]
  have hn' : (n : ℂ) ≠ 0 := by exact_mod_cast hn
  have hy : (n : ℂ) ^ z ≠ 0 := (cpow_ne_zero_iff_of_exponent_ne_zero (by
    intro h0'; rw [h0'] at h0; simp at h0)).mpr hn'
  have hpn : ((p n : ℝ) : ℂ) = ((n : ℂ) * (n : ℂ))⁻¹ := by
    rw [p]; push_cast; ring
  have harg : ((n : ℂ) * (n : ℂ)).arg ≠ Real.pi := by
    rw [← Nat.cast_mul, natCast_arg]; exact Real.pi_pos.ne
  rw [hpn, inv_cpow _ _ harg, natCast_mul_natCast_cpow, cpow_sub _ _ hn', cpow_ofNat_mul,
    a]
  have h3 : (n : ℂ) ^ (3 : ℂ) = (n : ℂ) ^ (3 : ℕ) := by exact_mod_cast cpow_natCast (n : ℂ) 3
  rw [h3]
  field_simp

lemma norm_uv_le {s : ℂ} {b x : ℝ} (hb0 : 0 < b) (hb1 : b ≤ 1) (hx : 0 < x) :
    ‖F x * ((x : ℂ) ^ s / s)‖ ≤ S b / ‖s‖ * x ^ (b + s.re) := by
  rw [norm_mul, norm_div, Complex.norm_cpow_eq_rpow_re_of_pos hx, F, Complex.norm_real,
    Real.norm_eq_abs, Real.rpow_add hx]
  have := abs_hlF_le hb0 hb1 hx.le
  calc |hlF x| * (x ^ s.re / ‖s‖) ≤ x ^ b * S b * (x ^ s.re / ‖s‖) := by gcongr
    _ = S b / ‖s‖ * (x ^ b * x ^ s.re) := by ring

/-- Integration by parts: M[hlF](s) = M[g](s+1)/s for -1 < Re s < -1/2. -/
lemma mellin_F_eq {s : ℂ} (hs1 : -1 < s.re) (hs2 : s.re < -1 / 2) :
    mellin F s = mellin G (s + 1) / s := by
  have hs0 : s ≠ 0 := by intro h0; rw [h0] at hs2; simp at hs2; linarith
  have hFc : IntegrableOn (fun t : ℝ => (t : ℂ) ^ (s - 1) • F t) (Ioi 0) :=
    mellinConvergent_of_isBigO_rpow (a := -(1 / 4)) (b := -1) locallyIntegrableOn_F
      (isBigO_F_atTop (θ := 1 / 4) ⟨S (1 / 4), fun x hx => by
        rw [mul_comm]; exact abs_hlF_le (by norm_num) (by norm_num) (by linarith)⟩)
      (by linarith) isBigO_F_zero hs1
  have hGc : IntegrableOn (fun t : ℝ => (t : ℂ) ^ (s + 1 - 1) • G t) (Ioi 0) :=
    mellinConvergent_of_isBigO_rpow (a := 1 / 2) (b := 0) locallyIntegrableOn_G
      isBigO_G_atTop (by simp; linarith) isBigO_G_zero (by simp; linarith)
  have hu : ∀ x ∈ Ioi (0 : ℝ), HasDerivAt F (-G x) x := by
    intro x hx
    have := (hasDerivAt_hlF hx).ofReal_comp
    rw [show -G x = ((-(g x) : ℝ) : ℂ) by simp [G]]
    exact this
  have hv : ∀ x ∈ Ioi (0 : ℝ),
      HasDerivAt (fun y : ℝ => (y : ℂ) ^ s / s) ((x : ℂ) ^ (s - 1)) x := by
    intro x hx
    have := hasDerivAt_ofReal_cpow_const' (ne_of_gt hx) (r := s - 1)
      (by intro h0; apply hs0; linear_combination h0)
    rwa [sub_add_cancel] at this
  have huv' : IntegrableOn (F * fun x : ℝ => (x : ℂ) ^ (s - 1)) (Ioi 0) :=
    hFc.congr_fun (fun x _ => by simp only [Pi.mul_apply, smul_eq_mul]; ring) measurableSet_Ioi
  have hu'v : IntegrableOn ((fun x => -G x) * fun y : ℝ => (y : ℂ) ^ s / s) (Ioi 0) :=
    IntegrableOn.congr_fun (Integrable.neg (Integrable.div_const hGc s)) (fun x _ => by
      simp only [Pi.mul_apply, Pi.neg_apply, smul_eq_mul, add_sub_cancel_right]; ring)
      measurableSet_Ioi
  have h_zero : Tendsto (F * fun y : ℝ => (y : ℂ) ^ s / s) (𝓝[>] 0) (𝓝 0) := by
    refine squeeze_zero_norm' (a := fun x => S 1 / ‖s‖ * x ^ (1 + s.re)) ?_ ?_
    · filter_upwards [self_mem_nhdsWithin] with x (hx : 0 < x)
      exact norm_uv_le one_pos le_rfl hx
    · have hc : 0 < 1 + s.re := by linarith
      have := ((Real.continuousAt_rpow_const 0 (1 + s.re) (Or.inr hc.le)).tendsto.mono_left
        (nhdsWithin_le_nhds (s := Ioi 0))).const_mul (S 1 / ‖s‖)
      rw [Real.zero_rpow hc.ne', mul_zero] at this
      exact this
  have h_infty : Tendsto (F * fun y : ℝ => (y : ℂ) ^ s / s) atTop (𝓝 0) := by
    refine squeeze_zero_norm' (a := fun x => S (1 / 4) / ‖s‖ * x ^ (1 / 4 + s.re)) ?_ ?_
    · filter_upwards [eventually_gt_atTop 0] with x hx
      exact norm_uv_le (by norm_num) (by norm_num) hx
    · have := (tendsto_rpow_neg_atTop (y := -(1 / 4 + s.re)) (by linarith)).const_mul
        (S (1 / 4) / ‖s‖)
      rw [neg_neg, mul_zero] at this
      exact this
  have ibp := integral_Ioi_mul_deriv_eq_deriv_mul hu hv huv' hu'v h_zero h_infty
  rw [sub_zero, zero_sub] at ibp
  rw [mellin, mellin]
  have e1 : (fun t : ℝ => (t : ℂ) ^ (s - 1) • F t) = fun t => F t * (t : ℂ) ^ (s - 1) := by
    funext t; rw [smul_eq_mul, mul_comm]
  rw [e1, ibp, ← integral_div, ← integral_neg]
  congr 1
  funext t
  simp only [smul_eq_mul, add_sub_cancel_right]
  ring

/-- M[hlF](s) = Γ(s) L(μ, 1 - 2s) = Γ(s)/ζ(1-2s) on -1 < Re s < -1/2. -/
lemma mellin_eq {s : ℂ} (hs1 : -1 < s.re) (hs2 : s.re < -1 / 2) :
    M s = Gamma s * LSeries (fun n => (ArithmeticFunction.moebius n : ℂ)) (1 - 2 * s) := by
  have hs0 : s ≠ 0 := by intro h0; rw [h0] at hs2; simp at hs2; linarith
  rw [M, mellin_F_eq hs1 hs2, mellin_G (by simp; linarith) (by simp; linarith),
    Complex.Gamma_add_one s hs0, show (3 : ℂ) - 2 * (s + 1) = 1 - 2 * s by ring]
  field_simp

lemma isLinearMap_im_sub_re : IsLinearMap ℝ (fun z : ℂ => z.im - z.re) :=
  ⟨fun x y => by simp; ring, fun c x => by simp; ring⟩

lemma isLinearMap_re_add_im : IsLinearMap ℝ (fun z : ℂ => z.re + z.im) :=
  ⟨fun x y => by simp; ring, fun c x => by simp; ring⟩

/-- Identity theorem: ζ(1-2z) M(z) = Γ(z) on any convex open D inside -1 < Re z < 1/4, off 0,
meeting Re z < -1/2. -/
lemma eq_on_convex
    (h : ∀ ε : ℝ, 0 < ε → ∃ C : ℝ, ∀ x : ℝ, 1 ≤ x → |hlF x| ≤ C * x ^ (-(1 / 4 : ℝ) + ε))
    {D : Set ℂ} (hDo : IsOpen D) (hDc : Convex ℝ D)
    (hD : ∀ z ∈ D, -1 < z.re ∧ z.re < 1 / 4 ∧ z ≠ 0)
    {c : ℂ} (hcD : c ∈ D) (hc : c.re < -1 / 2) {s : ℂ} (hs : s ∈ D) :
    riemannZeta (1 - 2 * s) * M s = Gamma s := by
  have hd : DifferentiableOn ℂ (fun z => riemannZeta (1 - 2 * z) * M z) D := by
    intro z hz
    obtain ⟨h1, h2, h3⟩ := hD z hz
    have h2z : 1 - 2 * z ≠ 1 := fun h0 => h3 (by linear_combination -h0 / 2)
    exact (((differentiableAt_riemannZeta h2z).comp z
      ((differentiableAt_const _).sub (differentiableAt_id.const_mul 2))).mul
      (differentiableAt_M h h1 h2)).differentiableWithinAt
  have hg : DifferentiableOn ℂ (fun z => Gamma z) D := by
    intro z hz
    obtain ⟨h1, _, h3⟩ := hD z hz
    have hG : ∀ m : ℕ, z ≠ -m := by
      intro m hm
      have hre := congrArg Complex.re hm
      simp at hre
      have hm1 : (m : ℝ) < 1 := by linarith
      have hm0 : m = 0 := by exact_mod_cast Nat.lt_one_iff.mp (by exact_mod_cast hm1)
      subst hm0
      simp at hm
      exact h3 hm
    exact (differentiableAt_Gamma _ hG).differentiableWithinAt
  have hev : (fun z => riemannZeta (1 - 2 * z) * M z) =ᶠ[𝓝 c] (fun z => Gamma z) := by
    have hmem : {z : ℂ | z.re < -1 / 2} ∩ D ∈ 𝓝 c :=
      ((isOpen_lt Complex.continuous_re continuous_const).inter hDo).mem_nhds ⟨hc, hcD⟩
    filter_upwards [hmem] with z hz
    obtain ⟨hz1, hzD⟩ := hz
    have hz1 : z.re < -1 / 2 := hz1
    have hz2 : 1 < (1 - 2 * z).re := by
      rw [show (1 - 2 * z).re = 1 - 2 * z.re by simp]; linarith
    have hL := ArithmeticFunction.LSeries_zeta_mul_Lseries_moebius hz2
    rw [mellin_eq (hD z hzD).1 hz1, ← ArithmeticFunction.LSeries_zeta_eq_riemannZeta hz2]
    linear_combination (Gamma z) * hL
  exact (hd.analyticOnNhd hDo).eqOn_of_preconnected_of_eventuallyEq (hg.analyticOnNhd hDo)
    hDc.isPreconnected hcD hev hs

lemma zeta_one_sub_two_mul_ne_zero
    (h : ∀ ε : ℝ, 0 < ε → ∃ C : ℝ, ∀ x : ℝ, 1 ≤ x → |hlF x| ≤ C * x ^ (-(1 / 4 : ℝ) + ε))
    {s : ℂ} (h1 : 0 < s.re) (h2 : s.re < 1 / 4) : riemannZeta (1 - 2 * s) ≠ 0 := by
  have hS : IsOpen {z : ℂ | -1 < z.re ∧ z.re < 1 / 4} :=
    (isOpen_lt continuous_const Complex.continuous_re).inter
      (isOpen_lt Complex.continuous_re continuous_const)
  have hSc : Convex ℝ {z : ℂ | -1 < z.re ∧ z.re < 1 / 4} :=
    (convex_halfSpace_re_gt _).inter (convex_halfSpace_re_lt _)
  have key : riemannZeta (1 - 2 * s) * M s = Gamma s := by
    by_cases hpos : 0 < s.re + s.im
    · refine eq_on_convex h (D := {z : ℂ | -1 < z.re ∧ z.re < 1 / 4} ∩
          {z : ℂ | 0 < z.re + z.im})
        (hS.inter (isOpen_lt continuous_const
          (Complex.continuous_re.add Complex.continuous_im)))
        (hSc.inter (convex_halfSpace_gt isLinearMap_re_add_im _)) ?_
        (c := -3 / 4 + I) ⟨⟨by simp; norm_num, by simp; norm_num⟩, by simp; norm_num⟩
        (by simp; norm_num) ⟨⟨by linarith, h2⟩, hpos⟩
      rintro z ⟨⟨hz1, hz2⟩, hz3⟩
      refine ⟨hz1, hz2, fun hz => ?_⟩
      subst hz; norm_num at hz3
    · push Not at hpos
      refine eq_on_convex h (D := {z : ℂ | -1 < z.re ∧ z.re < 1 / 4} ∩
          {z : ℂ | z.im - z.re < 0})
        (hS.inter (isOpen_lt (Complex.continuous_im.sub Complex.continuous_re)
          continuous_const))
        (hSc.inter (convex_halfSpace_lt isLinearMap_im_sub_re _)) ?_
        (c := -3 / 4 - I) ⟨⟨by simp; norm_num, by simp; norm_num⟩, by simp; norm_num⟩
        (by simp; norm_num) ⟨⟨by linarith, h2⟩, show s.im - s.re < 0 by linarith⟩
      rintro z ⟨⟨hz1, hz2⟩, hz3⟩
      refine ⟨hz1, hz2, fun hz => ?_⟩
      subst hz; norm_num at hz3
  intro h0
  rw [h0, zero_mul] at key
  exact Gamma_ne_zero_of_re_pos h1 key.symm

/-- **Hardy–Littlewood criterion (sufficiency).** If `hlF(x) = O(x^{-1/4+ε})` for every
`ε > 0`, then RH. -/
theorem RH_of_hlBound
    (h : ∀ ε : ℝ, 0 < ε → ∃ C : ℝ, ∀ x : ℝ, 1 ≤ x → |hlF x| ≤ C * x ^ (-(1 / 4 : ℝ) + ε)) :
    RiemannHypothesis := by
  rw [Fayed.R7.RH_iff_QuasiRH_half]
  intro w hw _
  rcases le_or_gt 1 w.re with h1 | h1
  · exact riemannZeta_ne_zero_of_one_le_re h1
  · have := zeta_one_sub_two_mul_ne_zero h (s := (1 - w) / 2) (by simp; linarith)
      (by simp; linarith)
    rwa [show 1 - 2 * ((1 - w) / 2) = w by ring] at this

end Fayed.HL

#print axioms Fayed.HL.RH_of_hlBound
