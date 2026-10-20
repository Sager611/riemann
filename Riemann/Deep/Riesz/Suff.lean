import Riemann.Tree.R007
import Riemann.Tree.R013

/-!
# Riesz's criterion, sufficiency: R(x) = O(x^{1/4+ε}) for all ε > 0 implies RH

Möbius form of the Riesz function: `rieszM x = Σ_n μ(n) (x/n²) e^{-x/n²}`.
(Classically `Σ_n μ(n)(x/n²)e^{-x/n²} = Σ_k (-x)^{k+1}/(k! ζ(2k+2))`, i.e. `rieszM` is
`Fayed.R13.rieszSeries` from `Riemann/Tree/R013.lean`; that power-series identity is not
formalized here — the criterion is stated directly for the Möbius form.)

Proof:
* `mellin_eq`: for 1/2 < Re s < 1, `∫_0^∞ rieszM(x) x^{-s-1} dx = Γ(1-s) L(μ, 2s)`
  (Mathlib's `hasSum_mellin` after writing `rieszM x = x · Σ (μ(n)/n²) e^{-x/n²}`).
* `differentiableAt_M`: the bound plus |rieszM x| ≤ ζ(2) x makes the left side holomorphic on
  1/4 < Re s < 1 (`mellin_differentiableAt_of_isBigO_rpow`).
* identity theorem on two convex domains avoiding s = 1/2: ζ(2s) M(s) = Γ(1-s) ≠ 0 on
  1/4 < Re s < 1/2, so ζ(w) ≠ 0 for 1/2 < Re w < 1; then `RH_iff_QuasiRH_half`.
-/

open Complex Asymptotics Filter Topology MeasureTheory Set

namespace Fayed.Riesz

/-- n-th term of the Möbius form of the Riesz function. -/
noncomputable def term (n : ℕ) (x : ℝ) : ℝ :=
  (ArithmeticFunction.moebius n : ℝ) * (x / (n : ℝ) ^ 2) * Real.exp (-(x / (n : ℝ) ^ 2))

/-- Riesz function, Möbius form: Σ μ(n) (x/n²) e^{-x/n²}. -/
noncomputable def rieszM (x : ℝ) : ℝ :=
  ∑' n : ℕ, (ArithmeticFunction.moebius n : ℝ) * (x / (n : ℝ) ^ 2) * Real.exp (-(x / (n : ℝ) ^ 2))

lemma rieszM_eq (x : ℝ) : rieszM x = ∑' n, term n x := rfl

/-- Complexified Riesz function. -/
noncomputable def F (x : ℝ) : ℂ := (rieszM x : ℂ)

/-- M(s) = ∫_0^∞ rieszM(x) x^{-s-1} dx. -/
noncomputable def M (s : ℂ) : ℂ := mellin F (-s)

/-- Dirichlet coefficients μ(n)/n². -/
noncomputable def a (n : ℕ) : ℂ := (ArithmeticFunction.moebius n : ℂ) / (n : ℂ) ^ 2

/-- Exponential rates 1/n². -/
noncomputable def p (n : ℕ) : ℝ := ((n : ℝ) ^ 2)⁻¹

/-- G(t) = Σ (μ(n)/n²) e^{-t/n²}, so rieszM t = t G(t). -/
noncomputable def G (t : ℝ) : ℂ := ∑' n, a n * (Real.exp (-p n * t) : ℂ)

lemma abs_mu_le (n : ℕ) : |(ArithmeticFunction.moebius n : ℝ)| ≤ 1 := by
  exact_mod_cast ArithmeticFunction.abs_moebius_le_one

lemma abs_term_le (n : ℕ) {x : ℝ} (hx : 0 ≤ x) : |term n x| ≤ x * (1 / (n : ℝ) ^ 2) := by
  have hq : 0 ≤ x / (n : ℝ) ^ 2 := by positivity
  have he : Real.exp (-(x / (n : ℝ) ^ 2)) ≤ 1 := Real.exp_le_one_iff.mpr (by linarith)
  rw [term, abs_mul, abs_mul, abs_of_nonneg hq, abs_of_pos (Real.exp_pos _)]
  calc |(ArithmeticFunction.moebius n : ℝ)| * (x / (n : ℝ) ^ 2) * Real.exp (-(x / (n : ℝ) ^ 2))
      ≤ 1 * (x / (n : ℝ) ^ 2) * 1 := by gcongr; exact abs_mu_le n
    _ = x * (1 / (n : ℝ) ^ 2) := by ring

lemma summable_inv_sq : Summable (fun n : ℕ => 1 / (n : ℝ) ^ 2) :=
  Real.summable_one_div_nat_pow.mpr one_lt_two

lemma abs_rieszM_le {x : ℝ} (hx : 0 ≤ x) :
    |rieszM x| ≤ x * ∑' n : ℕ, 1 / (n : ℝ) ^ 2 := by
  rw [← Real.norm_eq_abs, rieszM_eq]
  exact tsum_of_norm_bounded (summable_inv_sq.hasSum.mul_left x)
    (fun n => by rw [Real.norm_eq_abs]; exact abs_term_le n hx)

lemma continuousOn_rieszM (R : ℝ) : ContinuousOn rieszM (Icc 0 R) := by
  have : rieszM = fun x => ∑' n, term n x := rfl
  rw [this]
  refine continuousOn_tsum (u := fun n : ℕ => R * (1 / (n : ℝ) ^ 2)) (fun n => ?_)
    (summable_inv_sq.mul_left R) (fun n x hx => ?_)
  · have : Continuous (term n) := by unfold term; fun_prop
    exact this.continuousOn
  · rw [Real.norm_eq_abs]
    refine (abs_term_le n hx.1).trans ?_
    have : 0 ≤ 1 / (n : ℝ) ^ 2 := by positivity
    exact mul_le_mul_of_nonneg_right hx.2 this

lemma locallyIntegrableOn_F : LocallyIntegrableOn F (Ioi 0) := by
  refine ContinuousOn.locallyIntegrableOn (fun x hx => ?_) measurableSet_Ioi
  have hx : (0 : ℝ) < x := hx
  have hc : ContinuousAt rieszM x :=
    (continuousOn_rieszM (x + 1)).continuousAt (Icc_mem_nhds hx (by linarith))
  exact (continuous_ofReal.continuousAt.comp hc).continuousWithinAt

lemma isBigO_F_atTop {θ : ℝ} (h : ∃ C : ℝ, ∀ x : ℝ, 1 ≤ x → |rieszM x| ≤ C * x ^ θ) :
    F =O[atTop] (fun t : ℝ => t ^ (-(-θ))) := by
  obtain ⟨C, hC⟩ := h
  rw [neg_neg]
  refine IsBigO.of_bound C ?_
  filter_upwards [eventually_ge_atTop 1] with x hx
  rw [F, Complex.norm_real, Real.norm_eq_abs,
    Real.norm_of_nonneg (Real.rpow_nonneg (by linarith) _)]
  exact hC x hx

lemma isBigO_F_zero : F =O[𝓝[>] 0] (fun t : ℝ => t ^ (-(-1 : ℝ))) := by
  refine IsBigO.of_bound (∑' n : ℕ, 1 / (n : ℝ) ^ 2) ?_
  filter_upwards [self_mem_nhdsWithin] with x (hx : 0 < x)
  rw [neg_neg, Real.rpow_one, F, Complex.norm_real, Real.norm_eq_abs,
    Real.norm_of_nonneg hx.le, mul_comm]
  exact abs_rieszM_le hx.le

lemma differentiableAt_M
    (h : ∀ ε : ℝ, 0 < ε → ∃ C : ℝ, ∀ x : ℝ, 1 ≤ x → |rieszM x| ≤ C * x ^ (1 / 4 + ε))
    {s : ℂ} (h1 : 1 / 4 < s.re) (h2 : s.re < 1) : DifferentiableAt ℂ M s := by
  have hε : 0 < (s.re - 1 / 4) / 2 := by linarith
  have hm : DifferentiableAt ℂ (mellin F) (-s) :=
    mellin_differentiableAt_of_isBigO_rpow (a := -(1 / 4 + (s.re - 1 / 4) / 2)) (b := -1)
      locallyIntegrableOn_F (isBigO_F_atTop (h _ hε)) (by rw [neg_re]; linarith)
      isBigO_F_zero (by rw [neg_re]; linarith)
  exact hm.comp s differentiableAt_id.neg

lemma norm_a_le (n : ℕ) : ‖a n‖ ≤ 1 / (n : ℝ) ^ 2 := by
  rw [a, norm_div, norm_pow, Complex.norm_natCast, Complex.norm_intCast]
  exact div_le_div_of_nonneg_right (abs_mu_le n) (by positivity)

lemma hasSum_G {t : ℝ} (ht : 0 ≤ t) :
    HasSum (fun n => a n * (Real.exp (-p n * t) : ℂ)) (G t) := by
  refine (Summable.of_norm_bounded summable_inv_sq fun n => ?_).hasSum
  rw [norm_mul, Complex.norm_real, Real.norm_of_nonneg (Real.exp_pos _).le]
  have hp : 0 ≤ p n := by rw [p]; positivity
  have he : Real.exp (-p n * t) ≤ 1 := Real.exp_le_one_iff.mpr (by nlinarith)
  exact (mul_le_of_le_one_right (norm_nonneg _) he).trans (norm_a_le n)

lemma F_eq : F = fun t : ℝ => (t : ℂ) ^ (1 : ℂ) • G t := by
  funext t
  rw [cpow_one, smul_eq_mul, G, ← tsum_mul_left, F, rieszM, Complex.ofReal_tsum]
  congr 1
  funext n
  rw [show -p n * t = -(t / (n : ℝ) ^ 2) by rw [p]; ring]
  simp only [a]
  push_cast
  ring

lemma summable_coeff {s : ℂ} (h1 : 1 / 2 < s.re) :
    Summable fun n => ‖a n‖ / p n ^ (-s + 1).re := by
  have hre : (-s + 1).re = 1 - s.re := by simp; ring
  rw [hre]
  refine Summable.of_nonneg_of_le
    (fun n => div_nonneg (norm_nonneg _) (Real.rpow_nonneg (by rw [p]; positivity) _)) (fun n => ?_)
    (Real.summable_one_div_nat_rpow.mpr (show (1 : ℝ) < 2 * s.re by linarith))
  rcases eq_or_ne n 0 with rfl | hn
  · simp [a]; positivity
  have hn' : (0 : ℝ) < n := by exact_mod_cast Nat.pos_of_ne_zero hn
  rw [p, Real.inv_rpow (by positivity), div_inv_eq_mul]
  have h2 : ((n : ℝ) ^ 2) ^ (1 - s.re) = (n : ℝ) ^ (2 : ℝ) / (n : ℝ) ^ (2 * s.re) := by
    rw [← Real.rpow_natCast, ← Real.rpow_mul hn'.le, ← Real.rpow_sub hn']
    norm_num
    ring_nf
  rw [h2, Real.rpow_two]
  calc ‖a n‖ * ((n : ℝ) ^ 2 / (n : ℝ) ^ (2 * s.re))
      ≤ 1 / (n : ℝ) ^ 2 * ((n : ℝ) ^ 2 / (n : ℝ) ^ (2 * s.re)) := by
        gcongr; exact norm_a_le n
    _ = 1 / (n : ℝ) ^ (2 * s.re) := by
        have : (n : ℝ) ^ (2 * s.re) ≠ 0 := (Real.rpow_pos_of_pos hn' _).ne'
        field_simp

/-- Mellin transform of the Riesz function: Γ(1-s) L(μ, 2s) = Γ(1-s)/ζ(2s). -/
lemma mellin_eq {s : ℂ} (h1 : 1 / 2 < s.re) (h2 : s.re < 1) :
    M s = Gamma (1 - s) * LSeries (fun n => (ArithmeticFunction.moebius n : ℂ)) (2 * s) := by
  have hs' : 0 < (-s + 1).re := by simp; linarith
  have hp : ∀ n, a n = 0 ∨ 0 < p n := by
    intro n
    rcases eq_or_ne n 0 with rfl | hn
    · left; simp [a]
    · right; rw [p]; have : (0 : ℝ) < n := by exact_mod_cast Nat.pos_of_ne_zero hn
      positivity
  have key := hasSum_mellin hp hs' (fun t ht => hasSum_G (le_of_lt ht)) (summable_coeff h1)
  rw [M, F_eq, mellin_cpow_smul, ← key.tsum_eq, LSeries, ← tsum_mul_left]
  congr 1
  funext n
  rcases eq_or_ne n 0 with rfl | hn
  · simp [a, LSeries.term_zero]
  rw [LSeries.term_of_ne_zero hn]
  have hn' : (n : ℂ) ≠ 0 := by exact_mod_cast hn
  have hy : (n : ℂ) ^ s ≠ 0 := (cpow_ne_zero_iff_of_exponent_ne_zero (by
    intro h0; rw [h0] at h1; simp at h1; linarith)).mpr hn'
  have hpn : ((p n : ℝ) : ℂ) = ((n : ℂ) * (n : ℂ))⁻¹ := by
    rw [p]; push_cast; ring
  have harg : ((n : ℂ) * (n : ℂ)).arg ≠ Real.pi := by
    rw [← Nat.cast_mul, natCast_arg]; exact Real.pi_pos.ne
  rw [hpn, inv_cpow _ _ harg, natCast_mul_natCast_cpow, cpow_add _ _ hn', cpow_one, cpow_neg,
    cpow_ofNat_mul, show (1 : ℂ) - s = -s + 1 by ring, a]
  field_simp

lemma isLinearMap_im_sub_re : IsLinearMap ℝ (fun z : ℂ => z.im - z.re) :=
  ⟨fun x y => by simp; ring, fun c x => by simp; ring⟩

lemma isLinearMap_re_add_im : IsLinearMap ℝ (fun z : ℂ => z.re + z.im) :=
  ⟨fun x y => by simp; ring, fun c x => by simp; ring⟩

/-- Identity theorem: ζ(2z) M(z) = Γ(1-z) on any convex open D inside the strip, off 1/2,
meeting Re z > 1/2. -/
lemma eq_on_convex
    (h : ∀ ε : ℝ, 0 < ε → ∃ C : ℝ, ∀ x : ℝ, 1 ≤ x → |rieszM x| ≤ C * x ^ (1 / 4 + ε))
    {D : Set ℂ} (hDo : IsOpen D) (hDc : Convex ℝ D)
    (hD : ∀ z ∈ D, 1 / 4 < z.re ∧ z.re < 1 ∧ z ≠ 1 / 2)
    {c : ℂ} (hcD : c ∈ D) (hc : 1 / 2 < c.re) {s : ℂ} (hs : s ∈ D) :
    riemannZeta (2 * s) * M s = Gamma (1 - s) := by
  have hd : DifferentiableOn ℂ (fun z => riemannZeta (2 * z) * M z) D := by
    intro z hz
    obtain ⟨h1, h2, h3⟩ := hD z hz
    have h2z : 2 * z ≠ 1 := fun h0 => h3 (by linear_combination h0 / 2)
    exact (((differentiableAt_riemannZeta h2z).comp z
      (differentiableAt_id.const_mul 2)).mul (differentiableAt_M h h1 h2)).differentiableWithinAt
  have hg : DifferentiableOn ℂ (fun z => Gamma (1 - z)) D := by
    intro z hz
    obtain ⟨_, h2, _⟩ := hD z hz
    have hG : ∀ m : ℕ, 1 - z ≠ -m := by
      intro m hm
      have := congrArg Complex.re hm
      simp at this
      have : (0 : ℝ) ≤ m := Nat.cast_nonneg m
      linarith
    exact ((differentiableAt_Gamma _ hG).comp z
      (differentiableAt_id.const_sub 1)).differentiableWithinAt
  have hev : (fun z => riemannZeta (2 * z) * M z) =ᶠ[𝓝 c] (fun z => Gamma (1 - z)) := by
    have hmem : {z : ℂ | 1 / 2 < z.re} ∩ D ∈ 𝓝 c :=
      ((isOpen_lt continuous_const Complex.continuous_re).inter hDo).mem_nhds ⟨hc, hcD⟩
    filter_upwards [hmem] with z hz
    obtain ⟨hz1, hzD⟩ := hz
    have hz1 : 1 / 2 < z.re := hz1
    have hz2 : 1 < (2 * z).re := by
      rw [show (2 * z).re = 2 * z.re by simp]; linarith
    have hL := ArithmeticFunction.LSeries_zeta_mul_Lseries_moebius hz2
    rw [mellin_eq hz1 (hD z hzD).2.1, ← ArithmeticFunction.LSeries_zeta_eq_riemannZeta hz2]
    linear_combination (Gamma (1 - z)) * hL
  exact (hd.analyticOnNhd hDo).eqOn_of_preconnected_of_eventuallyEq (hg.analyticOnNhd hDo)
    hDc.isPreconnected hcD hev hs

lemma zeta_two_mul_ne_zero
    (h : ∀ ε : ℝ, 0 < ε → ∃ C : ℝ, ∀ x : ℝ, 1 ≤ x → |rieszM x| ≤ C * x ^ (1 / 4 + ε))
    {s : ℂ} (h1 : 1 / 4 < s.re) (h2 : s.re < 1 / 2) : riemannZeta (2 * s) ≠ 0 := by
  have hS : IsOpen {z : ℂ | 1 / 4 < z.re ∧ z.re < 1} :=
    (isOpen_lt continuous_const Complex.continuous_re).inter
      (isOpen_lt Complex.continuous_re continuous_const)
  have hSc : Convex ℝ {z : ℂ | 1 / 4 < z.re ∧ z.re < 1} :=
    (convex_halfSpace_re_gt _).inter (convex_halfSpace_re_lt _)
  have key : riemannZeta (2 * s) * M s = Gamma (1 - s) := by
    by_cases him : -(1 / 2) < s.im - s.re
    · refine eq_on_convex h (D := {z : ℂ | 1 / 4 < z.re ∧ z.re < 1} ∩
          {z : ℂ | -(1 / 2) < z.im - z.re})
        (hS.inter (isOpen_lt continuous_const
          (Complex.continuous_im.sub Complex.continuous_re)))
        (hSc.inter (convex_halfSpace_gt isLinearMap_im_sub_re _)) ?_
        (c := 3 / 4 + I) ⟨⟨by simp; norm_num, by simp; norm_num⟩, by simp; norm_num⟩
        (by simp; norm_num) ⟨⟨h1, by linarith⟩, him⟩
      rintro z ⟨⟨hz1, hz2⟩, hz3⟩
      refine ⟨hz1, hz2, fun hz => ?_⟩
      subst hz; norm_num at hz3
    · push Not at him
      refine eq_on_convex h (D := {z : ℂ | 1 / 4 < z.re ∧ z.re < 1} ∩
          {z : ℂ | z.re + z.im < 1 / 2})
        (hS.inter (isOpen_lt (Complex.continuous_re.add Complex.continuous_im)
          continuous_const))
        (hSc.inter (convex_halfSpace_lt isLinearMap_re_add_im _)) ?_
        (c := 3 / 4 - I) ⟨⟨by simp; norm_num, by simp; norm_num⟩, by simp; norm_num⟩
        (by simp; norm_num) ⟨⟨h1, by linarith⟩, show s.re + s.im < 1 / 2 by linarith⟩
      rintro z ⟨⟨hz1, hz2⟩, hz3⟩
      refine ⟨hz1, hz2, fun hz => ?_⟩
      subst hz; norm_num at hz3
  intro h0
  rw [h0, zero_mul] at key
  exact Gamma_ne_zero_of_re_pos (by simp; linarith) key.symm

/-- **Riesz's criterion (sufficiency).** If rieszM(x) = O(x^{1/4+ε}) for every ε > 0, then RH. -/
theorem RH_of_rieszBound
    (h : ∀ ε : ℝ, 0 < ε → ∃ C : ℝ, ∀ x : ℝ, 1 ≤ x → |rieszM x| ≤ C * x ^ (1 / 4 + ε)) :
    RiemannHypothesis := by
  rw [Fayed.R7.RH_iff_QuasiRH_half]
  intro w hw _
  rcases le_or_gt 1 w.re with h1 | h1
  · exact riemannZeta_ne_zero_of_one_le_re h1
  · have := zeta_two_mul_ne_zero h (s := w / 2) (by simp; linarith) (by simp; linarith)
    rwa [mul_div_cancel₀ w two_ne_zero] at this

end Fayed.Riesz

#print axioms Fayed.Riesz.RH_of_rieszBound
