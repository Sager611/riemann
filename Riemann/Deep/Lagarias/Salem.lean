import Riemann.Deep.Mertens.ZeroFree

/-!
# Salem's criterion, sufficiency direction

Salem (1953): RH ⟺ for every 1/2 < σ < 1 the integral equation
  ∫_ℝ e^{-σy} φ(y) / (exp(e^{x-y}) + 1) dy = 0   (all x)
has no bounded solution φ ≠ 0.  Direction formalized here: a zero ρ of ζ with
1/2 < Re ρ < 1 gives the bounded nonzero solution φ(y) = exp(-i·Im ρ·y).

Route:
* `K t = 1/(e^t+1)`; `mellin K` is holomorphic on Re s > 0 (exp decay, bounded at 0).
* Re s > 1: `mellin K s = Γ(s) η(s)`, η(s) = Σ (-1)^n/(n+1)^s = (1 - 2^{1-s}) ζ(s)
  (geometric series + `hasSum_mellin`, even/odd split of ζ's series).
* Identity theorem on {Re s > 0} \ {1} (`Fayed.Mertens.isPreconnected_U`): same for Re s > 0.
* Substitution t = e^{x-y}: the Salem integral is e^{-ρx} · mellin K ρ = 0.
-/

open Complex Asymptotics Filter Topology MeasureTheory Set

namespace Fayed.Lagarias.Salem

/-- Fermi–Dirac kernel K(t) = 1/(e^t + 1), complex-valued. -/
noncomputable def K (t : ℝ) : ℂ := ((1 / (Real.exp t + 1) : ℝ) : ℂ)

/-- Dirichlet η, written via ζ: η(s) = (1 - 2·2^{-s}) ζ(s). -/
noncomputable def eta (s : ℂ) : ℂ := (1 - 2 * (2 : ℂ) ^ (-s)) * riemannZeta s

lemma eta_eq (s : ℂ) : eta s = (1 - (2 : ℂ) ^ (1 - s)) * riemannZeta s := by
  rw [eta, sub_eq_add_neg (1 : ℂ) s, cpow_add _ _ two_ne_zero, cpow_one]

lemma K_continuous : Continuous K :=
  Complex.continuous_ofReal.comp (continuous_const.div (Real.continuous_exp.add continuous_const)
    (fun t => by positivity))

lemma norm_K (t : ℝ) : ‖K t‖ = 1 / (Real.exp t + 1) := by
  rw [K, Complex.norm_real, Real.norm_of_nonneg (by positivity)]

lemma K_locInt : LocallyIntegrableOn K (Ioi 0) :=
  K_continuous.continuousOn.locallyIntegrableOn measurableSet_Ioi

lemma K_top : K =O[atTop] fun t => Real.exp (-1 * t) := by
  refine IsBigO.of_bound 1 (Eventually.of_forall fun t => ?_)
  rw [norm_K, one_mul, Real.norm_of_nonneg (Real.exp_pos _).le, neg_one_mul, Real.exp_neg,
    ← one_div]
  exact one_div_le_one_div_of_le (Real.exp_pos t) (le_add_of_nonneg_right zero_le_one)

lemma K_bot : K =O[𝓝[>] 0] fun t : ℝ => t ^ (-(0 : ℝ)) := by
  refine IsBigO.of_bound 1 (Eventually.of_forall fun t => ?_)
  rw [norm_K, neg_zero, Real.rpow_zero, norm_one, one_mul]
  rw [div_le_one (by positivity)]
  linarith [Real.exp_pos t]

lemma mellinConvergent_K {s : ℂ} (hs : 0 < s.re) : MellinConvergent K s :=
  mellinConvergent_of_isBigO_rpow (a := s.re + 1) K_locInt
    (K_top.trans (isLittleO_exp_neg_mul_rpow_atTop one_pos _).isBigO) (by linarith) K_bot hs

lemma differentiableAt_mellin_K {s : ℂ} (hs : 0 < s.re) : DifferentiableAt ℂ (mellin K) s :=
  mellin_differentiableAt_of_isBigO_rpow_exp one_pos K_locInt K_top K_bot hs

/-- 1/(e^t+1) = Σ_{n≥0} (-1)^n e^{-(n+1)t}, t > 0. -/
lemma hasSum_K_real {t : ℝ} (ht : 0 < t) :
    HasSum (fun n : ℕ => (-1 : ℝ) ^ n * Real.exp (-((n : ℝ) + 1) * t)) (1 / (Real.exp t + 1)) := by
  have hq : ‖-Real.exp (-t)‖ < 1 := by
    rw [norm_neg, Real.norm_of_nonneg (Real.exp_pos _).le]
    exact Real.exp_lt_one_iff.mpr (by linarith)
  have h := (hasSum_geometric_of_norm_lt_one hq).mul_left (Real.exp (-t))
  convert h using 1
  · funext n
    rw [neg_pow (Real.exp (-t)) n, ← Real.exp_nat_mul, show -((n : ℝ) + 1) * t = -t + n * (-t) by ring,
      Real.exp_add]
    ring
  · rw [Real.exp_neg, sub_neg_eq_add]
    have h1 : Real.exp t ≠ 0 := (Real.exp_pos t).ne'
    have h2 : Real.exp t + 1 ≠ 0 := by positivity
    have h3 : 1 + (Real.exp t)⁻¹ ≠ 0 := by positivity
    rw [eq_mul_inv_iff_mul_eq₀ h3]
    field_simp

/-- ζ's terms 1/(n+1)^s. -/
noncomputable def zt (s : ℂ) (n : ℕ) : ℂ := 1 / ((n : ℂ) + 1) ^ s

lemma hasSum_zt {s : ℂ} (hs : 1 < s.re) : HasSum (zt s) (riemannZeta s) := by
  have hsum : Summable (zt s) := by
    have := (summable_nat_add_iff 1).mpr (Complex.summable_one_div_nat_cpow.mpr hs)
    show Summable fun n : ℕ => 1 / ((n : ℂ) + 1) ^ s
    simpa [Nat.cast_add, Nat.cast_one] using this
  rw [zeta_eq_tsum_one_div_nat_add_one_cpow hs]
  exact hsum.hasSum

lemma zt_odd (s : ℂ) (k : ℕ) : zt s (2 * k + 1) = (2 : ℂ) ^ (-s) * zt s k := by
  have h2 := natCast_mul_natCast_cpow 2 (k + 1) s
  push_cast at h2
  have e : ((2 * k + 1 : ℕ) : ℂ) + 1 = 2 * ((k : ℂ) + 1) := by push_cast; ring
  rw [zt, zt, e, h2, cpow_neg]
  have h2s : (2 : ℂ) ^ s ≠ 0 := cpow_ne_zero_iff.mpr (Or.inl two_ne_zero)
  field_simp

/-- Dirichlet series of η for Re s > 1. -/
lemma hasSum_eta {s : ℂ} (hs : 1 < s.re) :
    HasSum (fun n : ℕ => (-1 : ℂ) ^ n / ((n : ℂ) + 1) ^ s) (eta s) := by
  have hf := hasSum_zt hs
  have hE : Summable fun k : ℕ => zt s (2 * k) :=
    hf.summable.comp_injective (mul_right_injective₀ (two_ne_zero' ℕ))
  have hO : HasSum (fun k : ℕ => zt s (2 * k + 1)) ((2 : ℂ) ^ (-s) * riemannZeta s) := by
    simpa only [zt_odd] using hf.mul_left ((2 : ℂ) ^ (-s))
  have hEv : HasSum (fun k : ℕ => zt s (2 * k)) (riemannZeta s - (2 : ℂ) ^ (-s) * riemannZeta s) := by
    have := tsum_even_add_odd hE hO.summable
    rw [hf.tsum_eq, hO.tsum_eq] at this
    rw [← eq_sub_of_add_eq this]
    exact hE.hasSum
  have key : HasSum (fun n : ℕ => (-1 : ℂ) ^ n / ((n : ℂ) + 1) ^ s)
      ((riemannZeta s - (2 : ℂ) ^ (-s) * riemannZeta s) + -((2 : ℂ) ^ (-s) * riemannZeta s)) := by
    refine HasSum.even_add_odd ?_ ?_
    · convert hEv using 1; funext k; simp [zt, pow_mul]
    · convert hO.neg using 1; funext k; simp [zt, pow_succ, pow_mul, neg_div]
  convert key using 1
  rw [eta]; ring

/-- Re s > 1: mellin K s = Γ(s) η(s). -/
lemma mellin_K_of_one_lt {s : ℂ} (hs : 1 < s.re) : mellin K s = Gamma s * eta s := by
  have hm := hasSum_mellin (a := fun n : ℕ => (-1 : ℂ) ^ n) (p := fun n : ℕ => (n : ℝ) + 1)
    (F := K) (s := s) (fun n => Or.inr (by positivity)) (by linarith)
    (fun t ht => by
      convert (Complex.hasSum_ofReal.mpr (hasSum_K_real ht)) using 1
      · funext n; push_cast; ring
      · rfl)
    (by
      refine ((Real.summable_one_div_nat_add_rpow 1 s.re).mpr hs).congr fun n => ?_
      rw [norm_pow, norm_neg, norm_one, one_pow, abs_of_pos (by positivity)])
  refine hm.unique ?_
  convert (hasSum_eta hs).mul_left (Gamma s) using 1
  funext n; push_cast; ring

/-- Re s > 0, s ≠ 1: mellin K s = Γ(s) η(s), by the identity theorem. -/
lemma mellin_K_eq {s : ℂ} (hs : 0 < s.re) (hs1 : s ≠ 1) : mellin K s = Gamma s * eta s := by
  let U : Set ℂ := {s : ℂ | 0 < s.re} \ {1}
  have hUo : IsOpen U := (isOpen_lt continuous_const Complex.continuous_re).sdiff isClosed_singleton
  have hF : DifferentiableOn ℂ (mellin K) U :=
    fun z hz => (differentiableAt_mellin_K hz.1).differentiableWithinAt
  have hG : DifferentiableOn ℂ (fun z => Gamma z * eta z) U := by
    intro z hz
    have hz0 : 0 < z.re := hz.1
    have hz1 : z ≠ 1 := hz.2
    have hGa : DifferentiableAt ℂ Gamma z := differentiableAt_Gamma z fun m hm => by
      have := congrArg Complex.re hm
      simp at this
      linarith [(Nat.cast_nonneg m : (0 : ℝ) ≤ m)]
    have hc : DifferentiableAt ℂ (fun w : ℂ => (2 : ℂ) ^ (-w)) z :=
      differentiableAt_id.neg.const_cpow (Or.inl two_ne_zero)
    have he : DifferentiableAt ℂ eta z :=
      ((differentiableAt_const _).sub ((differentiableAt_const _).mul hc)).mul
        (differentiableAt_riemannZeta hz1)
    exact (hGa.mul he).differentiableWithinAt
  have hc : (2 : ℂ) ∈ U := ⟨show (0 : ℝ) < (2 : ℂ).re by simp, by norm_num⟩
  have hev : mellin K =ᶠ[𝓝 (2 : ℂ)] fun z => Gamma z * eta z := by
    have hmem : {z : ℂ | 1 < z.re} ∈ 𝓝 (2 : ℂ) :=
      (isOpen_lt continuous_const Complex.continuous_re).mem_nhds (show (1 : ℝ) < (2 : ℂ).re by simp)
    filter_upwards [hmem] with z hz
    exact mellin_K_of_one_lt hz
  exact (hF.analyticOnNhd hUo).eqOn_of_preconnected_of_eventuallyEq (hG.analyticOnNhd hUo)
    (Fayed.Mertens.isPreconnected_U 0) hc hev ⟨hs, hs1⟩

lemma mellin_K_eq_zero {ρ : ℂ} (hz : riemannZeta ρ = 0) (h0 : 0 < ρ.re) (h1 : ρ ≠ 1) :
    mellin K ρ = 0 := by
  rw [mellin_K_eq h0 h1, eta, hz, mul_zero, mul_zero]

lemma exp_smul_integrand (ρ : ℂ) (u : ℝ) :
    Real.exp u • ((Real.exp u : ℂ) ^ (ρ - 1) • K (Real.exp u)) =
      cexp (ρ * u) * K (Real.exp u) := by
  rw [smul_eq_mul, Complex.real_smul, cpow_def_of_ne_zero (ofReal_ne_zero.mpr (Real.exp_pos u).ne'),
    ← ofReal_log (Real.exp_pos u).le, Real.log_exp, ofReal_exp, ← mul_assoc, ← Complex.exp_add]
  congr 2
  ring

/-- Substitution t = e^{x-y}: ∫ e^{-ρy} K(e^{x-y}) dy = e^{-ρx} · mellin K ρ (Re ρ > 0). -/
lemma salem_integral (x : ℝ) {ρ : ℂ} (hρ : 0 < ρ.re) :
    Integrable (fun y : ℝ => cexp (-ρ * y) * K (Real.exp (x - y))) ∧
      ∫ y : ℝ, cexp (-ρ * y) * K (Real.exp (x - y)) = cexp (-ρ * x) * mellin K ρ := by
  set h : ℝ → ℂ := fun u => cexp (-ρ * x) * (cexp (ρ * u) * K (Real.exp u)) with hh
  have hint0 : Integrable fun u : ℝ => cexp (ρ * u) * K (Real.exp u) := by
    have := (integrable_comp_exp (fun t : ℝ => (t : ℂ) ^ (ρ - 1) • K t)).mpr (mellinConvergent_K hρ)
    simpa only [exp_smul_integrand] using this
  have hint : Integrable h := hint0.const_mul _
  have hfun : (fun y : ℝ => cexp (-ρ * y) * K (Real.exp (x - y))) = fun y => h (x - y) := by
    funext y
    simp only [hh]
    rw [← mul_assoc, ← Complex.exp_add]
    congr 2
    push_cast; ring
  rw [hfun]
  refine ⟨hint.comp_sub_left x, ?_⟩
  rw [integral_sub_left_eq_self h volume x, hh, integral_const_mul]
  congr 1
  rw [mellin, ← integral_comp_exp]
  simp only [exp_smul_integrand]

/-- φ(y) = exp(-i τ y). -/
noncomputable def phi (τ : ℝ) (y : ℝ) : ℂ := cexp (-(I * τ * y))

lemma norm_phi (τ y : ℝ) : ‖phi τ y‖ = 1 := by
  rw [phi, Complex.norm_exp]
  simp

/-- Salem, sufficiency direction: a zero ρ of ζ with 1/2 < Re ρ < 1 yields a bounded
(|φ| = 1) nonzero solution of ∫ e^{-σy} φ(y) K(e^{x-y}) dy = 0, σ = Re ρ, φ = exp(-i·Im ρ·y),
where K(t) = 1/(e^t+1). Hence "no bounded nonzero solution for all σ ∈ (1/2,1)" ⇒ RH. -/
theorem salem_solution {ρ : ℂ} (hz : riemannZeta ρ = 0) (h1 : 1 / 2 < ρ.re) (h2 : ρ.re < 1) :
    (∀ y, ‖phi ρ.im y‖ = 1) ∧ ∀ x : ℝ,
      Integrable (fun y : ℝ => (Real.exp (-ρ.re * y) : ℂ) * phi ρ.im y * K (Real.exp (x - y))) ∧
      ∫ y : ℝ, (Real.exp (-ρ.re * y) : ℂ) * phi ρ.im y * K (Real.exp (x - y)) = 0 := by
  have h0 : 0 < ρ.re := by linarith
  have hne : ρ ≠ 1 := fun h => by rw [h, one_re] at h2; exact lt_irrefl _ h2
  refine ⟨norm_phi ρ.im, fun x => ?_⟩
  have e : (fun y : ℝ => (Real.exp (-ρ.re * y) : ℂ) * phi ρ.im y * K (Real.exp (x - y))) =
      fun y : ℝ => cexp (-ρ * y) * K (Real.exp (x - y)) := by
    funext y
    rw [phi, ofReal_exp, ← Complex.exp_add]
    congr 2
    conv_rhs => rw [← re_add_im ρ]
    push_cast; ring
  obtain ⟨hi, hv⟩ := salem_integral x h0
  rw [e, hv, mellin_K_eq_zero hz h0 hne, mul_zero]
  exact ⟨hi, rfl⟩

end Fayed.Lagarias.Salem

#print axioms Fayed.Lagarias.Salem.salem_solution
