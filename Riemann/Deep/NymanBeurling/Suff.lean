import Riemann.Tree.R007
import Riemann.Deep.Mertens.ZeroFree

/-!
# Nyman-Beurling criterion, sufficiency half

If 1 is an L^2(0,1)-limit of finite combinations of rho_θ(x) = {θ/x} - θ{1/x} (0 < θ ≤ 1),
then RH holds.

Route: F(s) = ∫_1^∞ {t} t^{-s-1} dt (Mellin transform, holomorphic on Re s > 0);
ζ(s) = s/(s-1) - s F(s) for Re s > 1 (Abel summation, LSeries_eq_mul_integral'), extended to
{Re s > 0} \ {1} by the identity theorem. For 0 < Re s < 1,
M(s) := ∫_0^∞ {1/x} x^{s-1} dx = 1/(1-s) + F(s), so ζ(s) = 0 ⇒ M(s) = 0, hence every
∫_0^1 rho_θ x^{s-1} = θ^s M(s) - θ M(s) = 0. With ∫_0^1 x^{s-1} = 1/s, Cauchy-Schwarz (AM-GM form)
contradicts ‖1 - f‖₂ → 0 when Re s > 1/2.
-/

open Complex Asymptotics Filter Topology MeasureTheory Set

namespace Fayed.NB

noncomputable def rhoNB (θ x : ℝ) : ℝ := Int.fract (θ / x) - θ * Int.fract (1 / x)

def NBApprox : Prop := ∀ ε > 0, ∃ (n : ℕ) (c θ : Fin n → ℝ), (∀ k, 0 < θ k ∧ θ k ≤ 1) ∧
    ∫ x in Set.Ioc (0:ℝ) 1, (1 - ∑ k, c k * rhoNB (θ k) x) ^ 2 < ε

/-! ## A. F(s) = ∫_1^∞ {t} t^{-s-1} dt -/

/-- {t} on [1, ∞), zero below. -/
noncomputable def Q (t : ℝ) : ℂ := if 1 ≤ t then ((Int.fract t : ℝ) : ℂ) else 0

noncomputable def F (s : ℂ) : ℂ := mellin Q (-s)

lemma norm_fract_le (t : ℝ) : ‖((Int.fract t : ℝ) : ℂ)‖ ≤ 1 := by
  rw [Complex.norm_real, Real.norm_of_nonneg (Int.fract_nonneg t)]
  exact (Int.fract_lt_one t).le

lemma norm_Q_le (t : ℝ) : ‖Q t‖ ≤ 1 := by
  unfold Q; split_ifs
  · exact norm_fract_le t
  · simp

lemma Q_eq_zero {t : ℝ} (ht : t < 1) : Q t = 0 := by simp [Q, not_le.mpr ht]

lemma measurable_Q : Measurable Q :=
  Measurable.ite measurableSet_Ici (Complex.measurable_ofReal.comp measurable_fract)
    measurable_const

lemma locallyIntegrableOn_Q : LocallyIntegrableOn Q (Ioi 0) := by
  refine ((continuous_const (y := (1 : ℝ))).locallyIntegrable.locallyIntegrableOn (Ioi 0)).mono
    measurable_Q.aestronglyMeasurable (Eventually.of_forall fun t => ?_)
  simpa using norm_Q_le t

lemma isBigO_Q_atTop : Q =O[atTop] (fun t : ℝ => t ^ (-(0:ℝ))) := by
  refine IsBigO.of_bound 1 (Eventually.of_forall fun t => ?_)
  simpa using norm_Q_le t

lemma isBigO_Q_zero (b : ℝ) : Q =O[𝓝[>] 0] (fun t : ℝ => t ^ (-b)) := by
  have : Q =ᶠ[𝓝[>] (0 : ℝ)] (fun _ => (0 : ℂ)) := by
    filter_upwards [Ioo_mem_nhdsGT (zero_lt_one' ℝ)] with t ht using Q_eq_zero ht.2
  exact (isBigO_zero _ _).congr' this.symm EventuallyEq.rfl

lemma mellinConvergent_Q {s : ℂ} (hs : 0 < s.re) : MellinConvergent Q (-s) :=
  mellinConvergent_of_isBigO_rpow (a := 0) (b := -s.re - 1) locallyIntegrableOn_Q isBigO_Q_atTop
    (by rw [neg_re]; linarith) (isBigO_Q_zero _) (by rw [neg_re]; linarith)

lemma differentiableAt_F {s : ℂ} (hs : 0 < s.re) : DifferentiableAt ℂ F s := by
  have hm : DifferentiableAt ℂ (mellin Q) (-s) :=
    mellin_differentiableAt_of_isBigO_rpow (a := 0) (b := -s.re - 1) locallyIntegrableOn_Q
      isBigO_Q_atTop (by rw [neg_re]; linarith) (isBigO_Q_zero _) (by rw [neg_re]; linarith)
  show DifferentiableAt ℂ (fun z => mellin Q (-z)) s
  exact hm.comp s differentiableAt_id.neg

lemma F_eq_integral (s : ℂ) :
    F s = ∫ t in Ici (1 : ℝ), (t : ℂ) ^ (-s - 1) * ((Int.fract t : ℝ) : ℂ) := by
  rw [F, mellin, setIntegral_eq_of_subset_of_ae_sdiff_eq_zero measurableSet_Ioi.nullMeasurableSet
      (fun t (ht : 1 ≤ t) => (zero_lt_one.trans_le ht : (0 : ℝ) < t))
      (Eventually.of_forall fun t ht => ?_)]
  · refine setIntegral_congr_fun measurableSet_Ici fun t (ht : 1 ≤ t) => ?_
    simp [Q, ht]
  · rw [Q_eq_zero (not_le.mp ht.2), smul_zero]

lemma integrableOn_hi {s : ℂ} (hs : 0 < s.re) :
    IntegrableOn (fun t : ℝ => (t : ℂ) ^ (-s - 1) * ((Int.fract t : ℝ) : ℂ)) (Ici 1) := by
  refine (IntegrableOn.mono_set (mellinConvergent_Q hs)
    (fun t (ht : 1 ≤ t) => (zero_lt_one.trans_le ht : (0 : ℝ) < t))).congr_fun
    (fun t (ht : 1 ≤ t) => ?_) measurableSet_Ici
  simp [Q, ht]

/-! ## B. ζ(s) = s/(s-1) - s F(s) for Re s > 1 -/

lemma zeta_eq_of_one_lt {s : ℂ} (hs : 1 < s.re) :
    riemannZeta s = s / (s - 1) - s * F s := by
  have hO : (fun n : ℕ => ∑ k ∈ Finset.Icc 1 n, ‖(1 : ℕ → ℂ) k‖) =O[atTop]
      (fun n : ℕ => (n : ℝ) ^ (1 : ℝ)) := by
    refine IsBigO.of_bound 1 (Eventually.of_forall fun n => ?_)
    simp
  rw [← LSeries_one_eq_riemannZeta hs, LSeries_eq_mul_integral' _ zero_le_one (by linarith) hO]
  have h1 : ∫ t in Ioi (1:ℝ), (∑ k ∈ Finset.Icc 1 ⌊t⌋₊, (1 : ℕ → ℂ) k) * (t : ℂ) ^ (-(s + 1)) =
      ∫ t in Ioi (1:ℝ), ((t : ℂ) ^ (-s) - (t : ℂ) ^ (-s - 1) * ((Int.fract t : ℝ) : ℂ)) := by
    refine setIntegral_congr_fun measurableSet_Ioi fun t (ht : 1 < t) => ?_
    have h0 : (0:ℝ) ≤ t := zero_le_one.trans ht.le
    have ht0 : (t : ℂ) ≠ 0 := by exact_mod_cast (zero_lt_one.trans ht).ne'
    have hfl : ((⌊t⌋₊ : ℕ) : ℂ) = ((⌊t⌋ : ℤ) : ℂ) := by
      rw [← Int.natCast_floor_eq_floor h0, Int.cast_natCast]
    have hfr : ((Int.fract t : ℝ) : ℂ) = (t : ℂ) - ((⌊t⌋ : ℤ) : ℂ) := by
      rw [Int.fract]; push_cast; ring
    simp only [Pi.one_apply, Finset.sum_const, Nat.card_Icc, nsmul_eq_mul, mul_one,
      Nat.add_sub_cancel]
    rw [hfl, hfr, show -(s + 1) = -s - 1 by ring, cpow_sub _ _ ht0, cpow_one]
    field_simp
    ring
  have hi1 : IntegrableOn (fun t : ℝ => (t : ℂ) ^ (-s)) (Ioi 1) :=
    integrableOn_Ioi_cpow_of_lt (by rw [neg_re]; linarith) zero_lt_one
  have hi2 : IntegrableOn (fun t : ℝ => (t : ℂ) ^ (-s - 1) * ((Int.fract t : ℝ) : ℂ)) (Ioi 1) :=
    (integrableOn_hi (by linarith)).mono_set Ioi_subset_Ici_self
  have h3 : ∫ t in Ioi (1:ℝ), (t : ℂ) ^ (-s - 1) * ((Int.fract t : ℝ) : ℂ) = F s := by
    rw [F_eq_integral, integral_Ici_eq_integral_Ioi]
  have hs1 : s - 1 ≠ 0 := by
    intro h; have := congrArg Complex.re h; simp at this; linarith
  have hs1' : -s + 1 ≠ 0 := by
    intro h; apply hs1; linear_combination -h
  have h4 : ∫ t in Ioi (1:ℝ), (t : ℂ) ^ (-s) = 1 / (s - 1) := by
    rw [integral_Ioi_cpow_of_lt (by rw [neg_re]; linarith) zero_lt_one, ofReal_one, one_cpow]
    field_simp
    ring
  rw [h1, integral_sub hi1 hi2, h3, h4]
  ring

/-! ## C. Analytic continuation to {Re s > 0} \ {1} -/

lemma zeta_eq {s : ℂ} (hs : 0 < s.re) (hs1 : s ≠ 1) :
    riemannZeta s = s / (s - 1) - s * F s := by
  have hUo : IsOpen ({s : ℂ | 0 < s.re} \ {1}) :=
    (isOpen_lt continuous_const Complex.continuous_re).sdiff isClosed_singleton
  have hd1 : DifferentiableOn ℂ riemannZeta ({s : ℂ | 0 < s.re} \ {1}) :=
    fun z hz => (differentiableAt_riemannZeta hz.2).differentiableWithinAt
  have hd2 : DifferentiableOn ℂ (fun z => z / (z - 1) - z * F z) ({s : ℂ | 0 < s.re} \ {1}) := by
    intro z hz
    have hz1 : z - 1 ≠ 0 := sub_ne_zero.mpr hz.2
    have hF := differentiableAt_F hz.1
    have : DifferentiableAt ℂ (fun z => z / (z - 1) - z * F z) z := by
      fun_prop (disch := exact hz1)
    exact this.differentiableWithinAt
  have hc : (2 : ℂ) ∈ {s : ℂ | 0 < s.re} \ {1} :=
    ⟨show (0:ℝ) < (2:ℂ).re by norm_num, fun h => by norm_num at h⟩
  have hev : riemannZeta =ᶠ[𝓝 (2:ℂ)] (fun z => z / (z - 1) - z * F z) := by
    have hmem : {z : ℂ | 1 < z.re} ∈ 𝓝 (2:ℂ) :=
      (isOpen_lt continuous_const Complex.continuous_re).mem_nhds
        (show (1:ℝ) < (2:ℂ).re by norm_num)
    filter_upwards [hmem] with z hz using zeta_eq_of_one_lt hz
  exact (hd1.analyticOnNhd hUo).eqOn_of_preconnected_of_eventuallyEq (hd2.analyticOnNhd hUo)
    (Fayed.Mertens.isPreconnected_U 0) hc hev ⟨hs, hs1⟩

/-! ## D. M(s) = ∫_0^∞ {1/x} x^{s-1} dx = 1/(1-s) + F(s) on the strip -/

noncomputable def M (s : ℂ) : ℂ := mellin (fun x : ℝ => ((Int.fract (1 / x) : ℝ) : ℂ)) s

lemma lo_eq (s : ℂ) {t : ℝ} (ht : t ∈ Ioo (0:ℝ) 1) :
    (t : ℂ) ^ (-s - 1) * ((Int.fract t : ℝ) : ℂ) = (t : ℂ) ^ (-s) := by
  have ht0 : (t : ℂ) ≠ 0 := by exact_mod_cast ht.1.ne'
  rw [Int.fract_eq_self.mpr ⟨ht.1.le, ht.2⟩, cpow_sub _ _ ht0, cpow_one]
  field_simp

lemma integrableOn_lo {s : ℂ} (hs : s.re < 1) :
    IntegrableOn (fun t : ℝ => (t : ℂ) ^ (-s - 1) * ((Int.fract t : ℝ) : ℂ)) (Ioo 0 1) := by
  have h : IntegrableOn (fun t : ℝ => (t : ℂ) ^ (-s)) (Ioo 0 1) := by
    have := intervalIntegral.intervalIntegrable_cpow' (a := 0) (b := 1) (r := -s)
      (by rw [neg_re]; linarith)
    rw [intervalIntegrable_iff_integrableOn_Ioc_of_le zero_le_one] at this
    exact this.mono_set Ioo_subset_Ioc_self
  exact h.congr_fun (fun t ht => (lo_eq s ht).symm) measurableSet_Ioo

lemma disj : Disjoint (Ioo (0:ℝ) 1) (Ici 1) :=
  Set.disjoint_left.mpr fun _ ht ht' => absurd ht.2 (not_lt.mpr ht')

lemma mellinConvergent_fract {s : ℂ} (hs0 : 0 < s.re) (hs1 : s.re < 1) :
    MellinConvergent (fun t : ℝ => ((Int.fract t : ℝ) : ℂ)) (-s) := by
  unfold MellinConvergent
  rw [← Ioo_union_Ici_eq_Ioi zero_lt_one]
  simp_rw [smul_eq_mul]
  exact (integrableOn_lo hs1).union (integrableOn_hi hs0)

lemma mellin_fract {s : ℂ} (hs0 : 0 < s.re) (hs1 : s.re < 1) :
    mellin (fun t : ℝ => ((Int.fract t : ℝ) : ℂ)) (-s) = 1 / (1 - s) + F s := by
  have hs : -s + 1 ≠ 0 := by
    intro h; have := congrArg Complex.re h; simp at this; linarith
  rw [mellin, ← Ioo_union_Ici_eq_Ioi zero_lt_one]
  simp_rw [smul_eq_mul]
  rw [setIntegral_union disj measurableSet_Ici (integrableOn_lo hs1) (integrableOn_hi hs0),
    F_eq_integral]
  congr 1
  rw [setIntegral_congr_fun measurableSet_Ioo (fun t ht => lo_eq s ht),
    ← integral_Ioc_eq_integral_Ioo, ← intervalIntegral.integral_of_le zero_le_one,
    integral_cpow (Or.inl (by rw [neg_re]; linarith)), ofReal_one, ofReal_zero, one_cpow,
    zero_cpow hs, sub_zero, neg_add_eq_sub]

lemma M_eq {s : ℂ} (hs0 : 0 < s.re) (hs1 : s.re < 1) : M s = 1 / (1 - s) + F s := by
  rw [M, ← mellin_fract hs0 hs1, ← mellin_comp_inv]
  simp_rw [one_div]

lemma mellinConvergent_M {s : ℂ} (hs0 : 0 < s.re) (hs1 : s.re < 1) :
    MellinConvergent (fun x : ℝ => ((Int.fract (1 / x) : ℝ) : ℂ)) s := by
  have h := (MellinConvergent.comp_rpow (f := fun t : ℝ => ((Int.fract t : ℝ) : ℂ)) (s := s)
    (a := -1) (by norm_num)).mpr (by rw [ofReal_neg, ofReal_one, div_neg, div_one]; exact mellinConvergent_fract hs0 hs1)
  have e : (fun x : ℝ => ((Int.fract (1 / x) : ℝ) : ℂ)) =
      fun t => ((Int.fract (t ^ (-1:ℝ)) : ℝ) : ℂ) := by
    ext x; rw [Real.rpow_neg_one, one_div]
  rw [e]; exact h

lemma M_eq_zero {s : ℂ} (hs0 : 0 < s.re) (hs1 : s.re < 1) (hz : riemannZeta s = 0) : M s = 0 := by
  have hs : s ≠ 0 := by rintro rfl; simp at hs0
  have hne : s ≠ 1 := by rintro rfl; simp at hs1
  have hs1' : s - 1 ≠ 0 := sub_ne_zero.mpr hne
  have h := zeta_eq hs0 hne
  rw [hz] at h
  have h2 : s * (F s - 1 / (s - 1)) = 0 := by
    rw [mul_sub, mul_one_div]; linear_combination h
  have hF : F s = 1 / (s - 1) := sub_eq_zero.mp ((mul_eq_zero.mp h2).resolve_left hs)
  rw [M_eq hs0 hs1, hF, show (1:ℂ) - s = -(s - 1) by ring, one_div_neg_eq_neg_one_div]
  ring

/-! ## E/F. ∫_0^1 rho_θ x^{s-1} = θ^s M(s) - θ M(s) -/

lemma fract_div_eq (θ x : ℝ) : Int.fract (θ / x) = Int.fract (1 / (θ⁻¹ * x)) := by
  rw [one_div, mul_inv, inv_inv, div_eq_mul_inv]

lemma mellin_theta {θ : ℝ} (hθ : 0 < θ) (s : ℂ) :
    mellin (fun x : ℝ => ((Int.fract (θ / x) : ℝ) : ℂ)) s = ((θ⁻¹ : ℝ) : ℂ) ^ (-s) • M s := by
  simp_rw [fract_div_eq θ]
  exact mellin_comp_mul_left (fun x : ℝ => ((Int.fract (1 / x) : ℝ) : ℂ)) s (inv_pos.mpr hθ)

lemma mellinConvergent_theta {θ : ℝ} (hθ : 0 < θ) {s : ℂ} (hs0 : 0 < s.re) (hs1 : s.re < 1) :
    MellinConvergent (fun x : ℝ => ((Int.fract (θ / x) : ℝ) : ℂ)) s := by
  simp_rw [fract_div_eq θ]
  exact (MellinConvergent.comp_mul_left (f := fun x : ℝ => ((Int.fract (1 / x) : ℝ) : ℂ))
    (inv_pos.mpr hθ)).mpr (mellinConvergent_M hs0 hs1)

lemma hasMellin_rho {θ : ℝ} (hθ : 0 < θ) {s : ℂ} (hs0 : 0 < s.re) (hs1 : s.re < 1) :
    HasMellin (fun x : ℝ => ((rhoNB θ x : ℝ) : ℂ)) s
      (((θ⁻¹ : ℝ) : ℂ) ^ (-s) • M s - (θ : ℂ) • M s) := by
  have h1 := hasMellin_const_smul (mellinConvergent_M hs0 hs1) (θ : ℂ)
  have h := hasMellin_sub (mellinConvergent_theta hθ hs0 hs1) h1.1
  rw [h1.2, mellin_theta hθ] at h
  have e : (fun x : ℝ => ((rhoNB θ x : ℝ) : ℂ)) = fun x =>
      ((Int.fract (θ / x) : ℝ) : ℂ) - (θ : ℂ) • ((Int.fract (1 / x) : ℝ) : ℂ) := by
    ext x; simp [rhoNB]
  rw [e]; exact h

lemma rhoNB_eq_zero {θ x : ℝ} (hθ : 0 < θ) (hθ1 : θ ≤ 1) (hx : 1 < x) : rhoNB θ x = 0 := by
  have hx0 : 0 < x := zero_lt_one.trans hx
  unfold rhoNB
  rw [Int.fract_eq_self.mpr ⟨(div_pos hθ hx0).le, (div_lt_one hx0).mpr (hθ1.trans_lt hx)⟩,
    Int.fract_eq_self.mpr ⟨(div_pos one_pos hx0).le, (div_lt_one hx0).mpr hx⟩]
  ring

lemma integral_rho_eq_zero {θ : ℝ} (hθ : 0 < θ) (hθ1 : θ ≤ 1) {s : ℂ} (hs0 : 0 < s.re)
    (hs1 : s.re < 1) (hz : riemannZeta s = 0) :
    ∫ x in Ioc (0:ℝ) 1, ((rhoNB θ x : ℝ) : ℂ) * (x : ℂ) ^ (s - 1) = 0 := by
  have h := hasMellin_rho hθ hs0 hs1
  rw [M_eq_zero hs0 hs1 hz, smul_zero, smul_zero, sub_zero] at h
  rw [← h.2, mellin, setIntegral_eq_of_subset_of_ae_sdiff_eq_zero
    measurableSet_Ioi.nullMeasurableSet Ioc_subset_Ioi_self
    (Eventually.of_forall fun x hx => ?_)]
  · refine setIntegral_congr_fun measurableSet_Ioc fun x _ => ?_
    rw [smul_eq_mul, mul_comm]
  · rw [rhoNB_eq_zero hθ hθ1 (not_le.mp fun h => hx.2 ⟨hx.1, h⟩), ofReal_zero, smul_zero]

lemma integrableOn_rho {θ : ℝ} (hθ : 0 < θ) {s : ℂ} (hs0 : 0 < s.re) (hs1 : s.re < 1) :
    IntegrableOn (fun x : ℝ => ((rhoNB θ x : ℝ) : ℂ) * (x : ℂ) ^ (s - 1)) (Ioc 0 1) := by
  refine (IntegrableOn.mono_set (hasMellin_rho hθ hs0 hs1).1 Ioc_subset_Ioi_self).congr_fun
    (fun x _ => ?_) measurableSet_Ioc
  dsimp only
  rw [smul_eq_mul, mul_comm]

/-! ## G. The L^2 contradiction -/

lemma measurable_rhoNB (θ : ℝ) : Measurable (rhoNB θ) :=
  (measurable_fract.comp (measurable_const.div measurable_id)).sub
    (measurable_const.mul (measurable_fract.comp (measurable_const.div measurable_id)))

lemma abs_rhoNB_le {θ : ℝ} (hθ : 0 < θ) (x : ℝ) : |rhoNB θ x| ≤ 1 + θ := by
  unfold rhoNB
  have h1 := Int.fract_nonneg (θ / x)
  have h2 := Int.fract_lt_one (θ / x)
  have h3 := Int.fract_nonneg (1 / x)
  have h4 := Int.fract_lt_one (1 / x)
  rw [abs_le]; constructor <;> nlinarith

lemma amgm (lam a b : ℝ) (h : 0 < lam) : |a| * b ≤ (lam * a ^ 2 + b ^ 2 / lam) / 2 := by
  obtain ⟨t, rfl⟩ : ∃ t, b = lam * t := ⟨b / lam, by field_simp⟩
  have e : (lam * t) ^ 2 / lam = lam * t ^ 2 := by field_simp
  rw [e]
  nlinarith [mul_nonneg h.le (sq_nonneg (|a| - t)), sq_abs a]

theorem no_zero_of_NBApprox (h : NBApprox) {s : ℂ} (hs0 : 1 / 2 < s.re) (hs1 : s.re < 1) :
    riemannZeta s ≠ 0 := by
  intro hz
  have hσ0 : 0 < s.re := by linarith
  have hs : s ≠ 0 := by rintro rfl; simp at hσ0
  have hwint : IntegrableOn (fun x : ℝ => (x : ℂ) ^ (s - 1)) (Ioc 0 1) := by
    have := intervalIntegral.intervalIntegrable_cpow' (a := 0) (b := 1) (r := s - 1)
      (by simp; linarith)
    rwa [intervalIntegrable_iff_integrableOn_Ioc_of_le zero_le_one] at this
  have hwI : ∫ x in Ioc (0:ℝ) 1, (x : ℂ) ^ (s - 1) = 1 / s := by
    rw [← intervalIntegral.integral_of_le zero_le_one,
      integral_cpow (Or.inl (by simp; linarith)), sub_add_cancel, ofReal_one, ofReal_zero,
      one_cpow, zero_cpow hs, sub_zero]
  have hBint : IntegrableOn (fun x : ℝ => x ^ (2 * s.re - 2)) (Ioc 0 1) := by
    have := intervalIntegral.intervalIntegrable_rpow' (a := 0) (b := 1) (r := 2 * s.re - 2)
      (by linarith)
    rwa [intervalIntegrable_iff_integrableOn_Ioc_of_le zero_le_one] at this
  set B := ∫ x in Ioc (0:ℝ) 1, x ^ (2 * s.re - 2) with hB
  have hB0 : 0 ≤ B := setIntegral_nonneg measurableSet_Ioc fun x hx => Real.rpow_nonneg hx.1.le _
  set c := ‖1 / s‖ with hc
  have hc0 : 0 < c := norm_pos_iff.mpr (one_div_ne_zero hs)
  obtain ⟨n, a, θ, hθ, hA⟩ := h (c ^ 2 / (B + 1)) (by positivity)
  set f : ℝ → ℝ := fun x => 1 - ∑ k, a k * rhoNB (θ k) x with hf
  set A := ∫ x in Ioc (0:ℝ) 1, f x ^ 2 with hAdef
  have hfm : Measurable f := measurable_const.sub
    (Finset.measurable_sum _ fun k _ => measurable_const.mul (measurable_rhoNB _))
  set K := 1 + ∑ k, |a k| * (1 + θ k) with hK
  have hfb : ∀ x, |f x| ≤ K := by
    intro x
    have h1 : |f x| ≤ 1 + |∑ k, a k * rhoNB (θ k) x| := by
      have := abs_add_le (1 : ℝ) (-(∑ k, a k * rhoNB (θ k) x))
      simpa [hf, sub_eq_add_neg] using this
    have h2 : |∑ k, a k * rhoNB (θ k) x| ≤ ∑ k, |a k| * (1 + θ k) := by
      refine (Finset.abs_sum_le_sum_abs _ _).trans (Finset.sum_le_sum fun k _ => ?_)
      rw [abs_mul]
      exact mul_le_mul_of_nonneg_left (abs_rhoNB_le (hθ k).1 x) (abs_nonneg _)
    linarith
  have hf2int : Integrable (fun x => f x ^ 2) (volume.restrict (Ioc (0:ℝ) 1)) := by
    refine IntegrableOn.of_bound measure_Ioc_lt_top (hfm.pow_const 2).aestronglyMeasurable (K ^ 2)
      (Eventually.of_forall fun x => ?_)
    rw [Real.norm_eq_abs, abs_pow]
    exact pow_le_pow_left₀ (abs_nonneg _) (hfb x) 2
  have hrho : ∀ k, Integrable (fun x : ℝ => ((rhoNB (θ k) x : ℝ) : ℂ) * (x : ℂ) ^ (s - 1))
      (volume.restrict (Ioc (0:ℝ) 1)) := fun k => integrableOn_rho (hθ k).1 hσ0 hs1
  have hfI : ∫ x in Ioc (0:ℝ) 1, ((f x : ℝ) : ℂ) * (x : ℂ) ^ (s - 1) = 1 / s := by
    have e : ∀ x : ℝ, ((f x : ℝ) : ℂ) * (x : ℂ) ^ (s - 1) = (x : ℂ) ^ (s - 1) -
        ∑ k, (a k : ℂ) * (((rhoNB (θ k) x : ℝ) : ℂ) * (x : ℂ) ^ (s - 1)) := by
      intro x
      simp only [hf]
      push_cast
      rw [sub_mul, one_mul, Finset.sum_mul]
      congr 1
      exact Finset.sum_congr rfl fun k _ => by ring
    simp_rw [e]
    rw [integral_sub hwint (integrable_finsetSum _ fun k _ => (hrho k).const_mul _),
      integral_finsetSum _ fun k _ => (hrho k).const_mul _,
      Finset.sum_eq_zero fun k _ => by
        rw [integral_const_mul, integral_rho_eq_zero (hθ k).1 (hθ k).2 hσ0 hs1 hz, mul_zero],
      sub_zero, hwI]
  set lam := (B + 1) / c with hlam
  have hlam0 : 0 < lam := by positivity
  have hgint : Integrable (fun x : ℝ => (lam * f x ^ 2 + x ^ (2 * s.re - 2) / lam) / 2)
      (volume.restrict (Ioc (0:ℝ) 1)) :=
    ((hf2int.const_mul lam).add (Integrable.div_const hBint lam)).div_const 2
  have hbound : ‖∫ x in Ioc (0:ℝ) 1, ((f x : ℝ) : ℂ) * (x : ℂ) ^ (s - 1)‖ ≤
      ∫ x in Ioc (0:ℝ) 1, (lam * f x ^ 2 + x ^ (2 * s.re - 2) / lam) / 2 := by
    refine norm_integral_le_of_norm_le hgint ?_
    refine (ae_restrict_iff' measurableSet_Ioc).mpr (Eventually.of_forall fun x hx => ?_)
    rw [norm_mul, Complex.norm_real, Real.norm_eq_abs, norm_cpow_eq_rpow_re_of_pos hx.1]
    have hre : (s - 1).re = s.re - 1 := by simp
    have hsq : x ^ (2 * s.re - 2) = (x ^ (s.re - 1)) ^ 2 := by
      rw [← Real.rpow_natCast, ← Real.rpow_mul hx.1.le]
      congr 1; push_cast; ring
    rw [hre, hsq]
    exact amgm lam (f x) _ hlam0
  have hint_eq : ∫ x in Ioc (0:ℝ) 1, (lam * f x ^ 2 + x ^ (2 * s.re - 2) / lam) / 2 =
      (lam * A + B / lam) / 2 := by
    rw [integral_div, integral_add (hf2int.const_mul lam) (Integrable.div_const hBint lam),
      integral_const_mul, integral_div]
  rw [hfI, hint_eq, ← hc] at hbound
  have hA' : lam * A < c := by
    calc lam * A < lam * (c ^ 2 / (B + 1)) := mul_lt_mul_of_pos_left hA hlam0
      _ = c := by rw [hlam]; field_simp
  have hB' : B / lam ≤ c := by
    rw [hlam, div_div_eq_mul_div, div_le_iff₀ (by positivity)]
    nlinarith
  linarith

theorem RH_of_NBApprox (h : NBApprox) : RiemannHypothesis := by
  rw [Fayed.R7.RH_iff_QuasiRH_half]
  intro s hs _
  rcases le_or_gt 1 s.re with h1 | h1
  · exact riemannZeta_ne_zero_of_one_le_re h1
  · exact no_zero_of_NBApprox h (by linarith) h1

end Fayed.NB

#print axioms Fayed.NB.RH_of_NBApprox
