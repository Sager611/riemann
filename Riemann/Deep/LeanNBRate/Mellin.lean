import Riemann.Deep.LeanNBRate.LitInputsU1
import Mathlib.Analysis.MellinTransform
import Mathlib.Analysis.SpecialFunctions.ImproperIntegrals
import Mathlib.MeasureTheory.Integral.IntegralEqImproper
import Mathlib.MeasureTheory.Function.Floor

/-!
# Mellin transforms of the Nyman–Beurling functions and the Mellin–Plancherel identity
-/

open MeasureTheory Filter Topology Set Complex

noncomputable section

namespace Riemann.Deep.LeanNBRate

/-- The fractional part is Mellin-convergent at -s for 0 < Re s < 1. -/
theorem mellinConvergent_fract {s : ℂ} (h0 : 0 < s.re) (h1 : s.re < 1) :
    MellinConvergent (fun x : ℝ => ((Int.fract x : ℝ) : ℂ)) (-s) := by
  unfold MellinConvergent
  have hsplit : Ioi (0 : ℝ) = Ioc 0 1 ∪ Ioi 1 := (Ioc_union_Ioi_eq_Ioi zero_le_one).symm
  rw [hsplit]
  refine IntegrableOn.union ?_ ?_
  · -- on (0, 1]: fract x = x except at 1
    have hc : IntegrableOn (fun x : ℝ => (x : ℂ) ^ (-s)) (Ioc 0 1) := by
      have := intervalIntegral.intervalIntegrable_cpow' (a := 0) (b := 1) (r := -s)
        (by simp; linarith)
      exact (intervalIntegrable_iff_integrableOn_Ioc_of_le zero_le_one).mp this
    refine hc.congr_fun_ae ?_
    have hae : ∀ᵐ x ∂(volume.restrict (Ioc (0 : ℝ) 1)), x ∈ Ioo (0 : ℝ) 1 := by
      rw [ae_restrict_iff' measurableSet_Ioc]
      have : ∀ᵐ x ∂(volume : Measure ℝ), x ≠ 1 := by
        rw [ae_iff]; simp [measure_singleton]
      filter_upwards [this] with x hx hxI
      exact ⟨hxI.1, lt_of_le_of_ne hxI.2 hx⟩
    filter_upwards [hae] with x hx
    have hfr : Int.fract x = x := Int.fract_eq_self.mpr ⟨hx.1.le, hx.2⟩
    simp only [hfr, smul_eq_mul]
    have hx0 : (x : ℂ) ≠ 0 := by exact_mod_cast hx.1.ne'
    rw [show -s - 1 = -s + (-1) by ring, cpow_add _ _ hx0, cpow_neg_one]
    field_simp
  · -- on (1, ∞): |fract| ≤ 1
    have hc : IntegrableOn (fun x : ℝ => (x : ℂ) ^ (-s - 1)) (Ioi 1) :=
      integrableOn_Ioi_cpow_of_lt (by simp; linarith) one_pos
    refine hc.norm.mono' ?_ ?_
    · have hcont : ContinuousOn (fun x : ℝ => (x : ℂ) ^ (-s - 1)) (Ioi 1) := by
        intro x hx
        have hx1 : (1 : ℝ) < x := hx
        exact (continuousAt_ofReal_cpow_const _ _ (Or.inr (ne_of_gt (by linarith)))).continuousWithinAt
      have h1 : AEStronglyMeasurable (fun x : ℝ => (x : ℂ) ^ (-s - 1)) (volume.restrict (Ioi (1 : ℝ))) :=
        hcont.aestronglyMeasurable measurableSet_Ioi
      have h2 : AEStronglyMeasurable (fun x : ℝ => ((Int.fract x : ℝ) : ℂ)) (volume.restrict (Ioi (1 : ℝ))) :=
        (Complex.measurable_ofReal.comp measurable_fract).aestronglyMeasurable
      exact h1.smul h2
    · rw [ae_restrict_iff' measurableSet_Ioi]
      refine Filter.Eventually.of_forall fun x hx => ?_
      rw [norm_smul]
      have h1 : ‖((Int.fract x : ℝ) : ℂ)‖ ≤ 1 := by
        rw [Complex.norm_real, Real.norm_eq_abs, abs_of_nonneg (Int.fract_nonneg x)]
        exact (Int.fract_lt_one x).le
      calc ‖(x : ℂ) ^ (-s - 1)‖ * ‖((Int.fract x : ℝ) : ℂ)‖ ≤ ‖(x : ℂ) ^ (-s - 1)‖ * 1 := by gcongr
        _ = ‖(x : ℂ) ^ (-s - 1)‖ := mul_one _

/-- The Mellin transform of the fractional part at -s (Titchmarsh (2.1.5)). -/
theorem mellin_fract (h215 : Titchmarsh_2_1_5) {s : ℂ} (h0 : 0 < s.re) (h1 : s.re < 1) :
    mellin (fun x : ℝ => ((Int.fract x : ℝ) : ℂ)) (-s) = -riemannZeta s / s := by
  have hs0 : s ≠ 0 := fun h => by simp [h] at h0
  have hT := h215 s h0 h1
  unfold mellin
  have e : ∀ x : ℝ, ((⌊x⌋ : ℂ) - x) * (x : ℂ) ^ (-s - 1) = -((x : ℂ) ^ (-s - 1) • ((Int.fract x : ℝ) : ℂ)) := by
    intro x
    rw [smul_eq_mul, Int.fract]
    push_cast
    ring
  simp_rw [e, integral_neg] at hT
  rw [eq_div_iff hs0]
  linear_combination hT

/-- **Mellin transform of x ↦ {θ/x}** (θ > 0, 0 < Re s < 1): -θ^s ζ(s)/s. -/
theorem hasMellin_fract_div (h215 : Titchmarsh_2_1_5) {s : ℂ} (h0 : 0 < s.re) (h1 : s.re < 1)
    {θ : ℝ} (hθ : 0 < θ) :
    HasMellin (fun x : ℝ => ((Int.fract (θ / x) : ℝ) : ℂ)) s (-(θ : ℂ) ^ s * riemannZeta s / s) := by
  set f₁ : ℝ → ℂ := fun t => ((Int.fract (t⁻¹) : ℝ) : ℂ)
  have hf : (fun x : ℝ => ((Int.fract (θ / x) : ℝ) : ℂ)) = fun x => f₁ (x * θ⁻¹) := by
    funext x; simp only [f₁, mul_inv, inv_inv, div_eq_mul_inv, mul_comm]
  have hconv₁ : MellinConvergent f₁ s := by
    have := (MellinConvergent.comp_rpow (f := fun x : ℝ => ((Int.fract x : ℝ) : ℂ)) (s := s)
      (a := -1) (by norm_num)).mpr (by
        push_cast; rw [div_neg, div_one]; exact mellinConvergent_fract h0 h1)
    simpa [f₁, Real.rpow_neg_one] using this
  have hm₁ : mellin f₁ s = -riemannZeta s / s := by
    exact (mellin_comp_inv (fun x : ℝ => ((Int.fract x : ℝ) : ℂ)) s).trans (mellin_fract h215 h0 h1)
  refine ⟨?_, ?_⟩
  · rw [hf]
    have := (MellinConvergent.comp_mul_left (f := f₁) (s := s) (a := θ⁻¹) (inv_pos.mpr hθ)).mpr hconv₁
    simpa [mul_comm] using this
  · rw [hf, mellin_comp_mul_right _ _ (inv_pos.mpr hθ), hm₁]
    have hθc : ((θ⁻¹ : ℝ) : ℂ) ^ (-s) = (θ : ℂ) ^ s := by
      rw [ofReal_inv, inv_cpow _ _ (by rw [Complex.arg_ofReal_of_nonneg hθ.le]; exact Real.pi_ne_zero.symm),
        cpow_neg, inv_inv]
    rw [hθc, smul_eq_mul]
    ring


theorem integrable_exp_neg_abs_half : Integrable (fun v : ℝ => Real.exp (-|v| / 2)) := by
  have hI : IntegrableOn (fun v : ℝ => Real.exp (-|v| / 2)) (Ioi 0) := by
    refine (exp_neg_integrableOn_Ioi 0 (by norm_num : (0 : ℝ) < 1 / 2)).congr_fun (fun v hv => ?_)
      measurableSet_Ioi
    have hv' : (0 : ℝ) < v := hv
    simp only [abs_of_pos hv']
    ring_nf
  have int_Iic : IntegrableOn (fun v : ℝ => Real.exp (-|v| / 2)) (Iic 0) := by
    rw [← Measure.map_neg_eq_self (volume : Measure ℝ)]
    let m : MeasurableEmbedding fun x : ℝ => -x := (Homeomorph.neg ℝ).measurableEmbedding
    rw [m.integrableOn_map_iff]
    simp_rw [Function.comp_def, abs_neg, neg_preimage, neg_Iic, neg_zero]
    exact Iff.mpr integrableOn_Ici_iff_integrableOn_Ioi hI
  have := int_Iic.union hI
  rwa [Iic_union_Ioi, integrableOn_univ] at this

/-- Substitution x = e^{-v} for integrals over (0, ∞). -/
theorem integral_Ioi_eq_exp_neg {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E] (g : ℝ → E) :
    ∫ x in Ioi (0 : ℝ), g x = ∫ v : ℝ, Real.exp (-v) • g (Real.exp (-v)) := by
  rw [← integral_comp_exp g]
  exact (integral_neg_eq_self (fun x : ℝ => Real.exp x • g (Real.exp x)) volume).symm

theorem ofReal_exp_cpow (a : ℝ) (z : ℂ) : ((Real.exp a : ℝ) : ℂ) ^ z = Complex.exp (a * z) := by
  rw [Complex.ofReal_exp, Complex.cpow_def_of_ne_zero (Complex.exp_ne_zero _),
    Complex.log_exp (by simp [Real.pi_pos]) (by simp [Real.pi_pos.le])]

/-- **Mellin–Plancherel** on the critical line (from Rudin 9.13). -/
theorem mellin_plancherel (hP : Rudin1987_Thm9_13) (F : ℝ → ℂ)
    (hi : Integrable (fun v : ℝ => (Real.exp (-v / 2) : ℂ) * F (Real.exp (-v))))
    (h2 : MemLp (fun v : ℝ => (Real.exp (-v / 2) : ℂ) * F (Real.exp (-v))) 2) :
    ∫ t : ℝ, ‖mellin F (1 / 2 + t * I)‖ ^ 2 = 2 * Real.pi * ∫ x in Ioi (0 : ℝ), ‖F x‖ ^ 2 := by
  set h := fun v : ℝ => (Real.exp (-v / 2) : ℂ) * F (Real.exp (-v))
  have e1 : ∀ t : ℝ, mellin F (1 / 2 + t * I) = ∫ v : ℝ, h v * Complex.exp (-(v * t * I)) := by
    intro t
    unfold mellin
    rw [integral_Ioi_eq_exp_neg]
    congr 1; ext v
    simp only [h, Complex.real_smul, smul_eq_mul, ofReal_exp_cpow]
    rw [Complex.ofReal_exp, Complex.ofReal_exp]
    rw [show (Complex.exp ((-v : ℝ) : ℂ)) * (Complex.exp (((-v : ℝ) : ℂ) * (1 / 2 + t * I - 1)) *
        F (Real.exp (-v))) = (Complex.exp ((-v : ℝ) : ℂ) * Complex.exp (((-v : ℝ) : ℂ) * (1 / 2 + t * I - 1))) *
        F (Real.exp (-v)) by ring]
    rw [← Complex.exp_add]
    rw [show Complex.exp (((-v / 2 : ℝ) : ℂ)) * F (Real.exp (-v)) * Complex.exp (-(v * t * I)) =
        (Complex.exp (((-v / 2 : ℝ) : ℂ)) * Complex.exp (-(v * t * I))) * F (Real.exp (-v)) by ring]
    rw [← Complex.exp_add]
    congr 2
    push_cast
    ring
  have e2 : ∫ x in Ioi (0 : ℝ), ‖F x‖ ^ 2 = ∫ v : ℝ, ‖h v‖ ^ 2 := by
    rw [integral_Ioi_eq_exp_neg]
    congr 1; ext v
    simp only [h, smul_eq_mul, norm_mul, Complex.norm_real, Real.norm_eq_abs, mul_pow, sq_abs]
    have hx : Real.exp (-v / 2) ^ 2 = Real.exp (-v) := by rw [sq, ← Real.exp_add]; ring_nf
    rw [hx]
  simp_rw [e1, e2]
  exact hP h hi h2

end Riemann.Deep.LeanNBRate

