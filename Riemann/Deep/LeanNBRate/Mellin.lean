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
    MemLp (fun t : ℝ => mellin F (1 / 2 + t * I)) 2 ∧
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
  have he1 : (fun t : ℝ => mellin F (1 / 2 + t * I)) =
      fun t : ℝ => ∫ v : ℝ, h v * Complex.exp (-(v * t * I)) := funext e1
  refine ⟨?_, ?_⟩
  · rw [he1]; exact (hP h hi h2).1
  · simp_rw [e1]; rw [e2]; exact (hP h hi h2).2


/-! ## The Nyman–Beurling functions -/

/-- χ_{(0,1]} − Σ c_k {θ_k/x} (complex coefficients). -/
def nbF {n : ℕ} (c : Fin n → ℂ) (θ : Fin n → ℝ) (x : ℝ) : ℂ :=
  (Ioc (0 : ℝ) 1).indicator (fun _ => (1 : ℂ)) x - ∑ k, c k * ((Int.fract (θ k / x) : ℝ) : ℂ)

theorem fract_div_eq {θ x : ℝ} (hθ : 0 < θ) (hθx : θ < x) : Int.fract (θ / x) = θ / x := by
  have hx : 0 < x := by linarith
  exact Int.fract_eq_self.mpr ⟨by positivity, (div_lt_one hx).mpr hθx⟩

theorem measurable_nbF {n : ℕ} (c : Fin n → ℂ) (θ : Fin n → ℝ) : Measurable (nbF c θ) := by
  unfold nbF
  refine (measurable_const.indicator measurableSet_Ioc).sub ?_
  refine Finset.measurable_sum _ fun k _ => measurable_const.mul ?_
  exact Complex.measurable_ofReal.comp (measurable_fract.comp (measurable_const.div measurable_id))

theorem norm_nbF_exp_le {n : ℕ} (c : Fin n → ℂ) (θ : Fin n → ℝ) (hθ : ∀ k, θ k ∈ Ioc (0 : ℝ) 1)
    (v : ℝ) : ‖(Real.exp (-v / 2) : ℂ) * nbF c θ (Real.exp (-v))‖ ≤
      (1 + ∑ k, ‖c k‖) * Real.exp (-|v| / 2) := by
  have hS : 0 ≤ ∑ k, ‖c k‖ := Finset.sum_nonneg fun k _ => norm_nonneg _
  rw [norm_mul, Complex.norm_real, Real.norm_eq_abs, abs_of_pos (Real.exp_pos _)]
  rcases le_or_gt 0 v with hv | hv
  · rw [abs_of_nonneg hv]
    have hb : ‖nbF c θ (Real.exp (-v))‖ ≤ 1 + ∑ k, ‖c k‖ := by
      unfold nbF
      refine (norm_sub_le _ _).trans (add_le_add ?_ ?_)
      · by_cases hx : Real.exp (-v) ∈ Ioc (0 : ℝ) 1
        · rw [indicator_of_mem hx]; simp
        · rw [indicator_of_notMem hx]; simp
      · refine (norm_sum_le _ _).trans (Finset.sum_le_sum fun k _ => ?_)
        rw [norm_mul, Complex.norm_real, Real.norm_eq_abs, abs_of_nonneg (Int.fract_nonneg _)]
        calc ‖c k‖ * Int.fract (θ k / Real.exp (-v)) ≤ ‖c k‖ * 1 :=
              mul_le_mul_of_nonneg_left (Int.fract_lt_one _).le (norm_nonneg _)
          _ = ‖c k‖ := mul_one _
    calc Real.exp (-v / 2) * ‖nbF c θ (Real.exp (-v))‖ ≤ Real.exp (-v / 2) * (1 + ∑ k, ‖c k‖) :=
          mul_le_mul_of_nonneg_left hb (Real.exp_pos _).le
      _ = _ := by ring
  · rw [abs_of_neg hv]
    have hx1 : 1 < Real.exp (-v) := by rw [← Real.exp_zero]; exact Real.exp_lt_exp.mpr (by linarith)
    have hval : nbF c θ (Real.exp (-v)) = -((∑ k, c k * (θ k : ℂ)) * (Real.exp v : ℂ)) := by
      unfold nbF
      rw [indicator_of_notMem (fun h => by have := h.2; linarith)]
      have e : ∀ k, ((Int.fract (θ k / Real.exp (-v)) : ℝ) : ℂ) = (θ k : ℂ) * (Real.exp v : ℂ) := by
        intro k
        rw [fract_div_eq (hθ k).1 (by linarith [(hθ k).2]), div_eq_mul_inv, ← Real.exp_neg, neg_neg]
        push_cast; ring
      simp_rw [e, Finset.sum_mul]
      rw [zero_sub, neg_inj]
      refine Finset.sum_congr rfl fun k _ => by ring
    rw [hval, norm_neg, norm_mul, Complex.norm_real, Real.norm_eq_abs, abs_of_pos (Real.exp_pos _)]
    have hsum : ‖∑ k, c k * (θ k : ℂ)‖ ≤ ∑ k, ‖c k‖ := by
      refine (norm_sum_le _ _).trans (Finset.sum_le_sum fun k _ => ?_)
      rw [norm_mul, Complex.norm_real, Real.norm_eq_abs, abs_of_pos (hθ k).1]
      calc ‖c k‖ * θ k ≤ ‖c k‖ * 1 := mul_le_mul_of_nonneg_left (hθ k).2 (norm_nonneg _)
        _ = ‖c k‖ := mul_one _
    have e2 : Real.exp (-v / 2) * Real.exp v = Real.exp (- -v / 2) := by
      rw [← Real.exp_add]; ring_nf
    calc Real.exp (-v / 2) * (‖∑ k, c k * (θ k : ℂ)‖ * Real.exp v)
        = ‖∑ k, c k * (θ k : ℂ)‖ * (Real.exp (-v / 2) * Real.exp v) := by ring
      _ = ‖∑ k, c k * (θ k : ℂ)‖ * Real.exp (- -v / 2) := by rw [e2]
      _ ≤ (1 + ∑ k, ‖c k‖) * Real.exp (- -v / 2) := by
          apply mul_le_mul_of_nonneg_right _ (Real.exp_pos _).le; linarith

theorem integrable_nbF_exp {n : ℕ} (c : Fin n → ℂ) (θ : Fin n → ℝ) (hθ : ∀ k, θ k ∈ Ioc (0 : ℝ) 1) :
    Integrable (fun v : ℝ => (Real.exp (-v / 2) : ℂ) * nbF c θ (Real.exp (-v))) := by
  refine (integrable_exp_neg_abs_half.const_mul (1 + ∑ k, ‖c k‖)).mono' ?_ ?_
  · exact ((Complex.measurable_ofReal.comp (Real.measurable_exp.comp
      ((measurable_id.neg).div_const 2))).mul ((measurable_nbF c θ).comp
      (Real.measurable_exp.comp measurable_id.neg))).aestronglyMeasurable
  · exact Filter.Eventually.of_forall fun v => norm_nbF_exp_le c θ hθ v

theorem memLp_nbF_exp {n : ℕ} (c : Fin n → ℂ) (θ : Fin n → ℝ) (hθ : ∀ k, θ k ∈ Ioc (0 : ℝ) 1) :
    MemLp (fun v : ℝ => (Real.exp (-v / 2) : ℂ) * nbF c θ (Real.exp (-v))) 2 := by
  have hm := (integrable_nbF_exp c θ hθ).aestronglyMeasurable
  rw [memLp_two_iff_integrable_sq_norm hm]
  set C := 1 + ∑ k, ‖c k‖
  have hC : 0 ≤ C := by have := Finset.sum_nonneg fun k (_ : k ∈ Finset.univ) => norm_nonneg (c k); linarith
  refine (integrable_exp_neg_abs_half.const_mul (C ^ 2)).mono' (hm.norm.pow 2) ?_
  refine Filter.Eventually.of_forall fun v => ?_
  have h1 := norm_nbF_exp_le c θ hθ v
  have h0 : 0 ≤ ‖(Real.exp (-v / 2) : ℂ) * nbF c θ (Real.exp (-v))‖ := norm_nonneg _
  have h2 : Real.exp (-|v| / 2) ^ 2 ≤ Real.exp (-|v| / 2) := by
    have : Real.exp (-|v| / 2) ≤ 1 := Real.exp_le_one_iff.mpr (by have := abs_nonneg v; linarith)
    have : 0 ≤ Real.exp (-|v| / 2) := (Real.exp_pos _).le
    nlinarith
  rw [Real.norm_eq_abs, abs_of_nonneg (by positivity)]
  calc ‖(Real.exp (-v / 2) : ℂ) * nbF c θ (Real.exp (-v))‖ ^ 2 ≤ (C * Real.exp (-|v| / 2)) ^ 2 :=
        pow_le_pow_left₀ h0 h1 2
    _ = C ^ 2 * Real.exp (-|v| / 2) ^ 2 := by ring
    _ ≤ C ^ 2 * Real.exp (-|v| / 2) := mul_le_mul_of_nonneg_left h2 (sq_nonneg _)

theorem hasMellin_sum {n : ℕ} {s : ℂ} (c : Fin n → ℂ) (f : Fin n → ℝ → ℂ) (m : Fin n → ℂ)
    (hf : ∀ k, HasMellin (f k) s (m k)) :
    HasMellin (fun x => ∑ k, c k * f k x) s (∑ k, c k * m k) := by
  classical
  have key : ∀ S : Finset (Fin n), HasMellin (fun x => ∑ k ∈ S, c k * f k x) s (∑ k ∈ S, c k * m k) := by
    intro S
    induction S using Finset.induction_on with
    | empty =>
      refine ⟨?_, ?_⟩
      · simp [MellinConvergent]
      · simp [mellin]
    | insert a S ha ih =>
      simp only [Finset.sum_insert ha]
      have h1 := hasMellin_const_smul (hf a).1 (c a)
      have h2 := hasMellin_add h1.1 ih.1
      simp only [smul_eq_mul] at h1 h2
      refine ⟨h2.1, ?_⟩
      rw [h2.2, h1.2, (hf a).2, ih.2]
  simpa using key Finset.univ

/-- The Mellin transform of nbF on the critical strip. -/
theorem hasMellin_nbF (h215 : Titchmarsh_2_1_5) {n : ℕ} (c : Fin n → ℂ) (θ : Fin n → ℝ)
    (hθ : ∀ k, θ k ∈ Ioc (0 : ℝ) 1) {s : ℂ} (h0 : 0 < s.re) (h1 : s.re < 1) :
    HasMellin (nbF c θ) s ((1 + riemannZeta s * ∑ k, c k * (θ k : ℂ) ^ s) / s) := by
  have hs0 : s ≠ 0 := fun h => by simp [h] at h0
  have hI := hasMellin_one_Ioc h0
  have hS := hasMellin_sum (s := s) c (fun k x => ((Int.fract (θ k / x) : ℝ) : ℂ))
    (fun k => -((θ k : ℂ) ^ s) * riemannZeta s / s) (fun k => hasMellin_fract_div h215 h0 h1 (hθ k).1)
  have h := hasMellin_sub hI.1 hS.1
  refine ⟨h.1, ?_⟩
  unfold nbF
  rw [h.2, hI.2, hS.2, eq_div_iff hs0, sub_mul, Finset.sum_mul, Finset.mul_sum, one_div,
    inv_mul_cancel₀ hs0, sub_eq_add_neg, ← Finset.sum_neg_distrib]
  congr 1
  refine Finset.sum_congr rfl fun k _ => ?_
  field_simp

/-- **The Mellin–Plancherel identity for Nyman–Beurling combinations.** -/
theorem integral_nbF_sq (hP : Rudin1987_Thm9_13) (h215 : Titchmarsh_2_1_5) {n : ℕ} (c : Fin n → ℂ)
    (θ : Fin n → ℝ) (hθ : ∀ k, θ k ∈ Ioc (0 : ℝ) 1) :
    ∫ x in Ioi (0 : ℝ), ‖nbF c θ x‖ ^ 2 = (2 * Real.pi)⁻¹ * ∫ t : ℝ,
      ‖(1 + riemannZeta (1 / 2 + t * I) * ∑ k, c k * (θ k : ℂ) ^ (1 / 2 + t * I)) / (1 / 2 + t * I)‖ ^ 2 := by
  have h := (mellin_plancherel hP (nbF c θ) (integrable_nbF_exp c θ hθ) (memLp_nbF_exp c θ hθ)).2
  have e : ∀ t : ℝ, mellin (nbF c θ) (1 / 2 + t * I) =
      (1 + riemannZeta (1 / 2 + t * I) * ∑ k, c k * (θ k : ℂ) ^ (1 / 2 + t * I)) / (1 / 2 + t * I) := by
    intro t
    exact (hasMellin_nbF h215 c θ hθ (by simp) (by simp; norm_num)).2
  simp_rw [e] at h
  rw [h]
  field_simp

theorem memLp_line_nbF (hP : Rudin1987_Thm9_13) (h215 : Titchmarsh_2_1_5) {n : ℕ} (c : Fin n → ℂ)
    (θ : Fin n → ℝ) (hθ : ∀ k, θ k ∈ Ioc (0 : ℝ) 1) :
    MemLp (fun t : ℝ =>
      (1 + riemannZeta (1 / 2 + t * I) * ∑ k, c k * (θ k : ℂ) ^ (1 / 2 + t * I)) / (1 / 2 + t * I)) 2 := by
  have h := (mellin_plancherel hP (nbF c θ) (integrable_nbF_exp c θ hθ) (memLp_nbF_exp c θ hθ)).1
  have e : (fun t : ℝ => mellin (nbF c θ) (1 / 2 + t * I)) = fun t : ℝ =>
      (1 + riemannZeta (1 / 2 + t * I) * ∑ k, c k * (θ k : ℂ) ^ (1 / 2 + t * I)) / (1 / 2 + t * I) := by
    funext t
    exact (hasMellin_nbF h215 c θ hθ (by simp) (by simp; norm_num)).2
  rwa [e] at h

end Riemann.Deep.LeanNBRate

