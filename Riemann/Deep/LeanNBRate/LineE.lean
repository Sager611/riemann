import Riemann.Deep.LeanNBRate.Dual5
import Riemann.Deep.LeanNBRate.LitInputsU1
import Mathlib.Analysis.SpecialFunctions.ImproperIntegrals
import Mathlib.Analysis.InnerProductSpace.Calculus

/-!
# E_low on the critical line: the integral of |E_low(s)/s|² over the line is 2π N(p) (Plancherel)

E_low(1/2+it)/(1/2+it) is the Fourier transform of h = conj(p/2 - p') on (0, ℓ], continued by
conj(p ℓ) e^{-(x-ℓ)/2} on (ℓ, ∞); the L² norm of h is N(p).
-/

open MeasureTheory Filter Topology Set Complex ComplexConjugate
open scoped Real InnerProductSpace

noncomputable section

namespace Riemann.Deep.LeanNBRate

/-- The body q = conj(p/2 - p'). -/
def qE (p p' : ℝ → ℂ) (x : ℝ) : ℂ := conj ((1 / 2 : ℂ) * p x - p' x)

/-- The tail conj(p ℓ) e^{-(x-ℓ)/2}. -/
def tE (ℓ : ℝ) (p : ℝ → ℂ) (x : ℝ) : ℂ := conj (p ℓ) * cexp (-((x : ℂ) - ℓ) / 2)

/-- The function whose Fourier transform is E_low(s)/s on the line. -/
def hE (ℓ : ℝ) (p p' : ℝ → ℂ) (x : ℝ) : ℂ :=
  (Ioc 0 ℓ).indicator (qE p p') x + (Ioi ℓ).indicator (tE ℓ p) x

theorem tE_eq (ℓ : ℝ) (p : ℝ → ℂ) (x : ℝ) :
    tE ℓ p x = (conj (p ℓ) * cexp ((ℓ : ℂ) / 2)) * cexp ((-(1 / 2) : ℂ) * x) := by
  unfold tE
  rw [mul_assoc, ← Complex.exp_add]
  congr 2; ring

theorem integrableOn_tE (ℓ : ℝ) (p : ℝ → ℂ) : IntegrableOn (tE ℓ p) (Ioi ℓ) := by
  have h := (integrableOn_exp_mul_complex_Ioi (a := (-(1 / 2) : ℂ)) (by norm_num) ℓ).const_mul
    (conj (p ℓ) * cexp ((ℓ : ℂ) / 2))
  exact MeasureTheory.IntegrableOn.congr_fun h (fun x _ => (tE_eq ℓ p x).symm) measurableSet_Ioi

theorem norm_tE_sq (ℓ : ℝ) (p : ℝ → ℂ) (x : ℝ) : ‖tE ℓ p x‖ ^ 2 = ‖p ℓ‖ ^ 2 * Real.exp (-(x - ℓ)) := by
  unfold tE
  rw [norm_mul, Complex.norm_conj, Complex.norm_exp, mul_pow, sq (Real.exp _), ← Real.exp_add]
  congr 2
  simp

theorem integral_norm_tE_sq (ℓ : ℝ) (p : ℝ → ℂ) : ∫ x in Ioi ℓ, ‖tE ℓ p x‖ ^ 2 = ‖p ℓ‖ ^ 2 := by
  simp_rw [norm_tE_sq]
  rw [MeasureTheory.integral_const_mul]
  have e : ∀ x : ℝ, Real.exp (-(x - ℓ)) = Real.exp ℓ * Real.exp (-1 * x) := by
    intro x; rw [← Real.exp_add]; congr 1; ring
  simp_rw [e]
  rw [MeasureTheory.integral_const_mul, integral_exp_mul_Ioi (by norm_num : (-1 : ℝ) < 0)]
  rw [neg_div_neg_eq, div_one, neg_one_mul, ← Real.exp_add, add_neg_cancel, Real.exp_zero, mul_one]

theorem integrableOn_norm_tE_sq (ℓ : ℝ) (p : ℝ → ℂ) : IntegrableOn (fun x => ‖tE ℓ p x‖ ^ 2) (Ioi ℓ) := by
  have h := (integrableOn_exp_mul_Ioi (a := -1) (by norm_num) ℓ).const_mul (‖p ℓ‖ ^ 2 * Real.exp ℓ)
  refine MeasureTheory.IntegrableOn.congr_fun h (fun x _ => ?_) measurableSet_Ioi
  rw [norm_tE_sq, mul_assoc, ← Real.exp_add]
  congr 2; ring

theorem norm_qE_sq (p p' : ℝ → ℂ) (x : ℝ) :
    ‖qE p p' x‖ ^ 2 = ‖p' x‖ ^ 2 + 1 / 4 * ‖p x‖ ^ 2 - ⟪p x, p' x⟫_ℝ := by
  unfold qE
  rw [Complex.norm_conj, Complex.inner, Complex.sq_norm, Complex.sq_norm, Complex.sq_norm]
  simp [Complex.normSq_apply]
  ring

theorem integral_inner_deriv (ℓ : ℝ) (p p' : ℝ → ℂ) (hp : ∀ x, HasDerivAt p (p' x) x) (hcp : Continuous p)
    (hcp' : Continuous p') :
    ∫ x in (0 : ℝ)..ℓ, ⟪p x, p' x⟫_ℝ = (‖p ℓ‖ ^ 2 - ‖p 0‖ ^ 2) / 2 := by
  have h := intervalIntegral.integral_eq_sub_of_hasDerivAt (f := fun x => ‖p x‖ ^ 2)
    (f' := fun x => 2 * ⟪p x, p' x⟫_ℝ) (a := 0) (b := ℓ) (fun x _ => (hp x).norm_sq)
    ((continuous_const.mul (hcp.inner hcp')).intervalIntegrable _ _)
  rw [intervalIntegral.integral_const_mul] at h
  linarith

theorem integral_norm_qE_sq (ℓ : ℝ) (p p' : ℝ → ℂ) (hp : ∀ x, HasDerivAt p (p' x) x) (hcp : Continuous p)
    (hcp' : Continuous p') :
    ∫ x in (0 : ℝ)..ℓ, ‖qE p p' x‖ ^ 2 =
      (∫ x in (0 : ℝ)..ℓ, (‖p' x‖ ^ 2 + 1 / 4 * ‖p x‖ ^ 2)) - (‖p ℓ‖ ^ 2 - ‖p 0‖ ^ 2) / 2 := by
  simp_rw [norm_qE_sq]
  rw [intervalIntegral.integral_sub, integral_inner_deriv ℓ p p' hp hcp hcp']
  · exact ((hcp'.norm.pow 2).add (continuous_const.mul (hcp.norm.pow 2))).intervalIntegrable _ _
  · exact (hcp.inner hcp').intervalIntegrable _ _

theorem eS_line (t x : ℝ) : eS (1 / 2 + t * I) x = cexp (-(x * t * I)) := by
  unfold eS; congr 1; ring

theorem fourier_tE (ℓ : ℝ) (p : ℝ → ℂ) (t : ℝ) :
    ∫ x in Ioi ℓ, tE ℓ p x * eS (1 / 2 + t * I) x = conj (p ℓ) * eS (1 / 2 + t * I) ℓ / (1 / 2 + t * I) := by
  have ha : (-(1 / 2 + t * I) : ℂ).re < 0 := by norm_num
  have h := integral_exp_mul_complex_Ioi ha ℓ
  have e : ∀ x : ℝ, tE ℓ p x * eS (1 / 2 + t * I) x =
      (conj (p ℓ) * cexp ((ℓ : ℂ) / 2)) * cexp (-(1 / 2 + t * I) * x) := by
    intro x
    unfold tE eS
    rw [mul_assoc, ← Complex.exp_add, mul_assoc, ← Complex.exp_add]
    congr 2; ring
  simp_rw [e]
  rw [MeasureTheory.integral_const_mul, h]
  have e2 : cexp ((ℓ : ℂ) / 2) * cexp (-(1 / 2 + t * I) * ℓ) = eS (1 / 2 + t * I) ℓ := by
    unfold eS; rw [← Complex.exp_add]; congr 1; ring
  rw [neg_div_neg_eq, ← e2]
  ring

theorem integrableOn_tE_mul (ℓ : ℝ) (p : ℝ → ℂ) (t : ℝ) :
    IntegrableOn (fun x => tE ℓ p x * eS (1 / 2 + t * I) x) (Ioi ℓ) := by
  have ha : (-(1 / 2 + t * I) : ℂ).re < 0 := by norm_num
  have h := (integrableOn_exp_mul_complex_Ioi ha ℓ).const_mul (conj (p ℓ) * cexp ((ℓ : ℂ) / 2))
  refine MeasureTheory.IntegrableOn.congr_fun h (fun x _ => ?_) measurableSet_Ioi
  unfold tE eS
  rw [mul_assoc, ← Complex.exp_add, mul_assoc, ← Complex.exp_add]
  congr 2; ring

theorem fourier_qE (ℓ : ℝ) (p p' p'' : ℝ → ℂ) (hp : ∀ x, HasDerivAt p (p' x) x)
    (hp' : ∀ x, HasDerivAt p' (p'' x) x) (hcp : Continuous p) (hcp' : Continuous p') (hcp'' : Continuous p'')
    (s : ℂ) :
    s * ∫ x in (0 : ℝ)..ℓ, qE p p' x * eS s x =
      (∫ x in (0 : ℝ)..ℓ, cexp (-(s - 1 / 2) * x) * conj (gDens p p' p'' x)) + conj (gA p p') -
        qE p p' ℓ * eS s ℓ := by
  obtain ⟨q', hq'def⟩ : ∃ q' : ℝ → ℂ, q' = fun x => conj ((1 / 2 : ℂ) * p' x - p'' x) := ⟨_, rfl⟩
  have hq : ∀ x, HasDerivAt (qE p p') (q' x) x := fun x => by
    rw [hq'def]
    have h := ((hp x).const_mul (1 / 2 : ℂ)).sub (hp' x)
    exact hasDerivAt_conj (g := fun y => (1 / 2 : ℂ) * p y - p' y) (g' := fun y => (1 / 2 : ℂ) * p' y - p'' y) h
  have hcq : Continuous (qE p p') := by unfold qE; fun_prop
  have hcq' : Continuous q' := by rw [hq'def]; fun_prop
  have hce : Continuous (eS s) := continuous_eS s
  have hce' : Continuous (eS' s) := continuous_eS' s
  have hftc := intervalIntegral.integral_eq_sub_of_hasDerivAt (a := 0) (b := ℓ)
    (fun x _ => (hq x).mul (hasDerivAt_eS s x)) (((hcq'.mul hce).add (hcq.mul hce')).intervalIntegrable 0 ℓ)
  have hsplit : ∫ y in (0 : ℝ)..ℓ, (q' y * eS s y + qE p p' y * eS' s y) =
      (∫ y in (0 : ℝ)..ℓ, q' y * eS s y) + ∫ y in (0 : ℝ)..ℓ, qE p p' y * eS' s y :=
    intervalIntegral.integral_add ((hcq'.mul hce).intervalIntegrable _ _) ((hcq.mul hce').intervalIntegrable _ _)
  rw [hsplit] at hftc
  simp only [Pi.mul_apply] at hftc
  have h2 : ∫ x in (0 : ℝ)..ℓ, qE p p' x * eS' s x = -(s - 1 / 2) * ∫ x in (0 : ℝ)..ℓ, qE p p' x * eS s x := by
    rw [← intervalIntegral.integral_const_mul]
    refine intervalIntegral.integral_congr fun x _ => ?_
    simp only [eS', eS]; ring
  have h3 : ∫ x in (0 : ℝ)..ℓ, cexp (-(s - 1 / 2) * x) * conj (gDens p p' p'' x) =
      (1 / 2 : ℂ) * (∫ x in (0 : ℝ)..ℓ, qE p p' x * eS s x) + ∫ x in (0 : ℝ)..ℓ, q' x * eS s x := by
    have hi1 : IntervalIntegrable (fun x => (1 / 2 : ℂ) * (qE p p' x * eS s x)) volume 0 ℓ :=
      (continuous_const.mul (hcq.mul hce)).intervalIntegrable _ _
    have hi2 : IntervalIntegrable (fun x => q' x * eS s x) volume 0 ℓ := (hcq'.mul hce).intervalIntegrable _ _
    rw [← intervalIntegral.integral_const_mul, ← intervalIntegral.integral_add hi1 hi2]
    refine intervalIntegral.integral_congr fun x _ => ?_
    rw [hq'def]
    simp only [qE, gDens, eS, map_add, map_sub, map_mul, map_neg, map_div₀, map_one, map_ofNat]
    ring
  have h0 : eS s 0 = 1 := by simp [eS]
  have hA : qE p p' 0 = conj (gA p p') := by simp [qE, gA]
  rw [h2, h0, mul_one, hA] at hftc
  rw [h3]
  linear_combination (-1 : ℂ) * hftc

theorem line_ne_zero (t : ℝ) : (1 / 2 + t * I : ℂ) ≠ 0 := by
  intro h0; have := congrArg Complex.re h0; norm_num at this

theorem fourier_hE {ℓ : ℝ} (hℓ : 0 ≤ ℓ) (p p' p'' : ℝ → ℂ) (hp : ∀ x, HasDerivAt p (p' x) x)
    (hp' : ∀ x, HasDerivAt p' (p'' x) x) (hcp : Continuous p) (hcp' : Continuous p') (hcp'' : Continuous p'')
    (t : ℝ) : ∫ x : ℝ, hE ℓ p p' x * cexp (-(x * t * I)) = Elow ℓ p p' (1 / 2 + t * I) / (1 / 2 + t * I) := by
  obtain ⟨s, hsdef⟩ : ∃ s : ℂ, s = 1 / 2 + t * I := ⟨_, rfl⟩
  have hs : s ≠ 0 := by rw [hsdef]; exact line_ne_zero t
  have hcq : Continuous (qE p p') := by unfold qE; fun_prop
  have e1 : ∀ x : ℝ, hE ℓ p p' x * cexp (-(x * t * I)) = (Ioc 0 ℓ).indicator (fun x => qE p p' x * eS s x) x +
      (Ioi ℓ).indicator (fun x => tE ℓ p x * eS s x) x := by
    intro x
    unfold hE
    rw [add_mul, ← eS_line, ← hsdef, Set.indicator_mul_left, Set.indicator_mul_left]
  simp_rw [e1]
  have i1 : Integrable ((Ioc 0 ℓ).indicator (fun x => qE p p' x * eS s x)) :=
    (((hcq.mul (continuous_eS s)).integrableOn_Icc).mono_set Ioc_subset_Icc_self).integrable_indicator
      measurableSet_Ioc
  have i2 : Integrable ((Ioi ℓ).indicator (fun x => tE ℓ p x * eS s x)) := by
    rw [hsdef]; exact (integrableOn_tE_mul ℓ p t).integrable_indicator measurableSet_Ioi
  rw [integral_add i1 i2, integral_indicator measurableSet_Ioc, integral_indicator measurableSet_Ioi,
    ← intervalIntegral.integral_of_le hℓ, hsdef, fourier_tE, ← hsdef]
  rw [Elow_green p p' p'' hp' hcp hcp' hcp'' s]
  have hf := fourier_qE ℓ p p' p'' hp hp' hcp hcp' hcp'' s
  have hB : conj (gB ℓ p p') = conj (p ℓ) - qE p p' ℓ := by
    simp only [gB, qE, map_add, map_sub, map_mul, map_div₀, map_one, map_ofNat]; ring
  have heℓ : cexp (-(s - 1 / 2) * ℓ) = eS s ℓ := rfl
  rw [heℓ, hB, eq_div_iff hs, add_mul, div_mul_cancel₀ _ hs]
  linear_combination hf

theorem norm_hE_sq (ℓ : ℝ) (p p' : ℝ → ℂ) (x : ℝ) :
    ‖hE ℓ p p' x‖ ^ 2 = (Ioc 0 ℓ).indicator (fun x => ‖qE p p' x‖ ^ 2) x +
      (Ioi ℓ).indicator (fun x => ‖tE ℓ p x‖ ^ 2) x := by
  unfold hE
  by_cases h1 : x ∈ Ioc 0 ℓ
  · have h2 : x ∉ Ioi ℓ := fun h => by
      have a1 : x ≤ ℓ := h1.2
      have a2 : ℓ < x := h
      linarith
    rw [Set.indicator_of_mem h1, Set.indicator_of_notMem h2, Set.indicator_of_mem h1,
      Set.indicator_of_notMem h2, add_zero, add_zero]
  · by_cases h2 : x ∈ Ioi ℓ
    · rw [Set.indicator_of_mem h2, Set.indicator_of_notMem h1, Set.indicator_of_mem h2,
        Set.indicator_of_notMem h1, zero_add, zero_add]
    · rw [Set.indicator_of_notMem h1, Set.indicator_of_notMem h2, Set.indicator_of_notMem h1,
        Set.indicator_of_notMem h2]
      simp

/-- **Plancherel on the line**: E_low(s)/s is in L², with squared norm 2π N(p). -/
theorem line_Elow (hP : Rudin1987_Thm9_13) {ℓ : ℝ} (hℓ : 0 ≤ ℓ) (p p' p'' : ℝ → ℂ)
    (hp : ∀ x, HasDerivAt p (p' x) x) (hp' : ∀ x, HasDerivAt p' (p'' x) x) (hcp : Continuous p)
    (hcp' : Continuous p') (hcp'' : Continuous p'') :
    MemLp (fun t : ℝ => Elow ℓ p p' (1 / 2 + t * I) / (1 / 2 + t * I)) 2 ∧
    ∫ t : ℝ, ‖Elow ℓ p p' (1 / 2 + t * I) / (1 / 2 + t * I)‖ ^ 2 = 2 * π * Nself ℓ p p' := by
  have hcq : Continuous (qE p p') := by unfold qE; fun_prop
  have iq : IntegrableOn (qE p p') (Ioc 0 ℓ) := (hcq.integrableOn_Icc).mono_set Ioc_subset_Icc_self
  have iq2 : IntegrableOn (fun x => ‖qE p p' x‖ ^ 2) (Ioc 0 ℓ) :=
    ((hcq.norm.pow 2).integrableOn_Icc).mono_set Ioc_subset_Icc_self
  have hEdef : hE ℓ p p' = (Ioc 0 ℓ).indicator (qE p p') + (Ioi ℓ).indicator (tE ℓ p) := rfl
  have hi : Integrable (hE ℓ p p') := by
    rw [hEdef]
    exact (iq.integrable_indicator measurableSet_Ioc).add
      ((integrableOn_tE ℓ p).integrable_indicator measurableSet_Ioi)
  have e2 : (fun x => ‖hE ℓ p p' x‖ ^ 2) = (Ioc 0 ℓ).indicator (fun x => ‖qE p p' x‖ ^ 2) +
      (Ioi ℓ).indicator (fun x => ‖tE ℓ p x‖ ^ 2) := funext fun x => norm_hE_sq ℓ p p' x
  have i1 : Integrable ((Ioc 0 ℓ).indicator (fun x => ‖qE p p' x‖ ^ 2)) :=
    iq2.integrable_indicator measurableSet_Ioc
  have i2 : Integrable ((Ioi ℓ).indicator (fun x => ‖tE ℓ p x‖ ^ 2)) :=
    (integrableOn_norm_tE_sq ℓ p).integrable_indicator measurableSet_Ioi
  have hi2 : Integrable (fun x => ‖hE ℓ p p' x‖ ^ 2) := by rw [e2]; exact i1.add i2
  have hm : MemLp (hE ℓ p p') 2 := (memLp_two_iff_integrable_sq_norm hi.aestronglyMeasurable).mpr hi2
  obtain ⟨h1, h2⟩ := hP (hE ℓ p p') hi hm
  simp_rw [fourier_hE hℓ p p' p'' hp hp' hcp hcp' hcp''] at h1 h2
  refine ⟨h1, ?_⟩
  rw [h2]
  simp_rw [norm_hE_sq ℓ p p']
  rw [integral_add i1 i2, integral_indicator measurableSet_Ioc, integral_indicator measurableSet_Ioi,
    ← intervalIntegral.integral_of_le hℓ, integral_norm_qE_sq ℓ p p' hp hcp hcp', integral_norm_tE_sq]
  unfold Nself
  ring

end Riemann.Deep.LeanNBRate

