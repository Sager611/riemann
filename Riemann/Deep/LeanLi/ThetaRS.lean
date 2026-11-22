module

public import Riemann.Deep.LeanLi.CountBounds
public import Mathlib.Analysis.Complex.CauchyIntegral
public import Mathlib.Analysis.SpecialFunctions.Integrals.Basic

/-!
# LeanLi: `A(t) = θ_RS(t)/π + 1`, `A'(t) = (2π)⁻¹ log(t/2π) + r(t)`, `|r(t)| ≤ C/t²`, `∫₀^∞ r = -1/8`
-/

@[expose] public section

open Complex Real Filter Topology Set MeasureTheory

namespace LeanLi

lemma differentiableOn_Gamma_re_pos : DifferentiableOn ℂ Complex.Gamma {s : ℂ | 0 < s.re} := by
  intro s hs
  refine (differentiableAt_Gamma s fun m h => ?_).differentiableWithinAt
  have : s.re = -(m : ℝ) := by rw [h]; simp
  have hm : (0 : ℝ) ≤ m := Nat.cast_nonneg m
  simp only [Set.mem_setOf_eq] at hs
  linarith

lemma continuousOn_digamma : ContinuousOn digamma {s : ℂ | 0 < s.re} := by
  have hU : IsOpen {s : ℂ | 0 < s.re} := isOpen_lt continuous_const Complex.continuous_re
  have hd := (differentiableOn_Gamma_re_pos.deriv hU).continuousOn
  have hg := differentiableOn_Gamma_re_pos.continuousOn
  have hne : ∀ s ∈ {s : ℂ | 0 < s.re}, Complex.Gamma s ≠ 0 := fun s hs => Gamma_ne_zero_of_re_pos hs
  refine (hd.div hg hne).congr fun s _ => ?_
  rw [digamma_def, logDeriv_apply]; rfl

/-- The integrand `θ_RS'(u) = ½ Re ψ(1/4 + iu/2) - ½ log π`. -/
noncomputable def thetaD (u : ℝ) : ℝ := (digamma (1 / 4 + (u / 2 : ℝ) * I)).re / 2 - Real.log Real.pi / 2

lemma continuous_thetaD : Continuous thetaD := by
  have h1 : Continuous fun u : ℝ => (1 / 4 + (u / 2 : ℝ) * I : ℂ) := by fun_prop
  have h2 : Continuous fun u : ℝ => digamma (1 / 4 + (u / 2 : ℝ) * I) := by
    refine continuousOn_digamma.comp_continuous h1 fun u => ?_
    simp only [Set.mem_setOf_eq, add_re, mul_re, ofReal_re, I_re, mul_zero, ofReal_im, I_im,
      mul_one, sub_zero]
    norm_num
  unfold thetaD
  fun_prop

lemma rsTheta_eq (t : ℝ) : rsTheta t = ∫ u in (0 : ℝ)..t, thetaD u := rfl

/-- `A(t) = θ_RS(t)/π + 1`. -/
noncomputable def Afun (t : ℝ) : ℝ := rsTheta t / Real.pi + 1

lemma Afun_zero : Afun 0 = 1 := by simp [Afun, rsTheta]

/-- `A'(t) = θ_RS'(t)/π`. -/
noncomputable def Ad (t : ℝ) : ℝ := thetaD t / Real.pi

lemma hasDerivAt_Afun (t : ℝ) : HasDerivAt Afun (Ad t) t := by
  have h := (continuous_thetaD.integral_hasStrictDerivAt 0 t).hasDerivAt
  unfold Afun Ad
  exact (h.div_const Real.pi).add_const 1

lemma continuous_Ad : Continuous Ad := continuous_thetaD.div_const _

lemma integral_Ad (a b : ℝ) : ∫ t in a..b, Ad t = Afun b - Afun a :=
  intervalIntegral.integral_eq_sub_of_hasDerivAt (fun t _ => hasDerivAt_Afun t)
    (continuous_Ad.intervalIntegrable _ _)

/-- `r(t) = A'(t) - (2π)⁻¹ log(t/2π)`. -/
noncomputable def rfun (t : ℝ) : ℝ := Ad t - 1 / (2 * Real.pi) * Real.log (t / (2 * Real.pi))

/-- **`|r(t)| ≤ C/t²`** for large `t` (DLMF 5.11.2). -/
lemma rfun_bound (L : LiteratureInputs) : ∃ C t₀ : ℝ, 1 ≤ t₀ ∧ ∀ t, t₀ ≤ t → |rfun t| ≤ C / t ^ 2 := by
  obtain ⟨C, R₀, hC⟩ := L.digamma
  refine ⟨(|C| * 4 + 1) / Real.pi, max 1 (2 * |R₀| + 2), le_max_left _ _, fun t ht => ?_⟩
  have ht1 : 1 ≤ t := le_trans (le_max_left _ _) ht
  have htR : 2 * |R₀| + 2 ≤ t := le_trans (le_max_right _ _) ht
  have htpos : 0 < t := by linarith
  have hπ := Real.pi_pos
  set z : ℂ := 1 / 4 + (t / 2 : ℝ) * I with hz
  have hzre : z.re = 1 / 4 := by simp [hz]
  have hzim : z.im = t / 2 := by simp [hz]
  have hnsq : Complex.normSq z = 1 / 16 + t ^ 2 / 4 := by
    rw [Complex.normSq_apply, hzre, hzim]; ring
  have hnorm : ‖z‖ ^ 2 = 1 / 16 + t ^ 2 / 4 := by rw [← Complex.normSq_eq_norm_sq, hnsq]
  have hnorm_ge : t / 2 ≤ ‖z‖ := by
    have := Complex.abs_im_le_norm z
    rw [hzim, abs_of_pos (by positivity)] at this; exact this
  have hz0 : 0 < ‖z‖ := lt_of_lt_of_le (by positivity) hnorm_ge
  have hR : R₀ ≤ ‖z‖ := by linarith [le_abs_self R₀]
  have h := hC z hR (by rw [hzre]; norm_num)
  -- real parts
  have hre := (Complex.abs_re_le_norm (digamma z - (Complex.log z - 1 / (2 * z)))).trans h
  rw [Complex.sub_re, Complex.sub_re, Complex.log_re] at hre
  have hinv : (1 / (2 * z)).re = 1 / (8 * (1 / 16 + t ^ 2 / 4)) := by
    rw [one_div, mul_inv, Complex.mul_re, Complex.inv_re, Complex.inv_re, Complex.inv_im]
    simp [hzre, hnsq]
    ring
  rw [hinv] at hre
  -- log ‖z‖ = log(t/2) + small
  have hlog : |Real.log ‖z‖ - Real.log (t / 2)| ≤ 1 / (8 * t ^ 2) := by
    have e : Real.log ‖z‖ = Real.log (‖z‖ ^ 2) / 2 := by
      rw [Real.log_pow]; push_cast; ring
    have e2 : ‖z‖ ^ 2 = (t / 2) ^ 2 * (1 + 1 / (4 * t ^ 2)) := by rw [hnorm]; field_simp; ring
    rw [e, e2, Real.log_mul (by positivity) (by positivity), Real.log_pow]
    push_cast
    have h1 := Real.log_le_sub_one_of_pos (show 0 < 1 + 1 / (4 * t ^ 2) by positivity)
    have hx0 : 0 ≤ 1 / (4 * t ^ 2) := by positivity
    have h2 := Real.log_nonneg (show 1 ≤ 1 + 1 / (4 * t ^ 2) by linarith)
    have e3 : 1 / (8 * t ^ 2) = 1 / (4 * t ^ 2) / 2 := by field_simp; ring
    have e4 : (2 * Real.log (t / 2) + Real.log (1 + 1 / (4 * t ^ 2))) / 2 - Real.log (t / 2) =
        Real.log (1 + 1 / (4 * t ^ 2)) / 2 := by ring
    rw [e4, e3, abs_le]; constructor <;> linarith
  -- assemble: r = (Re ψ(z)/2 - log π/2)/π - (1/2π) log(t/2π)
  have hr : rfun t = ((digamma z).re - Real.log (t / 2)) / (2 * Real.pi) := by
    unfold rfun Ad thetaD
    rw [← hz]
    rw [Real.log_div (by positivity) (by positivity), Real.log_div (by positivity) (by norm_num),
      Real.log_mul (by norm_num) hπ.ne']
    field_simp
    ring
  rw [hr, abs_div, abs_of_pos (by positivity : (0 : ℝ) < 2 * Real.pi)]
  have hz2 : C / ‖z‖ ^ 2 ≤ |C| * 4 / t ^ 2 := by
    rw [hnorm]
    rw [div_le_div_iff₀ (by positivity) (by positivity)]
    have := le_abs_self C
    have h0 : 0 ≤ |C| := abs_nonneg C
    nlinarith
  have hsmall : 1 / (8 * (1 / 16 + t ^ 2 / 4)) ≤ 1 / (2 * t ^ 2) := by
    apply one_div_le_one_div_of_le (by positivity); nlinarith
  have hfin : |(digamma z).re - Real.log (t / 2)| ≤ (|C| * 4 + 1) / t ^ 2 := by
    have h1 := abs_le.1 hre
    have h2 := abs_le.1 hlog
    have hp : 0 ≤ 1 / (8 * (1 / 16 + t ^ 2 / 4)) := by positivity
    have e1 : 1 / (8 * t ^ 2) ≤ 1 / (2 * t ^ 2) := by
      apply one_div_le_one_div_of_le (by positivity); nlinarith
    have e3 : (|C| * 4 + 1) / t ^ 2 = |C| * 4 / t ^ 2 + 1 / (2 * t ^ 2) + 1 / (2 * t ^ 2) := by
      field_simp; ring
    rw [abs_le, e3]; constructor <;> linarith
  calc |(digamma z).re - Real.log (t / 2)| / (2 * Real.pi) ≤ (|C| * 4 + 1) / t ^ 2 / (2 * Real.pi) :=
        div_le_div_of_nonneg_right hfin (by positivity)
    _ ≤ (|C| * 4 + 1) / Real.pi / t ^ 2 := by
        rw [div_div, div_div, div_le_div_iff₀ (by positivity) (by positivity)]
        have : 0 ≤ (|C| * 4 + 1) * (π * t ^ 2) := by positivity
        linarith

end LeanLi
