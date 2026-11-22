import Mathlib.Analysis.Complex.PhragmenLindelof
import Mathlib.Analysis.SpecialFunctions.Gamma.Beta
import Mathlib.Analysis.SpecialFunctions.Gamma.BohrMollerup

/-!
# Bounds for the Gamma function on vertical strips (Phragmen-Lindelof from the reflection formula)

|Gamma(1/2+it)|^2 = pi/cosh(pi t), and on 1/2 <= Re z <= 3/2: |Gamma z| <= sqrt(2 pi) |z+1| exp(-pi |Im z|/2).
-/

open Complex Filter Topology Set MeasureTheory
open scoped Real

noncomputable section

namespace Riemann.Deep.LeanNBRate

theorem norm_Gamma_le_Gamma_re {z : ℂ} (hz : 0 < z.re) : ‖Complex.Gamma z‖ ≤ Real.Gamma z.re := by
  rw [Complex.Gamma_eq_integral hz, Real.Gamma_eq_integral hz, Complex.GammaIntegral]
  refine (norm_integral_le_integral_norm _).trans (le_of_eq ?_)
  refine setIntegral_congr_fun measurableSet_Ioi fun x hx => ?_
  simp only [norm_mul, Complex.norm_real, Real.norm_eq_abs, abs_of_pos (Real.exp_pos _)]
  rw [Complex.norm_cpow_eq_rpow_re_of_pos hx]
  simp

theorem Gamma_re_le_two {σ : ℝ} (h1 : 1 / 2 ≤ σ) (h2 : σ ≤ 3 / 2) : Real.Gamma σ ≤ 2 := by
  have h := Real.convexOn_Gamma.le_max_of_mem_Icc (x := 1 / 2) (y := 3 / 2) (z := σ)
    (by norm_num) (by norm_num) ⟨h1, h2⟩
  have ha : Real.Gamma (1 / 2) ≤ 2 := by
    rw [Real.Gamma_one_half_eq]
    rw [show (2 : ℝ) = Real.sqrt 4 by rw [show (4 : ℝ) = 2 ^ 2 by norm_num, Real.sqrt_sq (by norm_num)]]
    exact Real.sqrt_le_sqrt (by linarith [Real.pi_le_four])
  have hb : Real.Gamma (3 / 2) ≤ 2 := by linarith [Real.Gamma_three_div_two_lt_one]
  exact h.trans (max_le ha hb)

theorem norm_Gamma_le_two {z : ℂ} (h1 : 1 / 2 ≤ z.re) (h2 : z.re ≤ 3 / 2) : ‖Complex.Gamma z‖ ≤ 2 :=
  (norm_Gamma_le_Gamma_re (by linarith)).trans (Gamma_re_le_two h1 h2)

theorem sin_pi_half_add (t : ℝ) : Complex.sin (π * (1 / 2 + t * I)) = (Real.cosh (π * t) : ℂ) := by
  rw [show (π : ℂ) * (1 / 2 + t * I) = ((π * t : ℝ) : ℂ) * I + π / 2 by push_cast; ring,
    Complex.sin_add_pi_div_two, Complex.cos_mul_I, Complex.ofReal_cosh]

theorem norm_Gamma_half_sq (t : ℝ) : ‖Complex.Gamma (1 / 2 + t * I)‖ ^ 2 = π / Real.cosh (π * t) := by
  have h := Complex.Gamma_mul_Gamma_one_sub (1 / 2 + t * I)
  have h1 : (1 : ℂ) - (1 / 2 + t * I) = (starRingEnd ℂ) (1 / 2 + t * I) := by
    apply Complex.ext <;> simp <;> norm_num
  rw [h1, Complex.Gamma_conj, Complex.mul_conj, Complex.normSq_eq_norm_sq, sin_pi_half_add] at h
  have hc : (0 : ℝ) < Real.cosh (π * t) := Real.cosh_pos _
  have : ((‖Complex.Gamma (1 / 2 + t * I)‖ ^ 2 : ℝ) : ℂ) = ((π / Real.cosh (π * t) : ℝ) : ℂ) := by
    rw [h]; push_cast; ring
  exact_mod_cast this

theorem norm_Gamma_half_mul_exp_le (t : ℝ) :
    ‖Complex.Gamma (1 / 2 + t * I)‖ * Real.exp (π / 2 * t) ≤ Real.sqrt (2 * π) := by
  apply Real.le_sqrt_of_sq_le
  rw [mul_pow, norm_Gamma_half_sq, ← Real.exp_nat_mul]
  have hc : (0 : ℝ) < Real.cosh (π * t) := Real.cosh_pos _
  rw [div_mul_eq_mul_div, div_le_iff₀ hc, Real.cosh_eq]
  have e : ((2 : ℕ) : ℝ) * (π / 2 * t) = π * t := by push_cast; ring
  rw [e]
  have : 0 < Real.exp (-(π * t)) := Real.exp_pos _
  nlinarith [Real.pi_pos]

theorem norm_exp_neg_pi_half_I (z : ℂ) : ‖cexp (-(π / 2 : ℂ) * I * z)‖ = Real.exp (π / 2 * z.im) := by
  rw [Complex.norm_exp]
  congr 1
  simp [Complex.mul_re, Complex.mul_im]

/-- Phragmen-Lindelof for Gamma(z) e^{-i pi z/2}/(z+1) on the strip 1/2 <= Re z <= 3/2. -/
theorem norm_Gamma_strip_aux {z : ℂ} (h1 : 1 / 2 ≤ z.re) (h2 : z.re ≤ 3 / 2) :
    ‖Complex.Gamma z * cexp (-(π / 2 : ℂ) * I * z) / (z + 1)‖ ≤ Real.sqrt (2 * π) := by
  set f : ℂ → ℂ := fun z => Complex.Gamma z * cexp (-(π / 2 : ℂ) * I * z) / (z + 1) with hf
  have hne : ∀ w : ℂ, 0 < w.re → w + 1 ≠ 0 := fun w hw h => by
    have := congrArg Complex.re h; simp at this; linarith
  have hd : DifferentiableOn ℂ f (re ⁻¹' Ioi 0) := by
    intro w hw
    have hw' : (0 : ℝ) < w.re := hw
    have hG : DifferentiableAt ℂ Complex.Gamma w := by
      refine Complex.differentiableAt_Gamma _ fun m hm => ?_
      have := congrArg Complex.re hm; simp at this; linarith [(Nat.cast_nonneg m : (0 : ℝ) ≤ m)]
    exact ((hG.mul (by fun_prop)).div (by fun_prop) (hne w hw')).differentiableWithinAt
  have hsub : closure (re ⁻¹' Ioo (1 / 2 : ℝ) (3 / 2)) ⊆ re ⁻¹' Ioi 0 := by
    rw [closure_preimage_re, closure_Ioo (by norm_num)]
    intro w hw; simp only [mem_preimage, mem_Icc] at hw; simp only [mem_preimage, mem_Ioi]; linarith [hw.1]
  have hdc : DiffContOnCl ℂ f (re ⁻¹' Ioo (1 / 2 : ℝ) (3 / 2)) :=
    (hd.mono hsub).diffContOnCl
  have hbd : ∀ w : ℂ, 1 / 2 ≤ w.re → w.re ≤ 3 / 2 → ‖f w‖ ≤ 2 * Real.exp (π / 2 * |w.im|) := by
    intro w hw1 hw2
    simp only [hf, norm_div, norm_mul, norm_exp_neg_pi_half_I]
    have hn : 1 ≤ ‖w + 1‖ := by
      have := Complex.abs_re_le_norm (w + 1)
      simp at this
      have h3 : |w.re + 1| = w.re + 1 := abs_of_pos (by linarith)
      linarith
    rw [div_le_iff₀ (by linarith)]
    have e1 := norm_Gamma_le_two hw1 hw2
    have e2 : Real.exp (π / 2 * w.im) ≤ Real.exp (π / 2 * |w.im|) :=
      Real.exp_le_exp.mpr (mul_le_mul_of_nonneg_left (le_abs_self _) (by positivity))
    calc ‖Complex.Gamma w‖ * Real.exp (π / 2 * w.im) ≤ 2 * Real.exp (π / 2 * |w.im|) :=
          mul_le_mul e1 e2 (Real.exp_pos _).le (by norm_num)
      _ ≤ 2 * Real.exp (π / 2 * |w.im|) * ‖w + 1‖ := le_mul_of_one_le_right (by positivity) hn
  refine PhragmenLindelof.vertical_strip (a := 1 / 2) (b := 3 / 2) (f := f) hdc ?_ ?_ ?_ h1 h2
  · refine ⟨1, by rw [show (3 / 2 : ℝ) - 1 / 2 = 1 by norm_num, div_one]; linarith [Real.two_le_pi], π / 2, ?_⟩
    refine Asymptotics.IsBigO.of_bound 2 ?_
    refine Filter.eventually_inf_principal.2 (Filter.Eventually.of_forall fun w hw => ?_)
    simp only [mem_preimage, mem_Ioo] at hw
    refine (hbd w hw.1.le hw.2.le).trans ?_
    rw [Real.norm_eq_abs, abs_of_pos (Real.exp_pos _)]
    refine mul_le_mul_of_nonneg_left (Real.exp_le_exp.mpr ?_) (by norm_num)
    refine mul_le_mul_of_nonneg_left ?_ (by positivity)
    rw [one_mul]
    linarith [Real.add_one_le_exp (|w.im|), abs_nonneg w.im]
  · intro w hw
    have hw' : w = 1 / 2 + (w.im : ℂ) * I := by
      apply Complex.ext <;> simp [hw]
    rw [hw']
    simp only [hf, norm_div, norm_mul, norm_exp_neg_pi_half_I]
    have hn : 1 ≤ ‖(1 / 2 + (w.im : ℂ) * I) + 1‖ := by
      have := Complex.abs_re_le_norm ((1 / 2 + (w.im : ℂ) * I) + 1)
      simp at this; norm_num at this; linarith
    have him : (1 / 2 + (w.im : ℂ) * I).im = w.im := by simp
    rw [him, div_le_iff₀ (by linarith)]
    exact (norm_Gamma_half_mul_exp_le w.im).trans (le_mul_of_one_le_right (Real.sqrt_nonneg _) hn)
  · intro w hw
    have hw' : w = (1 / 2 + (w.im : ℂ) * I) + 1 := by
      apply Complex.ext <;> simp [hw] <;> norm_num
    rw [hw']
    have hz0 : (1 / 2 + (w.im : ℂ) * I) ≠ 0 := by
      intro h; have := congrArg Complex.re h; simp at this
    simp only [hf, Complex.Gamma_add_one _ hz0, norm_div, norm_mul, norm_exp_neg_pi_half_I]
    have him : ((1 / 2 + (w.im : ℂ) * I) + 1).im = w.im := by simp
    rw [him]
    have hle : ‖(1 / 2 + (w.im : ℂ) * I)‖ ≤ ‖(1 / 2 + (w.im : ℂ) * I) + 1 + 1‖ := by
      rw [show (1 / 2 + (w.im : ℂ) * I) + 1 + 1 = (5 / 2 : ℝ) + (w.im : ℂ) * I by push_cast; ring,
        show (1 / 2 + (w.im : ℂ) * I) = (1 / 2 : ℝ) + (w.im : ℂ) * I by push_cast; ring]
      rw [Complex.norm_add_mul_I, Complex.norm_add_mul_I]
      exact Real.sqrt_le_sqrt (by nlinarith)
    have hpos : 0 < ‖(1 / 2 + (w.im : ℂ) * I) + 1 + 1‖ := by
      refine norm_pos_iff.mpr fun h => ?_
      have := congrArg Complex.re h; simp at this; norm_num at this
    rw [div_le_iff₀ hpos]
    calc ‖(1 / 2 + (w.im : ℂ) * I)‖ * ‖Complex.Gamma (1 / 2 + (w.im : ℂ) * I)‖ * Real.exp (π / 2 * w.im)
        = ‖(1 / 2 + (w.im : ℂ) * I)‖ * (‖Complex.Gamma (1 / 2 + (w.im : ℂ) * I)‖ *
          Real.exp (π / 2 * w.im)) := by ring
      _ ≤ ‖(1 / 2 + (w.im : ℂ) * I) + 1 + 1‖ * Real.sqrt (2 * π) :=
          mul_le_mul hle (norm_Gamma_half_mul_exp_le w.im) (by positivity) (norm_nonneg _)
      _ = Real.sqrt (2 * π) * ‖(1 / 2 + (w.im : ℂ) * I) + 1 + 1‖ := by ring

/-- |Gamma z| <= sqrt(2 pi) |z+1| exp(-pi |Im z|/2) on 1/2 <= Re z <= 3/2. -/
theorem norm_Gamma_strip {z : ℂ} (h1 : 1 / 2 ≤ z.re) (h2 : z.re ≤ 3 / 2) :
    ‖Complex.Gamma z‖ ≤ Real.sqrt (2 * π) * ‖z + 1‖ * Real.exp (-(π / 2 * |z.im|)) := by
  have key : ∀ w : ℂ, 1 / 2 ≤ w.re → w.re ≤ 3 / 2 → 0 ≤ w.im →
      ‖Complex.Gamma w‖ ≤ Real.sqrt (2 * π) * ‖w + 1‖ * Real.exp (-(π / 2 * |w.im|)) := by
    intro w hw1 hw2 hw3
    have h := norm_Gamma_strip_aux hw1 hw2
    have hn : 0 < ‖w + 1‖ := norm_pos_iff.mpr fun h => by
      have := congrArg Complex.re h; simp at this; linarith
    rw [norm_div, norm_mul, norm_exp_neg_pi_half_I, div_le_iff₀ hn] at h
    rw [abs_of_nonneg hw3, Real.exp_neg]
    rw [← div_eq_mul_inv, le_div_iff₀ (Real.exp_pos _)]
    exact h
  rcases le_total 0 z.im with h | h
  · exact key z h1 h2 h
  · have h' := key ((starRingEnd ℂ) z) (by simpa using h1) (by simpa using h2) (by simpa using h)
    rw [Complex.Gamma_conj, Complex.norm_conj] at h'
    have e1 : ‖(starRingEnd ℂ) z + 1‖ = ‖z + 1‖ := by
      rw [show (starRingEnd ℂ) z + 1 = (starRingEnd ℂ) (z + 1) by simp, Complex.norm_conj]
    rw [e1, Complex.conj_im, abs_neg] at h'
    exact h'

end Riemann.Deep.LeanNBRate

