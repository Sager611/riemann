import Riemann.Deep.LeanNBRate.Dual5
import Mathlib.Analysis.SpecialFunctions.Trigonometric.Bounds
import Mathlib.Analysis.SpecialFunctions.Integrals.Basic

/-!
# The window: M_a(x) = int_0^1 e^{-axv} dv and Phi(s) = M_a(s - 1/2)^n
-/

open Complex Filter Topology Set MeasureTheory
open scoped Real

noncomputable section

namespace Riemann.Deep.LeanNBRate

/-- M_a(x) = int_0^1 e^{-a x v} dv (entire, = (1 - e^{-ax})/(ax) for ax != 0). -/
def Mw (a : ℝ) (x : ℂ) : ℂ := ∫ v in (0 : ℝ)..1, cexp (-(a * x) * v)

theorem differentiable_Mw (a : ℝ) : Differentiable ℂ (Mw a) := by
  have h := differentiable_laplace (ℓ := 1) zero_le_one (q := fun _ => (1 : ℂ)) continuous_const
  have e : Mw a = (fun s : ℂ => ∫ u in (0 : ℝ)..1, cexp (-(s - 1 / 2) * u) * (fun _ => (1 : ℂ)) u) ∘
      (fun x : ℂ => (a : ℂ) * x + 1 / 2) := by
    funext x; simp only [Function.comp, Mw, mul_one]
    congr 1; funext v; congr 1; ring
  rw [e]
  exact h.comp (by fun_prop)

theorem norm_Mw_le {a : ℝ} (ha : 0 ≤ a) (x : ℂ) : ‖Mw a x‖ ≤ Real.exp (a * max 0 (-x.re)) := by
  unfold Mw
  have h := intervalIntegral.norm_integral_le_of_norm_le_const (a := (0 : ℝ)) (b := 1)
    (f := fun v : ℝ => cexp (-(a * x) * v)) (C := Real.exp (a * max 0 (-x.re))) ?_
  · simpa using h
  intro v hv
  rw [Set.uIoc_of_le zero_le_one] at hv
  rw [Complex.norm_exp]
  apply Real.exp_le_exp.mpr
  have hre : (-((a : ℂ) * x) * (v : ℂ)).re = a * (-x.re) * v := by
    simp only [neg_mul, Complex.neg_re, Complex.mul_re, Complex.ofReal_re, Complex.ofReal_im]; ring
  rw [hre]
  have h1 : -x.re ≤ max 0 (-x.re) := le_max_right _ _
  have h2 : 0 ≤ max 0 (-x.re) := le_max_left _ _
  have hv0 : 0 ≤ v := hv.1.le
  have hv1 : v ≤ 1 := hv.2
  calc a * -x.re * v ≤ a * max 0 (-x.re) * v :=
        mul_le_mul_of_nonneg_right (mul_le_mul_of_nonneg_left h1 ha) hv0
    _ ≤ a * max 0 (-x.re) * 1 := mul_le_mul_of_nonneg_left hv1 (mul_nonneg ha h2)
    _ = a * max 0 (-x.re) := mul_one _

theorem Mw_eq {a : ℝ} {x : ℂ} (hx : (a : ℂ) * x ≠ 0) : Mw a x = (1 - cexp (-(a * x))) / (a * x) := by
  unfold Mw
  have hc : -((a : ℂ) * x) ≠ 0 := neg_ne_zero.mpr hx
  rw [integral_exp_mul_complex hc]
  simp only [Complex.ofReal_one, Complex.ofReal_zero, mul_one, mul_zero, Complex.exp_zero]
  field_simp
  ring

theorem norm_Mw_le_inv {a : ℝ} (ha : 0 < a) {x : ℂ} (hx : x ≠ 0) :
    ‖Mw a x‖ ≤ (1 + Real.exp (-(a * x.re))) / (a * ‖x‖) := by
  have hax : (a : ℂ) * x ≠ 0 := mul_ne_zero (by exact_mod_cast ha.ne') hx
  rw [Mw_eq hax, norm_div, norm_mul, Complex.norm_real, Real.norm_eq_abs, abs_of_pos ha]
  gcongr
  calc ‖1 - cexp (-(a * x))‖ ≤ ‖(1 : ℂ)‖ + ‖cexp (-(a * x))‖ := norm_sub_le _ _
    _ = 1 + Real.exp (-(a * x.re)) := by
        rw [norm_one, Complex.norm_exp]; congr 2; simp

theorem norm_one_sub_Mw_I_le {a : ℝ} (ha : 0 ≤ a) (t : ℝ) : ‖1 - Mw a (t * I)‖ ≤ a * |t| / 2 := by
  unfold Mw
  have h1 : (1 : ℂ) = ∫ v in (0 : ℝ)..1, (1 : ℂ) := by simp
  rw [h1, ← intervalIntegral.integral_sub intervalIntegrable_const
    ((by fun_prop : Continuous fun v : ℝ => cexp (-((a : ℂ) * (t * I)) * v)).intervalIntegrable _ _)]
  have hb : ∀ v ∈ Set.uIoc (0 : ℝ) 1, ‖(1 : ℂ) - cexp (-((a : ℂ) * (t * I)) * v)‖ ≤ a * |t| * v := by
    intro v hv
    rw [Set.uIoc_of_le zero_le_one] at hv
    have e : -((a : ℂ) * (t * I)) * v = I * ((-(a * t * v) : ℝ) : ℂ) := by push_cast; ring
    rw [e, norm_sub_rev]
    refine (Real.norm_exp_I_mul_ofReal_sub_one_le).trans ?_
    rw [Real.norm_eq_abs, abs_neg, abs_mul, abs_mul, abs_of_nonneg ha, abs_of_pos hv.1]
  have := intervalIntegral.norm_integral_le_of_norm_le (μ := volume) (a := (0 : ℝ)) (b := 1)
    (f := fun v : ℝ => (1 : ℂ) - cexp (-((a : ℂ) * (t * I)) * v)) (g := fun v : ℝ => a * |t| * v)
    zero_le_one (Filter.Eventually.of_forall
    (fun v hv => hb v (by rw [Set.uIoc_of_le zero_le_one]; exact ⟨hv.1, hv.2⟩)))
    ((continuous_const.mul continuous_id).intervalIntegrable 0 1)
  refine this.trans (le_of_eq ?_)
  show ∫ v in (0 : ℝ)..1, a * |t| * v = a * |t| / 2
  rw [intervalIntegral.integral_const_mul, integral_id]
  ring

/-- The window Phi(s) = M_a(s - 1/2)^n. -/
def Phi (a : ℝ) (n : ℕ) (s : ℂ) : ℂ := Mw a (s - 1 / 2) ^ n

theorem differentiable_Phi (a : ℝ) (n : ℕ) : Differentiable ℂ (Phi a n) := by
  unfold Phi
  exact ((differentiable_Mw a).comp (by fun_prop)).pow n

theorem norm_Phi_le {a : ℝ} (ha : 0 ≤ a) (n : ℕ) (s : ℂ) :
    ‖Phi a n s‖ ≤ Real.exp (n * a * max 0 (1 / 2 - s.re)) := by
  unfold Phi
  rw [norm_pow]
  have h := norm_Mw_le ha (s - 1 / 2)
  have hre : -(s - 1 / 2).re = 1 / 2 - s.re := by simp
  rw [hre] at h
  calc ‖Mw a (s - 1 / 2)‖ ^ n ≤ Real.exp (a * max 0 (1 / 2 - s.re)) ^ n :=
        pow_le_pow_left₀ (norm_nonneg _) h n
    _ = Real.exp (n * a * max 0 (1 / 2 - s.re)) := by rw [← Real.exp_nat_mul]; ring_nf

theorem norm_Phi_le_inv {a : ℝ} (ha : 0 < a) (n : ℕ) {s : ℂ} (hs : s ≠ 1 / 2) :
    ‖Phi a n s‖ ≤ ((1 + Real.exp (a * (1 / 2 - s.re))) / (a * ‖s - 1 / 2‖)) ^ n := by
  unfold Phi
  rw [norm_pow]
  have h := norm_Mw_le_inv ha (sub_ne_zero.mpr hs)
  have hre : -(a * (s - 1 / 2).re) = a * (1 / 2 - s.re) := by simp; ring
  rw [hre] at h
  exact pow_le_pow_left₀ (norm_nonneg _) h n

theorem norm_one_sub_pow_le {z : ℂ} (hz : ‖z‖ ≤ 1) (n : ℕ) : ‖1 - z ^ n‖ ≤ n * ‖1 - z‖ := by
  induction n with
  | zero => simp
  | succ n ih =>
    have e : 1 - z ^ (n + 1) = (1 - z ^ n) + z ^ n * (1 - z) := by ring
    rw [e]
    calc ‖(1 - z ^ n) + z ^ n * (1 - z)‖ ≤ ‖1 - z ^ n‖ + ‖z ^ n * (1 - z)‖ := norm_add_le _ _
      _ ≤ n * ‖1 - z‖ + 1 * ‖1 - z‖ := by
          gcongr
          rw [norm_mul, norm_pow]
          exact mul_le_mul_of_nonneg_right (pow_le_one₀ (norm_nonneg _) hz) (norm_nonneg _)
      _ = ((n + 1 : ℕ) : ℝ) * ‖1 - z‖ := by push_cast; ring

theorem Phi_line {a : ℝ} (n : ℕ) (t : ℝ) : Phi a n (1 / 2 + t * I) = Mw a (t * I) ^ n := by
  unfold Phi; congr 2; ring

theorem norm_Phi_line_le {a : ℝ} (ha : 0 ≤ a) (n : ℕ) (t : ℝ) : ‖Phi a n (1 / 2 + t * I)‖ ≤ 1 := by
  rw [Phi_line, norm_pow]
  have h := norm_Mw_le ha (t * I)
  simp at h
  exact pow_le_one₀ (norm_nonneg _) h

theorem norm_one_sub_Phi_line_le {a : ℝ} (ha : 0 ≤ a) (n : ℕ) (t : ℝ) :
    ‖1 - Phi a n (1 / 2 + t * I)‖ ≤ n * (a * |t| / 2) := by
  rw [Phi_line]
  have h := norm_Mw_le ha (t * I)
  simp at h
  exact (norm_one_sub_pow_le h n).trans (mul_le_mul_of_nonneg_left (norm_one_sub_Mw_I_le ha t)
    (Nat.cast_nonneg n))

end Riemann.Deep.LeanNBRate

