import Riemann.Deep.GL2Ceil.Lean.Crude
import Mathlib.Analysis.SpecialFunctions.Trigonometric.ArctanDeriv

/-!
# GL2Ceil (gl2lean): G(R) = e^R gd(R)/(4R), gd(R) = 2 arctan tanh(R/2) (the checker's crude formula)

ceil.py crude(): K = (C/ξ) e^R gd(R)/(4R) with gd = 2·atan(tanh(R/2)); here G(R) (Crude.lean) is proved
equal to e^R gd(R)/(4R), gd(R) = ∫₀^R sech = 2 arctan tanh(R/2), so CrudeCellCert and the crude tails
are the checker's inequalities.
-/

open Real MeasureTheory Set

noncomputable section

namespace GL2Ceil

/-- gd(R) = ∫₀^R sech t dt. -/
def gd (R : ℝ) : ℝ := ∫ t in (0:ℝ)..R, 1 / Real.cosh t

lemma G_eq {R : ℝ} (hR : R ≠ 0) : G R = Real.exp R * gd R / (4 * R) := by
  have hpt : ∀ y : ℝ, 1 / (Real.exp (2 * R * (y - 1)) + Real.exp (-(2 * R * y)))
      = Real.exp R / (2 * Real.cosh (R - 2 * R * y)) := by
    intro y
    have h1 : Real.exp R * Real.exp (2 * R * (y - 1)) = Real.exp (-(R - 2 * R * y)) := by
      rw [← Real.exp_add, show R + 2 * R * (y - 1) = -(R - 2 * R * y) by ring]
    have h2 : Real.exp R * Real.exp (-(2 * R * y)) = Real.exp (R - 2 * R * y) := by
      rw [← Real.exp_add, show R + -(2 * R * y) = R - 2 * R * y by ring]
    rw [Real.cosh_eq, ← h1, ← h2]
    have := Real.exp_pos R
    have := Real.exp_pos (2 * R * (y - 1))
    have := Real.exp_pos (-(2 * R * y))
    field_simp
    ring
  unfold G gd
  simp_rw [hpt]
  have h := intervalIntegral.integral_comp_sub_mul (fun t => Real.exp R / (2 * Real.cosh t))
    (a := 0) (b := 1 / 2) (mul_ne_zero two_ne_zero hR) R
  rw [h, show R - 2 * R * (1 / 2) = 0 by ring, show R - 2 * R * 0 = R by ring, smul_eq_mul]
  have h3 : ∫ x in (0:ℝ)..R, Real.exp R / (2 * Real.cosh x)
      = Real.exp R / 2 * ∫ x in (0:ℝ)..R, 1 / Real.cosh x := by
    rw [← intervalIntegral.integral_const_mul]
    congr 1
    funext x
    rw [div_mul_div_comm, mul_one]
  rw [h3]
  field_simp
  ring

lemma alg_gd (s c : ℝ) (hc : 0 < c) (h : c ^ 2 = s ^ 2 + 1) :
    1 / (c ^ 2 + s ^ 2) = 2 * (1 / (1 + (s / c) ^ 2) *
      ((c * (1 / 2) * c - s * (s * (1 / 2))) / c ^ 2)) := by
  have hc2 : c ^ 2 ≠ 0 := by positivity
  have e3 : 1 + (s / c) ^ 2 = (c ^ 2 + s ^ 2) / c ^ 2 := by
    rw [div_pow, eq_div_iff hc2, add_mul, div_mul_cancel₀ _ hc2]
    ring
  have e4 : (c * (1 / 2) * c - s * (s * (1 / 2))) / c ^ 2 = 1 / (2 * c ^ 2) := by
    rw [div_eq_div_iff hc2 (by positivity)]
    linear_combination c ^ 2 * h
  rw [e3, e4]
  field_simp

lemma gd_hasDerivAt (x : ℝ) :
    HasDerivAt (fun x => 2 * Real.arctan (Real.sinh (x / 2) / Real.cosh (x / 2)))
      (1 / Real.cosh x) x := by
  have hh : HasDerivAt (fun x : ℝ => x / 2) (1 / 2) x := (hasDerivAt_id' (x := x)).div_const 2
  have hcpos := Real.cosh_pos (x / 2)
  have ha := ((hh.sinh.fun_div hh.cosh hcpos.ne').arctan).const_mul 2
  convert ha using 1
  have e1 : Real.cosh x = Real.cosh (x / 2) ^ 2 + Real.sinh (x / 2) ^ 2 := by
    rw [← Real.cosh_two_mul]; ring_nf
  rw [e1]
  exact alg_gd _ _ hcpos (Real.cosh_sq (x / 2))

lemma gd_eq (R : ℝ) : gd R = 2 * Real.arctan (Real.tanh (R / 2)) := by
  unfold gd
  rw [intervalIntegral.integral_eq_sub_of_hasDerivAt (fun x _ => gd_hasDerivAt x)
    ((continuous_const.div₀ Real.continuous_cosh (fun t => (Real.cosh_pos t).ne')).intervalIntegrable _ _)]
  simp [Real.tanh_eq_sinh_div_cosh]

/-- G in the checker's closed form: G(R) = e^R · 2 arctan tanh(R/2) / (4R). -/
theorem G_closed {R : ℝ} (hR : R ≠ 0) :
    G R = Real.exp R * (2 * Real.arctan (Real.tanh (R / 2))) / (4 * R) := by
  rw [G_eq hR, gd_eq]

end GL2Ceil


