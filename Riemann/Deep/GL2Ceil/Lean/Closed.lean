import Riemann.Deep.GL2Ceil.Lean.InnerMin

/-!
# GL2Ceil (gl2lean): closed forms of C_λ, B_λ (DERIVATION §1.4)

The certificate (ceil.py, Arb) evaluates C_λ = λ coth λ/2 + λ²/(2 sinh²λ) and
B_λ = (sinh 2λ/(4λ) - 1/2)/sinh²λ (λ = 0: C = 1, B = 1/3); here they are proved equal to the
integrals Cl λ = ∫₀¹ P_λ'², Bl λ = ∫₀¹ P_λ² used in the Lean statements.
-/

open Real MeasureTheory Set

noncomputable section

namespace GL2Ceil

lemma int_sinh_sq {lam : ℝ} (h : lam ≠ 0) :
    ∫ x in (0:ℝ)..1, Real.sinh (lam * x) ^ 2 = Real.sinh (2 * lam) / (4 * lam) - 1 / 2 := by
  have hd : ∀ x ∈ uIcc (0:ℝ) 1, HasDerivAt (fun x => Real.sinh (2 * lam * x) / (4 * lam) - x / 2)
      (Real.sinh (lam * x) ^ 2) x := by
    intro x _
    have h1 := ((((hasDerivAt_id' (x := x)).const_mul (2 * lam)).sinh).div_const (4 * lam)).fun_sub
      ((hasDerivAt_id' (x := x)).div_const 2)
    convert h1 using 1
    rw [show 2 * lam * x = 2 * (lam * x) by ring, Real.cosh_two_mul, Real.cosh_sq]
    field_simp
    ring
  rw [intervalIntegral.integral_eq_sub_of_hasDerivAt hd
    ((by fun_prop : Continuous fun x => Real.sinh (lam * x) ^ 2).intervalIntegrable 0 1)]
  simp

theorem Bl_closed {lam : ℝ} (h : lam ≠ 0) :
    Bl lam = (Real.sinh (2 * lam) / (4 * lam) - 1 / 2) / Real.sinh lam ^ 2 := by
  unfold Bl
  have e : ∀ x, Pl lam x ^ 2 = (1 / Real.sinh lam ^ 2) * Real.sinh (lam * x) ^ 2 := by
    intro x
    rw [show Pl lam x = Real.sinh (lam * x) / Real.sinh lam by simp [Pl, h]]
    ring
  simp_rw [e]
  rw [intervalIntegral.integral_const_mul, int_sinh_sq h]
  ring

theorem Cl_closed {lam : ℝ} (h : lam ≠ 0) :
    Cl lam = lam * (Real.cosh lam / Real.sinh lam) / 2 + lam ^ 2 / (2 * Real.sinh lam ^ 2) := by
  have e : Cl lam = Pld lam 1 - lam ^ 2 * Bl lam := by linarith [Cl_add lam]
  rw [e, Bl_closed h, show Pld lam 1 = lam * Real.cosh lam / Real.sinh lam by simp [Pld, h],
    Real.sinh_two_mul]
  have hs : Real.sinh lam ≠ 0 := Real.sinh_ne_zero.2 h
  field_simp
  ring

theorem Bl_zero : Bl 0 = 1 / 3 := by
  have e : ∀ x : ℝ, Pl 0 x ^ 2 = x ^ 2 := fun x => by simp [Pl]
  unfold Bl
  simp_rw [e]
  have hd : ∀ x ∈ uIcc (0:ℝ) 1, HasDerivAt (fun x : ℝ => x ^ 3 / 3) (x ^ 2) x := by
    intro x _
    have := ((hasDerivAt_id' (x := x)).fun_pow 3).div_const 3
    convert this using 1
    ring
  rw [intervalIntegral.integral_eq_sub_of_hasDerivAt hd
    ((by fun_prop : Continuous fun x : ℝ => x ^ 2).intervalIntegrable 0 1)]
  norm_num

end GL2Ceil


