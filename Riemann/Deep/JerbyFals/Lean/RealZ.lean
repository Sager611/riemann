import Riemann.Deep.JerbyFals.Lean.Phase
import Mathlib.NumberTheory.Harmonic.ZetaAsymp

/-!
# Z is real: Z(t) = e^{iθ(t)} ζ(1/2 + it) exactly

From `gamma_phase`: Gammaℝ(1/2 + it) = ρ e^{iθ(t)} and Gammaℝ(1/2 − it) = ρ e^{−iθ(t)} with
ρ = π^{−1/4} |Γ(1/4 + it/2)| > 0. With ζ = Λ/Gammaℝ, Λ(1 − s) = Λ(s) and ζ(conj s) = conj ζ(s)
(Mathlib `riemannZeta_conj`), conj (e^{iθ} ζ(1/2 + it)) = e^{iθ} ζ(1/2 + it). So the real part in
`hardyZ` loses nothing (2406.06548 l. 182–184: "Z(t) = e^{iθ(t)} ζ(1/2+it) ... the real function").
-/

open Complex
open scoped Real ComplexConjugate

namespace JerbyLean

noncomputable section

/-- ρ(t) = π^{−1/4} |Γ(1/4 + it/2)|. -/
def rhoZ (t : ℝ) : ℝ := Real.exp (-(Real.log π) / 4) * ‖Gamma (1 / 4 + t / 2 * I)‖

lemma rhoZ_pos (t : ℝ) : 0 < rhoZ t :=
  mul_pos (Real.exp_pos _) (norm_pos_iff.2 (Gamma_ne_zero_of_re_pos (quarter_line_re_pos t)))

lemma pi_cpow (z : ℂ) : (π : ℂ) ^ z = exp ((Real.log π : ℂ) * z) := by
  rw [cpow_def_of_ne_zero (ofReal_ne_zero.2 Real.pi_ne_zero), ofReal_log Real.pi_pos.le]

lemma exp_mul_mul_exp (a b c : ℂ) : exp a * (b * exp c) = b * exp (a + c) := by
  rw [Complex.exp_add]; ring

lemma exp_mul_mul_exp' (a b c : ℂ) : exp a * b * exp c = b * exp (a + c) := by
  rw [Complex.exp_add]; ring

lemma Gammaℝ_line (t : ℝ) :
    Gammaℝ (1 / 2 + t * I) = (rhoZ t : ℂ) * exp (I * theta t) := by
  have h2 : (1 / 2 + t * I : ℂ) / 2 = 1 / 4 + t / 2 * I := by ring
  rw [Gammaℝ_def, h2, gamma_phase t, pi_cpow, rhoZ]
  push_cast
  rw [exp_mul_mul_exp, exp_mul_mul_exp']
  congr 2
  ring

lemma Gammaℝ_line_conj (t : ℝ) :
    Gammaℝ (1 - (1 / 2 + t * I)) = (rhoZ t : ℂ) * exp (-(I * theta t)) := by
  have h2 : (1 - (1 / 2 + t * I) : ℂ) / 2 = conj (1 / 4 + t / 2 * I) := by
    simp only [map_add, map_mul, map_div₀, map_one, map_ofNat, conj_ofReal, conj_I]; ring
  rw [Gammaℝ_def, h2, Gamma_conj, gamma_phase t, pi_cpow, rhoZ, map_mul, conj_ofReal,
    ← exp_conj, map_mul, conj_I, conj_ofReal]
  push_cast
  rw [exp_mul_mul_exp, exp_mul_mul_exp']
  congr 2
  ring

/-- **Z is real**: Z(t) = e^{iθ(t)} ζ(1/2 + it) as complex numbers. -/
theorem hardyZ_eq (t : ℝ) :
    (hardyZ t : ℂ) = exp (I * theta t) * riemannZeta (1 / 2 + t * I) := by
  have hs0 : (1 / 2 + t * I : ℂ) ≠ 0 := by
    intro h; have := congrArg re h; norm_num at this
  have hs1 : (1 - (1 / 2 + t * I) : ℂ) ≠ 0 := by
    intro h; have := congrArg re h; norm_num at this
  have hconj : conj (1 / 2 + t * I : ℂ) = 1 - (1 / 2 + t * I) := by
    simp only [map_add, map_mul, map_div₀, map_one, map_ofNat, conj_ofReal, conj_I]; ring
  have hρ : (rhoZ t : ℂ) ≠ 0 := ofReal_ne_zero.2 (rhoZ_pos t).ne'
  have key : ∀ e L : ℂ, e ≠ 0 → e * (L / ((rhoZ t : ℂ) * e)) = L / (rhoZ t : ℂ) := by
    intro e L he; field_simp
  have hw : conj (exp (I * theta t) * riemannZeta (1 / 2 + t * I)) =
      exp (I * theta t) * riemannZeta (1 / 2 + t * I) := by
    rw [map_mul, ← exp_conj, ← riemannZeta_conj, hconj, riemannZeta_def_of_ne_zero hs1,
      completedRiemannZeta_one_sub, riemannZeta_def_of_ne_zero hs0, Gammaℝ_line,
      Gammaℝ_line_conj]
    have hc : conj (I * (theta t : ℂ)) = -(I * theta t) := by simp
    rw [hc, key _ _ (exp_ne_zero _), key _ _ (exp_ne_zero _)]
  exact conj_eq_iff_re.1 hw

end

end JerbyLean
