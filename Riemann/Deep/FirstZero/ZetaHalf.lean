import Riemann.Deep.FirstZero.MinModel

/-!
# ζ(1/2) ∈ (-1.461, -1.460)

Order-8 Euler–Maclaurin at s = 1/2 with N = 2 (`em_err_le`): ζ(1/2) = 1 + 2^{-1/2}·EMKN 2 (1/2) + O(10⁻⁵),
EMKN 2 (1/2) = -72969743/20971520 exactly; √2 ∈ (1.41421, 1.41422). Im ζ(1/2) = 0 from `riemannZeta_conj`.
-/

open Complex

namespace Fayed.FirstZero.ZetaHalfCert
open Fayed.FirstZero Fayed.FirstZero.Min

lemma zeta_half_im : (riemannZeta (1 / 2 : ℂ)).im = 0 := by
  have h := riemannZeta_conj (1 / 2 : ℂ)
  have hc : (starRingEnd ℂ) (1 / 2 : ℂ) = 1 / 2 := by
    rw [show (1 / 2 : ℂ) = ((1 / 2 : ℝ) : ℂ) by push_cast; ring, Complex.conj_ofReal]
  rw [hc] at h
  exact Complex.conj_eq_iff_im.mp h.symm

lemma nps_half (n : ℕ) : nps (1 / 2 : ℂ) n = (((Real.sqrt n)⁻¹ : ℝ) : ℂ) := by
  rw [nps, Real.sqrt_eq_rpow, ← Real.rpow_neg (Nat.cast_nonneg n),
    Complex.ofReal_cpow (Nat.cast_nonneg n)]
  push_cast
  ring_nf

lemma K2_half : EMKN 2 (1 / 2 : ℂ) = ((-72969743 / 20971520 : ℝ) : ℂ) := by
  apply Complex.ext <;>
  simp [EMKN, gc, Complex.div_re, Complex.div_im, Complex.normSq] <;> norm_num

lemma Mf_half : Mf 2 (1 / 2 : ℂ) = ((1 + (Real.sqrt 2)⁻¹ * (-72969743 / 20971520) : ℝ) : ℂ) := by
  simp only [Mf, Nat.sub_self, Finset.range_zero, Finset.sum_empty, add_zero]
  rw [K2_half, nps_half]
  push_cast
  ring

lemma norm_half_add (x : ℝ) : ‖(1 / 2 : ℂ) + x‖ = |1 / 2 + x| := by
  rw [show (1 / 2 : ℂ) + x = ((1 / 2 + x : ℝ) : ℂ) by push_cast; ring, Complex.norm_real,
    Real.norm_eq_abs]

lemma em_half : ‖riemannZeta (1 / 2 : ℂ) - Mf 2 (1 / 2)‖ ≤ 1 / 100000 := by
  refine em_err_le (N := 2) (le_refl _) (σ0 := 1 / 2) (R := 1 / 2) (v := 1)
    (fun j => 3 / 2 + j) (by norm_num) ?_ (by norm_num) ?_ ?_ ?_ ?_
  · simp
  · have := norm_half_add 0
    simp only [Complex.ofReal_zero, add_zero] at this
    rw [this]; norm_num
  · intro j _
    have := norm_half_add (1 + j)
    push_cast at this
    rw [add_assoc, this, abs_of_nonneg (by positivity)]
    ring_nf; rfl
  · exact Real.rpow_le_one_of_one_le_of_nonpos (by norm_num) (by norm_num)
  · norm_num [Finset.prod_range_succ]

theorem zeta_half_re : -1.461 < (riemannZeta (1 / 2 : ℂ)).re ∧ (riemannZeta (1 / 2 : ℂ)).re < -1.460 := by
  have h := em_half
  rw [Mf_half] at h
  have hre := (Complex.abs_re_le_norm _).trans h
  rw [Complex.sub_re, Complex.ofReal_re, abs_le] at hre
  have s1 : (1.41421 : ℝ) < Real.sqrt 2 := by
    rw [Real.lt_sqrt (by norm_num)]; norm_num
  have s2 : Real.sqrt 2 < 1.41422 := by
    rw [Real.sqrt_lt' (by norm_num)]; norm_num
  have hs0 : 0 < Real.sqrt 2 := by linarith
  have i1 : (Real.sqrt 2)⁻¹ < 1 / 1.41421 := by
    rw [inv_eq_one_div]; exact one_div_lt_one_div_of_lt (by norm_num) s1
  have i2 : 1 / 1.41422 < (Real.sqrt 2)⁻¹ := by
    rw [inv_eq_one_div]; exact one_div_lt_one_div_of_lt hs0 s2
  constructor <;> nlinarith [hre.1, hre.2]

theorem zeta_half_neg : riemannZeta (1 / 2 : ℂ) ∈ {z : ℂ | z.im = 0 ∧ -1.461 < z.re ∧ z.re < -1.460} :=
  ⟨zeta_half_im, zeta_half_re⟩

end Fayed.FirstZero.ZetaHalfCert

namespace Fayed.FirstZero

/-- ζ(1/2) is real and lies in (-1.461, -1.460). -/
theorem zeta_half_neg : riemannZeta (1 / 2 : ℂ) ∈ {z : ℂ | z.im = 0 ∧ -1.461 < z.re ∧ z.re < -1.460} :=
  Fayed.FirstZero.ZetaHalfCert.zeta_half_neg

end Fayed.FirstZero

#print axioms Fayed.FirstZero.zeta_half_neg
