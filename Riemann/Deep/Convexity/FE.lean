import Riemann.Deep.Convexity.Partial

/-!
# Functional-equation inputs

* `left_edge`: on Re s = 0, ‖ζ(it)‖ ≤ √(4t/25)·‖ζ(1+it)‖ (exact |Γ(1-it)|² = πt/sinh πt from the
  reflection formula; |χ(it)|² = (t/2π)·tanh(πt/2)).
* `apriori`: ‖ζ(σ+it)‖ ≤ 24·e^{3t} for 0 ≤ σ ≤ 2, t ≥ 1 (crude; from Partial and the
  functional equation with |Γ(w)| ≤ Γ(Re w)).
-/

open Complex Filter Topology Set MeasureTheory Asymptotics
open scoped Real ComplexConjugate

namespace Fayed.Convexity

noncomputable section

lemma norm_Gamma_le {w : ℂ} (hw : 0 < w.re) : ‖Complex.Gamma w‖ ≤ Real.Gamma w.re := by
  rw [Complex.Gamma_eq_integral hw, Complex.GammaIntegral, Real.Gamma_eq_integral hw]
  refine (norm_integral_le_integral_norm _).trans (le_of_eq ?_)
  refine setIntegral_congr_fun measurableSet_Ioi (fun x hx => ?_)
  rw [norm_mul, Complex.norm_real, Real.norm_of_nonneg (Real.exp_pos _).le,
    Complex.norm_cpow_eq_rpow_re_of_pos hx]
  simp

lemma Gamma_le_two {x : ℝ} (h1 : 1 / 2 ≤ x) (h2 : x ≤ 1) : Real.Gamma x ≤ 2 := by
  have key := Real.convexOn_Gamma.2 (show (1 / 2 : ℝ) ∈ Ioi 0 by norm_num)
    (show (1 : ℝ) ∈ Ioi 0 by norm_num) (show 0 ≤ 2 - 2 * x by linarith)
    (show 0 ≤ 2 * x - 1 by linarith) (show (2 - 2 * x) + (2 * x - 1) = 1 by ring)
  have hx : (2 - 2 * x) • (1 / 2 : ℝ) + (2 * x - 1) • (1 : ℝ) = x := by
    simp only [smul_eq_mul]; ring
  rw [hx, Real.Gamma_one_half_eq, Real.Gamma_one] at key
  simp only [smul_eq_mul] at key
  have hpi : √π ≤ 2 := by
    refine (Real.sqrt_le_sqrt (by linarith [Real.pi_le_four] : π ≤ 4)).trans (le_of_eq ?_)
    rw [show (4 : ℝ) = 2 ^ 2 by norm_num, Real.sqrt_sq (by norm_num)]
  nlinarith

lemma norm_cos_le (z : ℂ) : ‖Complex.cos z‖ ≤ Real.exp |z.im| := by
  have h1 : ‖Complex.exp (z * I)‖ ≤ Real.exp |z.im| := by
    rw [Complex.norm_exp]; apply Real.exp_le_exp.mpr; simp [neg_le_abs]
  have h2 : ‖Complex.exp (-z * I)‖ ≤ Real.exp |z.im| := by
    rw [Complex.norm_exp]; apply Real.exp_le_exp.mpr; simp [le_abs_self]
  rw [Complex.cos, norm_div, Complex.norm_two]
  have := norm_add_le (Complex.exp (z * I)) (Complex.exp (-z * I))
  linarith

lemma norm_zeta_conj (s : ℂ) : ‖riemannZeta (conj s)‖ = ‖riemannZeta s‖ := by
  rw [riemannZeta_conj, Complex.norm_conj]

lemma zeta_bound_right {s : ℂ} (h1 : 1 / 2 ≤ s.re) (h2 : s.re ≤ 2) (h3 : 1 ≤ s.im) :
    ‖riemannZeta s‖ ≤ 6 * Real.exp s.im := by
  have hU : s ∈ U := ⟨by show 0 < s.re; linarith, by show 0 < s.im; linarith⟩
  refine (norm_zeta_le_on_U hU).trans ?_
  have hs1 : 1 ≤ ‖s - 1‖ := by
    refine le_trans ?_ (abs_im_le_norm _)
    simp only [sub_im, one_im, sub_zero]
    rw [abs_of_nonneg (by linarith)]; exact h3
  have hsn : ‖s‖ ≤ 2 + s.im := by
    refine (norm_le_abs_re_add_abs_im s).trans ?_
    rw [abs_of_nonneg (by linarith), abs_of_nonneg (by linarith)]; linarith
  have hexp : 2 + s.im ≤ 2 * Real.exp s.im := by linarith [Real.add_one_le_exp s.im]
  have e1 : ‖s‖ / ‖s - 1‖ ≤ ‖s‖ := div_le_self (norm_nonneg _) hs1
  have e2 : 1 / s.re ≤ 2 := by rw [div_le_iff₀ (by linarith)]; linarith
  have e3 : ‖s‖ * (1 / s.re) ≤ ‖s‖ * 2 := mul_le_mul_of_nonneg_left e2 (norm_nonneg _)
  linarith

lemma apriori {s : ℂ} (h0 : 0 ≤ s.re) (h2 : s.re ≤ 2) (h3 : 1 ≤ s.im) :
    ‖riemannZeta s‖ ≤ 24 * Real.exp (3 * s.im) := by
  have hmono : Real.exp s.im ≤ Real.exp (3 * s.im) := Real.exp_le_exp.mpr (by linarith)
  rcases le_or_gt (1 / 2) s.re with h | h
  · have := zeta_bound_right h h2 h3
    linarith [Real.exp_pos s.im]
  set w : ℂ := 1 - s with hw_def
  have hw : ∀ n : ℕ, w ≠ -n := by
    intro n hn
    have := congrArg Complex.re hn
    simp [w] at this
    linarith [(n.cast_nonneg : (0 : ℝ) ≤ n)]
  have hw1 : w ≠ 1 := by
    intro hn; have := congrArg Complex.im hn; simp [w] at this; linarith
  have hfe := riemannZeta_one_sub hw hw1
  rw [show 1 - w = s by simp [w]] at hfe
  rw [hfe, norm_mul, norm_mul, norm_mul, norm_mul, Complex.norm_two]
  have hpow : ‖(2 * (π : ℂ)) ^ (-w)‖ ≤ 1 := by
    rw [show (2 * (π : ℂ)) = ((2 * π : ℝ) : ℂ) by push_cast; ring,
      Complex.norm_cpow_eq_rpow_re_of_pos (by positivity)]
    refine Real.rpow_le_one_of_one_le_of_nonpos (by linarith [Real.pi_gt_three]) ?_
    simp [w]; linarith
  have hG : ‖Complex.Gamma w‖ ≤ 2 := by
    have hre : w.re = 1 - s.re := by simp [w]
    refine (norm_Gamma_le (by rw [hre]; linarith)).trans ?_
    rw [hre]; exact Gamma_le_two (by linarith) (by linarith)
  have hcos : ‖Complex.cos (π * w / 2)‖ ≤ Real.exp (π * s.im / 2) := by
    refine (norm_cos_le _).trans (le_of_eq ?_)
    congr 1
    have him : (↑π * w / 2).im = -(π * s.im / 2) := by simp [w]; ring
    rw [him, abs_neg, abs_of_nonneg (by have := Real.pi_pos; positivity)]
  have hz : ‖riemannZeta w‖ ≤ 6 * Real.exp s.im := by
    have hc : (conj w).im = s.im := by simp [w]
    rw [← norm_zeta_conj w, ← hc]
    refine zeta_bound_right ?_ ?_ ?_ <;> simp [w] <;> linarith
  have hexp : Real.exp (π * s.im / 2) * Real.exp s.im ≤ Real.exp (3 * s.im) := by
    rw [← Real.exp_add]; apply Real.exp_le_exp.mpr
    nlinarith [Real.pi_lt_four]
  have hpos1 := norm_nonneg ((2 * (π : ℂ)) ^ (-w))
  have hpos2 := norm_nonneg (Complex.Gamma w)
  have hpos3 := norm_nonneg (Complex.cos (π * w / 2))
  have hpos4 := norm_nonneg (riemannZeta w)
  calc 2 * ‖(2 * (π : ℂ)) ^ (-w)‖ * ‖Complex.Gamma w‖ * ‖Complex.cos (π * w / 2)‖ *
        ‖riemannZeta w‖
      ≤ 2 * 1 * 2 * Real.exp (π * s.im / 2) * (6 * Real.exp s.im) := by gcongr
    _ = 24 * (Real.exp (π * s.im / 2) * Real.exp s.im) := by ring
    _ ≤ 24 * Real.exp (3 * s.im) := by linarith

lemma normSq_Gamma_one_sub {t : ℝ} (ht : 0 < t) :
    ‖Complex.Gamma (1 - t * I)‖ ^ 2 = π * t / Real.sinh (π * t) := by
  have hconj : Complex.Gamma (1 - t * I) = conj (Complex.Gamma (1 + t * I)) := by
    rw [← Complex.Gamma_conj]; congr 1; apply Complex.ext <;> simp
  have hadd : Complex.Gamma (1 + t * I) = (t * I) * Complex.Gamma (t * I) := by
    rw [show (1 : ℂ) + t * I = t * I + 1 by ring]
    exact Complex.Gamma_add_one _ (by simp [ht.ne'])
  have hrefl := Complex.Gamma_mul_Gamma_one_sub (t * I)
  set S : ℂ := ((Real.sinh (π * t) : ℝ) : ℂ) with hS_def
  have hsin : Complex.sin (π * (t * I)) = S * I := by
    rw [show (π : ℂ) * (t * I) = ((π * t : ℝ) : ℂ) * I by push_cast; ring, Complex.sin_mul_I,
      hS_def, Complex.ofReal_sinh]
  have hsh : Real.sinh (π * t) ≠ 0 := (Real.sinh_pos_iff.mpr (by positivity)).ne'
  have hS : S ≠ 0 := by rw [hS_def]; exact_mod_cast hsh
  have key : Complex.Gamma (1 + t * I) * conj (Complex.Gamma (1 + t * I)) =
      ((π * t / Real.sinh (π * t) : ℝ) : ℂ) := by
    rw [← hconj]
    nth_rewrite 1 [hadd]
    rw [mul_assoc, hrefl, hsin, Complex.ofReal_div, Complex.ofReal_mul, ← hS_def]
    field_simp
  rw [Complex.mul_conj, Complex.normSq_eq_norm_sq] at key
  rw [hconj, Complex.norm_conj]
  exact_mod_cast key

lemma left_edge {t : ℝ} (ht : 0 < t) :
    ‖riemannZeta (t * I)‖ ≤ √((4 / 25) * t) * ‖riemannZeta (1 + t * I)‖ := by
  set s : ℂ := 1 - t * I with hs_def
  have hs : ∀ n : ℕ, s ≠ -n := by
    intro n h; have := congrArg Complex.re h; simp [s] at this
    linarith [(n.cast_nonneg : (0 : ℝ) ≤ n)]
  have hs1 : s ≠ 1 := by
    intro h; have := congrArg Complex.im h; simp [s] at this; exact ht.ne' this
  have hfe := riemannZeta_one_sub hs hs1
  rw [show 1 - s = t * I by simp [s]] at hfe
  have hz : ‖riemannZeta s‖ = ‖riemannZeta (1 + t * I)‖ := by
    rw [show s = conj (1 + t * I) by apply Complex.ext <;> simp [s], norm_zeta_conj]
  have h2pi : ‖(2 * (π : ℂ)) ^ (-s)‖ = 1 / (2 * π) := by
    rw [show (2 * (π : ℂ)) = ((2 * π : ℝ) : ℂ) by push_cast; ring,
      Complex.norm_cpow_eq_rpow_re_of_pos (by positivity)]
    have : (-s).re = -1 := by simp [s]
    rw [this, Real.rpow_neg_one, one_div]
  have hsinh : 0 < Real.sinh (π * t / 2) := Real.sinh_pos_iff.mpr (by positivity)
  have hcos : ‖Complex.cos (π * s / 2)‖ = Real.sinh (π * t / 2) := by
    rw [show (π : ℂ) * s / 2 = π / 2 - ((π * t / 2 : ℝ) : ℂ) * I by simp [s]; ring,
      Complex.cos_pi_div_two_sub, Complex.sin_mul_I, ← Complex.ofReal_sinh, norm_mul,
      Complex.norm_I, mul_one, Complex.norm_real, Real.norm_of_nonneg hsinh.le]
  have hG := normSq_Gamma_one_sub ht
  rw [hfe, norm_mul, norm_mul, norm_mul, norm_mul, hz, h2pi, hcos, Complex.norm_two]
  apply mul_le_mul_of_nonneg_right _ (norm_nonneg _)
  have hpos : 0 ≤ 2 * (1 / (2 * π)) * ‖Complex.Gamma s‖ * Real.sinh (π * t / 2) := by
    have := Real.pi_pos; positivity
  rw [← abs_of_nonneg hpos]
  apply Real.abs_le_sqrt
  have hpi := Real.pi_pos
  have hcosh := Real.cosh_pos (π * t / 2)
  have hlt := Real.sinh_lt_cosh (π * t / 2)
  have hsh2 : Real.sinh (π * t) = 2 * Real.sinh (π * t / 2) * Real.cosh (π * t / 2) := by
    rw [← Real.sinh_two_mul]; ring_nf
  have e : (2 * (1 / (2 * π)) * ‖Complex.Gamma s‖ * Real.sinh (π * t / 2)) ^ 2 =
      t * Real.sinh (π * t / 2) / (2 * π * Real.cosh (π * t / 2)) := by
    rw [mul_pow, mul_pow, mul_pow, hG, hsh2]
    field_simp
  rw [e, div_le_iff₀ (by positivity)]
  have hpi2 : (25 : ℝ) / 8 ≤ π := by linarith [Real.pi_gt_d2]
  nlinarith [mul_pos ht hcosh, mul_pos ht hsinh]

end

end Fayed.Convexity
