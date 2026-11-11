import Riemann.Deep.KimEff.Lean.Xi
import Mathlib.NumberTheory.Harmonic.ZetaAsymp
import Mathlib.Analysis.Complex.RealDeriv

/-!
# Explicit Kim–Lee wedge, part 9: γ(m) is the real Taylor coefficient

ξ is real on (0, ∞) (Mathlib's ζ(conj s) = conj ζ(s) and Γ real on the reals), so every
derivative of f(z) = 8ξ(1/2+z) at 0 is real: xiGammaC m = xiGamma m, and
Θ^{(n)}(w) = Σ_k γ(n+k) w^k/k! with the real γ of the Jensen polynomials.
-/

namespace KimEff

noncomputable section

open Complex

/-- ξ is real on the positive real axis. -/
lemma riemannXi_ofReal_im {x : ℝ} (hx : 0 < x) : (riemannXi x).im = 0 := by
  by_cases hx1 : x = 1
  · subst hx1
    simp [riemannXi]
  have hs0 : (x : ℂ) ≠ 0 := by exact_mod_cast hx.ne'
  have hs1 : (x : ℂ) ≠ 1 := by exact_mod_cast hx1
  have hre : 0 < (x : ℂ).re := by simpa using hx
  have hG : Gammaℝ (x : ℂ) ≠ 0 := Gammaℝ_ne_zero_of_re_pos hre
  have hΛ : completedRiemannZeta x = riemannZeta x * Gammaℝ x := by
    rw [riemannZeta_def_of_ne_zero hs0, div_mul_cancel₀ _ hG]
  have hz : (riemannZeta x).im = 0 := by
    have h := riemannZeta_conj (x : ℂ)
    rw [conj_ofReal] at h
    exact conj_eq_iff_im.mp h.symm
  have hg : (Gammaℝ x).im = 0 := by
    rw [Gammaℝ_def]
    have h1 : ((Real.pi : ℂ)) ^ (-(x : ℂ) / 2) = ((Real.pi ^ (-x / 2) : ℝ) : ℂ) := by
      rw [ofReal_cpow Real.pi_pos.le]
      push_cast
      ring_nf
    have h2 : Gamma ((x : ℂ) / 2) = ((Real.Gamma (x / 2) : ℝ) : ℂ) := by
      rw [← Gamma_ofReal]
      push_cast
      ring_nf
    rw [h1, h2, ← ofReal_mul, ofReal_im]
  rw [riemannXi_eq_completed hs0 hs1, hΛ]
  have ez : riemannZeta x = ((riemannZeta x).re : ℂ) := Complex.ext (by simp) (by simp [hz])
  have eg : Gammaℝ x = ((Gammaℝ x).re : ℂ) := Complex.ext (by simp) (by simp [hg])
  rw [ez, eg]
  have e : (x : ℂ) * (x - 1) * (((riemannZeta x).re : ℂ) * ((Gammaℝ x).re : ℂ)) / 2
      = ((x * (x - 1) * ((riemannZeta x).re * (Gammaℝ x).re) / 2 : ℝ) : ℂ) := by
    push_cast
    ring
  rw [e, ofReal_im]

/-- f(t) = 8ξ(1/2 + t) is real for real t > -1/2. -/
lemma xiF_ofReal_im {t : ℝ} (ht : -1 / 2 < t) : (xiF t).im = 0 := by
  unfold xiF
  have h : (1 / 2 + (t : ℂ)) = ((1 / 2 + t : ℝ) : ℂ) := by push_cast; ring
  rw [h]
  have h2 := riemannXi_ofReal_im (x := 1 / 2 + t) (by linarith)
  rw [Complex.mul_im, h2]
  simp

/-- A derivative of an entire function that is real on an open U ⊆ ℝ is real on U. -/
lemma deriv_im_eq_zero {g : ℂ → ℂ} (hg : Differentiable ℂ g) {U : Set ℝ} (hU : IsOpen U)
    (hreal : ∀ t ∈ U, (g t).im = 0) {t : ℝ} (ht : t ∈ U) : (deriv g t).im = 0 := by
  have h1 : HasDerivAt (fun z : ℂ => -I * g z) (-I * deriv g t) (t : ℂ) :=
    ((hg (t : ℂ)).hasDerivAt).const_mul (-I)
  have h2 : HasDerivAt (fun s : ℝ => (-I * g s).re) (-I * deriv g t).re t :=
    h1.real_of_complex
  have e1 : (fun s : ℝ => (-I * g s).re) = fun s : ℝ => (g s).im := by
    funext s
    simp [Complex.mul_re]
  have e2 : (-I * deriv g t).re = (deriv g t).im := by simp [Complex.mul_re]
  rw [e1, e2] at h2
  have h3 : HasDerivAt (fun s : ℝ => (g s).im) 0 t := by
    apply (hasDerivAt_const t (0 : ℝ)).congr_of_eventuallyEq
    filter_upwards [hU.mem_nhds ht] with s hs
    exact hreal s hs
  exact h2.unique h3

lemma differentiable_iteratedDeriv_of {g : ℂ → ℂ} (hg : Differentiable ℂ g) (k : ℕ) :
    Differentiable ℂ (iteratedDeriv k g) := by
  induction k with
  | zero => simpa using hg
  | succ k ih =>
    rw [iteratedDeriv_succ]
    exact ih.deriv

lemma iteratedDeriv_im_eq_zero {g : ℂ → ℂ} (hg : Differentiable ℂ g) {U : Set ℝ}
    (hU : IsOpen U) (hreal : ∀ t ∈ U, (g t).im = 0) (k : ℕ) :
    ∀ t ∈ U, (iteratedDeriv k g t).im = 0 := by
  induction k with
  | zero => simpa using hreal
  | succ k ih =>
    intro t ht
    rw [iteratedDeriv_succ]
    exact deriv_im_eq_zero (differentiable_iteratedDeriv_of hg k) hU ih ht

/-- The Taylor coefficients γ_C(m) are real. -/
theorem xiGammaC_im (m : ℕ) : (xiGammaC m).im = 0 := by
  have h := iteratedDeriv_im_eq_zero differentiable_xiF (isOpen_Ioi (a := (-1 / 2 : ℝ)))
    (fun t ht => xiF_ofReal_im ht) (2 * m) 0 (by norm_num)
  simp only [Complex.ofReal_zero] at h
  unfold xiGammaC
  have e : ((m.factorial : ℂ) / ((2 * m).factorial : ℂ))
      = (((m.factorial : ℝ) / ((2 * m).factorial : ℝ) : ℝ) : ℂ) := by push_cast; ring
  rw [e, im_ofReal_mul, h, mul_zero]

/-- γ_C(m) = γ(m): xiGamma is exactly the coefficient in 8ξ(1/2+z) = Σ γ(m) z^{2m}/m!. -/
theorem xiGammaC_eq (m : ℕ) : xiGammaC m = (xiGamma m : ℂ) :=
  Complex.ext (by simp [xiGamma]) (by simp [xiGamma, xiGammaC_im])

/-- 8ξ(1/2 + z) = Σ_m γ(m) z^{2m}/m! with the real γ of the Jensen polynomials. -/
theorem xiF_eq_tsum (z : ℂ) :
    xiF z = ∑' m : ℕ, (xiGamma m : ℂ) * z ^ (2 * m) / (m.factorial : ℂ) := by
  rw [← ThetaDeriv_zero_sq]
  unfold ThetaDeriv
  congr 1
  funext m
  rw [zero_add, xiGammaC_eq, pow_mul]

/-- Θ^{(n)}(w) = Σ_k γ(n+k) w^k/k! with the real γ. -/
theorem ThetaDeriv_eq (n : ℕ) (w : ℂ) :
    ThetaDeriv n w = ∑' k : ℕ, (xiGamma (n + k) : ℂ) * w ^ k / (k.factorial : ℂ) := by
  unfold ThetaDeriv
  congr 1
  funext k
  rw [xiGammaC_eq]

end

end KimEff

