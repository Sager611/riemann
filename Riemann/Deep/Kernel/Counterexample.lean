import Mathlib
import Riemann.Tree.R038

/-!
# R038 counterexample, fully machine-checked

K(t) = e^{-t²}(1 + t⁴/4) is positive, non-increasing on [0,∞), rapidly decreasing and strictly
log-concave, yet its two-sided Laplace transform
∫ K(t) e^{zt} dt = √π e^{z²/4} (19/16 + 3z²/16 + z⁴/64)
vanishes at some z with Re z ≠ 0 (z² = -6 + 2√10 i). This refutes "positive, non-increasing,
rapidly decreasing kernel ⇒ zeros of ∫ Ω(u) cosh(uz) du lie on z = iy" (Wünsche, arXiv 1607.04551).
-/

open MeasureTheory Complex Real

namespace Fayed.Kernel

/-- The kernel. -/
noncomputable def K (t : ℝ) : ℝ := Real.exp (-t ^ 2) * (1 + t ^ 4 / 4)

/-- (1) positivity. -/
theorem K_pos (t : ℝ) : 0 < K t := by
  unfold K; positivity

theorem continuous_K : Continuous K := by
  unfold K; fun_prop

theorem K_hasDerivAt (t : ℝ) :
    HasDerivAt K (-(t * Real.exp (-t ^ 2) * (2 - t ^ 2 + t ^ 4 / 2))) t := by
  have h := ((hasDerivAt_pow 2 t).neg).exp.mul (((hasDerivAt_pow 4 t).div_const 4).const_add 1)
  have hK : K = fun t => Real.exp (-t ^ 2) * (1 + t ^ 4 / 4) := rfl
  rw [hK]
  exact h.congr_deriv (by norm_num; ring)

/-- (2) K is non-increasing on [0, ∞). -/
theorem K_antitoneOn : AntitoneOn K (Set.Ici 0) := by
  refine antitoneOn_of_deriv_nonpos (convex_Ici 0) continuous_K.continuousOn
    (fun x _ => (K_hasDerivAt x).differentiableAt.differentiableWithinAt) ?_
  intro x hx
  rw [interior_Ici] at hx
  have hx : 0 < x := hx
  rw [(K_hasDerivAt x).deriv]
  have h2 : 0 < 2 - x ^ 2 + x ^ 4 / 2 := by nlinarith [sq_nonneg (x ^ 2 - 1)]
  have := Real.exp_pos (-x ^ 2)
  apply neg_nonpos.mpr
  positivity

theorem log_K (t : ℝ) : Real.log (K t) = -t ^ 2 + Real.log (1 + t ^ 4 / 4) := by
  unfold K
  rw [Real.log_mul (Real.exp_pos _).ne' (by positivity), Real.log_exp]

theorem logK_hasDerivAt (t : ℝ) :
    HasDerivAt (fun t => Real.log (K t)) (-2 * t + t ^ 3 / (1 + t ^ 4 / 4)) t := by
  have hu : (1 + t ^ 4 / 4) ≠ 0 := by positivity
  have h := ((hasDerivAt_pow 2 t).neg).add
    ((((hasDerivAt_pow 4 t).div_const 4).const_add 1).log hu)
  have hL : (fun t => Real.log (K t)) = fun t => -t ^ 2 + Real.log (1 + t ^ 4 / 4) :=
    funext log_K
  rw [hL]
  exact h.congr_deriv (by field_simp; ring)

theorem logK_deriv : deriv (fun t => Real.log (K t)) = fun t => -2 * t + t ^ 3 / (1 + t ^ 4 / 4) :=
  funext fun t => (logK_hasDerivAt t).deriv

theorem logK_deriv2_hasDerivAt (t : ℝ) :
    HasDerivAt (fun t : ℝ => -2 * t + t ^ 3 / (1 + t ^ 4 / 4))
      (-2 + (3 * t ^ 2 * (1 + t ^ 4 / 4) - t ^ 3 * t ^ 3) / (1 + t ^ 4 / 4) ^ 2) t := by
  have hu : (1 + t ^ 4 / 4) ≠ 0 := by positivity
  have h := ((hasDerivAt_id t).const_mul (-2)).add
    ((hasDerivAt_pow 3 t).div (((hasDerivAt_pow 4 t).div_const 4).const_add 1) hu)
  exact h.congr_deriv (by norm_num)

/-- (3) strict log-concavity of K on ℝ, via the polynomial certificate of R038. -/
theorem K_strictLogConcave : StrictConcaveOn ℝ Set.univ (fun t => Real.log (K t)) := by
  refine strictConcaveOn_of_deriv2_neg' convex_univ
    (continuous_K.log fun t => (K_pos t).ne').continuousOn ?_
  intro t _
  rw [Function.iterate_succ_apply, Function.iterate_one, logK_deriv,
    (logK_deriv2_hasDerivAt t).deriv]
  have hu : 0 < 1 + t ^ 4 / 4 := by positivity
  have key := Fayed.R38.logconcave_poly (t ^ 2) (sq_nonneg t)
  have : (3 * t ^ 2 * (1 + t ^ 4 / 4) - t ^ 3 * t ^ 3) / (1 + t ^ 4 / 4) ^ 2 < 2 := by
    rw [div_lt_iff₀ (by positivity)]
    nlinarith [key]
  linarith

/-! ## (4) Gaussian moments -/

/-- Gaussian with linear term. -/
noncomputable def g (z : ℂ) (t : ℝ) : ℂ := cexp (-(t : ℂ) ^ 2 + z * t)

/-- k-th moment of g. -/
noncomputable def M (z : ℂ) (k : ℕ) : ℂ := ∫ t : ℝ, (t : ℂ) ^ k * g z t

theorem integrable_mom (z : ℂ) (k : ℕ) : Integrable (fun t : ℝ => (t : ℂ) ^ k * g z t) := by
  have hk : (-1 : ℝ) < k := by have := k.cast_nonneg (α := ℝ); linarith
  have h1 : Integrable (fun t : ℝ => ((t ^ (k : ℝ) * Real.exp (-(1 / 2) * t ^ 2) : ℝ) : ℂ)) :=
    (integrable_rpow_mul_exp_neg_mul_sq (by norm_num) hk).ofReal
  have hb : ∀ t : ℝ, ‖cexp (((-t ^ 2 / 2 : ℝ) : ℂ) + z * t)‖ ≤ Real.exp (z.re ^ 2 / 2) := by
    intro t
    rw [Complex.norm_exp]
    apply Real.exp_le_exp.mpr
    simp only [add_re, ofReal_re, mul_re, ofReal_im, mul_zero, sub_zero]
    nlinarith [sq_nonneg (t - z.re)]
  have h2 := h1.mul_bdd (g := fun t : ℝ => cexp (((-t ^ 2 / 2 : ℝ) : ℂ) + z * t))
    (by fun_prop) (ae_of_all _ hb)
  refine h2.congr (ae_of_all _ fun t => ?_)
  simp only [g, Real.rpow_natCast]
  push_cast
  rw [mul_assoc, ← Complex.exp_add]
  congr 1
  congr 1
  ring

theorem M_zero (z : ℂ) : M z 0 = (Real.sqrt π : ℂ) * cexp (z ^ 2 / 4) := by
  have h := integral_cexp_quadratic (b := -1) (by simp) z 0
  have e : (fun t : ℝ => (t : ℂ) ^ 0 * g z t) = fun t : ℝ => cexp (-1 * (t : ℂ) ^ 2 + z * t + 0) := by
    funext t; simp only [g, pow_zero, one_mul]; ring_nf
  have hp : ((π : ℂ) / -(-1)) ^ (1 / 2 : ℂ) = (Real.sqrt π : ℂ) := by
    rw [neg_neg, div_one, Real.sqrt_eq_rpow, Complex.ofReal_cpow Real.pi_pos.le]
    norm_num
  rw [M, e, h, hp]
  congr 1
  congr 1
  ring

theorem g_hasDerivAt (z : ℂ) (t : ℝ) : HasDerivAt (g z) ((-2 * (t : ℂ) + z) * g z t) t := by
  have h1 : HasDerivAt (fun t : ℝ => (t : ℂ)) 1 t := by simpa using (hasDerivAt_id t).ofReal_comp
  have h2 := ((h1.pow 2).neg.add (h1.const_mul z)).cexp
  show HasDerivAt (fun t : ℝ => cexp (-(t : ℂ) ^ 2 + z * t))
    ((-2 * (t : ℂ) + z) * cexp (-(t : ℂ) ^ 2 + z * t)) t
  exact h2.congr_deriv (by norm_num; ring)

theorem M_one (z : ℂ) : 0 = z * M z 0 - 2 * M z 1 := by
  have h := integral_eq_zero_of_hasDerivAt_of_integrable (g_hasDerivAt z)
    ((((integrable_mom z 0).const_mul z).sub ((integrable_mom z 1).const_mul 2)).congr
      (ae_of_all _ fun t => by simp only [Pi.sub_apply, pow_zero, pow_one, one_mul]; ring))
    (by simpa using integrable_mom z 0)
  rw [← h, M, M, ← integral_const_mul, ← integral_const_mul,
    ← integral_sub ((integrable_mom z 0).const_mul z) ((integrable_mom z 1).const_mul 2)]
  congr 1
  funext t
  simp only [g, pow_zero, pow_one, one_mul]
  ring

theorem M_rec (z : ℂ) (k : ℕ) :
    0 = (k + 1 : ℂ) * M z k + z * M z (k + 1) - 2 * M z (k + 2) := by
  have hd : ∀ t : ℝ, HasDerivAt (fun t : ℝ => (t : ℂ) ^ (k + 1) * g z t)
      ((k + 1 : ℂ) * ((t : ℂ) ^ k * g z t) + z * ((t : ℂ) ^ (k + 1) * g z t)
        - 2 * ((t : ℂ) ^ (k + 2) * g z t)) t := by
    intro t
    have h1 : HasDerivAt (fun t : ℝ => (t : ℂ)) 1 t := by
      simpa using (hasDerivAt_id t).ofReal_comp
    convert (h1.pow (k + 1)).mul (g_hasDerivAt z t) using 1
    simp only [Pi.pow_apply, Nat.add_sub_cancel]
    push_cast
    ring
  have hint := (((integrable_mom z k).const_mul (k + 1 : ℂ)).add
    ((integrable_mom z (k + 1)).const_mul z)).sub ((integrable_mom z (k + 2)).const_mul 2)
  have h := integral_eq_zero_of_hasDerivAt_of_integrable hd hint (integrable_mom z (k + 1))
  rw [← h, integral_sub
      (f := fun x : ℝ => (k + 1 : ℂ) * ((x : ℂ) ^ k * g z x) + z * ((x : ℂ) ^ (k + 1) * g z x))
      (g := fun x : ℝ => 2 * ((x : ℂ) ^ (k + 2) * g z x)) (((integrable_mom z k).const_mul _).add
      ((integrable_mom z (k + 1)).const_mul z)) ((integrable_mom z (k + 2)).const_mul 2),
    integral_add ((integrable_mom z k).const_mul _) ((integrable_mom z (k + 1)).const_mul z),
    integral_const_mul, integral_const_mul, integral_const_mul, M, M, M]

/-- (4) Gaussian-moment identity for the two-sided Laplace transform of K. -/
theorem K_transform (z : ℂ) :
    ∫ t : ℝ, (K t : ℂ) * cexp (z * t) =
      (Real.sqrt π : ℂ) * cexp (z ^ 2 / 4) * (19 / 16 + 3 * z ^ 2 / 16 + z ^ 4 / 64) := by
  have e : (fun t : ℝ => (K t : ℂ) * cexp (z * t)) =
      fun t : ℝ => (t : ℂ) ^ 0 * g z t + (1 / 4 : ℂ) * ((t : ℂ) ^ 4 * g z t) := by
    funext t
    simp only [K, g]
    push_cast
    rw [Complex.exp_add]
    ring
  rw [e, integral_add (integrable_mom z 0) ((integrable_mom z 4).const_mul _),
    integral_const_mul]
  change M z 0 + 1 / 4 * M z 4 = _
  have e0 := M_zero z
  have e1 := M_one z
  have e2 : 0 = ((0 : ℕ) + 1 : ℂ) * M z 0 + z * M z 1 - 2 * M z 2 := M_rec z 0
  have e3 : 0 = ((1 : ℕ) + 1 : ℂ) * M z 1 + z * M z 2 - 2 * M z 3 := M_rec z 1
  have e4 : 0 = ((2 : ℕ) + 1 : ℂ) * M z 2 + z * M z 3 - 2 * M z 4 := M_rec z 2
  push_cast at e2 e3 e4
  set A := (Real.sqrt π : ℂ) * cexp (z ^ 2 / 4)
  rw [e0] at e1 e2 ⊢
  have h1 : M z 1 = z * A / 2 := by linear_combination (1 / 2 : ℂ) * e1
  have h2 : M z 2 = A / 2 + z ^ 2 * A / 4 := by
    linear_combination (1 / 2 : ℂ) * e2 + (z / 2) * h1
  have h3 : M z 3 = 3 * z * A / 4 + z ^ 3 * A / 8 := by
    linear_combination (1 / 2 : ℂ) * e3 + h1 + (z / 2) * h2
  have h4 : M z 4 = 3 * A / 4 + 3 * z ^ 2 * A / 4 + z ^ 4 * A / 16 := by
    linear_combination (1 / 2 : ℂ) * e4 + (3 / 2 : ℂ) * h2 + (z / 2) * h3
  rw [h4]
  ring

/-- (5) The transform of K has a zero off the imaginary axis. -/
theorem exists_zero_off_axis :
    ∃ z : ℂ, z.re ≠ 0 ∧ ∫ t : ℝ, (K t : ℂ) * cexp (z * t) = 0 := by
  obtain ⟨z, hz⟩ := IsAlgClosed.exists_pow_nat_eq (-6 + 2 * (Real.sqrt 10 : ℂ) * Complex.I) (n := 2) (by norm_num)
  have h10 : 0 < Real.sqrt 10 := by positivity
  refine ⟨z, ?_, ?_⟩
  · intro h0
    have := congrArg Complex.im hz
    simp [sq, h0] at this
  · have hs : (Real.sqrt 10 : ℂ) ^ 2 = 10 := by
      rw [← ofReal_pow, Real.sq_sqrt (by norm_num)]; norm_num
    rw [K_transform, show z ^ 4 = (z ^ 2) ^ 2 by ring, hz]
    apply mul_eq_zero_of_right
    linear_combination (Complex.I ^ 2 / 16) * hs + (5 / 8 : ℂ) * Complex.I_sq

end Fayed.Kernel

#print axioms Fayed.Kernel.K_pos
#print axioms Fayed.Kernel.K_antitoneOn
#print axioms Fayed.Kernel.K_strictLogConcave
#print axioms Fayed.Kernel.K_transform
#print axioms Fayed.Kernel.exists_zero_off_axis
