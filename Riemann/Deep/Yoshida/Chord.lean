/-
Riemann/Deep/Yoshida/Chord.lean — second-order (chord) lower bounds for the two parts of the weight.

* `mu_chord`: `2πμ(r)` is concave in `u = (r/2)²` (every digamma term `−c/(c² + u)` is), so on a cell
  it lies above the chord through its values at the ends, measured in `u`.
* `J_chord`: `J(r) = ∫_{−2a}^{2a} 2cosh(x/2) cos(rx) dx` and `y ↦ cos(yx) − x²y²/2` is concave, so
  `J(r) ≥ chord − (16/3) a³ cosh a · θ(1−θ)(hi − lo)²`.
Used by the affine-cell certificates (`Affine.lean`), whose loss per cell is second order.
-/
import Riemann.Deep.Yoshida.Hybrid
import Riemann.Deep.Yoshida.MuTail

noncomputable section

open Finset

namespace Fayed.Yoshida

open Zeta23

/-- Convexity of `1/X` at two points. -/
lemma inv_chord {X0 X1 θ : ℝ} (h0 : 0 < X0) (h1 : 0 < X1) (hθ0 : 0 ≤ θ) (hθ1 : θ ≤ 1) :
    1 / ((1 - θ) * X0 + θ * X1) ≤ (1 - θ) / X0 + θ / X1 := by
  have hX : 0 < (1 - θ) * X0 + θ * X1 := by
    rcases le_total X0 X1 with h | h <;> nlinarith
  rw [div_add_div _ _ h0.ne' h1.ne', div_le_div_iff₀ hX (by positivity)]
  nlinarith [mul_nonneg (mul_nonneg hθ0 (sub_nonneg.2 hθ1)) (sq_nonneg (X0 - X1))]

/-- One digamma term: `(1−θ) c/(c²+u₀) + θ c/(c²+u₁) ≥ c/(c² + u)`, `u = (1−θ)u₀ + θu₁`. -/
lemma term_chord {c u0 u1 θ : ℝ} (hc : 0 < c) (hu0 : 0 ≤ u0) (hu1 : 0 ≤ u1) (hθ0 : 0 ≤ θ)
    (hθ1 : θ ≤ 1) :
    c / (c ^ 2 + ((1 - θ) * u0 + θ * u1)) ≤ (1 - θ) * (c / (c ^ 2 + u0)) + θ * (c / (c ^ 2 + u1)) := by
  have h := inv_chord (X0 := c ^ 2 + u0) (X1 := c ^ 2 + u1) (by positivity) (by positivity) hθ0 hθ1
  have e : (1 - θ) * (c ^ 2 + u0) + θ * (c ^ 2 + u1) = c ^ 2 + ((1 - θ) * u0 + θ * u1) := by ring
  rw [e] at h
  have h2 := mul_le_mul_of_nonneg_left h hc.le
  have e1 : c * (1 / (c ^ 2 + ((1 - θ) * u0 + θ * u1))) = c / (c ^ 2 + ((1 - θ) * u0 + θ * u1)) := by
    ring
  have e2 : c * ((1 - θ) / (c ^ 2 + u0) + θ / (c ^ 2 + u1))
      = (1 - θ) * (c / (c ^ 2 + u0)) + θ * (c / (c ^ 2 + u1)) := by ring
  rw [e1, e2] at h2
  exact h2

/-- The exact series for `2πμ`. -/
lemma two_pi_mu_series (x : ℝ) : 2 * Real.pi * mu x = -Real.eulerMascheroniConstant
    - (1 / 4) / ((1 / 4) ^ 2 + (x / 2) ^ 2)
    + ∑' n : ℕ, (1 / ((n : ℝ) + 1) - ((n : ℝ) + 1 + 1 / 4) / (((n : ℝ) + 1 + 1 / 4) ^ 2 + (x / 2) ^ 2))
    - Real.log Real.pi := by
  have hπ := Real.pi_pos
  rw [MuFields.mu_eq, MuFields.re_digamma_vertical (a := 1 / 4) (by norm_num) (by norm_num)]
  field_simp

/-- **Chord bound for `μ`** (concavity in `(x/2)²`). -/
theorem mu_chord {x0 x1 x θ : ℝ} (hθ0 : 0 ≤ θ) (hθ1 : θ ≤ 1)
    (hx : (x / 2) ^ 2 = (1 - θ) * (x0 / 2) ^ 2 + θ * (x1 / 2) ^ 2) :
    (1 - θ) * (2 * Real.pi * mu x0) + θ * (2 * Real.pi * mu x1) ≤ 2 * Real.pi * mu x := by
  rw [two_pi_mu_series, two_pi_mu_series, two_pi_mu_series, hx]
  set u0 := (x0 / 2) ^ 2
  set u1 := (x1 / 2) ^ 2
  have hu0 : 0 ≤ u0 := by positivity
  have hu1 : 0 ≤ u1 := by positivity
  have hs0 := MuFields.summable_re_terms (a := 1 / 4) (by norm_num) (by norm_num) (x0 / 2)
  have hs1 := MuFields.summable_re_terms (a := 1 / 4) (by norm_num) (by norm_num) (x1 / 2)
  have hs := MuFields.summable_re_terms (a := 1 / 4) (by norm_num) (by norm_num) (x / 2)
  rw [hx] at hs
  have hterm : ∀ n : ℕ, (1 - θ) * (1 / ((n : ℝ) + 1) - ((n : ℝ) + 1 + 1 / 4) / (((n : ℝ) + 1 + 1 / 4) ^ 2 + u0))
      + θ * (1 / ((n : ℝ) + 1) - ((n : ℝ) + 1 + 1 / 4) / (((n : ℝ) + 1 + 1 / 4) ^ 2 + u1))
      ≤ 1 / ((n : ℝ) + 1) - ((n : ℝ) + 1 + 1 / 4)
          / (((n : ℝ) + 1 + 1 / 4) ^ 2 + ((1 - θ) * u0 + θ * u1)) := by
    intro n
    have := term_chord (c := (n : ℝ) + 1 + 1 / 4) (by positivity) hu0 hu1 hθ0 hθ1
    nlinarith
  have hsum := ((hs0.mul_left (1 - θ)).add (hs1.mul_left θ)).tsum_le_tsum hterm hs
  rw [Summable.tsum_add (hs0.mul_left (1 - θ)) (hs1.mul_left θ), tsum_mul_left, tsum_mul_left] at hsum
  have h0 := term_chord (c := 1 / 4) (by norm_num) hu0 hu1 hθ0 hθ1
  nlinarith

/-- `y ↦ cos(yx) − x²y²/2` is concave. -/
lemma concave_cos_sub (x : ℝ) : ConcaveOn ℝ Set.univ (fun y : ℝ => Real.cos (y * x) - x ^ 2 * y ^ 2 / 2) := by
  have hd : ∀ y, HasDerivAt (fun y : ℝ => Real.cos (y * x) - x ^ 2 * y ^ 2 / 2)
      (-Real.sin (y * x) * x - x ^ 2 * y) y := by
    intro y
    have h1 : HasDerivAt (fun y : ℝ => Real.cos (y * x)) (-Real.sin (y * x) * x) y :=
      (hasDerivAt_mul_const x).cos
    have h2 : HasDerivAt (fun y : ℝ => x ^ 2 * y ^ 2 / 2) (x ^ 2 * y) y := by
      have := ((hasDerivAt_pow 2 y).const_mul (x ^ 2)).div_const 2
      convert this using 1
      push_cast
      ring
    exact h1.sub h2
  have hd2 : ∀ y, HasDerivAt (fun y : ℝ => -Real.sin (y * x) * x - x ^ 2 * y)
      (-Real.cos (y * x) * x * x - x ^ 2) y := by
    intro y
    have h1 : HasDerivAt (fun y : ℝ => -Real.sin (y * x) * x) (-(Real.cos (y * x) * x) * x) y :=
      (hasDerivAt_mul_const x).sin.neg.mul_const x
    have h2 : HasDerivAt (fun y : ℝ => x ^ 2 * y) (x ^ 2) y := by
      simpa using (hasDerivAt_id y).const_mul (x ^ 2)
    convert h1.sub h2 using 1
    ring
  have hD : deriv (fun y : ℝ => Real.cos (y * x) - x ^ 2 * y ^ 2 / 2)
      = fun y => -Real.sin (y * x) * x - x ^ 2 * y := funext fun y => (hd y).deriv
  refine concaveOn_univ_of_deriv2_nonpos (fun y => (hd y).differentiableAt) ?_ ?_
  · rw [hD]; exact fun y => (hd2 y).differentiableAt
  · intro y
    rw [Function.iterate_succ, Function.iterate_one, Function.comp_apply, hD, (hd2 y).deriv]
    have := Real.neg_one_le_cos (y * x)
    nlinarith [sq_nonneg x]

/-- Pointwise chord bound for `cos(yx)`. -/
lemma cos_chord (x lo hi θ : ℝ) (hθ0 : 0 ≤ θ) (hθ1 : θ ≤ 1) :
    (1 - θ) * Real.cos (lo * x) + θ * Real.cos (hi * x) - x ^ 2 * (θ * (1 - θ) * (hi - lo) ^ 2) / 2
      ≤ Real.cos (((1 - θ) * lo + θ * hi) * x) := by
  have h := (concave_cos_sub x).2 (Set.mem_univ lo) (Set.mem_univ hi) (by linarith : 0 ≤ 1 - θ) hθ0
    (by ring)
  simp only [smul_eq_mul] at h
  nlinarith

/-- **Chord bound for `J`.** -/
theorem J_chord {a lo hi θ : ℝ} (ha : 0 ≤ a) (hθ0 : 0 ≤ θ) (hθ1 : θ ≤ 1) :
    (1 - θ) * J a lo + θ * J a hi - 16 / 3 * a ^ 3 * Real.cosh a * (θ * (1 - θ) * (hi - lo) ^ 2)
      ≤ J a ((1 - θ) * lo + θ * hi) := by
  rw [← integral_J, ← integral_J, ← integral_J]
  set q := θ * (1 - θ) * (hi - lo) ^ 2
  have hq : 0 ≤ q := by
    have : 0 ≤ θ * (1 - θ) := mul_nonneg hθ0 (by linarith)
    positivity
  have hab : -(2 * a) ≤ 2 * a := by linarith
  have c1 : Continuous fun x : ℝ => 2 * Real.cosh (x / 2) * Real.cos (lo * x) := by fun_prop
  have c2 : Continuous fun x : ℝ => 2 * Real.cosh (x / 2) * Real.cos (hi * x) := by fun_prop
  have c4 : Continuous fun x : ℝ => Real.cosh a * x ^ 2 * q := by fun_prop
  have hmono : ∫ x in (-(2 * a))..(2 * a), ((1 - θ) * (2 * Real.cosh (x / 2) * Real.cos (lo * x))
        + θ * (2 * Real.cosh (x / 2) * Real.cos (hi * x)) - Real.cosh a * x ^ 2 * q)
      ≤ ∫ x in (-(2 * a))..(2 * a), 2 * Real.cosh (x / 2) * Real.cos (((1 - θ) * lo + θ * hi) * x) := by
    refine intervalIntegral.integral_mono_on hab (Continuous.intervalIntegrable (by fun_prop) _ _)
      (Continuous.intervalIntegrable (by fun_prop) _ _) fun x hx => ?_
    have hc := cos_chord x lo hi θ hθ0 hθ1
    have hch0 : 0 < Real.cosh (x / 2) := Real.cosh_pos _
    have hcha : Real.cosh (x / 2) ≤ Real.cosh a := by
      rw [← Real.cosh_abs (x / 2)]
      apply Real.cosh_le_cosh.2
      rw [abs_abs, abs_of_nonneg ha, abs_le]
      constructor <;> linarith [hx.1, hx.2]
    have k1 := mul_le_mul_of_nonneg_left hc (by positivity : (0 : ℝ) ≤ 2 * Real.cosh (x / 2))
    have k2 : 2 * Real.cosh (x / 2) * (x ^ 2 * q / 2) ≤ Real.cosh a * x ^ 2 * q := by
      have : 0 ≤ x ^ 2 * q := by positivity
      nlinarith
    nlinarith
  refine le_trans (le_of_eq ?_) hmono
  have e2 : (fun x : ℝ => Real.cosh a * x ^ 2 * q) = fun x => (Real.cosh a * q) * x ^ 2 := by
    funext x; ring
  rw [intervalIntegral.integral_sub, intervalIntegral.integral_add, intervalIntegral.integral_const_mul,
    intervalIntegral.integral_const_mul, e2, intervalIntegral.integral_const_mul, integral_pow]
  · ring
  all_goals exact Continuous.intervalIntegrable (by fun_prop) _ _

end Fayed.Yoshida
