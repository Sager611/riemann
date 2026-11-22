module

public import Riemann.Deep.LeanLi.ChangeVar
public import Mathlib.MeasureTheory.Integral.IntervalIntegral.IntegrationByParts

/-!
# LeanLi: a Bessel inequality for the sines `sin(nθ(t))` on `[T, X]`

`Σ_{n∈s} (∫_T^X h sin nθ)² ≤ (π/2) ∫_T^X h² (t² + 1/4)` (orthogonality of `sin nv` on `(0, π)` after
`v = θ(t)`).
-/

@[expose] public section

open Real Set MeasureTheory Finset

namespace LeanLi

lemma integral_cos_int_mul (k : ℤ) (hk : k ≠ 0) : ∫ v in (0 : ℝ)..π, cos (k * v) = 0 := by
  have hk0 : (k : ℝ) ≠ 0 := by exact_mod_cast hk
  rw [intervalIntegral.integral_comp_mul_left (fun x => cos x) hk0, integral_cos]
  simp [Real.sin_int_mul_pi]

lemma integral_sin_mul_sin (n m : ℕ) (hn : 1 ≤ n) (hm : 1 ≤ m) :
    ∫ v in (0 : ℝ)..π, sin (n * v) * sin (m * v) = if n = m then π / 2 else 0 := by
  have e : ∀ v : ℝ, sin (n * v) * sin (m * v) =
      (cos (((n : ℤ) - m : ℤ) * v) - cos (((n : ℤ) + m : ℤ) * v)) / 2 := by
    intro v; push_cast
    rw [show ((n : ℝ) - m) * v = n * v - m * v by ring, show ((n : ℝ) + m) * v = n * v + m * v by ring,
      cos_sub, cos_add]; ring
  simp_rw [e]
  have hc : ∀ k : ℤ, IntervalIntegrable (fun v : ℝ => cos (k * v)) volume 0 π :=
    fun k => Continuous.intervalIntegrable (by fun_prop) _ _
  rw [intervalIntegral.integral_div, intervalIntegral.integral_sub (hc _) (hc _)]
  have hsum : ((n : ℤ) + m) ≠ 0 := by omega
  rw [integral_cos_int_mul _ hsum]
  split_ifs with h
  · subst h; simp
  · have hd : ((n : ℤ) - m) ≠ 0 := by omega
    rw [integral_cos_int_mul _ hd]; simp

lemma integral_sq_sin_sum (s : Finset ℕ) (hs : ∀ n ∈ s, 1 ≤ n) (a : ℕ → ℝ) :
    ∫ v in (0 : ℝ)..π, (∑ n ∈ s, a n * sin (n * v)) ^ 2 = π / 2 * ∑ n ∈ s, a n ^ 2 := by
  have e : ∀ v : ℝ, (∑ n ∈ s, a n * sin (n * v)) ^ 2 =
      ∑ n ∈ s, ∑ m ∈ s, a n * a m * (sin (n * v) * sin (m * v)) := by
    intro v; rw [sq, Finset.sum_mul_sum]
    refine Finset.sum_congr rfl fun n _ => Finset.sum_congr rfl fun m _ => by ring
  simp_rw [e]
  have hi : ∀ n m : ℕ, IntervalIntegrable (fun v : ℝ => a n * a m * (sin (n * v) * sin (m * v))) volume 0 π :=
    fun n m => Continuous.intervalIntegrable (by fun_prop) _ _
  rw [intervalIntegral.integral_finsetSum fun n _ => Continuous.intervalIntegrable (by fun_prop) _ _]
  · rw [Finset.mul_sum]
    refine Finset.sum_congr rfl fun n hn => ?_
    rw [intervalIntegral.integral_finsetSum fun m _ => hi n m]
    simp_rw [intervalIntegral.integral_const_mul]
    rw [Finset.sum_eq_single n]
    · rw [integral_sin_mul_sin n n (hs n hn) (hs n hn), if_pos rfl]; ring
    · intro m hm hmn
      rw [integral_sin_mul_sin n m (hs n hn) (hs m hm), if_neg (Ne.symm hmn), mul_zero]
    · intro h; exact absurd hn h

lemma continuous_inv_sq_add : Continuous fun t : ℝ => 1 / (t ^ 2 + 1 / 4) :=
  Continuous.div continuous_const (by fun_prop) fun t => by positivity

/-- **Bessel's inequality in `t`.** -/
theorem bessel_t {T X : ℝ} (hT : 0 ≤ T) (hTX : T ≤ X) (h : ℝ → ℝ)
    (hh : IntervalIntegrable h volume T X)
    (hh2 : IntervalIntegrable (fun t => h t ^ 2 * (t ^ 2 + 1 / 4)) volume T X)
    (s : Finset ℕ) (hs : ∀ n ∈ s, 1 ≤ n) :
    ∑ n ∈ s, (∫ t in T..X, h t * sin (n * liTheta t)) ^ 2 ≤
      π / 2 * ∫ t in T..X, h t ^ 2 * (t ^ 2 + 1 / 4) := by
  have hπ := pi_pos
  set a : ℕ → ℝ := fun n => ∫ t in T..X, h t * sin (n * liTheta t) with ha
  set A := ∑ n ∈ s, a n ^ 2 with hA
  set H := ∫ t in T..X, h t ^ 2 * (t ^ 2 + 1 / 4) with hH
  set P : ℝ → ℝ := fun t => ∑ n ∈ s, a n * sin (n * liTheta t) with hP
  have hcP : Continuous P := by
    have : Continuous liTheta := continuous_iff_continuousAt.2 fun t => (hasDerivAt_liTheta t).continuousAt
    simp only [hP]; fun_prop
  -- A = ∫ h P
  have hsin : ∀ n : ℕ, IntervalIntegrable (fun t => h t * sin (n * liTheta t)) volume T X := by
    intro n
    have : Continuous liTheta := continuous_iff_continuousAt.2 fun t => (hasDerivAt_liTheta t).continuousAt
    exact hh.mul_continuousOn (Continuous.continuousOn (by fun_prop))
  have hAP : A = ∫ t in T..X, h t * P t := by
    simp only [hA, hP, Finset.mul_sum]
    rw [intervalIntegral.integral_finsetSum fun n _ => by
      have := (hsin n).const_mul (a n)
      exact this.congr_uIoo fun t _ => by simp only; ring]
    refine Finset.sum_congr rfl fun n _ => ?_
    rw [sq]
    conv_rhs => rw [show (fun t => h t * (a n * sin (n * liTheta t))) =
      fun t => a n * (h t * sin (n * liTheta t)) from funext fun t => by ring]
    rw [intervalIntegral.integral_const_mul]
  -- ∫ |θ'| P² ≤ (π/2) A
  have hsub : ∫ t in T..X, (P t) ^ 2 * (1 / (t ^ 2 + 1 / 4)) ≤ π / 2 * A := by
    have hcv := intervalIntegral.integral_comp_mul_deriv (a := T) (b := X) (f := liTheta)
      (f' := fun t => -(1 / (t ^ 2 + 1 / 4))) (g := fun v => (∑ n ∈ s, a n * sin (n * v)) ^ 2)
      (fun t _ => hasDerivAt_liTheta t) (continuous_inv_sq_add.neg.continuousOn) (by fun_prop)
    have e1 : ∫ t in T..X, (P t) ^ 2 * (1 / (t ^ 2 + 1 / 4)) =
        -∫ t in T..X, ((fun v => (∑ n ∈ s, a n * sin (n * v)) ^ 2) ∘ liTheta) t * -(1 / (t ^ 2 + 1 / 4)) := by
      rw [← intervalIntegral.integral_neg]
      refine intervalIntegral.integral_congr fun t _ => ?_
      simp only [hP, Function.comp]; ring
    rw [e1, hcv, ← intervalIntegral.integral_symm]
    have hθ1 : liTheta X ≤ liTheta T := liTheta_strictAnti.antitone hTX
    have hθ0 : 0 ≤ liTheta X := by
      rcases eq_or_lt_of_le (le_trans hT hTX) with h0 | h0
      · rw [← h0]; simp [liTheta, pi_pos.le]
      · exact (liTheta_pos h0).le
    have hθπ : liTheta T ≤ π := by
      rcases eq_or_lt_of_le hT with h0 | h0
      · rw [← h0]; simp [liTheta]
      · exact (liTheta_lt_pi h0).le
    have hmono : ∫ v in liTheta X..liTheta T, (∑ n ∈ s, a n * sin (n * v)) ^ 2 ≤
        ∫ v in (0 : ℝ)..π, (∑ n ∈ s, a n * sin (n * v)) ^ 2 :=
      intervalIntegral.integral_mono_interval hθ0 hθ1 hθπ
        (Filter.Eventually.of_forall fun _ => sq_nonneg _) (Continuous.intervalIntegrable (by fun_prop) _ _)
    rw [integral_sq_sin_sum s hs a] at hmono
    exact hmono
  -- AM-GM
  have hpt : ∀ t ∈ Icc T X, h t * P t ≤ π / 4 * (h t ^ 2 * (t ^ 2 + 1 / 4)) +
      1 / π * ((P t) ^ 2 * (1 / (t ^ 2 + 1 / 4))) := by
    intro t _
    have hw : 0 < t ^ 2 + 1 / 4 := by positivity
    have key : 0 ≤ (π / 2 * h t * (t ^ 2 + 1 / 4) - P t) ^ 2 := sq_nonneg _
    have e : π / 4 * (h t ^ 2 * (t ^ 2 + 1 / 4)) + 1 / π * ((P t) ^ 2 * (1 / (t ^ 2 + 1 / 4))) - h t * P t =
        (π / 2 * h t * (t ^ 2 + 1 / 4) - P t) ^ 2 / (π * (t ^ 2 + 1 / 4)) := by
      field_simp; ring
    have : 0 ≤ (π / 2 * h t * (t ^ 2 + 1 / 4) - P t) ^ 2 / (π * (t ^ 2 + 1 / 4)) := by positivity
    linarith
  have hi1 : IntervalIntegrable (fun t => h t * P t) volume T X := hh.mul_continuousOn hcP.continuousOn
  have hi2 : IntervalIntegrable (fun t => (P t) ^ 2 * (1 / (t ^ 2 + 1 / 4))) volume T X :=
    Continuous.intervalIntegrable ((hcP.pow 2).mul continuous_inv_sq_add) _ _
  have hint := intervalIntegral.integral_mono_on hTX hi1 ((hh2.const_mul (π / 4)).add (hi2.const_mul (1 / π))) hpt
  rw [intervalIntegral.integral_add (hh2.const_mul _) (hi2.const_mul _), intervalIntegral.integral_const_mul,
    intervalIntegral.integral_const_mul, ← hAP] at hint
  have : 1 / π * ∫ t in T..X, (P t) ^ 2 * (1 / (t ^ 2 + 1 / 4)) ≤ A / 2 := by
    have := mul_le_mul_of_nonneg_left hsub (by positivity : (0 : ℝ) ≤ 1 / π)
    have e : 1 / π * (π / 2 * A) = A / 2 := by field_simp
    linarith
  have hA' : A ≤ π / 4 * H * 2 := by linarith
  show A ≤ π / 2 * H
  linarith

end LeanLi
