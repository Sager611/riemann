module

public import Riemann.Deep.LeanLi.TrendI
public import Riemann.Deep.LeanLi.LiTheta
public import Mathlib.MeasureTheory.Function.JacobianOneDim

/-!
# LeanLi: `I(n) = π⁻¹ ∫₀^∞ (1 - cos nθ(t)) log(t/2π) dt = Σ_{k odd<n} (n-k)/k - (n/2) log 4π`

The substitution `v = θ(t)`, `θ(t) = π - 2 arctan 2t`, maps `(0, ∞)` onto `(0, π)` with
`|θ'(t)| = 1/(t² + 1/4)` and `|θ'| a₀(θ(t)) = (2π)⁻¹ log(t/2π)`.
-/

@[expose] public section

open Real Set MeasureTheory

namespace LeanLi

lemma hasDerivAt_liTheta (t : ℝ) : HasDerivAt liTheta (-(1 / (t ^ 2 + 1 / 4))) t := by
  have h := ((hasDerivAt_id t).const_mul 2).arctan
  have := (h.const_mul 2).const_sub Real.pi
  unfold liTheta
  convert this using 1
  · funext x; simp
  · simp only [id]; field_simp; ring

lemma liTheta_strictAnti : StrictAnti liTheta := by
  intro a b hab
  unfold liTheta
  have := arctan_strictMono (show 2 * a < 2 * b by linarith)
  linarith

lemma liTheta_image : liTheta '' Ioi 0 = Ioo 0 π := by
  ext v
  constructor
  · rintro ⟨t, ht, rfl⟩
    exact ⟨liTheta_pos ht, liTheta_lt_pi ht⟩
  · rintro ⟨h0, h1⟩
    refine ⟨Real.cot (v / 2) / 2, ?_, ?_⟩
    · have hs : 0 < Real.sin (v / 2) := sin_pos_of_pos_of_lt_pi (by linarith) (by linarith)
      have hc : 0 < Real.cos (v / 2) := cos_pos_of_mem_Ioo ⟨by linarith, by linarith⟩
      rw [Real.cot_eq_cos_div_sin]
      show 0 < Real.cos (v / 2) / Real.sin (v / 2) / 2
      positivity
    · unfold liTheta
      have e : 2 * (Real.cot (v / 2) / 2) = Real.tan (π / 2 - v / 2) := by
        rw [Real.tan_pi_div_two_sub, Real.cot_eq_cos_div_sin, Real.tan_eq_sin_div_cos, inv_div]; ring
      rw [e, arctan_tan (by linarith) (by linarith)]
      ring

/-- `|θ'(t)| a₀(θ(t)) = (2π)⁻¹ log(t/2π)`. -/
lemma a0_liTheta {t : ℝ} (ht : 0 < t) :
    1 / (t ^ 2 + 1 / 4) * a0 (liTheta t) = 1 / (2 * π) * Real.log (t / (2 * π)) := by
  have hπ := pi_pos
  set y : ℝ := 1 / (2 * t) with hy
  have hy0 : 0 < y := by positivity
  have hθ : liTheta t / 2 = arctan y := by rw [liTheta_eq ht]; ring
  have hp : 0 < 1 + y ^ 2 := by positivity
  have hsq : √(1 + y ^ 2) ^ 2 = 1 + y ^ 2 := Real.sq_sqrt hp.le
  have hs0 : 0 < √(1 + y ^ 2) := Real.sqrt_pos.2 hp
  have hsin : Real.sin (liTheta t / 2) = y / √(1 + y ^ 2) := by rw [hθ, sin_arctan]
  have hcos : Real.cos (liTheta t / 2) = 1 / √(1 + y ^ 2) := by rw [hθ, cos_arctan]
  have hell : ell (liTheta t) = Real.log (t / (2 * π)) := by
    unfold ell
    rw [hsin, hcos, ← Real.log_div (by positivity) (by positivity), ← Real.log_div (by positivity) (by positivity)]
    congr 1
    rw [hy]; field_simp; ring
  unfold a0
  rw [hell, hsin, div_pow, hsq]
  rw [hy]
  field_simp
  ring

/-- **`I(n)` as a `t`-integral**: `∫_{t>0} π⁻¹(1 - cos nθ(t)) log(t/2π) dt = Iclosed n`. -/
theorem integral_I_eq (n : ℕ) :
    ∫ t in Ioi (0 : ℝ), 1 / π * ((1 - Real.cos (n * liTheta t)) * Real.log (t / (2 * π))) = Iclosed n := by
  have hπ := pi_pos
  have hd : ∀ t ∈ Ioi (0 : ℝ), HasDerivWithinAt liTheta (-(1 / (t ^ 2 + 1 / 4))) (Ioi 0) t :=
    fun t _ => (hasDerivAt_liTheta t).hasDerivWithinAt
  have hinj : InjOn liTheta (Ioi 0) := liTheta_strictAnti.injective.injOn
  have h := integral_image_eq_integral_abs_deriv_smul measurableSet_Ioi hd hinj
    (fun v => 2 * (1 - Real.cos (n * v)) * a0 v)
  rw [liTheta_image] at h
  rw [← Imod_eq_closed]
  unfold Imod
  rw [intervalIntegral.integral_of_le hπ.le, integral_Ioc_eq_integral_Ioo, h]
  refine setIntegral_congr_fun measurableSet_Ioi fun t ht => ?_
  simp only [smul_eq_mul]
  have ht0 : 0 < t := ht
  rw [abs_neg, abs_of_pos (by positivity : (0 : ℝ) < 1 / (t ^ 2 + 1 / 4))]
  have := a0_liTheta ht0
  symm
  calc 1 / (t ^ 2 + 1 / 4) * (2 * (1 - Real.cos (n * liTheta t)) * a0 (liTheta t))
      = 2 * (1 - Real.cos (n * liTheta t)) * (1 / (t ^ 2 + 1 / 4) * a0 (liTheta t)) := by ring
    _ = 2 * (1 - Real.cos (n * liTheta t)) * (1 / (2 * π) * Real.log (t / (2 * π))) := by rw [this]
    _ = 1 / π * ((1 - Real.cos (n * liTheta t)) * Real.log (t / (2 * π))) := by field_simp

end LeanLi
