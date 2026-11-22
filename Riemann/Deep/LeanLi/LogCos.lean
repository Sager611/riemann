module

public import Mathlib.Analysis.SpecialFunctions.Integrals.LogTrigonometric
public import Mathlib.Analysis.SpecialFunctions.Trigonometric.Deriv
public import Mathlib.Analysis.SpecialFunctions.Trigonometric.Bounds
public import Mathlib.Analysis.SpecialFunctions.Log.NegMulLog
public import Mathlib.MeasureTheory.Integral.IntervalIntegral.IntegrationByParts
public import Mathlib.MeasureTheory.Integral.DominatedConvergence

/-!
# LeanLi: the cosine–log integrals `∫_0^π cos(jv) log sin(v/2) dv = -π/(2j)`
-/

@[expose] public section

open Real Filter Topology Set MeasureTheory

namespace LeanLi

/-- Dirichlet kernel: `sin(v/2)(1 + 2Σ_{k=1}^{j-1} cos kv) = sin((2j-1)v/2)`. -/
lemma dirichlet_sin (v : ℝ) (j : ℕ) (hj : 1 ≤ j) :
    sin (v / 2) * (1 + 2 * ∑ k ∈ Finset.Ico 1 j, cos (k * v)) = sin ((2 * j - 1) * v / 2) := by
  induction j, hj using Nat.le_induction with
  | base => norm_num
  | succ j hj ih =>
    rw [Finset.sum_Ico_succ_top hj]
    have h2 : 2 * sin (v / 2) * cos (j * v) = sin ((2 * j + 1) * v / 2) - sin ((2 * j - 1) * v / 2) := by
      have := two_mul_sin_mul_cos (v / 2) (j * v)
      rw [show v / 2 - j * v = -((2 * j - 1) * v / 2) by ring, sin_neg,
        show v / 2 + j * v = (2 * j + 1) * v / 2 by ring] at this
      linarith
    push_cast
    rw [show (2 * ((j : ℝ) + 1) - 1) * v / 2 = (2 * j + 1) * v / 2 by ring]
    linear_combination ih + h2

lemma sin_mul_cos_half (v : ℝ) (j : ℕ) (hj : 1 ≤ j) :
    sin (j * v) * cos (v / 2) = sin (v / 2) * (1 + cos (j * v) + 2 * ∑ k ∈ Finset.Ico 1 j, cos (k * v)) := by
  have hD := dirichlet_sin v j hj
  have h1 : sin (v / 2) * cos (j * v) = (sin ((2 * j + 1) * v / 2) - sin ((2 * j - 1) * v / 2)) / 2 := by
    have := two_mul_sin_mul_cos (v / 2) (j * v)
    rw [show v / 2 - j * v = -((2 * j - 1) * v / 2) by ring, sin_neg,
      show v / 2 + j * v = (2 * j + 1) * v / 2 by ring] at this
    linarith
  have h2 : sin (j * v) * cos (v / 2) = (sin ((2 * j + 1) * v / 2) + sin ((2 * j - 1) * v / 2)) / 2 := by
    have := two_mul_sin_mul_cos (j * v) (v / 2)
    rw [show j * v - v / 2 = (2 * j - 1) * v / 2 by ring,
      show j * v + v / 2 = (2 * j + 1) * v / 2 by ring] at this
    linarith
  have e : sin (v / 2) * (1 + cos (j * v) + 2 * ∑ k ∈ Finset.Ico 1 j, cos (k * v)) =
      sin (v / 2) * (1 + 2 * ∑ k ∈ Finset.Ico 1 j, cos (k * v)) + sin (v / 2) * cos (j * v) := by ring
  rw [e, hD, h1, h2]; ring

lemma alg_aux (c s S D j : ℝ) (hs : s ≠ 0) (hj : j ≠ 0) (h : S * c = s * D) :
    c * (1 / 2) / s * (S / j) = 1 / (2 * j) * D := by
  field_simp
  linear_combination h

/-- `K_j = ∫_0^π cos(jv) log sin(v/2) dv = -π/(2j)` for `j ≥ 1`. -/
theorem integral_cos_mul_log_sin_half (j : ℕ) (hj : 1 ≤ j) :
    ∫ v in (0 : ℝ)..π, cos (j * v) * log (sin (v / 2)) = -π / (2 * j) := by
  set F : ℝ → ℝ := fun v => cos (j * v) * log (sin (v / 2)) with hF
  have hj0 : (0 : ℝ) < j := by exact_mod_cast hj
  have hlogint : ∀ a b : ℝ, IntervalIntegrable (fun v => log (sin (v / 2))) volume a b := by
    intro a b
    have := (intervalIntegrable_log_sin (a := a / 2) (b := b / 2)).comp_mul_left (c := 1 / 2)
    simpa [div_eq_mul_inv, mul_comm] using this
  have hFint : ∀ a b : ℝ, IntervalIntegrable F volume a b := fun a b =>
    (hlogint a b).continuousOn_mul (Continuous.continuousOn (by fun_prop))
  -- the integral over [ε, π] in closed form
  set G : ℝ → ℝ := fun ε => -(log (sin (ε / 2)) * (sin (j * ε) / j)) -
    (1 / (2 * j)) * ((π - ε) - sin (j * ε) / j - 2 * ∑ k ∈ Finset.Ico 1 j, sin (k * ε) / k) with hG
  have hIBP : ∀ ε ∈ Ioo 0 π, ∫ v in ε..π, F v = G ε := by
    intro ε hε
    have hpos : ∀ v ∈ uIcc ε π, 0 < sin (v / 2) := by
      intro v hv; rw [uIcc_of_le hε.2.le] at hv
      exact sin_pos_of_pos_of_lt_pi (by linarith [hv.1, hε.1]) (by linarith [hv.2, pi_pos])
    have hu : ∀ v ∈ uIcc ε π, HasDerivAt (fun v => log (sin (v / 2)))
        (cos (v / 2) * (1 / 2) / sin (v / 2)) v := by
      intro v hv
      have h1 := ((hasDerivAt_id' v).div_const 2).sin
      have h2 := h1.log (hpos v hv).ne'
      convert h2 using 1
    have hv : ∀ v ∈ uIcc ε π, HasDerivAt (fun v => sin (j * v) / j) (cos (j * v)) v := by
      intro v _
      have := (((hasDerivAt_id' v).const_mul (j : ℝ)).sin).div_const (j : ℝ)
      convert this using 1
      field_simp
    have hu' : IntervalIntegrable (fun v => cos (v / 2) * (1 / 2) / sin (v / 2)) volume ε π := by
      refine ContinuousOn.intervalIntegrable ?_
      refine ContinuousOn.div (Continuous.continuousOn (by fun_prop))
        (Continuous.continuousOn (by fun_prop)) fun v hv => (hpos v hv).ne'
    have hv' : IntervalIntegrable (fun v => cos (j * v)) volume ε π :=
      (Continuous.continuousOn (by fun_prop)).intervalIntegrable
    have h := intervalIntegral.integral_mul_deriv_eq_deriv_mul hu hv hu' hv'
    have hF' : ∫ v in ε..π, F v = ∫ v in ε..π, log (sin (v / 2)) * cos (j * v) := by
      refine intervalIntegral.integral_congr fun v _ => ?_
      simp only [hF]; ring
    rw [hF', h]
    have hD : ∫ v in ε..π, cos (v / 2) * (1 / 2) / sin (v / 2) * (sin (j * v) / j) =
        (1 / (2 * j)) * ∫ v in ε..π, (1 + cos (j * v) + 2 * ∑ k ∈ Finset.Ico 1 j, cos (k * v)) := by
      rw [← intervalIntegral.integral_const_mul]
      refine intervalIntegral.integral_congr fun v hv => ?_
      have hs := (hpos v hv).ne'
      exact alg_aux _ _ _ _ _ hs hj0.ne' (sin_mul_cos_half v j hj)
    have hcalc : ∫ v in ε..π, (1 + cos (j * v) + 2 * ∑ k ∈ Finset.Ico 1 j, cos (k * v)) =
        (π - ε) - sin (j * ε) / j - 2 * ∑ k ∈ Finset.Ico 1 j, sin (k * ε) / k := by
      have hderiv : ∀ v ∈ uIcc ε π, HasDerivAt
          (fun v => v + sin (j * v) / j + 2 * ∑ k ∈ Finset.Ico 1 j, sin (k * v) / k)
          (1 + cos (j * v) + 2 * ∑ k ∈ Finset.Ico 1 j, cos (k * v)) v := by
        intro v _
        have h1 : HasDerivAt (fun v => sin (j * v) / j) (cos (j * v)) v := by
          have := (((hasDerivAt_id' v).const_mul (j : ℝ)).sin).div_const (j : ℝ)
          convert this using 1; field_simp
        have h2 : HasDerivAt (fun v => ∑ k ∈ Finset.Ico 1 j, sin (k * v) / k)
            (∑ k ∈ Finset.Ico 1 j, cos (k * v)) v := by
          refine HasDerivAt.fun_sum fun k hk => ?_
          have hk0 : (k : ℝ) ≠ 0 := by
            have := (Finset.mem_Ico.1 hk).1; exact_mod_cast (by omega : k ≠ 0)
          have := (((hasDerivAt_id' v).const_mul (k : ℝ)).sin).div_const (k : ℝ)
          convert this using 1; field_simp
        convert ((hasDerivAt_id' v).add h1).add (h2.const_mul 2) using 1
      rw [intervalIntegral.integral_eq_sub_of_hasDerivAt hderiv
        ((Continuous.continuousOn (by fun_prop)).intervalIntegrable)]
      have hsinpi : sin (j * π) = 0 := by
        rw [show (j : ℝ) * π = (j : ℕ) * π by rfl]; exact sin_nat_mul_pi j
      have hsinpik : ∀ k ∈ Finset.Ico 1 j, sin (k * π) / k = 0 := fun k _ => by
        rw [show (k : ℝ) * π = (k : ℕ) * π by rfl, sin_nat_mul_pi]; simp
      rw [Finset.sum_congr rfl hsinpik, hsinpi]
      simp; ring
    rw [hD, hcalc]
    have hsinpi : sin (j * π) = 0 := by
      rw [show (j : ℝ) * π = (j : ℕ) * π by rfl]; exact sin_nat_mul_pi j
    simp only [hG, hsinpi, zero_div, mul_zero, zero_sub]
  -- limits as ε → 0⁺
  have hlim1 : Tendsto (fun ε => ∫ v in ε..π, F v) (𝓝[>] 0) (𝓝 (∫ v in (0 : ℝ)..π, F v)) := by
    have hc : ContinuousOn (fun b => ∫ v in (0 : ℝ)..b, F v) (uIcc 0 π) :=
      intervalIntegral.continuousOn_primitive_interval' (hFint 0 π) left_mem_uIcc
    have h0 : Tendsto (fun b => ∫ v in (0 : ℝ)..b, F v) (𝓝[>] 0) (𝓝 0) := by
      have := (hc 0 left_mem_uIcc).tendsto
      simp only [intervalIntegral.integral_same] at this
      refine this.mono_left ?_
      rw [uIcc_of_le pi_pos.le]
      exact nhdsWithin_le_of_mem (Icc_mem_nhdsGT pi_pos)
    have e : ∀ ε, ∫ v in ε..π, F v = (∫ v in (0 : ℝ)..π, F v) - ∫ v in (0 : ℝ)..ε, F v := by
      intro ε
      rw [← intervalIntegral.integral_add_adjacent_intervals (hFint 0 ε) (hFint ε π)]; ring
    rw [show (fun ε => ∫ v in ε..π, F v) = fun ε => (∫ v in (0 : ℝ)..π, F v) - ∫ v in (0 : ℝ)..ε, F v
      from funext e]
    simpa using tendsto_const_nhds.sub h0
  have hlim2 : Tendsto G (𝓝[>] 0) (𝓝 (-π / (2 * j))) := by
    set A : ℝ → ℝ := fun ε => -(log (sin (ε / 2)) * (sin (j * ε) / j)) with hA
    set B : ℝ → ℝ := fun ε => -(1 / (2 * j)) * ((π - ε) - sin (j * ε) / j -
      2 * ∑ k ∈ Finset.Ico 1 j, sin (k * ε) / k) with hB
    have hGAB : G = fun ε => A ε + B ε := by funext ε; simp only [hG, hA, hB]; ring
    have hBc : Continuous B := by fun_prop
    have hB0 : B 0 = -π / (2 * j) := by simp [hB]; field_simp
    have hBlim : Tendsto B (𝓝[>] 0) (𝓝 (-π / (2 * j))) := by
      rw [← hB0]; exact (hBc.tendsto 0).mono_left nhdsWithin_le_nhds
    have hAlim : Tendsto A (𝓝[>] 0) (𝓝 0) := by
      have hml : Tendsto (fun ε : ℝ => ε * log (ε / π)) (𝓝[>] 0) (𝓝 0) := by
        have h1 : Tendsto (fun ε : ℝ => ε / π) (𝓝[>] 0) (𝓝 0) := by
          have := ((continuous_id.div_const π).tendsto 0).mono_left (nhdsWithin_le_nhds (s := Ioi 0))
          simpa using this
        have h2 := (continuous_mul_log.tendsto 0).comp h1
        simp only [zero_mul, Real.log_zero] at h2
        have h3 := h2.const_mul π
        simp only [mul_zero] at h3
        refine h3.congr' (Eventually.of_forall fun ε => ?_)
        simp only [Function.comp]; field_simp
      refine squeeze_zero_norm' ?_ (by simpa using hml.neg)
      filter_upwards [Ioo_mem_nhdsGT pi_pos] with ε hε
      have hε0 := hε.1
      have hsin_lb : ε / π ≤ sin (ε / 2) := by
        have := Real.mul_le_sin (by linarith : (0 : ℝ) ≤ ε / 2) (by linarith [hε.2] : ε / 2 ≤ π / 2)
        have e : 2 / π * (ε / 2) = ε / π := by ring
        linarith
      have hsin_pos : 0 < sin (ε / 2) := lt_of_lt_of_le (by positivity) hsin_lb
      have hlog_le : log (sin (ε / 2)) ≤ 0 := Real.log_nonpos hsin_pos.le (sin_le_one _)
      have hlog_ge : log (ε / π) ≤ log (sin (ε / 2)) := Real.log_le_log (by positivity) hsin_lb
      have hsj : |sin (j * ε)| ≤ j * ε := by
        have := Real.abs_sin_le_abs (x := j * ε)
        rwa [abs_of_pos (mul_pos hj0 hε0)] at this
      rw [Real.norm_eq_abs, hA]
      simp only [abs_neg, abs_mul, abs_div, Nat.abs_cast]
      rw [abs_of_nonpos hlog_le]
      calc -log (sin (ε / 2)) * (|sin (j * ε)| / j) ≤ -log (ε / π) * ((j * ε) / j) := by
            apply mul_le_mul (by linarith) (div_le_div_of_nonneg_right hsj hj0.le) (by positivity)
            have : log (ε / π) ≤ 0 := Real.log_nonpos (by positivity) (by rw [div_le_one pi_pos]; linarith [hε.2])
            linarith
        _ = -(ε * log (ε / π)) := by field_simp
    rw [hGAB]
    simpa using hAlim.add hBlim
  have hev : (fun ε => ∫ v in ε..π, F v) =ᶠ[𝓝[>] 0] G := by
    filter_upwards [Ioo_mem_nhdsGT pi_pos] with ε hε using hIBP ε hε
  exact tendsto_nhds_unique (hlim1.congr' hev) hlim2

end LeanLi
