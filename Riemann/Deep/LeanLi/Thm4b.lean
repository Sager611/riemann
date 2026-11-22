module

public import Riemann.Deep.LeanLi.Thm4a

/-!
# LeanLi: Theorem 4 ingredients from the literature inputs:
`∫_0^X r = -1/8 + O(1/X)`, `∫ (2t+1)|r| ≪ log X`, `sup_{[0,X]} |S₁| ≪ log X`
-/

@[expose] public section

open Real Filter Set Finset MeasureTheory Topology

namespace LeanLi

lemma integral_log_div_two_pi (X : ℝ) :
    ∫ t in (0 : ℝ)..X, Real.log (t / (2 * Real.pi)) = X * Real.log (X / (2 * Real.pi)) - X := by
  have hc : (2 * Real.pi) ≠ 0 := by positivity
  rw [intervalIntegral.integral_comp_div (fun t => Real.log t) hc, integral_log, smul_eq_mul]
  simp only [zero_div, Real.log_zero, mul_zero, sub_zero, add_zero]
  field_simp

/-- `∫_0^X r = (θ(X) - (X/2)(log(X/2π) - 1))/π`. -/
lemma integral_rfun_eq (X : ℝ) (hX : 0 ≤ X) :
    ∫ t in (0 : ℝ)..X, rfun t = (rsTheta X - X / 2 * (Real.log (X / (2 * Real.pi)) - 1)) / Real.pi := by
  have hπ := Real.pi_pos
  unfold rfun
  rw [intervalIntegral.integral_sub (continuous_Ad.intervalIntegrable _ _)
    ((intervalIntegrable_log_div_two_pi le_rfl hX).const_mul _), integral_Ad, Afun_zero,
    intervalIntegral.integral_const_mul, integral_log_div_two_pi]
  unfold Afun
  field_simp
  ring

section
variable (L : LiteratureInputs)
include L

/-- **`∫_0^∞ r = -1/8`** (the `-π/8` of `θ_RS`), with rate. -/
lemma rfun_tail : ∃ C X₀ : ℝ, 1 ≤ X₀ ∧ ∀ X, X₀ ≤ X → |(∫ t in (0 : ℝ)..X, rfun t) + 1 / 8| ≤ C / X := by
  obtain ⟨C, t₀, hC⟩ := L.rs_theta
  refine ⟨|C|, max t₀ 1, le_max_right _ _, fun X hX => ?_⟩
  have hX1 : 1 ≤ X := le_trans (le_max_right _ _) hX
  have hX0 : 0 < X := by linarith
  have hπ := Real.pi_pos
  have hπ1 : 1 ≤ Real.pi := by linarith [Real.pi_gt_three]
  rw [integral_rfun_eq X hX0.le]
  have h := hC X (le_trans (le_max_left _ _) hX)
  have e : (rsTheta X - X / 2 * (Real.log (X / (2 * Real.pi)) - 1)) / Real.pi + 1 / 8 =
      (rsTheta X - (X / 2 * (Real.log (X / (2 * Real.pi)) - 1) - Real.pi / 8)) / Real.pi := by
    field_simp; ring
  rw [e, abs_div, abs_of_pos hπ, div_le_div_iff₀ hπ hX0]
  have h2 : C / X ≤ |C| / X := div_le_div_of_nonneg_right (le_abs_self C) hX0.le
  have h3 := h.trans h2
  rw [le_div_iff₀ hX0] at h3
  nlinarith [abs_nonneg C]

/-- `∫_{T₀}^X (2t+1)|r(t)| dt ≤ K₁ + K₂ log X`. -/
lemma int_weighted_r : ∃ K₁ K₂ X₀ : ℝ, 1 ≤ X₀ ∧ ∀ T₀ X : ℝ, 0 ≤ T₀ → T₀ ≤ 1 → X₀ ≤ X →
    ∫ t in T₀..X, (2 * t + 1) * |rfun t| ≤ K₁ + K₂ * Real.log X := by
  obtain ⟨C, t₀, ht₀, hC⟩ := rfun_bound L
  obtain ⟨Kr, hKr⟩ := rfun_L1 L
  refine ⟨(2 * t₀ + 1) * |Kr|, 3 * |C|, t₀, ht₀, fun T₀ X hT₀ hT₀1 hX => ?_⟩
  have ht0 : 0 < t₀ := by linarith
  have hX0 : 0 < X := by linarith
  have hc : Continuous fun t : ℝ => 2 * t + 1 := by fun_prop
  have hi : ∀ a b, 0 ≤ a → 0 ≤ b → IntervalIntegrable (fun t => (2 * t + 1) * |rfun t|) volume a b :=
    fun a b ha hb => (intervalIntegrable_rfun ha hb).abs.continuousOn_mul hc.continuousOn
  rw [← intervalIntegral.integral_add_adjacent_intervals (hi T₀ t₀ hT₀ ht0.le) (hi t₀ X ht0.le hX0.le)]
  have hTt : T₀ ≤ t₀ := by linarith
  have h1 : ∫ t in T₀..t₀, (2 * t + 1) * |rfun t| ≤ (2 * t₀ + 1) * |Kr| := by
    have hm : ∫ t in T₀..t₀, (2 * t + 1) * |rfun t| ≤ ∫ t in T₀..t₀, (2 * t₀ + 1) * |rfun t| := by
      apply intervalIntegral.integral_mono_on hTt (hi T₀ t₀ hT₀ ht0.le)
        ((intervalIntegrable_rfun hT₀ ht0.le).abs.const_mul _)
      intro t ht
      exact mul_le_mul_of_nonneg_right (by linarith [ht.2]) (abs_nonneg _)
    refine hm.trans ?_
    rw [intervalIntegral.integral_const_mul]
    apply mul_le_mul_of_nonneg_left _ (by linarith)
    have h0 : ∫ t in T₀..t₀, |rfun t| ≤ ∫ t in (0 : ℝ)..t₀, |rfun t| :=
      intervalIntegral.integral_mono_interval hT₀ hTt le_rfl
        (Filter.Eventually.of_forall fun t => abs_nonneg (rfun t)) (intervalIntegrable_rfun le_rfl ht0.le).abs
    exact h0.trans ((hKr t₀ ht0.le).trans (le_abs_self _))
  have h2 : ∫ t in t₀..X, (2 * t + 1) * |rfun t| ≤ 3 * |C| * Real.log X := by
    have hm : ∫ t in t₀..X, (2 * t + 1) * |rfun t| ≤ ∫ t in t₀..X, 3 * |C| * t⁻¹ := by
      apply intervalIntegral.integral_mono_on hX (hi t₀ X ht0.le hX0.le)
      · apply ContinuousOn.intervalIntegrable
        apply ContinuousOn.mul continuousOn_const
        apply ContinuousOn.inv₀ continuousOn_id
        intro t ht; rw [Set.uIcc_of_le hX] at ht; simp only [id]; linarith [ht.1]
      intro t ht
      have htp : 0 < t := by linarith [ht.1]
      have h3 := (hC t ht.1).trans (div_le_div_of_nonneg_right (le_abs_self C) (by positivity))
      have h4 : (2 * t + 1) * |rfun t| ≤ (3 * t) * (|C| / t ^ 2) :=
        mul_le_mul (by linarith [ht.1]) h3 (abs_nonneg _) (by linarith)
      have e : (3 * t) * (|C| / t ^ 2) = 3 * |C| * t⁻¹ := by field_simp
      linarith
    refine hm.trans ?_
    rw [intervalIntegral.integral_const_mul, integral_inv_of_pos ht0 hX0]
    apply mul_le_mul_of_nonneg_left _ (by positivity)
    rw [Real.log_div hX0.ne' ht0.ne']
    have := Real.log_nonneg (show (1 : ℝ) ≤ t₀ by linarith)
    linarith
  linarith

/-- `|S₁(t)| ≤ K₁ + K₂ log X` on `[0, X]`. -/
lemma S1_sup : ∃ K₁ K₂ X₀ : ℝ, 1 ≤ X₀ ∧ 0 ≤ K₁ ∧ 0 ≤ K₂ ∧ ∀ X, X₀ ≤ X → ∀ t, 0 ≤ t → t ≤ X →
    |S1 t| ≤ K₁ + K₂ * Real.log X := by
  obtain ⟨C, T₁, hC⟩ := L.s1
  set t₁ := max T₁ 1 with ht₁
  have ht₁1 : 1 ≤ t₁ := le_max_right _ _
  refine ⟨∫ u in (0 : ℝ)..t₁, |zetaS u|, |C|, t₁, ht₁1,
    intervalIntegral.integral_nonneg (by linarith) (fun u _ => abs_nonneg _), abs_nonneg _,
    fun X hX t ht0 htX => ?_⟩
  have hK : 0 ≤ ∫ u in (0 : ℝ)..t₁, |zetaS u| := intervalIntegral.integral_nonneg (by linarith) (fun u _ => abs_nonneg _)
  have hlogX : 0 ≤ Real.log X := Real.log_nonneg (by linarith)
  rcases le_or_gt t₁ t with h | h
  · have h1 := hC t (le_trans (le_max_left _ _) h)
    have hlt : Real.log t ≤ Real.log X := Real.log_le_log (by linarith) htX
    have hlt0 : 0 ≤ Real.log t := Real.log_nonneg (by linarith)
    have : C * Real.log t ≤ |C| * Real.log X :=
      (mul_le_mul_of_nonneg_right (le_abs_self C) hlt0).trans (mul_le_mul_of_nonneg_left hlt (abs_nonneg _))
    unfold S1
    linarith
  · have h1 : |S1 t| ≤ ∫ u in (0 : ℝ)..t, |zetaS u| := intervalIntegral.abs_integral_le_integral_abs ht0
    have h2 : ∫ u in (0 : ℝ)..t, |zetaS u| ≤ ∫ u in (0 : ℝ)..t₁, |zetaS u| :=
      intervalIntegral.integral_mono_interval le_rfl ht0 h.le
        (Filter.Eventually.of_forall fun u => abs_nonneg _) (intervalIntegrable_zetaS 0 t₁).abs
    have : 0 ≤ |C| * Real.log X := by positivity
    linarith

end

end LeanLi

