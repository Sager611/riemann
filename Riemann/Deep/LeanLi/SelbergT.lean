module

public import Riemann.Deep.LeanLi.AbelS

/-!
# LeanLi: `∫_T^X S(t)²/(t² + 1/4) dt ≤ C (log log T + 2)/T` (Selberg's mean value, dyadically)
-/

@[expose] public section

open Real Set MeasureTheory Finset

namespace LeanLi

lemma sum_geom_le (K : ℕ) : ∑ k ∈ range K, (1 / 2 : ℝ) ^ k ≤ 2 := by
  have h := geom_sum_eq (show (1 / 2 : ℝ) ≠ 1 by norm_num) K
  rw [h]
  have : 0 ≤ (1 / 2 : ℝ) ^ K := by positivity
  norm_num; linarith

lemma sum_k_geom_le (K : ℕ) : ∑ k ∈ range K, ((k : ℝ) + 1) * (1 / 2 : ℝ) ^ k ≤ 4 := by
  have key : ∀ K : ℕ, ∑ k ∈ range K, ((k : ℝ) + 1) * (1 / 2 : ℝ) ^ k = 4 - ((K : ℝ) + 2) * (1 / 2 : ℝ) ^ K / 1 * 2 := by
    intro K
    induction K with
    | zero => simp; norm_num
    | succ K ih =>
      rw [Finset.sum_range_succ, ih, pow_succ]; push_cast; ring
  rw [key]
  have : 0 ≤ ((K : ℝ) + 2) * (1 / 2 : ℝ) ^ K / 1 * 2 := by positivity
  linarith

lemma zeroCountAvg_nonneg (t : ℝ) : 0 ≤ zeroCountAvg t := by
  unfold zeroCountAvg; positivity

lemma intervalIntegrable_zetaS_sq (a b : ℝ) : IntervalIntegrable (fun t => zetaS t ^ 2) volume a b := by
  have hA : Continuous Afun := continuous_iff_continuousAt.2 fun t => (hasDerivAt_Afun t).continuousAt
  have hm2 : IntervalIntegrable (fun t => zeroCountAvg t ^ 2) volume a b := by
    have : Monotone fun t => zeroCountAvg t ^ 2 := fun s t h =>
      pow_le_pow_left₀ (zeroCountAvg_nonneg s) (zeroCountAvg_mono h) 2
    exact this.intervalIntegrable
  have hmA : IntervalIntegrable (fun t => zeroCountAvg t * Afun t) volume a b :=
    zeroCountAvg_mono.intervalIntegrable.mul_continuousOn hA.continuousOn
  have hA2 : IntervalIntegrable (fun t => Afun t ^ 2) volume a b := (hA.pow 2).intervalIntegrable _ _
  refine ((hm2.sub (hmA.const_mul 2)).add hA2).congr_uIoo fun t _ => ?_
  simp only [zetaS_eq]; ring

section
variable (L : LitCore)
include L

/-- **Selberg, dyadically.** -/
theorem selberg_tail : ∃ C T₀ : ℝ, 0 ≤ C ∧ exp 1 ≤ T₀ ∧ ∀ T X : ℝ, T₀ ≤ T → T ≤ X →
    ∫ t in T..X, zetaS t ^ 2 / (t ^ 2 + 1 / 4) ≤ C * (Real.log (Real.log T) + 2) / T := by
  obtain ⟨C, U₀, hC⟩ := L.selberg
  refine ⟨4 * |C|, max U₀ (exp 1), by positivity, le_max_right _ _, fun T X hT hTX => ?_⟩
  have hTe : exp 1 ≤ T := le_trans (le_max_right _ _) hT
  have hTU : U₀ ≤ T := le_trans (le_max_left _ _) hT
  have hTpos : 0 < T := lt_of_lt_of_le (exp_pos 1) hTe
  have hlogT : 1 ≤ Real.log T := by rw [← Real.log_exp 1]; exact Real.log_le_log (exp_pos 1) hTe
  have hll : 0 ≤ Real.log (Real.log T) := Real.log_nonneg hlogT
  set f : ℝ → ℝ := fun t => zetaS t ^ 2 / (t ^ 2 + 1 / 4) with hf
  have hcinv : Continuous fun t : ℝ => 1 / (t ^ 2 + 1 / 4) :=
    Continuous.div continuous_const (by fun_prop) fun t => by positivity
  have hint : ∀ a b, IntervalIntegrable f volume a b := by
    intro a b
    refine ((intervalIntegrable_zetaS_sq a b).mul_continuousOn hcinv.continuousOn).congr_uIoo fun t _ => ?_
    simp only [hf]; ring
  have hf0 : ∀ t, 0 ≤ f t := fun t => by simp only [hf]; positivity
  -- K with T 2^K ≥ X
  obtain ⟨K, hK⟩ := pow_unbounded_of_one_lt (X / T) (by norm_num : (1 : ℝ) < 2)
  have hXK : X ≤ T * 2 ^ K := by
    rw [div_lt_iff₀ hTpos] at hK; linarith
  have h1 : ∫ t in T..X, f t ≤ ∫ t in T..(T * 2 ^ K), f t :=
    intervalIntegral.integral_mono_interval le_rfl hTX hXK (Filter.Eventually.of_forall hf0) (hint _ _)
  have h2 : ∫ t in T..(T * 2 ^ K), f t = ∑ k ∈ range K, ∫ t in (T * 2 ^ k)..(T * 2 ^ (k + 1)), f t := by
    have := intervalIntegral.sum_integral_adjacent_intervals (a := fun k : ℕ => T * 2 ^ k) (n := K)
      (fun k _ => hint _ _)
    rw [this]; simp
  -- each piece
  have hpiece : ∀ k : ℕ, ∫ t in (T * 2 ^ k)..(T * 2 ^ (k + 1)), f t ≤
      2 * |C| / T * ((Real.log (Real.log T) + ((k : ℝ) + 1)) * (1 / 2 : ℝ) ^ k) := by
    intro k
    set a := T * 2 ^ k with ha
    set b := T * 2 ^ (k + 1) with hb
    have ha0 : 0 < a := by positivity
    have hab : a ≤ b := by rw [ha, hb, pow_succ]; nlinarith [pow_pos (two_pos : (0 : ℝ) < 2) k]
    have hb2 : b = 2 * a := by rw [ha, hb, pow_succ]; ring
    have hTa : T ≤ a := le_mul_of_one_le_right hTpos.le (one_le_pow₀ (by norm_num))
    have e1 : ∫ t in a..b, f t ≤ ∫ t in a..b, zetaS t ^ 2 / a ^ 2 := by
      apply intervalIntegral.integral_mono_on hab (hint _ _)
        ((intervalIntegrable_zetaS_sq a b).div_const _)
      intro t ht
      simp only [hf]
      apply div_le_div_of_nonneg_left (sq_nonneg _) (by positivity)
      have : a ^ 2 ≤ t ^ 2 := pow_le_pow_left₀ ha0.le ht.1 2
      linarith
    have e2 : ∫ t in a..b, zetaS t ^ 2 / a ^ 2 ≤ (∫ t in (0 : ℝ)..b, zetaS t ^ 2) / a ^ 2 := by
      rw [intervalIntegral.integral_div]
      apply div_le_div_of_nonneg_right _ (by positivity)
      exact intervalIntegral.integral_mono_interval ha0.le hab le_rfl
        (Filter.Eventually.of_forall fun _ => sq_nonneg _) (intervalIntegrable_zetaS_sq _ _)
    have hbU : U₀ ≤ b := le_trans hTU (le_trans hTa hab)
    have e3 := hC b hbU
    have hlb : Real.log (Real.log b) ≤ Real.log (Real.log T) + ((k : ℝ) + 1) := by
      have hbT : b = T * 2 ^ (k + 1) := rfl
      have hlb1 : Real.log b = Real.log T + (k + 1) * Real.log 2 := by
        rw [hbT, Real.log_mul hTpos.ne' (by positivity), Real.log_pow]; push_cast; ring
      have hl2 : Real.log 2 ≤ 1 := by
        have := Real.log_le_sub_one_of_pos (show (0 : ℝ) < 2 by norm_num); linarith
      have hl20 : 0 ≤ Real.log 2 := Real.log_nonneg (by norm_num)
      have hlb2 : Real.log b ≤ Real.log T * ((k : ℝ) + 2) := by
        rw [hlb1]
        have : ((k : ℝ) + 1) * Real.log 2 ≤ (k : ℝ) + 1 := by
          have h0 : (0 : ℝ) ≤ (k : ℝ) + 1 := by positivity
          nlinarith
        nlinarith
      have hlbpos : 0 < Real.log b := by rw [hlb1]; positivity
      have := Real.log_le_log hlbpos hlb2
      rw [Real.log_mul (show Real.log T ≠ 0 by linarith) (show ((k : ℝ) + 2) ≠ 0 by positivity)] at this
      have h3 : Real.log ((k : ℝ) + 2) ≤ (k : ℝ) + 1 := by
        have := Real.log_le_sub_one_of_pos (show (0 : ℝ) < (k : ℝ) + 2 by positivity); linarith
      linarith
    have hCabs : C * b * Real.log (Real.log b) ≤ |C| * b * (Real.log (Real.log T) + ((k : ℝ) + 1)) := by
      have hlbb : 0 ≤ Real.log (Real.log b) := by
        apply Real.log_nonneg
        rw [← Real.log_exp 1]; exact Real.log_le_log (exp_pos 1) (le_trans hTe (le_trans hTa hab))
      calc C * b * Real.log (Real.log b) ≤ |C| * b * Real.log (Real.log b) := by
            apply mul_le_mul_of_nonneg_right _ hlbb
            exact mul_le_mul_of_nonneg_right (le_abs_self C) (by positivity)
        _ ≤ |C| * b * (Real.log (Real.log T) + ((k : ℝ) + 1)) :=
            mul_le_mul_of_nonneg_left hlb (by positivity)
    calc ∫ t in a..b, f t ≤ (∫ t in (0 : ℝ)..b, zetaS t ^ 2) / a ^ 2 := e1.trans e2
      _ ≤ |C| * b * (Real.log (Real.log T) + ((k : ℝ) + 1)) / a ^ 2 := by
          apply div_le_div_of_nonneg_right (e3.trans hCabs) (by positivity)
      _ = 2 * |C| / T * ((Real.log (Real.log T) + ((k : ℝ) + 1)) * (1 / 2 : ℝ) ^ k) := by
          rw [hb2, ha]; field_simp
          rw [mul_assoc, ← mul_pow]; norm_num
  rw [h2] at h1
  have h3 := Finset.sum_le_sum fun k (_ : k ∈ range K) => hpiece k
  rw [← Finset.mul_sum] at h3
  have hs : ∑ k ∈ range K, (Real.log (Real.log T) + ((k : ℝ) + 1)) * (1 / 2 : ℝ) ^ k ≤
      2 * Real.log (Real.log T) + 4 := by
    have e : ∑ k ∈ range K, (Real.log (Real.log T) + ((k : ℝ) + 1)) * (1 / 2 : ℝ) ^ k =
        Real.log (Real.log T) * ∑ k ∈ range K, (1 / 2 : ℝ) ^ k + ∑ k ∈ range K, ((k : ℝ) + 1) * (1 / 2 : ℝ) ^ k := by
      rw [Finset.mul_sum, ← Finset.sum_add_distrib]
      exact Finset.sum_congr rfl fun k _ => by ring
    rw [e]
    have := sum_geom_le K
    have := sum_k_geom_le K
    nlinarith
  have hfin : 2 * |C| / T * (2 * Real.log (Real.log T) + 4) = 4 * |C| * (Real.log (Real.log T) + 2) / T := by
    field_simp; ring
  calc ∫ t in T..X, zetaS t ^ 2 / (t ^ 2 + 1 / 4) = ∫ t in T..X, f t := rfl
    _ ≤ 2 * |C| / T * ∑ k ∈ range K, (Real.log (Real.log T) + ((k : ℝ) + 1)) * (1 / 2 : ℝ) ^ k := h1.trans h3
    _ ≤ 2 * |C| / T * (2 * Real.log (Real.log T) + 4) := mul_le_mul_of_nonneg_left hs (by positivity)
    _ = 4 * |C| * (Real.log (Real.log T) + 2) / T := hfin

end

end LeanLi
